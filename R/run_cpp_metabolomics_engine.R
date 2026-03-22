run_cpp_metabolomics_engine <- function(
    dataA,
    rt_window       = 10,
    alpha           = 0.7,
    top_k           = 20,
    min_module_size = 10,
    soft_power      = 2      # additional soft threshold on top of C++ r^2
) {
  library(data.table)
  library(igraph)
  library(dynamicTreeCut)

  setDT(dataA)
  setorder(dataA, time)
  X <- as.matrix(dataA[, -(1:2)])


  # Clean X before variance check
  X[is.na(X)]      <- 0
  X[!is.finite(X)] <- 0

  row_vars <- apply(X, 1, var)
  keep     <- !is.na(row_vars) & row_vars > 1e-10

  if(any(!keep)){
    message("Removing ", sum(!keep), " zero-variance features")
  }
  X     <- X[keep, ]
  dataA <- dataA[keep, ]
  # AFTER — three-part keep condition:
  #   1. variance is not NA
  #   2. variance is above threshold
  #   3. no infinite values in the row
  row_vars <- apply(X, 1, var, na.rm = TRUE)
  keep     <- !is.na(row_vars) & row_vars > 1e-10 & apply(X, 1, function(r) all(is.finite(r)))

  if(any(!keep)){
    message("Removing ", sum(!keep), " zero-variance or non-finite features")
  }

  X     <- X[keep, ]
  dataA <- dataA[keep, ]

  X <- scale(X, center = TRUE, scale = FALSE)

  # ── Build sparse adjacency ──────────────────────────────────────────────────
  res <- build_sparse_graph_parallel(
    X         = X,
    rt        = dataA$time,
    rt_window = rt_window,
    alpha     = alpha,
    top_k     = top_k
  )

  # ── Convert to symmetric adjacency matrix ──────────────────────────────────
  # WGCNA works on a full symmetric matrix — we approximate TOM by
  # symmetrising: w_ij = max(w_ij, w_ji) so mutual edges get full weight
  n <- nrow(X)
  adj <- Matrix::sparseMatrix(
    i    = res$from,
    j    = res$to,
    x    = res$weight^soft_power,   # additional soft thresholding
    dims = c(n, n),
    symmetric = FALSE
  )
  # Symmetrise
  adj <- pmax(adj, Matrix::t(adj))

  # ── TOM-like similarity (optional but improves biological relevance) ────────
  # Full TOM is O(n^3) — this is a sparse approximation:
  # shared neighbor score between i and j
  # TOM_ij = (adj_ij + shared_neighbors) / (min(degree_i, degree_j) + 1 - adj_ij)
  # For large datasets skip this and use adj directly
  use_tom <- TRUE #n < 5000
  if(use_tom){
    message("Computing sparse TOM approximation...")
    # Shared neighbor count via matrix product
    shared <- as.matrix(adj %*% adj)
    deg    <- Matrix::rowSums(adj)
    tom    <- matrix(0, n, n)
    for(i in 1:n){
      for(j in i:n){
        if(adj[i,j] == 0 && shared[i,j] == 0) next
        denom      <- min(deg[i], deg[j]) - adj[i,j] + 1
        tom[i,j]   <- (adj[i,j] + shared[i,j]) / denom
        tom[j,i]   <- tom[i,j]
      }
    }
    dissim <- 1 - tom
  } else {
    message("Skipping TOM (n > 5000), using correlation distance directly")
    dissim <- 1 - as.matrix(adj)
  }

  # ── Hierarchical clustering + dynamic tree cut ─────────────────────────────
  # This is the WGCNA approach — no resolution parameter needed
  hc <- hclust(as.dist(dissim), method = "average")

  modules <- dynamicTreeCut::cutreeDynamic(
    dendro        = hc,
    distM         = dissim,
    deepSplit     = 2,          # 0-4: higher = more/smaller modules
    minClusterSize = min_module_size,
    method        = "hybrid"    # uses both tree shape and distance matrix
  )

  # Module 0 = unassigned — label them individually
  modules[modules == 0] <- max(modules) + seq_len(sum(modules == 0))

  dataA[, Module_RTclust := as.character(modules)]

  message(sprintf(
    "Found %d modules | median size: %.0f | range: %d–%d",
    length(unique(modules)),
    median(table(modules)),
    min(table(modules)),
    max(table(modules))
  ))

  return(dataA[, .(mz, time, Module_RTclust)])
}
