run_cpp_metabolomics_engine <- function(
    dataA,
    rt_window       = 10,
    alpha           = 0.70,
    top_k           = 15,
    soft_power      = 6,
    min_module_size = 10,
    use_igraph_infomap = FALSE   # TRUE = igraph (better), FALSE = C++ greedy
) {
  library(data.table)
  library(igraph)

  setDT(dataA); setorder(dataA, time)
  X <- as.matrix(dataA[, -(1:2)])
  X[!is.finite(X)] <- 0
  rv <- apply(X, 1, var)
  keep <- !is.na(rv) & rv > 1e-10
  X <- X[keep, ]; dataA <- dataA[keep, ]
  X <- scale(X, center=TRUE, scale=FALSE)
  n <- nrow(X)

  # ── Step 1: C++ sparse graph ─────────────────────────
  message("Building sparse graph...")
  res <- build_sparse_graph_parallel(
    X=X, rt=dataA$time, rt_window=rt_window,
    alpha=alpha, top_k=top_k
  )

  # ── Step 2: Soft threshold + C++ sparse TOM ──────────
  message("Computing TOM...")
  w_soft <- res$weight ^ soft_power

  tom_w  <- compute_sparse_tom(
    from   = as.integer(res$from),
    to     = as.integer(res$to),
    weight = w_soft,
    n      = n
  )

  # Remove near-zero TOM edges
  keep_edges <- tom_w > 1e-6
  from_k <- res$from[keep_edges]
  to_k   <- res$to[keep_edges]
  tom_k  <- tom_w[keep_edges]

  # ── Step 3: Cluster ───────────────────────────────────
  message("Clustering...")

  if(use_igraph_infomap){
    # Recommended: igraph's infomap (C-level, no R overhead)
    g <- igraph::graph_from_data_frame(
      data.frame(from=from_k, to=to_k, weight=tom_k),
      directed = FALSE,
      vertices = data.frame(name = seq_len(n))
    )
    g  <- igraph::simplify(g, edge.attr.comb="max")
    cl <- igraph::cluster_infomap(g, e.weights=igraph::E(g)$weight)
    modules <- igraph::membership(cl)

  } else {
    # Fallback: pure C++ greedy modularity
    modules <- cluster_sparse_greedy(
      from    = as.integer(from_k),
      to      = as.integer(to_k),
      weight  = tom_k,
      n       = n,
      max_iter = 20
    )
  }

  # Merge tiny modules
  mod_sizes <- table(modules)
  tiny      <- as.integer(names(mod_sizes[mod_sizes < min_module_size]))
  modules[modules %in% tiny] <- 0L

  dataA[, Module_RTclust := as.character(modules)]
  message(sprintf(
    "%d modules | median: %.0f | unassigned: %d",
    sum(as.integer(names(table(modules))) > 0),
    median(mod_sizes[mod_sizes >= min_module_size]),
    sum(modules == 0)
  ))

  return(dataA[, .(mz, time, Module_RTclust)])
}
