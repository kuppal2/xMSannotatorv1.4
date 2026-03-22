run_cpp_metabolomics_engine <- function(
    dataA,
    rt_window          = 60,
    alpha              = 0.70,
    top_k              = 100,
    soft_power         = 8,
    min_module_size    = 10,
    merge_cut_height   = 0.20,
    use_igraph_infomap = TRUE,
    seed               = 555L      # single seed controls all stochastic steps
) {
  library(data.table)
  library(igraph)

  # ============================================================
  # FIX 1 — Force single-threaded C++ execution
  # parallelReduce join() order is OS-scheduled — different thread
  # finish orders produce different edge-list orderings each run
  # even with identical input. Single thread = deterministic order.
  # ============================================================
  old_omp   <- Sys.getenv("OMP_NUM_THREADS")
  old_tbb   <- Sys.getenv("RCPP_PARALLEL_NUM_THREADS")
  on.exit({
    Sys.setenv(OMP_NUM_THREADS            = old_omp)
    Sys.setenv(RCPP_PARALLEL_NUM_THREADS  = old_tbb)
  }, add = TRUE)
  Sys.setenv(OMP_NUM_THREADS           = "1")
  Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "1")

  # ============================================================
  # MERGE CLOSE MODULES (eigengene-based, WGCNA-style)
  # Defined inside the function to keep it self-contained.
  # ============================================================
  merge_close_modules <- function(
    dataA,
    X,
    cut_height      = 0.25,
    min_module_size = 10,
    verbose         = TRUE
  ) {
    library(WGCNA)

    modules    <- as.integer(dataA$Module_RTclust)
    module_ids <- sort(unique(modules[modules > 0]))
    n_before   <- length(module_ids)

    if(n_before < 2){
      message("Only ", n_before, " module(s) — nothing to merge.")
      return(dataA)
    }

    if(verbose) message(sprintf("Merging close modules (cut_height=%.2f)...", cut_height))

    eigengenes    <- matrix(NA, nrow = ncol(X), ncol = length(module_ids))
    colnames(eigengenes) <- paste0("ME", module_ids)
    valid_modules <- integer(0)

    for(k in seq_along(module_ids)){
      mid <- module_ids[k]
      idx <- which(modules == mid)

      if(length(idx) < 3){
        eigengenes[, k] <- colMeans(X[idx, , drop = FALSE])
        valid_modules    <- c(valid_modules, mid)
        next
      }

      tryCatch({
        sv  <- svd(t(X[idx, ]), nu = 1, nv = 0)
        eg  <- sv$u[, 1]
        # FIX 4 — SVD sign: use complete.obs to avoid NaN order sensitivity
        avg_cor <- mean(cor(eg, t(X[idx, ]),
                            use = "complete.obs"), na.rm = TRUE)
        if(!is.na(avg_cor) && avg_cor < 0) eg <- -eg
        eigengenes[, k] <- eg
        valid_modules    <- c(valid_modules, mid)
      }, error = function(e){
        message(sprintf("  Eigengene failed for module %d: %s", mid, e$message))
      })
    }

    keep_cols  <- paste0("ME", valid_modules)
    eigengenes <- eigengenes[, keep_cols, drop = FALSE]
    module_ids <- valid_modules

    if(ncol(eigengenes) < 2){
      message("Too few valid eigengenes to merge.")
      return(dataA)
    }

    eg_cor  <- cor(eigengenes, use = "pairwise.complete.obs")
    eg_dist <- as.dist(1 - eg_cor)
    eg_hc   <- hclust(eg_dist, method = "average")

    merge_groups   <- cutree(eg_hc, h = cut_height)
    n_merge_groups <- length(unique(merge_groups))
    n_merges       <- n_before - n_merge_groups

    if(verbose){
      message(sprintf("  Modules before merge: %d", n_before))
      message(sprintf("  Merge groups formed:  %d", n_merge_groups))
      message(sprintf("  Modules merged:       %d", n_merges))
    }

    if(n_merges == 0){
      message("  No modules close enough to merge at cut_height=", cut_height)
      return(dataA)
    }

    merge_map <- integer(max(module_ids) + 1)
    for(g in unique(merge_groups)){
      members   <- module_ids[merge_groups == g]
      target_id <- min(members)
      for(m in members) merge_map[m] <- target_id
      if(verbose && length(members) > 1)
        message(sprintf("  Merging modules {%s} → module %d",
                        paste(members, collapse = ","), target_id))
    }
    merge_map[1] <- 0L   # module 0 stays unassigned

    new_modules <- ifelse(modules > 0, merge_map[modules], 0L)

    mod_sizes <- table(new_modules)
    tiny      <- as.integer(names(
      mod_sizes[mod_sizes < min_module_size & names(mod_sizes) != "0"]
    ))
    if(length(tiny) > 0){
      if(verbose) message(sprintf(
        "  Re-assigning %d newly-small modules to grey (0)", length(tiny)))
      new_modules[new_modules %in% tiny] <- 0L
    }

    unique_mods <- sort(unique(new_modules[new_modules > 0]))
    relabel     <- setNames(seq_along(unique_mods), unique_mods)
    final       <- ifelse(new_modules > 0,
                          relabel[as.character(new_modules)], 0L)

    dataA[, Module_RTclust := as.character(final)]

    if(verbose){
      mod_tab <- table(final)
      message(sprintf(
        "  Final: %d modules | median size: %.0f | unassigned: %d",
        sum(as.integer(names(mod_tab)) > 0),
        median(mod_tab[names(mod_tab) != "0"]),
        sum(final == 0)
      ))
    }

    return(dataA)
  }

  # ============================================================
  # DATA PREPARATION
  # ============================================================
  setDT(dataA)
  setorder(dataA, time)

  X <- as.matrix(dataA[, -(1:2)])
  X[!is.finite(X)] <- 0

  rv   <- apply(X, 1, var)
  keep <- !is.na(rv) & rv > 1e-10
  X     <- X[keep, ]
  dataA <- dataA[keep, ]
  X     <- scale(X, center = TRUE, scale = FALSE)
  n     <- nrow(X)

  # ============================================================
  # STEP 1 — Sparse graph (C++, now single-threaded = deterministic)
  # FIX 1+2: OMP_NUM_THREADS=1 set above ensures parallelReduce
  # join() always concatenates in the same worker order.
  # ============================================================
  message("Building sparse graph...")
  res <- build_sparse_graph_parallel(
    X         = X,
    rt        = dataA$time,
    rt_window = rt_window,
    alpha     = alpha,
    top_k     = top_k
  )

  # ============================================================
  # FIX 2 — Sort edges by (from, to) after parallel call
  # Belt-and-suspenders: even with single thread, sort guarantees
  # the edge list is in canonical order regardless of C++ internals.
  # ============================================================
  edge_order <- order(res$from, res$to)
  res$from   <- res$from[edge_order]
  res$to     <- res$to[edge_order]
  res$weight <- res$weight[edge_order]

  # ============================================================
  # STEP 2 — Sparse TOM (C++, single-threaded)
  # ============================================================
  message("Computing TOM...")
  w_soft <- res$weight ^ soft_power

  tom_w <- compute_sparse_tom(
    from   = as.integer(res$from),
    to     = as.integer(res$to),
    weight = w_soft,
    n      = n
  )

  # Sort TOM output to match canonical edge order
  tom_order <- order(res$from, res$to)
  from_k    <- res$from[tom_order]
  to_k      <- res$to[tom_order]
  tom_k     <- tom_w[tom_order]

  # Remove near-zero TOM edges
  keep_edges <- tom_k > 1e-6
  from_k     <- from_k[keep_edges]
  to_k       <- to_k[keep_edges]
  tom_k      <- tom_k[keep_edges]

  # ============================================================
  # STEP 3 — Cluster
  # FIX 3 — set.seed placed IMMEDIATELY before each stochastic call.
  # The seed placed before build_sparse_graph_parallel in the
  # original code was consumed by intermediate operations before
  # cluster_infomap ran.
  # ============================================================
  message("Clustering...")

  if(use_igraph_infomap){

    g <- igraph::graph_from_data_frame(
      data.frame(from = from_k, to = to_k, weight = tom_k),
      directed = FALSE,
      vertices = data.frame(name = seq_len(n))
    )
    g <- igraph::simplify(g, edge.attr.comb = "max")

    # FIX 3 — seed immediately before stochastic call
    set.seed(seed)
    cl      <- igraph::cluster_infomap(g, e.weights = igraph::E(g)$weight)
    modules <- igraph::membership(cl)

  } else {

    # FIX 3 — seed immediately before stochastic call
    set.seed(seed)
    modules <- cluster_sparse_greedy(
      from     = as.integer(from_k),
      to       = as.integer(to_k),
      weight   = tom_k,
      n        = n,
      max_iter = 20
    )
  }

  # ============================================================
  # STEP 4 — Remove tiny modules
  # ============================================================
  mod_sizes <- table(modules)
  tiny      <- as.integer(names(mod_sizes[mod_sizes < min_module_size]))
  modules[modules %in% tiny] <- 0L
  dataA[, Module_RTclust := as.character(modules)]

  # ============================================================
  # STEP 5 — Merge close modules (eigengene-based)
  # SVD inside is deterministic given identical input — no seed needed.
  # ============================================================
  if(merge_cut_height > 0){
    dataA <- merge_close_modules(
      dataA           = dataA,
      X               = X,
      cut_height      = merge_cut_height,
      min_module_size = min_module_size,
      verbose         = TRUE
    )
  }

  # Re-read modules from dataA after merge step
  modules_final <- as.integer(dataA$Module_RTclust)
  mod_tab_final <- table(modules_final)

  message(sprintf(
    "Done: %d modules | median: %.0f | unassigned: %d",
    sum(as.integer(names(mod_tab_final)) > 0),
    median(mod_tab_final[names(mod_tab_final) != "0"]),
    sum(modules_final == 0)
  ))

  return(dataA[, .(mz, time, Module_RTclust)])
}
