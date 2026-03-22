run_cpp_metabolomics_engine <- function(
    dataA,
    rt_window       = 1000,
    alpha           = 0.70,
    top_k           = 100,
    soft_power      = 6,
    min_module_size = 10,
    use_igraph_infomap = FALSE   # TRUE = igraph (better), FALSE = C++ greedy
) {
  library(data.table)
  library(igraph)

  merge_close_modules <- function(
    dataA,          # data.table with Module_RTclust column
    X,              # feature × sample matrix (mean-centered, same row order as dataA)
    cut_height = 0.25,   # merge if eigengene correlation > (1 - cut_height)
    # 0.25 = merge modules with r > 0.75 — WGCNA default
    min_module_size = 10,
    verbose = TRUE
  ) {
    library(data.table)
    library(WGCNA)

    modules     <- as.integer(dataA$Module_RTclust)
    module_ids  <- sort(unique(modules[modules > 0]))
    n_before    <- length(module_ids)

    if(n_before < 2){
      message("Only ", n_before, " module(s) — nothing to merge.")
      return(dataA)
    }

    if(verbose) message(sprintf("Merging close modules (cut_height=%.2f)...", cut_height))

    # ── Compute module eigengenes ──────────────────────────────
    # Eigengene = first principal component of each module's
    # feature expression matrix (samples × features within module)
    eigengenes <- matrix(NA, nrow = ncol(X), ncol = length(module_ids))
    colnames(eigengenes) <- paste0("ME", module_ids)
    valid_modules <- integer(0)

    for(k in seq_along(module_ids)){
      mid  <- module_ids[k]
      idx  <- which(modules == mid)

      if(length(idx) < 3){
        # Too few features for PCA — use mean as fallback eigengene
        eigengenes[, k] <- colMeans(X[idx, , drop=FALSE])
        valid_modules    <- c(valid_modules, mid)
        next
      }

      tryCatch({
        # SVD — first right singular vector = first PC across samples
        sv <- svd(t(X[idx, ]), nu=1, nv=0)
        eg <- sv$u[,1]

        # WGCNA sign convention: eigengene positively correlated
        # with majority of member features
        avg_cor <- mean(cor(eg, t(X[idx,])))
        if(avg_cor < 0) eg <- -eg

        eigengenes[, k] <- eg
        valid_modules    <- c(valid_modules, mid)

      }, error = function(e){
        message(sprintf("  Eigengene failed for module %d: %s", mid, e$message))
      })
    }

    # Keep only valid eigengenes
    keep_cols      <- paste0("ME", valid_modules)
    eigengenes     <- eigengenes[, keep_cols, drop=FALSE]
    module_ids     <- valid_modules

    if(ncol(eigengenes) < 2){
      message("Too few valid eigengenes to merge.")
      return(dataA)
    }

    # ── Eigengene correlation matrix ──────────────────────────
    eg_cor <- cor(eigengenes, use="pairwise.complete.obs")
    # 1 - correlation = dissimilarity (same scale as cut_height)
    eg_dist <- as.dist(1 - eg_cor)

    # ── Hierarchical clustering of eigengenes ─────────────────
    eg_hc <- hclust(eg_dist, method = "average")

    # ── Cut at cut_height to form merge groups ────────────────
    merge_groups <- cutree(eg_hc, h = cut_height)

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

    # ── Apply merge map ───────────────────────────────────────
    # Each merge group gets the ID of its smallest-numbered module
    # (matches WGCNA convention of keeping the "dominant" module)
    merge_map <- integer(max(module_ids) + 1)

    for(g in unique(merge_groups)){
      members    <- module_ids[merge_groups == g]
      target_id  <- min(members)   # keep lowest-numbered module ID
      for(m in members) merge_map[m] <- target_id

      if(verbose && length(members) > 1){
        message(sprintf("  Merging modules {%s} → module %d",
                        paste(members, collapse=","), target_id))
      }
    }

    # Unassigned features (module 0) stay at 0
    merge_map[0 + 1] <- 0   # index 1 = module 0

    # Apply
    new_modules <- ifelse(
      modules > 0,
      merge_map[modules],
      0L
    )

    # Re-check minimum size after merging
    mod_sizes <- table(new_modules)
    tiny      <- as.integer(names(mod_sizes[mod_sizes < min_module_size & names(mod_sizes) != "0"]))

    if(length(tiny) > 0){
      if(verbose) message(sprintf("  Re-assigning %d newly-small modules to grey (0)",
                                  length(tiny)))
      new_modules[new_modules %in% tiny] <- 0L
    }

    # Relabel to consecutive integers (cleaner output)
    unique_mods <- sort(unique(new_modules[new_modules > 0]))
    relabel     <- setNames(seq_along(unique_mods), unique_mods)
    final       <- ifelse(new_modules > 0, relabel[as.character(new_modules)], 0L)

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

  # ── Merge close modules ← THE MISSING STEP ────────────────
  if(merge_cut_height > 0){
    dataA <- merge_close_modules(
      dataA          = dataA,
      X              = X,
      cut_height     = merge_cut_height,
      min_module_size = min_module_size,
      verbose        = TRUE
    )
  }

  message(sprintf(
    "%d modules | median: %.0f | unassigned: %d",
    sum(as.integer(names(table(modules))) > 0),
    median(mod_sizes[mod_sizes >= min_module_size]),
    sum(modules == 0)
  ))

  return(dataA[, .(mz, time, Module_RTclust)])
}
