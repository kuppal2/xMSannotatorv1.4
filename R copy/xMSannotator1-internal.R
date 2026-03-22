.load_db <-
function(db_name, customDB, adduct_names, adduct_table,
                     customIDs, biofluid.location, origin, status,
                     HMDBselect) 
{
  
  if (db_name == "HMDB") {
    
    data(hmdbAllinf)
    hmdbAllinf$Name <- gsub(hmdbAllinf$Name, pattern="[\\\"\']", replacement="")
    hmdbAllinfv3.6  <- hmdbAllinf
    rm(hmdbAllinf)
    suppressWarnings(rm(hmdbAllinf, envir=.GlobalEnv))
    
    suppressWarnings({
      if (is.na(customIDs)[1L]) {
        customIDs <- hmdbAllinfv3.6[, c(1L, 20L)]
        
        .filter_hmdb <- function(df, bf, orig, stat, sel) {
          # FIX v2.1.2: gregexpr() returns a *list* so `gregexpr(...) == 1L`
          # compares the whole list object to a scalar — it evaluates to a
          # single logical, not a per-row vector, making which() return either
          # everything or nothing. This caused ~90% of HMDB entries to be
          # dropped (3249 vs 35186 DB matches). Fix: extract the first match
          # position from each list element via sapply(..., `[[`, 1), exactly
          # mirroring the old v2.0.1 pattern of `gres <- gregexpr(...)` then
          # `which(gres == 1)` which R auto-unlists during comparison.
          .greg_which <- function(vec, pat) {
            which(sapply(gregexpr(vec, pattern=pat, ignore.case=TRUE), `[[`, 1L) == 1L)
          }
          idx_bf   <- if (!is.na(bf))   .greg_which(df$BioFluidLocation, bf)   else seq_len(nrow(df))
          idx_orig <- if (!is.na(orig)) .greg_which(df$Origin,           orig)  else seq_len(nrow(df))
          idx_stat <- if (!is.na(stat)) .greg_which(df$HMDBStatus,       stat)  else seq_len(nrow(df))
          if (sel == "intersect")
            Reduce(intersect, list(idx_bf, idx_orig, idx_stat))
          else
            unique(c(idx_bf, idx_orig, idx_stat))
        }
        
        gres     <- .filter_hmdb(hmdbAllinfv3.6, biofluid.location, origin,
                                 status, HMDBselect)
        customIDs <- hmdbAllinfv3.6[gres, c(1L, 20L)]
      }
    })
    
    suppressWarnings(rm(hmdbAllinfv3.6, envir=.GlobalEnv))
    
    data(hmdbCompMZ)
    hmdbCompMZ$mz   <- round(as.numeric(as.character(hmdbCompMZ$mz)), 5L)
    hmdbCompMZ$Name <- gsub(hmdbCompMZ$Name, pattern="[\\\"\']", replacement="")
    
    suppressWarnings({
      if (!is.na(customIDs)[1L]) {
        customIDs  <- unique(customIDs)
        hmdbCompMZ <- hmdbCompMZ[hmdbCompMZ$HMDBID %in% customIDs[, 1L], ]
        hmdbCompMZ <- hmdbCompMZ[hmdbCompMZ$Adduct %in% adduct_names, ]
      }
    })
    
    hmdbCompMZ <- hmdbCompMZ[hmdbCompMZ$Adduct %in% adduct_names, ]
    chemCompMZ <- hmdbCompMZ
    print("Dimension of precomputed HMDB m/z database"); print(dim(chemCompMZ))
    try(rm(hmdbCompMZ), silent=TRUE)
    try(rm(hmdbCompMZ, envir=.GlobalEnv), silent=TRUE)
    try(rm(hmdbAllinf, envir=.GlobalEnv), silent=TRUE)
    
  } else if (db_name == "KEGG") {
    
    data(keggCompMZ)
    keggCompMZ$mz   <- round(as.numeric(as.character(keggCompMZ$mz)), 5L)
    keggCompMZ$Name <- gsub(keggCompMZ$Name, pattern="[\\\"\']", replacement="")
    suppressWarnings({
      if (!is.na(customIDs)[1L]) {
        customIDs  <- unique(customIDs)
        keggCompMZ <- keggCompMZ[keggCompMZ$KEGGID %in% customIDs[, 1L], ]
        keggCompMZ <- keggCompMZ[keggCompMZ$Adduct %in% adduct_names, ]
      }
    })
    keggCompMZ <- keggCompMZ[keggCompMZ$Adduct %in% adduct_names, ]
    chemCompMZ <- keggCompMZ
    print("Dimension of precomputed KEGG m/z database"); print(dim(chemCompMZ))
    try(rm(keggCompMZ), silent=TRUE)
    try(rm(keggCompMZ, envir=.GlobalEnv), silent=TRUE)
    
  } else if (db_name == "LipidMaps") {
    
    data(lipidmapsCompMZ)
    lipidmapsCompMZ <- lipidmapsCompMZ[lipidmapsCompMZ$Adduct %in% adduct_names, ]
    lipidmapsCompMZ$mz   <- round(as.numeric(as.character(lipidmapsCompMZ$mz)), 5L)
    lipidmapsCompMZ$Name <- gsub(lipidmapsCompMZ$Name, pattern="[\\\"\']", replacement="")
    chemCompMZ <- lipidmapsCompMZ
    print("Dimension of precomputed LipidMaps m/z database"); print(dim(chemCompMZ))
    try(rm(lipidmapsCompMZ), silent=TRUE)
    try(rm(lipidmapsCompMZ, envir=.GlobalEnv), silent=TRUE)
    
  } else if (db_name == "T3DB") {
    
    data(t3dbCompMZ)
    t3dbCompMZ <- t3dbCompMZ[t3dbCompMZ$Adduct %in% adduct_names, ]
    t3dbCompMZ$mz   <- round(as.numeric(as.character(t3dbCompMZ$mz)), 5L)
    t3dbCompMZ$Name <- gsub(t3dbCompMZ$Name, pattern="[\\\"\']", replacement="")
    chemCompMZ <- t3dbCompMZ
    print("Dimension of precomputed T3DB m/z database"); print(dim(chemCompMZ))
    try(rm(t3dbCompMZ), silent=TRUE)
    try(rm(t3dbCompMZ, envir=.GlobalEnv), silent=TRUE)
    
  } else if (db_name == "Custom") {
    
    data(adduct_table)
    inputmassmat <- customDB
    if (length(which(duplicated(inputmassmat[, 1L]))) > 0L)
      inputmassmat <- inputmassmat[-which(duplicated(inputmassmat[, 1L])), ]
    
    mz_search_list<-get_mz_by_monoisotopicmass(inputmassmat,queryadductlist=as.character(adduct_names),
                                               adduct_table=adduct_table)
    
    
    
    chemCompMZ <- as.data.frame(mz_search_list)
    save(chemCompMZ, file="mz_search_list.Rda")
    try(rm(customDB), silent=TRUE)
    try(rm(customDB, envir=.GlobalEnv), silent=TRUE)
    
  } else {
    stop("db_name should be: KEGG, HMDB, T3DB, LipidMaps, or Custom")
  }
  
  return(chemCompMZ)
}
.xms_fast_cluster <-
function(dataA, global_cor, mzid, corthresh,
                              max_diff_rt, minclustsize,
                              num_nodes, outloc) {
  
  n <- nrow(dataA)
  cat(sprintf("[FastClust] %d features, corthresh=%.2f, max_diff_rt=%.1f\n",
              n, corthresh, max_diff_rt))
  
  # ------------------------------------------------------------------
  # 1. BUILD SPARSE EDGE LIST
  # Upper-triangle pairs where cor >= corthresh AND |RT| <= max_diff_rt.
  # Row-by-row processing keeps peak RAM at O(n) per row, not O(n^2).
  # ------------------------------------------------------------------
  rt_vec    <- dataA$time
  edge_from <- integer(0)
  edge_to   <- integer(0)
  edge_w    <- numeric(0)
  
  for (i in seq_len(n - 1L)) {
    row_i   <- global_cor[i, (i + 1L):n]
    rt_diff <- abs(rt_vec[(i + 1L):n] - rt_vec[i])
    keep    <- which(row_i >= corthresh & rt_diff <= max_diff_rt)
    if (length(keep) == 0L) next
    j_idx       <- keep + i
    edge_from   <- c(edge_from, rep(i, length(keep)))
    edge_to     <- c(edge_to,   j_idx)
    edge_w      <- c(edge_w,    row_i[keep])
  }
  
  cat(sprintf("[FastClust] Sparse graph: %d edges (%.3f%% of all pairs)\n",
              length(edge_from),
              100 * length(edge_from) / max(1, n * (n - 1) / 2)))
  
  # ------------------------------------------------------------------
  # 2. BUILD igraph AND DETECT COMMUNITIES
  # ------------------------------------------------------------------
  if (length(edge_from) == 0L) {
    warning("[FastClust] No edges at corthresh=", corthresh,
            ". Assigning all features to singleton modules.")
    mod_list <- seq_len(n)
  } else {
    # ------------------------------------------------------------------
    # FIX: graph_from_edgelist() creates vertices 1..max(edge_to) which
    # can exceed n (e.g. 11838 > 11812) when the last few rows have no
    # edges and the highest index in an edge exceeds the feature count.
    # Solution: create an empty graph with exactly n vertices, then add
    # edges. Vertex i always corresponds to dataA row i — no reindexing.
    # ------------------------------------------------------------------
    g <- igraph::make_empty_graph(n=n, directed=FALSE)
    igraph::V(g)$name <- mzid          # vertex i == dataA row i == mzid[i]
    
    edge_mat <- matrix(c(edge_from, edge_to), ncol=2L)
    g <- igraph::add_edges(g, t(edge_mat), weight=edge_w)
    
    # Community detection: Leiden -> Louvain -> label propagation
    comm <- tryCatch({
      if (requireNamespace("leidenbase", quietly=TRUE)) {
        igraph::cluster_leiden(g, objective_function="modularity",
                               weights=igraph::E(g)$weight, n_iterations=5L)
      } else {
        igraph::cluster_louvain(g, weights=igraph::E(g)$weight)
      }
    }, error=function(e) {
      cat("[FastClust] Leiden/Louvain failed, using label propagation.\n")
      igraph::cluster_label_prop(g, weights=igraph::E(g)$weight)
    })
    
    # membership[i] is now directly the community for dataA row i
    mod_list <- as.integer(igraph::membership(comm))
    
    cat(sprintf("[FastClust] %d communities (before RT sub-grouping)\n",
                length(unique(mod_list))))
    
    # Absorb tiny communities (< floor(minclustsize/2)) into best neighbour
    min_size <- max(2L, minclustsize %/% 2L)
    comm_sz  <- table(mod_list)
    tiny     <- as.integer(names(comm_sz[comm_sz < min_size]))
    if (length(tiny) > 0L) {
      for (tc in tiny) {
        tc_rows <- which(mod_list == tc)
        for (fi in tc_rows) {
          nbr_v <- as.integer(igraph::neighbors(g, fi))  # vertex == row index
          nbr_v <- nbr_v[mod_list[nbr_v] != tc]
          if (length(nbr_v) > 0L)
            mod_list[fi] <- as.integer(
              names(which.max(table(mod_list[nbr_v]))))
        }
      }
    }
  }
  
  save(mod_list, file=file.path(outloc, "mod_list.Rda"))
  
  # ------------------------------------------------------------------
  # 3. RETURN COMMUNITIES AS Module_RTclust
  #
  # Graph edges are built with the constraint |RT| <= max_diff_rt, so
  # every feature within a Louvain community already co-elutes within
  # max_diff_rt by construction.  No further RT sub-grouping is needed
  # or beneficial — splitting communities further scatters adducts of
  # the same compound into different Module_RTclust groups and prevents
  # confidence-level-2 scoring.
  # ------------------------------------------------------------------
  levelA_res <- dataA[, c("mz", "time"), drop=FALSE]
  levelA_res$Module_RTclust <- as.character(mod_list)
  
  cat(sprintf("[FastClust] Final: %d Module_RTclust groups\n",
              length(unique(levelA_res$Module_RTclust))))
  levelA_res
}
.xms_fingerprint <-
function(...) {
  
  args <- list(...)
  
  normalize <- function(x) {
    if (is.null(x)) return("NULL")
    if (length(x) == 0) return("EMPTY")
    if (all(is.na(x))) return("NA")
    paste(as.character(x), collapse=",")
  }
  
  parts <- vapply(args, normalize, character(1))
  
  paste(parts, collapse="|")
}
.xms_read_fp <-
function(fp_file) {
  if (file.exists(fp_file)) readLines(fp_file, warn = FALSE)[1L] else ""
}
.xms_write_fp <-
function(fp_file, fp) {
  writeLines(as.character(fp), fp_file)
}
