get_peak_blocks_modulesvhclust <-
function(
    dataA, simmat, mycl_metabs,
    adjacencyfromsimilarity = FALSE,
    outloc,
    column.rm.index = NA,
    deepsplit    = 2,
    minclustsize = 20,
    cutheight    = 0.2,
    cormethod    = "spearman",
    networktype  = "unsigned",
    num_nodes    = 2,
    step1log2scale = TRUE)
{
  # ------------------------------------------------------------------
  # WGCNA-based module detection.
  # Returns dataA (mz + time + sample columns) with Module_RTclust set
  # to the raw integer module label.
  # ------------------------------------------------------------------
  
  dataA<-as.data.frame(dataA)
  cnames    <- colnames(dataA)
  cnames[1] <- "mz"
  cnames[2] <- "time"
  colnames(dataA) <- cnames
  
  data_mzrt <- dataA[, 1:2, drop=FALSE]
  
  if (!is.na(column.rm.index)[1])
    dataA <- dataA[, -column.rm.index, drop=FALSE]
  
  feat_inf <- paste(dataA[, 1], dataA[, 2], sep="_")
  dataA    <- dataA[, -c(1:2), drop=FALSE]
  data_m   <- t(dataA)
  
  colnames(data_m) <- feat_inf
  
  if (isTRUE(step1log2scale))
    data_m <- 2^data_m
  
  # ------------------------------------------------------------------
  # SOFT-THRESHOLD SELECTION
  # ------------------------------------------------------------------
  powers    <- c(1:10, seq(12, 20, 2))
  set.seed(555)
  sft       <- try(WGCNA::pickSoftThreshold(data=data_m, dataIsExpr=TRUE,
                                     powerVector=powers, verbose=0),
                   silent=TRUE)
  power_val <- if (inherits(sft, "try-error") || is.na(sft$powerEstimate))
    6L else sft$powerEstimate
  
  # ------------------------------------------------------------------
  # BLOCKWISE MODULES (primary path)
  # Falls back to stepwise TOM + flashClust if blockwiseModules fails.
  # ------------------------------------------------------------------
  set.seed(555)
  net <- try(WGCNA::blockwiseModules(
    datExpr              = data_m,
    checkMissingData     = FALSE,
    blocks               = mycl_metabs,
    maxBlockSize         = 500,
    blockSizePenaltyPower= 5,
    randomSeed           = 12345,
    loadTOM              = FALSE,
    corType              = "pearson",
    maxPOutliers         = 1,
    quickCor             = 1,
    pearsonFallback      = "individual",
    cosineCorrelation    = FALSE,
    power                = power_val,
    networkType          = "unsigned",
    TOMType              = "unsigned",
    TOMDenom             = "min",
    saveTOMs             = FALSE,
    deepSplit            = deepsplit,
    minModuleSize        = minclustsize,
    pamStage             = TRUE,
    pamRespectsDendro    = FALSE,
    reassignThreshold    = 1e-06,
    minCoreKME           = 0.5,
    minCoreKMESize       = minclustsize / 3,
    minKMEtoStay         = 0.3,
    mergeCutHeight       = cutheight,
    impute               = TRUE,
    trapErrors           = FALSE,
    numericLabels        = FALSE,
    nThreads             = num_nodes,
    verbose              = 0,
    indent               = 0),
    silent=TRUE)
  
  if (!inherits(net, "try-error")) {
    
    mod_list <- as.integer(as.factor(unlist(net$colors)))
    
  } else {
    
    # Stepwise fallback: adjacency -> TOM -> flashClust -> cutreeDynamic
    if (!adjacencyfromsimilarity) {
      set.seed(555)
      sft2 <- try(WGCNA::pickSoftThreshold(data=data_m, dataIsExpr=TRUE,
                                    powerVector=powers, verbose=0), silent=TRUE)
      power_val <- if (inherits(sft2, "try-error") || is.na(sft2$powerEstimate))
        6L else sft2$powerEstimate
      ADJ <- if (cormethod == "pearson")
        WGCNA::adjacency(datExpr=data_m, type=networktype,
                  power=power_val, corOptions="use = 'p'")
      else
        WGCNA::adjacency(datExpr=data_m, type=networktype,
                  power=power_val,
                  corOptions="use = 'p', method = 'spearman'")
    } else {
      sft2 <- try(WGCNA::pickSoftThreshold.fromSimilarity(
        similarity=simmat, powerVector=powers, verbose=0),
        silent=TRUE)
      power_val <- if (inherits(sft2, "try-error") || is.na(sft2$powerEstimate))
        6L else sft2$powerEstimate
      ADJ <- WGCNA::adjacency.fromSimilarity(similarity=simmat,
                                      power=power_val, type=networktype)
    }
    
    dup_idx <- which(duplicated(rownames(ADJ)))
    if (length(dup_idx) > 0)
      ADJ <- ADJ[-dup_idx, -dup_idx, drop=FALSE]
    
    TOM  <- WGCNA::TOMdist(ADJ)
    set.seed(555)
    hier <- flashClust::flashClust(as.dist(TOM), method="complete")
    set.seed(555)
    clust <- dynamicTreeCut::cutreeDynamic(hier, distM=TOM, deepSplit=deepsplit,
                           minClusterSize=minclustsize,
                           pamRespectsDendro=FALSE, pamStage=FALSE)
    
    l2 <- levels(as.factor(clust))
    if (length(l2) > 1) {
      set.seed(555)
      merged <- try(WGCNA::mergeCloseModules(data_m, colors=clust,
                                      cutHeight=cutheight), silent=TRUE)
      mod_list <- if (inherits(merged, "try-error"))
        as.integer(clust)
      else
        as.integer(merged$colors)
    } else {
      mod_list <- as.integer(clust)
    }
    
    rm(ADJ, TOM, hier)
  }
  
  rm(data_m)
  save(mod_list, file=file.path(outloc, "mod_list.Rda"))
  
  dataA <- cbind(data_mzrt, dataA)
  dataA$Module_RTclust <- as.character(mod_list)
  dataA
}
