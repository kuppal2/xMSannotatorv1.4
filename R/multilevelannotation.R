multilevelannotation <-
function(
    dataA, max.mz.diff=10, max.rt.diff=10, cormethod="pearson",
    num_nodes=2, queryadductlist=c("all"),
    gradienttype="Acetonitrile", mode="pos", outloc=getwd(), db_name="HMDB",
    adduct_weights=NA, num_sets=3000, allsteps=TRUE,
    corthresh=0.7, NOPS_check=TRUE, customIDs=NA,
    missing.value=NA, deepsplit=2, networktype="unsigned",
    minclustsize=10, module.merge.dissimilarity=0.2,
    filter.by=c("M+H"), redundancy_check=TRUE,
    min_ions_perchem=1, biofluid.location=NA, origin=NA,
    status=NA, boostIDs=NA, max_isp=5,
    MplusH.abundance.ratio.check=FALSE, customDB=NA,
    HMDBselect="union", mass_defect_window=0.01,
    mass_defect_mode="pos", dbAllinf=NA,
    pathwaycheckmode="pm", check_isp_abundance=TRUE,
    include_secondary_associations=FALSE, time_step=2,
    cor.type="full",
    isotopes_losses_fragments_transformations=NA,
    clustmethod="auto",
    require_primary_adduct=FALSE,annotation_mode=c("onepass","twopass_moduleschemical_with_primary_hits"),
    primary_adducts=c("M+H", "M+Na", "M+K", "M+NH4", "M+H-H2O","M-H", "M+FA-H", "M+Cl", "M-H-H2O", "M+Na-2H"),peakID_name=NA)   # v2.2.1: gate Confidence>0 on primary adduct presence
{
  options(warn=-1)

  .xms_read_fp <- function(fp_file) {
    if (file.exists(fp_file)) readLines(fp_file, warn = FALSE)[1L] else ""
  }

  .xms_write_fp <- function(fp_file, fp) {
    writeLines(as.character(fp), fp_file)
  }

  .xms_fingerprint <- function(...) {

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



  parallel:::setDefaultClusterOptions(setup_strategy="sequential")
  WGCNA::allowWGCNAThreads(nThreads=num_nodes)

  annotation_mode=annotation_mode[1]

  dataA      <- as.data.frame(dataA)
  #print(clustmethod)




  dataA$mz   <- round(dataA$mz, 5L)
  dataA$time <- round(dataA$time, 1L)

  check_id_col <- (match(tolower(colnames(dataA)[1]), tolower(c("Peak_ID","peakID","peak","feature","feature_ID",peakID_name))))

  if(!is.na(check_id_col)){

    peakID_mz_time<-dataA[,c(1:3)]
    dataA<-dataA[,-c(1)]
  }else{

    peakID_mz_time<-NA
  }

  dataA[, -c(1L:2L)] <- round(dataA[, -c(1L:2L)], 1L)

  s1 <- apply(dataA[, -c(1L:2L)], 1L, sd)

  bad_feats_allzeros <- which(s1 == 0)
  if (length(bad_feats_allzeros) > 0L) {
    print(paste("Removing", length(bad_feats_allzeros),
                "features with zero standard deviation."))
    dataA <- dataA[-bad_feats_allzeros, , drop=FALSE]
    print(dim(dataA))
  }

  if (is.na(module.merge.dissimilarity))
    module.merge.dissimilarity <- 1 - corthresh

  if (!is.na(isotopes_losses_fragments_transformations)[1L]) {
    rownames_isotopes_mat <- isotopes_losses_fragments_transformations$Adduct
    isotopes_losses_fragments_transformations <- as.numeric(as.character(
      isotopes_losses_fragments_transformations[, 2L]))
    isotopes_losses_fragments_transformations <- as.data.frame(
      isotopes_losses_fragments_transformations)
    rownames(isotopes_losses_fragments_transformations) <- rownames_isotopes_mat
  }

  if (!is.na(customIDs)[1L])
    customIDs <- as.data.frame(customIDs)

  data(adduct_table)

  print(paste("Annotating using", db_name, "database:"))

  max_diff_rt <- max.rt.diff
  cutheight   <- 1 - corthresh

  # --- normalise adduct query ---
  if (isTRUE(queryadductlist == "pos")) queryadductlist <- "positive"
  if (isTRUE(queryadductlist == "neg")) queryadductlist <- "negative"

  adduct_table$Adduct <- gsub(adduct_table$Adduct, pattern="2M\\+2H",  replacement="M+H")
  adduct_table$Adduct <- gsub(adduct_table$Adduct, pattern="3M\\+3H",  replacement="M+H")
  adduct_table$Adduct <- gsub(adduct_table$Adduct, pattern="2M\\+2Na", replacement="M+Na")
  adduct_table$Adduct <- gsub(adduct_table$Adduct, pattern="3M\\+3Na", replacement="M+Na")
  adduct_table$Adduct <- gsub(adduct_table$Adduct, pattern="2M\\+2K",  replacement="M+K")
  adduct_table$Adduct <- gsub(adduct_table$Adduct, pattern="3M\\+3K",  replacement="M+K")
  adduct_table <- unique(adduct_table)


  if (isTRUE((queryadductlist[1L] == "all" & mode == "pos")[1L])) {
    adduct_names <- adduct_table$Adduct[
      (adduct_table$Type == "S"           & adduct_table$Mode == "positive") |
        (adduct_table$Type == gradienttype  & adduct_table$Mode == "positive")]
    adduct_table <- adduct_table[adduct_table$Adduct %in% adduct_names, ]
  } else if (isTRUE((queryadductlist[1L] == "all" & mode == "neg")[1L])) {
    adduct_names <- adduct_table$Adduct[
      (adduct_table$Type == "S"           & adduct_table$Mode == "negative") |
        (adduct_table$Type == gradienttype  & adduct_table$Mode == "negative")]
    adduct_table <- adduct_table[adduct_table$Adduct %in% adduct_names, ]
  } else if (isTRUE((queryadductlist[1L] == "positive" & mode == "pos")[1L])) {
    adduct_names <- adduct_table$Adduct[
      (adduct_table$Type == "S"           & adduct_table$Mode == "positive") |
        (adduct_table$Type == gradienttype  & adduct_table$Mode == "positive")]
    adduct_table <- adduct_table[adduct_table$Adduct %in% adduct_names, ]
  } else if (isTRUE((queryadductlist[1L] == "negative" & mode == "neg")[1L])) {
    adduct_names <- adduct_table$Adduct[
      (adduct_table$Type == "S"           & adduct_table$Mode == "negative") |
        (adduct_table$Type == gradienttype  & adduct_table$Mode == "negative")]
    adduct_table <- adduct_table[adduct_table$Adduct %in% adduct_names, ]
  } else {
    # User passed a specific adduct list
    adduct_names <- adduct_table$Adduct[adduct_table$Adduct %in% queryadductlist]
    adduct_table <- adduct_table[adduct_table$Adduct %in% adduct_names, ]
    if (isTRUE((length(adduct_names) < 1L)[1L]))
      stop(paste("No valid adducts found in adduct_table for queryadductlist:",
                 paste(queryadductlist, collapse=", ")))
  }
  adduct_names <- unique(adduct_names)

  suppressWarnings(dir.create(outloc))
  #setwd(outloc)

  outloc_allres <- outloc

  l1            <- list.files(outloc_allres)
  check_step1   <- which(l1 == "step1_results.Rda")

  if (length(check_step1) < 1L) {


    check_levelA <- which(l1 == "xMSannotator_levelA_modules.Rda")


    if (!is.na(missing.value)) {
      for (j in seq_along(dataA)) {
        data.table::set(dataA,
            i = which(dataA[[j]] == missing.value),
            j = j,
            value = NA)
      }
    }
    data.table::setDT(dataA)
    dataA <- unique(dataA, by = c("mz","time"))
    mzid <- sprintf("%.5f_%.2f", dataA$mz, dataA$time)

    if (length(which(duplicated(mzid))) > 0L)
      dataA <- dataA[-which(duplicated(mzid)), , drop=FALSE]

    #mzid <- paste(dataA$mz, dataA$time, sep="_")
    dataA <- unique(dataA)
    mean_int_vec <- rowMeans(dataA[, -c(1L:2L)], na.rm=TRUE)
    print(dim(dataA))

    if (length(check_levelA) < 1L) {

      # ----------------------------------------------------------------
      # CLUSTERING + RT METHOD DISPATCH
      #
      # clustmethod accepts four compound strings:
      #   "Graph-RTLinear"  (default) Louvain graph + gap-based RT split
      #   "Graph-RTKDE"               Louvain graph + KDE turnpoint RT
      #   "WGCNA-RTLinear"            blockwiseModules + gap-based RT split
      #   "WGCNA-RTKDE"               blockwiseModules + KDE turnpoint RT
      #
      # Legacy values "WGCNA", "WGCNA1", "graph" are accepted and mapped
      # to their nearest equivalent for backward compatibility.
      # ----------------------------------------------------------------
      clustmethod_norm <- switch(tolower(clustmethod),
                                 "wgcna"        = "WGCNA-RTKDE",
                                 "wgcna1"       = "WGCNA-RTKDE",
                                 "graph"        = "Graph-RTLinear",
                                 "ultra"   = "Ultra-RTKDE",
                                 clustmethod
      )

      valid_methods <- c(
        "Graph-RTLinear", "Graph-RTKDE",
        "WGCNA-RTLinear", "WGCNA-RTKDE",
        "Ultra-RTKDE", "Ultra-RTLinear","Ultra","auto"
      )

      # normalize for comparison
      clustmethod_norm <- valid_methods[
        match(tolower(clustmethod_norm), tolower(valid_methods))
      ]

      if (is.na(clustmethod_norm)) {
        stop(
          "clustmethod must be one of: ",
          paste(valid_methods, collapse=", "),
          '. Got: "', clustmethod, '"'
        )
      }else{

        print(paste0("Status 1: Computing modules using ",clustmethod_norm," clustering"))
      }

      if(clustmethod=="auto"){

        if(nrow(dataA)>10000 & ncol(dataA)>30){
          use_ultra_clust=TRUE
          use_graph_clust=FALSE
          use_linear_rt   <- endsWith(tolower(clustmethod_norm),   "rtlinear")

        }else{
          clustmethod_norm<-"wgcna"
          use_wgcna_clust <- startsWith(tolower(clustmethod_norm), "wgcna")
          use_ultra_clust <- startsWith(tolower(clustmethod_norm), "ultra")
          use_graph_clust <- startsWith(tolower(clustmethod_norm), "graph")
          use_linear_rt   <- endsWith(tolower(clustmethod_norm),   "rtlinear")

        }
      }else{

        use_wgcna_clust <- startsWith(tolower(clustmethod_norm), "wgcna")
        use_ultra_clust <- startsWith(tolower(clustmethod_norm), "ultra")
        use_graph_clust <- startsWith(tolower(clustmethod_norm), "graph")
        use_linear_rt   <- endsWith(tolower(clustmethod_norm),   "rtlinear")

      }
      print(paste0("Status 1: Module clustering = \"", clustmethod_norm, "\""))


      if (use_graph_clust) {
        # ------------------------------------------------------------------
        # COMPUTE GLOBAL CORRELATION (needed by scorer in Steps 3+)
        # ------------------------------------------------------------------
        system.time(global_cor <- WGCNA::cor(
          t(dataA[, -c(1L:2L)]), nThreads=num_nodes,
          method=cormethod, use='p'))

        global_cor[(global_cor) < corthresh] <- 0
        global_cor <- drop0(Matrix(global_cor, sparse=TRUE))

        if (cor.type == "partial")
          global_cor <- cor2pcor(global_cor)

        global_cor <- round(global_cor, 2L)
        cor_na_index <- which(is.na(global_cor))
        if (length(cor_na_index) > 0L) global_cor[cor_na_index] <- 0



        setwd(outloc_allres)

        mzid <- paste(dataA$mz, dataA$time, sep="_")
        colnames(global_cor) <- mzid
        rownames(global_cor) <- mzid
        save(global_cor, file="global_cor.Rda")

        a1 <- rowSums(global_cor > corthresh)
        write.csv(a1, file=paste0("NumberOfCorrelationsPerFeature_cor", corthresh, ".csv"))

        levelA_res <- .xms_fast_cluster(
          dataA        = dataA,
          global_cor   = global_cor,
          mzid         = mzid,
          corthresh    = corthresh,
          max_diff_rt  = max_diff_rt,
          minclustsize = minclustsize,
          num_nodes    = num_nodes,
          outloc       = outloc)

      } else {

        if (use_ultra_clust) {

          # print("Running ultra one-shot RT-intensity clustering")
          levelA_res <- suppressWarnings(run_cpp_metabolomics_engine(dataA,alpha=corthresh))
          gc()
          ultra_mods <- sort(unique(levelA_res$Module_RTclust))
          # dataA<-dataA[,-which(colnames(dataA)=="Module_RTclust")]
          rt_grouped <- lapply(ultra_mods, function(gnum) {
            sub <- levelA_res[levelA_res$Module_RTclust == gnum,
                              c("mz","time"), drop=FALSE]
            if (nrow(sub) == 0L) return(NULL)
            if (use_linear_rt) {
              tryCatch(
                group_by_rt_linear(sub, time_step=time_step,
                                   max_diff_rt=max_diff_rt, groupnum=gnum),
                error=function(e) { sub$Module_RTclust <- paste0(gnum,"_1"); sub })
            } else {
              tryCatch(
                group_by_rt_histv1(sub, time_step=time_step,
                                   max_diff_rt=max_diff_rt, groupnum=gnum),
                error=function(e) { sub$Module_RTclust <- paste0(gnum,"_1"); sub })
            }
          })
          levelA_res <- as.data.frame(
            data.table::rbindlist(Filter(Negate(is.null), rt_grouped),
                                  use.names=TRUE, fill=TRUE))

        }else{
Annotation
          if (use_wgcna_clust) {
            # ------------------------------------------------------------------
            # COMPUTE GLOBAL CORRELATION (needed by scorer in Steps 3+)
            # ------------------------------------------------------------------
            system.time(global_cor <- WGCNA::cor(
              t(dataA[, -c(1L:2L)]), nThreads=num_nodes,
              method=cormethod, use='p'))

            # library(Matrix)
            #global_cor <- drop0(Matrix(global_cor * (abs(global_cor) > 0.6), sparse=TRUE))

            if (cor.type == "partial")
              global_cor <- cor2pcor(global_cor)

            global_cor <- round(global_cor, 2L)
            cor_na_index <- which(is.na(global_cor))
            if (length(cor_na_index) > 0L) global_cor[cor_na_index] <- 0


            dataA <- unique(dataA)
            setwd(outloc_allres)

            mzid <- paste(dataA$mz, dataA$time, sep="_")
            colnames(global_cor) <- mzid
            rownames(global_cor) <- mzid
            save(global_cor, file="global_cor.Rda")

            a1 <- rowSums(global_cor > corthresh)
            write.csv(a1, file=paste0("NumberOfCorrelationsPerFeature_cor", corthresh, ".csv"))

            mycl_metabs   <- NA
            set.seed(555)

            diss <- as.dist(1 - global_cor)
            rm(global_cor)
            gc()
            set.seed(555)
            hr <- fastcluster::hclust(diss, method="complete")
            rm(diss)
            gc()

            set.seed(555)
            mycl_metabs <- dynamicTreeCut::cutreeDynamic(
              hr,
              deepSplit=deepsplit,
              minClusterSize=round(0.1 * length(mzid)),
              pamRespectsDendro=FALSE,
              pamStage=FALSE,
              verbose=0
            )

            set.seed(555)
            levelA_res <- get_peak_blocks_modulesvhclust(
              dataA=dataA, simmat=global_cor, mycl_metabs=mycl_metabs,
              adjacencyfromsimilarity=FALSE, outloc=outloc,
              column.rm.index=NA, deepsplit=deepsplit,
              minclustsize=minclustsize, cutheight=cutheight,
              cormethod=cormethod, networktype=networktype,
              num_nodes=num_nodes, step1log2scale=FALSE)
            # RT sub-grouping applied here for the WGCNA path
            wgcna_mods <- sort(unique(levelA_res$Module_RTclust))
            rt_grouped <- lapply(wgcna_mods, function(gnum) {
              sub <- levelA_res[levelA_res$Module_RTclust == gnum,
                                c("mz","time"), drop=FALSE]
              if (nrow(sub) == 0L) return(NULL)
              if (use_linear_rt) {
                tryCatch(
                  group_by_rt_linear(sub, time_step=time_step,
                                     max_diff_rt=max_diff_rt, groupnum=gnum),
                  error=function(e) { sub$Module_RTclust <- paste0(gnum,"_1"); sub })
              } else {
                tryCatch(
                  group_by_rt_histv1(sub, time_step=time_step,
                                     max_diff_rt=max_diff_rt, groupnum=gnum),
                  error=function(e) { sub$Module_RTclust <- paste0(gnum,"_1"); sub })
              }
            })
            levelA_res <- as.data.frame(
              data.table::rbindlist(Filter(Negate(is.null), rt_grouped),
                                    use.names=TRUE, fill=TRUE))
          }
        }
      }
      #setwd(outloc_allres)
      save(levelA_res, file=file.path(outloc,"xMSannotator_levelA_modules.Rda"))
      rm(global_cor)



      write.csv(levelA_res, file=file.path(outloc,"Stage1.csv"), row.names=FALSE)


    } else {
      #setwd(outloc_allres)
      #load("xMSannotator_levelA_modules.Rda")
      load(file.path(outloc_allres, "xMSannotator_levelA_modules.Rda"))
    }

    #setwd(outloc)

    levelA_res <- levelA_res[order(levelA_res$mz, levelA_res$time), ]
    levelA_res <- levelA_res[, 1L:3L]
    #dataA      <- dataA[order(dataA$mz, dataA$time), ]

    levelA_res <- cbind(levelA_res[, 1L:3L], mean_int_vec)
    #dataA      <- dataA[, 1L:2L]
    #rm(dataA)
    #############################################################################
    # adduct_names and adduct_table already resolved correctly at function entry

    # ------------------------------------------------------------------
    # LOAD DATABASE
    # ------------------------------------------------------------------
    chemCompMZ <- load_db(db_name, customDB, adduct_names, adduct_table,
                           customIDs, biofluid.location, origin, status,
                           HMDBselect)

    data(adduct_table)

    if (is.na(adduct_weights)[1L]) {
      adduct_weights <- data.frame(
        Adduct=c("M+H","M-H"), Weight=c(1,1), stringsAsFactors=FALSE)
    }

    chemCompMZ$Name    <- gsub("([-:();])|[[:punct:]]","\\1", chemCompMZ$Name)
    chemCompMZ$Formula <- gsub("([-:();])|[[:punct:]]","\\1", chemCompMZ$Formula)

    cnames    <- colnames(chemCompMZ)
    cnames[2L] <- "chemical_ID"
    colnames(chemCompMZ) <- cnames

    chemCompMZ <- as.data.frame(chemCompMZ)

    #setwd(outloc)
    l1 <- list.files(outloc)
    check_levelB <- which(l1 == "xMSannotator_levelB.Rda")

    #save(chemCompMZ, file="chemCompMZ.Rda")

    save(chemCompMZ, file = file.path(outloc, "chemCompMZ.Rda"))

    primary_adducts   <- intersect(primary_adducts, adduct_names)
    secondary_adducts <- setdiff(adduct_names, primary_adducts)

    if (length(check_levelB) < 1L) {

      # ----------------------------------------------------------------
      # HELPER: run one annotation pass and return filtered levelB_res
      # ----------------------------------------------------------------
      .run_annotation_pass <- function(dataA_pass,
                                       chemCompMZ_pass,
                                       pass_adducts,
                                       max.mz.diff) {

        chemCompMZ_pass <- chemCompMZ_pass[
          chemCompMZ_pass$Adduct %in% pass_adducts, ]

        chemCompMZ_pass <- chemCompMZ_pass[
          order(chemCompMZ_pass$mz), ]

        res <- run_ultra_annotation_engine(
          data_mz        = as.numeric(dataA_pass$mz),
          data_time      = as.numeric(dataA_pass$time),
          db_mz_R          = as.numeric(chemCompMZ_pass$mz),
          db_id_R = chemCompMZ_pass$chemical_ID,
          db_name_R        = chemCompMZ_pass$Name,
          db_formula_R     = chemCompMZ_pass$Formula,
          db_mono_mass_R   = as.numeric(chemCompMZ_pass$MonoisotopicMass),
          db_adduct_R      = chemCompMZ_pass$Adduct,
          db_adduct_mass_R = as.numeric(chemCompMZ_pass$AdductMass),
          query_adducts = pass_adducts,
          ppm_tol        = max.mz.diff
        )

        if (nrow(res) == 0) return(NULL)

        res$MatchCategory <- "Multiple"

        return(res)
      }
      # ----------------------------------------------------------------
      # PASS 2A: PRIMARY ADDUCTS
      # ----------------------------------------------------------------
      print("=== Step 2A: Primary adduct annotation ===")
      chemCompMZ$mz <- suppressWarnings(as.numeric(chemCompMZ$mz))

      chemCompMZ <- chemCompMZ[is.finite(chemCompMZ$mz), ]

      chemCompMZ <- chemCompMZ[order(chemCompMZ$mz), ]

      levelB_primary <- .run_annotation_pass(
        dataA_pass      = levelA_res,
        chemCompMZ_pass = chemCompMZ,
        pass_adducts    = primary_adducts,
        max.mz.diff     = max.mz.diff)

      data.table::setDT(levelB_primary)
      data.table::setDT(levelA_res)
      levelB_primary <- merge(
        levelB_primary,
        levelA_res[, .(mz, time, Module_RTclust)],
        by = c("mz", "time"),
        all.x = TRUE,
        allow.cartesian = FALSE
      )


      if (is.null(levelB_primary) || nrow(levelB_primary) == 0L)
        stop("Step 2A produced no primary adduct matches. Check adduct list and mass tolerance.")

     # save(levelB_primary, file="xMSannotator_levelB_primary.Rda")

      save(levelB_primary, file=file.path(outloc, "xMSannotator_levelB_primary.Rda"))

      print(paste("Step 2A complete:", nrow(levelB_primary), "primary matches"))


      levelB_secondary <- NULL


      if (length(secondary_adducts) > 0L &&
          !is.null(levelB_primary) &&
          nrow(levelB_primary) > 0) {

        print("=== Step 2B: Secondary adduct annotation ===")

        primary_ids     <- unique(levelB_primary$chemical_ID)
        primary_modules <- unique(levelB_primary$Module_RTclust)

        dataA_secondary      <- levelA_res
        chemCompMZ_secondary <- chemCompMZ

        if (annotation_mode == "twopass_modules_with_primary_hits") {

          # OPTION 1:
          # Restrict features to modules that already have a primary hit.
          dataA_secondary <- levelA_res[
            dataA$Module_RTclust %in% primary_modules, ]

          chemCompMZ_secondary <- chemCompMZ[
            chemCompMZ$Adduct %in% secondary_adducts, ]

          print("Secondary mode: modules_with_primary_hits")

        } else if (annotation_mode == "twopass_compounds_with_primary_hits") {

          # OPTION 2:
          # Restrict DB to compounds that had a primary hit.
          chemCompMZ_secondary <- chemCompMZ[
            chemCompMZ$chemical_ID %in% primary_ids &
              chemCompMZ$Adduct %in% secondary_adducts, ]



          print("Secondary mode: compounds_with_primary_hits")

        } else if (annotation_mode == "twopass_moduleschemical_with_primary_hits") {

          # OPTION 3:
          # Restrict BOTH modules AND compounds.
          dataA_secondary <- levelA_res[
            levelA_res$Module_RTclust %in% primary_modules, ]

          chemCompMZ_secondary <- chemCompMZ[
            chemCompMZ$chemical_ID %in% primary_ids &
              chemCompMZ$Adduct %in% secondary_adducts, ]

          print("Secondary mode: moduleschemical_with_primary_hits")
        }


        #print(paste("dataA rows secondary:", nrow(dataA_secondary)))
        #print(paste("DB rows secondary:", nrow(chemCompMZ_secondary)))
        #print("Unique secondary adducts in DB:")
        #print(unique(chemCompMZ$Adduct))
        #print("Requested secondary adducts:")
        #print(secondary_adducts)

        if (nrow(dataA_secondary) > 0 &&
            nrow(chemCompMZ_secondary) > 0) {

          levelB_secondary <- .run_annotation_pass(
            dataA_pass      = dataA_secondary,
            chemCompMZ_pass = chemCompMZ_secondary,
            pass_adducts    = secondary_adducts,
            max.mz.diff     = max.mz.diff)

        }
        else
          {
          print("Secondary pass skipped. No eligible features or compounds")
        }
      }


      # ----------------------------------------------------------------
      # MERGE 2A + 2B into unified levelB_res
      # ----------------------------------------------------------------
      if (!is.null(levelB_secondary)) {
        #levelB_res <- unique(rbind(levelB_primary, levelB_secondary))
        data.table::setDT(levelB_secondary)
        data.table::setDT(levelA_res)
        levelB_secondary <- merge(
          levelB_secondary,
          levelA_res[, .(mz, time, Module_RTclust)],
          by = c("mz", "time"),
          all.x = TRUE,
          allow.cartesian = FALSE
        )
        # print(colnames(levelB_primary))
        # print(colnames(levelB_secondary))
        all_cols <- union(names(levelB_primary), names(levelB_secondary))

        data.table::setcolorder(levelB_primary, all_cols[all_cols %in% names(levelB_primary)])
        data.table::setcolorder(levelB_secondary, all_cols[all_cols %in% names(levelB_secondary)])

        levelB_res <- unique(
          data.table::rbindlist(
            list(levelB_primary, levelB_secondary),
            use.names = TRUE,
            fill = TRUE
          )
        )
      } else {
        levelB_res <- levelB_primary
      }
      rm(levelB_primary, levelB_secondary)

      print(paste("Combined DB matches (2A+2B):", nrow(levelB_res)))
      #save(levelB_res, file="xMSannotator_levelB.Rda")
      save(levelB_res, file=file.path(outloc, "xMSannotator_levelB.Rda"))

    } else {
      print("Status 2: Using existing m/z mapping results:")
      #load("xMSannotator_levelB.Rda")
       load(file.path(outloc, "xMSannotator_levelB.Rda"))
    }

    try(rm(hmdbAllinf, envir=.GlobalEnv), silent=TRUE)
    try(rm(hmdbAllinfv3.6), silent=TRUE)
    try(rm(hmdbCompMZ), silent=TRUE)

    levelA_res$mz   <- round(levelA_res$mz, 5L)
    levelB_res$mz   <- round(levelB_res$mz, 5L)

    levels_A <- names(table(levelA_res$Module_RTclust))
    save(levelA_res, file=file.path(outloc, "levelA_res.Rda"))
    save(levelB_res, file=file.path(outloc, "levelB_res.Rda"))

    levelA_res <- levelA_res[order(levelA_res$mz, levelA_res$time), ]

    dup_mzid <- which(duplicated(mzid))
    if (length(dup_mzid) > 0L) {
      levelA_res <- levelA_res[-dup_mzid, , drop=FALSE]
      dataA      <- dataA[-dup_mzid, , drop=FALSE]
    }

    levelA_res   <- levelA_res[order(levelA_res$Module_RTclust), ]
    module_num   <- gsub(levelA_res$Module_RTclust, pattern="_[0-9]{1,}", replacement="")

    levelA_res_all               <- levelA_res[, 1L:2L]
    levelA_res_all$Module_RTclust <- module_num

    levelA_res  <- levelA_res[order(levelA_res$mz, levelA_res$time), ]
    mzdefect    <- levelA_res$mz - round(levelA_res$mz)
    MD          <- mzdefect
    levelA_res1 <- cbind(levelA_res, MD)

    ISgroup      <- rep(1L, nrow(levelA_res1))
    #isop_res_md  <- cbind(levelA_res1[, c(2L:3L)], ISgroup,
    #                     levelA_res1[, c(1L, 4L, 5L)])
    #colnames(isop_res_md) <- c("mz","time","ISgroup","Module_RTclust","AvgIntensity","MD")
    levelA_res1 <- cbind(levelA_res1, ISgroup)
    rm(MD)

    multiresmat <- merge(levelB_res, levelA_res1, by=c("mz","time","Module_RTclust"))
    #print(head(multiresmat))
    multiresmat <- multiresmat[, c(
      "mz","time","MatchCategory","theoretical.mz","chemical_ID",
      "Name","Formula","MonoisotopicMass","Adduct","ISgroup",
      "Module_RTclust","mean_int_vec","MD")]
    # colnames(multiresmat) <- c(
    #  "mz","time","MatchCategory","theoretical.mz","chemical_ID",
    #  "Name","Formula","MonoisotopicMass","Adduct","ISgroup",
    #  "Module_RTclust","time.y","mean_int_vec","MD")

    rm(levelB_res)
    multiresmat <- multiresmat[order(multiresmat$Module_RTclust), ]

    t2        <- table(multiresmat$mz)
    uniquemz  <- names(which(t2 == 1L))
    multiresmat$MatchCategory <- rep("Multiple", nrow(multiresmat))
    multiresmat$MatchCategory[multiresmat$mz %in% uniquemz] <- "Unique"

    #setwd(outloc)
    multiresmat <- multiresmat[order(multiresmat$Module_RTclust,
                                     multiresmat$chemical_ID), ]
    dupmz <- multiresmat$mz[duplicated(multiresmat$mz)]
    multiresmat$MatchCategory <- rep("Multiple", nrow(multiresmat))
    multiresmat$MatchCategory[!multiresmat$mz %in% dupmz] <- "Unique"

    if (length(which(multiresmat$chemical_ID == "-")) > 0L)
      multiresmat <- multiresmat[-which(multiresmat$chemical_ID == "-"), ]

    cnames    <- colnames(dataA)
    cnames[2L] <- "time"
    colnames(dataA) <- cnames
    dataA      <- unique(dataA)

    cnames <- colnames(multiresmat)
    cnames[10L] <- "ISgroup"
    colnames(multiresmat) <- cnames

    #

    chem_ids <- names(which(table(multiresmat$chemical_ID) >= 1L))
    chemids  <- chem_ids

    cnames <- colnames(multiresmat)
    cnames[2L] <- "time"
    colnames(multiresmat) <- cnames

    rm(levelA_res); try(rm(levelB_res), silent=TRUE)
    mchemdata <- as.data.frame(multiresmat)
    rm(multiresmat)

    mchemdata$mz           <- as.numeric(as.character(mchemdata$mz))
    mchemdata$time         <- as.numeric(as.character(mchemdata$time))
    mchemdata$mz           <- round(mchemdata$mz, 5L)
    mchemdata$time         <- round(mchemdata$time, 1L)
    mchemdata$MD           <- round(as.numeric(as.character(mchemdata$MD)), 3L)
    mchemdata$mean_int_vec <- round(as.numeric(as.character(mchemdata$mean_int_vec)), 1L)

    level_module_isop_annot <- levelA_res1
    level_module_isop_annot$MD           <- round(as.numeric(as.character(level_module_isop_annot$MD)), 3L)
    level_module_isop_annot$mean_int_vec <- round(as.numeric(as.character(level_module_isop_annot$mean_int_vec)), 1L)
    #mchemdata$time.y       <- round(as.numeric(as.character(mchemdata$time.y)), 1L)

    bad_rows <- which(abs(mchemdata$time - mchemdata$time.y) > 0L)
    if (length(bad_rows) > 0L) mchemdata <- mchemdata[-bad_rows, ]

    chemids <- unique(mchemdata$chemical_ID)

    if (max_diff_rt >= 9999L) {
      module_num_mc <- gsub(mchemdata$Module_RTclust, pattern="_[0-9]{1,}", replacement="")
      mchemdata$Module_RTclust <- paste0(module_num_mc, "_0")
    }



    stage2_dir <- file.path(outloc, "stage2")
    dir.create(stage2_dir, showWarnings=FALSE)

    mzid_final <- paste(mchemdata$mz, mchemdata$time, sep="_")



    # --- partition chemids for parallel step 2 ---
    if (length(chemids) > 10L) {
      num_sets <- length(chemids) / 2L
    } else {
      num_sets <- 1L
    }

    list_winsize <- num_sets
    list_size    <- round(length(chemids) / list_winsize)

    if (length(chemids) > list_winsize) {
      g             <- seq(1L, length(chemids), list_size)
      chemids_split <- split(seq_along(chemids), f=factor(g))
      split_size    <- seq_len(list_winsize)
    } else {
      chemids_split <- split(seq_along(chemids), f=length(chemids))
      split_size    <- 1L
    }
    num_sets <- length(chemids_split)


    save(
      list=c("outloc","adduct_weights","boostIDs","pathwaycheckmode",
             "adduct_table","max_diff_rt","corthresh","filter.by",
             "max.rt.diff","max_isp","min_ions_perchem","max.mz.diff",
             "db_name","allsteps","redundancy_check","num_sets",
             "check_isp_abundance","include_secondary_associations",
             "isotopes_losses_fragments_transformations",
             "require_primary_adduct","check_id_col","peakID_mz_time"),
      file=file.path(outloc,"tempobjects.Rda"))

    #"isop_res_md",
    save(
      list=c("mchemdata","chemids","adduct_table","mzid","max_diff_rt",
             "corthresh","level_module_isop_annot",
             "chemids_split","max.mz.diff","outloc","num_sets","db_name",
             "num_nodes","adduct_weights","filter.by","max.rt.diff","max_isp",
             "MplusH.abundance.ratio.check","mass_defect_window",
             "mass_defect_mode","allsteps","check_isp_abundance",
             "include_secondary_associations",
             "isotopes_losses_fragments_transformations","primary_adducts"),
      file=file.path(outloc,"step1_results.Rda"))



    write.csv(mchemdata, file=file.path(outloc,"Stage2.csv"), row.names=FALSE)






    rm(mchemdata); rm(chemids); rm(mzid); try(rm(global_cor), silent=TRUE)
    rm(isop_res_md); rm(level_module_isop_annot); rm(dataA)
    try(rm(levelA_res1), silent=TRUE)
    try(rm(tall_unimzperchem), silent=TRUE)
    try(rm(tablemz), silent=TRUE)
    try(rm(levelB_res2), silent=TRUE)

    #rm(list=ls())
    load(file.path(outloc,"tempobjects.Rda"))

  } else {

    print("Status 1: Skipping step 1.")
    print("Status 2: Using existing step1_results.Rda file.")
    allsteps_temp <- allsteps
    load(file.path(outloc,"tempobjects.Rda"))
    load(file.path(outloc,"step1_results.Rda"))
    allsteps <- allsteps_temp
  }

  if (isTRUE(allsteps)) {

    # ==================================================================
    # STATUS 3 + 4: unified pre-check (v2.2.2.6)

    # mchemdata<-fread("Stage2.csv")


    .step2_sentinel <- file.path(outloc, "stage2", "step2_complete.txt")
    .stage3_rds     <- file.path(outloc, "Stage3.rds")
    .stage3_csv     <- file.path(outloc, "Stage3.csv")
    .fp3_file       <- file.path(outloc, "stage3_params.fp")
    .fp3_cached     <- .xms_read_fp(.fp3_file)

    # O(1) when sentinel exists; O(num_sets) only on the very first run
    # (or after a cache invalidation) when sentinel has not been written.
    .score_mtime <- if (file.exists(.step2_sentinel)) {
      readLines(.step2_sentinel, warn = FALSE)[1L]
    } else {
      .sf <- list.files(file.path(outloc, "stage2"),
                        pattern = "^chem_score.*\\.Rda$", full.names = TRUE)
      if (length(.sf) > 0L) as.character(max(file.mtime(.sf))) else "none"
    }

    .fp3_current  <- .xms_fingerprint(adduct_weights, boostIDs,
                                      pathwaycheckmode, .score_mtime)
    .stage34_skip <- (.fp3_current == .fp3_cached) &&
      (file.exists(.stage3_rds) || file.exists(.stage3_csv))

    if (.stage34_skip) {


      print("Status 3: Calculating scores. SKIPPED (Stage3 cache valid)")
      print("Status 4: Pathway evaluation.  SKIPPED (Stage3 cache valid)")


    } else {

      # ---- Slow path: run step2 workers then step3 ----
      print("Status 3: Calculating scores for individual chemicals/metabolites")
      load(file.path(outloc,"step1_results.Rda"))


      prep <- prepare_stage3_hybrid(
        mchemdata,
        adduct_table,
        adduct_weights,
        isotopes_losses_fragments_transformations,
        level_module_isop_annot
      )

      DT  <- prep$DT
      ISP <- prep$ISP
      trans_vec <- prep$trans_vec


      # ---------------------------
      # Sort main table
      # ---------------------------
      DT <- DT[order(DT$chemical_ID, DT$Module_RTclust), ]

      group_key <- paste(DT$chemical_ID, DT$Module_RTclust)

      group_factor <- as.integer(factor(group_key))

      group_start <- tapply(
        seq_along(group_factor) - 1,
        group_factor,
        min
      )

      group_end <- tapply(
        seq_along(group_factor) - 1,
        group_factor,
        max
      ) + 1

      group_start_cpp <- as.integer(group_start)
      group_end_cpp   <- as.integer(group_end)

      # ---------------------------
      # Sort isotope table
      # ---------------------------
      ISP <- ISP[order(ISP$Module_RTclust), ]


      expanded_new <- suppressWarnings(expand_isotopes_ultra_fast_cpp(
        DT$mz,
        DT$time,
        DT$theoretical.mz,
        DT$module_id,
        group_start_cpp,
        group_end_cpp,
        ISP$mz,
        ISP$time,
        ISP$module_id,
        DT$Formula,
        max_diff_rt,
        max_isp,
        max.mz.diff
      ))

      gc()
      parent_rows <- DT[expanded_new$parent_index + 1]
      parent_rows <- DT[expanded_new$parent_index + 1]

      parent_rows$mean_int_vec <- DT$mean_int_vec[expanded_new$parent_index + 1]
      parent_rows$mz            <- expanded_new$mz
      parent_rows$time          <- expanded_new$time

      parent_rows$Formula       <- expanded_new$Formula
      parent_rows$Adduct        <- expanded_new$Adduct

      expanded_full <- rbind(DT, parent_rows)

      rm(DT)
      rm(parent_rows)
      # ---------------------------
      # Re-sort + re-index
      # ---------------------------
      expanded_full <- expanded_full[
        order(expanded_full$chemical_ID,
              expanded_full$Module_RTclust),
      ]

      group_key <- paste(
        expanded_full$chemical_ID,
        expanded_full$Module_RTclust
      )

      group_factor <- as.integer(factor(group_key))

      idx <- seq_len(nrow(expanded_full)) - 1

      group_start_cpp <- tapply(idx, group_factor, min)
      group_end_cpp   <- tapply(idx, group_factor, max) + 1

      group_start_cpp <- as.integer(group_start_cpp)
      group_end_cpp   <- as.integer(group_end_cpp)

      stopifnot(
        max(group_end_cpp) <= nrow(expanded_full),
        min(group_start_cpp) >= 0,
        all(group_start_cpp < group_end_cpp)
      )

      #call the scoring function
      result <- suppressWarnings(stage3_score_engine_cpp(
        expanded_full,
        adduct_weights,
        primary_adducts,
        max_diff_rt,
        group_start_cpp,
        group_end_cpp
      ))
      gc()

      #print(summary(result$score))

      # keep only rows that survived scoring
      DT <- as.data.table(result)


      # summarize per chemical
      diag <- DT[, .(
        n_rows = .N,
        n_adducts = uniqueN(Adduct),
        adduct_list = paste(unique(Adduct), collapse=", "),
        rt_range = IQR(time),
        max_int = max(mean_int_vec),
        score = max(score,na.rm=TRUE)
      ), by = .(chemical_ID,Module_RTclust)]

      # sort by adduct diversity
      setorder(diag, -n_adducts, -score)


      write.csv(diag,file=file.path(outloc,"Stage3A_scoring_summary.csv"))

      write.csv(result,file=file.path(outloc,"Stage3A.csv"))






      # Write sentinel: future runs use this single file for O(1) mtime
      # lookup instead of stat-ing every chem_score*.Rda.
      writeLines(as.character(Sys.time()),
                 file.path(outloc, "stage2", "step2_complete.txt"))

      #setwd(outloc)

      #rm(list=ls())

     # try(rm(hmdbCompMZ), silent=TRUE)
      load(file.path(outloc,"tempobjects.Rda"))
      if (!exists("require_primary_adduct")) require_primary_adduct <- FALSE
      if (!exists("max.mz.diff"))           max.mz.diff            <- 10
      if (!exists("min_ions_perchem"))      min_ions_perchem        <- 1
      if (!exists("redundancy_check"))      redundancy_check        <- TRUE
      if (!exists("check_isp_abundance"))   check_isp_abundance     <- TRUE
      if (!exists("include_secondary_associations")) include_secondary_associations <- FALSE
      if (!exists("isotopes_losses_fragments_transformations")) isotopes_losses_fragments_transformations <- NA

      # Recompute fp3 after rm/reload (sentinel now guaranteed to exist)
      .step2_sentinel <- file.path(outloc, "stage2", "step2_complete.txt")
      .score_mtime    <- readLines(.step2_sentinel, warn = FALSE)[1L]
      .fp3_current    <- .xms_fingerprint(adduct_weights, boostIDs,
                                          pathwaycheckmode, .score_mtime)
      .fp3_file       <- file.path(outloc, "stage3_params.fp")

      print("Status 4: Pathway evaluation")
      annotres <- multilevelannotationstep3(
        outloc=outloc, adduct_weights=adduct_weights,
        boostIDs=boostIDs, pathwaycheckmode=pathwaycheckmode,
        require_primary_adduct=require_primary_adduct)
      .xms_write_fp(.fp3_file, .fp3_current)

    }
   # setwd(outloc)
    #rm(list=ls())
    try(rm(hmdbCompMZ), silent=TRUE)
    try(rm(hmdbCompMZ, env=.GlobalEnv), silent=TRUE)
    load(file.path(outloc,"tempobjects.Rda"))
    # Defensive defaults: older tempobjects.Rda files (pre-v2.2.1) may be missing
    # parameters added in later patches. Provide safe fallbacks so cached runs
    # don't fail with "object not found".
    if (!exists("require_primary_adduct")) require_primary_adduct <- FALSE
    if (!exists("max.mz.diff"))           max.mz.diff            <- 10
    if (!exists("min_ions_perchem"))      min_ions_perchem        <- 1
    if (!exists("redundancy_check"))      redundancy_check        <- TRUE
    if (!exists("check_isp_abundance"))   check_isp_abundance     <- TRUE
    if (!exists("include_secondary_associations")) include_secondary_associations <- FALSE
    if (!exists("isotopes_losses_fragments_transformations")) isotopes_losses_fragments_transformations <- NA

    ###Assign confidence levels###################
    stage3_mtime <- if (file.exists(file.path(outloc, "Stage3.rds")))
      as.character(file.mtime(file.path(outloc, "Stage3.rds")))
    else if (file.exists(file.path(outloc, "Stage3.csv")))
      as.character(file.mtime(file.path(outloc, "Stage3.csv")))
    else "none"

    fp4_current <- .xms_fingerprint(
      max_diff_rt, filter.by, adduct_weights, min_ions_perchem,
      boostIDs, max_isp, max.mz.diff,
      require_primary_adduct,
      stage3_mtime)
    fp4_file   <- file.path(outloc, "stage4_params.fp")
    fp4_cached <- .xms_read_fp(fp4_file)
    stage4_csv <- file.path(outloc, "Stage4.csv")

    if (fp4_current == fp4_cached && file.exists(stage4_csv)) {

      print("Status 5: Assigning confidence levels. SKIPPED (Stage4 cache valid)")
      annotresstage4 <- data.table::fread(stage4_csv, data.table=FALSE)

    } else {

      print("Status 5: Assigning confidence levels")
      stage3_csv <- file.path(outloc, "Stage3B.csv")
      annotresstage3<- data.table::fread(stage3_csv, data.table=FALSE)
      annotresstage4<-simple_xms_conf_fast(stage3_results=annotresstage3)

      write.csv(annotresstage4$cluster_summary,file=file.path(outloc,"Stage4_cluster_summary.csv"),row.names=FALSE)



      write.csv(annotresstage4$conf_mat%>%as.data.frame(),file=file.path(outloc,"Stage4.csv"),row.names=FALSE)




      .xms_write_fp(fp4_file, fp4_current)

    }

   # rm(list=ls())
    load(file.path(outloc,"tempobjects.Rda"))
    try(rm(hmdbAllinf, env=.GlobalEnv), silent=TRUE)

    if (!exists("require_primary_adduct")) require_primary_adduct <- FALSE
    if (!exists("max.mz.diff"))           max.mz.diff            <- 10
    if (!exists("min_ions_perchem"))      min_ions_perchem        <- 1
    if (!exists("redundancy_check"))      redundancy_check        <- TRUE
    if (!exists("check_isp_abundance"))   check_isp_abundance     <- TRUE
    if (!exists("include_secondary_associations")) include_secondary_associations <- FALSE
    if (!exists("isotopes_losses_fragments_transformations")) isotopes_losses_fragments_transformations <- NA

    if (isTRUE(redundancy_check)) {

      # ------------------------------------------------------------------
      # Status 6: Redundancy filtering (step5)
      # ------------------------------------------------------------------
      stage4_mtime <- if (file.exists(file.path(outloc, "Stage4.csv")))
        as.character(file.mtime(file.path(outloc, "Stage4.csv")))
      else "none"

      fp5_current <- .xms_fingerprint(
        max_diff_rt, filter.by, adduct_weights, min_ions_perchem,
        boostIDs, max_isp, max.mz.diff, db_name,
        require_primary_adduct,
        stage4_mtime)
      fp5_file   <- file.path(outloc, "stage5_params.fp")
      fp5_cached <- .xms_read_fp(fp5_file)
      stage5_csv <- file.path(outloc, "Stage5.csv")

      if (fp5_current == fp5_cached && file.exists(stage5_csv)) {

        print("Status 6: Redundancy filtering. SKIPPED (Stage5 cache valid)")
        #rm(list=ls())

        annotresstage5 <- data.table::fread(stage5_csv, data.table=FALSE)

      } else {

        print("Status 6: Redundancy filtering")
        annotresstage5 <- multilevelannotationstep5(
          outloc=outloc)


        # write_fp BEFORE rm. fp5_file/fp5_current are destroyed by rm(list=ls()).

        .xms_write_fp(fp5_file, fp5_current)
        #rm(list=ls())

      }

    } else {
      annotresstage5 <- NULL
    }

  } else {
    annotres <- NULL
  }


  if(!is.na(check_id_col)){
    print("Adding Peak ID column")

    stage1<-fread(file.path(outloc,"Stage1.csv"))

    annotresstage1<-merge(peakID_mz_time, stage1,by=c("mz","time"))
    write.csv(annotresstage1,file="Stage1.csv",row.names=FALSE)

    stage2<-fread(file.path(outloc,"Stage2.csv"))

    annotresstage2<-merge(peakID_mz_time, stage2,by=c("mz","time"))
    write.csv(annotresstage2,file="Stage2.csv",row.names=FALSE)

    stage3a<-fread(file.path(outloc,"Stage3A.csv"))

    annotresstage3A<-merge(peakID_mz_time, stage3a,by=c("mz","time"))
    write.csv(annotresstage3A,file="Stage3A.csv",row.names=FALSE)

    stage3b<-fread(file.path(outloc,"Stage3B.csv"))

    annotresstage3B<-merge(peakID_mz_time, stage3b,by=c("mz","time"))
    write.csv(annotresstage3B,file="Stage3B.csv",row.names=FALSE)

    stage4<-fread(file.path(outloc,"Stage4.csv"))

    annotresstage4<-merge(peakID_mz_time, stage4,by=c("mz","time"))
    write.csv(annotresstage4,file="Stage4_curated.csv",row.names=FALSE)

    stage5<-fread("Stage5.csv")

    annotresstage5<-merge(peakID_mz_time, stage5,by=c("mz","time"))
    write.csv(annotresstage5,file="Stage5.csv",row.names=FALSE)

  }

  print("################")
  print("Final: Processing complete")
  print("Output files description:")
  print("Stage 1 includes clustering of features based on intensity and retention time without annotations")
  print("Stage 2 includes clustering results along with simple m/z based database matching")

  if (isTRUE(allsteps)) {
    print("Stage 3 includes scores for annotations assigned in stage 2 based on multiple criteria")
    print("Stages 4 and 5 include the confidence levels before and after redundancy (multiple matches) filtering, respectively")
  }

  suppressWarnings(sink(file=NULL))
  #setwd(outloc)
  sink(file=file.path(outloc,"Readme.txt"))
  print("Output files description:")
  print("Stage1.csv: clustering without annotations")
  print("Stage2.csv: clustering + m/z database matching")
  if (isTRUE(allsteps)) {
    print("Stage3.csv: annotation scores")
    print("Stage4_allmatches.csv: confidence levels for all chemical_ID, Module_RTClust, mz, time groups")
    print("Stage4_curated.csv: confidence levels for best scoring chemical_ID, Module_RTClust, mz, time groups")
    print("Stage5.csv: confidence levels for best scoring chemical_ID, mz, time groups")
  }
  suppressWarnings(sink(file=NULL))
  return(list(results.stage4=annotresstage4,results.stage5=annotresstage5))
}
