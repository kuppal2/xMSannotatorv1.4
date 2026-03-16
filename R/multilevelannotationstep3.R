multilevelannotationstep3 <-
function(
    outloc1,
    adduct_weights = NA,
    boostIDs = NA,
    pathwaycheckmode = "p",
    scorethresh = 0.1,
    require_primary_adduct = TRUE
) {
  
  suppressWarnings(library(data.table))
  setwd(outloc1)
  
  load("step1_results.Rda")
  load("chemCompMZ.Rda")
  
  if (!exists("adduct_weights") || isTRUE(is.na(adduct_weights)[1])) {
    adduct_weights <- data.frame(
      Adduct = c("M+H","M-H"),
      Weight = c(1,1)
    )
  }
  
  if (!exists("max_diff_rt")) max_diff_rt <- 10
  
  DT <- fread("Stage3A.csv")
  setDT(DT)
  
  hmdbbad <- c("HMDB29244","HMDB29245","HMDB29246")
  DT <- DT[!chemical_ID %in% hmdbbad]
  
  DT[, module_num := sub("_[0-9]*", "", Module_RTclust)]
  
  primary_adducts <- as.character(adduct_weights[,1])
  
  DTf <- DT[
    score >= scorethresh &
      Adduct %in% primary_adducts
  ]
  
  total_chems <- uniqueN(DTf$chemical_ID)
  pthresh <- 0.05
  
  ########################################################
  # --------- FUNCTION FOR PATHWAY ENRICHMENT ----------
  ########################################################
  
  run_pathway_enrichment <- function(mapDT) {
    
    DT_join <- DTf[mapDT, on="chemical_ID", nomatch=0]
    if (nrow(DT_join) == 0) return(NULL)
    
    # counts per pathway + module
    pm <- DT_join[
      , .(a = uniqueN(chemical_ID)),
      by = .(pathway, module_num)
    ]
    
    # pathway totals
    pt <- DT_join[
      , .(path_total = uniqueN(chemical_ID)),
      by = pathway
    ]
    
    # module totals
    mt <- DTf[
      , .(module_total = uniqueN(chemical_ID)),
      by = module_num
    ]
    
    pm[pt, path_total := i.path_total, on="pathway"]
    pm[mt, module_total := i.module_total, on="module_num"]
    
    # VECTORISED hypergeometric test
    pm[, p_module :=
         phyper(
           q = a - 1,
           m = path_total,
           n = total_chems - path_total,
           k = module_total,
           lower.tail = FALSE
         )
    ]
    
    sig_paths <- unique(pm[p_module <= pthresh]$pathway)
    
    if (length(sig_paths) == 0) return(NULL)
    
    DT_join[pathway %in% sig_paths,
            unique(chemical_ID)]
  }
  
  ########################################################
  # ---------------- KEGG -------------------------------
  ########################################################
  
  if (!is.na(pathwaycheckmode) &&
      exists("db_name") &&
      db_name == "KEGG") {
    
    data(keggotherinf)
    keggDT <- as.data.table(keggotherinf)
    
    chem_col <- names(keggDT)[1]
    path_col <- names(keggDT)[4]
    
    kegg_map <- keggDT[
      , .(pathway = unlist(strsplit(as.character(get(path_col)), ";"))),
      by = .(chemical_ID = get(chem_col))
    ]
    
    kegg_map <- kegg_map[
      pathway != "-" & pathway != "map01100"
    ]
    
    boost_chems <- run_pathway_enrichment(kegg_map)
    
    if (!is.null(boost_chems) && length(boost_chems) > 0) {
      
      n_ids <- length(unique(boost_chems))
      n_rows_before <- DT[chemical_ID %in% boost_chems, .N]
      
      DT[chemical_ID %in% boost_chems, score := score + 10]
      
      message(sprintf(
        "Pathway boost applied to %d unique chemical_IDs (%d rows updated).",
        n_ids,
        n_rows_before
      ))
      
    } else {
      message("No chemicals boosted during pathway evaluation.")
    }
  }
  
  ########################################################
  # ---------------- HMDB -------------------------------
  ########################################################
  
  if (!is.na(pathwaycheckmode) &&
      exists("db_name") &&
      db_name == "HMDB") {
    
    data(hmdbAllinf)
    hmdbDT <- as.data.table(hmdbAllinf)
    
    if (ncol(hmdbDT) >= 27)
      hmdbDT <- hmdbDT[, -c(26:27), with=FALSE]
    
    chem_col <- names(hmdbDT)[1]
    path_col <- names(hmdbDT)[14]
    
    hmdb_map <- hmdbDT[
      , .(pathway = unlist(strsplit(as.character(get(path_col)), ";"))),
      by = .(chemical_ID = get(chem_col))
    ]
    
    hmdb_map <- hmdb_map[pathway != "-"]
    
    boost_chems <- run_pathway_enrichment(hmdb_map)
    
    
    if (!is.null(boost_chems) && length(boost_chems) > 0) {
      
      n_ids <- length(unique(boost_chems))
      n_rows_before <- DT[chemical_ID %in% boost_chems, .N]
      
      DT[chemical_ID %in% boost_chems, score := score + 10]
      
      message(sprintf(
        "Pathway boost applied to %d unique chemical_IDs (%d rows updated).",
        n_ids,
        n_rows_before
      ))
      
    } else {
      message("No chemicals boosted during pathway evaluation.")
    }
  }
  
  ########################################################
  # ---------------- boostIDs ---------------------------
  ########################################################
  
  if (length(boostIDs) > 0 && !all(is.na(boostIDs))) {
    boost_idx <- DT$chemical_ID %in% boostIDs
    if (any(boost_idx))
      DT[boost_idx, score :=
           max(DT$score, na.rm=TRUE) * 100]
  }
  
  ########################################################
  # ---------------- OUTPUT -----------------------------
  ########################################################
  
  core_cols <- c(
    "chemical_ID","score","Module_RTclust",
    "mz","time","MatchCategory",
    "theoretical.mz","Name","Formula",
    "MonoisotopicMass","Adduct","ISgroup",
    "mean_int_vec","MD"
  )
  
  avail_cols <- intersect(core_cols, names(DT))
  DT_out <- DT[, ..avail_cols]
  
  fwrite(DT_out, "Stage3B.csv")
  saveRDS(DT_out, "Stage3.rds")
  
  gc()
  
  return(as.data.frame(DT_out))
}
