load_db <- function(db_name, customDB, adduct_names, adduct_table,
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
