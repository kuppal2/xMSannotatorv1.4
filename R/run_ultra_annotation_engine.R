run_ultra_annotation_engine <-
function(
    data_mz,
    data_time,
    db_mz_R,
    db_id_R,
    db_name_R,
    db_formula_R,
    db_mono_mass_R,
    db_adduct_R,
    db_adduct_mass_R,
    query_adducts,
    ppm_tol = 10
) {
  
  n_features <- length(data_mz)
  n_db       <- length(db_mz_R)
  
  if (n_db == 0 || n_features == 0)
    return(data.frame())
  
  # Ensure DB sorted
  ord <- order(db_mz_R)
  db_mz_R        <- db_mz_R[ord]
  db_id_R        <- db_id_R[ord]
  db_name_R      <- db_name_R[ord]
  db_formula_R   <- db_formula_R[ord]
  db_mono_mass_R <- db_mono_mass_R[ord]
  db_adduct_R    <- db_adduct_R[ord]
  db_adduct_mass_R <- db_adduct_mass_R[ord]
  
  results <- vector("list", n_features)
  res_i <- 1
  
  for (i in seq_len(n_features)) {
    
    mz <- data_mz[i]
    
    ppm_window <- mz * ppm_tol / 1e6
    
    lower <- mz - ppm_window
    upper <- mz + ppm_window
    
    # Binary search region
    left  <- findInterval(lower, db_mz_R)
    right <- findInterval(upper, db_mz_R)
    
    if (right > left) {
      
      idx <- seq.int(left + 1, right)
      
      results[[res_i]] <- data.frame(
        mz               = rep(mz, length(idx)),
        time             = rep(data_time[i], length(idx)),
        theoretical.mz   = db_mz_R[idx],
        chemical_ID      = db_id_R[idx],
        Name             = db_name_R[idx],
        Formula          = db_formula_R[idx],
        MonoisotopicMass = db_mono_mass_R[idx],
        Adduct           = db_adduct_R[idx],
        AdductMass       = db_adduct_mass_R[idx],
        stringsAsFactors = FALSE
      )
      
      res_i <- res_i + 1
    }
  }
  
  if (res_i == 1)
    return(data.frame())
  
  res <- do.call(rbind, results[seq_len(res_i - 1)])
  
  return(res)
}
