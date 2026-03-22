prepare_stage3_hybrid <-
function(
    mchemdata,
    adduct_table,
    adduct_weights,
    isotopes_losses_fragments_transformations = NA,
    level_module_isop_annot = NULL
) {
  
  library(data.table)
  
  DT <- as.data.table(unique(mchemdata))
  setorder(DT, chemical_ID, Module_RTclust)
  
  
  # --- Adduct mapping
  adduct_levels <- sort(unique(DT$Adduct))
  DT[, adduct_id := match(Adduct, adduct_levels) - 1L]
  
  adduct_weight_vec <- adduct_weights[
    match(adduct_levels, adduct_weights[,1]), 2
  ]
  adduct_weight_vec[is.na(adduct_weight_vec)] <- 0
  adduct_weight_vec <- as.numeric(adduct_weight_vec)
  
  # --- Charge mapping
  adduct_lookup <- as.data.table(adduct_table)[
    , .(charge = charge[1]), by = Adduct
  ]
  
  adduct_charge_vec <- adduct_lookup[
    match(adduct_levels, Adduct),
    charge
  ]
  
  adduct_charge_vec[is.na(adduct_charge_vec)] <- 1L
  adduct_charge_vec <- as.integer(adduct_charge_vec)
  
  DT[, charge := adduct_charge_vec[adduct_id + 1L]]
  
  # --- Group index
  DT[, group_id := .GRP, by=.(chemical_ID, Module_RTclust)]
  group_sizes <- DT[, .N, by=group_id]$N
  group_index <- cumsum(c(0L, group_sizes))
  group_index <- as.integer(group_index)
  
  # --- Transform vector
  if (!is.null(isotopes_losses_fragments_transformations) &&
      !is.na(isotopes_losses_fragments_transformations)[1]) {
    trans_vec <- as.numeric(
      isotopes_losses_fragments_transformations[,1]
    )
  } else {
    trans_vec <- numeric(0)
  }
  
  # --- ISP table (required)
  if (!is.null(level_module_isop_annot)) {
    ISP <- as.data.table(level_module_isop_annot)
    setorder(ISP, Module_RTclust)
    
    # ---- Convert Module_RTclust to integer IDs
    all_modules <- unique(
      c(DT$Module_RTclust, ISP$Module_RTclust)
    )
    
    module_map <- setNames(
      seq_along(all_modules),
      all_modules
    )
    
    DT$module_id  <- module_map[DT$Module_RTclust]
    ISP$module_id <- module_map[ISP$Module_RTclust]
    
  } else {
    ISP <- NULL
  }
  
  list(
    DT = DT,
    ISP = ISP,
    group_index = group_index,
    adduct_weight_vec = adduct_weight_vec,
    trans_vec = trans_vec
  )
}
