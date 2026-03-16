simpleAnnotation <-
function(
    dataA,
    max.mz.diff = 10,
    num_nodes = 2,
    queryadductlist = c("M+2H","M+H+NH4","M+ACN+2H","M+2ACN+2H","M+H",
                        "M+NH4","M+Na","M+ACN+H","M+ACN+Na",
                        "M+2ACN+H","2M+H","2M+Na","2M+ACN+H"),
    gradienttype = "Acetonitrile",
    mode = "pos",
    outloc,
    db_name = "HMDB",
    customDB = NA,
    NOPS_check = TRUE
) {
  
  library(data.table)
  
  setDT(dataA)
  dataA <- dataA[, .(mz = as.numeric(mz), time = as.numeric(time))]
  
  data(adduct_table)
  setDT(adduct_table)
  
  # Normalize adducts (vectorized)
  adduct_table[, Adduct := gsub("2M\\+2H|3M\\+3H", "M+H", Adduct)]
  adduct_table[, Adduct := gsub("2M\\+2Na|3M\\+3Na", "M+Na", Adduct)]
  adduct_table[, Adduct := gsub("2M\\+2K|3M\\+3K", "M+K", Adduct)]
  adduct_table <- unique(adduct_table)
  
  # Filter adducts (vectorized)
  if ((queryadductlist[1] == "all" & mode == "pos") ||
      queryadductlist[1] == "positive") {
    
    adduct_names <- adduct_table[
      (Type == "S" & Mode == "positive") |
        (Type == gradienttype & Mode == "positive"),
      unique(Adduct)
    ]
    
  } else if ((queryadductlist[1] == "all" & mode == "neg") ||
             queryadductlist[1] == "negative") {
    
    adduct_names <- adduct_table[
      (Type == "S" & Mode == "negative") |
        (Type == gradienttype & Mode == "negative"),
      unique(Adduct)
    ]
    
  } else {
    
    adduct_names <- intersect(unique(adduct_table$Adduct), queryadductlist)
    
    if (length(adduct_names) == 0)
      stop("Invalid adducts selected.")
  }
  
  # Load database
  chemCompMZ <- switch(
    db_name,
    "KEGG" = { data(keggCompMZ); as.data.table(keggCompMZ) },
    "HMDB" = { data(hmdbCompMZ); as.data.table(hmdbCompMZ) },
    "T3DB" = { data(t3dbCompMZ); as.data.table(t3dbCompMZ) },
    "LipidMaps" = { data(lipidmapsCompMZ); as.data.table(lipidmapsCompMZ) },
    "Phytochelatin" = { data(phytochelatinCompMZ); as.data.table(phytochelatinCompMZ) },
    "Custom" = as.data.table(customDB),
    stop("Unsupported database")
  )
  
  chemCompMZ <- chemCompMZ[Adduct %in% adduct_names]
  
  setDT(dataA)
  setDT(chemCompMZ)
  
  dataA[, mz := as.numeric(mz)]
  chemCompMZ[, mz := as.numeric(mz)]
  
  # Compute ppm window
  dataA[, mz_tol := mz * max.mz.diff / 1e6]
  dataA[, mz_lower := mz - mz_tol]
  dataA[, mz_upper := mz + mz_tol]
  
  setkey(chemCompMZ, mz)
  
  levelB_res <- dataA[
    chemCompMZ,
    on = .(mz_lower <= mz, mz_upper >= mz),
    allow.cartesian = TRUE,
    nomatch = 0
  ]
  
  if (nrow(levelB_res) == 0) {
    message("No matches found.")
    return(-1)
  }
  
  # Golden rule filtering (vectorized)
  uniq_formula <- unique(levelB_res$Formula)
  golden_check <- rbindlist(
    lapply(uniq_formula, check_golden_rules, NOPS_check = NOPS_check)
  )
  valid_formula <- golden_check[bool_check == 1, curformula]
  
  levelB_res <- levelB_res[Formula %in% valid_formula]
  
  # Water filtering (vectorized)
  water_adducts <- c("M+H-H2O","M+H-2H2O","M-H2O-H")
  
  water_res <- levelB_res[Adduct %in% water_adducts]
  non_water_res <- levelB_res[!Adduct %in% water_adducts]
  
  if (nrow(water_res) > 0) {
    water_res[, numO := sapply(Formula, check_element, element="O")]
    water_res <- water_res[numO > 0]
  }
  
  levelB_res <- rbind(non_water_res, water_res, fill = TRUE)
  
  # MatchCategory (vectorized)
  levelB_res[, MatchCategory :=
               ifelse(.N > 1, "Multiple", "Unique"),
             by = mz
  ]
  
  return(as.data.frame(levelB_res))
}
