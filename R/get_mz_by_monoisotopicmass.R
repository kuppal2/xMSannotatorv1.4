get_mz_by_monoisotopicmass <-
function(
    inputmassmat,
    queryadductlist = c("M+H"),
    adduct_table
) {
  
  cnames <- c("mz","ID","Name","Formula",
              "MonoisotopicMass","Adduct","AdductMass")
  
  # ---- Ensure data.frame ----
  inputmassmat <- as.data.frame(inputmassmat)
  adduct_table <- as.data.frame(adduct_table)
  
  # ---- Extract compound columns safely ----
  exact_mass <- as.numeric(inputmassmat$Monoisotopic_Mass)
  dbid       <- as.character(inputmassmat$Formula_ID)
  name       <- as.character(inputmassmat$Compound_names)
  formula    <- as.character(inputmassmat$Molecular_Formula)
  
  # ---- Extract adduct parameters ----
  adduct_names  <- as.character(adduct_table$Adduct)
  adduct_mass   <- as.numeric(adduct_table$adductMass)
  adduct_charge <- as.numeric(adduct_table$charge)
  adduct_nmol   <- as.numeric(adduct_table$num_molecules)
  adduct_mode   <- as.character(adduct_table$Mode)
  
  names(adduct_mass)   <- adduct_names
  names(adduct_charge) <- adduct_names
  names(adduct_nmol)   <- adduct_names
  
  # ---- Resolve adduct selection ----
  if (queryadductlist[1] == "positive") {
    queryadductlist <- adduct_names[adduct_mode == "positive"]
  } else if (queryadductlist[1] == "negative") {
    queryadductlist <- adduct_names[adduct_mode == "negative"]
  } else if (queryadductlist[1] == "all") {
    queryadductlist <- adduct_names
  } else {
    if (any(!queryadductlist %in% adduct_names)) {
      stop("Invalid adduct in queryadductlist.")
    }
  }
  
  n_comp <- length(exact_mass)
  n_add  <- length(queryadductlist)
  
  # ---- Cartesian expansion ----
  comp_index <- rep(seq_len(n_comp), each = n_add)
  add_index  <- rep(seq_len(n_add), times = n_comp)
  
  sel_mass   <- adduct_mass[queryadductlist]
  sel_charge <- adduct_charge[queryadductlist]
  sel_nmol   <- adduct_nmol[queryadductlist]
  
  # ---- Vectorized m/z computation ----
  mz_vals <- ((exact_mass[comp_index] * sel_nmol[add_index]) +
                sel_mass[add_index]) /
    sel_charge[add_index]
  
  # ---- Build result ----
  result <- data.frame(
    mz = mz_vals,
    ID = dbid[comp_index],
    Name = name[comp_index],
    Formula = formula[comp_index],
    MonoisotopicMass = round(exact_mass[comp_index],5),
    Adduct = queryadductlist[add_index],
    AdductMass = sel_mass[add_index],
    stringsAsFactors = FALSE
  )
  
  colnames(result) <- cnames
  
  unique(result)
}
