simple_xms_conf_fast <-
function(stage3_results,
                                 filter.by = "M+H",
                                 max_diff_rt = 10,boostIDs=NA){


  DT <- as.data.table(stage3_results)

  DT[, `:=`(
    mz = as.numeric(mz),
    time = as.numeric(time),
    score = as.numeric(score),
    mean_int_vec = as.numeric(mean_int_vec)
  )]

  if(!"pathway_boosted" %in% names(DT)){
    DT[, pathway_boosted := FALSE]
  } else {
    # Coerce to logical in case fread read it as integer (0/1)
    DT[, pathway_boosted := as.logical(pathway_boosted)]
    DT[is.na(pathway_boosted), pathway_boosted := FALSE]
  }
  ############################################################
  # 1️⃣ CLUSTER SUMMARY
  ############################################################

  cluster_summary <- DT[, {

    #valid <- score > 0
    valid <- !is.na(score)
    if(!any(valid)){

      list(
        best_score      = 0.0,
        n_adducts       = 0L,
        has_primary     = FALSE,
        rt_range        = Inf,
        multimer_valid  = FALSE,
        charge_valid    = FALSE,
        pathway_boosted = FALSE
      )

    } else {

      mono_ind  <- !grepl("^2M|^3M", Adduct) & valid
      multi_ind <-  grepl("^2M|^3M", Adduct) & valid

      mono_int  <- if(any(mono_ind)) max(mean_int_vec[mono_ind], na.rm=TRUE) else NA_real_
      multi_int <- if(any(multi_ind)) max(mean_int_vec[multi_ind], na.rm=TRUE) else NA_real_

      charge_valid <- TRUE
      if("charge" %in% names(.SD)){
        mono_charge <- max(mean_int_vec[valid & charge==1], na.rm=TRUE)
        high_charge <- max(mean_int_vec[valid & charge>1], na.rm=TRUE)
        charge_valid <- is.na(high_charge) | high_charge <= mono_charge
      }

      list(
        best_score      = max(score[valid], na.rm=TRUE),
        n_adducts       = as.integer(uniqueN(Adduct[valid])),
        has_primary     = any(Adduct[valid] %in% filter.by),
        rt_range        = IQR(time[valid]),
        multimer_valid  = is.na(multi_int) | multi_int <= mono_int,
        charge_valid    = charge_valid,
        pathway_boosted  = any(pathway_boosted, na.rm = TRUE)   # ← new
      )
    }

  }, by = .(chemical_ID, Module_RTclust)]

  cluster_summary[, rt_valid := rt_range <= max_diff_rt]

  ############################################################
  # 2️⃣ CONFIDENCE RULES
  ############################################################

  cluster_summary[, Confidence :=
                    fifelse(best_score<=0,0L,
                            fifelse(
                              n_adducts>=2 &
                                has_primary &
                                best_score>=10 &
                                rt_valid &
                                multimer_valid &
                                charge_valid,
                              3L,
                              fifelse(
                                has_primary &
                                  best_score>0 &
                                  pathway_boosted,
                                2L,
                                fifelse(
                                  n_adducts>=2 &
                                    best_score>0,
                                  1L,
                                0L
                              ))))
  ]

  ############################################################
  # 3️⃣ ATTACH CONFIDENCE BACK
  ############################################################

  curated_res <- merge(
    DT,
    cluster_summary[,.(chemical_ID,Module_RTclust,Confidence)],
    by=c("chemical_ID","Module_RTclust"),
    all.x=TRUE
  )


  if (length(boostIDs) > 0 && !all(is.na(boostIDs))) {
    curated_res[
      chemical_ID %in% boostIDs,
      Confidence := 4L
    ]
  }
  ############################################################
  # 4️⃣ REPORTING
  ############################################################

  #cat("Stage 4 confidence distribution (chemical IDs)\n")
  # print(table(unique(curated_res[,.(chemical_ID,Confidence)])$Confidence))
  # Sort rows by confidence (highest first)
  setorder(curated_res, -Confidence)




  return(list(
    cluster_summary = cluster_summary,
    conf_mat = curated_res
  ))
}
