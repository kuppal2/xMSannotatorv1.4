multilevelannotationstep5 <-
function(outloc) {




  DT <- data.table::fread(file.path(outloc,"Stage4.csv"))
  setDT(DT)

  setorder(DT, -Confidence, -score)

  DT[, best_conf := max(Confidence, na.rm=TRUE),
     by=.(mz, time, Module_RTclust)]

  DT <- DT[Confidence == best_conf]

  DT[, best_score := max(score, na.rm=TRUE),
     by=.(mz, time, Module_RTclust)]

  DT <- DT[score == best_score]

  ############################################################
  ############################################################

  final <- unique(DT)

  ############################################################
  ############################################################

  summary <- final[, .(
    max_score = max(score, na.rm=TRUE),
    max_confidence = max(Confidence, na.rm=TRUE),
    n_features = .N
  ), by = chemical_ID]

  setorder(summary, -max_confidence, -max_score)
  chem_conf <-DT[
    , list(Confidence = max(Confidence, na.rm = TRUE)),
    by = chemical_ID
  ]

  cat("Stage 5 confidence distribution (unique chemical IDs)\n")

  print(
    chem_conf[
      , list(N = .N),
      by = Confidence
    ][order(-Confidence)]
  )

  ############################################################
  ############################################################

  fwrite(final, file.path(outloc,"Stage5.csv"))
  fwrite(summary, file.path(outloc,"Stage5_summary.csv"))

  print("Stage5 redundancy filtering complete")

  return(list(
    feature_annotations = as.data.frame(final),
    chemical_summary = as.data.frame(summary)
  ))
}
