multilevelannotationstep5 <-
function(outloc) {
  
  
  setwd(outloc)
  
  DT <- data.table::fread("Stage4_curated.csv")
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
  
  ############################################################
  ############################################################
  
  fwrite(final, "Stage5_feature_annotations.csv")
  fwrite(summary, "Stage5_summary.csv")
  
  print("Stage5 redundancy filtering complete")
  
  return(list(
    feature_annotations = as.data.frame(final),
    chemical_summary = as.data.frame(summary)
  ))
}
