group_by_rt_linear <-
function(mchemicaldata, time_step=2, max_diff_rt=10,
                               groupnum) {
  mcd <- as.data.frame(mchemicaldata)
  if (nrow(mcd) == 0L) {
    mcd$Module_RTclust <- character(0)
    return(mcd[, c("Module_RTclust", setdiff(colnames(mcd), "Module_RTclust"))])
  }
  ord    <- order(mcd$time)
  mcd_s  <- mcd[ord, , drop=FALSE]
  clust  <- integer(nrow(mcd_s))
  clust[1L] <- 1L
  if (nrow(mcd_s) > 1L) {
    cur <- 1L
    gaps <- diff(mcd_s$time)
    for (i in seq_along(gaps)) {
      if (gaps[i] > time_step) cur <- cur + 1L
      clust[i + 1L] <- cur
    }
  }
  mcd_s$Module_RTclust <- paste0(groupnum, "_", clust)
  mcd_s <- mcd_s[order(ord), , drop=FALSE]
  other <- setdiff(colnames(mcd_s), "Module_RTclust")
  mcd_s[, c("Module_RTclust", other), drop=FALSE]
}
