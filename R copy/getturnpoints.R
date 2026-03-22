getturnpoints <-
function (x, calc.proba = TRUE) 
{
  data <- deparse(substitute(x))
  if (isTRUE((is.null(ncol(x)) == FALSE)[1])) 
    stop("Only one series can be treated at a time")
  x <- as.vector(x)
  n <- length(x)
  diffs <- c(x[1] - 1, x[1:(n - 1)]) != x
  uniques <- x[diffs]
  n2 <- length(uniques)
  poss <- (1:n)[diffs]
  exaequos <- c(poss[2:n2], n + 1) - poss - 1
  if (isTRUE((n2 < 3)[1])) {
    warning("Less than 3 unique values, no calculation!")
    nturns <- NA
    firstispeak <- FALSE
    peaks <- rep(FALSE, n2)
    pits <- rep(FALSE, n2)
    tppos <- NA
    proba <- NA
    info <- NA
  }
  else {
    m <- n2 - 2
    ex <- matrix(uniques[1:m + rep(3:1, rep(m, 3)) - 1], 
                 m)
    peaks <- c(FALSE, apply(ex, 1, max, na.rm = TRUE) == 
                 ex[, 2], FALSE)
    pits <- c(FALSE, apply(ex, 1, min, na.rm = TRUE) == ex[, 
                                                           2], FALSE)
    tpts <- peaks | pits
    if (isTRUE((sum(tpts) == 0)[1])) {
      nturns <- 0
      firstispeak <- FALSE
      peaks <- rep(FALSE, n2)
      pits <- rep(FALSE, n2)
      tppos <- NA
      proba <- NA
      info <- NA
    }
    else {
      tppos <- (poss + exaequos)[tpts]
      tptspos <- (1:n2)[tpts]
      firstispeak <- tptspos[1] == (1:n2)[peaks][1]
      nturns <- length(tptspos)
      if (isTRUE((isTRUE(calc.proba))[1])) {
        if (isTRUE((nturns < 2)[1])) {
          inter <- n2 + 1
          posinter1 <- tptspos[1]
        }
        else {
          inter <- c(tptspos[2:nturns], n2) - c(1, tptspos[1:(nturns - 
                                                                1)]) + 1
          posinter1 <- tptspos - c(1, tptspos[1:(nturns - 
                                                   1)])
        }
        posinter2 <- inter - posinter1
        posinter <- pmax(posinter1, posinter2)
        proba <- 2/(inter * gamma(posinter) * gamma(inter - 
                                                      posinter + 1))
        info <- -log(proba, base = 2)
      }
      else {
        proba = NULL
        info <- NULL
      }
    }
  }
  res <- list(data = data, n = n, points = uniques, pos = (poss + 
                                                             exaequos), exaequos = exaequos, nturns = nturns, firstispeak = firstispeak, 
              peaks = peaks, pits = pits, tppos = tppos, proba = proba, 
              info = info)
  class(res) <- "turnpoints"
  res
}
