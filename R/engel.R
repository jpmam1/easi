engel <- function(object = object, file = FALSE, sd = FALSE, lim.y = FALSE) {
  # Compute and plot the Engel curves for an estimated EASI model.
  #
  # Args:
  #   object: the results of easi::easi().
  #   file: optionally, the base name (without file extension) of a file to
  #         contain plotted Engel curves. The file is written in PDF format in
  #         the current directory. If FALSE (default), no output is generated.
  #   sd: if TRUE, compute, plot and return the 95% confidence interval on the
  #       budget shares.
  #   lim.y: if given, the Y-limits for plots.
  #
  # Returns: a list containing:
  #  w: data frame containing the budget share for every observation.
  #  w.pctile: data frame containing median budget share.
  #  w.pctile.upper: confidence interval, upper bound.
  #  w.pctile.lower: confidence interval, lower bound.

  WDELTA <- ifelse(sd != FALSE & !is.na(sd), TRUE, FALSE)
  fit3sls <- object$fit3sls
  var.soc <- object$var.soc
  log.price <- object$log.price
  neq <- object$neq
  y.power <- object$y.power
  nsoc <- object$nsoc
  interact <- object$interact
  py.inter <- object$py.inter
  zy.inter <- object$zy.inter
  pz.inter <- object$pz.inter
  interpz <- object$interpz
  labels.share <- object$labels.share
  y <- object$y

  n <- length(object$log.exp)

  temp <- intermediate.blocs(object)
  my.array <- temp$my.array
  bjk <- temp$bjk
  P <- temp$P
  w <- temp$w
  Z <- temp$Z
  bjr <- temp$bjr
  gjt <- temp$gjt
  hjt <- temp$hjt
  cc <- temp$cc
  noms <- temp$noms
  lnx <- object$log.exp
  y <- temp$y

  result <- list()

  # Calculation of w_j
  W <- matrix(0, n, neq)

  ajk <- my.array
  for (i in 1:neq) {
    tot3 <- tot4 <- tot5 <- tot6 <- tot7 <- 0

    for (j in 1:y.power) {
      tempo <- bjr[j, i] * y^j
      tot3 <- tot3 + tempo
    }

    if (nsoc) {
      for (j in 1:nsoc) {
        tempo <- gjt[j, i] * Z[, j + 1]
        tot4 <- tot4 + tempo
      }
    }

    if (zy.inter) {
      for (j in 1:nsoc) {
        tempo <- hjt[j, i] * Z[, j + 1] * y
        tot5 <- tot5 + tempo
      }
    }

    if (pz.inter) {
      for (k in 1:neq) {
        for (t in (1:(nsoc + 1))) {
          tempo <- ajk[t, k, i] * Z[, t] * P[, k]
          tot6 <- tot6 + tempo
        }
      }
    }

    if (py.inter) {
      for (k in 1:neq) {
        tempo <- bjk[k, i] * P[, k] * y
        tot7 <- tot7 + tempo
      }
    }

    W[, i] <- cc[i] + tot3 + tot4 + tot5 + tot6 + tot7
  }

  # Labels of W matrix
  W <- cbind(W, 1 - apply(W, 1, sum))
  colnames(W) <- labels.share
  result$w <- W

  # Calculation of standard deviations of the fitted budget shares (delta
  # method)
  if (WDELTA) {
    MAT <- rep(1, n)

    for (i in 1:y.power) {
      MAT <- cbind(MAT, y^i)
    }

    if (nsoc) {
      for (i in 1:nsoc) {
        MAT <- cbind(MAT, Z[, i + 1])
      }
    }

    if (zy.inter) {
      for (i in 1:nsoc) MAT <- cbind(MAT, y * Z[, i + 1])
    }

    for (i in 1:neq) {
      MAT <- cbind(MAT, P[, i])
    }

    if (py.inter) {
      for (i in 1:neq) {
        MAT <- cbind(MAT, y * P[, i])
      }
    }

    if (pz.inter) {
      for (i in interpz) {
        for (j in 1:neq) {
          MAT <- cbind(MAT, Z[, i + 1] * P[, j])
        }
      }
    }

    W_ecart <- matrix(0, n, neq)
    D <- summary(fit3sls)$coefCov

    # index D in chunks of dim_varlist
    d <- seq(0, neq * object$dim_varlist, object$dim_varlist)

    for (i in 1:neq) {
      idx <- (d[i]+1):d[i+1]
      W_ecart[, i] <- sqrt(rowSums((MAT %*% D[idx, idx]) * MAT))
    }

    rm(D, MAT)
  }

  # Engel Curves
  ee <- cut(lnx, breaks = quantile(lnx, seq(0, 1, 0.01)),
            include.lowest = TRUE, labels = 1:100)

  Wm <- matrix(0, 100, neq)
  for (i in 1:100) {
    for (j in 1:neq) {
      Wm[i, j] <- median(W[ee == i, j])
    }
  }

  Wm <- cbind(Wm, 1 - apply(Wm, 1, sum))
  colnames(Wm) <- labels.share
  result$w.pctile <- Wm

  # Confidence intervals for fitted budget shares
  if (WDELTA) {
    Wme <- matrix(0, 100, neq + 1)
    for (i in 1:100) {
      for (j in 1:neq) {
        Wme[i, j] <- median(W_ecart[ee == i, j])
      }
    }

    rm(W_ecart)

    for (i in 1:100) {
      Wme[i, neq + 1] <- sqrt(sum(Wme[i, 1:neq]^2))
    }

    Wmep <- Wm + 1.96 * Wme
    Wmem <- Wm - 1.96 * Wme

    Wmep <- Wmep[(1:20) * 5, ]
    Wmem <- Wmem[(1:20) * 5, ]

    result$w.pctile.upper <- Wmep
    result$w.pctile.lower <- Wmem
  }

  write.file <- ifelse(file != FALSE & !is.na(file), TRUE, FALSE)
  if (write.file) {
    # management of labels.share
    if (length(labels.share) < 2) {
      labels.share <- noms
    }
    limYY <- c()
    if (length(lim.y) < 2) {
      for (i in 1:neq) {
        limYY <- c(limYY, c(0, summary(w[, i])[5]))
      }
    } else {
      limYY <- lim.y
    }

    ss <- seq(1, neq * 2, by = 2)
    # Export of Engel curves in the parent folder under the name 'file'.
    
    pdf(paste0(file, ".pdf"))

    xx <- seq(1, 100, len = 20)

    for (i in 1:neq + 1) {
      # smoothing cubic
      sp <- smooth.spline(c(1:100), Wm[, i], spar = 0.9)
      y.loess <- loess(Wm[, i] ~ c(1:100), span = 0.75,
                       data.frame(xxx = c(1:100), yyy = Wm[, i]))
      y.predict <- predict(y.loess, data.frame(xxx = c(1:100)))

      plot(c(1:100), Wm[, i], xlab = "Percentiles of total expenditure",
           ylab = "Budget shares", col = "green",
           ylim = c(limYY[ss[i]], limYY[ss[i] + 1]))

      if (i <= neq) {
        title(main = labels.share[i])
      } else {
        title(main = "Other goods")
      }

      # plot of the adjustment curve
      lines(predict(sp, xx), col = "red")
      lines(c(1:100), y.predict, col = "blue")
      lines(ksmooth(c(1:100), Wm[, i], "normal", bandwidth = 10), col = "black")

      if (!is.na(WDELTA)) {
        points(c((1:20) * 5), Wmep[, i], pch = "+", cex = 1, col = "violet")
        points(c((1:20) * 5), Wmem[, i], pch = "+", cex = 1, col = "violet")
      }
    }
    
    dev.off()
  }

  return(result)
}
