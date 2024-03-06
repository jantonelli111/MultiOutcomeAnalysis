###########################################################
## Code for analysis of multiple treatments              ##
## and multiple outcomes simultaneously                  ##
##                                                       ##
##                                                       ##
##                                                       ##
## Inputs to the function:                               ##
##                                                       ##
## 1: Y: an N by Q matrix of outcomes                    ##
## 2: Tr: an N by P matrix of exposures                  ##
## 3: X: an N by C matrix of confounders                 ##
## 4: b: J by Q matrix for NC contrasts                  ##
## 5: t1: exposure level 1 for estimand                  ##
## 6: t2: exposure level 2 for estimand                  ##
## 7: t1NC: list of exposure 1 values for NC             ##
## 8: t2NC: list of exposure 2 values for NC             ##
###########################################################

# Our main function
multiFunc <- function(Y, Tr, X, b, t1, t2, t1NC, t2NC, maxM){
  
  # Check no. of observations, outcomes, exposures, and observed confounders
  N = nrow(Y)
  Q = ncol(Y)
  P = ncol(Tr)
  C = ncol(X)
  
  # Check no. of NC outcomes
  J = ncol(b)
  
  # Estimate no. of unobserved confounders for Tr and Y models
  Mest = matrix(NA, 2, 8)
  rownames(Mest) = c("Tr", "Y")
  colnames(Mest) = c("vss.comp1", "vss.comp2", "MAP", "BIC", "adjBIC",
                     "parallel.fa", "parallel.pc", "eigen")
  
  naivef = earth::earth(x = X, y = Tr)
  TrResidual = Tr - predict(naivef)
  naiveg = earth::earth(x = cbind(Tr, X), y = Y)
  YResidual = Y - predict(naiveg)
  
  ## Store residual densities
  densYlist = list()
  for (qq in 1 : Q) {
    densYlist[[qq]] = density(YResidual[,qq])
  }
  
  densTrlist = list()
  for (pp in 1 : P) {
    densTrlist[[pp]] = density(TrResidual[,pp])
  }
  
  invisible(capture.output(vssTr <- psych::vss(TrResidual, plot = FALSE)))
  invisible(capture.output(parallelTr <- psych::fa.parallel(TrResidual,
                                                            plot = FALSE)))
  Mest[1, ] = c(which.max(vssTr$cfit.1), which.max(vssTr$cfit.2),
                which.min(vssTr$map),
                which.min(vssTr$vss.stats[, "BIC"]),
                which.min(vssTr$vss.stats[, "SABIC"]),
                parallelTr$nfact, parallelTr$ncomp,
                sum(eigen(cor(TrResidual))$values > 1))
  
  invisible(capture.output(vssY <- psych::vss(YResidual, plot = FALSE)))
  invisible(capture.output(parallelY <- psych::fa.parallel(YResidual,
                                                           plot = FALSE)))
  Mest[2, ] = c(which.max(vssY$cfit.1), which.max(vssY$cfit.2),
                which.min(vssY$map),
                which.min(vssY$vss.stats[, "BIC"]),
                which.min(vssY$vss.stats[, "SABIC"]),
                parallelY$nfact, parallelY$ncomp,
                sum(eigen(cor(YResidual))$values > 1))
  
  # Save p-values from the test if the given M is sufficient
  # to capture the full dimensionality of data.
  Mpval = matrix(NA, 2, maxM)
  rownames(Mpval) = c("Tr", "Y")
  colnames(Mpval) = as.factor(1:maxM)
  for(mm in 1:maxM){
    Mpval[1, mm] = as.numeric(factanal(TrResidual, factors = mm, 
                                       nstart=100, lower = 0.01)$PVAL)
    Mpval[2, mm] = as.numeric(factanal(YResidual, factors = mm, 
                                       nstart=100, lower = 0.01)$PVAL)
  }
  
  # For screeplots
  TrScree = data.frame(M = rep(1:P, 2),
                       Type = rep(c("PC", "FA"), each = P),
                       VarExp = c(parallelTr$pc.values, parallelTr$fa.values))
  TrScree$Type = factor(TrScree$Type, levels = c("PC", "FA"))
  YScree = data.frame(M = rep(1:Q, 2),
                      Type = rep(c("PC", "FA"), each = Q),
                      VarExp = c(parallelY$pc.values, parallelY$fa.values))
  YScree$Type = factor(YScree$Type, levels = c("PC", "FA"))
  
  # Determine M - We try multiple values of m and save the results.
  saveResult = array(list(), maxM)
  for(mm in 1:maxM){
    
    M = mm
    
    # Print what simulation number we're on
    #cat(paste("M =", mm, "; outcome of interest =", aa), "\r")
    cat(paste("M =", mm), "\r")
    
    #--------------------------------------------------------------------------#
    # Estimation of B.hat, Gamma.hat, Sigma_u.tx.hat, Sigma_y.tx.hat
    #--------------------------------------------------------------------------#
    # Perform factor analysis on residuals of Tr model
    TrResidualFA = factanal(TrResidual, factors = M,
                            nstart=100, lower = 0.01)
    B.hat = sqrt(diag(var(TrResidual))) * TrResidualFA$loadings[]
    sigma2_t.xu.hat = mean(TrResidualFA$uniquenesses * diag(var(TrResidual)))
    Sigma_t.x.hat = B.hat %*% t(B.hat) + diag(sigma2_t.xu.hat, P)
    
    coef_mu_u.tx.hat = t(B.hat) %*% solve(Sigma_t.x.hat)
    mu_u.deltatx.hat <- coef_mu_u.tx.hat %*% (t1 - t2)
    Sigma_u.tx.hat = diag(M) - t(B.hat) %*% solve(Sigma_t.x.hat) %*% B.hat
    Sigma_u.tx.hat_chol = chol(solve(Sigma_u.tx.hat))
    
    # Perform factor analysis on residuals of Y model
    YResidualFA = factanal(YResidual, factors = M,
                           nstart=100, lower = 0.01)
    Gamma.hat = sqrt(diag(var(YResidual))) * YResidualFA$loadings[]
    sigma2_y.txu.hat = mean(YResidualFA$uniquenesses * diag(var(YResidual)))
    Sigma_y.tx.hat = Gamma.hat %*% t(Gamma.hat) + diag(sigma2_y.txu.hat, Q)
    
    G = as.numeric(colMeans(predict(naiveg, cbind(t(replicate(N, t1)), X)) -
                              predict(naiveg, cbind(t(replicate(N, t2)), X))))
    
    #--------------------------------------------------------------------------#
    # Negative Controls
    #--------------------------------------------------------------------------#
    D = NULL
    GhatList = list()
    MncList = list()
    MncInvList = list()
    for(jj in 1:J){
      # For jj-th NC outcome
      t1NC_jj = t1NC[[jj]]
      t2NC_jj = t2NC[[jj]]
      
      Ghat_j = Mnc_j = NULL
      for(cj in 1:nrow(t1NC_jj)){
        t1NC_cj = t1NC_jj[cj, ]
        t2NC_cj = t2NC_jj[cj, ]
        tmp1 = colMeans(predict(naiveg, cbind(t(replicate(N, t1NC_cj)), X)))
        tmp2 = colMeans(predict(naiveg, cbind(t(replicate(N, t2NC_cj)), X)))
        Ghat_j_onecol = as.numeric(tmp1 - tmp2)
        Ghat_j = cbind(Ghat_j, Ghat_j_onecol)
        Mnc_j_onecol = Sigma_u.tx.hat_chol %*% coef_mu_u.tx.hat %*%
          (t1NC_cj - t2NC_cj)
        Mnc_j = cbind(Mnc_j, Mnc_j_onecol)
        MncInv_j = MASS::ginv(Mnc_j)
        D_j = t(b[, jj, drop = FALSE]) %*% Ghat_j %*% MncInv_j
      }
      D = rbind(D, D_j)
      GhatList[[jj]] = Ghat_j
      MncList[[jj]] = Mnc_j
      MncInvList[[jj]] = MncInv_j
    }
    
    out = list(B.hat = B.hat, Gamma.hat = Gamma.hat,
               mu_u.deltatx.hat = mu_u.deltatx.hat,
               Sigma_u.tx.hat = Sigma_u.tx.hat,
               Sigma_u.tx.hat_chol = Sigma_u.tx.hat_chol,
               Sigma_t.x.hat = Sigma_t.x.hat,
               Sigma_y.tx.hat = Sigma_y.tx.hat,
               G = G,
               D = D, GhatList = GhatList,
               MncList = MncList, MncInvList = MncInvList)
    
    saveResult[[mm]] = out
    
  }
  
  allResults = list(Dimensions = list(N = N, Q = Q, P = P, C = C, J = J),
                    Settings = list(t1 = t1, t2 = t2,
                                    b = b, t1NC = t1NC, t2NC = t2NC),
                    Estimates = saveResult,
                    Mdetermination = list(Mest = Mest, Mpval = Mpval,
                                          TrScree = TrScree, YScree = YScree),
                    Densities = list(densTrlist = densTrlist,
                                     densYlist = densYlist))
  return(allResults)
  
}
