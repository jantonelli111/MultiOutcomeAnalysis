# Our main function
multiFunc <- function(Y, Tr, X, b, t1, t2, t1NC, t2NC, maxM, scaleData = FALSE, nB = 100){
  
  # Convert to data.table if not already
  Y = as.data.table(Y)
  Tr = as.data.table(Tr)
  X = as.data.table(X)
  
  # Check no. of observations, outcomes, exposures, and observed confounders
  N = nrow(Y)
  Q = ncol(Y)
  P = ncol(Tr)
  C = ncol(X)
  # Check no. of NC outcomes
  J = ncol(b)
  
  # Scale data
  if (scaleData) {
    Tr[] = lapply(Tr, scale, center = FALSE, scale = TRUE)
    Y[] = lapply(Y, scale, center = FALSE, scale = TRUE)
  }
  
  # Setup for matrix of estimation results
  Mest = matrix(NA, 2, 8,
                 dimnames = list(c("Tr", "Y"),
                                 c("vss.comp1", "vss.comp2",
                                   "MAP", "BIC", "adjBIC",
                                   "parallel.fa", "parallel.pc",
                                   "eigen")))
  
  # Modeling and residual calculation
  naivef = earth::earth(x = X, y = Tr)
  TrResidual = Tr - predict(naivef) #Tr -  X[, predict(naivef, newdata = .SD), .SDcols = names(Tr)]
  
  for(gMod in 1:2){
    
    if(gMod == 1){  # gLinear
      naiveg = lm(as.matrix(Y) ~ as.matrix(Tr) + as.matrix(X))
      YResidual = Y - data.table(predict(naiveg, newdata = cbind(Tr, X)))
    } else if(gMod == 2){  # gNonLinear
      naiveg = earth::earth(x = cbind(Tr, X), y = Y)
      YResidual = Y - data.table(predict(naiveg, newdata = cbind(Tr, X)))
    }
    
    ## Store residual densities
    densYlist = lapply(1:Q, function(qq) density(YResidual[[qq]]))
    densTrlist = lapply(1:P, function(pp) density(TrResidual[[pp]]))
    
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
    Mpval = matrix(NA, 2, maxM, dimnames = list(c("Tr", "Y"), as.character(1:maxM)))
    for(mm in 1:maxM){
      Mpval[1, mm] = as.numeric(factanal(TrResidual, factors = mm, 
                                         nstart=100, lower = 0.01)$PVAL)
      Mpval[2, mm] = as.numeric(factanal(YResidual, factors = mm, 
                                         nstart=100, lower = 0.01)$PVAL)
    }
    
    # For screeplots
    TrScree = data.table(M = rep(1:P, each = 2),
                         Type = rep(c("PC", "FA"), each = P),
                         VarExp = c(parallelTr$pc.values, parallelTr$fa.values))
    TrScree[, Type := factor(Type, levels = c("PC", "FA"))]
    
    YScree = data.table(M = rep(1:Q, each = 2),
                        Type = rep(c("PC", "FA"), each = Q),
                        VarExp = c(parallelY$pc.values, parallelY$fa.values))
    YScree[, Type := factor(Type, levels = c("PC", "FA"))]
    
    # Determine M - We try multiple values of m and save the results.
    saveResult = array(list(), maxM)
    Gboot = array(list(), maxM)
    for(mm in 1:maxM){
      
      M = mm
      
      # Print what simulation number we're on
      cat(paste(ifelse(gMod == 1, "g is linear", "g is nonlinear"), "; M =", mm), "\r")
      #cat(paste("M =", mm), "\r")
      
      #--------------------------------------------------------------------------#
      # Estimation of B.hat, Gamma.hat, Sigma_u.tx.hat, Sigma_y.tx.hat
      #--------------------------------------------------------------------------#
      # Perform factor analysis on residuals of Tr model
      TrResidualFA = factanal(TrResidual, factors = M,
                              nstart=100, lower = 0.01)
      B.hat = sqrt(diag(var(TrResidual))) * TrResidualFA$loadings[]
      FATrRes = TrResidualFA$uniquenesses * diag(var(TrResidual))
      sigma2_t.xu.hat = mean(FATrRes)
      Sigma_t.x.hat = B.hat %*% t(B.hat) + diag(sigma2_t.xu.hat, P)
      
      coef_mu_u.tx.hat = t(B.hat) %*% solve(Sigma_t.x.hat)
      #mu_u.deltatx.hat <- coef_mu_u.tx.hat %*% (t1 - t2)
      Sigma_u.tx.hat = diag(M) - t(B.hat) %*% solve(Sigma_t.x.hat) %*% B.hat
      Sigma_u.tx.hat_chol = chol(solve(Sigma_u.tx.hat))
      
      # Perform factor analysis on residuals of Y model
      YResidualFA = factanal(YResidual, factors = M,
                             nstart=100, lower = 0.01)
      Gamma.hat = sqrt(diag(var(YResidual))) * YResidualFA$loadings[]
      FAYRes = YResidualFA$uniquenesses * diag(var(YResidual))
      sigma2_y.txu.hat = mean(FAYRes)
      Sigma_y.tx.hat = Gamma.hat %*% t(Gamma.hat) + diag(sigma2_y.txu.hat, Q)
      
      #G = as.numeric(colMeans(predict(naiveg, cbind(t(replicate(N, t1)), X)) -
      #                          predict(naiveg, cbind(t(replicate(N, t2)), X))))
      
      mu_u.deltatx.hat = list()
      G = list()
      for(jj in 1:length(t1)){
        t1_jj = t1[[jj]]
        t2_jj = t2[[jj]]
        mu_u.deltatx.hat[[jj]] = coef_mu_u.tx.hat %*% (t1_jj - t2_jj)
        G1 = predict(naiveg, as.data.frame(cbind(matrix(t1_jj, nrow = N, ncol = P, byrow = TRUE), X)))
        G2 = predict(naiveg, as.data.frame(cbind(matrix(t2_jj, nrow = N, ncol = P, byrow = TRUE), X)))
        G[[jj]] = as.numeric(colMeans(G1 - G2))
      }
      
      #--------------------------------------------------------------------------#
      # Negative Controls
      #--------------------------------------------------------------------------#
      #D = NULL
      D = matrix(NA, J, M)
      GhatList = list()
      MncList = list()
      MncInvList = list()
      for(jj in 1:J){
        # For jj-th NC outcome
        Ghat_j = Mnc_j = NULL
        for(cj in seq_len(nrow(t1NC[[jj]]))){
          #tmp1 = colMeans(predict(naiveg, cbind(t(replicate(N, t1NC_cj)), X)))
          #tmp2 = colMeans(predict(naiveg, cbind(t(replicate(N, t2NC_cj)), X)))
          t1NC_cj_data = data.table(matrix(t1NC[[jj]][cj, ],
                                           nrow = N, ncol = P, byrow = TRUE), X)
          t2NC_cj_data = data.table(matrix(t2NC[[jj]][cj, ],
                                           nrow = N, ncol = P, byrow = TRUE), X)
          tmp1 = colMeans(predict(naiveg, newdata = t1NC_cj_data))
          tmp2 = colMeans(predict(naiveg, newdata = t2NC_cj_data))
          Ghat_j = cbind(Ghat_j, tmp1 - tmp2)
          Mnc_j = cbind(Mnc_j, Sigma_u.tx.hat_chol %*%
                          coef_mu_u.tx.hat %*%
                          (t1NC[[jj]][cj, ] - t2NC[[jj]][cj, ]))
        }
        #D = rbind(D, D_j)
        D[jj, ] = t(b[, jj, drop = FALSE]) %*% Ghat_j %*% MASS::ginv(Mnc_j)
        GhatList[[jj]] = Ghat_j
        MncList[[jj]] = Mnc_j
        MncInvList[[jj]] = MASS::ginv(Mnc_j)
      }
      
      saveResult[[mm]] = list(B.hat = B.hat, Gamma.hat = Gamma.hat,
                              mu_u.deltatx.hat = mu_u.deltatx.hat,
                              Sigma_u.tx.hat = Sigma_u.tx.hat,
                              Sigma_u.tx.hat_chol = Sigma_u.tx.hat_chol,
                              Sigma_t.x.hat = Sigma_t.x.hat,
                              Sigma_y.tx.hat = Sigma_y.tx.hat,
                              FATrRes = FATrRes, FAYRes = FAYRes,
                              G = G,
                              D = D, GhatList = GhatList,
                              MncList = MncList, MncInvList = MncInvList)
      
      ## Measuring uncertainty of PATE under NUC
      Gboot[[mm]] = array(NA, dim = c(nB, Q, length(t1)))
      
      system.time(for(bt in 1:nB){
        
        # Generate bootstrapped data
        set.seed(bt)
        indB = sample(N, replace = TRUE)
        Trbt = Tr[indB, ]
        Ybt = Y[indB, ]
        Xbt = X[indB, ]
        
        if(gMod == 1){ # gLinear
          naiveg = lm(as.matrix(Ybt) ~ as.matrix(Trbt) + as.matrix(Xbt))
          #YResidual = Y - predict(naiveg)
        } else if(gMod == 2){ # gNonLinear
          naiveg = earth::earth(x = cbind(Trbt, Xbt), y = Ybt)
          #YResidual = Y - predict(naiveg)
        }
        
        for(jj in 1:length(t1)){
          t1_jj_data = data.table(matrix(t1[[jj]], nrow = N, ncol = P, byrow = TRUE), Xbt)
          t2_jj_data = data.table(matrix(t2[[jj]], nrow = N, ncol = P, byrow = TRUE), Xbt)
          G1 = predict(naiveg, newdata = t1_jj_data)
          G2 = predict(naiveg, newdata = t2_jj_data)
          Gboot[[mm]][bt, , jj] = as.numeric(colMeans(G1 - G2))
        }
      })

      if(gMod == 1){
        gLinear = list(Estimates = saveResult,
                       Mdetermination = list(Mest = Mest, Mpval = Mpval,
                                             TrScree = TrScree, YScree = YScree),
                       Densities = list(densTrlist = densTrlist,
                                        densYlist = densYlist),
                       Corr = list(CorrTr = cor(Tr), CorrTrResidual = cor(TrResidual),
                                   CorrY = cor(Y), CorrYResidual = cor(YResidual)),
                       gSummary = summary(naiveg),
                       gvcov = vcov(naiveg),
                       Gboot = Gboot)
      } else if(gMod == 2){
        gNonLinear = list(Estimates = saveResult,
                       Mdetermination = list(Mest = Mest, Mpval = Mpval,
                                             TrScree = TrScree, YScree = YScree),
                       Densities = list(densTrlist = densTrlist,
                                        densYlist = densYlist),
                       Corr = list(CorrTr = cor(Tr), CorrTrResidual = cor(TrResidual),
                                   CorrY = cor(Y), CorrYResidual = cor(YResidual)),
                       Gboot = Gboot)
      }

    }
    
  }
  
  allResults = list(Dimensions = list(N = N, Q = Q, P = P, C = C, J = J),
                    Settings = list(t1 = t1, t2 = t2,
                                    b = b, t1NC = t1NC, t2NC = t2NC),
                    gLinear = gLinear,
                    gNonLinear = gNonLinear)

  return(allResults)
  
}
