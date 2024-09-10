
#' model15
#'
#' power analysis of model 15 in Introduction to Mediation, Moderation, and Conditional Process Analysis. Powers are obtained through either the percentile bootstrap method or the Monte Carlo method. The conditional indirect effect value is a1(b1 + b2w); the conditional direct effect value is cp + d2w; the index of moderated mediation is c2b1.
#'
#' @param a1 regression coefficient of mediator (m) on predictor (x)
#' @param cp regression coefficient of outcome (y) on predictor (x)
#' @param b1 regression coefficient of outcome (y) on mediator (m)
#' @param b2 regression coefficient of outcome (y) on the product (mw)
#' @param d1 regression coefficient of outcome (y) on moderator (w)
#' @param d2 regression coefficient of outcome (y) on the product (xw)
#' @param sigx2 variance of predictor (x)
#' @param sigw2 variance of moderator (w)
#' @param sige12 variance of error in the first regression equation
#' @param sige22 variance of error in the second regression equation
#' @param sigx_w covariance between predictor (x) and moderator (w)
#' @param n sample size
#' @param nrep_power number of replications for finding power
#' @param alpha type 1 error rate
#' @param b number of bootstrap iterations 
#' @param nb bootstrap sample size, default to n, used when simulation method is "percentile"
#' @param w_value moderator level, value of w
#' @param power_method "product" for using the indirect effect value in power calculation, or "joint" for using joint significance in power calculation
#' @param simulation_method "percentile" for using percentile bootstrap CI in finding significance of mediation, or "MC" for using Monte Carlo CI in finding significance of mediation
#' @param ncore number of cores to use for the percentile bootstrap method, default is 1, when ncore > 1, parallel is used
#' @param pop.cov covariance matrix, default to NULL if using the regression coefficient approach
#' @param mu mean vector, default to NULL if using the regression coefficient approach
#' @param varnames name of variables for the covariance matrix
#' @return power of indirect effect, direct effect, moderation, and the index of moderated mediation
#' @export
#' @examples
#' test = wp.modmed.m15(a1 = 0.6, cp = 0.2, b1 = 0.3, b2 = 0.2, d1 = 0.2, d2 = 0.1,
#'                      sigx2 = 1, sigw2 = 1, sige12 = 1, sige22 = 1, sigx_w = 0.4,
#'                      w_value = 0.3, simulation_method = "MC",
#'                      n = 50, nrep_power = 1000, b = 1000, alpha = 0.05, ncore = 1)
#'print(test)
wp.modmed.m15 <- function(a1, cp, b1, b2, d1, d2, sige12, sige22, sigx_w, n,
                   sigx2 = 1, sigw2 = 1, 
                   nrep_power = 1000, alpha = 0.05, b = 1000, nb = n,
                   w_value = 0, power_method = "product",
                   simulation_method = "percentile", ncore = 1,
                   pop.cov = NULL, mu = NULL, 
                   varnames = c('y','x','w','m','xw','mw')){
  mmindex_theoretical = a1*b2
  if (is.null(pop.cov) || is.null(mu)){
    sigxw2 = sigx2*sigw2 + sigx_w^2
    sigm_xw = 0
    sigy_xw = (d2 + a1*b2)*sigxw2
    sigm_x = a1*sigx2
    sigm_w = a1*sigx_w
    sigy_x = d1*sigx_w + (cp + a1*b1)*sigx2
    sigy_w = d1*sigw2 + (cp + a1*b1)*sigx_w
    sige1w2 = sige12*sigw2
    sigy2 = (cp + a1*b1)^2*sigx2 + (d2 + a1*b2)^2*sigxw2 + b2 ^
      2*sige1w2 + d1^2*sigw2 + 2*d1*(cp + a1*b1)*sigx_w + b1^2 *
      sige12 + sige22
    sigm2 = a1^2*sigx2 + sige12
    sigy_m = a1*(cp + a1*b1)*sigx2 + a1*d1*sigx_w + b1*sige12
    sigy_mw = a1*(d2 + a1*b2)*sigxw2 + b2*sige1w2
    sigxw_mw = a1*sigxw2
    sigmw2 = a1^2*sigxw2 + sige1w2

    pop.cov=array(c(sigy2, sigy_x, sigy_w, sigy_m, sigy_xw, sigy_mw, b1*sige12, sige22,
                    sigy_x, sigx2, sigx_w, sigm_x, 0, 0, 0, 0,
                    sigy_w, sigx_w, sigw2, sigm_w, 0, 0, 0, 0,
                    sigy_m, sigm_x, sigm_w, sigm2, sigm_xw, 0, sige12, 0,
                    sigy_xw, 0, 0, sigm_xw, sigxw2, sigxw_mw, 0, 0,
                    sigy_mw, 0, 0, 0, sigxw_mw, sigmw2, 0, 0,
                    b1*sige12, 0, 0, sige12, 0, 0, sige12, 0,
                    sige22, 0, 0, 0, 0, 0, 0, sige22),
                    dim = c(8, 8))

    pop.cov = pop.cov[1:6, 1:6]
    rownames(pop.cov) = colnames(pop.cov) = c('y', 'x', 'w', 'm', 'xw', 'mw')
    u_xw = sigx_w
    u_mw = sigm_w
    u_y = d2*u_xw + b2*u_mw
    mu = c(u_y, 0, 0, 0, u_xw, u_mw)
  }else{
    pop.cov = pop.cov
    mu = mu
    colnames(pop.cov) = varnames
  }


  runonce <- function(i){
    if (simulation_method == "percentile"){
      simdata <- MASS::mvrnorm(n, mu = mu, Sigma = pop.cov)
      simdata <- as.data.frame(simdata)
      test_a <- lm(m ~ x, data = simdata)
      test_b <- lm(y ~ x + m + w + xw + mw, data = simdata)
      
      bootstrap=function(i){
        boot_dataint = sample.int(n, nb, replace = T)
        boot_data = simdata[boot_dataint, ]
        test_boot1 = lm(m ~ x, data = boot_data)
        test_boot2 = lm(y ~ x + m + w + xw + mw, data = boot_data)
        boot_CI = (test_boot2$coefficients[3] + test_boot2$coefficients[6]*w_value)*test_boot1$coefficients[2]
        boot_CD = test_boot2$coefficients[2] + test_boot2$coefficients[5]*w_value
        boot_d2 = test_boot2$coefficients[5]
        boot_b2 = test_boot2$coefficients[6]
        boot_CI1 = test_boot1$coefficients[2]
        boot_CI2 = (test_boot2$coefficients[3] + test_boot2$coefficients[6]*w_value)
        boot_index = test_boot2$coefficients[6]*test_boot1$coefficients[2]
        return(list(boot_CI, boot_CD, boot_d2, boot_b2, boot_CI1, boot_CI2, boot_index))
      }
      boot_effect = lapply(1:b, bootstrap)
      boot_CI = matrix(0, ncol = 1, nrow = b)
      boot_CI1 = matrix(0, ncol = 1, nrow = b)
      boot_CI2 = matrix(0, ncol = 1, nrow = b)
      boot_CD = matrix(0, ncol = 1, nrow = b)
      boot_d2 = matrix(0, ncol = 1, nrow = b)
      boot_b2 = matrix(0, ncol = 1, nrow = b)
      boot_index = matrix(0, ncol = 1, nrow = b)
      
      boot_CI=t(sapply(1:b,function(i) unlist(boot_effect[[i]][1])))
      boot_CD=t(sapply(1:b,function(i) unlist(boot_effect[[i]][2])))
      boot_d2=t(sapply(1:b,function(i) unlist(boot_effect[[i]][3])))
      boot_b2=t(sapply(1:b,function(i) unlist(boot_effect[[i]][4])))
      boot_CI1=t(sapply(1:b,function(i) unlist(boot_effect[[i]][5])))
      boot_CI2=t(sapply(1:b,function(i) unlist(boot_effect[[i]][6])))
      boot_index=t(sapply(1:b,function(i) unlist(boot_effect[[i]][7])))
      
      interval_CI=matrix(0,ncol=1,nrow=2)
      interval_CI1=matrix(0,ncol=1,nrow=2)
      interval_CI2=matrix(0,ncol=1,nrow=2)
      interval_CD=matrix(0,ncol=1,nrow=2)
      interval_d2=matrix(0,ncol=1,nrow=2)
      interval_b2=matrix(0,ncol=1,nrow=2)
      interval_index=matrix(0,ncol=1,nrow=2)
      
      interval_CI[, 1] = quantile(boot_CI,
                                  probs = c(alpha / 2, 1 - alpha / 2),
                                  names = T)
      interval_CI1[, 1] = quantile(boot_CI1,
                                   probs = c(alpha / 2, 1 - alpha / 2),
                                   names = T)
      interval_CI2[, 1] = quantile(boot_CI2,
                                   probs = c(alpha / 2, 1 - alpha / 2),
                                   names = T)
      interval_CD[, 1] = quantile(boot_CD,
                                  probs = c(alpha / 2, 1 - alpha / 2),
                                  names = T)
      interval_d2[, 1] = quantile(boot_d2,
                                  probs = c(alpha / 2, 1 - alpha / 2),
                                  names = T)
      interval_b2[, 1] = quantile(boot_b2,
                                  probs = c(alpha / 2, 1 - alpha / 2),
                                  names = T)
      interval_index[, 1] = quantile(boot_index,
                                  probs = c(alpha / 2, 1 - alpha / 2),
                                  names = T)
      
      r_CI=as.numeric(!sapply(1,function(i) dplyr::between(0,interval_CI[1,i],interval_CI[2,i])))
      r_DI=as.numeric(!sapply(1,function(i) dplyr::between(0,interval_CD[1,i],interval_CD[2,i])))
      r_d2=as.numeric(!sapply(1:1,function(i) dplyr::between(0,interval_d2[1,i],interval_d2[2,i])))
      r_b2=as.numeric(!sapply(1:1,function(i) dplyr::between(0,interval_b2[1,i],interval_b2[2,i])))
      r_index=as.numeric(!sapply(1:1,function(i) dplyr::between(0,interval_index[1,i],interval_index[2,i])))
      if (power_method =="joint"){
        r_CI=as.numeric(!dplyr::between(0,interval_CI1[1,1],interval_CI1[2,1]))*as.numeric(!dplyr::between(0,interval_CI2[1,1],interval_CI2[2,1]))
      }
    }else if (simulation_method == "MC"){
      ncore = 1
      simdata <- MASS::mvrnorm(n, mu = mu, Sigma = pop.cov)
      simdata <- as.data.frame(simdata)
      
      model <- '
        m ~ x
        y ~ x + m + w + xw + mw
      '
      
      fit <- lavaan::sem(model = model, data = simdata)
      covm <- lavaan::vcov(fit)
      means <- lavaan::coef(fit)
      
      simmc <- MASS::mvrnorm(b, mu = means, Sigma = covm)
      
      path1_dist <- simmc[,1]
      path2_dist <- simmc[,3] + simmc[,6]*w_value
      med_dist <- path1_dist*path2_dist
      b2_dist <- simmc[,6]
      d2_dist <- simmc[,5]
      cp_dist <- simmc[,2] + w_value*simmc[,5]
      index_dist <-simmc[,6]*simmc[,1]
      path1_interval <- quantile(path1_dist, probs = c(alpha / 2, 1 - alpha / 2))
      path2_interval <- quantile(path2_dist, probs = c(alpha / 2, 1 - alpha / 2))
      med_interval <- quantile(med_dist, probs = c(alpha / 2, 1 - alpha / 2))
      d2_interval <- quantile(d2_dist, probs = c(alpha / 2, 1 - alpha / 2))
      b2_interval <- quantile(b2_dist, probs = c(alpha / 2, 1 - alpha / 2))
      cp_interval <- quantile(cp_dist, probs = c(alpha / 2, 1 - alpha / 2))
      index_interval <- quantile(index_dist, probs = c(alpha / 2, 1 - alpha / 2))
      
      r_CI = as.numeric(!dplyr::between(0, med_interval[1], med_interval[2]))
      r_DI = as.numeric(!dplyr::between(0, cp_interval[1], cp_interval[2]))
      r_d2 = as.numeric(!dplyr::between(0, d2_interval[1], d2_interval[2]))
      r_b2 = as.numeric(!dplyr::between(0, b2_interval[1], b2_interval[2]))
      r_index = as.numeric(!dplyr::between(0, index_interval[1], index_interval[2]))
      if (power_method == "joint") {
        r_CI = as.numeric(!dplyr::between(0, path1_interval[1], path1_interval[2]))*as.numeric(!dplyr::between(0, path2_interval[1], path2_interval[2]))
      }
    }
    power=c(r_CI,r_DI, r_d2, r_b2, r_index)
    return(power)
  }
  if (ncore > 1){
    CL1=parallel::makeCluster(ncore)
    parallel::clusterExport(CL1,c('a1','c','b1','b2','d1','d2',
                                  'sigx2','sigw2','sige12','sige22','sigx_w',
                                  'n','nrep','alpha','b','nb','pop.cov','u_y','u_xw',
                                  'u_mw', 'mu', 'simulation_method', 'w_value',
                                  'power_method'),envir = environment())
    
    allsim <- parallel::parLapply(CL1,1:nrep_power, runonce)
    parallel::clusterExport(CL1,'allsim',envir=environment())
    allsim1=t(parallel::parSapply(CL1,1:nrep_power,function(i) unlist(allsim[[i]])))
    power <- colMeans(allsim1)
    parallel::stopCluster(CL1)
  }else{
    allsim <- sapply(1:nrep_power, runonce)
    power <- colMeans(t(allsim))
  }
 

  power.structure=structure(list(n = n,
                                 alpha = alpha,
                                 samples = nrep_power,
                                 w = w_value,
                                 power1 = power[1],
                                 power2 = power[2],
                                 power3 = power[3],
                                 power4 = power[4],
                                 power5 = power[5],
                                 indirect = a1*(b1 + b2*w_value),
                                 direct = cp + d2*w_value,
                                 index = a1*b2,
                 method="moderated mediation model 15",
                 url = "https://webpower.psychstat.org/models/modmed15/",
                 note = "power1 is the power value of the conditional indirect effect of x on y.
power2 is the power value of the conditional direct effect of x on y.
power3 is the power of moderation on the path x to y.
power4 is the power of moderation on the path m to y.
power5 is the power of the index of moderated mediation.
indirect is the value of the conditional indirect effect.
direct is the value of the conditional direct effect.
index is the value of the index of moderated mediation."), class = "webpower")
  return(power.structure)
}

# test = wp.modmed.m15(a1 = 0.6, cp = 0.2, b1 = 0.3, b2 = 0.2, d1 = 0.2, d2 = 0.1,
#                      sigx2 = 1, sigw2 = 1, sige12 = 1, sige22 = 1, sigx_w = 0.4,
#                      w_value = 0.3, simulation_method = "MC",
#                      n = 50, nrep_power = 100, b = 100, alpha = 0.05, ncore = 1)
# print(test)

