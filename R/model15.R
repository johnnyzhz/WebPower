
#' model15
#'
#' power analysis of model 15 in Introduction to Mediation, Moderation, and Conditional Process Analysis
#'
#' @param a1 regression coefficient of mediator (m) on predictor (x)
#' @param c regression coefficient of outcome (y) on predictor (x)
#' @param b1 regression coefficient of outcome (y) on mediator (m)
#' @param b2 regression coefficient of outcome (y) on the product (mw)
#' @param c1 regression coefficient of outcome (y) on moderator (w)
#' @param c2 regression coefficient of outcome (y) on the product (xw)
#' @param sigx2 variance of predictor (x)
#' @param sigw2 variance of moderator (w)
#' @param sige12 variance of error in the first regression equation
#' @param sige22 variance of error in the second regression equation
#' @param sigx_w covariance between predictor (x) and moderator (w)
#' @param n sample size
#' @param nrep number of replications
#' @param alpha type 1 error rate
#' @param b number of bootstrap iterations
#' @param nb bootstrap sample size, default to n
#' @param w_value moderator level
#' @param method "value" for using the indirect effect value in power calculation, or "joint" for using joint significance in power calculation
#' @param ncore number of cores to use, default is 1, when ncore > 1, parallel is used
#' @param pop.cov covariance matrix, default to NULL if using the regression coefficient approach
#' @param mu mean vector, default to NULL if using the regression coefficient approach
#' @param varnames name of variables for the covariance matrix
#' @return power of indirect effect, direct effect, and moderation
#' @export
#' @examples
#' test = wp.modmed.m15(a1 = 0.6, c = 0.2, b1 = 0.3, b2 = 0.2, c1 = 0.2, c2 = 0.1,
#'                      sigx2 = 1, sigw2 = 1, sige12 = 1, sige22 = 1, sigx_w = 0.4,
#'                      n = 50, nrep = 100, alpha = 0.05, b = 1000, ncore = 1)
#'print(test)
wp.modmed.m15 <- function(a1 = 0.6, c = 0.2, b1 = 0.3, b2 = 0.5, c1 = 0.2, c2 = 0.1,
                   sigx2 = 1, sigw2 = 1, sige12 = 1, sige22 = 1, sigx_w = 0.4,
                   n = 500, nrep = 1000, alpha = 0.05, b = 1000, nb = n,
                   w_value = 0, method = "value", ncore = 1,
                   pop.cov = NULL, mu = NULL,
                   varnames = c('y','x','w','m','xw','mw'))
{
  if (is.null(pop.cov) || is.null(mu)){
    sigxw2 = sigx2*sigw2 + sigx_w^2
    sigm_xw = 0
    sigy_xw = (c2 + a1*b2)*sigxw2
    sigm_x = a1*sigx2
    sigm_w = a1*sigx_w
    sigy_x = c1*sigx_w + (c + a1*b1)*sigx2
    sigy_w = c1*sigw2 + (c + a1*b1)*sigx_w
    sige1w2 = sige12*sigw2
    sigy2 = (c + a1*b1)^2*sigx2 + (c2 + a1*b2)^2*sigxw2 + b2 ^
      2*sige1w2 + c1^2*sigw2 + 2*c1*(c + a1*b1)*sigx_w + b1^2 *
      sige12 + sige22
    sigm2 = a1^2*sigx2 + sige12
    sigy_m = a1*(c + a1*b1)*sigx2 + a1*c1*sigx_w + b1*sige12
    sigy_mw = a1*(c2 + a1*b2)*sigxw2 + b2*sige1w2
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
    u_y = c2*u_xw + b2*u_mw
    mu = c(u_y, 0, 0, 0, u_xw, u_mw)
  }else{
    pop.cov = pop.cov
    mu = mu
    colnames(pop.cov) = varnames
  }


  runonce <- function(i){
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
      boot_c2 = test_boot2$coefficients[5]
      boot_b2 = test_boot2$coefficients[6]
      boot_CI1 = test_boot1$coefficients[2]
      boot_CI2 = (test_boot2$coefficients[3] + test_boot2$coefficients[6]*w_value)
      return(list(boot_CI, boot_CD, boot_c2, boot_b2, boot_CI1, boot_CI2))
    }
    boot_effect = lapply(1:b, bootstrap)
    boot_CI = matrix(0, ncol = 1, nrow = b)
    boot_CI1 = matrix(0, ncol = 1, nrow = b)
    boot_CI2 = matrix(0, ncol = 1, nrow = b)
    boot_CD = matrix(0, ncol = 1, nrow = b)
    boot_c2 = matrix(0, ncol = 1, nrow = b)
    boot_b2 = matrix(0, ncol = 1, nrow = b)

    boot_CI=t(sapply(1:b,function(i) unlist(boot_effect[[i]][1])))
    boot_CD=t(sapply(1:b,function(i) unlist(boot_effect[[i]][2])))
    boot_c2=t(sapply(1:b,function(i) unlist(boot_effect[[i]][3])))
    boot_b2=t(sapply(1:b,function(i) unlist(boot_effect[[i]][4])))
    boot_CI1=t(sapply(1:b,function(i) unlist(boot_effect[[i]][5])))
    boot_CI2=t(sapply(1:b,function(i) unlist(boot_effect[[i]][6])))

    interval_CI=matrix(0,ncol=1,nrow=2)
    interval_CI1=matrix(0,ncol=1,nrow=2)
    interval_CI2=matrix(0,ncol=1,nrow=2)
    interval_CD=matrix(0,ncol=1,nrow=2)
    interval_c2=matrix(0,ncol=1,nrow=2)
    interval_b2=matrix(0,ncol=1,nrow=2)

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
    interval_c2[, 1] = quantile(boot_c2,
                                probs = c(alpha / 2, 1 - alpha / 2),
                                names = T)
    interval_b2[, 1] = quantile(boot_b2,
                                probs = c(alpha / 2, 1 - alpha / 2),
                                names = T)


    r_CI=as.numeric(!sapply(1,function(i) dplyr::between(0,interval_CI[1,i],interval_CI[2,i])))
    r_DI=as.numeric(!sapply(1,function(i) dplyr::between(0,interval_CD[1,i],interval_CD[2,i])))
    r_c2=as.numeric(!sapply(1:1,function(i) dplyr::between(0,interval_c2[1,i],interval_c2[2,i])))
    r_b2=as.numeric(!sapply(1:1,function(i) dplyr::between(0,interval_b2[1,i],interval_b2[2,i])))
    if (method =="joint"){
      r_CI=as.numeric(!dplyr::between(0,interval_CI1[1,1],interval_CI1[2,1]))*as.numeric(!dplyr::between(0,interval_CI2[1,1],interval_CI2[2,1]))
    }
    power=c(r_CI,r_DI, r_c2, r_b2)
    return(power)
  }
  if (ncore > 1){
    CL1=parallel::makeCluster(ncore)
    parallel::clusterExport(CL1,c('a1','c','b1','b2','c1','c2',
                                  'sigx2','sigw2','sige12','sige22','sigx_w',
                                  'n','nrep','alpha','b','nb','pop.cov','u_y','u_xw',
                                  'u_mw', 'mu', 'method', 'w_value'),envir = environment())
    
    allsim <- parallel::parLapply(CL1,1:nrep, runonce)
    parallel::clusterExport(CL1,'allsim',envir=environment())
    allsim1=t(parallel::parSapply(CL1,1:nrep,function(i) unlist(allsim[[i]])))
    power <- colMeans(allsim1)
    parallel::stopCluster(CL1)
  }else{
    allsim <- sapply(1:nrep, runonce)
    power <- colMeans(t(allsim))
  }
 

  power.structure=structure(list(n = n,
                                 alpha = alpha,
                                 samples = nrep,
                                 w = w_value,
                                 power1 = power[1],
                                 power2 = power[2],
                                 power3 = power[3],
                                 power4 = power[4],
                 method="moderated mediation model 15",
                 url = "https://webpower.psychstat.org/models/modmed15/",
                 note = "power1 is the power value of the conditional indirect effect of x on y.
power2 is the power value of the conditional direct effect of x on y.
power3 is the power of moderation on the path x to y.
power4 is the power of moderation on the path m to y."), class = "webpower")
  return(power.structure)
}


