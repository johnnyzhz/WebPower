

#' model8
#'
#' power analysis of model 8 in Introduction to Mediation, Moderation, and Conditional Process Analysis
#'
#' @param beta1 regression coefficient of mediator (m) on predictor (x)
#' @param beta2 regression coefficient of outcome (y) on predictor (m)
#' @param beta3 regression coefficient of outcome (y) on mediator (x)
#' @param a1 regression coefficient of mediator (m) on moderator (w)
#' @param a2 regression coefficient of mediator (y) on moderator (w)
#' @param c1 regression coefficient of mediator (m) on the product (xw)
#' @param c2 regression coefficient of mediator (y) on the product (xw)
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
#' # usage of wp.modmed.m8
#' test = wp.modmed.m8(beta1 = 0.2, beta2 = 0.2, beta3 = 0.2,
#'                    a1 = 0.2, a2=0.2, c1 = 0.2, c2 = 0.2,
#'                    sigx2 = 1, sigw2 = 1, sige12 = 1, sige22 = 1, sigx_w = 0.5,
#'                    n = 50, nrep = 100, alpha = 0.05, b = 1000, ncore = 1)
#' print(test)
wp.modmed.m8 <- function (beta1 = 0.2, beta2 = 0.75, beta3 = 0.8,
                  a1 = 0.6, a2 = 0.3, c1 = 0.8, c2 = 0.79,
                  sigx2 = 1, sigw2 = 1, sige12 = 1, sige22 = 1, sigx_w = 0.5,
                  n = 100, nrep = 1000, alpha = 0.05, b = 1000, nb = n,
                  w_value = 0, method = "value", ncore = 1,
                  pop.cov = NULL, mu = NULL, varnames = c('y','x','w','m','xw'))
{
  if (is.null(pop.cov) || is.null(mu)){
    sigxw2 = sigx2*sigw2 + sigx_w^2
    sigm_xw = c1*sigxw2
    sigy_xw = (beta2*c1 + c2)*sigxw2
    sigm_x = a1*sigx_w + beta1*sigx2
    sigm_w = beta1*sigx_w+a1*sigw2
    sigy_x = (a2 + a1*beta2)*sigx_w + (beta3 + beta1*beta2)*sigx2
    sigy_w = (a2+a1*beta2)*sigw2 + (beta3 + beta1*beta2)*sigx_w
    sigy2 = (a2 + a1*beta2)^2*sigw2 + (beta3 + beta1*beta2)^2*sigx2 + beta2^2*sige12 +
      (beta2*c1 + c2)^2*sigxw2 + sige22 + 2*(a2 + a1*beta2)*(beta3 + beta1*beta2)*sigx_w
    sigm2 = beta1^2*sigx2 + a1^2*sigw2 + c1^2*sigxw2 + sige12 + 2*beta1*a1*sigx_w
    sigy_m = a2*sigm_w + beta3*sigm_x + beta2*sigm2 + c2*sigm_xw
    pop.cov = array(c(sigy2, sigy_x, sigy_w, sigy_m, sigy_xw, beta2*sige12, sige22,
                    sigy_x, sigx2, sigx_w, sigm_x, 0, 0, 0,
                    sigy_w, sigx_w, sigw2, sigm_w, 0, 0, 0,
                    sigy_m, sigm_x, sigm_w, sigm2, sigm_xw, sige12, 0,
                    sigy_xw, 0, 0, sigm_xw, sigxw2, 0, 0,
                    beta2*sige12, 0, 0, sige12, 0, sige12, 0,
                    sige22, 0, 0, 0, 0, 0, sige22),
                    dim = c(7, 7))
    pop.cov = pop.cov[1:5,1:5]
    u_xw = sigx_w
    u_m = c1*u_xw
    u_y = beta2*u_m + c2*u_xw
    colnames(pop.cov) = rownames(pop.cov) = c('y','x','w','m','xw')
    mu = c(u_y, 0, 0, u_m, u_xw)
  }else{
    pop.cov = pop.cov
    mu = mu
    colnames(pop.cov) = varnames
  }


  runonce <- function(i){
    simdata <- MASS::mvrnorm(n, mu = mu, Sigma = pop.cov)
    simdata <- as.data.frame(simdata)
    test_a <- lm(m ~ x + w + xw, data = simdata)
    test_b <- lm(y ~ x + m+w+xw, data = simdata)

    bootstrap=function(i){
      boot_dataint = sample.int(n, nb, replace = T)
      boot_data = simdata[boot_dataint,]
      test_boot1 = lm(m ~ x + w + xw,data = boot_data)
      test_boot2 = lm(y ~ x + m + w + xw,data = boot_data)
      boot_CI = (test_boot1$coefficients[2] + test_boot1$coefficients[4]*w_value)*test_boot2$coefficients[3]
      boot_CD = test_boot2$coefficients[2] + test_boot2$coefficients[5]*w_value
      boot_c1 = as.numeric(test_boot1$coefficients[4])
      boot_c2 = as.numeric(test_boot2$coefficients[5])
      boot_CI1 = (test_boot1$coefficients[2] + test_boot1$coefficients[4]*w_value)
      boot_CI2 = test_boot2$coefficients[3]
      return(list(boot_CI, boot_CD, boot_c1, boot_c2, boot_CI1, boot_CI2))
    }
    boot_effect = lapply(1:b, bootstrap)
    boot_CI = matrix(0, ncol = 1, nrow = b)
    boot_CD = matrix(0, ncol = 1, nrow = b)
    boot_c1 = matrix(0, ncol = 1, nrow = b)
    boot_c2 = matrix(0, ncol = 1, nrow = b)
    boot_CI1 = matrix(0, ncol = 1, nrow = b)
    boot_CI2 = matrix(0, ncol = 1, nrow = b)

    boot_CI = t(sapply(1:b, function(i) unlist(boot_effect[[i]][1])))
    boot_CD = t(sapply(1:b, function(i) unlist(boot_effect[[i]][2])))
    boot_c1 = t(sapply(1:b, function(i) unlist(boot_effect[[i]][3])))
    boot_c2 = t(sapply(1:b, function(i) unlist(boot_effect[[i]][4])))
    boot_CI1 = t(sapply(1:b, function(i) unlist(boot_effect[[i]][5])))
    boot_CI2 = t(sapply(1:b, function(i) unlist(boot_effect[[i]][6])))


    interval_CI = matrix(0, ncol = 1, nrow = 2)
    interval_CD = matrix(0, ncol = 1, nrow = 2)
    interval_c1 = matrix(0, ncol = 1, nrow = 2)
    interval_c2 = matrix(0, ncol = 1, nrow = 2)
    interval_CI1 = matrix(0, ncol = 1, nrow = 2)
    interval_CI2 = matrix(0, ncol = 1, nrow = 2)


    interval_CI[, 1] = quantile(boot_CI,
                                probs = c(alpha / 2, 1 - alpha / 2),
                                names = T)
    interval_CD[, 1] = quantile(boot_CD,
                                probs = c(alpha / 2, 1 - alpha / 2),
                                names = T)
    interval_c1[, 1] = quantile(boot_c1,
                                probs = c(alpha / 2, 1 - alpha / 2),
                                names = T)
    interval_c2[, 1] = quantile(boot_c2,
                                probs = c(alpha / 2, 1 - alpha / 2),
                                names = T)
    interval_CI1[, 1] = quantile(boot_CI1,
                                 probs = c(alpha / 2, 1 - alpha / 2),
                                 names = T)
    interval_CI2[, 1] = quantile(boot_CI2,
                                 probs = c(alpha / 2, 1 - alpha / 2),
                                 names = T)


    r_CI = as.numeric(!sapply(1,function(i) dplyr::between(0,interval_CI[1,i], interval_CI[2,i])))
    r_DI = as.numeric(!sapply(1,function(i) dplyr::between(0,interval_CD[1,i], interval_CD[2,i])))
    r_c1 = as.numeric(!sapply(1, function(i) dplyr::between(0,interval_c1[1,i], interval_c1[2,i])))
    r_c2 = as.numeric(!sapply(1, function(i) dplyr::between(0,interval_c2[1,i], interval_c2[2,i])))
    if (method == "joint"){
      r_CI = as.numeric(!dplyr::between(0, interval_CI1[1,1], interval_CI1[2,1]))*as.numeric(!dplyr::between(0, interval_CI2[1,1], interval_CI2[2,1]))
    }
    power=c(r_CI, r_DI, r_c1, r_c2)
    return(power)
  }
  if (ncore > 1){
    CL1 = parallel::makeCluster(ncore)
    parallel::clusterExport(CL1, c('beta1', 'beta2', 'beta3', 'a1', 'a2', 'c1','c2',
                                   'sigx2', 'sigw2', 'sige12', 'sige22', 'sigx_w',
                                   'n', 'nrep', 'alpha', 'b', 'nb', 'pop.cov',
                                   'mu', 'w_value', 'method'), envir = environment())
    allsim <- parallel::parLapply(CL1, 1:nrep, runonce)
    parallel::clusterExport(CL1, 'allsim', envir=environment())
    allsim1 = t(parallel::parSapply(CL1, 1:nrep, function(i) unlist(allsim[[i]])))
    power <- colMeans(allsim1)
    parallel::stopCluster(CL1)
  }else{
    allsim <- sapply(1:nrep, runonce)
    power <- colMeans(t(allsim))
  }

  power.structure = structure(list(n = n,
                                 alpha = alpha,
                                 samples = nrep,
                                 w = w_value,
                                 power1 = power[1],
                                 power2 = power[2],
                                 power3 = power[3],
                                 power4 = power[4],
                 method = "moderated mediation model 8",
                 url = "https://webpower.psychstat.org/models/modmed8/",
                 note="power1 is  the power of the conditional indirect effect of x on y through m.
power2 is the power value of the conditional direct effect of x on y.
power3 is the power of moderation on the path x to m.
power4 is the power of moderation on the path x to y."), class = "webpower")
  return(power.structure)

}
