

#' model7
#'
#' power analysis of model 7 in Introduction to Mediation, Moderation, and Conditional Process Analysis
#'
#' @param beta1 regression coefficient of mediator (m) on predictor (x)
#' @param beta2 regression coefficient of outcome (y) on predictor (x)
#' @param beta3 regression coefficient of outcome (y) on mediator (m)
#' @param a1 regression coefficient of mediator (m) on moderator (w)
#' @param c1 regression coefficient of mediator (m) on the product (xw)
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
#' @param ncore number of cores to use in parallel
#' @param pop.cov covariance matrix, default to NULL if using the regression coefficient approach
#' @param mu mean vector, default to NULL if using the regression coefficient approach
#' @param varnames name of variables for the covariance matrix
#' @return power of indirect effect, direct effect, and moderation
#' @export
#' @examples
#' # usage of wp.modmed.m7
#' test = wp.modmed.m7(beta1 = 0.39, beta2 = 0.2, beta3 = 0.3, a1 = 0.39,
#'           c1 = 0.2, sigx2 = 1, sigw2 = 1, sige12 = 1,
#'          sige22 = 1, sigx_w = 0.5, n = 500, nrep = 200,
#'          alpha = 0.05, b = 1000, ncore = 5)
#' print(test)
wp.modmed.m7 <- function (beta1 = 0.39, beta2 = 0.2, beta3 = 0.3, a1 = 0.39,
                  c1 = 0.2, sigx2 = 1, sigw2 = 1, sige12 = 1,
                  sige22 = 1, sigx_w = 0.5, n = 500, nrep = 1000,
                  alpha = 0.05, b = 1000, nb = n, w_value = 0, method="value",
                  ncore = 5, pop.cov = NULL, mu = NULL,
                  varnames = c('y', 'x', 'w', 'm', 'xw'))
{

  if (is.null(pop.cov) || is.null(mu)){
    sigxw2 = sigx2*sigw2 + sigx_w^2
    sigm_x = beta1*sigx2 + a1*sigx_w
    sigm_w = beta1*sigx_w + a1*sigw2
    sigm_xw = c1*sigxw2
    sigm2 = beta1^2*sigx2 + a1^2*sigw2 + c1^2*sigxw2 + sige12 + 2*beta1*a1*sigx_w
    sigy_x = (beta2 + beta3*beta1)*sigx2 + beta3*a1*sigx_w
    sigy_w = (beta2 + beta3*beta1)*sigx_w + beta3*a1*sigw2
    sigy_xw = beta3*c1*sigxw2
    sigy_m = beta2*sigm_x + beta3*sigm2
    sigy2 = (beta2 + beta3*beta1)^2*sigx2 + beta3^2*a1^2*sigw2 + beta3^2*c1^2*sigxw2 + beta3^2*sige12 + sige22 + 2*beta3*a1*(beta2 + beta3*beta1)*sigx_w

    pop.cov = array(c(sigy2, sigy_x, sigy_w, sigy_m, sigy_xw, beta3*sige12,sige22,
                    sigy_x, sigx2, sigx_w, sigm_x, 0, 0, 0,
                    sigy_w, sigx_w, sigw2, sigm_w, 0, 0, 0,
                    sigy_m, sigm_x, sigm_w, sigm2, sigm_xw, sige12, 0,
                    sigy_xw, 0, 0, sigm_xw, sigxw2, 0, 0,
                    beta3*sige12, 0, 0, sige12, 0, sige12, 0,
                    sige22, 0, 0, 0, 0, 0, sige22),
                    dim = c(7, 7))
    pop.cov = pop.cov[1:5, 1:5]
    colnames(pop.cov) = rownames(pop.cov)=c('y','x','w','m','xw')
    # means
    u_xw = sigx_w
    u_m = c1*u_xw
    u_y = beta3*u_m
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
    test_b <- lm(y ~ x + m, data = simdata)

    bootstrap = function(i) {
      boot_dataint = sample.int(n, nb, replace = T)
      boot_data = simdata[boot_dataint, ]
      test_boot1 = lm(m ~ x + w + xw, data = boot_data)
      test_boot2 = lm(y ~ x + m, data = boot_data)
      boot_CI = (test_boot1$coefficients[2] + test_boot1$coefficients[4] *
                   w_value) * test_boot2$coefficients[3]
      boot_DI = as.numeric(test_boot2$coefficients[2])
      boot_mod = as.numeric(test_boot1$coefficients[4])
      boot_CI1 = (test_boot1$coefficients[2] + test_boot1$coefficients[4] *
                    w_value)
      boot_CI2 = test_boot2$coefficients[3]
      return(list(boot_CI, boot_DI, boot_mod, boot_CI1, boot_CI2))
    }
    boot_effect = lapply(1:b, bootstrap)
    boot_CI = matrix(0, ncol = 1, nrow = b)
    boot_CI1 = matrix(0, ncol = 1, nrow = b)
    boot_CI2 = matrix(0, ncol = 1, nrow = b)
    boot_DI = matrix(0, ncol = 1, nrow = b)
    boot_mod = matrix(0, ncol = 1, nrow = b)

    boot_CI = sapply(1:b,function(i) unlist(boot_effect[[i]][1]))
    boot_DI = sapply(1:b,function(i) unlist(boot_effect[[i]][2]))
    boot_mod = sapply(1:b,function(i) unlist(boot_effect[[i]][3]))
    boot_CI1 = sapply(1:b,function(i) unlist(boot_effect[[i]][4]))
    boot_CI2 = sapply(1:b,function(i) unlist(boot_effect[[i]][5]))

    interval_CI = matrix(0, ncol = 1, nrow = 2)
    interval_CI1 = matrix(0, ncol = 1, nrow = 2)
    interval_CI2 = matrix(0, ncol = 1, nrow = 2)
    interval_DI = matrix(0, ncol = 1, nrow = 2)
    interval_mod = matrix(0, ncol = 1, nrow = 2)

    interval_CI[, 1] = quantile(boot_CI,
                                probs = c(alpha / 2, 1 - alpha / 2),
                                names = T)
    interval_CI1[, 1] = quantile(boot_CI1,
                                 probs = c(alpha / 2, 1 - alpha / 2),
                                 names = T)
    interval_CI2[, 1] = quantile(boot_CI2,
                                 probs = c(alpha / 2, 1 - alpha / 2),
                                 names = T)
    interval_DI[, 1] = quantile(boot_DI,
                                probs = c(alpha / 2, 1 - alpha / 2),
                                names = T)
    interval_mod[, 1] = quantile(boot_mod,
                                 probs = c(alpha / 2, 1 - alpha / 2),
                                 names = T)

    r_CI = as.numeric(!dplyr::between(0, interval_CI[1, 1], interval_CI[2, 1]))
    r_DI = as.numeric(!dplyr::between(0, interval_DI[1, 1], interval_DI[2, 1]))
    r_mod = as.numeric(!dplyr::between(0, interval_mod[1, 1], interval_mod[2, 1]))
    if (method == "joint") {
      r_CI = as.numeric(!dplyr::between(0, interval_CI1[1, 1], interval_CI1[2, 1]))*as.numeric(!dplyr::between(0, interval_CI2[1, 1], interval_CI2[2, 1]))
    }
    power = c(r_CI, r_DI, r_mod)
    return(power)
  }
  CL1 = parallel::makeCluster(ncore)
  parallel::clusterExport(CL1, c('beta1', 'beta2', 'beta3', 'a1', 'c1',
                      'sigx2', 'sigw2', 'sige12', 'sige22', 'sigx_w',
                      'n', 'nrep', 'alpha', 'b', 'nb', 'pop.cov',
                      'method', 'mu', 'w_value'), envir = environment())
  allsim <- parallel::parLapply(CL1,1:nrep, runonce)
  parallel::clusterExport(CL1, 'allsim', envir = environment())
  allsim1 = t(parallel::parSapply(CL1, 1:nrep, function(i) unlist(allsim[[i]])))
  power <- colMeans(allsim1)
  parallel::stopCluster(CL1)
  power.structure = structure(list(n = n,
                                 alpha = alpha,
                                 samples = nrep,
                                 w = w_value,
                                 power1 = power[1],
                                 power2 = power[2],
                                 power3 = power[3],
                                 method = "moderated mediation model 7",
                                 url = "https://webpower.psychstat.org/models/modmed7/",
                 note = "power1 is the power of the conditional indirect effect of x on y through m.
power2 is the power of the direct effect of x on  y.
power3 is the power of moderation on the path x to m."), class = "webpower")
  return(power.structure)
}






