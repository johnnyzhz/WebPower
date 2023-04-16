
#' model58
#'
#' power analysis of model 58 in Introduction to Mediation, Moderation, and Conditional Process Analysis
#'
#' @param a1 regression coefficient of outcome (m) on moderator (w)
#' @param b1 regression coefficient of mediator (m) on predictor (x)
#' @param c1 regression coefficient of outcome (m) on the product (xw)
#' @param a2 regression coefficient of outcome (y) on moderator (w)
#' @param b2 regression coefficient of outcome (y) on mediator (m)
#' @param c2 regression coefficient of outcome (y) on the product (mw)
#' @param b3 regression coefficient of outcome (y) on predictor (x)
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
#' test = wp.modmed.m58(a1 = 0.2, b1 = 0.2, c1 = 0.1, c2 = 0.1,
#'      b2 = 0.2, b3 = 0.2, a2 = 0.2,
#'      sigx2 = 1, sigw2 = 1, sige12 = 1, sige22 = 1, sigx_w = 0.5,
#'      n = 100, nrep = 50, alpha = 0.05, b = 1000, ncore = 5)
#' print(test)
wp.modmed.m58 <- function(a1 = 0.5, b1 = 0.75, c1 = 0.1, a2 = 0.4, b2 = 0.6, c2 = 0.1,
                   b3 = 0.5, sigx2 = 1, sigw2 = 1, sige12 = 1, sige22 = 1, sigx_w = 0.4,
                   n = 500, nrep = 1000, alpha = 0.05, b=1000, nb = n,
                   w_value = 0, method = "value", ncore=5,
                   pop.cov = NULL, mu = NULL)
{
  if (is.null(pop.cov) || is.null(mu)){
    sigx_m = b1*sigx2 + a1*sigx_w
    sigx_mw = c1*(1 + 2*sigx_w^2 / sigx2 / sigw2)*sigx2*sigw2
    sigx_y = a2*sigx_w + b3*sigx2 + b2*sigx_m + c2*sigx_mw
    sigx_xw = sigw_xw = 0

    sigw2 = sigw2
    sigw_m = b1*sigx_w + a1*sigw2
    sigw_mw = 3*c1*sigx_w / sqrt(sigx2) / sqrt(sigw2)*sqrt(sigx2)*(sqrt(sigw2))^3
    sigw_y = a2*sigw2 + b3*sigx_w + b2*sigw_m + c2*sigw_mw

    sigxw2 = sigx2*sigw2 + sigx_w^2
    sigxw_w2 = 3*sigx_w*sigw2 - sigx_w*sigw2

    sigm2 = b1^2*sigx2 + a1^2*sigw2 + c1^2*sigxw2 + sige12 +
      2*b1*a1*sigx_w
    sigm_xw = c1*sigxw2
    sigm_mw = b1*c1*(sigx2*sigw2 + 2*sigx_w^2) + 3*a1*c1*sigx_w*sigw2 + c1*b1*sigxw2 + c1*a1*sigxw_w2
    sigm_y = (a2 + a1*b2)*sigw_m + (b3 + b1*b2)*sigx_m + b2*c1*sigm_xw + c2*sigm_mw + b2*sige12

    sigxw_mw = b1*sigxw2 + a1*(3*sigx_w*sigw2 - sigx_w*sigw2)
    sigxw_y = (a2 + a1*b2)*sigw_xw + (b3 + b1*b2)*sigx_xw + b2*c1*sigxw2 + c2*sigxw_mw

    sigxw22 = 3*sigx2*sigw2^2 - 3*sigx_w^2*sigw2 + 15*sigx_w^2*sigw2
    sige1w2 = sige12*sigw2

    sigmw2 = b1^2*sigxw2 + 2*a1^2*sigw2^2 + c1^2*sigxw22 + sige1w2 + 2*b1*a1*sigxw_w2
    sigmw_y = (a2 + a1*b2)*sigw_mw + (b3 + b1*b2)*sigx_mw + b2*c1*sigxw_mw + c2*sigmw2

    sigy2 = (a2 + a1*b2)^2*sigw2 + (b3 + b1*b2)^2*sigx2 + (b2*c1)^2*sigxw2 + c2^2*sigmw2 + b2^2*sige12 + sige22 + 2*(a2 + a1*b2)*(b3 + b1*b2)*sigx_w + 2*(a2 + a1*b2)*b2*c1*sigw_xw + 2*(a2 + a1*b2)*c2*sigw_mw + 2*(b3 + b1*b2)*c2*sigx_mw + 2*b2*c1*c2*sigxw_mw


    pop.cov=array(c(sigx2, sigx_w, sigx_m, sigx_xw, sigx_mw, sigx_y, 0, 0,
                    sigx_w, sigw2, sigw_m, sigw_xw, sigw_mw, sigw_y, 0, 0,
                    sigx_m, sigw_m, sigm2, sigm_xw, sigm_mw, sigm_y, sige12, 0,
                    sigx_xw, sigw_xw, sigm_xw, sigxw2, sigxw_mw, sigxw_y, 0, 0,
                    sigx_mw, sigw_mw, sigm_mw, sigxw_mw, sigmw2, sigmw_y, 0, 0,
                    sigx_y, sigw_y, sigm_y, sigxw_y, sigmw_y, sigy2, b2*sige12, sige22,
                    0, 0, sige12, 0, 0, b2*sige12, sige12, 0,
                    0, 0, 0, 0, 0, sige22, 0, sige22), dim=c(8, 8))

    pop.cov = pop.cov[1:6, 1:6]
    u_xw = sigx_w
    u_m = c1*u_xw
    u_mw = b1*u_xw + a1*sigw2
    u_y = b2*u_m + c2*u_mw
    mu = c(0, 0, u_m, u_xw, u_mw, u_y)
    rownames(pop.cov) = colnames(pop.cov) = c('x', 'w', 'm', 'xw', 'mw', 'y')
  }else{
    pop.cov = pop.cov
    mu = mu
    colnames(pop.cov) = varnames
  }

  ## conduct the analysis once
  ##bootstrap sampling

  runonce <- function(i){
    simdata <- MASS::mvrnorm(n, mu = mu, Sigma = pop.cov)
    simdata <- as.data.frame(simdata)
    test_a <- lm(m ~ x + w + xw, data = simdata)
    test_b <- lm(y ~ x + m + w + mw, data = simdata)

    bootstrap=function(i){
      boot_dataint = sample.int(n, nb, replace = T)
      boot_data = simdata[boot_dataint, ]
      test_boot1 = lm(m ~ x + w + xw, data = boot_data)
      test_boot2 = lm(y ~ x + m + w + mw, data = boot_data)
      boot_CI = (test_boot1$coefficients[2] + test_boot1$coefficients[4]*w_value)*(test_boot2$coefficients[3] + test_boot2$coefficients[5]*w_value)
      boot_CI1 = (test_boot1$coefficients[2] + test_boot1$coefficients[4]*w_value)
      boot_CI2 = (test_boot2$coefficients[3] + test_boot2$coefficients[5]*w_value)
      boot_DI = test_boot2$coefficients[2]
      boot_c1 = test_boot1$coefficients[4]
      boot_c2 = test_boot2$coefficients[5]
      return(list(boot_CI, boot_DI, boot_c1, boot_c2, boot_CI1, boot_CI2))
    }
    boot_effect = lapply(1:b, bootstrap)
    boot_CI = matrix(0, ncol = 1, nrow = b)
    boot_CI1 = matrix(0, ncol = 1, nrow = b)
    boot_CI2 = matrix(0, ncol = 1, nrow = b)
    boot_DI = matrix(0, ncol = 1, nrow = b)
    boot_c1 = matrix(0, ncol = 1, nrow = b)
    boot_c2 = matrix(0, ncol = 1, nrow = b)

    boot_CI = t(sapply(1:b,function(i) unlist(boot_effect[[i]][1])))
    boot_DI = as.matrix(sapply(1:b,function(i) unlist(boot_effect[[i]][2])))
    boot_c1 = as.matrix(sapply(1:b,function(i) unlist(boot_effect[[i]][3])))
    boot_c2 = as.matrix(sapply(1:b,function(i) unlist(boot_effect[[i]][4])))
    boot_CI1 = t(sapply(1:b,function(i) unlist(boot_effect[[i]][5])))
    boot_CI2 = t(sapply(1:b,function(i) unlist(boot_effect[[i]][6])))

    interval_CI = matrix(0, ncol = 1, nrow = 2)
    interval_DI = matrix(0, ncol = 1, nrow = 2)
    interval_c1 = matrix(0, ncol = 1, nrow = 2)
    interval_c2 = matrix(0, ncol = 1, nrow = 2)
    interval_CI1 = matrix(0, ncol = 1, nrow = 2)
    interval_CI2 = matrix(0, ncol = 1, nrow = 2)

    interval_CI[, 1] = quantile(boot_CI,
                                probs = c(alpha / 2, 1 - alpha / 2),
                                names = T)
    interval_CI1[, 1] = quantile(boot_CI1,
                                 probs = c(alpha / 2, 1 - alpha / 2),
                                 names = T)
    interval_CI2[, 1] = quantile(boot_CI2,
                                 probs = c(alpha / 2, 1 - alpha / 2),
                                 names = T)
    interval_DI[, 1] = quantile(boot_DI[, 1],
                                probs = c(alpha / 2, 1 - alpha / 2),
                                names = T)
    interval_c1[, 1] = quantile(boot_c1[, 1],
                                probs = c(alpha / 2, 1 - alpha / 2),
                                names = T)
    interval_c2[, 1] = quantile(boot_c2[, 1],
                                probs = c(alpha / 2, 1 - alpha / 2),
                                names = T)

    r_CI = as.numeric(!sapply(1,function(i) dplyr::between(0,interval_CI[1,i],interval_CI[2,i])))
    r_DI = as.numeric(!sapply(1,function(i) dplyr::between(0,interval_DI[1,i],interval_DI[2,i])))
    r_c1 = as.numeric(!sapply(1,function(i) dplyr::between(0,interval_c1[1,i],interval_c1[2,i])))
    r_c2 = as.numeric(!sapply(1,function(i) dplyr::between(0,interval_c2[1,i],interval_c2[2,i])))
    power = c(r_CI, r_DI, r_c1, r_c2)
    if (method == "joint"){
      r_CI = as.numeric(!dplyr::between(0,interval_CI1[1,1], interval_CI1[2,1]))*as.numeric(!dplyr::between(0, interval_CI2[1,1], interval_CI2[2,1]))
    }
    return(power)
  }
  CL1 = parallel::makeCluster(ncore)
  parallel::clusterExport(CL1,c('a1', 'b1', 'c1', 'c2', 'b2', 'b3', 'a2',
                      'sigx2', 'sigw2', 'sige12', 'sige22', 'sigx_w',
                      'n', 'nrep', 'alpha','b','nb','pop.cov',
                      'mu', 'method', 'w_value'),envir = environment())

  allsim <- parallel::parLapply(CL1, 1:nrep, runonce)
  parallel::clusterExport(CL1, 'allsim', envir = environment())
  allsim1 = t(parallel::parSapply(CL1, 1:nrep, function(i) unlist(allsim[[i]])))
  power <- colMeans(allsim1)
  parallel::stopCluster(CL1)

  power.structure=structure(list(n = n,
                                 alpha = alpha,
                                 samples = nrep,
                                 w = w_value,
                                 power1 = power[1],
                                 power2 = power[2],
                                 power3 = power[3],
                                 power4 = power[4],
                 method = "moderated mediation model 58",
                 url = "https://webpower.psychstat.org/models/modmed58/",
                 note="power1 is the power of the conditional indirect effect of x on y through m.
power2 is  the power of the direct effect of x on y.
power3 is the power of moderation on the path x to m.
power4 is the powe of moderation on the path m to y."), class = "webpower")
  return(power.structure)

}




