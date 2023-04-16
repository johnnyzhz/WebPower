

#' model14
#'
#' power analysis of model 14 in Introduction to Mediation, Moderation, and Conditional Process Analysis
#'
#' @param a1 regression coefficient of mediator (m) on predictor (x)
#' @param cc regression coefficient of outcome (y) on predictor (x)
#' @param b1 regression coefficient of outcome (y) on mediator (m)
#' @param c1 regression coefficient of outcome (y) on moderator (w)
#' @param c2 regression coefficient of outcome (y) on the product (mw)
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
#' test = wp.modmed.m14(a1 = 0.2, cc = 0.2, b1 = 0.5, c1 = 0.5, c2 = 0.2, sigx2 = 1,
#'                     sigw2 = 1, sige12 = 1, sige22 = 1, sigx_w = 0.5, n = 100,
#'                     nrep = 50, alpha = 0.05, b = 1000, ncore = 5)
#' print(test)
wp.modmed.m14 <- function (a1 = 0.2, cc = 0.2, b1 = 0.5, c1 = 0.5, c2 = 0.2,
                           sigx2 = 1, sigw2 = 1, sige12 = 1, sige22 = 1,
                           sigx_w = 0.5, n = 100, nrep = 1000, alpha = 0.05,
                           b = 1000, nb = n, w_value = 0, method = "value",
                           ncore = 5, pop.cov = NULL, mu = NULL,
                           varnames =  c('y', 'x', 'w', 'm', 'mw'))
{

  if (is.null(pop.cov) || is.null(mu)) {
    sigm2 = a1^2*sigx2 + sige12
    sigxw2 = sigx2*sigw2 + sigx_w^2
    sigy2 = (cc + a1*b1)^2*sigx2 + c1^2*sigw2 + c2^2*sige12 *
      sigw2 + (a1*c2)^2*sigxw2 + b1^2*sige12 + sige22 + 2*(cc + a1*b1)*c1*sigx_w
    sigmw2 = a1^2*sigxw2 + sige12*sigw2
    sigx_m = a1*sigx2
    sigx_y = (cc + a1*b1)*sigx2 + c1*sigx_w
    sigm_y = a1*(cc + a1*b1)*sigx2 + a1*c1*sigx_w + b1*sige12
    sigm_w = a1*sigx_w
    sigy_w = (cc + a1*b1)*sigx_w + c1*sigw2
    sigy_mw = a1^2*c2*sigxw2 + c2*sige12*sigw2

    # y, x, w, m, mw
    pop.cov = array(
      c(sigy2, sigx_y, sigy_w, sigm_y, sigy_mw,
        sigx_y, sigx2, sigx_w, sigx_m, 0,
        sigy_w, sigx_w, sigw2, sigm_w, 0,
        sigm_y, sigx_m, sigm_w, sigm2, 0,
        sigy_mw, 0, 0, 0, sigmw2),
        dim = c(5, 5)
    )

    u_mw = a1*sigx_w
    u_y = a1*c2*sigx_w
    colnames(pop.cov) = rownames(pop.cov) = c('y', 'x', 'w', 'm', 'mw')
    mu = c(u_y, 0, 0, 0, u_mw)
  }else{
    pop.cov = pop.cov
    mu = mu
    colnames(pop.cov) = varnames
  }

  runonce <- function(i) {
    simdata <- MASS::mvrnorm(n, mu = mu, Sigma = pop.cov)
    simdata <- as.data.frame(simdata)
    test_a <- lm(m ~ x, data = simdata)
    test_b <- lm(y ~ x + m + w + mw, data = simdata)

    bootstrap = function(i) {
      boot_dataint = sample.int(n, nb, replace = T)
      boot_data = simdata[boot_dataint, ]
      test_boot1 = lm(m ~ x, data = boot_data)
      test_boot2 = lm(y ~ x + m + w + mw, data = boot_data)
      boot_CI = test_boot1$coefficients[2]*(test_boot2$coefficients[3] + test_boot2$coefficients[5]*w_value)
      boot_CD = test_boot2$coefficients[2]
      boot_c2 = as.numeric(test_boot2$coefficients[5])
      boot_CI1 = test_boot1$coefficients[2]
      boot_CI2 = (test_boot2$coefficients[3] + test_boot2$coefficients[5]*w_value)
      return(list(boot_CI, boot_CD, boot_c2, boot_CI1, boot_CI2))
    }
    boot_effect = lapply(1:b, bootstrap)
    boot_CI = matrix(0, ncol = 1, nrow = b)
    boot_CD = matrix(0, ncol = 1, nrow = b)
    boot_c2 = matrix(0, ncol = 1, nrow = b)
    boot_CI1 = matrix(0, ncol = 1, nrow = b)
    boot_CI2 = matrix(0, ncol = 1, nrow = b)

    boot_CI = t(sapply(1:b, function(i) unlist(boot_effect[[i]][1])))
    boot_CD = t(sapply(1:b, function(i) unlist(boot_effect[[i]][2])))
    boot_c2 = t(sapply(1:b, function(i) unlist(boot_effect[[i]][3])))
    boot_CI1 = t(sapply(1:b, function(i) unlist(boot_effect[[i]][4])))
    boot_CI2 = t(sapply(1:b, function(i) unlist(boot_effect[[i]][5])))


    interval_CI = matrix(0, ncol = 1, nrow = 2)
    interval_CI1 = matrix(0, ncol = 1, nrow = 2)
    interval_CI2 = matrix(0, ncol = 1, nrow = 2)
    interval_CD = matrix(0, ncol = 1, nrow = 2)
    interval_c2 = matrix(0, ncol = 1, nrow = 2)

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


    r_CI = as.numeric(!sapply(1, function(i) dplyr::between(0, interval_CI[1, i], interval_CI[2, i])))
    r_CD = as.numeric(!sapply(1, function(i) dplyr::between(0, interval_CD[1, i], interval_CD[2, i])))
    r_c2 = as.numeric(!sapply(1, function(i) dplyr::between(0, interval_c2[1, i], interval_c2[2, i])))
    if (method == "joint") {
      r_CI = as.numeric(!dplyr::between(0, interval_CI1[1, 1], interval_CI1[2, 1]))*as.numeric(!dplyr::between(0, interval_CI2[1, 1], interval_CI2[2, 1]))
    }
    power = c(r_CI, r_CD, r_c2)
    return(power)
  }


  CL1 = parallel::makeCluster(ncore)

  parallel::clusterExport(
    CL1,
    c('a1', 'cc', 'b1', 'c1', 'c2', 'sigx2', 'sigw2',
      'sige12', 'sige22', 'sigx_w', 'n', 'nrep',
      'alpha', 'b','nb', 'pop.cov', 'mu', 'method'
    ),
    envir = environment()
  )
  allsim <- parallel::parLapply(CL1, 1:nrep, runonce)
  parallel::clusterExport(CL1, 'allsim', envir = environment())
  allsim1 = t(parallel::parSapply(CL1, 1:nrep, function(i)
    unlist(allsim[[i]])))
  power <- colMeans(allsim1)
  parallel::stopCluster(CL1)
  power.structure = structure(
    list(
      n = n,
      alpha = alpha,
      samples = nrep,
      w = w_value,
      power1 = power[1],
      power2 = power[2],
      power3 = power[3],
      method = "moderated mediation model 14",
      url = "https://webpower.psychstat.org/models/modmed14/",
      note = "power1 is the power of the conditional indirect effect of x on y through m.
power2 is the power value of the direct effect of x on y.
power3 is the power of moderation on the path m to y."
    ),
    class = "webpower"
  )
  return(power.structure)

}



