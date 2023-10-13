
#' model58
#'
#' power analysis of model 58 in Introduction to Mediation, Moderation, and Conditional Process Analysis
#'
#' @param c1 regression coefficient of outcome (m) on moderator (w)
#' @param a1 regression coefficient of mediator (m) on predictor (x)
#' @param c2 regression coefficient of outcome (m) on the product (xw)
#' @param d1 regression coefficient of outcome (y) on moderator (w)
#' @param b1 regression coefficient of outcome (y) on mediator (m)
#' @param b2 regression coefficient of outcome (y) on the product (mw)
#' @param cp regression coefficient of outcome (y) on predictor (x)
#' @param sigx2 variance of predictor (x)
#' @param sigw2 variance of moderator (w)
#' @param sige12 variance of error in the first regression equation
#' @param sige22 variance of error in the second regression equation
#' @param sigx_w covariance between predictor (x) and moderator (w)
#' @param n sample size
#' @param nrep number of replications for finding power
#' @param alpha type 1 error rate
#' @param b number of bootstrap iterations used when simulation method is "percentile"
#' @param MCrep number of repetitions used for finding distribution when simulation method is "MC"
#' @param nb bootstrap sample size, default to n, used when simulation method is "percentile"
#' @param w_value moderator level
#' @param power_method "product" for using the indirect effect value in power calculation, or "joint" for using joint significance in power calculation
#' @param simulation_method "percentile" for using percentile bootstrap CI in finding significance of mediation, or "MC" for using Monte Carlo CI in finding significance of mediation
#' @param ncore number of cores to use, default is 1, when ncore > 1, parallel is used
#' @param pop.cov covariance matrix, default to NULL if using the regression coefficient approach
#' @param mu mean vector, default to NULL if using the regression coefficient approach
#' @param varnames name of variables for the covariance matrix
#' @return power of indirect effect, direct effect, and moderation
#' @export
#' @examples
#' test = wp.modmed.m58(c1 = 0.2, a1 = 0.2, c2 = 0.1, b2 = 0.1,
#'      b1 = 0.2, cp = 0.2, d1 = 0.2, w_value = 0.3, simulation_method = "MC",
#'      sigx2 = 1, sigw2 = 1, sige12 = 1, sige22 = 1, sigx_w = 0.5,
#'      n = 50, nrep = 1000, alpha = 0.05, ncore = 1)
#' print(test)
wp.modmed.m58 <- function(c1, a1, c2, d1, b1, b2, cp, sige12, sige22, sigx_w, n,
                   sigx2 = 1, sigw2 = 1, 
                   nrep = 1000, alpha = 0.05, b = 1000, nb = n,
                   w_value = 0, power_method = "product", MCrep = 1000,
                   ncore = 1, simulation_method = "percentile",
                   pop.cov = NULL, mu = NULL, varnames =  c('x', 'w', 'm', 'xw', 'mw', 'y'))
{
  if (is.null(pop.cov) || is.null(mu)){
    sigx_m = a1*sigx2 + c1*sigx_w
    sigx_mw = c2*(1 + 2*sigx_w^2 / sigx2 / sigw2)*sigx2*sigw2
    sigx_y = d1*sigx_w + cp*sigx2 + b1*sigx_m + b2*sigx_mw
    sigx_xw = sigw_xw = 0

    sigw2 = sigw2
    sigw_m = a1*sigx_w + c1*sigw2
    sigw_mw = 3*c2*sigx_w / sqrt(sigx2) / sqrt(sigw2)*sqrt(sigx2)*(sqrt(sigw2))^3
    sigw_y = d1*sigw2 + cp*sigx_w + b1*sigw_m + b2*sigw_mw

    sigxw2 = sigx2*sigw2 + sigx_w^2
    sigxw_w2 = 3*sigx_w*sigw2 - sigx_w*sigw2

    sigm2 = a1^2*sigx2 + c1^2*sigw2 + c2^2*sigxw2 + sige12 +
      2*a1*c1*sigx_w
    sigm_xw = c2*sigxw2
    sigm_mw = a1*c2*(sigx2*sigw2 + 2*sigx_w^2) + 3*c1*c2*sigx_w*sigw2 + c2*a1*sigxw2 + c2*c1*sigxw_w2
    sigm_y = (d1 + c1*b1)*sigw_m + (cp + a1*b1)*sigx_m + b1*c2*sigm_xw + b2*sigm_mw + b1*sige12

    sigxw_mw = a1*sigxw2 + c1*(3*sigx_w*sigw2 - sigx_w*sigw2)
    sigxw_y = (d1 + c1*b1)*sigw_xw + (cp + a1*b1)*sigx_xw + b1*c2*sigxw2 + b2*sigxw_mw

    sigxw22 = 3*sigx2*sigw2^2 - 3*sigx_w^2*sigw2 + 15*sigx_w^2*sigw2
    sige1w2 = sige12*sigw2

    sigmw2 = a1^2*sigxw2 + 2*c1^2*sigw2^2 + c2^2*sigxw22 + sige1w2 + 2*a1*c1*sigxw_w2
    sigmw_y = (d1 + c1*b1)*sigw_mw + (cp + a1*b1)*sigx_mw + b1*c2*sigxw_mw + b2*sigmw2

    sigy2 = (d1 + c1*b1)^2*sigw2 + (cp + a1*b1)^2*sigx2 + (b1*c2)^2*sigxw2 + b2^2*sigmw2 + b1^2*sige12 + sige22 + 2*(d1 + c1*b1)*(cp + a1*b1)*sigx_w + 2*(d1 + c1*b1)*b1*c2*sigw_xw + 2*(d1 + c1*b1)*b2*sigw_mw + 2*(cp + a1*b1)*b2*sigx_mw + 2*b1*c2*b2*sigxw_mw


    pop.cov=array(c(sigx2, sigx_w, sigx_m, sigx_xw, sigx_mw, sigx_y, 0, 0,
                    sigx_w, sigw2, sigw_m, sigw_xw, sigw_mw, sigw_y, 0, 0,
                    sigx_m, sigw_m, sigm2, sigm_xw, sigm_mw, sigm_y, sige12, 0,
                    sigx_xw, sigw_xw, sigm_xw, sigxw2, sigxw_mw, sigxw_y, 0, 0,
                    sigx_mw, sigw_mw, sigm_mw, sigxw_mw, sigmw2, sigmw_y, 0, 0,
                    sigx_y, sigw_y, sigm_y, sigxw_y, sigmw_y, sigy2, b1*sige12, sige22,
                    0, 0, sige12, 0, 0, b1*sige12, sige12, 0,
                    0, 0, 0, 0, 0, sige22, 0, sige22), dim=c(8, 8))

    pop.cov = pop.cov[1:6, 1:6]
    u_xw = sigx_w
    u_m = c2*u_xw
    u_mw = a1*u_xw + c1*sigw2
    u_y = b1*u_m + b2*u_mw
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
    if (simulation_method == "percentile"){
      simdata <- MASS::mvrnorm(n, mu = mu, Sigma = pop.cov)
      simdata <- as.data.frame(simdata)
      test_a <- lm(m ~ x + w + xw, data = simdata)
      test_b <- lm(y ~ x + m + w + mw, data = simdata)
      
      bootstrap = function(i){
        boot_dataint = sample.int(n, nb, replace = T)
        boot_data = simdata[boot_dataint, ]
        test_boot1 = lm(m ~ x + w + xw, data = boot_data)
        test_boot2 = lm(y ~ x + m + w + mw, data = boot_data)
        boot_CI = (test_boot1$coefficients[2] + test_boot1$coefficients[4]*w_value)*(test_boot2$coefficients[3] + test_boot2$coefficients[5]*w_value)
        boot_CI1 = (test_boot1$coefficients[2] + test_boot1$coefficients[4]*w_value)
        boot_CI2 = (test_boot2$coefficients[3] + test_boot2$coefficients[5]*w_value)
        boot_DI = test_boot2$coefficients[2]
        boot_c2 = test_boot1$coefficients[4]
        boot_b2 = test_boot2$coefficients[5]
        return(list(boot_CI, boot_DI, boot_c2, boot_b2, boot_CI1, boot_CI2))
      }
      boot_effect = lapply(1:b, bootstrap)
      boot_CI = matrix(0, ncol = 1, nrow = b)
      boot_CI1 = matrix(0, ncol = 1, nrow = b)
      boot_CI2 = matrix(0, ncol = 1, nrow = b)
      boot_DI = matrix(0, ncol = 1, nrow = b)
      boot_c2 = matrix(0, ncol = 1, nrow = b)
      boot_b2 = matrix(0, ncol = 1, nrow = b)
      
      boot_CI = t(sapply(1:b,function(i) unlist(boot_effect[[i]][1])))
      boot_DI = as.matrix(sapply(1:b,function(i) unlist(boot_effect[[i]][2])))
      boot_c2 = as.matrix(sapply(1:b,function(i) unlist(boot_effect[[i]][3])))
      boot_b2 = as.matrix(sapply(1:b,function(i) unlist(boot_effect[[i]][4])))
      boot_CI1 = t(sapply(1:b,function(i) unlist(boot_effect[[i]][5])))
      boot_CI2 = t(sapply(1:b,function(i) unlist(boot_effect[[i]][6])))
      
      interval_CI = matrix(0, ncol = 1, nrow = 2)
      interval_DI = matrix(0, ncol = 1, nrow = 2)
      interval_c2 = matrix(0, ncol = 1, nrow = 2)
      interval_b2 = matrix(0, ncol = 1, nrow = 2)
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
      interval_c2[, 1] = quantile(boot_c2[, 1],
                                  probs = c(alpha / 2, 1 - alpha / 2),
                                  names = T)
      interval_b2[, 1] = quantile(boot_b2[, 1],
                                  probs = c(alpha / 2, 1 - alpha / 2),
                                  names = T)
      
      r_CI = as.numeric(!sapply(1,function(i) dplyr::between(0,interval_CI[1,i],interval_CI[2,i])))
      r_DI = as.numeric(!sapply(1,function(i) dplyr::between(0,interval_DI[1,i],interval_DI[2,i])))
      r_c2 = as.numeric(!sapply(1,function(i) dplyr::between(0,interval_c2[1,i],interval_c2[2,i])))
      r_b2 = as.numeric(!sapply(1,function(i) dplyr::between(0,interval_b2[1,i],interval_b2[2,i])))
      if (power_method == "joint"){
        r_CI = as.numeric(!dplyr::between(0,interval_CI1[1,1], interval_CI1[2,1]))*as.numeric(!dplyr::between(0, interval_CI2[1,1], interval_CI2[2,1]))
      }
    }else if (simulation_method == "MC"){
      simdata <- MASS::mvrnorm(n, mu = mu, Sigma = pop.cov)
      simdata <- as.data.frame(simdata)
      test_a <- lm(m ~ x + w + xw, data = simdata)
      test_b <- lm(y ~ x + m + w + mw, data = simdata)
      
      a1_mean <- summary(test_a)$coefficients[2, 1]
      c1_mean <- summary(test_a)$coefficients[3, 1]
      c2_mean <- summary(test_a)$coefficients[4, 1]
      cp_mean <- summary(test_b)$coefficients[2, 1]
      b1_mean <- summary(test_b)$coefficients[3, 1]
      d1_mean <- summary(test_b)$coefficients[4, 1]
      b2_mean <- summary(test_b)$coefficients[5, 1]
      
      a1_se <- summary(test_a)$coefficients[2, 2]
      c1_se <- summary(test_a)$coefficients[3, 2]
      c2_se <- summary(test_a)$coefficients[4, 2]
      cp_se <- summary(test_b)$coefficients[2, 2]
      b1_se <- summary(test_b)$coefficients[3, 2]
      d1_se <- summary(test_b)$coefficients[4, 2]
      b2_se <- summary(test_b)$coefficients[5, 2]
      
      path1_dist <- rnorm(MCrep, a1_mean, a1_se) + rnorm(MCrep, c2_mean, c2_se)*w_value
      path2_dist <- rnorm(MCrep, b1_mean, b1_se) + rnorm(MCrep, b2_mean, b2_se)*w_value
      med_dist <- path1_dist*path2_dist
      c2_dist <- rnorm(MCrep, c2_mean, c2_se)
      b2_dist <- rnorm(MCrep, b2_mean, b2_se)
      cp_dist <- rnorm(MCrep, cp_mean, cp_se)
      path1_interval <- quantile(path1_dist, probs = c(alpha / 2, 1 - alpha / 2))
      path2_interval <- quantile(path2_dist, probs = c(alpha / 2, 1 - alpha / 2))
      med_interval <- quantile(med_dist, probs = c(alpha / 2, 1 - alpha / 2))
      c2_interval <- quantile(c2_dist, probs = c(alpha / 2, 1 - alpha / 2))
      b2_interval <- quantile(b2_dist, probs = c(alpha / 2, 1 - alpha / 2))
      cp_interval <- quantile(cp_dist, probs = c(alpha / 2, 1 - alpha / 2))
      
      r_CI = as.numeric(!dplyr::between(0, med_interval[1], med_interval[2]))
      r_DI = as.numeric(!dplyr::between(0, cp_interval[1], cp_interval[2]))
      r_c2 = as.numeric(!dplyr::between(0, c2_interval[1], c2_interval[2]))
      r_b2 = as.numeric(!dplyr::between(0, b2_interval[1], b2_interval[2]))
      if (power_method == "joint") {
        r_CI = as.numeric(!dplyr::between(0, path1_interval[1], path1_interval[2]))*as.numeric(!dplyr::between(0, path2_interval[1], path2_interval[2]))
      }
    }
    
    power = c(r_CI, r_DI, r_c2, r_b2)
   
    return(power)
  }
  if (ncore > 1){
    CL1 = parallel::makeCluster(ncore)
    parallel::clusterExport(CL1,c('c1', 'a1', 'c2', 'b2', 'b1', 'cp', 'd1',
                                  'sigx2', 'sigw2', 'sige12', 'sige22', 'sigx_w',
                                  'n', 'nrep', 'alpha','b','nb','pop.cov',
                                  'mu', 'method', 'w_value'),envir = environment())
    
    allsim <- parallel::parLapply(CL1, 1:nrep, runonce)
    parallel::clusterExport(CL1, 'allsim', envir = environment())
    allsim1 = t(parallel::parSapply(CL1, 1:nrep, function(i) unlist(allsim[[i]])))
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
                 method = "moderated mediation model 58",
                 url = "https://webpower.psychstat.org/models/modmed58/",
                 note="power1 is the power of the conditional indirect effect of x on y through m.
power2 is  the power of the direct effect of x on y.
power3 is the power of moderation on the path x to m.
power4 is the power of moderation on the path m to y."), class = "webpower")
  return(power.structure)
  

}

# test = wp.modmed.m58(c1 = 0.2, a1 = 0.2, c2 = 0.1, b2 = 0.1,
#      b1 = 0.2, cp = 0.2, d1 = 0.2, w_value = 0.3, simulation_method = "MC",
#      sigx2 = 1, sigw2 = 1, sige12 = 1, sige22 = 1, sigx_w = 0.5,
#      n = 50, nrep = 1000, alpha = 0.05, ncore = 1)
# print(test)
