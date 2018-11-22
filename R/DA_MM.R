#' Estimation of single Michaelis-Menten constant using diffusion approximation
#'
#' The function estimates single Michaelis-Menten constant using the
#' likelihood function with diffusion approximation method.
#' @param method method selection: T=TQ model, F=SQ model(default = T)
#' @param dat observed dataset ( time & trajectory columns)
#' @param enz enzyme concentrate
#' @param subs substrate concentrate
#' @param MM initial value of MM constant
#' @param catal true value of catalytic constant
#' @param nrepeat total number of iteration (default=10000)
#' @param jump length of distance (default = 1)
#' @param burning length of burning period (default=0)
#' @param MM_m_v MM prior gamma mean, variance(default=c(1,10000))
#' @param sig standard deviation of univariate Normal proposal distribution
#' @param scale_tun scale tunning constant for stochastic simulation
#' @return A vector of posterior samples of catalytic constant
#' @details The function DA.MM generates a set of MCMC simulation samples
#' from the conditional posterior distribution of Michaelis-Menten constant
#' of enzyme kinetics model. As the MM constant is only parameter to be
#' estimated in the function the user should assign catalytic constant
#' as well as initial enzyme concentration and substrate concentration.
#' The prior information for the parameter can be given. The turning
#' constant (scale_tun) and standard deviation (sig) can be set to
#' controlled proper mixing and acceptance ratio of the parameter from
#' the conditional posterior distribution. The posterior samples are only
#' stored with fixed interval according to set "jump" to reduce serial
#' correlation. The initial iterations are removed for convergence.
#' The “burning” is set the length of initial iterations. The diffusion
#' approximation method is used for construction of the likelihood.
#' @importFrom stats dgamma rgamma rnorm pnorm runif acf sd dnorm
#' @importFrom graphics par plot lines
#' @importFrom utils write.csv
#' @importFrom numDeriv hessian
#'
#' @examples
#' \dontrun{
#' data('Chymo_low')
#' sm_DA=DA.MM(method=TRUE,dat=Chymo_low,enz=4.4e+7,subs = 4.4e+7,MM=4.4e+8,catal=0.05
#'                   ,nrepeat=10000,jump=1,burning = 0,MM_m_v=c(1,1e+10),sig=500)
#'}
#'
#' @export DA.MM
#'

DA.MM = function(method = T, dat, enz, subs, MM, catal, nrepeat = 10000,
                 jump = 1, burning = 0, MM_m_v = c(1, 10000), sig, scale_tun = 80) {
  # basic setting # prior alpha & beta setting #

  MM_m = MM_m_v[1]
  MM_v = MM_m_v[2]

  b_MM = MM_m/MM_v
  a_MM = MM_m * b_MM

  # scale tunning #

  scale = scale_tun/subs
  dat[, 2] = dat[, 2] * scale
  subs = subs * scale
  enz = enz * scale

  # trajectory scaling #

  init = dat[1, ]
  for (i in 2:nrow(dat)) if (dat[i, 2] <= max(dat[, 2])) {
    init = rbind(init, dat[i, ])
    if (i == which.max(dat[, 2]))
      break
  }
  n = nrow(init) - 1
  t_diff = rep(1, n)
  n_i = rep(1, n)
  S_i = rep(1, n)
  t_diff[1] = init[1, 1]
  n_i[1] = init[1, 2]
  S_i[1] = subs
  for (i in 1:n) {
    t_diff[i] = init[i + 1, 1] - init[i, 1]
    n_i[i] = init[i + 1, 2] - init[i, 2]
    S_i[i] = subs - init[i, 2]
  }

  {
    if (method == T) {
      M_st = "tQ"
    } else {
      M_st = "sQ"
    }
  }

  ############################################ MH step #
  nrepeat = nrepeat + burning
  x = rep(0, (nrepeat * jump))
  Y = rep(0, nrepeat)
  x[1] = MM * scale
  Y[1] = MM * scale
  count = 0
  for (i in 2:(nrepeat * jump)) {
    # random proposal gen #
    while (1) {
      u = rnorm(1, 0, sig)
      x_s = x[i - 1] + u
      if (x_s >= 0)
        break
    }
    # posterior probability #

    {
      if (M_st == "tQ") {
        lambda1 = (catal/2 * ((enz + x_s + S_i) - sqrt(((enz +
                                                           x_s + S_i)^2) - (4 * enz * S_i))))
        lambda2 = (catal/2 * ((enz + x[i - 1] + S_i) - sqrt(((enz +
                                                                x[i - 1] + S_i)^2) - (4 * enz * S_i))))
      } else {
        lambda1 = (catal * (enz * S_i)/(x_s + S_i))
        lambda2 = (catal * (enz * S_i)/(x[i - 1] + S_i))
      }
    }

    loglike1 = 0
    loglike2 = 0

    for (j in 1:length(S_i)) {
      loglike1[j] = log(dnorm(n_i[j], lambda1[j] * t_diff[j], sqrt(lambda1[j] *
                                                                     t_diff[j])))
      loglike2[j] = log(dnorm(n_i[j], lambda2[j] * t_diff[j], sqrt(lambda2[j] *
                                                                     t_diff[j])))
    }
    post = exp(dgamma(x_s, a_MM, b_MM/scale, log = T) + sum(loglike1) +
                 pnorm(x[i - 1], 0, sig, log.p = T) - dgamma(x[i - 1], a_MM,
                                                             b_MM/scale, log = T) - sum(loglike2) - pnorm(x_s, 0, sig, log.p = T))

    # MH random walk step #

    accept = min(1, post)
    U = runif(1)
    if (U < accept) {
      x[i] = x_s
      count = count + 1
    } else {
      x[i] = x[i - 1]
    }
    # if(i%%1000==0){print(i);print(x[i])}
    if (i%%jump == 0) {
      rep1 = i/jump
      Y[rep1] = x[i]
    }
  }
  theta = Y[(burning + 1):nrepeat]/scale

  {
    if (M_st == "tQ") {
      main = "TQ model: Michaelis-Menten constant"
    } else {
      main = "SQ model: Michaelis-Menten constant"
    }
  }
  par(oma = c(0, 0, 4, 0), mar = c(4, 4, 1, 1))
  mat = matrix(c(1, 1, 2, 3), 2, 2, byrow = T)
  layout(mat)
  plot(theta, type = "l", main = "", xlab = "iteration", ylab = "MM constant")
  acf(theta, main = "")
  plot(density(theta), main = "", xlab = "MM constant")
  mtext(side = 3, line = 1, outer = T, text = main, cex = 1.5)

  cat("MCMC simulation summary", "\n")
  cat("Posterior mean:       ", format(mean(theta), digits = 4, justify = "right",
                                       scientific = TRUE), "\n")
  cat("Posterior sd:         ", format(sd(theta), digits = 4, justify = "right",
                                       scientific = TRUE), "\n")
  cat("Credible interval(l): ", format(quantile(theta, probs = 0.025),
                                       digits = 4, justify = "right", scientific = TRUE), "\n")
  cat("Credible interval(U): ", format(quantile(theta, probs = 0.975),
                                       digits = 4, justify = "right", scientific = TRUE), "\n")
  cat("Relative CV:          ", format((sd(theta)/mean(theta))/(sqrt(MM_v)/MM_m),
                                       digits = 4, justify = "right", scientific = TRUE), "\n")
  cat("Acceptance ratio:     ", format(count/(nrepeat * jump), digits = 4,
                                       justify = "right", scientific = TRUE), "\n")

  theta = as.data.frame(theta)
  names(theta) = c("MM")
  write.csv(theta, "MM_DA.csv")
  return(theta)
}
