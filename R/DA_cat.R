#' Estimation of single catalytic constant using the diffusion approximation
#'
#' The function estimates single catalytic constant using single data set with
#' an initial enzyme concentrations and substrate concentration.
#' The diffusion approximation is utilized for the likelihood function.
#' @param method method selection: T=TQ model, F=SQ model(default = T)
#' @param dat observed dataset ( time & trajectory columns)
#' @param enz enzyme concentration
#' @param subs substrate concentration
#' @param MM true value of MM constant
#' @param catal initial value of catalytic constant
#' @param nrepeat total number of iteration (default=10000)
#' @param jump length of distance (default =1)
#' @param burning lenth of burning period (default =0)
#' @param catal_m_v Catalytic prior gamma mean, variance(default=c(1,10000))
#' @param sig standard deviation of univariate Normal proposal distribution
#' @param scale_tun scale tunning constant for stochastic simulation
#' @return A vector of posterior samples of catalytic constant
#' @details The function DA.cat generates a set of MCMC simulation samples from
#' the conditional posterior distribution of catalytic constant of enzyme kinetics model.
#' As the catalytic constant is only parameter to be estimated in the function the user
#' should assign MM constant as well as initial enzyme concentration and
#' substrate concentration. The prior information for the parameter can be given.
#' The turning constant (scale_tun) and standard deviation (sig) can be set to controlled
#' proper mixing and acceptance ratio of the parameter from the conditional posterior
#' distribution. The posterior samples are only stored with fixed interval according to
#' set "jump" to reduce serial correlation. The initial iterations are removed for
#' convergence. The “burning” is set the length of initial iterations. The diffusion
#' approximation method is used for construction of the likelihood.
#' @importFrom stats dgamma rgamma rnorm pnorm runif acf sd dnorm density quantile
#' @importFrom graphics par layout mtext plot lines
#' @importFrom utils write.csv
#'
#' @examples
#' \dontrun{
#' data('Chymo_low')
#' sk_DA=DA.cat(method=TRUE,dat=Chymo_low,enz=4.4e+7,subs = 4.4e+7,MM=4.4e+8,catal=0.05
#'                   ,nrepeat=10000,jump=1,burning = 0,catal_m_v=c(1,10000),sig=0.005)
#'}
#'
#' @export DA.cat
#'
#'

DA.cat <- function(method = T, dat, enz, subs, MM, catal, nrepeat = 10000,
                   jump = 1, burning = 0, catal_m_v = c(1, 10000), sig, scale_tun = 80) {
  # basic setting # prior alpha & beta setting #

  catal_m = catal_m_v[1]
  catal_v = catal_m_v[2]

  b_catal = catal_m/catal_v
  a_catal = catal_m * b_catal

  # scale tunning #
  scale = scale_tun/subs
  dat[, 2] = dat[, 2] * scale
  subs = subs * scale
  enz = enz * scale
  MM = MM * scale

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
  ###################################################### MH step #
  nrepeat = nrepeat + burning
  x = rep(0, nrepeat * jump)
  Y = rep(0, nrepeat)
  x[1] = catal
  Y[1] = catal
  count = 0
  for (i in 2:(nrepeat * jump)) {
    # proposal gen #
    while (1) {
      u = rnorm(1, 0, sig)
      x_s = x[i - 1] + u
      if (x_s >= 0)
        break
    }

    # posterior probability gen #

    {
      if (M_st == "sQ") {
        lambda1 = (x_s * (enz * S_i)/(MM + S_i))
        lambda2 = (x[i - 1] * (enz * S_i)/(MM + S_i))
      } else (M_st == "tQ")
      {
        lambda1 = (x_s/2 * ((enz + MM + S_i) - sqrt(((enz + MM +
                                                        S_i)^2) - (4 * enz * S_i))))
        lambda2 = (x[i - 1]/2 * ((enz + MM + S_i) - sqrt(((enz +
                                                             MM + S_i)^2) - (4 * enz * S_i))))
      }
    }

    loglike1 = 0
    loglike2 = 0

    for (j in 1:length(S_i)) {
      loglike1[j] = dnorm(n_i[j], lambda1[j] * t_diff[j], sqrt(lambda1[j] *
                                                                 t_diff[j]), log = T)
      loglike2[j] = dnorm(n_i[j], lambda2[j] * t_diff[j], sqrt(lambda2[j] *
                                                                 t_diff[j]), log = T)
    }

    post = exp(dgamma(x_s, a_catal, b_catal, log = T) + sum(loglike1) -
                 dgamma(x[i - 1], a_catal, b_catal, log = T) - sum(loglike2))

    accept = min(1, post)

    # MH random walk step #
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

  # summary step #
  theta = Y[(burning + 1):nrepeat]
  {
    if (M_st == "tQ") {
      main = "TQ model: Catalytic constant"
    } else {
      main = "SQ model: Catalytic constant"
    }
  }
  par(oma = c(0, 0, 4, 0), mar = c(4, 4, 1, 1))
  mat = matrix(c(1, 1, 2, 3), 2, 2, byrow = T)
  layout(mat)
  plot(theta, type = "l", main = "", xlab = "iteration", ylab = "Catalytic constant")
  acf(theta, main = "")
  plot(density(theta), main = "", xlab = "Catalytic constant")
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
  cat("Relative CV:          ", format((sd(theta)/mean(theta))/(sqrt(catal_v)/catal_m),
                                       digits = 4, justify = "right", scientific = TRUE), "\n")
  theta = as.data.frame(theta)
  names(theta) = c("Catalytic")
  write.csv(theta, "catalytic_DA.csv")
  return(theta)
}
