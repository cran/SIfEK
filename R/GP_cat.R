#' Estimation of single catalytic constant using Gaussian processes
#'
#'The function estimates single catalytic constant using single data set
#'with an initial enzyme concentrations and substrate concentration.
#'The Gaussian processes is utilized for the likelihood function.
#' @param method method selection: T=TQ model, F=SQ model(default = T)
#' @param RAM Robust Adaptive MCMC options (default = F)
#' @param time total time of data
#' @param dat observed dataset (trajectory column)
#' @param enz enzyme concentrate
#' @param subs substrate concentrate
#' @param MM true value of MM constant
#' @param catal initial value of catalytic constant
#' @param nrepeat total number of iteration (default=10000)
#' @param jump length of distance (default =1)
#' @param burning length of burning period ( default =0)
#' @param catal_m_v Catalytic prior gamma mean, variance(default=c(1,10000))
#' @param sig standard deviation of univariate Normal proposal distribution
#' @param va variance of dataset
#' @return A vector of posterior samples of catalytic constant
#' @details The function GP.cat generates a set of MCMC simulation samples
#' from the conditional posterior distribution of catalytic constant of
#' enzyme kinetics model. As the catalytic constant is only parameter to
#' be estimated in the function the user should assign MM constant as well as
#' initial enzyme concentration and substrate concentration. The prior
#' information for the parameter can be given. The GP.cat function can
#' select Robust Adaptive Metropolis (RAM) algorithm as well as
#' Metropolis-Hastings algorithm with random walk chain for MCMC procedure.
#' When “RAM” is assigned T then the function use RAM method and the “sig”
#' is used as initial standard deviation of normal proposal distribution.
#' When “RAM” is F, the function use Metropolis-Hastings algorithm with
#' random walk chain and the “sig” can be set to controlled proper mixing
#' and acceptance ratio of the parameter from the conditional posterior
#' distribution. The “va” is the variance of the Gaussian process. The
#' posterior samples are only stored with fixed interval according to
#' set "jump" to reduce serial correlation. The initial iterations are
#' removed for convergence. The “burning” is set the length of initial
#' iterations. The Gaussian process method is used for construction of
#' the likelihood.
#' @importFrom stats dgamma rgamma rnorm pnorm runif acf sd dnorm
#' @importFrom graphics par plot lines
#' @importFrom utils write.csv
#' @importFrom MASS mvrnorm
#' @importFrom smfsb simpleEuler
#'
#' @examples
#' \dontrun{
#' data('Chymo_low')
#' time1=max(Chymo_low[,1])*1.01
#' sk_GPMH=GP.cat(method=TRUE,time=time1,dat=Chymo_low[,2],enz=4.4e+7,subs=4.4e+7
#'                       ,MM=4.4e+8,catal=0.05,nrepeat=10000,jump=1,burning=0
#'                       ,catal_m_v=c(1,10000),sig=0.016,va=var(Chymo_low[,2]))
#' # use RAM algorithm #
#' sk_GPRAM=GP.cat(method=TRUE,RAM=TRUE,time=time1,dat=Chymo_low[,2],enz=4.4e+7,subs=4.4e+7
#'                       ,MM=4.4e+8,catal=0.05,nrepeat=10000,jump=1,burning=0
#'                       ,catal_m_v=c(1,10000),sig=0.1,va=var(Chymo_low[,2]))
#'}
#'
#' @export GP.cat
#'
GP.cat = function(method = T, RAM = F, time, dat, enz, subs, MM, catal,
                  nrepeat = 10000, jump = 1, burning = 0, catal_m_v = c(1, 10000), sig,
                  va) {
  catal_m = catal_m_v[1]
  catal_v = catal_m_v[2]

  b_catal = catal_m/catal_v
  a_catal = catal_m * b_catal

  # MM ODE From sQ model#
  MM_sQ <- function(x, t, k = c(k_m = 0, k_p = 0, ET = 0)) {
    with(as.list(c(x, k)), {
      c(-k_p * (ET * ST)/(k_m + ST), k_p * (ET * ST)/(k_m + ST))
    })
  }

  # MM ODE From tQ model #
  MM_tQ <- function(x, t, k = c(k_m = 0, k_p = 0, ET = 0)) {
    with(as.list(c(x, k)), {
      c(-k_p * (1/2) * ((ET + k_m + ST) - sqrt(((ET + k_m + ST)^2) -
                                                 4 * ET * ST)), k_p * (1/2) * ((ET + k_m + ST) - sqrt(((ET +
                                                                                                          k_m + ST)^2) - 4 * ET * ST)))
    })
  }
  if (method == T) {
    M_st = MM_tQ
  } else {
    M_st = MM_sQ
  }

  nrepeat = nrepeat + burning
  x = rep(0, nrepeat * jump)
  Y = rep(0, nrepeat)
  x[1] = catal
  Y[1] = catal

  count = 0
  S = 1
  for (i in 2:(nrepeat * jump)) {
    while (1) {
      u = rnorm(1, 0, sig)
      x_s = (x[i - 1]) + (u * S)
      if (x_s >= 0)
        break
    }

    SE = as.data.frame(smfsb::simpleEuler(t = time, fun = M_st, k = c(k_m = MM,
                                                                      k_p = x_s, ET = enz), ic = c(ST = subs, P = 0), dt = time/101))

    SE2 = as.data.frame(smfsb::simpleEuler(t = time, fun = M_st, k = c(k_m = MM,
                                                                       k_p = x[i - 1], ET = enz), ic = c(ST = subs, P = 0), dt = time/101))

    posterior = (log(dgamma(x_s, a_catal, b_catal)) - sum((dat - SE[,
                                                                    2])^2)/(2 * va) - log(dgamma(x[i - 1], a_catal, b_catal)) +
                   sum((dat - SE2[, 2])^2)/(2 * va))


    accept = min(1, exp(posterior))
    U = runif(1)
    if (U < accept) {
      x[i] = x_s
      count = count + 1
    } else {
      x[i] = x[i - 1]
    }
    if (i%%jump == 0) {
      rep1 = i/jump
      Y[rep1] = x[i]
    }
    {
      if (RAM == T) {
        S = ramcmc::adapt_S(S, u, accept, i, target = 0.44, gamma = min(1,
                                                                        (i)^(-2/3)))
      } else {
        S = 1
      }
    }

  }
  theta = Y[(burning + 1):nrepeat]
  {
    if (method == T) {
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
  write.csv(theta, "catalytic_GPMH.csv")
  return(theta)
}
