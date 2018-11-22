#' Estimation of single Michaelis-Menten constant using the Gaussian process method

#' The function estimates single Michaelis-Menten constant using the likelihood
#' function with the Gaussian process method.
#' @param method method selection: T=TQ model, F=SQ model(default = T)
#' @param RAM Robust Adaptive MCMC options (default = F)
#' @param time total time of data
#' @param dat observed dataset (trajectory column)
#' @param enz enzyme concentrate
#' @param subs substrate concentrate
#' @param MM initial value of MM constant
#' @param catal true value of catalytic constant
#' @param nrepeat total number of iteration (default=10000)
#' @param jump length of distance (default = 1)
#' @param burning length of burning period (default=0)
#' @param MM_m_v MM prior gamma mean, variance(default=c(1,10000))
#' @param sig standard deviation of univariate Normal proposal distribution
#' @param va variance of dataset
#' @return A vector of posterior samples of catalytic constant
#' @details The function GP.MM generates a set of MCMC simulation samples from the
#' conditional posterior distribution of Michaelis-Menten constant of enzyme kinetics
#' model. As the MM constant is only parameter to be estimated in the function the user
#' should assign catalytic constant as well as initial enzyme concentration and substrate
#' concentration. The prior information for the parameter can be given. The GP.MM function
#' can select Robust Adaptive Metropolis (RAM) algorithm as well as Metropolis-Hastings
#' algorithm with random walk chain for MCMC procedure. When “RAM” is assigned T then the
#' function use RAM method and the “sig” is used as initial standard deviation of normal
#' proposal distribution. When “RAM” is F, the function use Metropolis-Hastings algorithm
#' with random walk chain and the “sig” can be set to controlled proper mixing and acceptance
#' ratio of the parameter from the conditional posterior distribution. The “va” is the
#' variance of the Gaussian process. The posterior samples are only stored with fixed
#' interval according to set "jump" to reduce serial correlation. The initial iterations
#' are removed for convergence. The “burning” is set the length of initial iterations.
#' The diffusion approximation method is used for construction of the likelihood.
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
#' sm_GPMH=GP.MM(method=TRUE,time=time1,dat=Chymo_low[,2],enz=4.4e+7,subs=4.4e+7
#'                       ,MM=4.4e+8,catal=0.05,nrepeat=10000,jump=1,burning=0
#'                       ,MM_m_v=c(1,1e+10),sig=8e+7,va=var(Chymo_low[,2]))
#' # use RAM algorithm #
#' sm_GPRAM=GP.MM(method=TRUE,RAM=TRUE,time=time1,dat=Chymo_low[,2],enz=4.4e+7,subs=4.4e+7
#'                         ,MM=4.4e+8,catal=0.05,nrepeat=10000,jump=1,burning=0
#'                         ,MM_m_v=c(1,1e+10),sig=500,va=var(Chymo_low[,2]))
#'}
#'
#' @export GP.MM
#'
#'

GP.MM = function(method = T, RAM = F, time, dat, enz, subs, MM, catal,
                 nrepeat = 10000, jump = 1, burning = 0, MM_m_v = c(1, 10000), sig,
                 va) {
  MM_m = MM_m_v[1]
  MM_v = MM_m_v[2]

  b_MM = MM_m/MM_v
  a_MM = MM_m * b_MM

  # MM ODE from sQ model #
  MM_sQ <- function(x, t, k = c(k_m = 0, k_p = 0, ET = 0)) {
    with(as.list(c(x, k)), {
      c(-k_p * (ET * ST)/(k_m + ST), k_p * (ET * ST)/(k_m + ST))
    })
  }

  # MM ODE from tQ model #
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
  x = rep(0, (nrepeat * jump))
  Y = rep(0, nrepeat)
  x[1] = MM
  Y[1] = MM
  count = 0
  S = 1
  for (i in 2:(nrepeat * jump)) {
    while (1) {
      u = rnorm(1, 0, sig)
      x_s = x[i - 1] + (u * S)
      if (x_s >= 0)
        break
    }

    SE = as.data.frame(simpleEuler(t = time, fun = M_st, k = c(k_m = x_s,
                                                               k_p = catal, ET = enz), ic = c(ST = subs, P = 0), dt = time/101))
    SE2 = as.data.frame(simpleEuler(t = time, fun = M_st, k = c(k_m = x[i -
                                                                          1], k_p = catal, ET = enz), ic = c(ST = subs, P = 0), dt = time/101))

    posterior = (log(dgamma(x_s, a_MM, b_MM)) - sum((dat - SE[, 2])^2)/(2 *
                                                                          va) - log(dgamma(x[i - 1], a_MM, b_MM)) + sum((dat - SE2[,
                                                                                                                                   2])^2)/(2 * va))

    accept = min(1, exp(posterior))
    U = runif(1)
    if (U < accept) {
      x[i] = x_s
      count = count + 1
    } else {
      x[i] = x[i - 1]
    }

    if (RAM == T) {
      S = ramcmc::adapt_S(S, u, accept, i, target = 0.44, gamma = min(1,
                                                                      (i)^(-2/3)))
    } else {
      S = 1
    }
    if (i%%jump == 0) {
      rep1 = i/jump
      Y[rep1] = x[i]
    }
  }

  theta = Y[(burning + 1):nrepeat]

  {
    if (method == T) {
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
  write.csv(theta, "MM_GPMH.csv")
  return(theta)
}
