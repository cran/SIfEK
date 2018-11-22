#' Simultaneous estimation of Michaelis-Menten constant and catalytic constant
#' using the likelihood function with the Gaussian process method.
#'
#' The function estimates both catalytic constant and Michaelis-Menten constant
#' simultaneously using single data set with an initial enzyme concentrations and
#' substrate concentration. the Gaussian process is utilized for the likelihood function.
#' @param method method selection: T=TQ model, F=SQ model(default = T)
#' @param RAM Robust Adaptive MCMC options (default = F)
#' @param time total time of data
#' @param dat observed dataset (trajectory column)
#' @param enz enzyme concentrate
#' @param subs substrate concentrate
#' @param MM initial value of MM constant
#' @param catal initial value of cataldattic constant
#' @param nrepeat total number of iteration (default=10000)
#' @param jump length of distance (default =1)
#' @param burning length of burning period (default =0)
#' @param catal_m_v MM prior gamma mean, variance(default=c(1,10000))
#' @param MM_m_v MM prior gamma mean, variance(default=c(1,10000))
#' @param sig variance of bivariate Normal proposal distribution
#' @param va variance of dataset
#' @return A vector of posterior samples of catalytic constant
#' @details The function GP.multi generates a set of MCMC simulation samples from the posterior
#' distribution of catalytic constant and MM constant of enzyme kinetics model. As the function
#' estimates both two constants the user should input the enzyme and substrate initial
#' concentration. The prior information for both two parameters can be given. The function
#' can select Robust Adaptive Metropolis (RAM) algorithm as well as Metropolis-Hastings
#' algorithm with random walk chain for MCMC procedure. When “RAM” is assigned T then the
#' function use RAM method and the “sig” is used as initial variances of normal proposal
#' distribution for catalytic and MM constant. When “RAM” is F, the function use
#' Metropolis-Hastings algorithm with random walk chain and the “sig” can be set to controlled
#' proper mixing and acceptance ratio of the parameter for updating two parameters simultaneously.
#' The “va” is the variance of the Gaussian process. The posterior samples are only stored with
#' fixed interval according to set "jump" to reduce serial correlation. The initial iterations
#' are removed for convergence. The “burning” is set the length of initial iterations.
#' The diffusion approximation method is used for construction of the likelihood
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
#' dou_GPMH=GP.multi(method=TRUE,time=time1,dat=Chymo_low[,2],enz=4.4e+7,subs=4.4e+7
#'                     ,MM=4.4e+8,catal=0.05,nrepeat=10000,jump=1,burning=0,catal_m_v=c(1,1e+10)
#'                     ,MM_m_v=c(1e+9,1e+18),sig=c(0.05,4.4e+8)^2,va=var(Chymo_low[,2]))
#' # use RAM algorithm #
#' dou_GPRAM=GP.multi(method=TRUE,RAM=TRUE,time=time1,dat=Chymo_low[,2],enz=4.4e+7,subs=4.4e+7
#'                       ,MM=4.4e+8,catal=0.05,nrepeat=10000,jump=1,burning=0,catal_m_v=c(1,1e+10)
#'                       ,MM_m_v=c(1e+9,1e+18),sig=c(1,1),va=var(Chymo_low[,2]))
#'}
#'
#' @export GP.multi
#'
GP.multi = function(method = T, RAM = F, time, dat, enz, subs, MM, catal,
                    nrepeat = 10000, jump = 1, burning = 0, catal_m_v = c(1, 10000), MM_m_v = c(1,
                                                                                                10000), sig, va) {
  catal_m = catal_m_v[1]
  catal_v = catal_m_v[2]
  MM_m = MM_m_v[1]
  MM_v = MM_m_v[2]

  b_catal = catal_m/catal_v
  a_catal = catal_m * b_catal
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
  x = matrix(rep(0, 2 * (nrepeat * jump)), ncol = 2)
  Y = matrix(rep(0, 2 * nrepeat), ncol = 2)
  x[1, ] = c(catal, MM)
  Y[1, ] = c(catal, MM)
  if (RAM == T) {
    S = diag(c((2.4)/sqrt(2) * c((catal), (MM))))
  } else {
    S = diag(2)
  }
  count = 0

  for (i in 2:(nrepeat * jump)) {
    while (1) {
      u = as.matrix(mvrnorm(1, c(0, 0), diag(c(sig[1], sig[2]))))
      x_s = (x[i - 1, ] + t(S %*% u))
      if (all(x_s >= 0))
        break
    }

    SE = as.data.frame(simpleEuler(t = time, fun = M_st, k = c(k_m = x_s[2],
                                                               k_p = x_s[1], ET = enz), ic = c(ST = subs, P = 0), dt = time/101))
    SE2 = as.data.frame(simpleEuler(t = time, fun = M_st, k = c(k_m = x[i -
                                                                          1, 2], k_p = x[i - 1, 1], ET = enz), ic = c(ST = subs, P = 0),
                                    dt = time/101))

    posterior = (log(dgamma(x_s[1], a_catal, b_catal)) + log(dgamma(x_s[2],
                                                                    a_MM, b_MM)) - sum((dat - SE[, 2])^2)/(2 * va) - log(dgamma(x[i -
                                                                                                                                    1, 1], a_catal, b_catal)) - log(dgamma(x[i - 1, 2], a_MM, b_MM)) +
                   sum((dat - SE2[, 2])^2)/(2 * va))


    accept = min(1, exp(posterior))

    U = runif(1)
    if (U < accept) {
      x[i, ] = x_s
      count = count + 1
    } else {
      x[i, ] = x[i - 1, ]
    }
    if (RAM == T) {
      S = ramcmc::adapt_S(S, u, accept, i, gamma = min(1, (2 * i)^(-2/3)))
    } else {
      S = diag(2)
    }
    if (i%%jump == 0) {
      rep1 = i/jump
      Y[rep1, ] = x[i, ]
    }
  }
  theta = Y[(burning + 1):nrepeat, ]

  if (method == T) {
    main1 = "TQ model: Catalytic constant"
    main2 = "TQ model: Michaelis-Menten constant"
  } else {
    main1 = "SQ model: Catalytic constant"
    main2 = "SQ model: Michaelis-Menten constant"
  }
  theta1 = theta[, 1]
  theta2 = theta[, 2]

  par(oma = c(0, 0, 4, 0), mar = c(4, 4, 1, 1))
  mat = matrix(c(1, 1, 2, 3), 2, 2, byrow = T)
  layout(mat)
  plot(theta1, type = "l", main = "", xlab = "iteration", ylab = "Catalytic constant")
  acf(theta1, main = "")
  plot(density(theta1), main = "", xlab = "Catalytic constant")
  mtext(side = 3, line = 1, outer = T, text = main1, cex = 1.5)

  cat("MCMC simulation summary of catalytic constanat", "\n")
  cat("Posterior mean:       ", format(mean(theta1), digits = 4, justify = "right",
                                       scientific = TRUE), "\n")
  cat("Posterior sd:         ", format(sd(theta1), digits = 4, justify = "right",
                                       scientific = TRUE), "\n")
  cat("Credible interval(l): ", format(quantile(theta1, probs = 0.025),
                                       digits = 4, justify = "right", scientific = TRUE), "\n")
  cat("Credible interval(U): ", format(quantile(theta1, probs = 0.975),
                                       digits = 4, justify = "right", scientific = TRUE), "\n")
  cat("Relative CV:          ", format((sd(theta1)/mean(theta1))/(sqrt(catal_v)/catal_m),
                                       digits = 4, justify = "right", scientific = TRUE), "\n\n")

  par(oma = c(0, 0, 4, 0), mar = c(4, 4, 1, 1))
  mat = matrix(c(1, 1, 2, 3), 2, 2, byrow = T)
  layout(mat)
  plot(theta2, type = "l", main = "", xlab = "iteration", ylab = "MM constant")
  acf(theta2, main = "")
  plot(density(theta2), main = "", xlab = "MM constant")
  mtext(side = 3, line = 1, outer = T, text = main2, cex = 1.5)

  cat("MCMC simulation summary of Michaelis-Menten constanat", "\n")
  cat("Posterior mean:       ", format(mean(theta2), digits = 4, justify = "right",
                                       scientific = TRUE), "\n")
  cat("Posterior sd:         ", format(sd(theta2), digits = 4, justify = "right",
                                       scientific = TRUE), "\n")
  cat("Credible interval(l): ", format(quantile(theta2, probs = 0.025),
                                       digits = 4, justify = "right", scientific = TRUE), "\n")
  cat("Credible interval(U): ", format(quantile(theta2, probs = 0.975),
                                       digits = 4, justify = "right", scientific = TRUE), "\n")
  cat("Relative CV:          ", format((sd(theta2)/mean(theta2))/(sqrt(MM_v)/MM_m),
                                       digits = 4, justify = "right", scientific = TRUE), "\n")
  cat("Acceptance ratio:     ", format(count/(nrepeat * jump), digits = 4,
                                       justify = "right", scientific = TRUE), "\n")

  par(mfrow = c(1, 1))
  plot(log(theta[, 2]), log(theta[, 1]), type = "p", pch = 19, cex = 0.5,
       main = "Scatter plot", ylab = "Catalytic constant", xlab = "MM constant")


  theta = as.data.frame(theta)
  names(theta) = c("Catalytic", "MM")
  write.csv(theta, "Simul_GPMH.csv")

  return(theta)
}

