#' Simultaneous estimation of Michaelis-Menten constant and catalytic constant
#' using the likelihood function with the stochastic simulation approximation method
#'
#' The function estimates both catalytic constant and Michaelis-Menten constant simultaneously
#' using single data set with an initial enzyme concentrations and substrate concentration.
#' The stochastic simulation approximation is utilized for the likelihood function
#' @param method method selection: T=TQ model, F=SQ model(default = T)
#' @param time observed time interval
#' @param species observed trajectory of product
#' @param enz enzyme concentration
#' @param subs substrate concentration
#' @param MM true value of MM constant
#' @param catal initial value of catalytic constant
#' @param tun tunning constant of MH algorithm (default=2.4)
#' @param std standard deviation of proposal distribution (if =0, caclulated by Opt. function)
#' @param nrepeat total number of iteration (default=10000)
#' @param jump length of distance (default =1)
#' @param burning lenth of burning period (default =0)
#' @param catal_m prior mean of gamma prior (default =1)
#' @param catal_v prior variance of gamma prior (default =10000)
#' @param MM_m prior mean of gamma prior (default =1)
#' @param MM_v prior variance of gamma prior (default =10000)
#' @return A n*2 matrix of postrior samples of catalytic constant and MM constant
#' @details The function DA.multi generates a set of MCMC simulation samples from
#' the posterior distribution of catalytic constant and MM constant of enzyme
#' kinetics model. As the function estimates both two constants the user should
#' input the enzyme and substrate initial concentration. The prior information
#' for both two parameters can be given. The function utilizes the Gibbs sampler
#' to update two parameters iteratively from conditional posterior distribution.
#' Updating catalytic constant is conducted using conditional gamma distribution.
#' The posterior samples of MM constant are drawn vis Metropolis-Hasting algorithm
#' with random walk chain. The turning constant (scale_tun) and standard deviation
#' of proposal normal distribution (sig) can be set to controlled proper mixing and
#' acceptance ratio of the parameter from the conditional posterior distribution.
#' The posterior samples are only stored with fixed interval according to set "jump"
#' to reduce serial correlation. The initial iterations are removed for convergence.
#' The “burning” is set the length of initial iterations. The stochastic simulation
#' approximation method is used for construction of the likelihood.
#' @importFrom numDeriv hessian
#' @importFrom stats dgamma rgamma rnorm pnorm runif acf density quantile sd rexp
#' @importFrom graphics par layout mtext plot lines
#' @importFrom utils write.csv
#' @examples
#' data("Chymo_low")
#' time1=Chymo_low[,1]
#' species1=Chymo_low[,2]
#' Chymotrypsin.low<-SSA.multi(method=TRUE, time=time1,species=species1,enz=4.4e+7
#' ,subs=4.4e+7,MM=1e+9,catal=0.01,tun=2.4,std=8e+7,nrepeat=10000,jump=1,
#' burning=0,catal_m=1,catal_v=1e+10, MM_m=1e+9,MM_v=1e+18)
#' @export

SSA.multi <- function(method = T, time, species, enz, subs, MM, catal, tun = 2.4, std, nrepeat, jump = 1, burning = 0,
                         catal_m = 1, catal_v = 10000, MM_m = 1, MM_v = 10000) {
  b_catal = catal_m/catal_v
  a_catal = catal_m * b_catal
  nrepeat = nrepeat + burning
  scale = (length(time) - 1)/subs
  b_MM = MM_m/MM_v
  a_MM = MM_m * b_MM
  dat = cbind(time, species * scale)
  enz = enz * scale
  subs = subs * scale

  init = dat[1, ]
  for (i in 2:nrow(dat)) if (dat[i, 2] > max(dat[1:i - 1, 2])) {
    init = rbind(init, dat[i, ])
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

  # log-likelihood function
  if (method == T) {
    L.posterior <- function(kd, E = enz, S = S_i, t = t_diff, ni = n_i, k1 = catal, a = a_MM, b = b_MM, sc=scale) {
      l_lik = 0
      n = length(t)
      for (i in 1:n) {
        lambda_i = k1 * ((E + S[i] + kd * sc)/2 - sqrt((E + S[i] + kd * sc)^2 - 4 * E * S[i])/2)
        l_f_i = log(dgamma(t[i], shape = ni[i], rate = lambda_i) + 1e-300)
        l_lik = l_lik + l_f_i
      }
      l_lik = l_lik + log(dgamma(kd, a, b) + 1e-300)
      return(l_lik)
    }
  } else {
    L.posterior <- function(kd, E = enz, S = S_i, t = t_diff, ni = n_i, k1 = catal, a = a_MM, b = b_MM, sc=scale) {
      l_lik = 0
      n = length(t)
      for (i in 1:n) {
        lambda_i = (k1 * E * S[i]/(kd  * sc + S[i]))
        l_f_i = log(dgamma(t[i], shape = ni[i], rate = lambda_i) + 1e-300)
        l_lik = l_lik + l_f_i
      }
      l_lik = l_lik + log(dgamma(kd, a, b) + 1e-300)
      return(l_lik)
    }
  }

  # random walk chain sampler
  proposal <- function(kd, sd) {
    ans = rnorm(1, mean = kd, sd)
    while (ans < 0) {
      ans = rnorm(1, mean = kd, sd)
    }
    return(ans)
  }

  MH <- function(l.star, l.m, kd.star, kd.m, sd) {
    logMH = l.star - l.m
    logMH = logMH + log(pnorm(kd.m, mean = 0, sd)) - log(pnorm(kd.star, mean = 0, sd))
    if (logMH < 0) {
      MH = exp(logMH)
      uni = runif(1)
      if (uni < MH) {
        kd = kd.star
      } else {
        kd = kd.m
      }
    } else {
      kd = kd.star
    }
    return(kd)
  }



  # Gibbs sampler of recovery parameter:
  Gibbs_catal <- function(E, S, kd, t, ni, a, b, sc) {
    sum_x_t = 0
    n = length(t)
    if (method == T) {
      for (i in 1:n) {
        sum_x_t = sum_x_t + ((E + S[i] + kd * sc)/2 - sqrt((E + S[i] + kd * sc)^2 - 4 * E * S[i])/2) * t[i]
      }
    } else {
      for (i in 1:n) {
        sum_x_t = sum_x_t + (E * S[i]/(kd * sc + S[i])) * t[i]
      }
    }
    alpha = sum(ni) + a
    beta = sum_x_t + b
    k1 = rgamma(1, shape = alpha, rate = beta)
    return(k1)
  }

  ################################################################ MCMC ITERATION START
  MCMC_gen = matrix(rep(0, nrepeat * 2), nrow = nrepeat, ncol = 2)
  kd = MM
  MCMC_gen[1, ] = cbind(catal, kd)
  count.kd = 0
  for (rep in 2:(nrepeat * jump)) {
    if (std > 0) {
      sd = sqrt(tun) * std
    } else {
      hes <- hessian(func = L.posterior, kd, E = enz, S = S_i, t = t_diff, ni = n_i, k1 = catal, a = a_MM,
                     b = b_MM, sc=scale)
      sd = sqrt(tun) * sqrt(abs(1/hes))
    }
    kd.star = proposal(kd, sd)
    L.kd.m = L.posterior(kd, enz, S_i, t_diff, n_i, catal, a_MM, b_MM, scale)
    L.kd.star = L.posterior(kd.star, enz, S_i, t_diff, n_i, catal, a_MM, b_MM, scale)
    kd = MH(L.kd.star, L.kd.m, kd.star, kd, sd)
    if (kd == kd.star)
      count.kd = count.kd + 1
    catal = Gibbs_catal(enz, S_i, kd, t_diff, n_i, a_catal, b_catal, scale)
    if (rep%%jump == 0) {
      rep1 = rep/jump
      MCMC_gen[rep1, 1] = catal
      MCMC_gen[rep1, 2] = kd
    }
  }

  theta = MCMC_gen[(burning + 1):nrepeat, ]
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
  cat("Posterior mean:       ", format(mean(theta1), digits = 4, justify = "right", scientific = TRUE), "\n")
  cat("Posterior sd:         ", format(sd(theta1), digits = 4, justify = "right", scientific = TRUE), "\n")
  cat("Credible interval(l): ", format(quantile(theta1, probs = 0.025), digits = 4, justify = "right", scientific = TRUE),
      "\n")
  cat("Credible interval(U): ", format(quantile(theta1, probs = 0.975), digits = 4, justify = "right", scientific = TRUE),
      "\n")
  cat("Relative CV:          ", format((sd(theta1)/mean(theta1))/(sqrt(catal_v)/catal_m), digits = 4, justify = "right",
                                       scientific = TRUE), "\n\n")

  par(oma = c(0, 0, 4, 0), mar = c(4, 4, 1, 1))
  mat = matrix(c(1, 1, 2, 3), 2, 2, byrow = T)
  layout(mat)
  plot(theta2, type = "l", main = "", xlab = "iteration", ylab = "MM constant")
  acf(theta2, main = "")
  plot(density(theta2), main = "", xlab = "MM constant")
  mtext(side = 3, line = 1, outer = T, text = main2, cex = 1.5)

  cat("MCMC simulation summary of Michaelis-Menten constanat", "\n")
  cat("Posterior mean:       ", format(mean(theta2), digits = 4, justify = "right", scientific = TRUE), "\n")
  cat("Posterior sd:         ", format(sd(theta2), digits = 4, justify = "right", scientific = TRUE), "\n")
  cat("Credible interval(l): ", format(quantile(theta2, probs = 0.025), digits = 4, justify = "right", scientific = TRUE),
      "\n")
  cat("Credible interval(U): ", format(quantile(theta2, probs = 0.975), digits = 4, justify = "right", scientific = TRUE),
      "\n")
  cat("Relative CV:          ", format((sd(theta2)/mean(theta2))/(sqrt(MM_v)/MM_m), digits = 4, justify = "right",
                                       scientific = TRUE), "\n")
  cat("Acceptance ratio:     ", format(count.kd/(nrepeat * jump), digits = 4, justify = "right", scientific = TRUE),
      "\n")

  par(mfrow = c(1, 1))
  plot(log(theta[, 2]), log(theta[, 1]), type = "p", pch = 19, cex = 0.5, main = "Scatter plot", ylab = "Catalytic constant",
       xlab = "MM constant")


  theta = as.data.frame(theta)
  names(theta) = c("Catalytic", "MM")
  write.csv(theta, "Simul.csv")
  return(theta)
}
