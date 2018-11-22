#' Simultaneous estimation of Michaelis-Menten constant and catalytic constant
#' using the likelihood function with diffusion approximation method
#'
#' The function estimates both catalytic constant and Michaelis-Menten constant
#' simultaneously using single data set with an initial enzyme concentrations
#' and substrate concentration. The diffusion approximation is utilized for the
#' likelihood function.
#' @param method method selection: T=TQ model, F=SQ model(default = T)
#' @param dat observed dataset (time & trajectory columns)
#' @param enz enzyme concentrate
#' @param subs substrate concentrate
#' @param MM initial value of MM constant
#' @param catal initial value of catalytic constant
#' @param nrepeat total number of iteration (default=10000)
#' @param jump length of distance (default =1)
#' @param burning length of burning period (default =0)
#' @param catal_m_v catalytic prior gamma mean, variance(default=c(1,10000))
#' @param MM_m_v MM prior gamma mean, variance(default=c(1,10000))
#' @param sig variance of bivariate Normal proposal distribution
#' @param scale_tun scale tunning constant for stochastic simulation
#' @return A n*2 matrix of postrior samples of catalytic constant and MM constant
#' @details The function DA.multi generates a set of MCMC simulation samples from
#' the posterior distribution of catalytic constant and MM constant of enzyme
#' kinetics model. As the function estimates both two constants the user should
#' input the enzyme and substrate initial concentration. The prior information for
#' both two parameters can be given. The turning constant (scale_tun) and variances
#' for two constants (sig) can be set to controlled proper mixing and acceptance
#' ratio for updating two parameters simultaneously. The posterior samples are only
#' stored with fixed interval according to set "jump" to reduce serial correlation.
#' The initial iterations are removed for convergence. The “burning” is set the
#' length of initial iterations. The diffusion approximation method is used for
#' construction of the likelihood.
#' @importFrom stats dgamma rgamma rnorm pnorm runif acf sd dnorm
#' @importFrom graphics par plot lines
#' @importFrom utils write.csv
#' @importFrom MASS mvrnorm
#'
#' @examples
#' \dontrun{
#' data('Chymo_low')
#' dou_DA=DA.multi(method=TRUE,dat=Chymo_low,enz=4.4e+7,subs=4.4e+7,MM=4.4e+8,catal=0.05
#'                 ,nrepeat=10000,jump=1,burning=1,catal_m_v=c(1,1e+10),MM_m_v=c(1e+9,1e+18)
#'                 ,sig=2.4*0.001*c(0.05,4.4e+8),scale_tun=80)
#'}
#'
#' @export DA.multi
#'

DA.multi = function(method = T, dat, enz, subs, MM, catal, nrepeat = 10000,
                    jump = 1, burning, catal_m_v = c(1, 10000), MM_m_v = c(1, 10000), sig,
                    scale_tun = 80) {
  catal_m = catal_m_v[1]
  catal_v = catal_m_v[2]
  MM_m = MM_m_v[1]
  MM_v = MM_m_v[2]

  b_catal = catal_m/catal_v
  a_catal = catal_m * b_catal
  b_MM = MM_m/MM_v
  a_MM = MM_m * b_MM

  scale = scale_tun/subs
  dat[, 2] = dat[, 2] * scale
  subs = subs * scale
  enz = enz * scale

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

  {
    if (method == 1) {
      M_st = "tQ"
    } else {
      M_st = "sQ"
    }
  }

  nrepeat = nrepeat + burning
  x = matrix(rep(0, (nrepeat * jump) * 2), ncol = 2)
  Y = matrix(rep(0, (nrepeat) * 2), ncol = 2)
  x[1, ] = c(catal, MM)
  x[1, 2] = x[1, 2] * scale
  Y[1,] = x[1,]
  count = 0
  for (i in 2:(nrepeat * jump)) {
    # random proposal gen #
    while (1) {
      u = mvrnorm(1, c(0, 0), diag(c(sig[1], sig[2])))
      x_s = x[i - 1, ] + u
      if (all(x_s >= 0))
        break
    }
    # posterior probability #

    {
      if (M_st == "tQ") {
        lambda1 = (x_s[1]/2 * ((enz + x_s[2] + S_i) - sqrt(((enz +
                                                               x_s[2] + S_i)^2) - (4 * enz * S_i))))
        lambda2 = (x[i - 1, 1]/2 * ((enz + x[i - 1, 2] + S_i) -
                                      sqrt(((enz + x[i - 1, 2] + S_i)^2) - (4 * enz * S_i))))
      } else {
        lambda1 = (x_s[1] * (enz * S_i)/(x_s[2] + S_i))
        lambda2 = (x[i - 1, 1] * (enz * S_i)/(x[i - 1, 2] + S_i))
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
    post = exp(dgamma(x_s[1], a_catal, b_catal, log = T) + dgamma(x_s[2],
                                                                  a_MM, b_MM/scale, log = T) + sum(loglike1) - dgamma(x[i - 1,
                                                                                                                        1], a_catal, b_catal, log = T) - dgamma(x[i - 1, 2], a_MM,
                                                                                                                                                                b_MM/scale, log = T) - sum(loglike2))



    # MH random walk step #

    accept = min(1, post)
    U = runif(1)
    {
      if (U < accept) {
        x[i, ] = x_s
        count = count + 1
      } else {
        x[i, ] = x[i - 1, ]
      }
    }
    # if(i%%1000==0){print(i);print(x[i,])}
    if (i%%jump == 0) {
      rep1 = i/jump
      Y[rep1, ] = x[i, ]
    }
  }

  theta = cbind(Y[(burning + 1):nrepeat, 1], Y[(burning + 1):nrepeat,
                                               2]/scale)
  if (M_st == "tQ") {
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
  write.csv(theta, "Simul_DA.csv")

  return(theta)
}
