#' Simultaneous estimation of Michaelis-Menten constant and catalytic constant using
#' combined data and the likelihood function with the diffusion approximation method.
#'
#'The function estimates both catalytic constant and Michaelis-Menten constant
#'simultaneously using combined data sets with different enzyme concentrations
#'or substrate concentrations. The diffusion approximation is utilized for
#'the likelihood function.
#' @param method method selection: T=TQ model, F=SQ model(default = T)
#' @param dat1 observed dataset1 ( time & trajectory columns)
#' @param dat2 observed dataset2 ( time & trajectory columns)
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
#' @return  A n*2 matrix of postrior samples of catalytic constant and MM constant
#' @details The function DA.combi generates a set of MCMC simulation samples from
#' the posterior distribution of catalytic constant and MM constant of enzyme
#' kinetics model. As the function uses combined data set with different initial
#' concentration of enzyme or substrate concentration the user should input two
#' values of enzyme and substrate initial concentration. The prior information
#' for both two parameters can be given. The turning constant (scale_tun) and
#' variances for two constants (sig) can be set to controlled proper mixing and
#' acceptance ratio for updating two parameters simultaneously. The posterior
#' samples are only stored with fixed interval according to set "jump" to reduce
#' serial correlation. The initial iterations are removed for convergence.
#' The “burning” is set the length of initial iterations. The diffusion
#' approximation method is used for construction of the likelihood.
#' @importFrom stats dgamma rgamma rnorm pnorm runif acf sd dnorm
#' @importFrom graphics par plot lines
#' @importFrom utils write.csv
#' @importFrom MASS mvrnorm
#' @examples
#' \dontrun{
#' data('Chymo_low')
#' data('Chymo_high')
#' comb_DA=DA.combi(method=TRUE,dat1=Chymo_low,dat2=Chymo_high,enz=c(4.4e+7,4.4e+9)
#'                ,subs=c(4.4e+7,4.4e+7),MM=4.4e+8,catal=0.05,sig=2.0*(c(0.005,8e+7))^2
#'                ,nrepeat=10000,jump=1,burning=0,catal_m_v = c(1, 1e+6)
#'                ,MM_m_v = c(1, 1e+10),scale_tun=100)
#'}
#'
#' @export DA.combi
#'

DA.combi = function(method = T, dat1, dat2, enz, subs, MM, catal, sig,
                    nrepeat = 11000, jump = 1, burning = 1000, catal_m_v = c(1, 10000),
                    MM_m_v = c(1, 10000), scale_tun = 80) {
  catal_m = catal_m_v[1]
  catal_v = catal_m_v[2]
  MM_m = MM_m_v[1]
  MM_v = MM_m_v[2]

  b_catal = catal_m/catal_v
  a_catal = catal_m * b_catal
  b_MM = MM_m/MM_v
  a_MM = MM_m * b_MM

  scale = scale_tun/subs
  {
    if ((subs[1] == 80 & enz[1] == 0.2)) {
      scale[1] = scale[1] * 5
    } else if ((subs[2] == 80 & enz[2] == 0.2)) {
      scale[2] = scale[2] * 5
    }
  }

  dat1[, 2] = dat1[, 2] * scale[1]
  dat2[, 2] = dat2[, 2] * scale[2]
  subs = subs * scale
  enz = enz * scale

  init1 = dat1[1, ]
  for (i in 2:nrow(dat1)) if (dat1[i, 2] > max(dat1[1:i - 1, 2])) {
    init1 = rbind(init1, dat1[i, ])
  }
  n1 = nrow(init1) - 1
  t_diff1 = rep(1, n1)
  n_i1 = rep(1, n1)
  S_i1 = rep(1, n1)
  t_diff1[1] = init1[1, 1]
  n_i1[1] = init1[1, 2]
  S_i1[1] = subs[1]
  for (i in 1:n1) {
    t_diff1[i] = init1[i + 1, 1] - init1[i, 1]
    n_i1[i] = init1[i + 1, 2] - init1[i, 2]
    S_i1[i] = subs[1] - init1[i, 2]
  }

  init2 = dat2[1, ]
  for (i in 2:nrow(dat2)) if (dat2[i, 2] > max(dat2[1:i - 1, 2])) {
    init2 = rbind(init2, dat2[i, ])
  }
  n2 = nrow(init2) - 1
  t_diff2 = rep(1, n2)
  n_i2 = rep(1, n2)
  S_i2 = rep(1, n2)
  t_diff2[1] = init2[1, 1]
  n_i2[1] = init2[1, 2]
  S_i2[1] = subs[2]
  for (i in 1:n2) {
    t_diff2[i] = init2[i + 1, 1] - init2[i, 1]
    n_i2[i] = init2[i + 1, 2] - init2[i, 2]
    S_i2[i] = subs[2] - init2[i, 2]
  }

  {
    if (method == T) {
      M_st = "tQ"
    } else {
      M_st = "sQ"
    }
  }

  nrepeat = nrepeat + burning
  x = matrix(rep(0, (nrepeat * jump) * 2), ncol = 2)
  Y = matrix(rep(0, nrepeat * 2), ncol = 2)
  x[1, ] = c(catal, MM)
  Y[1, ] = c(catal, MM)

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
        lambda1 = (x_s[1]/2 * ((enz[1] + x_s[2] * scale[1] + S_i1) -
                                 sqrt(((enz[1] + x_s[2] * scale[1] + S_i1)^2) - (4 * enz[1] *
                                                                                   S_i1))))
        lambda2 = (x[i - 1, 1]/2 * ((enz[1] + x[i - 1, 2] * scale[1] +
                                       S_i1) - sqrt(((enz[1] + x[i - 1, 2] * scale[1] + S_i1)^2) -
                                                      (4 * enz[1] * S_i1))))
        lambda3 = (x_s[1]/2 * ((enz[2] + x_s[2] * scale[2] + S_i2) -
                                 sqrt(((enz[2] + x_s[2] * scale[2] + S_i2)^2) - (4 * enz[2] *
                                                                                   S_i2))))
        lambda4 = (x[i - 1, 1]/2 * ((enz[2] + x[i - 1, 2] * scale[2] +
                                       S_i2) - sqrt(((enz[2] + x[i - 1, 2] * scale[2] + S_i2)^2) -
                                                      (4 * enz[2] * S_i2))))
      } else {
        lambda1 = (x_s[1] * (enz[1] * S_i1)/(x_s[2] * scale[1] +
                                               S_i1))
        lambda2 = (x[i - 1, 1] * (enz[1] * S_i1)/(x[i - 1, 2] *
                                                    scale[1] + S_i1))
        lambda3 = (x_s[1] * (enz[2] * S_i2)/(x_s[2] * scale[2] +
                                               S_i2))
        lambda4 = (x[i - 1, 1] * (enz[2] * S_i2)/(x[i - 1, 2] *
                                                    scale[2] + S_i2))
      }
    }

    loglike1 = 0
    loglike2 = 0
    loglike3 = 0
    loglike4 = 0

    for (j in 1:length(S_i1)) {
      loglike1[j] = dnorm(n_i1[j], lambda1[j] * t_diff1[j], sqrt(lambda1[j] *
                                                                   t_diff1[j]), log = T)
      loglike2[j] = dnorm(n_i1[j], lambda2[j] * t_diff1[j], sqrt(lambda2[j] *
                                                                   t_diff1[j]), log = T)
    }
    for (j in 1:length(S_i2)) {
      loglike3[j] = dnorm(n_i2[j], lambda3[j] * t_diff2[j], sqrt(lambda3[j] *
                                                                   t_diff2[j]), log = T)
      loglike4[j] = dnorm(n_i2[j], lambda4[j] * t_diff2[j], sqrt(lambda4[j] *
                                                                   t_diff2[j]), log = T)
    }

    post = exp(dgamma(x_s[1], a_catal, b_catal, log = T) + dgamma(x_s[2],
                                                                  a_MM, b_MM, log = T) + sum(loglike1) + sum(loglike3) - dgamma(x[i -
                                                                                                                                    1, 1], a_catal, b_catal, log = T) - dgamma(x[i - 1, 2], a_MM,
                                                                                                                                                                               b_MM, log = T) - sum(loglike2) - sum(loglike4))



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
                                               2])

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
  write.csv(theta, "combined_DA.csv")

  return(theta)
}
