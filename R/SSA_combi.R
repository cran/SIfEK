#' Simultaneous estimation of Michaelis-Menten constant and catalytic constant
#' using combined data and the likelihood function with the Stochastic
#' Simulation Approximation method.
#'
#' The function estimates both catalytic constant and Michaelis-Menten constant
#' simultaneously using combined data sets with different enzyme concentrations
#' or substrate concentrations. The diffusion approximation is utilized for the
#' likelihood function.
#' @param method method selection: T=TQ model, F=SQ model(default = T)
#' @param time1 observed time interval for data1
#' @param time2 observed time interval for data2
#' @param species1 observed trajectory of product for data1
#' @param species2 observed trajectory of product for data2
#' @param enz1 enzyme concentration for data1
#' @param enz2 enzyme concentration for data2
#' @param subs1 substrate concentration for data1
#' @param subs2 substrate concentration for data2
#' @param MM initial value of MM constant
#' @param catal initial value of catalytic constant
#' @param tun tunning constant of MH algorithm (default =2.4)
#' @param std standard deviation of proposal distribution
#' @param nrepeat total number of iteration
#' @param jump length of distance (default =1)
#' @param burning lenth of burning period (default =0)
#' @param catal_m prior mean of gamma prior (default =1)
#' @param catal_v prior variance of gamma prior (default =10000)
#' @param MM_m prior mean of gamma prior (default =1)
#' @param MM_v prior variance of gamma prior (default =10000)
#' @return  A n*2 matrix of postrior samples of catalytic constant and MM constant
#' @details The function DA.combi generates a set of MCMC simulation samples from
#' the posterior distribution of catalytic constant and MM constant of enzyme kinetics
#' model. As the function uses combined data set with different initial concentration
#' of enzyme or substrate concentration the user should input two values of enzyme and
#' substrate initial concentration. The prior information for both two parameters can be
#' given. The function utilizes the Gibbs sampler to update two parameters iteratively
#' from conditional posterior distribution. Updating catalytic constant is conducted
#' using conditional gamma distribution. The posterior samples of MM constant are drawn
#' vis Metropolis-Hasting algorithm with random walk chain.
#' The turning constant (scale_tun) and standard deviation of proposal normal
#' distribution (sig) can be set to controlled proper mixing and acceptance ratio
#' of the parameter from the conditional posterior distribution. The posterior samples
#' are only stored with fixed interval according to set "jump" to reduce serial correlation.
#' The initial iterations are removed for convergence. The “burning” is set the length of
#' initial iterations. The stochastic simulation approximation method is used for
#' construction of the likelihood.
#' @examples
#' \dontrun{
#' data("Chymo_low")
#' time1=Chymo_low[,1]
#' species1=Chymo_low[,2]
#' data("Chymo_high")
#' time2=Chymo_high[,1]
#' species2=Chymo_high[,2]
#' enz.Chymotrypsin<-SSA.combi(method=TRUE, time1=time1 ,time2=time2 ,species1=species1
#'                                ,species2=species2,enz1=4.4e+7,enz2=4.4e+9
#'                                ,subs1=4.4e+7,subs2=4.4e+7,MM=1e+9,catal=0.01,
#'                                tun=2.0,std=8e+7,nrepeat=10000,jump=1,burning=0
#'                                ,catal_m=1,catal_v=1e+6, MM_m=1,MM_v=1e+10)
#' }
#' @export

SSA.combi <- function(method = T, time1, time2, species1, species2, enz1
                         , enz2, subs1, subs2, MM, catal, tun = 2.4, std
                         , nrepeat, jump = 1, burning = 0, catal_m = 1
                         , catal_v = 1e+05, MM_m = 1, MM_v = 1e+05) {
    b_catal = catal_m/catal_v
    a_catal = catal_m * b_catal
    nrepeat = nrepeat + burning
    scale1 = (length(time1) - 1)/subs1
    scale2 = (length(time2) - 1)/subs2
    b_MM = MM_m/MM_v
    a_MM = MM_m * b_MM

    # data re-scaling
    dat1 = cbind(time1, species1 * scale1)
    enz1 = enz1 * scale1
    subs1 = subs1 * scale1
    dat2 = cbind(time2, species2 * scale2)
    enz2 = enz2 * scale2
    subs2 = subs2 * scale2

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
    S_i1[1] = subs1
    for (i in 1:n1) {
        t_diff1[i] = init1[i + 1, 1] - init1[i, 1]
        n_i1[i] = init1[i + 1, 2] - init1[i, 2]
        S_i1[i] = subs1 - init1[i, 2]
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
    S_i2[1] = subs2
    for (i in 1:n2) {
        t_diff2[i] = init2[i + 1, 1] - init2[i, 1]
        n_i2[i] = init2[i + 1, 2] - init2[i, 2]
        S_i2[i] = subs2 - init2[i, 2]
    }

    # log posterior:
    if (method == T) {
      L.posterior <- function(kd, E, S, t, ni, k1, sc) {
        l_lik = 0
        n = length(t)
        for (i in 1:n) {
          lambda_i = k1 * ((E + S[i] + kd * sc)/2 - sqrt((E + S[i] + kd * sc)^2 - 4 * E * S[i])/2)
          l_f_i = log(dgamma(t[i], shape = ni[i], rate = lambda_i) + 1e-300)
          l_lik = l_lik + l_f_i
        }
        return(l_lik)
      }
    } else {
      L.posterior <- function(kd, E, S, t, ni, k1, sc) {
        l_lik = 0
        n = length(t)
        for (i in 1:n) {
          lambda_i = (k1 * E * S[i]/(kd * sc + S[i]))
          l_f_i = log(dgamma(t[i], shape = ni[i], rate = lambda_i) + 1e-300)
          l_lik = l_lik + l_f_i
        }
        return(l_lik)
      }
    }


    if (method == T) {
      L.posterior_com <- function(kd, k1, E1, E2, S1, S2, t1, t2, ni1, ni2, a, b, sc1, sc2) {
        l_lik1 = 0
        n1 = length(t1); kd1 <- kd * sc1;
        l_lik2 = 0
        n2 = length(t2); kd2 <- kd * sc2;
        for (i in 1:n1) {
          lambda_i = k1 * ((E1 + S1[i] + kd1)/2 - sqrt((E1 + S1[i] + kd1)^2 - 4 * E1 * S1[i])/2)
          l_f_i = log(dgamma(t1[i], shape = ni1[i], rate = lambda_i) + 1e-300)
          l_lik1 = l_lik1 + l_f_i
        }
        for (i in 1:n2) {
          lambda_i = k1 * ((E2 + S2[i] + kd2)/2 - sqrt((E2 + S2[i] + kd2)^2 - 4 * E2 * S2[i])/2)
          l_f_i = log(dgamma(t2[i], shape = ni2[i], rate = lambda_i) + 1e-300)
          l_lik2 = l_lik2 + l_f_i
        }

        l_lik = l_lik1 + l_lik2 + log(dgamma(kd, a, b)+ 1e-300)
        return(l_lik)
      }
    } else {
      L.posterior_com <- function(kd, k1, E1, E2, S1, S2, t1, t2, ni1, ni2, a, b, sc1, sc2) {
        l_lik1 = 0
        n1 = length(t1); kd1 <- kd * sc1;
        l_lik2 = 0
        n2 = length(t2); kd2 <- kd * sc2;
        for (i in 1:n1) {
          lambda_i = (k1 * E1 * S1[i]/(kd1 + S1[i]))
          l_f_i = log(dgamma(t1[i], shape = ni1[i], rate = lambda_i) + 1e-300)
          l_lik1 = l_lik1 + l_f_i
        }
        for (i in 1:n2) {
          lambda_i = (k1 * E2 * S2[i]/(kd2 + S2[i]))
          l_f_i = log(dgamma(t2[i], shape = ni2[i], rate = lambda_i) + 1e-300)
          l_lik2 = l_lik2 + l_f_i
        }
        l_lik = l_lik1 + l_lik2 + log(dgamma(kd, a, b)+ 1e-300)
        return(l_lik)
      }
    }


    # Gibbs sampler of recovery parameter:
    Gibbs_catal <- function(E1, E2, S1, S2, kd, t1, t2, ni1, ni2, a, b, sc1, sc2) {
      sum_x_t1 = 0
      n1 = length(t1); kd1 <- kd * sc1;
      sum_x_t2 = 0
      n2 = length(t2); kd2 <- kd * sc2;
      if (method == T) {
        for (i in 1:n1) {
          sum_x_t1 = sum_x_t1 + ((E1 + S1[i] + kd1)/2 - sqrt((E1 + S1[i] + kd1)^2 - 4 * E1 * S1[i])/2) * t1[i]
        }
        for (i in 1:n2) {
          sum_x_t2 = sum_x_t2 + ((E2 + S2[i] + kd2)/2 - sqrt((E2 + S2[i] + kd2)^2 - 4 * E2 * S2[i])/2) * t2[i]
        }
      } else {
        for (i in 1:n1) {
          sum_x_t1 = sum_x_t1 + (E1 * S1[i]/(kd1 + S1[i])) * t1[i]
        }
        for (i in 1:n2) {
          sum_x_t2 = sum_x_t2 + (E2 * S2[i]/(kd2 + S2[i])) * t2[i]
        }
      }
      alpha = sum(ni1) + sum(ni2) + a
      beta = sum_x_t1 + sum_x_t2 + b
      k1 = rgamma(1, shape = alpha, rate = beta)
      return(k1)
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

    ################################################################ MCMC ITERATION START
    MCMC_gen = matrix(rep(0, nrepeat * 2), nrow = nrepeat, ncol = 2)
    kd = MM
    MCMC_gen[1, ] = cbind(catal, kd)
    count.kd = 0

    for (rep in 2:(nrepeat * jump)) {
      if (std > 0) {
        sd = sqrt(tun) * std
      } else {
        hes <- hessian(func = L.posterior_com, kd, k1 = catal, E1 = enz1
                       , E2 = enz2, S1 = S_i1, S2 = S_i2, t1 = t_diff1
                       , t2 = t_diff2, ni1 = n_i1, ni2 = n_i2, a = a_MM
                       , b = b_MM, sc1 =scale1, sc2=scale2)
        if(is.nan(hes)==F) sd = sqrt(abs(1/hes)) * sqrt(tun)
      }
      kd.star = proposal(kd, sd)
      L.kd.m1 <- L.posterior(kd, enz1, S_i1, t_diff1, n_i1, catal, scale1)
      L.kd.star1 <- L.posterior(kd.star, enz1, S_i1, t_diff1, n_i1, catal, scale1)
      L.kd.m2 <- L.posterior(kd, enz2, S_i2, t_diff2, n_i2, catal, scale2)
      L.kd.star2 <- L.posterior(kd.star, enz2, S_i2, t_diff2, n_i2, catal, scale2)
      L.kd.m <- L.kd.m1 + L.kd.m2 + log(dgamma(kd, a_MM, b_MM) + 1e-300)
      L.kd.star <- L.kd.star1 + L.kd.star2 + log(dgamma(kd.star, a_MM, b_MM) + 1e-300)
      #cat('iter=',rep, ', lm=',L.kd.m,', lstar=',L.kd.star,', m=',kd,', m*=',kd.star,'\n')
      kd = MH(L.kd.star, L.kd.m, kd.star, kd, sd)
      catal = Gibbs_catal(enz1, enz2, S_i1, S_i2, kd, t_diff1, t_diff2, n_i1, n_i2
                          , a_catal, b_catal, scale1, scale2)
      if (kd == kd.star)
        count.kd = count.kd + 1
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
    write.csv(theta, "combined.csv")
    return(theta)
}
