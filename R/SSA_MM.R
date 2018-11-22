#' Estimation of single Michaelis-Menten constant using the stochastic simulation approximation
#'
#' The function estimates single Michaelis-Menten constant using the likelihood function
#' with the stochastic simulation approximation method.
#' @param method method selection: T=TQ model, F=SQ model(default = T)
#' @param time observed time interval
#' @param species observed trajectory of product
#' @param enz enzyme concentration
#' @param subs substrate concentration
#' @param MM initial value of MM constant
#' @param catal true value of catalytic constant
#' @param tun tunning constant of MH algorithm (default=2.4)
#' @param std standard deviation of proposal distribution (if =0, caclulated by Opt. function)
#' @param nrepeat total number of iteration (default=10000)
#' @param jump length of distance (default =1)
#' @param burning lenth of burning period (default =0)
#' @param MM_m prior mean of gamma prior (default =1)
#' @param MM_v prior variance of gamma prior (default =10000)
#' @return A vector of posterior samples of Michaelis-Menten constant
#' @details The function SSA.MM generates a set of MCMC simulation samples from the
#' conditional posterior distribution of Michaelis-Menten constant of enzyme kinetics
#' model. As the MM constant is only parameter to be estimated in the function the user
#' should assign catalytic constant as well as initial enzyme concentration and substrate
#' concentration. The prior information for the parameter can be given.
#' The turning constant (scale_tun) and standard deviation of proposal normal
#' distribution (sig) can be set to controlled proper mixing and acceptance ratio of
#' the parameter from the conditional posterior distribution. The posterior samples
#' are only stored with fixed interval according to set "jump" to reduce serial correlation.
#' The initial iterations are removed for convergence. The “burning” is set the length of
#' initial iterations. The diffusion approximation method is used for construction of
#' the likelihood.
#' @examples
#' data("Chymo_low")
#' time1=Chymo_low[,1]
#' species1=Chymo_low[,2]
#'  Chymotrypsin.mm<-SSA.MM(method=TRUE,time=time1,species=species1,enz=4.4e+7,subs=4.4e+7
#'  ,MM=10000,catal=0.051,tun=2.4,std=8e+7 ,nrepeat=10000,jump=1
#'  ,burning=0,MM_m=1,MM_v=1e+10)
#' @export
SSA.MM <- function(method = T, time, species, enz, subs, MM, catal, tun = 2.4
                   , std, nrepeat, jump = 1, burning = 0, MM_m = 1, MM_v = 1e+6) {
    nrepeat = nrepeat + burning
    scale = (length(time) - 1)/subs
    #if (MM_m != 1) {
    #    MM_m = MM_m * scale
    #    MM_v = MM_v * (scale^2)
    #}
    b_k = MM_m/MM_v
    a_k = MM_m * b_k
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
        L.posterior <- function(kd, E = enz, S = S_i, t = t_diff, ni = n_i, k1 = catal
                                , a = a_k, b = b_k, sc = scale) {
            l_lik = 0
            n = length(t)
            for (i in 1:n) {
                lambda_i = k1 * ((E + S[i] + kd *sc)/2 - sqrt((E + S[i] + kd *sc)^2 - 4 * E * S[i])/2)
                l_f_i = log(dgamma(t[i], shape = ni[i], rate = lambda_i) + 1e-300)
                l_lik = l_lik + l_f_i
            }
            l_lik = l_lik + log(dgamma(kd, a, b))
            return(l_lik)
        }
    } else {
        L.posterior <- function(kd, E = enz, S = S_i, t = t_diff, ni = n_i, k1 = catal
                                , a = a_k, b = b_k, sc = scale) {
            l_lik = 0
            n = length(t)
            for (i in 1:n) {
                lambda_i = (k1 * E * S[i]/(kd *sc + S[i]))
                l_f_i = log(dgamma(t[i], shape = ni[i], rate = lambda_i) + 1e-300)
                l_lik = l_lik + l_f_i
            }
            l_lik = l_lik + log(dgamma(kd, a, b))
            return(l_lik)
        }
    }


    # Trancated normal proposal
    proposal <- function(kd, sd) {
        ans = rnorm(1, mean = kd, sd)
        while (ans < 0) {
            ans = rnorm(1, mean = kd, sd)
        }
        return(ans)
    }
    # Metropolis-Hastings step
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
    MCMC_gen = matrix(rep(0, nrepeat), nrow = nrepeat, ncol = 1)
    init = MM
    MCMC_gen[1, ] = init
    kd = init
    count.kd = 0
    for (rep in 2:(nrepeat * jump)) {
        if (std > 0) {
            sd = sqrt(tun) * std
        } else {
            hes <- hessian(func = L.posterior, kd, E = enz, S = S_i, t = t_diff, ni = n_i
                           , k1 = catal, a = a_k, b = b_k, sc=scale)
            sd = sqrt(abs(1/hes)) * tun
        }
        kd.star = proposal(kd, sd)
        L.kd.m <- L.posterior(kd, enz, S_i, t_diff, n_i, catal, a_k, b_k, scale)
        L.kd.star = L.posterior(kd.star, enz, S_i, t_diff, n_i, catal, a_k, b_k, scale)
        kd = MH(L.kd.star, L.kd.m, kd.star, kd, sd)
        if (kd == kd.star)
            count.kd = count.kd + 1
        if (rep%%jump == 0) {
            rep1 = rep/jump
            MCMC_gen[rep1, 1] = kd
        }
    }
    theta = MCMC_gen[(burning + 1):nrepeat, ]
    if (method == T) {
        main = "TQ model: Michaelis-Menten constant"
    } else {
        main = "SQ model: Michaelis-Menten constant"
    }
    par(oma = c(0, 0, 4, 0), mar = c(4, 4, 1, 1))
    mat = matrix(c(1, 1, 2, 3), 2, 2, byrow = T)
    layout(mat)
    plot(theta, type = "l", main = "", xlab = "iteration", ylab = "MM constant")
    acf(theta, main = "")
    plot(density(theta), main = "", xlab = "MM constant")
    mtext(side = 3, line = 1, outer = T, text = main, cex = 1.5)

    cat("MCMC simulation summary", "\n")
    cat("Posterior mean:       ", format(mean(theta), digits = 4, justify = "right", scientific = TRUE), "\n")
    cat("Posterior sd:         ", format(sd(theta), digits = 4, justify = "right", scientific = TRUE), "\n")
    cat("Credible interval(l): ", format(quantile(theta, probs = 0.025), digits = 4, justify = "right", scientific = TRUE),
        "\n")
    cat("Credible interval(U): ", format(quantile(theta, probs = 0.975), digits = 4, justify = "right", scientific = TRUE),
        "\n")
    cat("Relative CV:          ", format((sd(theta)/mean(theta))/(sqrt(MM_v)/MM_m), digits = 4, justify = "right",
        scientific = TRUE), "\n")
    cat("Acceptance ratio:     ", format(count.kd/(nrepeat * jump), digits = 4, justify = "right", scientific = TRUE),
        "\n")

    theta = as.data.frame(theta)
    names(theta) = c("MM")
    write.csv(theta, "MM.csv")
    return(theta)
}



