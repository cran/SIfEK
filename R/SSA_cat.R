#' Estimation of single catalytic constant using the stochastic simulation approximation method
#'
#' The function estimates single catalytic constant using single data set with an initial
#' enzyme concentrations and substrate concentration. The stochastic simulation approximation
#' method is utilized for the likelihood function.
#' @param method method selection: T=TQ model, F=SQ model(default = T)
#' @param time observed time interval
#' @param species observed trajectory of product
#' @param enz enzyme concentration
#' @param subs substrate concentration
#' @param MM true value of MM constant
#' @param catal initial value of catalytic constant
#' @param nrepeat total number of iteration (default=10000)
#' @param jump length of distance (default =1)
#' @param burning lenth of burning period (default =0)
#' @param catal_m prior mean of gamma prior (default =1)
#' @param catal_v prior variance of gamma prior (default =1e+6)
#' @return A vector of posterior samples of catalytic constant
#' @details The function SSA.cat generates a set of Monte Carlo simulation samples from
#' the conditional posterior distribution of catalytic constant of enzyme kinetics model.
#' As the catalytic constant is only parameter to be estimated in the function the user
#' should assign MM constant as well as initial enzyme concentration and substrate
#' concentration. The prior information for the parameter can be given.
#' @examples
#' data("Chymo_low")
#' time1=Chymo_low[,1]
#' species1=Chymo_low[,2]
#' Chymotrypsin.low<-SSA.cat(method=TRUE,time=time1,species=species1,enz=4.4e+7,subs=4.4e+7
#'                  ,MM=4.4e+8, catal=0.1,nrepeat = 10000)
#' @export
SSA.cat <- function(method = T, time, species, enz, subs, MM, catal
                  , nrepeat = 10000, jump = 1, burning = 0, catal_m = 1, catal_v = 1e+06) {
    b_k = catal_m/catal_v
    a_k = catal_m * b_k
    nrepeat = nrepeat + burning
    scale = (length(time) - 1)/subs
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

    ################################################################ MCMC ITERATION START
    MCMC_gen = matrix(rep(0, nrepeat), nrow = nrepeat, ncol = 1)
    init = catal
    MCMC_gen[1, ] = init
    catal = init[1]
    MM = MM * scale

    sum_x_t = 0
    n = length(t_diff)
    if (method == T) {
        for (i in 1:n) {
            sum_x_t = sum_x_t + ((enz + S_i[i] + MM)/2 - sqrt((enz + S_i[i] + MM)^2 - 4 * enz * S_i[i])/2) * t_diff[i]
        }
    } else {
        for (i in 1:n) {
            sum_x_t = sum_x_t + (enz * S_i[i]/(MM + S_i[i])) * t_diff[i]
        }
    }
    alpha = sum(n_i) + a_k
    beta = sum_x_t + b_k
    for (rep in 2:nrepeat * jump) {
        catal = rgamma(1, shape = alpha, rate = beta)
        if (rep%%jump == 0) {
            rep1 = rep/jump
            MCMC_gen[rep1, 1] = catal
        }
    }
    theta = MCMC_gen[(burning + 1):nrepeat, ]

    if (method == T) {
        main = "TQ model: Catalytic constant"
    } else {
        main = "SQ model: Catalytic constant"
    }
    par(oma = c(0, 0, 4, 0), mar = c(4, 4, 1, 1))
    mat = matrix(c(1, 1, 2, 3), 2, 2, byrow = T)
    layout(mat)
    plot(theta, type = "l", main = "", xlab = "iteration", ylab = "Catalytic constant")
    acf(theta, main = "")
    plot(density(theta), main = "", xlab = "Catalytic constant")
    mtext(side = 3, line = 1, outer = T, text = main, cex = 1.5)

    cat("MCMC simulation summary", "\n")
    cat("Posterior mean:       ", format(mean(theta), digits = 4, justify = "right", scientific = TRUE), "\n")
    cat("Posterior sd:         ", format(sd(theta), digits = 4, justify = "right", scientific = TRUE), "\n")
    cat("Credible interval(l): ", format(quantile(theta, probs = 0.025), digits = 4, justify = "right", scientific = TRUE),
        "\n")
    cat("Credible interval(U): ", format(quantile(theta, probs = 0.975), digits = 4, justify = "right", scientific = TRUE),
        "\n")
    cat("Relative CV:          ", format((sd(theta)/mean(theta))/(sqrt(catal_v)/catal_m), digits = 4, justify = "right",
        scientific = TRUE), "\n")

    theta = as.data.frame(theta)
    names(theta) = c("Catalytic")
    write.csv(theta, "catalytic.csv")
    return(theta)
}

