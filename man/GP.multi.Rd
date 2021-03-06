% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GP_multi.R
\name{GP.multi}
\alias{GP.multi}
\title{Simultaneous estimation of Michaelis-Menten constant and catalytic constant
using the likelihood function with the Gaussian process method.}
\usage{
GP.multi(method = T, RAM = F, time, dat, enz, subs, MM, catal,
  nrepeat = 10000, jump = 1, burning = 0, catal_m_v = c(1, 10000),
  MM_m_v = c(1, 10000), sig, va)
}
\arguments{
\item{method}{method selection: T=TQ model, F=SQ model(default = T)}

\item{RAM}{Robust Adaptive MCMC options (default = F)}

\item{time}{total time of data}

\item{dat}{observed dataset (trajectory column)}

\item{enz}{enzyme concentrate}

\item{subs}{substrate concentrate}

\item{MM}{initial value of MM constant}

\item{catal}{initial value of cataldattic constant}

\item{nrepeat}{total number of iteration (default=10000)}

\item{jump}{length of distance (default =1)}

\item{burning}{length of burning period (default =0)}

\item{catal_m_v}{MM prior gamma mean, variance(default=c(1,10000))}

\item{MM_m_v}{MM prior gamma mean, variance(default=c(1,10000))}

\item{sig}{variance of bivariate Normal proposal distribution}

\item{va}{variance of dataset}
}
\value{
A vector of posterior samples of catalytic constant
}
\description{
The function estimates both catalytic constant and Michaelis-Menten constant
simultaneously using single data set with an initial enzyme concentrations and
substrate concentration. the Gaussian process is utilized for the likelihood function.
}
\details{
The function GP.multi generates a set of MCMC simulation samples from the posterior
distribution of catalytic constant and MM constant of enzyme kinetics model. As the function
estimates both two constants the user should input the enzyme and substrate initial
concentration. The prior information for both two parameters can be given. The function
can select Robust Adaptive Metropolis (RAM) algorithm as well as Metropolis-Hastings
algorithm with random walk chain for MCMC procedure. When “RAM” is assigned T then the
function use RAM method and the “sig” is used as initial variances of normal proposal
distribution for catalytic and MM constant. When “RAM” is F, the function use
Metropolis-Hastings algorithm with random walk chain and the “sig” can be set to controlled
proper mixing and acceptance ratio of the parameter for updating two parameters simultaneously.
The “va” is the variance of the Gaussian process. The posterior samples are only stored with
fixed interval according to set "jump" to reduce serial correlation. The initial iterations
are removed for convergence. The “burning” is set the length of initial iterations.
The diffusion approximation method is used for construction of the likelihood
}
\examples{
\dontrun{
data('Chymo_low')
time1=max(Chymo_low[,1])*1.01
dou_GPMH=GP.multi(method=TRUE,time=time1,dat=Chymo_low[,2],enz=4.4e+7,subs=4.4e+7
                    ,MM=4.4e+8,catal=0.05,nrepeat=10000,jump=1,burning=0,catal_m_v=c(1,1e+10)
                    ,MM_m_v=c(1e+9,1e+18),sig=c(0.05,4.4e+8)^2,va=var(Chymo_low[,2]))
# use RAM algorithm #
dou_GPRAM=GP.multi(method=TRUE,RAM=TRUE,time=time1,dat=Chymo_low[,2],enz=4.4e+7,subs=4.4e+7
                      ,MM=4.4e+8,catal=0.05,nrepeat=10000,jump=1,burning=0,catal_m_v=c(1,1e+10)
                      ,MM_m_v=c(1e+9,1e+18),sig=c(1,1),va=var(Chymo_low[,2]))
}

}
