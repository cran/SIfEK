% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DA_combi.R
\name{DA.combi}
\alias{DA.combi}
\title{Simultaneous estimation of Michaelis-Menten constant and catalytic constant using
combined data and the likelihood function with the diffusion approximation method.}
\usage{
DA.combi(method = T, dat1, dat2, enz, subs, MM, catal, sig,
  nrepeat = 11000, jump = 1, burning = 1000, catal_m_v = c(1,
  10000), MM_m_v = c(1, 10000), scale_tun = 80)
}
\arguments{
\item{method}{method selection: T=TQ model, F=SQ model(default = T)}

\item{dat1}{observed dataset1 ( time & trajectory columns)}

\item{dat2}{observed dataset2 ( time & trajectory columns)}

\item{enz}{enzyme concentrate}

\item{subs}{substrate concentrate}

\item{MM}{initial value of MM constant}

\item{catal}{initial value of catalytic constant}

\item{sig}{variance of bivariate Normal proposal distribution}

\item{nrepeat}{total number of iteration (default=10000)}

\item{jump}{length of distance (default =1)}

\item{burning}{length of burning period (default =0)}

\item{catal_m_v}{catalytic prior gamma mean, variance(default=c(1,10000))}

\item{MM_m_v}{MM prior gamma mean, variance(default=c(1,10000))}

\item{scale_tun}{scale tunning constant for stochastic simulation}
}
\value{
A n*2 matrix of postrior samples of catalytic constant and MM constant
}
\description{
The function estimates both catalytic constant and Michaelis-Menten constant
simultaneously using combined data sets with different enzyme concentrations
or substrate concentrations. The diffusion approximation is utilized for
the likelihood function.
}
\details{
The function DA.combi generates a set of MCMC simulation samples from
the posterior distribution of catalytic constant and MM constant of enzyme
kinetics model. As the function uses combined data set with different initial
concentration of enzyme or substrate concentration the user should input two
values of enzyme and substrate initial concentration. The prior information
for both two parameters can be given. The turning constant (scale_tun) and
variances for two constants (sig) can be set to controlled proper mixing and
acceptance ratio for updating two parameters simultaneously. The posterior
samples are only stored with fixed interval according to set "jump" to reduce
serial correlation. The initial iterations are removed for convergence.
The “burning” is set the length of initial iterations. The diffusion
approximation method is used for construction of the likelihood.
}
\examples{
\dontrun{
data('Chymo_low')
data('Chymo_high')
comb_DA=DA.combi(method=TRUE,dat1=Chymo_low,dat2=Chymo_high,enz=c(4.4e+7,4.4e+9)
               ,subs=c(4.4e+7,4.4e+7),MM=4.4e+8,catal=0.05,sig=2.0*(c(0.005,8e+7))^2
               ,nrepeat=10000,jump=1,burning=0,catal_m_v = c(1, 1e+6)
               ,MM_m_v = c(1, 1e+10),scale_tun=100)
}

}
