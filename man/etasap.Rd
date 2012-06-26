\name{etasap}
\alias{etasap}
\title{Maximum Likelihood Estimates of The ETAS Model}
\description{
  Compute the maximum likelihood estimates of five parameters of ETAS model. This function consists of two (exact and approximated) versions of the calculation algorithm for the maximization of likelihood.
}
\usage{
  etasap(time, mag, threshold=0.0, reference=0.0, parami,
         zts=0.0, tstart, zte, approx=2, tmpfile=NULL, nlmax=1000, plot=TRUE)
}
\arguments{
  \item{time}{the time measured from the main shock(t=0).}
  \item{mag}{magnitude.}
  \item{threshold}{threshold magnitude.}
  \item{reference}{reference magnitude.}
  \item{parami}{initial estimates of five parameters \eqn{\mu}{mu}, \eqn{K}, \eqn{c}, \eqn{\alpha}{alpha} and \eqn{p}.}
  \item{zts}{the start of the precursory period.}
  \item{tstart}{the start of the target period.}
  \item{zte}{the end of the target period.}
  \item{approx}{>0 : the level for approximation version, which is one of the five levels 1, 2, 4, 8 and 16.
 The higher level means faster processing but lower accuracy. \cr
                =0 : the exact version.}
  \item{tmpfile}{write the process of maximum likelihood procedure to \eqn{tmpfile}.\cr
                 If "" print the process to the standard output and if NULL (default) no report.}
  \item{nlmax}{the maximum number of steps in the process of minimizing.} 
  \item{plot}{logical. If TRUE (default) the graph of cumulative number and magnitude of earthquakes
 against the ordinary time is plotted.}
}
\details{
  The ETAS model is a point-process model representing the activity of earthquakes of magnitude \eqn{M_z} and larger occurring in a certain region during a certain interval of time.
 The total number of such earthquakes is denoted by \eqn{N}. The seismic activity includes primary activity of constant
 occurrence rate \eqn{\mu}{mu} in time (Poisson process). Each earthquake ( including aftershock of another earthquake)
 is followed by its aftershock activity, though only aftershocks of magnitude \eqn{M_z} and larger are included in the data.
 The aftershock activity is represented by the Omori-Usu formula in the time domain. The rate of aftershock occurrence
 at time \eqn{t} following the \eqn{i}th earthquake (time: \eqn{t_i}, magnitude: \eqn{M_i}) is given by 

  \deqn{n_i(t) = K exp[\alpha(M_i-M_z)]/(t-t_i+c)^p,}{%
        n_i(t) = K exp[alpha(M_i-M_z)]/(t-t_i+c)^p,}

 for \eqn{ t>t_i } where \eqn{K}, \eqn{\alpha}{alpha}, \eqn{c}, and \eqn{p} are constants, which are common to all aftershock sequences
 in the region. The rate of occurrence of the whole earthquake series at time \eqn{t} becomes

  \deqn{\lambda(t) = \mu + \Sigma_i  n_i(t).}{%
        lambda(t) = mu + Sigma_i  n_i(t).}

  The summation is done for all \eqn{i} satisfying \eqn{t_i < t}. Five parameters \eqn{\mu}{mu}, \eqn{K}, \eqn{c}, \eqn{\alpha}{alpha} and \eqn{p}
 represent characteristics of seismic activity of the region.
}
\value{
  \item{ngmle}{negative max log-likelihood.}
  \item{param}{list of maximum likelihood estimates of five parameters \eqn{\mu}{mu}, \eqn{K}, \eqn{c}, \eqn{\alpha}{alpha} and \eqn{p}.}
  \item{aic2}{AIC/2.}
}
\references{
  Y.Ogata (2006)
  \emph{Computer Science Monographs, No.33, Statistical Analysis of Seismicity - updated version (SASeies2006).}
  The Institute of Statistical Mathematics.
}
\examples{
  data(main2003JUL26)  # The aftershock data of 26th July 2003 earthquake of M6.2 
  x <- main2003JUL26
  etasap(x$time, x$magnitude, 2.5, 6.2, 
         c(0, 0.63348E+02, 0.38209E-01, 0.26423E+01, 0.10169E+01),, 0.01, 18.68)
}

\keyword{spatial}
