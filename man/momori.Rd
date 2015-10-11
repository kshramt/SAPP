\name{momori}
\alias{momori}
\title{Maximum Likelihood Estimates of Parameters in the Omori-Utsu (Modified Omori) Formula}
\description{
  Compute the maximum likelihood estimates (MLEs) of parameters in the Omori-Utsu (modified Omori) formula representing for the decay of occurrence rate of aftershocks with time.
}
\usage{
momori(data,mag=NULL,threshold=0.0,tstart,tend,parami,tmpfile=NULL, nlmax=1000)
}
\arguments{
  \item{data}{point process data.}
  \item{mag}{magnitude.}
  \item{threshold}{threshold magnitude.}
  \item{tstart}{the start of the target period.}
  \item{tend}{the end of the target period.}
  \item{parami}{the initial estimates of the four parameters \eqn{B}, \eqn{K}, \eqn{c} and \eqn{p}.}
  \item{tmpfile}{write the process of minimizing to \eqn{tmpfile}.\cr
                 If "" print the process to the standard output and if NULL (default) no report.}
  \item{nlmax}{the maximum number of steps in the process of minimizing.} 
}
\details{
  The modified Omori formula represent the delay law of aftershock activity in time.
  In this equation, \eqn{f(t)} represents the rate of aftershock occurrence at time \eqn{t}, where \eqn{t} is the time measured from the origin time of the main shock.
  \eqn{B}, \eqn{K}, \eqn{c} and \eqn{p} are non-negative constants.
  \eqn{B} represents constant-rate background seismicity which may included in the aftershock data.

     \deqn{f(t) = B + K/(t+c)^p}

  In this function the negative log-likelihood function is minimized by the Davidon-Fletcher-Powell algorithm.
  Starting from a given set of initial guess of the parameters \code{parai}, \code{momori()} repeats calculations of function values and its gradients at each step of parameter vector.
  At each cycle of iteration, the linearly searched step (\eqn{lambda}), negative log-likelihood value (\eqn{-LL}),
  and two estimates of square sum of gradients are shown (\eqn{process=1}).

  The cumulative number of earthquakes at time \eqn{t} since \eqn{t_0} is given by the integration of \eqn{f(t)} with respect to the time \eqn{t},

     \deqn{F(t) = B(t-t_0) + K\{c^{1-p}-(t-t_i+c)^{1-p}\} / (p-1)}{%
           F(t) = B(t-t_0) + K{c^{1-p}-(t-t_i+c)^{1-p}} / (p-1)}

  where the summation of \eqn{i} is taken for all data event.
}
\value{
  \item{param}{the final estimates of the four parameters \eqn{B}, \eqn{K}, \eqn{c} and \eqn{p}.}
  \item{ngmle}{negative max likelihood.}
  \item{aic}{AIC = -2\eqn{LL} + 2*(number of variables), and the number=4 in this case.}
  \item{plist}{list of parameters \eqn{t_i}, \eqn{K}, \eqn{c}, \eqn{p} and \eqn{cls}.}
}
\references{
  Ogata, Y. (2006)
  \emph{Computer Science Monographs, No.33, Statistical Analysis of Seismicity - updated version (SASeies2006).}
  The Institute of Statistical Mathematics.
}
\examples{
data(main2003JUL26)  # The aftershock data of 26th July 2003 earthquake of M6.2 
x <- main2003JUL26
momori(x$time, x$magnitude, 2.5, 0.01, 18.68,
       c(0,0.96021E+02,0.58563E-01,0.96611E+00))
}

\keyword{spatial}
