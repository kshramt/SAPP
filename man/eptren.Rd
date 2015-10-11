\name{eptren}
\alias{eptren}
\title{Maximum Likelihood Estimates of Intensity Rates}
\description{
  Compute the maximum likelihood estimates of intensity rates of either exponential polynomial or exponential Fourier series of non-stationary Poisson process models.
}
\usage{
eptren(data, mag=NULL, threshold=0.0, nparam, nsub, cycle=0, tmpfile=NULL,
       nlmax=1000, plot=TRUE)
}
\arguments{
  \item{data}{point process data.}
  \item{mag}{magnitude.}
  \item{threshold}{threshold magnitude.}
  \item{nparam}{maximum number of parameters.}
  \item{nsub}{number of subdivisions in either (0,\eqn{t}) or (0,\eqn{cycle}), where \eqn{t} is the length of observed time interval of points.}
  \item{cycle}{periodicity to be investigated days in a Poisson process model. If zero (default) fit an exponential polynomial model.}
  \item{tmpfile}{file name to write the process of minimizing by davidon-fletcher-powell procedure.
                 If "" print the process to the standard output and if NULL (default) no report.}
  \item{nlmax}{the maximum number of steps in the process of minimizing.} 
  \item{plot}{logical. If TRUE (default) intensity rates are plotted.}
}
\details{
  This function computes the maximum likelihood estimates (MLEs) of the coefficients
  \eqn{A_1, A_2,\ldots A_n}{A_1, A_2,..., A_n} is an exponential polynomial

  \deqn{ f(t) = exp(A_1 + A_2t + A_3t^2 + ... ) }

  or \eqn{A_1, A_2, B_2, ..., A_n, B_n} in a Poisson process model with an intensity taking the form of an
  exponential Fourier series

  \deqn{f(t) = exp\{ A_1 + A_2cos(2\pi t/p) + B_2sin(2\pi t/p) + A_3cos(4\pi t/p) + B_3sin(4\pi t/p) +... \}}{%
        f(t) = exp{ A_1 + A_2cos(2\pi t/p) + B_2sin(2\pi t/p) + A_3cos(4\pi t/p) + B_3sin(4\pi t/p) +... }}

  which represents the time varying rate of occurrence (intensity function) of earthquakes in a region.

  These two models belong to the family of non-stationary Poisson process. The optimal order \eqn{n} can be determined by minimize the value of the Akaike Information Criterion (AIC).
}
\value{
  \item{aic}{AIC.}
  \item{param}{parameters.}
  \item{aicmin}{minimum AIC.}
  \item{maice.order}{number of parameters of minimum AIC.}
  \item{time}{time ( \code{cycle}=0 ) or superposed occurrence time ( \code{cycle}>0 ).}
  \item{intensity}{intensity rates.}
}
\references{
  Ogata, Y., Katsura, K. and Zhuang, J. (2006) \emph{Computer Science Monographs, No.32,
  TIMSAC84: STATISTICAL ANALYSIS OF SERIES OF EVENTS (TIMSAC84-SASE) VERSION 2}.
   The Institute of Statistical Mathematics.

  Ogata, Y. (2006)
  \emph{Computer Science Monographs, No.33, Statistical Analysis of Seismicity - updated version (SASeies2006).}
  The Institute of Statistical Mathematics.
}
\examples{
## The Occurrence Times Data of 627 Blastings
data(Brastings)
eptren(Brastings,,, 10, 1000)      # exponential polynomial trend fitting

eptren(Brastings,,, 10, 1000, 1)   # exponential Fourier series fitting

## Poisson Process data
data(PoissonData)
eptren(PoissonData,,, 10, 1000)      # exponential polynomial trend fitting

eptren(PoissonData,,, 10, 1000, 1)   # exponential Fourier series fitting

## The aftershock data of 26th July 2003 earthquake of M6.2
data(main2003JUL26)
x <- main2003JUL26
eptren(x$time, x$magnitude,, 10, 1000)      # exponential polynomial trend fitting

eptren(x$time, x$magnitude,, 10, 1000, 1)   # exponential Fourier series fitting
}

\keyword{spatial}
