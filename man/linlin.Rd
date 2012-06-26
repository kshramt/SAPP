\name{linlin}
\alias{linlin}
\title{Maximum Likelihood Estimates of Linear Intensity Models}
\description{
  Perform the maximum likelihood estimates of linear intensity models of self-exciting point process with another point
 process input, cyclic and trend components.
}

\usage{
  linlin(external, self.excit, interval, c, d, ax=NULL, ay=NULL, ac=NULL,
         at=NULL, opt=0, tmpfile=NULL, nlmax=1000)
}

\arguments{
  \item{external}{another point process data.}
  \item{self.excit}{self-exciting data.}
  \item{interval}{length of observed time interval of event.}
  \item{c}{exponential coefficient of lgp in self-exciting part.}
  \item{d}{exponential coefficient of lgp in input part.}
  \item{ax}{coefficients of self-exciting response function.}
  \item{ay}{coefficients of input response function.}
  \item{ac}{coefficients of cycle.}
  \item{at}{coefficients of trend.}
  \item{opt}{0 : minimize the likelihood with fixed exponential coefficient \eqn{c}\cr
             1 : not fixed \eqn{d}.}
  \item{tmpfile}{write the process of minimizing to \eqn{tmpfile}.\cr
                 If "" print the process to the standard output and if NULL (default) no report.}
  \item{nlmax}{the maximum number of steps in the process of minimizing.} 
}

\value{
  \item{c1}{initial estimate of exponential coefficient of lgp in self-exciting part.}
  \item{d1}{initial estimate of exponential coefficient of lgp in input part.}
  \item{ax1}{initial estimates of lgp coefficients in self-exciting part.}
  \item{ay1}{initial estimates of lgp coefficients in the input part.}
  \item{ac1}{initial estimates of coefficients of Fourier series.}
  \item{at1}{initial estimates of coefficients of the polynomial trend.}
  \item{c2}{final estimate of exponential coefficient of lgp in self-exciting part.}
  \item{d2}{final estimate of exponential coefficient of lgp in input part.}
  \item{ax2}{final estimates of lgp coefficients in self-exciting part.}
  \item{ay2}{final estimates of lgp coefficients in the input part.}
  \item{ac2}{final estimates of coefficients of Fourier series.}
  \item{at2}{final estimates of coefficients of the polynomial trend.}
  \item{aic2}{AIC/2.}
  \item{ngmle}{negative max likelihood.}
  \item{rayleigh.prob}{Rayleigh probability.}
  \item{distance}{= \eqn{\sqrt(rwx^2+rwy^2)}{sqrt(rwx^2+rwy^2)}.}
  \item{phase}{phase.}
}

\details{
 The cyclic part is given by the Fourier series, the trend is given by usual polynomial.
The response functions of the self-exciting and the input are given by the Laguerre type polynomials (lgp), 
where the scaling parameters in the exponential function, say \eqn{c} and \eqn{d}, can be different. 
However it is advised to estimate \eqn{c} first without the input component, and then to estimate \eqn{d} with the fixed \eqn{c} 
(this means that the gradient corresponding to the \eqn{c} is set to keep \eqn{0}), which are good initial estimates for 
the \eqn{c} and \eqn{d} of the mixed self-exciting and input model.
 
 Note that estimated intensity sometimes happen to be negative on some part of time interval 
outside the neighborhood of events.  this take place more easily the larger the number of parameters. 
This causes some difficulty in getting the m.l.e., because the negativity of the intensity contributes to 
the seeming increase of the likelihood.

 Note that for the initial estimates of \eqn{ax(1)}, \eqn{ay(1)} and \eqn{at(1)}, some positive value are necessary. 
Especially 0.0 is not suitable.
}

\references{
  Y.Ogata, K.Katsura and J.Zhuang (2006) \emph{Computer Science Monographs, No.32,
  TIMSAC84: STATISTICAL ANALYSIS OF SERIES OF EVENTS (TIMSAC84-SASE) VERSION 2}.
  The Institute of Statistical Mathematics.

  Y.Ogata and H.Akaike (1982) \emph{On linear intensity models for mixed doubly stochastic poisson and self-exciting point processes}.
  J. royal statist. soc. b, vol. 44, pp. 102-107.

  Y.Ogata, H.Akaike and K.Katsura (1982) \emph{The application of linear intensity models to the investigation of causal relations between a point process and another stochastic process}.
  Ann. inst. statist. math., vol. 34. pp. 373-387.
}

\examples{
  data(PProcess)  # point process data 
  data(SelfExcit) # self-exciting point process data
  linlin( PProcess[1:69], SelfExcit, 20000, 0.13, 0.026,
          c(0.035,-0.0048), c(0.0,0.00017),, c(0.007,-.00000029) )
}

\keyword{spatial}
