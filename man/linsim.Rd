\name{linsim}
\alias{linsim}
\title{Similation of a Self-exciting Point Process}
\description{
  Perform simulation of a self-exciting point process whose intensity also includes a component triggered by another given point process data
and a non-stationary Poisson trend.
}
\usage{
  linsim(data, interval, c, d, ax, ay, at, ptmax)
}

\arguments{
  \item{data}{point process data.}
  \item{interval}{length of time interval in which events take place.}
  \item{c}{exponential coefficient of lgp corresponding to simulated data.}
  \item{d}{exponential coefficient of lgp corresponding to input data.}
  \item{ax}{lgp coefficients in self-exciting part.}
  \item{ay}{lgp coefficients in the input part.}
  \item{at}{coefficients of the polynomial trend.}
  \item{ptmax}{an upper bound of trend polynomial.}
} 

\value{
  \item{in.data}{input data for \code{sim.data}.}
  \item{sim.data}{self-exciting simulated data.}
}

\details{
  This function performs simulation of a self-exciting point process whose intensity also includes 
  a component triggered by another given point process data and non-stationary Poisson trend. 
  The trend is given by usual polynomial, and the response functions to the self-exciting and 
  the external inputs are given the Laguerre-type polynomials (lgp), where the scaling parameters 
  in the exponential functions, say \eqn{c} and \eqn{d}, can be different.
}
\references{
  Y.Ogata, K.Katsura and J.Zhuang (2006) \emph{Computer Science Monographs, No.32,
  TIMSAC84: STATISTICAL ANALYSIS OF SERIES OF EVENTS (TIMSAC84-SASE) VERSION 2}.
  The Institute of Statistical Mathematics.

  Y.Ogata (1981) \emph{On lewis' simulation method for point processes}.
  IEEE information theory, vol. it-27, pp. 23-31.

  Y.Ogata and H.Akaike (1982) \emph{On linear intensity models for mixed doubly stochastic poisson and self-exciting point processes}.
  J. royal statist. soc. b, vol. 44, pp. 102-107.

  Y.Ogata, H.Akaike and K.Katsura (1982) \emph{The application of linear intensity models to the investigatio of causal relations between a point process and another stochastic process}.
  Ann. inst. statist math., vol. 34. pp. 373-387.
}

\examples{
  data(PProcess)   ## The point process data
  linsim( PProcess, 20000, 0.13, 0.026,
          c(0.035,-0.0048), c(0.0,0.00017), c(0.007,-0.00000029), 0.007 )
}

\keyword{spatial}
