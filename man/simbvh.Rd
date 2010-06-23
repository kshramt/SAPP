\name{simbvh}
\alias{simbvh}
\title{Simulation of Bi-variate Hawkes' Mutually Exciting Point Processes
}
\description{
  Perform the simulation of bi-variate Hawkes' mutually exciting point processes.
  The response functions are parameterized by the Laguerre-type polynomials.
}
\usage{
  simbvh(interval,axx=NULL,axy=NULL,axz=NULL,ayx=NULL,ayy=NULL,ayz=NULL,
         c,d,c2,d2,ptxmax,ptymax)
}

\arguments{
  \item{interval}{length of time interval in which events take place.}
  \item{axx}{coefficients of Laguerre polynomial (lgp) of the transfer function (= response function) from the data events
x to x (trf; x --> x).}
  \item{axy}{coefficients of lgp (trf; y --> x).}
  \item{ayx}{coefficients of lgp (trf; x --> y).}
  \item{ayy}{coefficients of lgp (trf; y --> y).}
  \item{axz}{coefficients of polynomial for x data.}
  \item{ayz}{coefficients of polynomial for y data.}
  \item{c}{exponential coefficient of lgp corresponding to xx.}
  \item{d}{exponential coefficient of lgp corresponding to xy.}
  \item{c2}{exponential coefficient of lgp corresponding to yx.}
  \item{d2}{exponential coefficient of lgp corresponding to yy.}
  \item{ptxmax}{an upper bound of trend polynomial corresponding to xz.}
  \item{ptymax}{an upper bound of trend polynomial corresponding to yz.}
} 
\value{
  \item{x}{simulated data X.}
  \item{y}{simulated data Y.}
}
\references{
  Y.Ogata, K.Katsura and J.Zhuang (2006) \emph{Computer Science Monographs, No.32,
  TIMSAC84: STATISTICAL ANALYSIS OF SERIES OF EVENTS (TIMSAC84-SASE) VERSION 2}.
  The Institute of Statistical Mathematics.

  Y.Ogata (1981) \emph{On Lewis' simulation method for point processes}.
  IEEE Information Theory, IT-27, pp.23-31.
}
\examples{
  simbvh(interval=20000,
         axx=0.01623,
         axy=0.007306,
         axz=c(0.006187, -0.00000023),
         ayz=c(0.0046786, -0.00000048, 0.2557e-10),
         c=0.4032,d=0.0219,c2=1.0,d2=1.0,
         ptxmax=0.0062,ptymax=0.08)
}

\keyword{spatial}
