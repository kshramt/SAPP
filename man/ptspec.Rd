\name{ptspec}
\alias{ptspec}
\title{The Periodogram of Point Process Data}
\description{
  Provide the periodogram of point process data with the significant band (0.90, 0.95 and 0.99) of the maximum power
 in searching a cyclic component, for stationary Poisson Process.
}
\usage{
  ptspec( data, nfre, prdmin, prd, nsmooth=1, pprd, interval, plot=TRUE )
}
\arguments{
  \item{data}{data of events.}
  \item{nfre}{number of sampling frequencies of spectra.}
  \item{prdmin}{the minimum periodicity of the sampling.}
  \item{prd}{a periodicity for calculating the Rayleigh probability.}
  \item{nsmooth}{number for smoothing of periodgram.}
  \item{pprd}{particular periodicities to be investigated among others.}
  \item{interval}{length of observed time interval of events.}
  \item{plot}{logical. If TRUE (default) the periodogram is plotted.}
}
\value{
  \item{f}{frequency.}
  \item{db}{D.B.}
  \item{power}{power.}
  \item{rayleigh.prob}{the probability of Rayleigh.}
  \item{distance}{= \eqn{\sqrt(rwx^2+rwy^2)}{sqrt(rwx^2+rwy^2)}.}
  \item{phase}{phase.}
}
\references{
  Y.Ogata, K.Katsura and J.Zhuang (2006) \emph{Computer Science Monographs, No.32,
  TIMSAC84: STATISTICAL ANALYSIS OF SERIES OF EVENTS (TIMSAC84-SASE) VERSION 2}.
  The Institute of Statistical Mathematics.
}
\examples{
  data(Brastings)  # The Occurrence Times Data of 627 Brastings
  ptspec( Brastings, 1000, 0.5, 1.0,, c(2.0, 1.0, 0.5), 4600 )

  data(PoissonData)  # to see the contrasting difference
  ptspec( PoissonData, 1000, 0.5, 1.0,, c(2.0, 1.0, 0.5), 5000 )

}
\keyword{spatial}
