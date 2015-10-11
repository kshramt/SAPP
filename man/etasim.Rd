\name{etasim}
\alias{etasim1}
\alias{etasim2}
\title{Simulation of Earthquake Dataset Based on the ETAS Model}
\description{
  Produce simulated dataset for given sets of parameters in the point process model used in ETAS.
}
\usage{
etasim1(bvalue, nd, threshold=0.0, reference=0.0, param)

etasim2(etas, tstart, threshold=0.0, reference=0.0, param)
}
\arguments{
  \item{bvalue}{\eqn{b}-value of G-R law if \code{etasim1}.}
  \item{nd}{the number of the simulated events if \code{etasim1}.}
  \item{etas}{a etas-format dataset on 9 variables
             (no., longitude, latitude, magnitude, time, depth, year, month and days).}
  \item{tstart}{the end of precursory period if \code{etasim2}.}
  \item{threshold}{threshold magnitude.}
  \item{reference}{reference magnitude.}
  \item{param}{five parameters \eqn{\mu}{mu}, \eqn{K}, \eqn{c}, \eqn{\alpha}{alpha} and \eqn{p}.}
}
\details{
  There are two versions; either simulating magnitude by Gutenberg-Richter's Law \code{etasim1}
 or using magnitudes from \eqn{etas} dataset \code{etasim2}.
 For \code{etasim1}, \eqn{b}-value of G-R law and number of events to be simulated are provided.
 \code{stasim2} simulates the same number of events that are not less than threshold magnitude in the dataset \eqn{etas},
 and simulation starts after a precursory period depending on the same history of events in \eqn{etas} in the period.
}
\value{
  \code{etasim1} and \code{etasim2} generate a etas-format dataset given values of 'no.', 'magnitude' and 'time'.
}
\references{
  Ogata, Y. (2006)
  \emph{Computer Science Monographs, No.33, Statistical Analysis of Seismicity - updated version (SASeies2006).}
  The Institute of Statistical Mathematics.
}
\examples{
## by Gutenberg-Richter's Law 
etasim1(1.0, 999, 3.5, 3.5, c(0.2e-02, 0.4e-02, 0.3e-02, 0.24e+01, 0.13e+01))

## from a etas-format dataset
data(main2003JUL26)  # The aftershock data of 26th July 2003 earthquake of M6.2
etasim2(main2003JUL26, 0.01, 2.5, 6.2,
        c(0, 0.63348e+02, 0.38209e-01, 0.26423e+01, 0.10169e+01))
}

\keyword{spatial}
