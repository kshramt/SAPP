\name{etarpp}
\alias{etarpp}
\alias{etarpp2}
\title{Residual Point Process of The ETAS Model}
\description{
  Compute the residual data using the ETAS model with MLEs.
}
\usage{
  etarpp(time, mag, threshold=0.0, reference=0.0, parami,
         zts=0.0, tstart, zte, ztend=NULL, plot=TRUE)

  etarpp2(etas, threshold=0.0, reference=0.0, parami,
          zts=0.0, tstart, zte, ztend=NULL, plot=TRUE)
}
\arguments{
  \item{time}{the time measured from the main shock(t=0).}
  \item{mag}{magnitude.}
  \item{etas}{a etas-format dataset on 9 variables
             (no., longitude, latitude, magnitude, time, depth, year, month and days).}
  \item{threshold}{threshold magnitude.}
  \item{reference}{reference magnitude.}
  \item{parami}{initial estimates of five parameters \eqn{\mu}{mu}, \eqn{K}, \eqn{c}, \eqn{\alpha}{alpha} and \eqn{p}.}
  \item{zts}{the start of the precursory period.}
  \item{tstart}{the start of the target period.}
  \item{zte}{the end of the target period.}
  \item{ztend}{the end of the prediction period. If NULL (default) the last time of available data is set.}
  \item{plot}{logical. If TRUE (default) the graphs of cumulative number and magnitude against the ordinary time and transformed time are plotted.}
}
\details{
  The cumulative number of earthquakes at time \eqn{t} since \eqn{t_0} is given by the integration of \eqn{\lambda(t)}{lambda(t)} ( see \code{\link{etasap}} ) with respect to the time \eqn{t},

   \deqn{\Lambda(t) = \mu(t-t_0) + K \Sigma_i \exp[\alpha(M_i-M_z)]\{c^{(1-p)}-(t-t_i+c)^{(1-p)}\}/(p-1),}{%
         Lambda(t) = mu(t-t_0) + K*Sigma_i exp[alpha(M_i-M_z)]{c^{(1-p)}-(t-t_i+c)^{(1-p)}}/(p-1),}

  where the summation of \eqn{i} is taken for all data event. 
  The output of \code{etarpp2} is given in a res-format dataset which includes the column of \eqn{\{\Lambda(t_i),i=1,2,...,N\}}{{Lambda(t_i),i=1,2,...,N}}.
}
\value{
  \item{trans.time}{transformed time \eqn{\Lambda(t_i), i=1,2,...,N}{Lambda(t_i), i=1,2,...,N}.}
  \item{no.tstart}{data number of the start of the target period.}
  \item{resData}{a res-format dataset on 7 variables
                (no., longitude, latitude, magnitude, time, depth and transformed time).}
}
\references{
  Y.Ogata (2006)
  \emph{Computer Science Monographs, No.33, Statistical Analysis of Seismicity - updated version (SASeies2006).}
  The Institute of Statistical Mathematics.
}
\examples{
  data(main2003JUL26)  # The aftershock data of 26th July 2003 earthquake of M6.2

## output transformed times and cumulative numbers
  x <- main2003JUL26
  etarpp(x$time, x$magnitude, 2.5, 6.2,
         c(0, 0.63348E+02, 0.38209E-01, 0.26423E+01, 0.10169E+01),, 0.01, 7, 18.68)

## output a res-format dataset
  etarpp2(main2003JUL26, 2.5, 6.2,
         c(0, 0.63348E+02, 0.38209E-01, 0.26423E+01, 0.10169E+01),, 0.01, 7, 18.68)

}

\keyword{spatial}
