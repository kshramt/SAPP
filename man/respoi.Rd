\name{respoi}
\alias{respoi}
\alias{respoi2}
\title{The residual point process of the ETAS model}
\description{
  Compute the residual of modified Omori Poisson process and display the cumulative curve and magnitude v.s. transformed time.
}
\usage{
  respoi(time,mag,param,zts,tstart,zte,threshold=0.0,plot=TRUE)

  respoi2(etas,param,zts,tstart,zte,threshold=0.0,plot=TRUE)
}
\arguments{
  \item{time}{the time measured from the main shock(t=0).}
  \item{mag}{magnitude.}
  \item{etas}{a etas-format dataset on 9 variables

             (no., longitude, latitude, magnitude, time, depth, year, month and days).}
  \item{param}{the four parameters \eqn{B}, \eqn{K}, \eqn{c} and \eqn{p}.}
  \item{zts}{the start of the precursory period.}
  \item{tstart}{the start of the target period.}
  \item{zte}{the end of the target period.}
  \item{threshold}{threshold magnitude.}
  \item{plot}{logical. If TRUE (default) cumulative curve and magnitude v.s. transformed time \eqn{F(t_i)} are plotted.}
}
\details{
  The function \code{respoi} and \code{respoi2} compute the following output for displaying the goodness-of-fit of Omori-Utsu model to the data.
  The cumulative number of earthquakes at time \eqn{t} since \eqn{t_0} is given by the integration of \eqn{f(t)} with respect to the time \eqn{t},

   \deqn{F(t) = B(t-t_0) + K\{ c^{(1-p)}-(t-t_i+c)^{(1-p)} \}/(p-1)}{
         F(t) = B(t-t_0) + K{ c^{(1-p)}-(t-t_i+c)^{(1-p)} }/(p-1) }

  where the summation of \eqn{i} is taken for all data event. 

  \code{respoi2} is equivalent to \code{respoi} except that input and output forms are different. When a etas-format dataset is given, \code{respoi2} returns the dataset with the format as described below.
}
\value{
  \item{trans.time}{transformed time \eqn{F(t_i), i=1,2,...,N.}}
  \item{cnum}{cumulative number of events.}
  \item{resData}{a res-format dataset on 7 variables

                (no., longitude, latitude, magnitude, time, depth and trans.time) }
}
\references{
  Y.Ogata (2006)
  \emph{Computer Science Monographs, No.33, Statistical Analysis of Seismicity - updated version (SASeies2006).}
  The Institute of Statistical Mathematics.
}
\examples{
  data(main2003JUL26)  # The aftershock data of 26th July 2003 earthquake of M6.2 

# output transformed times and cumulative numbers
  x <- main2003JUL26
  respoi(x$time, x$magnitude, c(0,0.96021E+02,0.58563E-01,0.96611E+00),
         0.0, 0.01, 18.68, 2.5)

# output a res-format dataset
  respoi2(main2003JUL26, c(0,0.96021E+02,0.58563E-01,0.96611E+00),
          0.0, 0.01, 18.68, 2.5)
}

\keyword{spatial}
