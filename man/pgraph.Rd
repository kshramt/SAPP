\name{pgraph}
\alias{pgraph}
\title{Graphical Outputs for the Point Process Data Set}
\description{
  Provide the several graphical outputs for the point process data set.
}
\usage{
pgraph(data, mag, threshold=0.0, h, npoint, days, delta=0.0, dmax=0.0,
       separate.graphics=FALSE)
}
\arguments{
  \item{data}{point process data.}
  \item{mag}{magnitude.}
  \item{threshold}{threshold magnitude.}
  \item{h}{time length of the moving interval in which points are counted to show the graph.}
  \item{npoint}{number of subintervals in (0,days) to estimate a non parametric intensity under the palm probability measure.}
  \item{days}{length of interval to display the intensity estimate under the palm probability.}
  \item{delta}{length of a subinterval unit in (0,dmax) to compute the variance time curve.}
  \item{dmax}{time length of a interval to display the variance time curve; \cr
  this is less than (length of whole interval)/4.
    As the default setting of either delta=0.0 or dmax=0.0, set dmax = (length of whole interval)/4 and delta = dmax/100.}
  \item{separate.graphics}{logical. If TRUE a graphic device is opened for each graphics display.}
}
\value{
  \item{cnum}{cumulative numbers of events time.}
  \item{lintv}{interval length.}
  \item{tau}{=time*(total number of events)/(time end).}
  \item{nevent}{number of events in [tau, tau+h].}
  \item{survivor}{log survivor curve with i*(standard error), i=1,2,3.}
  \item{deviation}{deviation of survivor function from the Poisson.}
  \item{nomal.cnum}{normalized cumulative number.}
  \item{nomal.lintv}{U(i)=-exp(-(normalized interval length)).}
  \item{success.intv}{successive pair of intervals.}
  \item{occur}{occurrence rate.}
  \item{time}{time assuming the stationary Poisson process.}
  \item{variance}{Var(N(0,time)).}
  \item{error}{the 0.95 and 0.99 error lines assuming the stationary Poisson process.}
}
\references{
  Ogata, Y., Katsura, K. and Zhuang, J. (2006) \emph{Computer Science Monographs, No.32,
  TIMSAC84: STATISTICAL ANALYSIS OF SERIES OF EVENTS (TIMSAC84-SASE) VERSION 2}.
   The Institute of Statistical Mathematics.

  Ogata, Y. (2006)
  \emph{Computer Science Monographs, No.33, Statistical Analysis of Seismicity - updated version (SASeies2006).}
  The Institute of Statistical Mathematics.

  Ogata, Y. and Shimazaki, K. (1984) \emph{Transition from aftershock to normal activity: The 1965 Rat islands earthquake aftershock sequence}.
  Bulletin of the seismological society of America, vol. 74, no. 5, pp. 1757-1765.
}
\examples{
## The aftershock data of 26th July 2003 earthquake of M6.2
data(main2003JUL26)
x <- main2003JUL26
pgraph(data=x$time, mag=x$magnitude, h=6, npoint=100, days=10)

## The residual point process data of 26th July 2003 earthquake of M6.2
data(res2003JUL26)
y <- res2003JUL26
pgraph(data=y$trans.time, mag=y$magnitude, h=6, npoint=100, days=10)
}

\keyword{spatial}
