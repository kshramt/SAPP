.packageName <- "SAPP"

eptren <- function( data,mag=NULL,threshold=0.0,nparam,nsub,cycle=0,tmp.file=NULL,plot=TRUE )
{

# data                   #  point process data
  n <- length(data)      #  total number of variable
# mag                    #  magnitude
# threshold              #  threshold magnitude
# nparam                 #  maximum number of parameters
# nsub                   #  number of subdivisions in either (0,t) or (0,cycle)
                         #  for the numerical integration of intensity.
# cycle             　　　　      # > 0 : periodicity to be investigated in days
  if( cycle == 0 ) nfunct = 1     #  trend  (exponential polynomial trend fitting)
  if( cycle > 0 )  nfunct = 2     #  cycle  (exponential fourier series fitting)

  nn <- 0
  xx <- NULL
  if( is.null(mag) ) {
    xx <- data
    magnitude <- data
    nn <- n
    t <- data[nn]
    if( nfunct == 2 )
      for( i in 1:nn )  magnitude[i] <- magnitude[i]-as.integer(data[i]/t)*t
  } else {
    magnitude <- NULL
    for( i in 1:n ) 
      if( (mag[i]>=threshold) && (data[i]>=0.0) )  {
        nn <- nn+1
        xx <- c(xx,data[i])
        magnitude <- c(magnitude,mag[i])
      }
    t <- data[nn]
  }

  nmax <- nparam
  if( nfunct==2 ) nmax <- 2*nparam-1
  tmp.file <- paste(tmp.file, " ")
    
  xa <- array(0, dim=c(nmax,nparam))
  aic <- rep(0, nparam)
  aicmin <- 0
  min_order <- 0

  np <- 101
  xval <- rep(0,np)
  fval <- rep(0,np)

  z <- .Fortran("eptrenf",
	as.double(xx),
	as.double(t),
	as.integer(nn),
	as.integer(nfunct),
	as.integer(nparam),
	as.integer(nsub),
	as.double(cycle),
	xa = as.double(xa),
	aic = as.double(aic),
	aicmin = as.double(aicmin),
	min_order = as.integer(min_order),
       xval = as.double(xval),
       fval = as.double(fval),
	as.integer(nmax),
	as.integer(np),
	as.character(tmp.file) )

  xa <- array(z$xa, dim=c(nmax,nparam))
  param <- list()
  for( i in 1:nparam ) {
    m <- i
    if( nfunct == 2 ) m <- m*2-1
    param[[i]] <- xa[1:m,i]
  }
 
  if( plot == TRUE ) {
    par(mfrow=c(2,1))
    data1 <-matrix(, nrow=nn, ncol=2)  
    data1[,1] <- xx
    data1[,2] <- magnitude
    if( is.null(mag) ) data1[1:nn,2] <- rep(1.0, nn)
    data2 <- matrix(, nrow=np, ncol=2)
    data2[,1] <- z$xval
    data2[,2] <- z$fval
    if( cycle == 0 ) {
      plot(data2, type="l", main="EPTREN-Trend", xlab="Time", ylab="Intensity Rate" )
      plot(data1, type="h", main="", xlab="Time", ylab="Magnitude" )
    } else if( cycle > 0 ) {
      for( i in 1:nn ) data1[i,1] <- data1[i,1]-as.integer(data1[i,1]/cycle)*cycle
      plot(data2, type='l', main="EPTREN > EXP(CYCLE)", xlab="Superposed Occurrence Times", ylab="Intensity Rate" )
      plot(data1, type="h", main="", xlab="Superposed Occurrence Times", ylab="Magnitude"  ) 
    }
    par(mfrow=c(1,1))
  }

  eptren.out <- list( aic=z$aic,param=param,aicmin=z$aicmin,maice.order=z$min_order,time=z$xval,intensity=z$fval )
  class( eptren.out ) <- "eptren"
  return( eptren.out )
}

print.eptren <- function(x, ...)
{
# np <- length(x$time)
# cat("Time\t\tIntensity Rate\n")
# for( i in 1:np )
# cat(sprintf("%f\t%f\n", x$time, x$intensity))

  nparam <- length(x$aic)
  for( i in 1:nparam ) {
    cat(sprintf("\n AIC\t%f\n", x$aic[i]))
    cat(" parameters\n")
    ii <- length(x$param[[i]])
    for( j in 1:ii ) {
      cat(sprintf("\t%e",x$param[[i]][j]))
      if( (j%%5 == 0) & j>1 ) cat("\n")
    }
    cat("\n")
  }

  cat(sprintf("\n minimum AIC\t%f\n", x$aicmin))
  cat(" parameters\n")
  i <- x$maice.order
  ii <- length(x$param[[i]])
  for( j in 1:ii ) {
    cat(sprintf("\t%e",x$param[[i]][j]))
    if( (j%%5 == 0) & j>1 ) cat("\n")
  }
  cat("\n")
}


pgraph <- function( data,mag,threshold=0.0,h,npoint,days,delta=0.0,dmax=0.0,separate.graphics=FALSE )
{
  nfunct <- 0
  isi <- 0
  n1 <- length(data)     #  total number of variable

  if(delta==0.0 || dmax==0.0) {
    dmax <- data[n1]/4
    delta <- dmax/100
  }
  td <- data[n1]/delta
  kmax <- as.integer(td/16.0 + 2)

  zd <- rep(0,n1)
  xmg <- rep(0,n1)
  xmg1 <- 0.0
  xmg2 <- 10.0
  nn <- 0
  for( i in 1:n1 )
    if( mag[i]==threshold || mag[i]>threshold )
    if( mag[i]==xmg2 || mag[i]<xmg2 )
    if( mag[i]==xmg1 || mag[i]>xmg1 )
    if( data[i]==0.0 || data[i]>0.0 ) {
        nn <- nn+1
        zd[nn] <- data[i]
        xmg[nn] <- mag[i]
#        if( xmg[nn]==0.0 ) xmg[nn]=6.0
    }
  zd <- zd[1:nn]
  xmg <- xmg[1:nn]

  nn1 <- nn-1
  xtau <- rep(0,2*nn)
  y <- rep(0,2*nn)
  kn <- 0
  xl <- rep(0,nn1)
  xx <- array(0, dim=c(nn1,6))
  ydev <- rep(0,nn1)
  ui <- rep(0,nn1)
  ncn <- rep(0,nn1)
  sui <- rep(0,nn1)
  xp <- rep(0,4)
  xrate <- rep(0,npoint)
  dlt <- 0
  xtime <- rep(0,kmax)
  yvar <- array(0, dim=c(5,kmax))
  sigma <- rep(0,kmax)
  k <- 0
  ier <- 0

  z <- .Fortran("pgraphf",
	as.integer(nfunct),
	as.integer(isi),
	as.double(zd),
	as.integer(nn),
	as.integer(npoint),
	as.double(days),
	as.double(h),
	as.double(delta),
	as.double(dmax),
	as.integer(kmax),
	xtau = as.double(xtau),
	y = as.double(y),
	kn = as.integer(kn),
	xl = as.double(xl),
	xx = as.double(xx),
	ydev = as.double(ydev),
	ui = as.double(ui),
	ncn = as.double(ncn),
	sui = as.double(sui),
	xp = as.double(xp),
	xrate = as.double(xrate),
	dlt = as.double(dlt),
	xtime = as.double(xtime),
	yvar = as.double(yvar),
	sigma = as.double(sigma),
	k = as.integer(k),
	ier = as.integer(ier) )

  kn <- z$kn
  k <- z$k

  tn <- data[n1]/as.double(nn)
  xtau <- z$xtau[1:kn]
  y <- z$y[1:kn]
  xl <- z$xl[1:nn1] * tn
  xx <- array(z$xx, dim=c(nn1,6)) * tn
  xx <- xx[1:nn1,1:6]
  ydev <- z$ydev[1:nn1]
  ui <- z$ui[1:nn1]
  ncn <- z$ncn[1:nn1]
  sui <- z$sui[1:nn1]
  xp <- z$xp
  xrate <- z$xrate
  dlt <- z$dlt
  xtime <- z$xtime[1:k]
  yvar <- array(z$yvar, dim=c(5,kmax))
  yvar <- yvar[1:5,1:k]
  sigma <- z$sigma[1:k]

  plot.pgraph( zd,xmg,h,kn,xtau,y,xl,xx,ydev,ui,ncn,sui,xp,npoint,dlt,xrate,k,xtime,sigma,yvar,separate.graphics )

  pgraph.out <- list( cnum=zd,lintv=xl,tau=xtau,nevent=y,survivor=xx,deviation=ydev,
normal.cnum=ncn,normal.lintv=ui,success.intv=sui,occur=xrate,time=xtime,variance=sigma,error=yvar )
  return( pgraph.out )
}


plot.pgraph <- function( zd,xmg,h,kn,xtau,y,xl,xx,ydev,ui,ncn,sui,xp,npoint,dlt,xrate,k,xtime,sigma,yvar,separate.graphics )
{
  nn <- length(zd)
  nn1 <- nn-1
  rnn <- as.double(nn)
  ymax <- nn
  xmax <- zd[nn]
  nband1 <- 1.35810*sqrt(rnn)
  nband2 <- 1.62762*sqrt(rnn)

# r.pgCumMT
  par(mfrow=c(2,1))
  plot(x=zd, y=c(1:nn), xlab="Time", ylab="Cumulative Number", type="s", main="Cumulative Curve & M-T plot", xlim=c(0,xmax), ylim=c(0,ymax))
  points(c(0,xmax),c(0,ymax), lty=2,  type="l")
  points(c(0,xmax),c(0,ymax)+nband1, col=2, lty=2, type="l")
  points(c(0,xmax),c(0,ymax)+nband2, col=2, lty=2, type="l")
  points(c(0,xmax),c(0,ymax)-nband1, col=2, lty=2, type="l")
  points(c(0,xmax),c(0,ymax)-nband2, col=2, lty=2, type="l")
  
  plot(x=zd, y=xmg, type="h", xlab="Time", ylab="Magnitude", xlim=c(0,xmax))

# r.pgPTnum
  par(ask=TRUE)
  if( separate.graphics == TRUE ) X11()
  par(mfrow=c(2,1))
  plot(x=xtau, y=y, type="l", ylab="Number of Events in [tau, tau+h]", xlab="tau = Time x (Total Number of Events) / (Time End)",xlim=c(0,nn))
  abline(h=seq(-3,3,1),lty=2, col=2)
  title(main=paste("Normalized number of point in [t, t+h], where h=",h,sep=""))
  plot(x=zd, y=xmg, type="h", ylab="Magnitude", xlab="Ordinary or Transformed Time",xlim=c(0,xmax))

# r.pgSurviv
  if( separate.graphics == TRUE ) X11()
  par(mfrow=c(2,1))
  yy <- log10(c(nn1:1))
  plot(x=xl, y=yy, type="p", pch="+", cex=0.5, xlab="Interval Length",ylab="Cumulative Number", axes=F, main="Survivor Function")
  points(xx[,1], yy, type="l", col=2)
  points(xx[,2], yy, type="l", col=2)
  points(xx[,3], yy, type="l", col=2)
  points(xx[,4], yy, type="l", col=2)
  points(xx[,5], yy, type="l", col=2)
  points(xx[,6], yy, type="l", col=2)
  box()
  axis(1)
  axis(2, at=log10(c(1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,20000,30000,40000,50000,60000,70000,8000,90000,100000)), labels=rep("", 46))
  axis(2, at=log10(c(1,10, 100, 1000,10000, 10000)), labels=c(expression(10^0),
  expression(10^1), expression(10^2), expression(10^3), 
  expression(10^4), expression(10^5)))

  xx <- c(1:nn1)
  plot(x=xx, y=ydev, type="n", xlab="Interval Length", ylab="Deviations", main="Deviation of Survivor Function from the Poisson")
  points(xx, ydev, type="p",pch="+",cex=0.5)
  abline(h=seq(-3,3,1),lty=2, col=2)
  abline(h=0, lty=1)

# r.pgInterP
# Inter1
  if( separate.graphics == TRUE ) X11()
  par(mfrow=c(2,1))
  xmax <- ui[1]
  ymax <- ncn[1]
  nband1 <- nband1/nn1
  nband2 <- nband2/nn1
  plot(x=ui, y=ncn, type="l", xlab="U(i) = exp{-(Normalized Interval Length)}", ylab="Normalized Cumulative Number", xlim=c(0,xmax), ylim=c(0,ymax), xaxs="i", yaxs="i", main="Interval-Length Distribution")
  points(c(0,xmax), c(0,ymax)+nband1, type="l", lty=2, col=2)
  points(c(0,xmax), c(0,ymax)+nband2, type="l", lty=2, col=2)
  points(c(0,xmax), c(0,ymax)-nband1, type="l", lty=2, col=2)
  points(c(0,xmax), c(0,ymax)-nband2, type="l", lty=2, col=2)
  points(c(0,xmax), c(0,ymax), type="l", lty=2, col=2)

# Inter2
  data<- matrix(, ncol=2, nrow=nn1-1)
  data[,1] <- sui[1:(nn1-1)]
  data[,2] <- sui[2:nn1]
  plot(data, type="p", pch="+", xlab="U(i)", ylab="U(i+1)", xaxs="i", yaxs="i", main="Successive Pair of Intervals")

# r.pgPalm
  if( separate.graphics == TRUE ) X11()
  par(mfrow=c(1,1))
  band <- xp
  data <- matrix(, ncol=2, nrow=npoint)
  data[,1] <- c(1:npoint)*dlt
  data[,2] <- xrate
  plot(data,type="p", pch="*", xlab="Elapsed Time After Each Event", ylab="Occurrence Rate",main="Palm Intensity")
  points(data,type="l")
  abline(h=band,lty=2)
  abline(h=mean(band))

# r.pgVTC
  if( separate.graphics == TRUE ) X11()
  plot(x=xtime, y=sigma, type="p", pch="+", cex=0.5, xlab="Time", ylab="Var{N(0,Time)}", 
       ylim=range(pretty(c(sigma,yvar))),main="Variance-Time Curve")
  points(xtime, yvar[1,], type="l", col=2)
  points(xtime, yvar[2,], type="l", lty=2, col=2)
  points(xtime, yvar[3,], type="l", lty=2, col=2)
  points(xtime, yvar[4,], type="l", lty=2, col=2)
  points(xtime, yvar[5,], type="l", lty=2, col=2)
  abline(h=0)
  abline(v=0)
  par(ask=FALSE)
}


ptspec <- function( data,nfre,prdmin,prd,nsmooth=1,pprd,interval,plot=TRUE )
{
# data                   #  data of events
  n <- length(data)      #  number of events of data set
# nfre                   #  number of sampling frequencies of spectra
  nh1 <- nfre+1
# prdmin                 #  the minimum periodicity of the sampling
# prd                    #  a periodicity for calculating the rayleigh probability  
# smooth                 #  number for smoothing of periodgram
  is <- nsmooth
# pprd                   #  particular periodcities to be investigated among others
  nt <- length(pprd)     #  number of particular periodcities to be investigated
# interval               #  length of observed time interval of events

  prb <- 0.0
  r1 <- 0.0
  rwx <- 0.0
  rwy <- 0.0
  phs <- 0.0
  wt <- rep(0,nt)
  ht <- rep(0,nt)
  w <- rep(0,nh1)
  h <- rep(0,nh1)
  g <- rep(0,nh1)

  z <- .Fortran("ptspecf",
	as.double(data),
	as.integer(n),
	as.double(interval),
	as.double(pprd),
	as.double(prdmin),
	as.double(prd),
	as.integer(nfre),
	as.integer(nt),
	as.integer(is),
	prb = as.double(prb),
	r1 = as.double(r1),
	rwx = as.double(rwx),
	rwy = as.double(rwy),
	phs = as.double(phs),
	wt = as.double(wt),
	ht = as.double(ht),
	w = as.double(w),
	h = as.double(h),
	g = as.double(g) )

  wt <- z$wt
  ht <- z$ht
  w <- z$w
  h <- z$h
  g <- z$g

  if( plot == TRUE ) {
    par(mfrow=c(2,1))
    pt.spec <- matrix(, nrow=nh1-1, ncol=2)
    pt.spec[,1] <- w[2:nh1]
    pt.spec[,2] <- h[2:nh1]
    ymin <- min(0.0, pt.spec[,2])
#    plot(pt.spec, type="n", xlab="Frequency", ylab="D.B.")
#    for( i in 1:nfre ) points(x=c(pt.spec[i,1],pt.spec[i,1]), y=c(pt.spec[i,2],ymin), type='l')
#    abline(h=0,lty=2)

    ram <- nh1/interval
    sgm <- log(interval/prdmin/pi)-log(0.05)
    sgm1 <- log10(ram*sgm)*10.e0
    sgm <- log(interval/prdmin/pi)-log(0.1)
    sgm2 <- log10(ram*sgm)*10.e0
    sgm <- log(interval/prdmin/pi)-log(0.01)
    sgm3 <- log10(ram*sgm)*10.e0
 
    confidence.bar <- c(sgm1,sgm2,sgm3)
    vertical.bar <- matrix(, nrow=nt, ncol=2)
    vertical.bar[,1] <- wt/(2*pi)
    vertical.bar[,2] <- ht
    pt.spec[,1] <- pt.spec[,1]/(2*pi)
    plot(pt.spec, pch=1, cex=0.5, xlab="Frequency", ylab="DB")
    points(vertical.bar, pch="+", cex=2,col=2)
    points(vertical.bar, lty=2, type="h")
    abline(h=confidence.bar, lty=3)

    plot(pt.spec[,1], 10^(pt.spec[,2]/10), pch=1, cex=0.5, xlab="Frequency", ylab="")
    points(vertical.bar[,1], 10^(vertical.bar[,2]/10), pch="+", cex=2,col=2)
    points(vertical.bar[,1], 10^(vertical.bar[,2]/10), lty=2, type="h")
    abline(h=10^(confidence.bar/10), lty=3)
    par(mfrow=c(1,1))
  }

  ptspec.out <- list( f=w,db=h,power=g,rayleigh.prob=z$prb,distance=z$r1,rwx=z$rwx,rwy=z$rwy,phase=z$phs,
                      n=n,nfre=nfre,prdmin=prdmin,nsmooth=nsmooth,interval=interval )

  class( ptspec.out ) <- "ptspec"
  return( ptspec.out )
}


print.ptspec <- function(x, ...)
{
  om <- 2*pi/x$prdmin
  cat(sprintf("\n\n maximum frequency= %e divided into %d points\n", om,x$nfre))

  ave <- x$n/x$interval
  cat(sprintf(" average= %e\n", ave))

  cat(sprintf("\n\n rayleigh probability = %f\n", x$rayleigh.prob))
  cat(sprintf(" distance = %f\n", x$distance))
  cat(sprintf(" rwx = %f\trwy = %f\n", x$rwx, x$rwy))
  cat(sprintf(" phase = %f\n", x$phase))

  nh1 <- length(x$f)
  cat("\n\n frequency\tD.B.\t\tpower\n")
  for( i in 2:nh1 )
    cat(sprintf(" %f\t%f\t%f\n", x$f[i], x$db[i], x$power[i]))

  if( x$nsmooth == 1 ) {
    xi <- x$f/(2*pi)
    lt4 <- 4/x$interval
    cat("\n\n list of high powers(>4) with period < t/4\n")
    cat(" frequency\tperiods\tpowers\n")
    for( i in 1:nh1 )
      if( xi[i]>=lt4 && x$power[i]>=4 ) {
        prd4 <- 1/xi[i]
        cat(sprintf(" %f\t%f\t%f\n", xi[i],prd4,x$power[i]))
      }
    cat("\n\n")
  }
}


linsim <- function( data,interval,c,d,ax,ay,at,ptmax )
{
  kxx <- length(ax)      # kxx-1 is the order of lgp trf; xx --> xx
  kxy <- length(ay)      # kxy-1 is the oeder of lgp trf; yy --> xx
  kxz <- length(at)      # kxz-1 is the order of polynomial for xx data
  t <- interval          # length of time interval in which events take place
# c,d                    # exponential coefficients in lgp corresponding to xx and xy, respectively
# ax,ay,at               # coefficients of the corresponding polynomials
  mm <- length(data)
  yy <- c(data,t)
# ptmax                  # ptxmax: an upper bound of trend polynomial
  kmax <- max(kxx,kxy,kxz)

  xx <- rep(0,2*mm)
  i1 <- 0
  j1 <- 0
  err <- 0.

  z <- .Fortran("linsimf",
	as.integer(kxx),
	as.integer(kxy),
	as.integer(kxz),
	as.double(t),
	as.double(c),
	as.double(d),
	as.double(ax),
	as.double(ay),
	as.double(at),
	as.double(yy),
	as.integer(mm),
	as.double(ptmax),
	as.integer(kmax),
	xx = as.double(xx),
	i1 = as.integer(i1),
	j1 = as.integer(j1),
	err = as.double(err) )

  simul <- z$xx[1:z$i1]
  input <- data[1:z$j1]

  if( z$err!=0. )
    cat(sprintf(" warning: is ptmax correct?  prob=%f\n",z$err))

  linsim.out <- list( in.data=input,sim.data=simul )
  return( linsim.out )
}


linlin <- function( external, self.excit, interval, c, d, ax=NULL, ay=NULL, ac=NULL, at=NULL, opt=0, tmp.file=NULL )
{
  yy <- external       # input data set
  xx <- self.excit     # self-exciting data set
  t <- interval        # length of time interval in which events take place

  x <- c(c,d)
  nn <- length(xx)
  mm <- length(yy)
  if( is.null(ax) ) {
    kkx <- 0
    ax <- 0.0
  }  else {
    x <- c(x,ax)
    kkx <- length(ax)  # kxx-1 is the order of lgp trf; xx --> xx
  }
  if( is.null(ay) ) {
    kky <- 0
    ay <- 0.0
  } else {
    x <- c(x,ay)
    kky <- length(ay)  # kxy-1 is the oeder of lgp trf; yy --> xx
  }
  if( is.null(ac) ) {
    kkc <- 0
    ac <- 0.0
  } else {
    x <- c(x,ac)
    kkc <- length(ac)  # kxz-1 is the order of polynomial for xx data
  }
  if( is.null(at)) {
    kkt <- 0
    at <- 0.0
  } else {
    x <- c(x,at)
    kkt <- length(at)  # kxz-1 is the order of polynomial for xx data
  }
  n <- length(x)

# opt                  # =0 : minimize the likelihood with fixed exponential coefficient c
                       # =1 :  not fixed d
  tmp.file <- paste(tmp.file, " ")

  x1 <- rep(0,n)
  x2 <- rep(0,n)
  aic <- 0.0
  f <- 0.0
  prb <- 0.0
  r1 <- 0.0
  rwx <- 0.0
  rwy <- 0.0
  phs <- 0.0

  z <- .Fortran("linlinf",
	as.integer(n),
	as.double(x),
	as.integer(opt),
	as.double(t),
	as.integer(nn),
	as.integer(mm),
	as.double(xx),
	as.double(yy),
	as.integer(kkx),
	as.integer(kky),
	as.integer(kkc),
	as.integer(kkt),
	x1 = as.double(x1),
	x2 = as.double(x2),
	aic = as.double(aic),
	f = as.double(f),
	prb = as.double(prb),
	r1 = as.double(r1),
	rwx = as.double(rwx),
	rwy = as.double(rwy),
	phs = as.double(phs),
	as.character(tmp.file) )

# initial estimates
  x1 <- z$x1
  c1 <- x1[1]
  d1 <- x1[2]
  ax1 <- NULL
  ay1 <- NULL
  ac1 <- NULL
  at1 <- NULL
  if( kkx>0 ) ax1 <- x1[3:(kkx+2)]
  if( kky>0 ) ay1 <- x1[(kkx+3):(kkx+kky+2)]
  if( kkc>0 ) ac1 <- x1[(kkx+kky+3):(kkx+kky+kkc+2)]
  if( kkt>0 ) at1 <- x1[(kkx+kky+kkc+3):n]

# final estimates
  x2 <- z$x2
  c2 <- x2[1]
  d2 <- x2[2]
  ax2 <- NULL
  ay2 <- NULL
  ac2 <- NULL
  at2 <- NULL
  if( kkx>0 ) ax2 <- x2[3:(kkx+2)]
  if( kky>0 ) ay2 <- x2[(kkx+3):(kkx+kky+2)]
  if( kkc>0 ) ac2 <- x2[(kkx+kky+3):(kkx+kky+kkc+2)]
  if( kkt>0 ) at2 <- x2[(kkx+kky+kkc+3):n]

  linlin.out <- list( c1=c1,d1=d1,ax1=ax1,ay1=ay1,ac1=ac1,at1=at1,c2=c2,d2=d2,ax2=ax2,ay2=ay2,ac2=ac2,at2=at2,
  aic2=z$aic,ngmle=z$f,rayleigh.prob=z$prb,distance=z$r1,rwx=z$rwx,rwy=z$rwy,phase=z$phs )

  class( linlin.out ) <- "linlin"
  return( linlin.out )
}

print.linlin <- function(x, ...)
{
  cat(sprintf("\n rayleigh probability = %f\n", x$rayleigh.prob))
  cat(sprintf(" distance = %f\n", x$distance))
  cat(sprintf(" rwx = %f\t\trwy = %f\n", x$rwx, x$rwy))
  cat(sprintf(" phase = %f\n", x$phase))

  kkx <- length(x$ax1)
  kky <- length(x$ay1)
  kkc <- length(x$ac1)
  kkt <- length(x$at1)

  cat("\n initial estimates\n")
  cat(sprintf(" c, d\t%f\t%f\n",x$c1,x$d1))
  if( kkx > 0 ) {
    cat(" ax")
    for( i in 1:kkx )  cat(sprintf("\t%e",x$ax1[i]))
    cat("\n")
  }
  if( kky > 0 ) {
    cat(" ay")
    for( i in 1:kky )  cat(sprintf("\t%e",x$ay1[i]))
    cat("\n")
  }
  if( kkc > 0 ) {
    cat(" ac")
    for( i in 1:kkc )  cat(sprintf("\t%e",x$ac1[i]))
    cat("\n")
  }
  if( kkt > 0 ) {
    cat(" at")
    for( i in 1:kkt )  cat(sprintf("\t%e",x$at1[i]))
    cat("\n")
  }

  cat("\n final outputs\n")
  cat(sprintf(" c, d\t%f\t%f\n",x$c2,x$d2))
  if( kkx > 0 ) {
    cat(" ax")
    for( i in 1:kkx )  cat(sprintf("\t%e",x$ax2[i]))
    cat("\n")
  }
  if( kky > 0 ) {
    cat(" ay")
    for( i in 1:kky )  cat(sprintf("\t%e",x$ay2[i]))
    cat("\n")
  }
  if( kkc > 0 ) {
    cat(" ac")
    for( i in 1:kkc )  cat(sprintf("\t%e",x$ac2[i]))
    cat("\n")
  }
  if( kkt > 0 ) {
    cat(" at")
    for( i in 1:kkt )  cat(sprintf("\t%e",x$at2[i]))
    cat("\n")
  }

  cat(sprintf("\n negative max likelihood = %f\n",x$ngmle))
  cat(sprintf(" AIC/2 = %f\n\n",x$aic2))
}

simbvh <- function( interval,axx=NULL,axy=NULL,axz=NULL,ayx=NULL,ayy=NULL,ayz=NULL,c,d,c2,d2,ptxmax,ptymax )
{
  t <- interval          # length of time interval in which events take place

# axx(i),axy(i),ayx(i),ayy(i),axz(i),ayz(i) # coefficients of the corresponding polynomials
  if( is.null(axx) ) {
    axx <- 0
    kxx <- 0
  } else {
    kxx <- length(axx)     # kxx-1 is the order of lgp trf; xx --> xx
  }
  if( is.null(axy) ) {
    axy <- 0
    kxy <- 0
  } else {
    kxy <- length(axy)     # kxy-1 is the oeder of lgp trf; yy --> xx
  }
  if( is.null(ayx) ) {
    ayx <- 0
    kyx <- 0
  } else {
    kyx <- length(ayx)     # kyx-1 is the oeder of lgp trf; xx --> yy
  }
  if( is.null(ayy) ) {
    ayy <- 0
    kyy <- 0
  } else {
    kyy <- length(ayy)     # kyy-1 is the oeder of lgp trf; yy --> yy
  }
  if( is.null(axz) ) {
    axz <- 0
    kxz <- 0
  } else {
    kxz <- length(axz)     # kxz-1 is the order of polynomial for xx data
  }
  if( is.null(ayz) ) {
    ayz <- 0
    kyz <- 0
  } else {
    kyz <- length(ayz)     # kyz-1 is the order of polynomial for yy data
  }

# c,d,c2,d2              # exponential coefficients in lgp corresponding to
                         # xx,xy,yx and yy, respectively

# ptxmax                 # upper bounds of trend polynomials corresponding to xz
# ptymax                 # upper bounds of trend polynomials corresponding to yz
  kmax <- max(kxx,kxy,kyx,kyy)

  nnmax <- 10000
  mmmax <- 10000
  xx <- rep(0,nnmax)
  yy <- rep(0,mmmax)
  i1 <- 0
  j1 <- 0
  err <- 0.0

  z <- .Fortran("simbvhf",
	as.integer(kxx),
	as.integer(kxy),
	as.integer(kxz),
	as.integer(kyx),
	as.integer(kyy),
	as.integer(kyz),
	as.double(t),
	as.double(c),
	as.double(d),
	as.double(c2),
	as.double(d2),
	as.double(axx),
	as.double(axy),
	as.double(axz),
	as.double(ayx),
	as.double(ayy),
	as.double(ayz),
	as.double(ptxmax),
	as.double(ptymax),
	as.integer(kmax),
	xx = as.double(xx),
	yy = as.double(yy),
	i1 = as.integer(i1),
	j1 = as.integer(j1),
	err = as.double(err),
	as.integer(nnmax),
	as.integer(mmmax) )

  if( z$err!=0. )
    cat(sprintf(" warning: are ptxmax & ptymax correct? prob=%f\n",z$err))

  simbvh.out <- list( x=z$xx[1:z$i1], y=z$yy[1:z$j1] )
  return( simbvh.out )
}


momori <- function( data,mag=NULL,threshold=0.0,tstart,tend,parami,tmp.file=NULL )
{

# data                   #  point process data
  n <- length(data)      #  total number of variable
# mag                    #  magnitude
# threshold              #  threshold magnitude
# tstart                 #  the start of the target perio
  zts=tstart
  zte=tend
# tend                   #  the end of the target period
# parami                 #  initial estimates of parameters
  np <- length(parami)
  np1 <- np+1
  if( np1 != 5 ) {
    cat(" warning: number of parameters is worng.")
    return()
  }
  parami <- c(parami[1:3],0,parami[4])
  tmp.file <- paste(tmp.file, " ")

  nn <- 0
  xx <- NULL
  magnitude <- NULL

  if( is.null(mag) ) {
    xx <- data
    nn <- n
  } else {
    for( i in 1:n ) 
      if( mag[i] >= threshold ) 
        if( (data[i]>=zts) && (data[i]<=zte) ) {
          nn <- nn+1
          xx <- c(xx,data[i])
          magnitude <- c(magnitude,mag[i])
        }
  }

  nfunct <- 6
  ncount <- 1
  kn <- (np-1)/3
  ff <- 0.0
  x1 <- rep(0,np)
  x <- rep(0,np)
  ahaic <- rep(0,ncount)
#  kn <- (np-1)/3
  kn <- np/3
  if( kn==0 ) kn <- 1
  t0 <- 0
  ti <- rep(0,kn)
  ak <- rep(0,kn)
  c <- rep(0,kn)
  p <- rep(0,kn)
  cls <- rep(0,kn)

  z <- .Fortran("momorif",
	as.double(xx),
	as.integer(nn),
	as.double(parami),
	as.integer(np),
	as.double(zts),
	as.double(zte),
	as.double(tstart),
	as.integer(ncount),
	as.integer(nfunct),
	ff=as.double(ff),
	x1=as.double(x1),
	x=as.double(x),
	ahaic=as.double(ahaic),
	t0=as.double(t0),
	ti=as.double(ti),
	ak=as.double(ak),
	c=as.double(c),
	p=as.double(p),
	cls=as.double(cls),
	as.character(tmp.file) )

#  x0 <- 0
#  param2 <- c(t=z$x[1], k=z$x[2], c=z$x[3], p=x0, cls=z$x[4])
  pa2 <- z$x

  ti <- c(z$t0,z$ti[1:kn-1])
  pa <- list(t_i=ti, K=z$ak, c=z$c, p=z$p, cls=z$cls)
  momori.out <- list( param1=z$x1, param2=pa2, ngmle=z$ff, aic=z$ahaic, param=pa )
  return( momori.out )
}


respoi <- function( time,mag,param,zts,tstart,zte,threshold=0.0,plot=TRUE )
{
# time                   #  time from the main shock in days
  nd <- length(time)
# mag                    #  magnitude

# param                  #  estimates of parameters
  np <- length(param) + 1
  if( np != 5 ) {
    cat(" warning: number of parameters is worng.")
    return()
  }
  param <- c(param[1:3],0,param[4])

# zts                    #  the start of the precursory period
# zstart                 #  the start of the target period
# zte                    #  the end of the target period

  dep <- rep(0,nd)
  xp <- rep(0,nd)
  yp <- rep(0,nd)
  ntstar <- 0
  xx <- rep(0,nd)
  x <- rep(0,nd)
  nn <- 0

  z <- .Fortran("respoif",
	as.double(time),
	mag=as.double(mag),
	as.double(dep),
	as.double(xp),
	as.double(yp),
	as.integer(nd),
	as.double(param),
	as.integer(np),
	as.double(zts),
	as.double(zte),
	as.double(tstart),
	as.double(threshold),
	ntstar=as.integer(ntstar),
	xx=as.double(xx),
	x=as.double(x),
	nn=as.integer(nn))

  n <- z$nn
  cn <- 1:n - z$ntstar
  ti <- z$x[1:n]
  mag <- z$mag[1:n]

  if( plot == TRUE ) {
    mag1 <- threshold
    level <- min(0, cn[1])
    xrange <- c(min(ti), max(ti))
    mgrange <- max(cn)/4
    bottom <- min(cn) - mgrange
    yrange <- c(bottom, max(max(cn),max(ti)))
    plot(xrange, yrange, type = "n", 
       main = "Omori-Utsu Residual", xlab="Transformed Time", 
       ylab = "Cumulative Number of Events", lwd=1, pty="s",xaxs='r')
    points(ti, 1:length(ti)+cn[1]+1, type='s') 
    mgmax <- max(mag - mag1 + 1)
    mag <- mag - mag1 + 0.5
    segments(ti, bottom, ti, mag/mgmax * mgrange + bottom)	
    abline(h=bottom)
    abline(h=0)
    abline(v=0,lty=2)
    abline(0,1,lty=1,col='red')
    s <- tstart
  }

  respoi.out <- list( trans.time=ti, cnum=cn )
  return( respoi.out )
}


respoi2 <- function( etas,param,zts,tstart,zte,threshold=0.0,plot=TRUE )
{
  if( dim(etas)[2] != 9 ) {
    cat(" warning: etas is a sequential data with 9 columns-format.")
    return()
  }

  nd <- dim(etas)[1]
  xp <- etas[,2]         #  longitude
  yp <- etas[,3]         #  latitude
  mag <- etas[,4]        #  magnitude
  time <- etas[,5]       #  time from the main shock in days
  dep <- etas[,6]        #  depth

# param                  #  estimates of parameters
  np <- length(param) + 1
  if( np != 5 ) {
    cat(" warning: number of parameters is worng.")
    return()
  }
  param <- c(param[1:3],0,param[4])

# zts                    #  the start of the precursory period
# zstart                 #  the start of the target period
# zte                    #  the end of the target period

  ntstar <- 0
  xx <- rep(0,nd)
  x <- rep(0,nd)
  nn <- 0

  z <- .Fortran("respoif",
	as.double(time),
	mag=as.double(mag),
	dep=as.double(dep),
	xp=as.double(xp),
	yp=as.double(yp),
	as.integer(nd),
	as.double(param),
	as.integer(np),
	as.double(zts),
	as.double(zte),
	as.double(tstart),
	as.double(threshold),
	ntstar=as.integer(ntstar),
	xx=as.double(xx),
	x=as.double(x),
	nn=as.integer(nn))

  n <- z$nn
  cn <- 1:n - z$ntstar
  ti <- z$x[1:n]
  mag <- z$mag[1:n]

  if( plot == TRUE ) {
    mag1 <- threshold
    level <- min(0, cn[1])
    xrange <- c(min(ti), max(ti))
    mgrange <- max(cn)/4
    bottom <- min(cn) - mgrange
    yrange <- c(bottom, max(max(cn),max(ti)))
    plot(xrange, yrange, type = "n", 
       main = "Omori-Utsu Residual", xlab="Transformed Time", 
       ylab = "Cumulative Number of Events", lwd=1, pty="s", xaxs='r')
    points(ti, 1:length(ti)+cn[1]+1, type='s') 
    mgmax <- max(mag - mag1 + 1)
    mag <- mag - mag1 + 0.5
    segments(ti, bottom, ti, mag/mgmax * mgrange + bottom)	
    abline(h=bottom)
    abline(h=0)
    abline(v=0,lty=2)
    abline(0,1,lty=1,col='red')
    s <- tstart
  }

  res <- matrix(, ncol=7, nrow=n)
  res[,1] <- cn
  res[,2] <- z$xp[1:n]
  res[,3] <- z$yp[1:n]
  res[,4] <- z$mag[1:n]
  res[,5] <- z$xx[1:n]
  res[,6] <- z$dep[1:n]
  res[,7] <- ti
  res <- data.frame(res)
  names(res) <- c("no.", "longitude", "latitude", "magnitude", "time", "depth", "trans.time")

  respoi2.out <- list(resData=res)
  return( respoi2.out )

}


etasap <- function( time,mag,threshold=0.0,reference=0.0,parami,zts=0.0,tstart,zte,approx=2,tmp.file=NULL,plot=TRUE )
{
# time                   #  time from the main shock in days
  n <- length(time)
# mag                    #  magnitude
# threshold              #  threshold magnitude
# reference              #  reference magnitude
# parami                 #  initial estimates of parameters
  np <- length(parami)
  if( np != 5 ) {
    cat(" warning: number of parameters is worng.")
    return()
  }
# zts                    #  the start of the precursory period
# tstart                 #  the start of the target period
# zte                    #  the end of the target period
# approx                 #  the level for approximation version (1, 2, 4, 8, 16)
                         #  =0 : the exact version
  tmp.file <- paste(tmp.file, " ")

  nn <- 0
  xx <- NULL
  magnitude <- NULL
  amx1 <- threshold
  amx2 <- 10.0

#  if( is.null(mag) ) {
#    for( i in 1:n )
#      if( time[i]>=zts && time[i]<=zte ) {
#          nn <- nn+1
#          xx <- c(xx,time[i])
#      }
#  } else {
    for( i in 1:n ) 
      if( mag[i]>=amx1 && mag[i]<=amx2 )
      if( time[i]>=zts && time[i]<=zte ) {
          nn <- nn+1
          xx <- c(xx,time[i])
          magnitude <- c(magnitude,mag[i])
     }
#  }

  nfunct <- 9
  if( approx == 0 ) nfunct <- 2

  f <- 0.0
  x <- rep(0,np)
  aic2 <- 0.0


  z <- .Fortran("etasapf",
	as.double(xx),
	as.double(magnitude),
	as.integer(nn),
	as.double(reference),
	as.double(threshold),
	as.double(parami),
	as.integer(np),
	as.double(zts),
	as.double(zte),
	as.double(tstart),
	as.integer(nfunct),
	as.integer(approx),
	f=as.double(f),
	x=as.double(x),
	aic2=as.double(aic2),
	as.character(tmp.file) )

  if( plot == TRUE ) {
    mag1 <- min(magnitude)
    ti <- xx
    zero <- rep(0,nn)
    cn <- 1:nn
    xrange <- c(min(xx), max(xx))
    mgrange <- max(cn)/4
    bottom <- min(cn) - mgrange
    yrange <- c(bottom, max(cn))
     plot(xrange, yrange, type = "n", main = "Seismicity CUM & M-T plot", 
         xlab = "Ordinary Time (Days)", 
         ylab = "Cumulative Number of Events", lwd=1, pty="s", xaxs='r', yaxs='i')
    lines(ti, cn, type='S')
    mgmax <- max(magnitude - mag1 + 1)
    magnitude1 <- magnitude - mag1 + 0.5
    segments(xx, bottom, xx, magnitude1/mgmax * mgrange + bottom)	
    abline(h = 0)
    abline(h = bottom)
  }

  pa <- list(B=z$x[1], K=z$x[2], c=z$x[3], p=z$x[4], cls=z$x[5])
  etasap.out <- list( ngmle=z$f, aic2=z$aic2, param=pa )
  return( etasap.out )
}


etarpp <- function( time,mag,threshold=0.0,reference=0.0,parami,zts=0.0,tstart,zte,ztend=NULL, tmp.file=NULL,plot=TRUE )
{
# time                   #  time from the main shock in days
  n <- length(time)
# mag                    #  magnitude
# threshold              #  threshold magnitude
# reference              #  reference magnitude
# parami                 #  initial estimates of parameters
  np <- length(parami)
  if( np != 5 ) {
    cat(" warning: number of parameters is worng.")
    return()
  }
# zts                    #  the start of the precursory period
# tstart                 #  the start of the target period
# zte                    #  the end of the target period
# ztend                  #  the end of the predictionperiod
  if( is.null(ztend) ) ztend <- time[n]
  tmp.file <- paste(tmp.file, " ")

  mm <- 0
  nn <- 0
  time0 <- NULL
  mag0 <- NULL
  time1 <- NULL
  mag1 <- NULL
  for( i in 1:n ) 
    if( mag[i] >= threshold ) {
      mm <- mm+1
      time0 <- c(time0,time[i])
      mag0 <- c(mag0,mag[i])
      if( time[i]>=zts && time[i]<=ztend ) {
        nn <- nn+1
        time1 <- c(time1,time[i])
        mag1 <- c(mag1,mag[i])
      }
    }

  x <- rep(0,nn)
  ntstar <- 0

  z <- .Fortran("etarppf",
	as.double(time1),
	as.double(mag1),
	as.double(reference),
	as.integer(nn),
	as.double(parami),
	as.integer(np),
	as.double(zts),
	as.double(ztend),
	as.double(tstart),
	x=as.double(x),
	ntstar=as.integer(ntstar),
	as.character(tmp.file) )

  cnum <- 1:nn - z$ntstar

  if( plot == TRUE ) {
# r.seisetas
## scan("work.etas")
    mgmin <- min(mag0)
    zero <- rep(0,mm)
    cn <- 1:mm
    cna <- append(cn,0,after=0)
    cn1 <- cna[1:mm]
    tia <- append(time0,0,after=0)
    ti1 <- tia[1:mm]
    xrange <- c(min(time0),ztend)
    mgrange <- max(cn)/4
    bottom <- min(cn)-mgrange
    yrange <- c(bottom,max(max(cn),max(z$x-z$ntstar)))
#    plot(xrange,yrange,type="n",main="ETAS Fit and Prediction",
    plot(xrange,yrange,type="n", main=paste("ETAS Fit and Prediction\nM>=",threshold,"S=",zts,"T=",zte,"Tend=",ztend),
         xlab="Ordinary Time (Days)", ylab="Cumulative Number of Events", lwd=1, xlim=c(xrange), pty="s", xaxs="r",yaxs="r")
## scan("work.res")
    lines(time1,z$x+z$ntstar,type='l',col='red')
    segments(ti1,cn1,time0,cn1)
    segments(time0,cn1,time0,cn)
    mgmax <- max(mag0-mgmin+1)
    mag2 <- mag0-mgmin+0.5
    segments(time0,bottom,time0,mag2/mgmax*mgrange+bottom)
    abline(h=0)
    abline(h=bottom)
    mark0 <- tstart; abline(v=mark0,lty=2)
#   mark1 <- ngmle ; abline(v=mark1,lty=2)
    mark2 <- ztend;  abline(v=mark2,lty=2)
    abline(v=tstart,lty=2)
    abline(v=zte,lty=2)
#    text(max(time0)*0.3, max(cn)*0.95, paste('M>=',mgmin,'S=',zts,'T=',zte,'Tend=',ztend))

# r.retas
## scan("work.res")
#    X11()
    par(ask=TRUE)
    level <- min(0, cnum[1])
    cn <- 1:mm+cnum[1]
    ti <- z$x
    xrange <- c(min(ti),max(ti))
    mgrange <- max(cn/4)
    bottom <- min(cn)-mgrange
    yrange <- c(bottom, max(max(cn),max(ti)))
#    plot(xrange,yrange,type="n",main="ETAS Residual",
    plot(xrange,yrange,type="n",main=paste("ETAS Residual\nM>=",threshold,"S=",zts,"T=",zte,"Tend=",ztend),
         xlab="Transformed Time",ylab="Cumulative Number of Events",lwd=1, pty="s", xaxs="r")
    points(z$x,1:nn+cnum[1]+1,type='s')
    mgmax <- max(mag1-threshold+1)
    mag3 <- mag1-threshold+0.5
    segments(ti,bottom,ti,mag3/mgmax*mgrange+bottom)
    abline(h=bottom)
    abline(h=0)
    abline(v=0,lty=2)
    abline(0,1,lty=1,col='red')
    timax <- max(ti[time1<=zte])
    abline(v=timax,lty=2)
#    text(max(ti)*0.4,max(cn)*0.9,paste('M>=',threshold,'S=',zts,'T=',zte,'Tend=',ztend))
    par(ask=FALSE)
  }

  etarpp.out <- list( trans.time=z$x, no.tstart=z$ntstar )
  return( etarpp.out )
}


etarpp2 <- function( etas,threshold=0.0,reference=0.0,parami,zts=0.0,tstart,zte,ztend=NULL,tmp.file=NULL,plot=TRUE )
{
  n <- dim(etas)[1]
  xp <- etas[,2]        #  longitude
  yp <- etas[,3]        #  latitude
  mag <- etas[,4]       #  magnitude
  time <- etas[,5]      #  time from the main shock in days
  dep <- etas[,6]       #  depth

# threshold              #  threshold magnitude
# reference              #  reference magnitude
# parami                 #  initial estimates of parameters
  np <- length(parami)
  if( np != 5 ) {
    cat(" warning: number of parameters is worng.")
    return()
  }
# zts                    #  the start of the precursory period
# tstart                 #  the start of the target period
# zte                    #  the end of the target period
# ztend                  #  the end of the predictionperiod
  if( is.null(ztend) ) ztend <- time[n]
  tmp.file <- paste(tmp.file, " ")

  mm <- 0
  nn <- 0
  time0 <- NULL
  mag0 <- NULL
  xp1 <- NULL
  yp1 <- NULL
  mag1 <- NULL
  time1 <- NULL
  dep1 <- NULL
  for( i in 1:n ) 
    if( mag[i] >= threshold ) {
      mm <- mm+1
      time0 <- c(time0,time[i])
      mag0 <- c(mag0,mag[i])
      if( time[i]>=zts && time[i]<=ztend ) {
        nn <- nn+1
        xp1 <- c(xp1,xp[i])
        yp1 <- c(yp1,yp[i])
        mag1 <- c(mag1,mag[i])
        time1 <- c(time1,time[i])
        dep1 <- c(dep1,dep[i])
      }
    }

  x <- rep(0,nn)
  ntstar <- 0

  z <- .Fortran("etarppf",
	as.double(time1),
	as.double(mag1),
	as.double(reference),
	as.integer(nn),
	as.double(parami),
	as.integer(np),
	as.double(zts),
	as.double(ztend),
	as.double(tstart),
	x=as.double(x),
	ntstar=as.integer(ntstar),
	as.character(tmp.file) )

  cnum <- 1:nn - z$ntstar

  if( plot == TRUE ) {
# r.seisetas
## scan("work.etas")
    mgmin <- min(mag0)
    zero <- rep(0,mm)
    cn <- 1:mm
    cna <- append(cn,0,after=0)
    cn1 <- cna[1:mm]
    tia <- append(time0,0,after=0)
    ti1 <- tia[1:mm]
    par(pty="s", xaxs="r",yaxs="r")
    xrange <- c(min(time0),ztend)
    mgrange <- max(cn)/4
    bottom <- min(cn)-mgrange
    yrange <- c(bottom,max(max(cn),max(z$x-z$ntstar)))
#    plot(xrange,yrange,type="n",main="ETAS Fit and Prediction",
    plot(xrange,yrange,type="n",main=paste("ETAS Fit and Prediction\nM>=",threshold,"S=",zts,"T=",zte,"Tend=",ztend),
         xlab="Ordinary Time (Days)",ylab="Cumulative Number of Events",lwd=1,xlim=c(xrange))
## scan("work.res")
    lines(time1,z$x+z$ntstar,type='l',col='red')
    segments(ti1,cn1,time0,cn1)
    segments(time0,cn1,time0,cn)
    mgmax <- max(mag0-mgmin+1)
    mag2 <- mag0-mgmin+0.5
    segments(time0,bottom,time0,mag2/mgmax*mgrange+bottom)
    abline(h=0)
    abline(h=bottom)
    mark0 <- tstart; abline(v=mark0,lty=2)
#   mark1 <- ngmle ; abline(v=mark1,lty=2)
    mark2 <- ztend;  abline(v=mark2,lty=2)
    abline(v=tstart,lty=2)
    abline(v=zte,lty=2)
#    text(max(time0)*0.3, max(cn)*0.95, paste('M>=',mgmin,'S=',zts,'T=',zte,'Tend=',ztend))

# r.retas
## scan("work.res")
#    X11()
    par(ask=TRUE)
    level <- min(0, cnum[1])
    cn <- 1:mm+cnum[1]
    ti <- z$x
    xrange <- c(min(ti),max(ti))
    mgrange <- max(cn/4)
    bottom <- min(cn)-mgrange
    yrange <- c(bottom, max(max(cn),max(ti)))
    par(pty="s",xaxs="r")
#    plot(xrange,yrange,type="n",main="ETAS Residual",
    plot(xrange,yrange,type="n",main=paste("ETAS Residual\nM>=",threshold,"S=",zts,"T=",zte,"Tend=",ztend),
         xlab="Transformed Time",ylab="Cumulative Number of Events",lwd=1)
    points(z$x,1:nn+cnum[1]+1,type='s')
    mgmax <- max(mag1-threshold+1)
    mag3 <- mag1-threshold+0.5
    segments(ti,bottom,ti,mag3/mgmax*mgrange+bottom)
    abline(h=bottom)
    abline(h=0)
    abline(v=0,lty=2)
    abline(0,1,lty=1,col='red')
    timax <- max(ti[time1<=zte])
    abline(v=timax,lty=2)
#    text(max(ti)*0.4,max(cn)*0.9,paste('M>=',threshold,'S=',zts,'T=',zte,'Tend=',ztend))
    par(ask=FALSE)
  }

  res <- matrix(, ncol=7, nrow=nn)
  res[,1] <- cnum
  res[,2] <- xp1
  res[,3] <- yp1
  res[,4] <- mag1
  res[,5] <- time1
  res[,6] <- dep1
  res[,7] <- z$x
  res <- data.frame(res)
  names(res) <- c("no.", "longitude", "latitude", "magnitude", "time", "depth", "trans.time")

  etarpp2.out <- list(resData=res)
  return( etarpp2.out )
}


etasim1 <- function( bvalue,nd,threshold=0.0,reference=0.0,param )
{
# bvlalue                #  b-value of G-R law
# nd                     #  the number of the simulated events
# threshold              #  threshold magnitude
# reference              #  reference magnitude
# param                  #  initial estimates of parameters
  np <- length(param)
  if( np != 5 ) {
    cat(" warning: number of parameters is worng.")
    return()
  }
  a <- param[1]
  b <- param[2]
  c <- param[3]
  d <- param[4]
  p <- param[5]

  ic <- 0
  tstart <- 0
  xm <- rep(0,nd)
  zz <- rep(0,nd)
  xx <- rep(0,nd)
  probx <- 0

  z <- .Fortran("etasimf",
	as.integer(ic),
	as.double(bvalue),
	as.double(tstart),
	as.integer(nd),
	as.double(threshold),
	as.double(reference),
	as.double(a),
	as.double(b),
	as.double(c),
	as.double(d),
	as.double(p),
	xm=as.double(xm),
	as.double(zz),
	xx=as.double(xx),
	probx=as.double(probx) )

  sim <- matrix(, ncol=9, nrow=nd)
  sim[,1] <- 1:nd
  sim[,2] <- rep(0,nd)
  sim[,3] <- rep(0,nd)
  sim[,4] <- z$xm
  sim[,5] <- z$xx
  sim[,6] <- rep(0,nd)
  sim[,7] <- rep(0,nd)
  sim[,8] <- rep(0,nd)
  sim[,9] <- rep(0,nd)
  sim <- data.frame(sim)
  names(sim) <- c("no.", "longitude", "latitude", "magnitude", "time", "depth", "year", "month", "days")

  return( etasim=sim )
}


etasim2 <- function( etas,tstart,threshold=0.0,reference=0.0,param )
{
  n <- dim(etas)[1]
  mag <- etas[,4]       #  magnitude
  time <- etas[,5]      #  time from the main shock in days
# tstart                #  
# threshold             #  threshold magnitude
# reference             #  reference magnitude
# param                 #  initial estimates of parameters
  np <- length(param)
  if( np != 5 ) {
    cat(" warning: number of parameters is worng.")
    return()
  }
  a <- param[1]
  b <- param[2]
  c <- param[3]
  d <- param[4]
  p <- param[5]

  nd <- 0
  mag1 <- NULL
  time1 <- NULL
  for( i in 1:n )
      if( mag[i] >= threshold ) {
        nd <- nd+1
        mag1 <- c(mag1,mag[i])
        time1 <- c(time1,time[i])
      }

  ic <- 1
  bvalue <- 0
#  time1 <- rep(0,nd)
  xx <- rep(0,nd)
  probx <- 0

  z <- .Fortran("etasimf",
	as.integer(ic),
	as.double(bvalue),
	as.double(tstart),
	as.integer(nd),
	as.double(threshold),
	as.double(reference),
	as.double(a),
	as.double(b),
	as.double(c),
	as.double(d),
	as.double(p),
	as.double(mag1),
	as.double(time1),
	xx=as.double(xx),
	probx=as.double(probx) )

  if( z$probx > 1 ) {
    cat(sprintf("\n prob>1%f\n", z$probx))
    return()
  } else {
    sim <- matrix(, ncol=9, nrow=nd)
    sim[,1] <- 1:nd
    sim[,2] <- rep(0,nd)
    sim[,3] <- rep(0,nd)
    sim[,4] <- mag1
    sim[,5] <- z$xx
    sim[,6] <- rep(0,nd)
    sim[,7] <- rep(0,nd)
    sim[,8] <- rep(0,nd)
    sim[,9] <- rep(0,nd)
    sim <- data.frame(sim)
    names(sim) <- c("no.", "longitude", "latitude", "magnitude", "time", "depth", "year", "month", "days")

    return( etasim2=sim )
  }
}

.noGenerics <- TRUE

options(warn.FPU=FALSE)
