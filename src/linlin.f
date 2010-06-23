cc     program linlin
      subroutine linlinf(n,x,iopt,t,nn,mm,xx,yy,kkx,kky,kkc,kkt,x1,x2,
     & aic,f,prb,r1,rwx,rwy,phs,tmp)
c
      include 'sapp.h'
c
c     this program performs the maximum likelihood estimates of linear
c  intensity models of self-exciting point process with another point
c  process input, cyclic and trend components.  the cyclic part is
c  given by the fourier series, the trend is given by usual polynomial.
c  the response functions of the self-exciting and the input are given
c  by the laguerre type polynomials, where the scaling parameters in the
c  exponential function, say c and d, can be different.  however it is
c  adviced to estimate c first without the input component, and then to
c  estimate d with the fixed c (this means taht the gradient
c  corresponding to the c is set to keep 0), which are good initial
c  estimates for the c and d of the mixed self-exciting and input model.
c     it is recommended to use akaike information criterion (aic) for
c  selecting optimal orders of respective components.
c     finally it should be noted that estimated intensitysometimes
c  happen to be negative on some part of time interval outside the
c  neighbourhood of events.  this take place more easily the larger the
c  number of parameters.  this causes some difficulty in getting the
c  m.l.e., because the negativity of the intensity contributes to the
c  seeming increase of the likelihood.  see the references below.
c
c     structure of the program
c
c          linlin
c             i---input
c             i---comfac
c             i---cycle
c             +---dav-------davidn----------funct1 or 2
c                              i---hesian
c                              +---linear---funct1 or 2
c
c     this program is designed by y. ogata, and programmed by
c  y. ogata and k. katsura, inst. statist. math., tokyo.  (31/01/85)
c
c  references
c  y. ogata & h. akaike (1982). "on linear intensity models for mixed
c     doubly stochastic poisson and self-exciting point processes." j.
c     royal statist. soc. b, vol. 44, pp. 102-107.
c  y. ogata, h. akaike and k. katsura (1982). "the application of linear
c     intensity models to the investigation of causal relations between
c     a point process and another stochastic process."  ann. inst.
c     statist. math., vol. 34. pp. 373-387.
c  y. ogata and k. katsura (1984). "point process model with linearly
c     parametrized intensity for the application to earthquake data."
c     essays in time series and allied processes (festscrift for prof.
c     e. j. hannan), j. gani and m. b. priestley eds., j. appl. probab.
c     vol. 23a, to appear.
c
      implicit real * 8 (a-h,o-z)
cc      dimension x(50)
      dimension x(n), x1(n), x2(n)
      dimension xx(nn),yy(nn)
      dimension lf(51,51)
c
      integer*1  tmp(1)
      character cname*80
c
      lu=3
      do 5 i = 1,80
    5 cname(i:i) = ' '
      ic = 0
      ifg = 0
      if ( tmp(1).ne.ichar(' ') ) then
      do 6 i = 1,80
         if ( tmp(i).ne.ichar(' ') ) then
            cname(i:i) = char(tmp(i))
            ic = ic+1
         end if
    6 continue
      end if
      if ( ic.gt.0 ) then
         ifg = 1
         open (lu,file=cname,iostat=ivar)
         if (ivar .ne. 0) then
            write(*,*) ' ***  linlin temp FILE OPEN ERROR :',cname,ivar
            ifg=0
         end if
      end if
c
cc      call input(n,x)
cc      call comfac
      call comfac(lf)
cc      call cycle
      prd=365.25d0
      call cycle(xx,nn,prd,prb,r1,rwx,rwy,phs)
cc      call dav(n,x)
      do 10 i=1,n
      x2(i)=x(i)
   10 continue
      call dav(n,x1,xx,yy,nn,kkx,kky,kkc,kkt,t,mm,iopt,lf,x2,aic,f,xm,
     & lu,ifg)
      if (ifg.eq.1) close(lu)
cc      stop
      return
      end
cc      subroutine dav(n,x)
      subroutine dav(n,x1,xx,yy,nn,kkx,kky,kkc,kkt,t,mm,iopt,lf,
     & x,aic,f,xm,lu,ifg)
      implicit real * 8 (a-h,o-z)
cc      external funct
cc      common t,nn,mm,iopt
cc      common /ddd/r,f,aic,sd
cc      common /kkxy/kkx,kky
cc      common /ct/kkc,kkt
cc      dimension x(50),r(31,31)
      dimension x(n),x1(n),x2(n)
      dimension xx(nn),yy(nn)
      if(n.eq.1) go to 100
      x(1)=sqrt(x(1))
      x(2)=sqrt(x(2))
      k2=2
      if(kkx.ne.0) x(k2+1)=sqrt(x(k2+1))
      if(kky.ne.0) x(kkx+k2+1)=sqrt(x(kkx+k2+1))
      if(kkt.ne.0) x(kkx+kky+kkc+k2+1)=sqrt(x(kkx+kky+kkc+k2+1))
c
c     polynomial(trend) coefficient's transformation
c
      if(kkt.eq.1) go to 42
      do 43 k=2,kkt
      x(k+kkx+kky+kkc+k2)=x(k+kkx+kky+kkc+k2)*t**(k-1)
   43 continue
   42 continue
cc      write(6,1020) n,(x(i),i=1,n)
      do 40 i=1,n
      x1(i)=x(i)
   40 continue
      do 30 ii=1,5
c----------------------------------
cc      call davidn(x,n,0,funct)
      call davidn( x,n,0,xx,yy,nn,kkx,kky,kkc,kkt,mm,iopt,lf,
     & t,f,xm,lu,ifg )
c----------------------------------
   30 continue
c
      x(1) = x(1)**2
      x(2) = x(2)**2
      if(kkx.ne.0) x(k2+1)=x(k2+1)**2
      if(kky.ne.0) x(kkx+k2+1)=x(kkx+k2+1)**2
      if(kkt.ne.0) x(kkx+kky+kkc+k2+1)=x(kkx+kky+kkc+k2+1)**2
c
c     polynomial(trend) coefficient's transformation
c
      if(kkt.eq.1) go to 52
      do 53 k=2,kkt
      x(k+kkx+kky+kkc+k2)=x(k+kkx+kky+kkc+k2)/t**(k-1)
   53 continue
   52 continue
c
cc      write(6,1040) f,(x(i),i=1,n)
cc      write(6,1080) x(1),x(2)
cc      if(kkx.ne.0) write(6,1100) (x(k+k2),k=1,kkx)
cc      if(kky.ne.0) write(6,1110) (x(k+kkx+k2),k=1,kky)
cc      if(kkc.ne.0) write(6,1120) (x(k+kkx+kky+k2),k=1,kkc)
cc      if(kkt.ne.0) write(6,1130) (x(k+kkx+kky+kkc+k2),k=1,kkt)
      aic=f+n
      if(kkx.eq.0) aic=aic-1
      if(kky.eq.0) aic=aic-1
      if(kkx.ne.0.and.iopt.eq.1) aic=aic-1
 1005 format(e20.10,3i10)
cc      write(6,1001) aic
 1001 format(/1h ,'aic/2 =',e20.10)
   50 continue
   20 continue
      return
  100 aic=-nn*log(nn/t)+nn+1
cc      write(6,1001) aic
      return
 1000 format(3i10,2f10.6)
 1010 format(8f10.4)
 1020 format(1h //1h ,' n =',i5//1h ,' initial estimates (x(i),i=1,n)'/
     &            (1h ,5d13.5))
 1040 format(/1h ,'neg max lklhd=',1 e16.7
     1      /1h ,'max lklhd est.=',4d14.6/('                ',4d14.6))
 1050 format(4d20.13)
 1060 format(e25.15)
 1070 format(/1h ,'  c = ',e25.15)
 1071 format(/1h ,'  d = ',e25.15)
 1080 format(////1h ,' final outputs '/1h ,' c , d       ',2d14.6)
 1100 format(1h ,'ax(i),i=1,kkx',4d14.6/(1h ,13x,4d14.6))
 1110 format(1h ,'ay(i),i=1,kky',4d14.6/(1h ,13x,4d14.6))
 1120 format(1h ,'ac(i),i=1,kkc',4d14.6/(1h ,13x,4d14.6))
 1130 format(1h ,'at(i),i=1,kkt',4d14.6/(1h ,13x,4d14.6))
      end
cc      subroutine  linear( x,h,ram,ee,k,ig,funct )
      subroutine  linear( x,h,ram,f,ee,k,ig,xx,yy,t,nn,mm,iopt,ff,kkx,
     & kky,kkc,kkt,lf,lu,jfg )
c
c     this subroutine performs the linear search along the direction
c     specified by the vector h
c
c     this subroutine is copied from timsac 78.
c     ---------------------------------------------------------------
c     the following subroutine is directly called by this subroutine:
c           funct
c     ---------------------------------------------------------------
c
c     inputs:
c        x:       vector of position
c        h:       search direction
c        k:       dimension of vector x
c
c     outputs:
c        ram:     optimal step width
c        e2:      minimum function value
c        ig:      error code
c
      implicit  real  *8 ( a-h,o-z )
      integer  return,sub
cc      dimension  x(1) , h(1) , x1(50)
cc      dimension  g(50)
        common     / ccc /  isw , ipr
cc      external funct
      dimension  x(k) , h(k) , x1(k)
      dimension  g(k)
      dimension  xx(nn),yy(nn)
      dimension  lf(51,51)
c
      isw = 1
      ipr = 7
      if( ram .le. 1.0d-30 )  ram = 0.01d0
      const2 = 1.0d-60
      hnorm = 0.d0
      do 10  i=1,k
   10 hnorm = hnorm + h(i)**2
      hnorm = dsqrt( hnorm )
c
      ram2 = ram
      e1 =ee
      ram1 = 0.d0
c
      do 20  i=1,k
   20 x1(i) = x(i) + ram2*h(i)
cc      call  funct( k,x1,e2,g,ig )
      call funct(k,x1,e2,g,ig,xx,yy,t,nn,mm,iopt,ff,kkx,kky,
     & kkc,kkt,lf)
cc      if(ipr.ge.7)  write(6,2)  ram2,e2
      if( jfg.eq.1 )  write(lu,2) ram2,e2
c
      if( ig .eq. 1 )  go to  50
      if( e2 .gt. e1 )  go to 50
   30 ram3 = ram2*2.d0
      do 40  i=1,k
   40 x1(i) = x(i) + ram3*h(i)
cc      call  funct( k,x1,e3,g,ig )
      call funct(k,x1,e3,g,ig,xx,yy,t,nn,mm,iopt,ff,kkx,kky,
     & kkc,kkt,lf)
      if( ig.eq.1 )  go to  500
cc      if( ipr.ge.7 )  write(6,3)  ram3,e3
      if( jfg.eq.1 )  write(lu,3)  ram3,e3
      if( e3 .gt. e2 )  go to 70
      ram1 = ram2
      ram2 = ram3
      e1 = e2
      e2 = e3
      go to 30
c
   50 ram3 = ram2
      e3 = e2
      ram2 = ram3*0.1d0
      if( ram2*hnorm .lt. const2 )  go to  400
      do 60  i=1,k
   60 x1(i) = x(i) + ram2*h(i)
cc      call  funct( k,x1,e2,g,ig )
      call funct(k,x1,e2,g,ig,xx,yy,t,nn,mm,iopt,ff,kkx,kky,
     & kkc,kkt,lf)
cc      if(ipr.ge.7)  write(6,4)  ram2,e2
      if( jfg.eq.1 )  write(lu,4)  ram2,e2
      if( e2.gt.e1 )  go to 50
c
cc   70 assign 80 to return
   70 return = 80
      go to 200
c
   80 do 90  i=1,k
   90 x1(i) = x(i) + ram*h(i)
cc      call  funct( k,x1,ee,g,ig )
      call funct(k,x1,ee,g,ig,xx,yy,t,nn,mm,iopt,ff,kkx,kky,
     & kkc,kkt,lf)
cc      if(ipr.ge.7)  write(6,5)  ram,ee
      if(jfg.eq.1)  write(lu,5)  ram,ee
c
      ifg = 0
cc      assign  300 to  sub
cc      assign 200 to sub
cc   95 assign 130 to return
      sub = 200
   95 return = 130
      if( ram .gt. ram2 )  go to 110
      if( ee .ge. e2 )  go to 100
      ram3 = ram2
      ram2 = ram
      e3 =e2
      e2 =ee
cc      go to  sub,( 200,300 )
      if( sub .eq. 200 ) go to 200
      if( sub .eq. 300 ) go to 300
c
  100 ram1 = ram
      e1 = ee
cc      go to  sub,( 200,300 )
      if( sub .eq. 200 ) go to 200
      if( sub .eq. 300 ) go to 300
c
  110 if( ee .le. e2 )  go to 120
      ram3 = ram
      e3 = ee
cc      go to  sub,( 200,300 )
      if( sub .eq. 200 ) go to 200
      if( sub .eq. 300 ) go to 300
c
  120 ram1 = ram2
      ram2 = ram
      e1 = e2
      e2 = ee
cc      go to  sub,( 200,300 )
      if( sub .eq. 200 ) go to 200
      if( sub .eq. 300 ) go to 300
c
  130 do 140  i=1,k
  140 x1(i) = x(i) + ram*h(i)
cc      call  funct( k,x1,ee,g,ig )
      call funct(k,x1,ee,g,ig,xx,yy,t,nn,mm,iopt,ff,kkx,kky,
     & kkc,kkt,lf)
cc      if( ipr.ge.7 )  write(6,6)  ram,ee
      if( jfg.eq.1 )  write(lu,6)  ram,ee
cc      assign 200 to sub
      sub = 200
      ifg = ifg+1
      ifg = 0
      if( ifg .eq. 1 )  go to 95
c
      if( e2 .lt. ee )  ram = ram2
      return
c
c      -------  internal subroutine sub1  -------
  200 a1 = (ram3-ram2)*e1
      a2 = (ram1-ram3)*e2
      a3 = (ram2-ram1)*e3
      b2 = (a1+a2+a3)*2.d0
      b1 = a1*(ram3+ram2) + a2*(ram1+ram3) + a3*(ram2+ram1)
      if( b2 .eq. 0.d0 )  go to 210
      ram = b1 /b2
cc      go to return ,( 80,130 )
      if( return .eq. 80 ) go to 80
      if( return .eq. 130 ) go to 130
c
  210 ig = 1
      ram = ram2
      return
c
c      -------  internal subroutine sub2  -------
c
  300 if( ram3-ram2 .gt. ram2-ram1 )  go to 310
      ram = (ram1+ram2)*0.5d0
cc      go to return ,( 80,130 )
      if( return .eq. 80 ) go to 80
      if( return .eq. 130 ) go to 130
c
  310 ram = (ram2+ram3)*0.5d0
cc      go to return ,( 80,130 )
      if( return .eq. 80 ) go to 80
      if( return .eq. 130 ) go to 130
c ------------------------------------------------------------
c
  400 ram = 0.d0
      return
c ------------------------------------------------------------
c
  500 ram = (ram2+ram3)*0.5d0
  510 do 520  i=1,k
  520 x1(i) = x(i) + ram*h(i)
cc      call  funct( k,x1,e3,g,ig )
      call funct(k,x1,e3,g,ig,xx,yy,t,nn,mm,iopt,ff,kkx,kky,
     & kkc,kkt,lf)
cc      if( ipr.ge.7 )  write(6,7)  ram,e3
      if( jfg.eq.1 )  write(lu,7)  ram3,e3
      if( ig.eq.1 )  go to 540
      if( e3.gt.e2 )  go to 530
      ram1 = ram2
      ram2 = ram
      e1 = e2
      e2 = e3
      go to 500
c
  530 ram3 = ram
      go to 70
c
  540 ram = (ram2+ram)*0.5d0
      go to 510
c
c ------------------------------------------------------------
    1 format( 1h ,'lambda =',d18.10, 10x,'e1 =',d25.17 )
    2 format( 1h ,'lambda =',d18.10, 10x,'e2 =',d25.17 )
    3 format( 1h ,'lambda =',d18.10, 10x,'e3 =',d25.17 )
    4 format( 1h ,'lambda =',d18.10, 10x,'e4 =',d25.17 )
    5 format( 1h ,'lambda =',d18.10, 10x,'e5 =',d25.17 )
    6 format( 1h ,'lambda =',d18.10, 10x,'e6 =',d25.17 )
    7 format( 1h ,'lambda =',d18.10, 10x,'e7 =',d25.17 )
      e n d
cc      subroutine  davidn( x,n,ihes,funct )
      subroutine  davidn( x,n,ihes,xx,yy,nn,kkx,kky,kkc,kkt,mm,iopt,lf,
     & t,f,xm,lu,ifg )
c
c     minimization by davidon-fletcher-powell procedure
c     this subroutine is copied from timsac 78.
c
c     ---------------------------------------------------------------
c     the following subroutines are directly called by this subroutine
c           funct
c           hesian
c           linear
c     ---------------------------------------------------------------
c     inputs:
c        x:       vector of initial values
c        k:       dimension of the vector x
c        ihes:    =0   inverse of hessian matrix is not available
c                 =1   inverse of hessian matrix is available
c
c     output:
c        x:       vector of minimizing solution
c
      implicit  real * 8  ( a-h , o-z )
cc      dimension  x(50) , dx(50) , g(50) , g0(50) , y(50)
cc      dimension  h(50,50) , wrk(50) , s(50)
      dimension  x(n) , dx(n) , g(n) , g0(n) , y(n)
      dimension  h(n,n) , wrk(n) , s(n)
      dimension  xx(nn),yy(nn)
cc      dimension  r(31,31)
      common     / ccc /  isw,ipr
cc      common     / ddd /  r , f , aic , sd
cc      external funct
      data  tau1 , tau2  /  1.0d-5 , 1.0d-5  /
      data  eps1 , eps2  / 1.0d-5 , 1.0d-5  /
      ramda = 0.5d0
      const1 = 1.0d-70
c
c          initial estimate of inverse of hessian
c
      do  20   i=1,n
      do  10   j=1,n
   10 h(i,j) = 0.0d00
      s(i) = 0.0d00
      dx(i) = 0.0d00
   20 h(i,i) = 1.0d00
cc      isw = 0
c
cc      call  funct( n,x,xm,g,ig )
      call funct(n,x,xm,g,ig,xx,yy,t,nn,mm,iopt,f,kkx,kky,
     & kkc,kkt,lf)
c
cc      write( 6,340 )     xm
      if( ifg.eq.1 )  write( lu,340 )     xm
c
c          inverse of hessian computation (if available)
c
c     if( ihes .eq. 1 )   call  hesian( x,n,h )
c
      icc = 0
c      iteration
11110 continue
      icc = icc + 1
      do  11111   ic=1,n
      if( ic .eq. 1 .and. icc .eq. 1 )     go to 120
c
      do  40   i=1,n
   40 y(i) = g(i) - g0(i)
      do  60   i=1,n
      sum = 0.0d00
      do  50   j=1,n
   50 sum = sum + y(j) * h(i,j)
   60 wrk(i) = sum
      s1 = 0.0d00
      s2 = 0.0d00
      do  70   i=1,n
      s1 = s1 + wrk(i) * y(i)
   70 s2 = s2 + dx(i) * y(i)
      if( s1.le.const1 .or. s2.le.const1 )  go to 900
      if( s1 .le. s2 )     go to 100
c
c          update the inverse of hessian matrix
c
c               ---  davidon-fletcher-powell type correction  ---
c
      do  90   i=1,n
      do  90   j=i,n
      h(i,j) = h(i,j) + dx(i)*dx(j)/s2 - wrk(i)*wrk(j)/s1
   90 h(j,i) = h(i,j)
      go to  120
c
c               ---  fletcher type correction  ---
c
  100 continue
      stem = s1 / s2 + 1.0d00
      do  110   i=1,n
      do  110   j=i,n
      h(i,j) = h(i,j)- (dx(i)*wrk(j)+wrk(i)*dx(j)-dx(i)*dx(j)*stem)/s2
  110 h(j,i) = h(i,j)
c
c
c
  120 continue
      ss = 0.0d00
      do  150   i=1,n
      sum = 0.0d00
      do  140   j=1,n
  140 sum = sum + h(i,j)*g(j)
      ss = ss + sum * sum
  150 s(i) = -sum
c
c
      s1 = 0.0d00
      s2 = 0.0d00
      do  170   i=1,n
      s1 = s1 + s(i)*g(i)
  170 s2 = s2 + g(i)*g(i)
      ds2 = dsqrt(s2)
      gtem = dabs(s1) / ds2
c     write(6,610)gtem,ds2
      if( gtem .le. tau1  .and.  ds2 .le. tau2 )     go to  900
      if( s1 .lt. 0.0d00 )     go to  200
      do  190   i=1,n
      do  180   j=1,n
  180 h(i,j) = 0.0d00
      h(i,i) = 1.0d00
  190 s(i) = -s(i)
  200 continue
c
      ed = xm
c
c          linear  search
c
cc      call  linear( x,s,ramda,ed,n,ig,funct )
      call  linear( x,s,ramda,f,ed,n,ig,xx,yy,t,nn,mm,iopt,f,kkx,
     & kky,kkc,kkt,lf,lu,ifg )
c
cc      write( 6,330 )     ramda , f
      if( ifg.eq.1 )  write( lu,330 )     ramda , f
c
      s1 = 0.0d00
      do  210   i=1,n
      dx(i) = s(i) * ramda
      s1 = s1 + dx(i) * dx(i)
      g0(i) = g(i)
  210 x(i) = x(i) + dx(i)
      xmb = xm
cc      isw = 0
c
cc      call  funct( n,x,xm,g,ig )
      call funct(n,x,xm,g,ig,xx,yy,t,nn,mm,iopt,f,kkx,kky,
     & kkc,kkt,lf)
c
      s2 = 0.d0
      do  220     i=1,n
  220 s2 = s2 + g(i)*g(i)
      if( dsqrt(s2) .gt. tau2 )   go to  11111
      if( xmb/xm-1.d0 .lt. eps1  .and.  dsqrt(s1) .lt. eps2 )  go to 900
11111 continue
      if( icc .ge. 5 )     go to 900
      go to 11110
  900 continue
cc      write( 6,600 )
cc      write( 6,610 )     (x(i),i=1,n)
cc      write( 6,601 )
cc      write( 6,610 )     (g(i),i=1,n)
      if( ifg.eq.1 )  then
         write( lu,600 )
         write( lu,610 )     (x(i),i=1,n)
         write( lu,601 )
         write( lu,610 )     (g(i),i=1,n)
      end if
      return
  330 format( 1h ,'lambda =',d15.7,5x,'log likelihood =',d18.10)
  340 format( 1h ,28x,'log-likelihood =',d18.10)
  600 format( /1h ,'-----  x  -----' )
  601 format( /1h ,'***  gradient  ***' )
  610 format( 1h ,5d13.5 )
      end
cc      subroutine funct(kk2,b,f,h,ifg)
      subroutine funct(kk2,b,f,h,ifg,xx,yy,t,nn,mm,iopt,ff,
     & kkx,kky,kkc,kkt,lf)
c
c     likelihood function of l-g hawkes' point process
c
      implicit real * 8 (a-h,o-z)
cc      common/xyod/xx(10000),yy(10000)
cc      common t,nn,mm,iopt
cc      common /ddd/r,ff,aic,sd
cc      common /kkxy/kkx,kky
cc      common /ct/kkc,kkt
cc      common /comf/lf(51,51)
cc      dimension r(31,31)
cc      dimension b(50),h(50)
cc      dimension a(50),ei(50),si(50),s(50),ds(50),g(50),gg(50)
cc      dimension ay(50),eyi(50),syi(50),sy(50),dsy(50),gy(50),ggy(50)
cc      dimension ac(50),at(50),gc(50),gt(50),ggc(50),ggt(50)
cc      dimension dei(50)
      dimension xx(nn),yy(nn)
      dimension lf(51,51)
      dimension b(kk2),h(kk2)
      dimension a(kkx),ei(kkx+1),si(kkx+1),s(kkx+1),ds(kkx+1)
      dimension g(kkx),gg(kkx),ay(kky),eyi(kky+1),syi(kky+1)
      dimension sy(kky+1),dsy(kky+1),gy(kky),ggy(kky)
      dimension ac(kkc),at(kkt),gc(kkc),gt(kkt),ggc(kkc),ggt(kkt)
      dimension dei(kkx+kky+1)
      pi=3.14159265358979d0
      t0=365.25d0
      if(kk2.eq.kkx+kky+kkc+kkt+2) go to 5
cc      write(6,4) kk2,kkx,kky,kkc,kkt
    4 format(1h ,'n or kkx or kky kkc kkt error',5i10)
      ifg=-1
      stop
    5 ifg=0
      ix=1
      wop=b(kkx+kky+kkc+3)**2
      c=b(1)**2
      d=b(2)**2
c
      if(kkx.eq.0) go to 201
      do 200 k=1,kkx
  200 a(k)=b(k+2)
      a(1)=b(3)**2
  201 continue
c
      if(kky.eq.0) go to 202
      do 205 k=1,kky
  205 ay(k)=b(k+kkx+2)
      ay(1)=b(kkx+3)**2
  202 continue
      if(kky.eq.0) ay(1)=0.0
c
      if(kkc.eq.0) go to 203
      do 206 k=1,kkc
  206 ac(k)=b(k+kkx+kky+2)
  203 continue
c
      if(kkt.eq.1) go to 204
      do 207 k=2,kkt
  207 at(k)=b(k+kkx+kky+kkc+2)
  204 continue
c
      f1=0.0
c
      if(kkx.eq.0) go to 11
      do 10 k=1,kkx
   10 gg(k)=0.0
   11 continue
c
      if(kky.eq.0) go to 12
      do 15 k=1,kky
   15 ggy(k)=0.0
   12 continue
c
      if(kkc.eq.0) go to 13
      do 16 k=1,kkc
   16 ggc(k)=0.0
   13 continue
c
      if(kkt.eq.0) go to 14
      do 17 k=1,kkt
   17 ggt(k)=0.0
   14 continue
c
      ggwop=0.0
      ggoc=0.0
      ggod=0.0
      do 30 k=1,kkx+1
   30 ei(k)=0.0
      do 35 k=1,kky+1
   35 eyi(k)=0.0
c
c
c
      do 20 i=1,nn
      if(i.eq.1) go to 310
      dxxi=xx(i)-xx(i-1)
      if(c.le.0.0) go to 50
      ecdxxi=0.0d00
      if(-c*dxxi.ge.-100.0) ecdxxi=exp(-c*dxxi)
      do 25 j=1,kkx+1
      jj=kkx+2-j
      eii=ei(1)+1.d0
      if(jj.eq.1) go to 27
      do 26 k=2,jj
      if(lf(jj,k).gt.1) eii=lf(jj,k)*ei(k)+dxxi*eii
      if(lf(jj,k).eq.1) eii=ei(k)+dxxi*eii
   26 continue
   27 continue
      ei(jj)=ecdxxi*eii
   25 continue
  310 continue
c
      if(d.le.0.0) go to 50
      if(i.eq.1) dxxi=xx(i)
      eddxxi=0.0d00
      if(-d*dxxi.ge.-100.0) eddxxi=exp(-d*dxxi)
      do 45 j=1,kkx+1
      dei(j)=0.0
   45 continue
c
      if(mm.eq.0) go to 320
      do 40 j=ix,mm
      if(yy(j).gt.xx(i)) go to 300
      dxyij=xx(i)-yy(j)
      if(-d*dxyij.lt.-100.0) go to 40
      deij=exp(-d*dxyij)
      dei(1)=dei(1)+deij
      if(kky.eq.0) go to 40
      do 41 k=2,kky+1
      deij=dxyij*deij
      dei(k)=dei(k)+deij
   41 continue
   40 continue
  300 ix=j
  320 continue
c
c
      do 125 j=1,kky+1
      jj=kky+2-j
      eyii=eyi(1)
      if(jj.eq.1) go to 127
      do 126 k=2,jj
      if(lf(jj,k).gt.1) eyii=lf(jj,k)*eyi(k)+dxxi*eyii
      if(lf(jj,k).eq.1) eyii=eyi(k)+dxxi*eyii
  126 continue
  127 continue
      eyi(jj)=eddxxi*eyii+dei(jj)
  125 continue
c
      rmdi1=0.0
      rmdi2=0.0
      if(kkx.eq.0) go to 81
      do 80 k=1,kkx
      rmdi1=rmdi1+a(k)*ei(k)
      rmdi2=rmdi2-a(k)*ei(k+1)
c     write(6,*) k,rmdi1,a(k),ei(k)
   80 continue
   81 continue
c
      rmdyi1=0.0
      rmdyi2=0.0
      rmdci1=0.0
      if(kky.eq.0) go to 82
      do 60 k=1,kky
      rmdyi1=rmdyi1+ay(k)*eyi(k)
      rmdyi2=rmdyi2-ay(k)*eyi(k+1)
   60 continue
   82 continue
c
      if(kkc.eq.0) go to 83
      do 400 k=1,kkc
      if(mod(k,2).eq.1) rmdci1=rmdci1+ac(k)*cos(2*((k+1)/2)*pi*xx(i)/t0)
      if(mod(k,2).eq.0) rmdci1=rmdci1+ac(k)*sin(2*((k+1)/2)*pi*xx(i)/t0)
  400 continue
   83 continue
c
      rmdti1=0.0
      if(kkt.eq.1) go to 84
      do 410 k=2,kkt
  410 rmdti1=rmdti1+at(k)*(xx(i)/t)**(k-1)
   84 continue
c
c
      ramdai=wop+rmdi1+rmdyi1+rmdci1+rmdti1
      if(ramdai.le.0.0) go to 50
      f1=f1+log(ramdai)
      ggwop=ggwop+1.0/ramdai
      ggoc=ggoc+rmdi2/ramdai
      ggod=ggod+rmdyi2/ramdai
c
      if(kkx.eq.0) go to 71
      do 70 k=1,kkx
      gg(k)=gg(k)+ei(k)/ramdai
   70 continue
   71 continue
c
      if(kky.eq.0) go to 91
      do 90 k=1,kky
      ggy(k)=ggy(k)+eyi(k)/ramdai
   90 continue
   91 continue
c
      if(kkc.eq.0) go to 391
      do 390 k=1,kkc
      if(mod(k,2).eq.1)ggc(k)=ggc(k)+cos(2*pi*((k+1)/2)*xx(i)/t0)/ramdai
      if(mod(k,2).eq.0)ggc(k)=ggc(k)+sin(2*pi*((k+1)/2)*xx(i)/t0)/ramdai
  390 continue
  391 continue
c
      if(kkt.eq.1) go to 381
      do 380 k=2,kkt
      ggt(k)=ggt(k)+(xx(i)/t)**(k-1)/ramdai
  380 continue
  381 continue
c
c
   20 continue
c
c
c
      do 100 k=1,kkx+1
      s(k)=0.0
      ds(k)=0.0
  100 continue
c
      do 110 i=1,nn
      fatxxi=-c*(t-xx(nn-i+1))
      if(fatxxi.gt.70.0) go to 50
      si(1)=1.0d00/c
      if(fatxxi.ge.-100.0) si(1)=(1.0d00-exp(fatxxi))/c
      s(1)=s(1)+si(1)
      if(kkx.eq.0) go to 121
      do 120 k=2,kkx+1
      si(k)=(t-xx(nn-i+1))**(k-1)*(si(1)-1.0/c)+(k-1)/c*si(k-1)
      s(k)=s(k)+si(k)
  120 continue
  121 continue
  110 continue
c
      if(kkx.eq.0) go to 141
      do 140 k=1,kkx
      ds(k)=-s(k+1)
  140 continue
c
  141 continue
      sasum=0.0
      dsasum=0.0
      if(kkx.eq.0) go to 151
      do 150 k=1,kkx
      sasum=sasum+a(k)*s(k)
      dsasum=dsasum+a(k)*ds(k)
  150 continue
  151 continue
c
c
      do 180 k=1,kky+1
      sy(k)=0.0
      dsy(k)=0.0
  180 continue
c
      if(mm.eq.0) go to 195
      do 190 i=1,mm
      fatyyi=-d*(t-yy(mm-i+1))
      if(fatyyi.gt.70.0) go to 50
      syi(1)=1.0d00/d
      if(fatyyi.ge.-100.0) syi(1)=(1.0d00-exp(fatyyi))/d
      sy(1)=sy(1)+syi(1)
      if(kky.eq.0) go to 211
      do 210 k=2,kky+1
      syi(k)=(t-yy(mm-i+1))**(k-1)*(syi(1)-1.0/d)+(k-1)/d*syi(k-1)
      sy(k)=sy(k)+syi(k)
  210 continue
  211 continue
  190 continue
  195 continue
c
      if(kky.eq.0) go to 221
      do 220 k=1,kky
      dsy(k)=-sy(k+1)
  220 continue
  221 continue
c
      sysum=0.0
      dsysum=0.0
      if(kky.eq.0) go to 235
      do 230 k=1,kky
      sysum=sysum+ay(k)*sy(k)
      dsysum=dsysum+ay(k)*dsy(k)
  230 continue
  235 continue
c
      scsum=0.0
      if(kkc.eq.0) go to 375
      do 370 k=1,kkc
      if(mod(k,2).eq.1) scsum=scsum+ac(k)*sin(2*((k+1)/2)*pi*t/t0)
     &                                    *t0/pi/(2*((k+1)/2))
      if(mod(k,2).eq.0)scsum=scsum-ac(k)*(cos(2*((k+1)/2)*pi*t/t0)-1.d0)
     &                                    *t0/pi/(2*((k+1)/2))
  370 continue
  375 continue
c
      stsum=0.0
      if(kkt.eq.1) go to 365
      do 360 k=2,kkt
  360 stsum=stsum+at(k)*t*(t/t)**(k)/(k)
  365 continue
c
      f=f1-wop*t-sasum-sysum-scsum-stsum
      gwop=ggwop-t
      goc=ggoc-dsasum
      god=ggod-dsysum
c
      if(kkx.eq.0) go to 251
      do 250 k=1,kkx
      g(k)=gg(k)-s(k)
  250 continue
  251 continue
c
      if(kky.eq.0) go to 254
      do 255 k=1,kky
      gy(k)=ggy(k)-sy(k)
  255 continue
  254 continue
c
      if(kkc.eq.0) go to 253
      do 256 k=1,kkc
      if(mod(k,2).eq.1) gc(k)=ggc(k)-sin(2*pi*((k+1)/2)*t/t0)
     &                               *t0/(2*pi)/((k+1)/2)
      if(mod(k,2).eq.0) gc(k)=ggc(k)+(cos(2*pi*((k+1)/2)*t/t0)-1.d0)
     &                               *t0/(2*pi)/((k+1)/2)
  256 continue
  253 continue
c
      if(kkt.eq.1) go to 355
      do 350 k=2,kkt
  350 gt(k)=ggt(k)-t*(t/t)**(k)/(k)
  355 continue
c
c
      f=-f
      h(kkx+kky+kkc+3)=-gwop*2.d0*b(kkx+kky+kkc+3)
      h(1)=-goc*2.d0*b(1)
      h(2)=-god*2.d0*b(2)
      if(iopt.eq.1) h(1)=0.0
c
      if(kkx.eq.0) go to 261
      do 260 k=1,kkx
  260 h(k+2)=-g(k)
      h(3)=-g(1)*2.0d00*b(3)
  261 continue
c
      if(kky.eq.0) go to 264
      do 265 k=1,kky
  265 h(k+kkx+2)=-gy(k)
      h(kkx+3)=-gy(1)*2.0d00*b(kkx+3)
  264 continue
c
      if(kkc.eq.0) go to 263
      do 266 k=1,kkc
  266 h(k+kkx+kky+2)=-gc(k)
  263 continue
c
      if(kkt.eq.1) go to 268
      do 267 k=2,kkt
  267 h(k+kkx+kky+kkc+2)=-gt(k)
  268 continue
c
      ff=f
    3 format(1h ,110x,d18.10)
    1 format(1h ,7d18.10)
      return
c
   50 continue
      ifg=1
      f=1.0d30
      return
c
      end
