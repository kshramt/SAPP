cc      program linsim
      subroutine linsimf(kxx,kxy,kxz,t,c,d,axx,axy,axz,yy,mm,ptxmax,
     & kmax,xx,i1,j1,err)
c
      include 'sapp.h'
c
c     this program performs similation of a self-exciting point process
c  which is also stimulated by another given point process data; non
c  stationary poisson trend is also included in the intensity
c  function.
c
c     structure of the program
c
c          linsim
c             i---input
c             i---initl
c             i---simxx----comfac
c             i     i------duf
c             i     i------unifor
c             i     +------fx--------ptrend
c             +---output
c
c     this program is designed by y. ogata, and programmed by y. ogata
c  and k. katsura, inst. statist. math., tokyo. (31/01/85)
c
c     references
c
c  y. ogata (1981). "on lewis' simulation method for point processes."
c     ieee information theory, vol. it-27, pp. 23-31.
c  y. ogata & h. akaike (1982). "on linear intensity models for mixed
c     doubly stochastic poisson and self-exciting point processes." j.
c     royal statist. soc. b, vol. 44, pp. 102-107.
c  y. ogata, h. akaike and k. katsura (1982). "the application of linear
c     intensity models to the investigatio of causal relations between a
c     point process and another stochastic process."  ann. inst.
c     statist math., vol. 34. pp. 373-387.
c  y. ogata and k. katsura (1985). "point process model with linearly
c     parametrized intensity for the application to earthquake data."
c     essays in time series and allied processes (festscrift for prof.
c     e. j. hannan), j. gani and m. b. priestley eds., j.appl. probab.
c     vol. 23a, to appear.
c
      implicit real*8(a-h,o-z)
cc      dimension xx(10000),yy(10000)
cc      dimension axx(100),axy(100),axz(100)
cc      dimension ei(100),ej(100),fi(100),fj(100),lf(51,51)
      dimension xx(2*mm),yy(mm+1)
      dimension axx(kxx),axy(kxy),axz(kxz)
c
cc      call input(kxx,kxy,kxz,t,c,d,axx,axy,axz,yy,mm,ptxmax,kmax)
c
      call initl(kxx,axx,c,fxxmax)
      call initl(kxy,axy,d,fxymax)
c
      call simxx(kxx,kxy,kxz,t,c,d,axx,axy,axz,fxxmax,fxymax,
cc     &           xx,yy,ei,ej,fi,fj,kmax,ptxmax,lf,i1,j1)
     &           xx,yy,kmax,ptxmax,i1,j1,err)
c
cc      call output(xx,yy,i1,j1,t)
c
cc      stop
      return
      end
      subroutine simxx(kxx,kxy,kxz,t,c,d,axx,axy,axz,fxxmax,fxymax,
cc     &                 xx,yy,ei,ej,fi,fj,kmax,ptxmax,lf,i1,j1)
     &                 xx,yy,kmax,ptxmax,i1,j1,err)
      implicit real*8(a-h,o-z)
cc      dimension axx(1),axy(1),axz(1),xx(1),yy(1),ei(1),ej(1),fi(1),fj(1)
      dimension axx(1),axy(1),axz(1),xx(1),yy(1)
      dimension ei(kmax),ej(kmax),fi(kmax),fj(kmax)
      dimension lf(51,51)
      real*4r
c
      err=0.0
      ir=584287
c
      do 100 k=1,kmax
      ei(k)=0.0
      ej(k)=0.0
      fi(k)=0.0
      fj(k)=0.0
  100 continue
c-------------------------
      call comfac(lf)
c-------------------------
      x=0.0
      xx(1)=0.0
      uity=ptxmax
      i=0
      j=0
      ij=0
   30 continue
c----------------------------------------------------------------
      call duf(i,j,x,duity,xx,yy,axx,axy,kxx,kxy,c,d,ei,ej,fi,fj,
     &         ptxmax)
c----------------------------------------------------------------
c
      if(duity.lt.uity) uity=duity
   20 continue
c--------------------
cc  130 call unifor(r)
  130 call unifor(r,ir)
c--------------------
   40 continue
      e=-alog(r)/uity
      x=x+e
      if(x.gt.t) go to 80
    6 format(2i3,7f10.5)
      if(kxy.eq.0) go to 110
      if(x.le.yy(j+1)) go to 110
      x=yy(j+1)
c-----------------------------------------------------------------
      call fx(i,j,yy(j+1),dmx,axx,axy,axz,kxx,kxy,kxz,c,d,lf,ei,ej,
     &        fi,fj,xx,yy)
c-----------------------------------------------------------------
      do 120 k=1,kxy
  120 ej(k)=fj(k)
      j=j+1
      uity=uity+fxymax
      go to 130
c
  110 continue
c
c------------------------------------------------------------
      call fx(i,j,x,xity,axx,axy,axz,kxx,kxy,kxz,c,d,lf,ei,ej,
     &        fi,fj,xx,yy)
c------------------------------------------------------------
      probx=xity/uity
      prob=xity/uity
      if(prob.le.1.0d00)go to 50
cc      write(6,1) prob
      err=prob
      stop
   50 continue
   60 continue
c-----------------------
cc      call unifor(r)
      call unifor(r,ir)
c-----------------------
      ij=ij+1
      if(probx.le.r)go to 30
      i=i+1
      xx(i)=x
      do 10 k=1,kxx
      ei(k)=fi(k)
   10 continue
      uity=uity+fxxmax
      go to 20
   80 continue
      i1=i
      j1=j
    1 format(1h ,'warning: is ptmax correct?  prob=',f10.5)
      return
      end
      subroutine duf(i,j,x,duity,xx,yy,axx,axy,kxx,kxy,c,d,ei,ej,fi,
     &               fj,ptxmax)
c
c     decreasing process except jumps, which is always greater than
c     the intensity process in subroutine fx.
c
      implicit real * 8 (a-h,o-z)
cc      dimension bxx(100),bxy(100)
      dimension bxx(kxx),bxy(kxy)
      dimension axx(1),axy(1),ei(1),ej(1),fi(1),fj(1)
      dimension xx(1),yy(1)
      ixf=1
      iyf=1
      cxp=0.0
      cyp=0.0
      do 30 ix=1,kxx
      bxx(ix)=max(axx(ix),0.d0)
      cx=bxx(ix)/(c/2)**(ix-1) *ixf
      cxp=max(cx,cxp)
      ixf=ixf*ix
   30 continue
      do 50 iy=1,kxy
      bxy(iy)=max(axy(iy),0.d0)
      cy=bxy(iy)/(d/2)**(iy-1) *iyf
      cyp=max(cy,cyp)
      iyf=iyf*iy
   50 continue
      if(i.eq.0) go to 10
      dxxi=x-xx(i)
      ecdxxi=0.0d00
      if(-c/2*dxxi.ge.-20.0) ecdxxi=exp(-c/2*dxxi)
      fi(1)=ecdxxi*(ei(1)+1.d0)
   10 continue
      if(j.eq.0) go to 20
      dyyj=x-yy(j)
      ecdyyj=0.0d00
      if(-d/2*dyyj.ge.-20.0) ecdyyj=exp(-d/2*dyyj)
      fj(1)=ecdyyj*(ej(1)+1.d0)
   20 continue
      xity=ptxmax+cxp*fi(1)+cyp*fj(1)
      duity=xity
      return
      end
      subroutine fx(i,j,x,xity,axx,axy,axz,kxx,kxy,kxz,c,d,lf,ei,
     &              ej,fi,fj,xx,yy)
c
c     intensity processes
c
      implicit real * 8 (a-h,o-z)
      dimension axx(1),axy(1),axz(1),ei(1),ej(1),fi(1),fj(1),lf(51,51)
      dimension xx(1),yy(1)
      if(i.eq.0) go to 30
      dxxi=x-xx(i)
      ecdxxi=0.0d00
      if(-c*dxxi.ge.-20.0) ecdxxi=exp(-c*dxxi)
      do 70 ix=1,kxx
      ff=ei(1)*lf(ix,1)+1.d0
      if(ix.eq.1) go to 85
      do 80 jx=2,ix
      ff=ff*dxxi+ei(jx)*lf(ix,jx)
   80 continue
   85 continue
      fi(ix)=ecdxxi*ff
   70 continue
   30 continue
      if(j.eq.0) go to 40
      dyyj=x-yy(j)
      ecdyyj=0.0d00
      if(-d*dyyj.ge.-20.0) ecdyyj=exp(-d*dyyj)
      do 90 ix=1,kxy
      ff=ej(1)*lf(ix,1)+1.d0
      if(ix.eq.1) go to 105
      do 100 jx=2,ix
      ff=ff*dyyj+ej(jx)*lf(ix,jx)
  100 continue
  105 continue
      fj(ix)=ecdyyj*ff
   90 continue
   40 continue
      call ptrend(x,ptx,axz,kxz)
      xity=ptx
      if(i.eq.1) go to 15
      do 10 k=1,kxx
      xity=xity+axx(k)*fi(k)
   10 continue
   15 continue
      if(j.eq.1) go to 25
      do 20 k=1,kxy
      xity=xity+axy(k)*fj(k)
   20 continue
   25 continue
      return
      end
