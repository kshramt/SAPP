cc      program etarpp
      subroutine etarppf(xx,xmg,xmag0,nn,xini,n,zts,zte,tstart0,
     &                   x,ntstar)
c-----------------------------------------------------------------------
c Subroutine FUNC4 corresponds to the exact version and FUNC9 to the
c approximate version: see the references below.
c 
c     References:
c     Ogata, Y. (1988). J. Amer. Statist. Soc., 83, pp. 9-27.
c       ---     (1989). Tectonophysics, 169, pp. 159-174.
c       ---     (1992). J. Geophys. Res. 97, pp. 19845-19871.
c     Ogata, Y., Matsu'ura S.R., Katsura, K. (1993). submitted to 
c                Geophys. Res. Letters.
c-----------------------------------------------------------------------
c
      include 'sapp_f.h'
c
      implicit real * 8 (a-h,o-z)
cc      parameter(ldata=17777, npara=5)
cc      common/xyod/xx(ldata),xmg(ldata),xmag0
cc      common/param/xini(npara),n
cc      common /range/tstart,ntstar
cc      common t,nn,mm,iappr,nfunct
      dimension xx(nn),xmg(nn),xini(n)
      dimension x(nn)
c
cc      call input
      t=zte-zts
      tstart=tstart0-zts
      ntstar=0
      do 10 i=1,nn
      if(xx(i).lt.tstart) ntstar=i
      xx(i)=xx(i)-zts
   10 continue
c
cc      do 10 i=2,nn
cc      if(xx(i).ge.xx(i-1)) go to 10
cc      write(6,*) 'reverse occurrence time'
cc      write(6,*) i,xx(i),xx(i-1),xmg(i),xmg(i-1)
cc   10 continue
cc      write(6,6) nfunct
cc      write(6,4) t,nn,mm
cc      write(6,5) xmag0
cc      write(6,*)
    6 format(1h ,' funct = ',i5)
    5 format(1h ,'reference magnitudes; xmag0',5x,f10.4)
    4 format(1h ,'(T,nn,mm) =',5x,f10.4,2i6)
    3 format(1h ,10f12.4/(1h ,10f12.4))
    2 format(f10.2,2i10)
    1 format(8f10.2)
c
cc      call residual
      call eresidual(xx,xmg,xmag0,nn,xini,n,t,tstart,ntstar,x)
c
   20 continue
      return
      end
c***********************************************************************
cc      subroutine residual
      subroutine eresidual(xx,xmg,xmag0,nn,a,n,t,tstart,ntstar,x)
      implicit real * 8 (a-h,o-z)
cc      parameter(ldata=17777, npara=5)
cc      common/param/a(npara),n
cc      common /range/tstart,ntstar
cc      common t,nn,mm,iappr,nfunct
cc      common/xyod/xx(ldata),xmg(ldata),xmag0
c     common/epi/ xp(ldata),yp(ldata)
cc      common/hyp/ xp(ldata),yp(ldata),dep(ldata)
cc      dimension x(ldata),xmg0(ldata)
      dimension xx(nn),xmg(nn),a(n)
      dimension x(nn),xmg0(nn)
      func41(t,tx,xm,a3,a4)=(log(t-tx+a3)-log(a3))*exp(a4*xm)
      func4p(t,tx,xm,a3,a4,a5)=1.d0/(1.d0-a5)*((t-tx+a3)**(1.d0-a5)
     &      -a3**(1.d0-a5))*exp(a4*xm)
c
      do 40 i=1,nn
   40 xmg0(i)=xmg(i)-xmag0
c
      chtsta=a(1)*tstart
      ft=0.0
      do 30 j=1,ntstar
      if(a(5).eq.1.d0) ft=ft+func41(tstart,xx(j),xmg0(j),a(3),a(4))
      if(a(5).ne.1.d0) ft=ft+func4p(tstart,xx(j),xmg0(j),a(3),a(4),a(5))
   30 continue
      chtsta=chtsta+a(2)*ft
c
      x(1)=a(1)*xx(1)-chtsta
      do 10 i=2,nn
      ft=0.0
      do 20 j=1,i-1
      if(a(5).eq.1.d0) ft=ft+func41(xx(i),xx(j),xmg0(j),a(3),a(4))
      if(a(5).ne.1.d0) ft=ft+func4p(xx(i),xx(j),xmg0(j),a(3),a(4),a(5))
   20 continue
      x(i)=a(1)*xx(i)+a(2)*ft-chtsta
   10 continue
c
c     write(6,1002) (i-ntstar,xmg(i),x(i),x(i)-x(i-1),i=1,nn)
cc      write(6,1002) (i-ntstar,xmg(i),x(i),x(i)-x(i-1),xx(i),i=1,nn)
cc      open(unit=1,file='work.res')
*     write(1,1001) (i-ntstar,xmg(i),x(i),i=1,nn)
c     write(1,1003) (i-ntstar,xp(i),yp(i),xmg(i),x(i),i=1,nn)
*     write(1,1003) (i-ntstar,xp(i),yp(i),xmg(i),x(i),xx(i),i=1,nn)
c     write(1,1004) (i-ntstar,xp(i),yp(i),xmg(i),x(i),xx(i),i=1,nn)
cc      write(1,1005) (i-ntstar,xp(i),yp(i),xmg(i),xx(i),
cc     &                                         dep(i),x(i),i=1,nn)
c1003 format(i5,4f12.5,5x)
 1003 format(i5,5f12.5)
 1004 format(i6,2f12.5,f6.2,2f15.5)
 1005 format(i5,2f12.5,f12.1,f12.5,f8.2,2x,f12.5)
cc      close(unit=1)
 1001 format(i5,24x,2f12.5,5x)
c1002 format(i5,24x,3f12.5)
 1002 format(i5,4f12.5)
      return
      end
