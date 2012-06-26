cc      program main
      subroutine etasimf(ic,bvalue,tstart,nd,ctmg,rfmg,a,b,c,d,p,
     &                   xm,zz,xx,probx)
c
      include 'sapp_f.h'
c
      implicit real*8(a-h,o-z)
      real*4 r
cc      parameter(ndata=19999)
cc      data bvalue,ctmg   / 1.0, 3.5 /
cc      data a,b,c,d,p / 0.00200,.00400,.00300,2.40000,1.30000 /
c     data a,b,c,d,p / 0.02000,.04000,.00300,0.40000,1.30000 /
cc      data nd / 1000 /
c     data xm / ndata*0.0 /
cc      dimension xx(ndata),xm(ndata),zz(ndata)
      dimension xx(nd),xm(nd),zz(nd)
cc      open(unit=2,file='etasim.open')
cc      read(2,*) ic,bvalue
cc      read(2,*) tstart,nd
cc      read(2,*) ctmg,rfmg
cc      read(2,*) a,b,c,d,p
cc      close(2)
c---
      ix=1992
      iy=1111
      iz=1151
c---
      if(ic.eq.0) then
        do 40 i=1,nd
cc          call pseudo(r)
          call pseudo(r,ix,iy,iz)
          xm(i)=-log10(r)/bvalue +ctmg
          xx(i)=0
   40   continue
cc      else
cc        open(unit=10,file='work.etasim0')
cc        read(10,1009) fmt
cc        write(6,1009) fmt
 1009   format(a)
cc        i=1
cc   50   continue
cc        read(10,*, end=60) idum,e,en,xmi,zzi,depi
cc 1008   format(i5,4f12.5,f5.0)
cc          if(xmi.lt.ctmg) go to 50
cc          write(6,*) idum,e,en,xmi,zzi,depi
cc          xm(i)=xmi
cc          zz(i)=zzi
cc          i=i+1
cc          go to 50
cc   60   continue
cc        close(3)
cc        nd=i-1
      endif
cc      write(6,1002) nd
cc      write(6,1001) (xm(i),i=1,nd)
 1001 format(8f10.5)
 1002 format(i10)
cc      write(6,*) 'b-value =',bvalue, 'cutoff magnitude=', ctmg
cc      write(6,*) 'a,b,c,d,p,bvalue ='
cc      write(6,1001) a,b,c,d,p,bvalue
      if(ic.eq.0) then
        x=0.0
        xity=a
        i=1
      else
c       if(xity.eq.0.0) then
        i=1
   70   continue
        xx(i)=zz(i)
        i=i+1
        if(zz(i).lt.tstart) go to 70
        x=xx(i-1)
cc        if(i.gt.1) call fx(i-1,x,a,b,c,d,p,rfmg,xx,xm,xity)
        if(i.gt.1) call fx1(i-1,x,a,b,c,d,p,rfmg,xx,xm,xity)
cc        write(6,*) i,x,xity
      endif
   10 continue
      uity=xity
   20 continue
cc      call pseudo(r)
      call pseudo(r,ix,iy,iz)
      e=-log(r)/uity
      x=x+e
cc      if(i.gt.1) call fx(i-1,x,a,b,c,d,p,rfmg,xx,xm,xity)
      if(i.gt.1) call fx1(i-1,x,a,b,c,d,p,rfmg,xx,xm,xity)
      probx=xity/uity
cc      if(probx.gt.1.d0) write(6,*) 'prob>1', probx
cc      if(probx.gt.1.d0) stop 119
      if(probx.gt.1.d0) return
cc      call pseudo(r)
      call pseudo(r,ix,iy,iz)
      if(r.gt.probx) go to 10
      xx(i)=x
cc      if(i.eq.1) write(6,*) 'event#  occurrence time'
cc      if(mod(i,10).eq.0) write(6,*) i,xx(i)
      uity=xity+b/c**p*exp(d*(xm(i)-rfmg))
      if(i.ge.nd) go to 30
      i=i+1
      go to 20
c
   30 continue
      t=xx(i)+e
      nn=i
cc      open(unit=1,file='work.etas')
c     write(1,1004) nn,t,ctmg,bvalue
c     write(1,1006) a,b,c,d,p
c     write(1,*) '(5x,4f12.5,f5.0)'
*     write(1,1003) (i,xm(i),xx(i),i=1,nn)
cc      if(ic.eq.0) 
cc     & write(1,*) 'ETAS simlated data: MAGs generated bt G-R'
cc      if(ic.eq.1) 
cc     & write(1,*) 'ETAS simlated data: MAGS from [work.etas]'
      zr=0.0
      izr=0
cc      write(1,1007) (i,zr,zr,xm(i),xx(i),zr,izr,izr,izr, i=1,nn)
cc      close(unit=1)
 1003 format(i5,24x,2f12.5)
 1007 format(i5,5f12.5,i4,2i3)
 1004 format('simulated data',' N,T=',i5,f9.3,';   M>=',f4.2,
     &       ';   b-value=',f4.2)
 1006 format(' ETAS parameter; mu,K_0,c,alpha,p =',5f8.5)
 1005 format(1h ,6f13.4)
cc      stop
      return
      end
c
cc      subroutine fx(i,x,a,b,c,d,p,rfmg,xx,xm,xity)
      subroutine fx1(i,x,a,b,c,d,p,rfmg,xx,xm,xity)
      implicit real*8(a-h,o-z)
      dimension xx(1),xm(1)
      xity=a
ctren xity=a-x*0.2821d-6
      do 10 j=1,i
   10 xity=xity+b/(x-xx(j)+c)**p*exp(d*(xm(j)-rfmg))
      return
      end
      subroutine pseud0(r)
c     generation of pseudo-random numbers
c     data ir/584287/
      data ir/574289/
      ir=ir*48828125
      if(ir) 10,20,20
   10 ir=(ir+2147483647)+1
   20 r=float(ir)*0.4656613e-9
      return
      end
cc      subroutine pseudo(random)
      subroutine pseudo(random,ix,iy,iz)
c     wichmann+hill (1982) Appl. Statist 31
cc      data ix,iy,iz /1992,1111,1151/
      ix=171*mod(ix,177)-2*(ix/177)
      iy=172*mod(iy,176)-35*(iy/176)
      iz=170*mod(iz,178)-63*(iz/178)
      if (ix.lt.0) ix=ix+30269
      if (iy.lt.0) iy=iy+30307
      if (iz.lt.0) iz=iz+30323
      random=mod(float(ix)/30269.0+float(iy)/30307.0+
     &            float(iz)/30323.0,1.0)
      return
      end
