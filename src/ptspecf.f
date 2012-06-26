cc      program ptspec
      subroutine ptspecf(t,n,t0,tmpr,tmp,prd,nh,nt,is,prb,r1,rwx,rwy,
     & phs,wt,ht,w,h,g)
c
      include 'sapp_f.h'
c
c     this program provides perodgrams of point process data with
c     the significant bands ( 0.90, 0.95 and 0.99 ) assuming the
c     stationary poisson process.   powers of interested frequencies
c     or periods are plotted by a particular sign.
c
c     structure of the program
c
c          ptspec
c             i---input
c             i---cycle
c             i---period
c             i---printr
c             +---smooth
c
c     this is designed by y. ogata and t. ozaki, and programmed by
c     k. katsura and t. ozaki, inst. statist. math., tokyo. (31/01/85)
c
c     reference
c     d. vere-jones and t. ozaki (1982). "some examples of statistical
c        estimation applied to earthquake data, 1. cyclic poisson and
c        self-exciting models."  ann. inst. statist. math., vol. 34,
c        no. 1, pp. 189-207.
c
      implicit real*8 (a-h,o-z)
cc      dimension t(5000)
cc      dimension h(5000),g(5000),s(5000)
cc      dimension w(5000)
cc      dimension tmpr(100),wt(200),ht(200),gt(200)
      dimension t(n)
      dimension h(nh+1),g(nh+1),s(nh+1)
      dimension w(nh+1)
      dimension tmpr(nt),wt(nt),ht(nt),gt(nt)
cc      real*4 widthx,widthy
cc      call input(t,n,t0,pi2,rpt,tmpr,tmp,prd,nh1,nt,is,ipl)
      nh1=nh+1
      pi=3.1415926536d 00
      pi2=pi*2.d 00
      f=1.d0/tmp
      om=pi2*f
      rpt=om/nh
cc      call cycle(t,n,prd)
      call cycle(t,n,prd,prb,r1,rwx,rwy,phs)
      call period(h,g,w,n,t,nh1,ht,gt,wt,nt,rpt,t0,pi2,tmpr)
      if(is.gt.1) go to 10
c
cc      isw1=1
cc      call printr(nh1,w,h,g,nt,wt,ht,t0,tmp,isw1)
c
      go to 20
   10 call smooth(s,h,nh1,is)
c-----
      do 11 i=1,nh1
      h(i)=s(i)
   11 continue
c-----
cc      isw1=2
cc      call printr(nh1,w,s,g,nt,wt,ht,t0,tmp,isw1)
   20 continue
      return
      end
      subroutine period(h,g,w,n,t,nh1,ht,gt,wt,nt,rpt,t0,pi2,tmpr)
      implicit real*8 (a-h,o-z)
cc      dimension t(5000)
cc      dimension h(5000),g(5000)
cc      dimension w(5000)
cc      dimension tmpr(100),wt(200),ht(200),gt(200)
      dimension t(n)
      dimension h(nh1),g(nh1)
      dimension w(nh1)
      dimension tmpr(nt),wt(nt),ht(nt),gt(nt)
      do 120 i=1,nh1
      om=(i-1)*rpt
      w(i)=(i-1)*rpt
      a=0.d 00
      b=0.d 00
      do 110 j=1,n
      a=a+dcos(t(j)*om)
      b=b+dsin(t(j)*om)
110   continue
      h(i)=(a*a+b*b)/t0
      g(i)=h(i)/pi2
      ram=n/t0
      g(i)=g(i)/(ram/pi2)
      h(i)=dlog10(h(i))*10.d 00
120   continue
      do 1200 i=1,nt
      om=pi2/tmpr(i)
      wt(i)=om
      a=0.d 00
      b=0.d 00
      do 1100 j=1,n
      a=a+dcos(t(j)*om)
      b=b+dsin(t(j)*om)
1100  continue
      ht(i)=(a*a+b*b)/t0
      gt(i)=ht(i)/pi2
      ram=n/t0
      gt(i)=gt(i)/(ram/pi2)
      ht(i)=dlog10(ht(i))*10.d 00
1200  continue
      return
      end
      subroutine smooth(s,h,nh1,is)
      implicit real*8 (a-h,o-z)
      dimension h(nh1),s(nh1)
      is1=2*is-1
      do 140 i=1,nh1
      ij=0
      s(i)=0.0
      do 130 j=1,is1
      ij1=i-is+j-1
      if(ij1.lt.1.or.ij1.gt.nh1) go to 130
      ij=ij+1
      s(i)=s(i)+h(ij1)
130   continue
      s(i)=s(i)/ij
140   continue
      return
      end
