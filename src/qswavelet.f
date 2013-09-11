      subroutine qswavelet(wvf,mm2)
      implicit none
      integer mm2
      double complex wvf(mm2)
c
      include 'qsglobal.h'
c
      integer l,n
      double precision f,omi,x,dt0
      double complex alfa,beta,gamma,eta
c
      double precision pi,pi2,eps
      data pi,pi2,eps/3.14159265358979d0,6.28318530717959d0,1.0d-04/
c
      if(wdeg.ne.0)then
c
c       for wavelet: normalized square half-sinus
c
        do l=1,mm2
          f=df*dble(l-1)
          x=f*tau
          if(x.eq.0.d0)then
            wvf(l)=(1.d0,0.d0)
          else if(x.ge.1.d0-eps.and.x.le.1.d0+eps)then
            wvf(l)=dcmplx(-1.d0/x/(1+x),0.d0)
          else
            wvf(l)=dcmplx(0.d0,1.d0/(pi2*x*(1.d0+x)*(1.d0-x)))
     &               *(cdexp(dcmplx(0.d0,-pi2*x))-(1.d0,0.d0))
          endif
        enddo
      else
c
c       user's own wavelet function
c
        dt0=tau/dble(nn0-1)
c
        wvf(1)=dcmplx(0.5d0*(wv0(1)+wv0(nn0)),0.d0)
        do n=2,nn0-1
          wvf(1)=wvf(1)+dcmplx(wv0(n),0.d0)
        enddo
        wvf(1)=wvf(1)*dcmplx(dt0,0.d0)
c
        do l=2,mm2
          wvf(l)=(0.d0,0.d0)
          omi=2.d0*pi*df*dble(l-1)
          alfa=cdexp(dcmplx(0.d0,-omi*dt0))
          beta=(alfa-(1.d0,0.d0))*dcmplx(0.d0,1.d0/omi)
          gamma=alfa*dcmplx(0.d0,1.d0/omi)
     &           -beta*dcmplx(0.d0,1.d0/omi/dt0)
          eta=(1.d0,0.d0)
          do n=1,nn0-1
            wvf(l)=wvf(l)+eta*(dcmplx(wv0(n),0.d0)*(beta-gamma)
     &            +dcmplx(wv0(n+1),0.d0)*gamma)
            eta=eta*alfa
          enddo
        enddo
      endif
      if(iexist.ne.1)iexist=1
      return
      end
