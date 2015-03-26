      subroutine qsmatrix0(a,ck,z,n)
      implicit none
c
      integer n
      double precision z
      double complex ck
      double complex a(4,4)
c
      include 'qsglobal.h'
c
      double complex c0,c1,c2
      data c0,c1,c2/(0.d0,0.d0),(1.d0,0.d0),(2.d0,0.d0)/
c
      integer i
      double complex ck2,cz,cxi,cfac,b(4,4)
      double complex cxps,cdkp,cdks,cps,csp,cdeps,cdesp
c
      ck2=ck*ck
c
      a(1,1)=kp(n)
      a(2,1)=c2*cmu(n)*ck2-acc(n)
      a(3,1)=ck
      a(4,1)=c2*ck*cmu(n)*kp(n)
      a(1,2)=-kp(n)
      a(2,2)=a(2,1)
      a(3,2)=ck
      a(4,2)=-a(4,1)
c
      b(1,3)=ck
      b(2,3)=c2*ck*cmu(n)*ks(n)
      b(3,3)=ks(n)
      b(4,3)=c2*cmu(n)*ck2-acc(n)
      b(1,4)=ck
      b(2,4)=-b(2,3)
      b(3,4)=-ks(n)
      b(4,4)=b(4,3)
c
      cz=dcmplx(z,0.d0)
      cxi=cla(n)+c2*cmu(n)
c
c     cxps=(kp-ks)*z
c     cdkp=(k-kp)*k^2*vp^2/omiga^2
c     cdks=(k-ks)*k^2*vp^2/omiga^2
c     cdeps=(1-exp(cxps))*k^2*vp^2/omiga^2
c     cdesp=(1-exp(-cxps))*k^2*vp^2/omiga^2
c
      cxps=acc(n)*(cla(n)+cmu(n))*cz/(cmu(n)*cxi*(kp(n)+ks(n)))
      cdkp=ck2/(ck+kp(n))
      cdks=ck2*cxi/((ck+ks(n))*cmu(n))
      if(dabs(z).gt.0.d0)then
        i=1
        cps=(1.d0,0.d0)
        csp=(1.d0,0.d0)
        cdeps=(1.d0,0.d0)
        cdesp=(1.d0,0.d0)
10      i=i+1
        cfac=cxps/dcmplx(dble(i),0.d0)
        cps=cps*cfac
        csp=-csp*cfac
        cdeps=cdeps+cps
        cdesp=cdesp+csp
        if(cdabs(cps).gt.1.0d-12)goto 10
        cfac=(cla(n)+cmu(n))*ck2*cz/(cmu(n)*(kp(n)+ks(n)))
        cdeps=-cdeps*cfac
        cdesp=cdesp*cfac
      else
        cdeps=0.d0
        cdesp=0.d0
      endif
c
c     define stable propagator for omiga -> 0:
c
c     a(i,3)=[a(i,3)-a(i,1)*exp(cxps)]*k^2*vp^2/omiga^2
c     a(i,4)=[a(i,4)+a(i,2)*exp(-cxps)]*k^2*vp^2/omiga^2
c
      a(1,3)=cdkp+kp(n)*cdeps
      a(2,3)=c2*cmu(n)*ck*(ck*cdeps-cdks)+cxi*ck2*cdexp(cxps)
      a(3,3)=-cdks+ck*cdeps
      a(4,3)=c2*cmu(n)*ck*(cdkp+kp(n)*cdeps)-cxi*ck2
      a(1,4)=cdkp+kp(n)*cdesp
      a(2,4)=c2*cmu(n)*ck*(cdks-ck*cdesp)-cxi*ck2*cdexp(-cxps)
      a(3,4)=cdks-ck*cdesp
      a(4,4)=c2*cmu(n)*ck*(cdkp+kp(n)*cdesp)-cxi*ck2
c
      return
      end