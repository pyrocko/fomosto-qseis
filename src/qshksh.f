      subroutine qshksh(hk,z,n,rsite)
      implicit none
c
      integer n
      double precision z
      double complex hk(2,2)
      logical rsite
c
      include 'qsglobal.h'
c
      double complex cx,cem,cch,csh
c
      if(rsite)then
        cx=ksrs(n)*dcmplx(2.d0*z,0.d0)
        if(z.gt.0.d0)then
          cem=cdexp(-cx)
          cch=(0.5d0,0.d0)*((1.d0,0.d0)+cem)
          csh=(0.5d0,0.d0)*((1.d0,0.d0)-cem)
        else
          cem=cdexp(cx)
          cch=(0.5d0,0.d0)*((1.d0,0.d0)+cem)
          csh=-(0.5d0,0.d0)*((1.d0,0.d0)-cem)
        endif
c
c       propagator matrix for SH waves
c
        hk(1,1)=cch
        hk(1,2)=csh/(cmurs(n)*ksrs(n))
        hk(2,1)=csh*cmurs(n)*ksrs(n)
        hk(2,2)=cch
      else
        cx=ks(n)*dcmplx(2.d0*z,0.d0)
        if(z.gt.0.d0)then
          cem=cdexp(-cx)
          cch=(0.5d0,0.d0)*((1.d0,0.d0)+cem)
          csh=(0.5d0,0.d0)*((1.d0,0.d0)-cem)
        else
          cem=cdexp(cx)
          cch=(0.5d0,0.d0)*((1.d0,0.d0)+cem)
          csh=-(0.5d0,0.d0)*((1.d0,0.d0)-cem)
        endif
c
c       propagator matrix for SH waves
c
        hk(1,1)=cch
        hk(1,2)=csh/(cmu(n)*ks(n))
        hk(2,1)=csh*cmu(n)*ks(n)
        hk(2,2)=cch
      endif
      return
      end
