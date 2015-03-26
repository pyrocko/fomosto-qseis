      subroutine qsfftinv(icmp,istp)
      implicit none
c
      integer icmp,istp
c
      include 'qsglobal.h'
c
      integer nn,nn2,nch,nch2
      parameter (nn=2*nfmax,nn2=nfmax,nch=2*nrmax,nch2=nrmax)
      integer lf,mf,ir,j,it
      double precision t,pi,pi2,slw0,omi
      double precision f(nn),y0(nch),y(nch)
      double complex s
      double complex cy(nn,nch2)
c
      double complex wvf(nn2)
      save wvf
c
      pi=4.d0*datan(1.d0)
      pi2=2.d0*pi
c
      if(v0.gt.0.d0)then
        slw0=1.d0/v0
      else
        slw0=0.d0
      endif
c
      print *,' Output: '//outfile(icmp,istp)(1:flen(icmp,istp))
      do lf=1,nf
        f(lf)=dble(lf-1)*df
        do ir=1,nr
          cy(lf,ir)=grns(lf,icmp,ir,istp)
        enddo
      enddo
c
      if(v0.gt.0.d0.or.dabs(tstart).gt.0.d0)then
c
c       for time reduction
c
        do lf=1,nf
          do ir=1,nr
            s=dcmplx(-fi,f(lf))
     &         *dcmplx(2.d0*pi*(tstart+r(ir)*slw0),0.d0)
            cy(lf,ir)=cy(lf,ir)*cdexp(s)
          enddo
        enddo
      endif
c
c     seismometer filtering
c
      if(cdabs(asm).ne.1.d0)then
        do lf=1,nf
          do ir=1,nr
            cy(lf,ir)=asm*cy(lf,ir)
          enddo
        enddo
      endif
c
      if(nroot+npole.gt.0)then
        do lf=1,nf
          s=dcmplx(-2.d0*pi*fi,2.d0*pi*f(lf))
          do ir=1,nr
            do j=1,nroot
              cy(lf,ir)=cy(lf,ir)*(s-root(j))
            enddo
            do j=1,npole
              cy(lf,ir)=cy(lf,ir)/(s-pole(j))
            enddo
          enddo
        enddo
      endif
c
c     muliplication with wavelet spectrum
c
      if(iexist.ne.1)then
        call qswavelet(wvf,nf)
        iexist=1
      endif
      do lf=1,nf
        do ir=1,nr
          cy(lf,ir)=cy(lf,ir)*wvf(lf)
        enddo
      enddo
c
      mf=1
      do lf=2*nf,nf+2,-1
        mf=mf+1
        do ir=1,nr
          cy(lf,ir)=dconjg(cy(mf,ir))
        enddo
      enddo
      do ir=1,nr
        cy(nf+1,ir)=(0.d0,0.d0)
      enddo
c
c     convention for Fourier transform:
c     f(t)=\int F(f) exp(i2\pi f t) df
c
      do ir=1,nr
        call four1(cy(1,ir),2*nf,+1)
      enddo
c
      open(20,file=outfile(icmp,istp),status='unknown')
      write(20,'(a,$)')'   T_sec    '
      do ir=1,nr-1
        write(20,'(a4,a1,a1,2a4,$)')'   ',
     &         varbtxt,comptxt(icmp),rcvtxt(ir),'   '
      enddo
      write(20,'(a4,a1,a1,2a4)')'   ',
     &         varbtxt,comptxt(icmp),rcvtxt(nr),'   '
c
      do ir=1,nr
        y0(ir)=0.d0
      enddo
      omi=2.d0*pi*fi
      do it=1,nt
        t=dble(it-1)*dt
        do ir=1,nr
          y(ir)=df*dreal(cy(it,ir))*dexp(-omi*t)
        enddo
        if(wdeg.eq.2)then
          do ir=1,nr
            y(ir)=y0(ir)+y(ir)*dt
            y0(ir)=y(ir)
          enddo
        endif
        write(20,1003)t+tstart
        do ir=1,nr-1
          write(20,1004)y(ir)
        enddo
        write(20,1005)y(nr)
      enddo
      close(20)
1003  format(f12.5,$)
1004  format(E12.4,$)
1005  format(E12.4)
      return
      end
