      subroutine qskern(y,f,k)
      implicit none
c
c     calculation of response function in frequency-wavelength domain
c     y(6,6): solution vector (complex)
c     k: wave number
c     f: frequency. 2*pi*f = cyclar frequency
c     fi: imaginary part of the frequency
c
      double precision f,k
      double complex y(6,6)
c
      include 'qsglobal.h'
c
      integer i,istp,lup,llw
c
      call qswaveno(f,k)
c
      do istp=1,6
        do i=1,6
          y(i,istp)=(0.d0,0.d0)
        enddo
      enddo
c
c     determination of starting upper sublayer
c
      lup=1
      llw=lp
      call qspsv(y,k,lup,llw)
      if(calsh)call qssh(y,k,lup,llw)
c
      return
      end