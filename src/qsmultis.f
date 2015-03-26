      subroutine qsmultis(grnexist)
      implicit none
      logical grnexist
      include 'qsglobal.h'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     this program calculates convolution integral (summation)         c
c     of discrete sources to model the seismic ground motion of        c
c     an earthquake with arbitrary moment tensor distribution          c
c                                                                      c
c     The inputs are:                                                  c
c                                                                      c
c     1) Green's functions for surface movement (displacement          c
c        or velocity);                                                 c
c                                                                      c
c     2) source parameters including the source coordinates, the       c
c        strike, dip and rake angles as well as the released           c
c        energy (moment);                                              c
c                                                                      c
c     3) receiver parameters including the receiver coordinates.       c
c                                                                      c
c     The outputs are the surface motion at receiver sites.            c
c                                                                      c
c                                                                      c
c     ~~~~~~~~~~~~~~~~~~                                               c
c     Coordinate system:                                               c
c     ~~~~~~~~~~~~~~~~~~                                               c
c     cylindrical (z,r,t) with z = downward,                           c
c                t = azmuth angle from north;                          c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer ir,lf,icmp,istp,unit
      double precision az1,az2
      double precision weight(0:5)
      double precision t(2*nfmax),y1(nrmax),y2(nrmax)
      character*80 textline
c
      double precision deg2rad
      data deg2rad/1.745329251994328d-02/
c
      if(grnexist)then
        do istp=1,4
          do icmp=1,4
            if(ms(istp).eq.0.and.icmp.eq.3)then
              do lf=1,nf
                do ir=1,nr
                  grns(lf,icmp,ir,istp)=(0.d0,0.d0)
                enddo
              enddo
            else
              unit=istp*10+icmp
              open(unit,file=outfile(icmp,istp),status='old')
              read(unit,'(a1)')textline
              do lf=1,nf
                read(unit,*,end=100)t(2*lf-1),(y1(ir),ir=1,nr)
                read(unit,*,end=100)t(2*lf),(y2(ir),ir=1,nr)
                do ir=1,nr
                  grns(lf,icmp,ir,istp)=dcmplx(y1(ir),y2(ir))
                enddo
              enddo
100           close(unit)
            endif
          enddo
        enddo
      endif
c
      weight(0)=(mtensor(1)+mtensor(2)+mtensor(3))/3.d0
      weight(1)=mtensor(4)
      weight(2)=mtensor(6)
      weight(3)=mtensor(3)-weight(0)
      weight(4)=0.5d0*(mtensor(1)-mtensor(2))
      weight(5)=mtensor(5)
c
      do ir=1,nr
        az1=azimuth(ir)*deg2rad
        az2=2.d0*az1
        do lf=1,nf
         grns(lf,1,ir,7)=grns(lf,1,ir,1)*dcmplx(weight(0),0.d0)
     &                  +grns(lf,1,ir,2)
     &     *dcmplx(weight(1)*dsin(az2)+weight(4)*dcos(az2),0.d0)
     &                  +grns(lf,1,ir,3)
     &     *dcmplx(weight(2)*dcos(az1)+weight(5)*dsin(az1),0.d0)
     &                  +grns(lf,1,ir,4)*dcmplx(weight(3),0.d0)
c
         grns(lf,2,ir,7)=grns(lf,2,ir,1)*dcmplx(weight(0),0.d0)
     &                  +grns(lf,2,ir,2)
     &     *dcmplx(weight(1)*dsin(az2)+weight(4)*dcos(az2),0.d0)
     &                  +grns(lf,2,ir,3)
     &     *dcmplx(weight(2)*dcos(az1)+weight(5)*dsin(az1),0.d0)
     &                  +grns(lf,2,ir,4)*dcmplx(weight(3),0.d0)
c
         grns(lf,3,ir,7)=grns(lf,3,ir,2)
     &     *dcmplx(weight(1)*dcos(az2)-weight(4)*dsin(az2),0.d0)
     &                  +grns(lf,3,ir,3)
     &     *dcmplx(weight(2)*dsin(az1)-weight(5)*dcos(az1),0.d0)
c
         grns(lf,4,ir,7)=grns(lf,4,ir,1)*dcmplx(weight(0),0.d0)
     &                  +grns(lf,4,ir,2)
     &     *dcmplx(weight(1)*dsin(az2)+weight(4)*dcos(az2),0.d0)
     &                  +grns(lf,4,ir,3)
     &     *dcmplx(weight(2)*dcos(az1)+weight(5)*dsin(az1),0.d0)
     &                  +grns(lf,4,ir,4)*dcmplx(weight(3),0.d0)
        enddo
      enddo
      if(grnexist)then
        do icmp=1,4
          if(fsel(icmp,7).eq.1)then
            unit=40+icmp
            open(unit,file=outfile(icmp,7),status='unknown')
            write(unit,'(a,$)')'   T_sec    '
            do ir=1,nr-1
              write(unit,'(a4,a1,a1,2a4,$)')'   ',
     &             varbtxt,comptxt(icmp),rcvtxt(ir),'   '
            enddo
            write(unit,'(a4,a1,a1,2a4)')'   ',
     &             varbtxt,comptxt(icmp),rcvtxt(nr),'   '
            do lf=1,nf
              write(unit,'(f12.5,$)')t(2*lf-1)
              do ir=1,nr-1
                write(unit,'(E12.4,$)')dreal(grns(lf,icmp,ir,7))
              enddo
              write(unit,'(E12.4)')dreal(grns(lf,icmp,nr,7))
c
              write(unit,'(f12.5,$)')t(2*lf)
              do ir=1,nr-1
                write(unit,'(E12.4,$)')dimag(grns(lf,icmp,ir,7))
              enddo
              write(unit,'(E12.4)')dimag(grns(lf,icmp,nr,7))
            enddo
            close(unit)
          endif
        enddo
      endif
c
      return
      end
