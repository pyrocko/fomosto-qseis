      subroutine qssource(ros,vps,vss)
      implicit none
c
      double precision ros,vps,vss
c
      include 'qsglobal.h'
c
      integer i,istp
      double precision pi,pi2
c
      do istp=1,6
        do i=1,6
          sfct0(i,istp)=0.d0
          sfct1(i,istp)=0.d0
        enddo
      enddo
c
      pi=4.d0*datan(1.d0)
      pi2=2.d0*pi
c
c     istp = 1
c     explosion source (m11=m22=m33=1)
c
      ms(1)=0
      ics(1)=1
      sfct0(1,1)=-1.d0/(pi2*ros*vps*vps)
      sfct1(4,1)=-(vss/vps)**2/pi
c
c     istype = 2
c     strike-slip (m12=m21=1)
c
      ms(2)=2
      ics(2)=-1
      sfct1(4,2)=1.d0/pi2
      sfct1(6,2)=-sfct1(4,2)
c
c     istype = 3
c     dip-slip (m13=m31=1)
c
      ms(3)=1
      ics(3)=1
      sfct0(3,3)=-1.d0/(pi2*ros*vss*vss)
      sfct0(5,3)=sfct0(3,3)
c
c     istp = 4
c     compensated linear vector dipole (CLVD) (m11=m22=-1/2, M33=1)
c
      ms(4)=0
      ics(4)=1
      sfct0(1,4)=-1.d0/(pi2*ros*vps*vps)
      sfct1(4,4)=(3.d0-4.d0*(vss/vps)**2)/(2.d0*pi2)
c
c     istp = 5
c     vertical-single-force (fz=1)
c
      ms(5)=0
      ics(5)=1
      sfct0(2,5)=1.d0/pi2
c
c     istp = 6
c     horizontal-single-force (fx=1)
c
      ms(6)=1
      ics(6)=1
      sfct0(4,6)=1.d0/pi2
      sfct0(6,6)=sfct0(4,6)
c
      return
      end
