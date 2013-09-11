      subroutine qsbsj(dk,nk)
      implicit none
      integer nk
      double precision dk
c
      include 'qsglobal.h'
c
      integer i,ir,ik
      double precision k,x
      double precision bessj0,bessj1,bessj
c
      do ik=1,nk
        k=dble(ik)*dk
        do ir=1,nr
          x=k*r(ir)
          bsj(ik,0,ir)=bessj0(x)
          bsj(ik,1,ir)=bessj1(x)
          if(x.gt.2.d0)then
            bsj(ik,2,ir)=bsj(ik,1,ir)*2.d0/x-bsj(ik,0,ir)
          else
            bsj(ik,2,ir)=bessj(2,x)
          endif
          if(x.gt.3.d0)then
            bsj(ik,3,ir)=bsj(ik,2,ir)*4.d0/x-bsj(ik,1,ir)
          else
            bsj(ik,3,ir)=bessj(3,x)
          endif
          bsj(ik,-1,ir)=-bsj(ik,1,ir)
          do i=-1,3
            bsj(ik,i,ir)=bsj(ik,i,ir)*geospr(ir)
          enddo
        enddo
      enddo
c
      return
      end