      subroutine qsam2ve(n,ck,y,c,nj,rsite)
      implicit none
c
c     converting wave amplitudes to displacement-stress vectors
c     by Langer block-diagonal decomposition technique for com-
c     putational efficiency
c
      integer n,nj
      double complex ck
      double complex y(4,nj),c(4,nj)
      logical rsite
c
      include 'qsglobal.h'
c
      integer j
      double complex bp,bs,swap(4)
c
      if(rsite)then
        bp=wbrs(n)*kprs(n)
        bs=wbrs(n)*ksrs(n)
        do j=1,nj
          swap(1)= c(1,j)-c(2,j)
          swap(2)= c(1,j)+c(2,j)
          swap(3)=-c(3,j)+c(4,j)
          swap(4)= c(3,j)+c(4,j)
c
          y(1,j)=kprs(n)*swap(1)+ck*swap(4)
          y(2,j)=wars(n)*swap(2)-bs*swap(3)
          y(3,j)=ck*swap(2)-ksrs(n)*swap(3)
          y(4,j)=bp*swap(1)+wars(n)*swap(4)
        enddo
      else
        bp=wb(n)*kp(n)
        bs=wb(n)*ks(n)
        do j=1,nj
          swap(1)= c(1,j)-c(2,j)
          swap(2)= c(1,j)+c(2,j)
          swap(3)=-c(3,j)+c(4,j)
          swap(4)= c(3,j)+c(4,j)
c
          y(1,j)=kp(n)*swap(1)+ck*swap(4)
          y(2,j)=wa(n)*swap(2)-bs*swap(3)
          y(3,j)=ck*swap(2)-ks(n)*swap(3)
          y(4,j)=bp*swap(1)+wa(n)*swap(4)
        enddo
      endif
c
      return
      end
