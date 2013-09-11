      subroutine qspsv(y,k,lup,llw)
      implicit none
c
c     calculation of response to p-sv source
c     y(6,6): solution vector (complex)
c     k: wave number
c
      integer lup,llw
      double precision k
      double complex y(6,6)
c
      include 'qsglobal.h'
c
c     work space
c
      integer i,istp,j,l,n,key
      double complex cfac,ck,ch0,pwave,swave
      double complex y0(4,2),c0(4,2),c1(4,2),b(4,6),b0(4,6)
      double complex cinc(4,6)
      double complex y1(4,2),yup(4,2),ylw(4,2),orth(2,2)
      double complex coef(4,4),cnorm(2),coefrs(2,2),brs(2,6)
      logical rsite
c
      double complex c2
	data c2/(2.d0,0.d0)/
c
      ck=dcmplx(k,0.d0)
c
c===============================================================================
c
c     matrix propagation from surface to source
c
c     determination of starting upper sublayer
c
      rsite=.false.
c
      if(lup.eq.1.and.isurf.eq.0)then
        do j=1,2
          do i=1,4
            yup(i,j)=(0.d0,0.d0)
          enddo
        enddo
        yup(1,1)=(1.d0,0.d0)
        if(vpair.gt.0.d0)yup(2,1)=-accair/kpair
        yup(3,2)=(1.d0,0.d0)
      else
        n=nno(lup)
c
        yup(1,1)=kp(n)
        yup(2,1)=wa(n)
        yup(3,1)=ck
        yup(4,1)=wb(n)*kp(n)
c
        yup(1,2)=ck
        yup(2,2)=wb(n)*ks(n)
        yup(3,2)=ks(n)
        yup(4,2)=wa(n)
      endif
      if(lup.eq.lzr)call cmemcpy(yup,y0,8)
c
      do l=lup+1,ls
        ch0=dcmplx(hp(l-1),0.d0)
        n=nno(l-1)
c
c       determination of propagation matrix
c
        call qsve2am(n,ck,yup,c0,2,rsite)
        pwave=cdexp(-kp(n)*ch0)
        swave=cdexp(-ks(n)*ch0)
c
c       orthonormalization of the p-sv modes
c
        cfac=(1.d0,0.d0)/(c0(3,2)*c0(1,1)-c0(1,2)*c0(3,1))
        orth(1,1)=c0(3,2)*cfac
        orth(1,2)=-c0(1,2)*cfac
        orth(2,1)=-c0(3,1)*cfac
        orth(2,2)=c0(1,1)*cfac
        call caxcb(c0,orth,4,2,2,c1)
c
        if(l.gt.lzr)then
c
c         additional normalization to avoid overflow
c
          do i=1,2
            orth(i,1)=orth(i,1)*pwave
            orth(i,2)=orth(i,2)*swave
          enddo
          call caxcb(y0,orth,4,2,2,y1)
          call cmemcpy(y1,y0,8)
        endif
c
c       for partial solution only!
c
        if(.not.pup(l-1))then
          c1(2,1)=(0.d0,0.d0)
          c1(4,1)=(0.d0,0.d0)
          if(l.gt.lzr)then
            do j=1,4
              y0(j,1)=(0.d0,0.d0)
            enddo
          endif
        endif
        if(.not.svup(l-1))then
          c1(2,2)=(0.d0,0.d0)
          c1(4,2)=(0.d0,0.d0)
          if(l.gt.lzr)then
            do j=1,4
              y0(j,2)=(0.d0,0.d0)
            enddo
          endif
        endif
        if(.not.pdw(l-1))then
          c1(2,1)=(0.d0,0.d0)
          c1(2,2)=(0.d0,0.d0)
          if(l-1.eq.lzr)then
            call qsve2am(n,ck,y0,c0,2,rsite)
            c0(2,1)=(0.d0,0.d0)
            c0(2,2)=(0.d0,0.d0)
            call qsam2ve(n,ck,y0,c0,2,rsite)
          endif
        endif
        if(.not.svdw(l-1))then
          c1(4,1)=(0.d0,0.d0)
          c1(4,2)=(0.d0,0.d0)
          if(l-1.eq.lzr)then
            call qsve2am(n,ck,y0,c0,2,rsite)
            c0(4,1)=(0.d0,0.d0)
            c0(4,2)=(0.d0,0.d0)
            call qsam2ve(n,ck,y0,c0,2,rsite)
          endif
        endif
c
c       end partial solution procedure!
c
c        c1(1,1)=c1(1,1)
        c1(2,1)=c1(2,1)*pwave*pwave
        c1(3,1)=(0.d0,0.d0)
        c1(4,1)=c1(4,1)*pwave*swave
c
        c1(1,2)=(0.d0,0.d0)
        c1(2,2)=c1(2,2)*swave*pwave
c        c1(3,2)=c1(3,2)
        c1(4,2)=c1(4,2)*swave*swave
c
        call qsam2ve(n,ck,yup,c1,2,rsite)
c
        if(l.eq.lzr)call cmemcpy(yup,y0,8)
      enddo
c
c===============================================================================
c
c     matrix propagation from half-space to source
c
      n=nno(llw)
      ylw(1,1)=-kp(n)
      ylw(2,1)=wa(n)
      ylw(3,1)=ck
      ylw(4,1)=-wb(n)*kp(n)
c
      ylw(1,2)=ck
      ylw(2,2)=-wb(n)*ks(n)
      ylw(3,2)=-ks(n)
      ylw(4,2)=wa(n)
c
      if(llw.gt.ls.and.llw.eq.lzr)call cmemcpy(ylw,y0,8)
c
      do l=llw-1,ls,-1
        ch0=dcmplx(hp(l),0.d0)
        n=nno(l)
c
c       determination of propagation matrix
c
        call qsve2am(n,ck,ylw,c0,2,rsite)
        pwave=cdexp(-kp(n)*ch0)
        swave=cdexp(-ks(n)*ch0)
c
c       orthonormalization of the p-sv modes
c
        cfac=(1.d0,0.d0)/(c0(4,2)*c0(2,1)-c0(2,2)*c0(4,1))
        orth(1,1)=c0(4,2)*cfac
        orth(1,2)=-c0(2,2)*cfac
        orth(2,1)=-c0(4,1)*cfac
        orth(2,2)=c0(2,1)*cfac
        call caxcb(c0,orth,4,2,2,c1)
c
        if(l.lt.lzr)then
c
c         additional normalization to avoid overflow
c
          do i=1,2
            orth(i,1)=orth(i,1)*pwave
            orth(i,2)=orth(i,2)*swave
          enddo
          call caxcb(y0,orth,4,2,2,y1)
          call cmemcpy(y1,y0,8)
        endif
c
c       for partial solution only!
c
        if(.not.pup(l))then
          c1(1,1)=(0.d0,0.d0)
          c1(1,2)=(0.d0,0.d0)
          if(l+1.eq.lzr)then
            call qsve2am(n,ck,y0,c0,2,rsite)
            c0(1,1)=(0.d0,0.d0)
            c0(1,2)=(0.d0,0.d0)
            call qsam2ve(n,ck,y0,c0,2,rsite)
          endif
        endif
        if(.not.svup(l))then
          c1(3,1)=(0.d0,0.d0)
          c1(3,2)=(0.d0,0.d0)
          if(l+1.eq.lzr)then
            call qsve2am(n,ck,y0,c0,2,rsite)
            c0(3,1)=(0.d0,0.d0)
            c0(3,2)=(0.d0,0.d0)
            call qsam2ve(n,ck,y0,c0,2,rsite)
          endif
        endif
        if(.not.pdw(l))then
          c1(1,1)=(0.d0,0.d0)
          c1(3,1)=(0.d0,0.d0)
          if(l.lt.lzr)then
            do j=1,4
              y0(j,1)=(0.d0,0.d0)
            enddo
          endif
        endif
        if(.not.svdw(l))then
          c1(1,2)=(0.d0,0.d0)
          c1(3,2)=(0.d0,0.d0)
          if(l.lt.lzr)then
            do j=1,4
              y0(j,2)=(0.d0,0.d0)
            enddo
          endif
        endif
c
c       end partial solution procedure!
c
        c1(1,1)=c1(1,1)*pwave*pwave
c        c1(2,1)=c1(2,1)
        c1(3,1)=c1(3,1)*pwave*swave
        c1(4,1)=(0.d0,0.d0)
c
        c1(1,2)=c1(1,2)*swave*pwave
        c1(2,2)=(0.d0,0.d0)
        c1(3,2)=c1(3,2)*swave*swave
c        c1(4,2)=c1(4,2)
c
        call qsam2ve(n,ck,ylw,c1,2,rsite)
        if(l.gt.ls.and.l.eq.lzr)call cmemcpy(ylw,y0,8)
      enddo
c
c===============================================================================
c     source function
c===============================================================================
c
      do istp=1,6
        do i=1,4
          b(i,istp)=dcmplx(sfct0(i,istp)+k*sfct1(i,istp),0.d0)
        enddo
      enddo
      do i=1,4
        do j=1,2
          coef(i,j)=yup(i,j)
          coef(i,j+2)=-ylw(i,j)
        enddo
      enddo
      if(ipath.eq.1)call cmemcpy(b,b0,24)
      key=0
      call cdgemp(coef,b,4,6,0.d0,key)
      if(key.eq.0)then
        print *,'warning in qspsv: anormal exit from cdgemp!'
        return
      endif
      if(lzr.le.ls)then
        do istp=1,6
          do i=1,4
            y(i,istp)=(0.d0,0.d0)
            do j=1,2
              y(i,istp)=y(i,istp)+b(j,istp)*y0(i,j)
            enddo
          enddo
        enddo
      else
        do istp=1,6
          do i=1,4
            y(i,istp)=(0.d0,0.d0)
            do j=1,2
              y(i,istp)=y(i,istp)+b(j+2,istp)*y0(i,j)
            enddo
          enddo
        enddo
      endif
c
      if(ipath.eq.1)then
        n=nno(lpath)
        ylw(1,1)=-kp(n)
        ylw(2,1)=wa(n)
        ylw(3,1)=ck
        ylw(4,1)=-wb(n)*kp(n)
c
        ylw(1,2)=ck
        ylw(2,2)=-wb(n)*ks(n)
        ylw(3,2)=-ks(n)
        ylw(4,2)=wa(n)
c
        if(lpath.gt.ls.and.lpath.eq.lzr)call cmemcpy(ylw,y0,8)
        do l=lpath-1,ls,-1
          ch0=dcmplx(hp(l),0.d0)
          n=nno(l)
c
c         determination of propagation matrix
c
          call qsve2am(n,ck,ylw,c0,2,rsite)
          pwave=cdexp(-kp(n)*ch0)
          swave=cdexp(-ks(n)*ch0)
c
c         orthonormalization of the p-sv modes
c
          cfac=(1.d0,0.d0)/(c0(4,2)*c0(2,1)-c0(2,2)*c0(4,1))
          orth(1,1)=c0(4,2)*cfac
          orth(1,2)=-c0(2,2)*cfac
          orth(2,1)=-c0(4,1)*cfac
          orth(2,2)=c0(2,1)*cfac
c
          call caxcb(c0,orth,4,2,2,c1)
          if(l.lt.lzr)then
c
c           additional normalization to avoid overflow
c
            do i=1,2
              orth(i,1)=orth(i,1)*pwave
              orth(i,2)=orth(i,2)*swave
            enddo
            call caxcb(y0,orth,4,2,2,y1)
            call cmemcpy(y1,y0,8)
          endif
c
c         for partial solution only!
c
          if(.not.pup(l))then
            c1(1,1)=(0.d0,0.d0)
            c1(1,2)=(0.d0,0.d0)
            if(l+1.eq.lzr)then
              call qsve2am(n,ck,y0,c0,2,rsite)
              c0(1,1)=(0.d0,0.d0)
              c0(1,2)=(0.d0,0.d0)
              call qsam2ve(n,ck,y0,c0,2,rsite)
            endif
          endif
          if(.not.svup(l))then
            c1(3,1)=(0.d0,0.d0)
            c1(3,2)=(0.d0,0.d0)
            if(l+1.eq.lzr)then
              call qsve2am(n,ck,y0,c0,2,rsite)
              c0(3,1)=(0.d0,0.d0)
              c0(3,2)=(0.d0,0.d0)
              call qsam2ve(n,ck,y0,c0,2,rsite)
            endif
          endif
          if(.not.pdw(l))then
            c1(1,1)=(0.d0,0.d0)
            c1(3,1)=(0.d0,0.d0)
            if(l.lt.lzr)then
              do j=1,4
                  y0(j,1)=(0.d0,0.d0)
                enddo
              endif
            endif
          if(.not.svdw(l))then
            c1(1,2)=(0.d0,0.d0)
            c1(3,2)=(0.d0,0.d0)
            if(l.lt.lzr)then
              do j=1,4
                y0(j,2)=(0.d0,0.d0)
             enddo
            endif
          endif
c
c         end partial solution procedure!
c
          c1(1,1)=c1(1,1)*pwave*pwave
c          c1(2,1)=c1(2,1)
          c1(3,1)=c1(3,1)*pwave*swave
          c1(4,1)=(0.d0,0.d0)
c
          c1(1,2)=c1(1,2)*swave*pwave
          c1(2,2)=(0.d0,0.d0)
          c1(3,2)=c1(3,2)*swave*swave
c          c1(4,2)=c1(4,2)
c
          call qsam2ve(n,ck,ylw,c1,2,rsite)
          if(l.gt.ls.and.l.eq.lzr)call cmemcpy(ylw,y0,8)
        enddo
        do i=1,4
          do j=1,2
            coef(i,j)=yup(i,j)
            coef(i,j+2)=-ylw(i,j)
          enddo
        enddo
        key=0
        call cdgemp(coef,b0,4,6,1.d-30,key)
        if(key.eq.0)then
          print *,'warning in qspsv: anormal exit from cdgemp!'
          return
        endif
        if(lzr.le.ls)then
          do istp=1,6
            do i=1,4
              do j=1,2
                y(i,istp)=y(i,istp)-b0(j,istp)*y0(i,j)
              enddo
            enddo
          enddo
        else
          do istp=1,6
            do i=1,4
              do j=1,2
                y(i,istp)=y(i,istp)-b0(j+2,istp)*y0(i,j)
              enddo
            enddo
          enddo
        endif
      endif
c
      if(n0rs.gt.0)then
c
c       for receiver-site structure different from source-site structure
c
        n=nno(lzr)
        do istp=1,6
          call qsve2am(n,ck,y(1,istp),cinc(1,istp),1,rsite)
        enddo
c
c       determination of starting upper sublayer
c
        rsite=.true.
c
        if(isurf.eq.0)then
          do j=1,2
            do i=1,4
              yup(i,j)=(0.d0,0.d0)
            enddo
          enddo
          yup(1,1)=(1.d0,0.d0)
          if(vpair.gt.0.d0)yup(2,1)=-accair/kpair
          yup(3,2)=(1.d0,0.d0)
        else
          n=nnors(lup)
c
          yup(1,1)=kprs(n)
          yup(2,1)=wars(n)
          yup(3,1)=ck
          yup(4,1)=wbrs(n)*kprs(n)
c
          yup(1,2)=ck
          yup(2,2)=wbrs(n)*ksrs(n)
          yup(3,2)=ksrs(n)
          yup(4,2)=wars(n)
        endif
        if(lzrrs.eq.1)call cmemcpy(yup,y0,8)
c
        do l=2,lprs
          ch0=dcmplx(hprs(l-1),0.d0)
          n=nnors(l-1)
c
c
c         determination of propagation matrix
c
          call qsve2am(n,ck,yup,c0,2,rsite)
          pwave=cdexp(-kprs(n)*ch0)
          swave=cdexp(-ksrs(n)*ch0)
c
c         orthonormalization of the p-sv modes
c
          cfac=(1.d0,0.d0)/(c0(3,2)*c0(1,1)-c0(1,2)*c0(3,1))
          orth(1,1)=c0(3,2)*cfac
          orth(1,2)=-c0(1,2)*cfac
          orth(2,1)=-c0(3,1)*cfac
          orth(2,2)=c0(1,1)*cfac
c
          call caxcb(c0,orth,4,2,2,c1)
          if(l.gt.lzrrs)then
c
c           additional normalization to avoid overflow
c
            do i=1,2
              orth(i,1)=orth(i,1)*pwave
              orth(i,2)=orth(i,2)*swave
            enddo
            call caxcb(y0,orth,4,2,2,y1)
            call cmemcpy(y1,y0,8)
          endif
c          c1(1,1)=c1(1,1)
          c1(2,1)=c1(2,1)*pwave*pwave
          c1(3,1)=(0.d0,0.d0)
          c1(4,1)=c1(4,1)*pwave*swave
c
          c1(1,2)=(0.d0,0.d0)
          c1(2,2)=c1(2,2)*swave*pwave
c          c1(3,2)=c1(3,2)
          c1(4,2)=c1(4,2)*swave*swave
c
          call qsam2ve(n,ck,yup,c1,2,rsite)
          if(l.eq.lzrrs)call cmemcpy(yup,y0,8)
        enddo
        n=nnors(lprs)
        call qsve2am(n,ck,yup,c0,2,rsite)
        do istp=1,6
          brs(1,istp)=cinc(1,istp)
          brs(2,istp)=cinc(3,istp)
        enddo
        coefrs(1,1)=c0(1,1)
        coefrs(1,2)=c0(1,2)
        coefrs(2,1)=c0(3,1)
        coefrs(2,2)=c0(3,2)
        key=0
        call cdgemp(coefrs,brs,2,6,0.d0,key)
        if(key.eq.0)then
          print *,'warning in qspsv: anormal exit from cdgemp!'
          return
        endif
        do istp=1,6
          do i=1,4
            y(i,istp)=(0.d0,0.d0)
            do j=1,2
              y(i,istp)=y(i,istp)+brs(j,istp)*y0(i,j)
            enddo
          enddo
        enddo
      endif
      if(isurf.eq.2)then
c
c       free surface corrections
c
        n=1
        if(rsite)then
          yup(1,1)=kprs(n)
          yup(2,1)=wars(n)
          yup(3,1)=ck
          yup(4,1)=wbrs(n)*kprs(n)
c
          yup(1,2)=ck
          yup(2,2)=wbrs(n)*ksrs(n)
          yup(3,2)=ksrs(n)
          yup(4,2)=wars(n)
c
          ylw(1,1)=-kprs(n)
          ylw(2,1)=wars(n)
          ylw(3,1)=ck
          ylw(4,1)=-wbrs(n)*kprs(n)
c
          ylw(1,2)=ck
          ylw(2,2)=-wbrs(n)*ksrs(n)
          ylw(3,2)=-ksrs(n)
          ylw(4,2)=wars(n)
        else
          yup(1,1)=kp(n)
          yup(2,1)=wa(n)
          yup(3,1)=ck
          yup(4,1)=wb(n)*kp(n)
c
          yup(1,2)=ck
          yup(2,2)=wb(n)*ks(n)
          yup(3,2)=ks(n)
          yup(4,2)=wa(n)
c
          ylw(1,1)=-kp(n)
          ylw(2,1)=wa(n)
          ylw(3,1)=ck
          ylw(4,1)=-wb(n)*kp(n)
c
          ylw(1,2)=ck
          ylw(2,2)=-wb(n)*ks(n)
          ylw(3,2)=-ks(n)
          ylw(4,2)=wa(n)
        endif
c
        do istp=1,6
          call qsve2am(1,ck,y(1,istp),b0(1,istp),1,rsite)

          b(1,istp)=-yup(2,1)*b0(1,istp)-yup(2,2)*b0(3,istp)
          b(2,istp)=-yup(4,1)*b0(1,istp)-yup(4,2)*b0(3,istp)
          orth(1,1)=ylw(2,1)
          orth(1,2)=ylw(2,2)
          orth(2,1)=ylw(4,1)
          orth(2,2)=ylw(4,2)
          key=0
          call cdgemp(orth,b(1,istp),2,1,0.d0,key)
          if(key.eq.0)then
            print *,'warning in qspsv: anormal exit from cdgemp!'
            return
          endif
          b0(2,istp)=b(1,istp)
          b0(4,istp)=b(2,istp)
          call qsam2ve(1,ck,y(1,istp),b0(1,istp),1,rsite)
        enddo
      endif
      return
      end
