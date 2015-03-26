      subroutine qssh(y,k,lup,llw)
      implicit none
c
c     calculation of response to sh source
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
      integer i,istp,l,n,key
      double complex cnorm
      double complex y0(2),b(2,6),b0(2,6)
      double complex y1(2),yup(2),ylw(2)
      double complex hk(2,2,nzmax),coef(2,2),cinc(6),c0
      logical rsite
c
c===============================================================================
c
c     matrix propagation from surface to source
c
      rsite=.false.
      do l=lup,ls-1
        n=nno(l)
        call qshksh(hk(1,1,l),hp(l),n,rsite)
      enddo
      do l=ls,llw-1
        n=nno(l)
        call qshksh(hk(1,1,l),-hp(l),n,rsite)
      enddo
c
c     yup: the starting solution vector
c
      yup(1)=(1.d0,0.d0)
      if(lup.eq.1.and.isurf.eq.0)then
        yup(2)=(0.d0,0.d0)
      else
        n=nno(lup)
        yup(2)=cmu(n)*ks(n)
      endif
      if(lup.eq.lzr)call cmemcpy(yup,y0,2)
c
      do l=lup+1,ls
        n=nno(l-1)
c
c       determination of propagation matrix
c
        call caxcb(hk(1,1,l-1),yup,2,2,1,y1)
        call cmemcpy(y1,yup,2)
c
c       normalization
c
        if(l.gt.lzr)then
          cnorm=cdexp(-ks(n)*dcmplx(hp(l-1),0.d0))
          y0(1)=y0(1)*cnorm
          y0(2)=y0(2)*cnorm
        else if(l.eq.lzr)then
          call cmemcpy(yup,y0,2)
        endif
      enddo
c
c===============================================================================
c
c     matrix propagation from half-space to source
c
c     ylw: the starting solution vector
c
      n=nno(llw)
      ylw(1)=(1.d0,0.d0)
      if(sh(llw))then
        ylw(2)=-cmu(n)*ks(n)
      else
        ylw(2)=(0.d0,0.d0)
      endif
      if(llw.gt.ls.and.llw.eq.lzr)call cmemcpy(ylw,y0,2)
c
      do l=llw-1,ls,-1
        n=nno(l)
c
c       determination of propagation matrix
c
        call caxcb(hk(1,1,l),ylw,2,2,1,y1)
        call cmemcpy(y1,ylw,2)
        if(l.lt.lzr)then
          cnorm=cdexp(-ks(n)*dcmplx(hp(l),0.d0))
          y0(1)=y0(1)*cnorm
          y0(2)=y0(2)*cnorm
        else if(l.gt.ls.and.l.eq.lzr)then
          call cmemcpy(ylw,y0,2)
        endif
      enddo
c
c===============================================================================
c     source function
c===============================================================================
c
      do istp=1,6
        do i=1,2
          b(i,istp)=sfct0(i+4,istp)+dcmplx(k,0.d0)*sfct1(i+4,istp)
        enddo
      enddo
      do i=1,2
        coef(i,1)=yup(i)
        coef(i,2)=-ylw(i)
      enddo
      if(ipath.eq.1)call cmemcpy(b,b0,12)
      key=0
      call cdgemp(coef,b,2,6,0.d0,key)
      if(key.eq.0)then
        print *,'warning in qssh: anormal exit from cdgemp!'
        return
      endif
      if(lzr.le.ls)then
        do istp=1,6
          do i=1,2
            y(i+4,istp)=b(1,istp)*y0(i)
          enddo
        enddo
      else
        do istp=1,6
          do i=1,2
            y(i+4,istp)=b(2,istp)*y0(i)
          enddo
        enddo
      endif
      if(ipath.eq.1)then
        n=nno(lpath)
        ylw(1)=(1.d0,0.d0)
        if(sh(lpath))then
          ylw(2)=-cmu(n)*ks(n)
        else
          ylw(2)=(0.d0,0.d0)
        endif
        if(lpath.gt.ls.and.lpath.eq.lzr)call cmemcpy(ylw,y0,2)
        do l=lpath-1,ls,-1
          n=nno(l)
c
c         determination of propagation matrix
c
          call caxcb(hk(1,1,l),ylw,2,2,1,y1)
          call cmemcpy(y1,ylw,2)
c
          if(l.lt.lzr)then
            cnorm=cdexp(-ks(n)*dcmplx(hp(l),0.d0))
            y0(1)=y0(1)*cnorm
            y0(2)=y0(2)*cnorm
          else if(l.gt.ls.and.l.eq.lzr)then
            call cmemcpy(ylw,y0,2)
          endif
        enddo
        do i=1,2
          coef(i,1)=yup(i)
          coef(i,2)=-ylw(i)
        enddo
        key=0
        call cdgemp(coef,b0,2,6,0.d0,key)
        if(key.eq.0)then
          print *,'warning in qssh: anormal exit from cdgemp!'
          return
        endif
        if(lzr.le.ls)then
          do istp=1,6
            do i=1,2
              y(i+4,istp)=y(i+4,istp)-b0(1,istp)*y0(i)
            enddo
          enddo
        else
          do istp=1,6
            do i=1,2
              y(i+4,istp)=y(i+4,istp)-b0(2,istp)*y0(i)
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
          cinc(istp)=(0.5d0,0.d0)*(y(5,istp)+y(6,istp)/(cmu(n)*ks(n)))
        enddo
        rsite=.true.
        do l=1,lprs-1
          n=nnors(l)
          call qshksh(hk(1,1,l),hprs(l),n,rsite)
        enddo
c
c       yup: the starting solution vector
c
        yup(1)=(1.d0,0.d0)
        if(isurf.eq.0)then
          yup(2)=(0.d0,0.d0)
        else
          yup(2)=cmurs(1)*ksrs(1)
        endif
        if(lzrrs.eq.1)call cmemcpy(yup,y0,2)
c
        do l=2,lprs
          n=nnors(l-1)
c
c         determination of propagation matrix
c
          call caxcb(hk(1,1,l-1),yup,2,2,1,y1)
          call cmemcpy(y1,yup,2)
c
c         normalization
c
          if(l.gt.lzrrs)then
            cnorm=cdexp(-ksrs(n)*dcmplx(hprs(l-1),0.d0))
            y0(1)=y0(1)*cnorm
            y0(2)=y0(2)*cnorm
          else if(l.eq.lzrrs)then
            call cmemcpy(yup,y0,2)
          endif
        enddo
        n=nnors(lprs)
        c0=(0.5d0,0.d0)*(yup(1)+yup(2)/(cmurs(n)*ksrs(n)))
        do istp=1,6
          y(5,istp)=y0(1)*cinc(istp)/c0
          y(6,istp)=y0(2)*cinc(istp)/c0
        enddo
      endif
      if(isurf.eq.2)then
c
c       free surface corrections
c
        do istp=1,6
          y(5,istp)=(2.d0,0.d0)*y(5,istp)
          y(6,istp)=(0.d0,0.d0)
        enddo
      endif
      return
      end
