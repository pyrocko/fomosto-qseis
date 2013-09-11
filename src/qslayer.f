      subroutine qslayer(ierr)
      implicit none
      integer ierr
c
      include 'qsglobal.h'
c
      integer l,n,li,lp0
      double precision zswap
      double precision z(nzmax),z0(nzmax)
c
      lp0=1
      z0(lp0)=0.d0
      do n=1,n0-1
        lp0=lp0+1
        z0(lp0)=z0(lp0-1)+h(n)
      enddo
      lp0=lp0+1
      z0(lp0)=zs
      lp0=lp0+1
      z0(lp0)=zr
c
c     sort the z0-profile
c
      do l=1,lp0-1
        do li=l+1,lp0
          if(z0(li).lt.z0(l))then
            zswap=z0(l)
            z0(l)=z0(li)
            z0(li)=zswap
          endif
        enddo
      enddo
c
c     delete duplicates
c
      lp=1
      z(lp)=0.d0
      do l=2,lp0
        if(z0(l).gt.z(lp))then
          hp(lp)=z0(l)-z(lp)
          lp=lp+1
          z(lp)=z0(l)
        endif
      enddo
      hp(lp)=0.d0
c
c     determine ls,lzr
c
      do l=1,lp
        if(z(l).eq.zs)ls=l
        if(z(l).eq.zr)lzr=l
      enddo
c
c     determine layer no of each depth
c
      li=1
      zswap=h(1)
      nno(1)=1
      do l=2,lp
        if(z(l).ge.zswap.and.li.lt.n0)then
          li=li+1
          zswap=zswap+h(li)
        endif
        nno(l)=li
      enddo
c
c     for the receiver site
c
      if(l0rs.le.0)then
        n0rs=0
      else
        lp0=1
        z0(lp0)=0.d0
        do n=1,n0rs-1
          lp0=lp0+1
          z0(lp0)=z0(lp0-1)+hrs(n)
        enddo
        lp0=lp0+1
        z0(lp0)=zrrs
c
c       sort the z0-profile
c
        do l=1,lp0-1
          do li=l+1,lp0
            if(z0(li).lt.z0(l))then
              zswap=z0(l)
              z0(l)=z0(li)
              z0(li)=zswap
            endif
          enddo
        enddo
c
c       delete duplicates
c
        lprs=1
        z(lprs)=0.d0
        do l=2,lp0
          if(z0(l).gt.z(lprs))then
            hprs(lprs)=z0(l)-z(lprs)
            lprs=lprs+1
            z(lprs)=z0(l)
          endif
        enddo
        hprs(lprs)=0.d0
c
c       determine lzrrs
c
        do l=1,lprs
          if(z(l).eq.zrrs)lzrrs=l
        enddo
c
c       determine layer no of each depth
c
        li=1
        zswap=h(1)
        nnors(1)=1
        do l=2,lprs
          if(z(l).ge.zswap.and.li.lt.n0rs)then
            li=li+1
            zswap=zswap+hrs(li)
          endif
          nnors(l)=li
        enddo
      endif
c
      ierr=0
      return
      end
