      subroutine qsgetinp(unit,srate,nssel)
      implicit none
      integer unit,nssel
      double precision srate
c
      include 'qsglobal.h'
c
c     work space
c
      integer i,ir,istp,j,l,l1,lrs,lcut,n,ierr,iazi
      integer irlast,irnow,is,ns
      integer ieqdis,kmordeg,iv0,flen0,nup,nlw
      double precision rr,z,pi,taunorm,rnow,rlast,tnow,tlast
      double precision v00,depth,hpmin,vsliquid
      double precision s1,s2,s,ds,smin,shead,twinmin
      double precision r1,r2,dr,dm,t1,rdis
      double precision mis,mcl,mdc,st,di,ra,deg2rad
      double precision suppress,ros,vps,vss,fcut
      double precision rot(3,3),sm(3,3),swap(3,3)
      double precision resolut(3),t0(nrmax)
      character*80 outfile0(7)
c
c     source parameters
c     =================
c
      pi=4.d0*datan(1.d0)
      deg2rad=pi/180.d0
      call skip_comments(unit)
      read(unit,*)zs
      zs=dmax1(0.d0,km2m*zs)
c
c     receiver parameters
c     ===================
c
      call skip_comments(unit)
      read(unit,*)zr
      zr=km2m*zr
      call skip_comments(unit)
      read(unit,*)ieqdis,kmordeg
      call skip_comments(unit)
      read(unit,*)nr
      if(nr.gt.nrmax)then
        stop 'Error in input: nr > nrmax!'  
      endif
      if(kmordeg.eq.1)then
        dm=km2m
      else
        dm=rr0*pi/180.d0
      endif
      if(ieqdis.eq.1.and.nr.gt.1)then
        read(unit,*)r1,r2
        if(nr.eq.1)then
          dr=0.d0
        else
          dr=(r2-r1)/dble(nr-1)
        endif
        do i=1,nr
          r(i)=dm*(r1+dr*dble(i-1))
        enddo
      else
        read(unit,*)(r(i),i=1,nr)
        do i=1,nr
          r(i)=dm*r(i)
        enddo
      endif
c
      call skip_comments(unit)
      read(unit,*)tstart,twindow,nt
      if(twindow.le.0.d0.or.nt.le.0)then
        stop 'Error in input: time window or sampling no <= 0!'  
      endif
c
      call skip_comments(unit)
      read(unit,*)iv0,v0
      if(iv0.eq.1)then
        v0=km2m*v0
      else if(v0.gt.0.d0)then
        v0=rr0*pi/180.d0/v0
        write(*,'(a,f10.4,a)')' Velocity for time reduction: ',
     &                        v0/km2m,' km/s'
      else
        v0=0.d0
      endif
c
c     wavenumber integration parameters
c     =================================
c
      call skip_comments(unit)
      read(unit,*)ndtrans
      if(ndtrans.lt.0.or.ndtrans.gt.ndtransmax)then
        stop 'Error in input: wrong select of integration algorithm!'
      endif
      call skip_comments(unit)
      read(unit,*)(slw(j),j=1,4)
      do j=1,4
        slw(j)=slw(j)/km2m
      enddo
      if(slw(1).lt.0.d0.or.slw(2).lt.0.d0.or.
     +   slw(3).le.0.d0.or.slw(4).le.0.d0.or.
     +   slw(2).lt.slw(1).or.slw(3).lt.slw(2).or.
     +   slw(4).lt.slw(3))then
        fullwave=.true.
      else
        fullwave=.false.
      endif
      call skip_comments(unit)
      read(unit,*)srate
      if(srate.lt.1.d0)srate=1.d0
c
      call skip_comments(unit)
      read(unit,*)suppress
      if(suppress.le.0.d0.or.suppress.ge.1.d0)then
        suppress=dexp(-1.d0)
        print *,'warning in qsmain: aliasing suppression'
        print *,'factor is replaced by the default value of 1/e.'
      endif
      fi=dlog(suppress)/(2.d0*pi*twindow)
c
c     partial solution parameters
c     ===========================
c
      call skip_comments(unit)
      read(unit,*)isurf
      if(isurf.lt.0.or.isurf.gt.2)then
        stop 'Error: wrong switch for filtering surface reflection!'
      else if(isurf.eq.2.and.zr.gt.0.d0)then
        stop 'Error: filtering surface multiples for zr > 0!'
      endif
      call skip_comments(unit)
      read(unit,*)ipath,pathdepth
      pathdepth=pathdepth*km2m
      if(ipath.eq.1.and.(pathdepth.lt.zs.or.pathdepth.lt.zr))then
        print *,'warning: condition for path filter is not satisfied,'
        print *,'==> path filter will not be selected!'
        ipath=0
      endif
      call skip_comments(unit)
      read(unit,*)npar
      if(npar.ge.1)then
        ipartial=1
        do i=1,npar
          call skip_comments(unit)
          read(unit,*)zup(i),zlow(i),ipsv(i)
          if(ipsv(i).le.0.or.ipsv(i).ge.5)then
            stop ' Error in qsmain: wrong partial solution selection!'
          endif
          zup(i)=zup(i)*km2m
          zlow(i)=zlow(i)*km2m
        enddo
      endif
c
c     wavelet parameters
c     ==================
c
      call skip_comments(unit)
      read(unit,*)taunorm,wdeg
      if(wdeg.lt.0.or.wdeg.gt.2)then
        stop ' Error in qsmain: wrong wavelet selection!'
      else if(wdeg.eq.0)then
        call skip_comments(unit)
        read(unit,*)nn0
        read(unit,*)(wv0(i),i=1,nn0)
      endif
c
c     seimometer parameters
c     =====================
c
      call skip_comments(unit)
      read(unit,*)asm
      call skip_comments(unit)
      read(unit,*)nroot
      read(unit,*)(root(i),i=1,nroot)
      call skip_comments(unit)
      read(unit,*)npole
      read(unit,*)(pole(i),i=1,npole)
c
c     output files
c     ============
c
      varbtxt='U'
      call skip_comments(unit)
      read(unit,*)(ssel(istp),istp=1,6)
      call skip_comments(unit)
      read(unit,*)(outfile0(istp),istp=1,6)
      do istp=1,6
        if(ssel(istp).ne.1)ssel(istp)=0
        do flen0=80,1,-1
          if(outfile0(istp)(flen0:flen0).ne.' ')goto 100
        enddo
100     continue
        outfile(1,istp)=outfile0(istp)(1:flen0)//'.tz'
        outfile(2,istp)=outfile0(istp)(1:flen0)//'.tr'
        outfile(3,istp)=outfile0(istp)(1:flen0)//'.tt'
        outfile(4,istp)=outfile0(istp)(1:flen0)//'.tv'
        do i=1,4
          flen(i,istp)=flen0+3
        enddo
      enddo
      call skip_comments(unit)
      read(unit,*)ssel(7)
      if (ssel(7) .ne. 0) then
          backspace(unit)
      endif
      if(ssel(7).eq.1)then
        read(unit,*)ssel(7),(mtensor(i),i=1,6),outfile0(7)
      else if(ssel(7).eq.2)then
        read(unit,*)ssel(7),mis,mcl,mdc,st,di,ra,outfile0(7)
        st=st*deg2rad
        di=di*deg2rad
        ra=ra*deg2rad
c
c       use principal stress coordinates:
c       x: along T-axis
c       y: along N-axis
c       z: along P-axis (symmetry axis of CLVD)
c
        do i=1,3
          do j=1,3
            sm(i,j)=0.d0
          enddo
        enddo
        sm(1,1)=mis-0.5d0*mcl+mdc
        sm(2,2)=mis-0.5d0*mcl
        sm(3,3)=mis+mcl-mdc
c
c       construct the rotation matrix:
c       1. around x by -45 deg;
c       2. around z by rake angle
c       3. around x -dip angle;
c       4. around z by -strike angle.
c
        rot(1,1)=(dcos(st)*dcos(ra)
     &           +dsin(st)*(dcos(di)*dsin(ra)-dsin(di)))/dsqrt(2.d0)
        rot(1,2)=dcos(st)*dsin(ra)-dsin(st)*dcos(di)*dcos(ra)
        rot(1,3)=(dcos(st)*dcos(ra)
     &           +dsin(st)*(dcos(di)*dsin(ra)+dsin(di)))/dsqrt(2.d0)
        rot(2,1)=(dsin(st)*dcos(ra)
     &           -dcos(st)*(dcos(di)*dsin(ra)-dsin(di)))/dsqrt(2.d0)
        rot(2,2)=dsin(st)*dsin(ra)+dcos(st)*dcos(di)*dcos(ra)
        rot(2,3)=(dsin(st)*dcos(ra)
     &           -dcos(st)*(dcos(di)*dsin(ra)+dsin(di)))/dsqrt(2.d0)
        rot(3,1)=(-dsin(di)*dsin(ra)-dcos(di))/dsqrt(2.d0)
        rot(3,2)=dsin(di)*dcos(ra)
        rot(3,3)=(-dsin(di)*dsin(ra)+dcos(di))/dsqrt(2.d0)
c
        do i=1,3
          do j=1,3
            swap(i,j)=0.d0
            do l=1,3
              swap(i,j)=swap(i,j)+rot(i,l)*sm(l,j)
            enddo
          enddo
        enddo
        do i=1,3
          do j=1,3
            sm(i,j)=0.d0
            do l=1,3
              sm(i,j)=sm(i,j)+swap(i,l)*rot(j,l)
            enddo
          enddo
        enddo
        mtensor(1)=sm(1,1)
        mtensor(2)=sm(2,2)
        mtensor(3)=sm(3,3)
        mtensor(4)=sm(1,2)
        mtensor(5)=sm(2,3)
        mtensor(6)=sm(3,1)
      else
        do i=1,6
          mtensor(i)=0.d0
        enddo
        ssel(7)=0
      endif
      call skip_comments(unit)
      read(unit,*)iazi
      if(iazi.eq.0)then
        read(unit,*)azimuth(1)
        do i=2,nr
          azimuth(i)=azimuth(1)
        enddo
      else
        read(unit,*)(azimuth(i),i=1,nr)
      endif
      do flen0=80,1,-1
        if(outfile0(7)(flen0:flen0).ne.' ')goto 200
      enddo
200   continue
      outfile(1,7)=outfile0(7)(1:flen0)//'.tz'
      outfile(2,7)=outfile0(7)(1:flen0)//'.tr'
      outfile(3,7)=outfile0(7)(1:flen0)//'.tt'
      outfile(4,7)=outfile0(7)(1:flen0)//'.tv'
      do i=1,4
        flen(i,7)=flen0+3
      enddo
c
      nssel=0
      do istp=1,7
        nssel=nssel+ssel(istp)
      enddo
      if(nssel.le.0)then
        stop ' Error in input file: no outputs selected!'
      endif
c
c     global model parameters
c     =======================
c
      call skip_comments(unit)
      read(unit,*)iflat
      call skip_comments(unit)
      read(unit,*)(resolut(i),i=1,3)
      do i=1,3
        if(resolut(i).le.0.d0)resolut(i)=0.1d0
        resolut(i)=1.d-02*resolut(i)
      enddo
      call skip_comments(unit)
      read(unit,*)l
      if(l.gt.lmax)then
        stop ' Error in input: to large number of layers!'
      endif
c
c     multilayered model parameters
c     =============================
c
      do i=1,l
        call skip_comments(unit)
        read(unit,*)j,h(i),vp(i),vs(i),ro(i),qp(i),qs(i)
c
c       input units:    -,km,  km/s, km/s, g/cm^3,-,-
c
        h(i)=km2m*h(i)
        vp(i)=km2m*vp(i)
        vs(i)=km2m*vs(i)
        ro(i)=km2m*ro(i)
        if(vs(i).le.vspmin*vp(i))vs(i)=0.9d0*vspmin*vp(i)
      enddo
c
      call skip_comments(unit)
      read(unit,*)lrs
      if(lrs.gt.lmax)then
        stop ' Error in input: to large number of layers!'
      endif
c
c     multilayered model parameters
c     =============================
c
      do i=1,lrs
        call skip_comments(unit)
        read(unit,*)j,hrs(i),vprs(i),vsrs(i),rors(i),qprs(i),qsrs(i)
c
c       input units:    -,km,  km/s, km/s, g/cm^3,-,-
c
        hrs(i)=km2m*hrs(i)
        vprs(i)=km2m*vprs(i)
        vsrs(i)=km2m*vsrs(i)
        rors(i)=km2m*rors(i)
        if(vsrs(i).le.vspmin*vprs(i))vsrs(i)=0.9d0*vspmin*vprs(i)
      enddo
      zrrs=zr
      if(lrs.gt.0)zr=hrs(lrs)
c
c     end of inputs
c     =============
c
      if(iflat.eq.1)then
c
c       flat earth transformation (Mueller, 1985)
c
        zs=rr0*dlog(rr0/(rr0-zs))
        zr=rr0*dlog(rr0/(rr0-zr))
        zrrs=rr0*dlog(rr0/(rr0-zrrs))
        if(ipartial.eq.1)then
          do i=1,npar
            zup(i)=rr0*dlog(rr0/(rr0-zup(i)))
            zlow(i)=rr0*dlog(rr0/(rr0-zlow(i)))
          enddo
        endif
        if(ipath.eq.1)then
          pathdepth=rr0*dlog(rr0/(rr0-pathdepth))
        endif
c
        do i=1,l
          rr=rr0-h(i)
          h(i)=rr0*dlog(rr0/rr)
          vp(i)=vp(i)*rr0/rr
          vs(i)=vs(i)*rr0/rr
          ro(i)=ro(i)*(rr/rr0)**ndens
        enddo
c
        do i=1,lrs
          rr=rr0-hrs(i)
          hrs(i)=rr0*dlog(rr0/rr)
          vprs(i)=vprs(i)*rr0/rr
          vsrs(i)=vsrs(i)*rr0/rr
          rors(i)=rors(i)*(rr/rr0)**ndens
        enddo
      endif
c
c     end of the flat earth transformation
c
      dt=twindow/dble(nt-1)
      nf=1
300   nf=2*nf
      if(nf.lt.nt)goto 300
      nf=nf/2
      if(nf.gt.nfmax)then
        print *,'Error in input: time sampling no > ',2*nfmax,'!'
        stop
      endif
      df=1.d0/(dble(2*nf)*dt)
      fcut=0.5d0/dt
      tau=taunorm*dt
      if(taunorm.le.0.d0)tau=2.d0*dt
c
      comptxt(1)='z'
      comptxt(2)='r'
      comptxt(3)='t'
      comptxt(4)='v'
      do j=1,nr
        i=j/1000
        rcvtxt(j)(1:1)=char(ichar('0')+i)
        i=mod(j,1000)/100
        rcvtxt(j)(2:2)=char(ichar('0')+i)
        i=mod(j,100)/10
        rcvtxt(j)(3:3)=char(ichar('0')+i)
        i=mod(j,10)
        rcvtxt(j)(4:4)=char(ichar('0')+i)
      enddo
c
c     determine upper und lower parameter values of each layer
c
      l0=1
      z1(l0)=0.d0
      do i=2,l
        if(h(i).gt.h(i-1))then
          z1(l0)=h(i-1)
          vp1(l0)=vp(i-1)
          vs1(l0)=vs(i-1)
          ro1(l0)=ro(i-1)
          qp1(l0)=qp(i-1)
          qs1(l0)=qs(i-1)
c
          z2(l0)=h(i)
          vp2(l0)=vp(i)
          vs2(l0)=vs(i)
          ro2(l0)=ro(i)
          qp2(l0)=qp(i)
          qs2(l0)=qs(i)
          l0=l0+1
        else
          z1(l0)=h(i)
          vp1(l0)=vp(i)
          vs1(l0)=vs(i)
          ro1(l0)=ro(i)
          qp1(l0)=qp(i)
          qs1(l0)=qs(i)
        endif
      enddo
      z1(l0)=h(l)
      vp1(l0)=vp(l)
      vs1(l0)=vs(l)
      ro1(l0)=ro(l)
      qp1(l0)=qp(l)
      qs1(l0)=qs(l)
c
c     determine upper und lower parameter values of each layer (receiver site)
c
      if(lrs.le.0)then
        l0rs=0
      else
        l0rs=1
        z1rs(l0rs)=0.d0
        do i=2,lrs
          if(hrs(i).gt.hrs(i-1))then
            z1rs(l0rs)=hrs(i-1)
            vp1rs(l0rs)=vprs(i-1)
            vs1rs(l0rs)=vsrs(i-1)
            ro1rs(l0rs)=rors(i-1)
            qp1rs(l0rs)=qprs(i-1)
            qs1rs(l0rs)=qsrs(i-1)
c
            z2rs(l0rs)=hrs(i)
            vp2rs(l0rs)=vprs(i)
            vs2rs(l0rs)=vsrs(i)
            ro2rs(l0rs)=rors(i)
            qp2rs(l0rs)=qprs(i)
            qs2rs(l0rs)=qsrs(i)
            l0rs=l0rs+1
          else
            z1rs(l0rs)=hrs(i)
            vp1rs(l0rs)=vprs(i)
            vs1rs(l0rs)=vsrs(i)
            ro1rs(l0rs)=rors(i)
            qp1rs(l0rs)=qprs(i)
            qs1rs(l0rs)=qsrs(i)
          endif
        enddo
        z1rs(l0rs)=hrs(lrs)
        vp1rs(l0rs)=vprs(lrs)
        vs1rs(l0rs)=vsrs(lrs)
        ro1rs(l0rs)=rors(lrs)
        qp1rs(l0rs)=qprs(lrs)
        qs1rs(l0rs)=qsrs(lrs)
c
        ipath=1
        pathdepth=dmax1(pathdepth,zr)
      endif
c
c     construction of sublayers at the cutoff frequency
c
      call qssublay(resolut,fcut)
      write(*,*)' The layered model of source site:'
      write(*,'(7a)')'    no ','  z(km)  ',
     &               '  vp(km/s) ','  vs(km/s) ',' ro(g/cm^3)',
     &               '    qp   ','    qs'
      depth=0.d0
      do i=1,n0
        if(vs(i).le.vspmin*vp(i))then
          vsliquid=0.d0
        else
          vsliquid=vs(i)
        endif
        write(*,1000)i,depth/km2m,vp(i)/km2m,
     &               vsliquid/km2m,ro(i)/km2m,qp(i),qs(i)
        depth=depth+h(i)
        if(i.lt.n0)then
          write(*,1000)i,depth/km2m,vp(i)/km2m,
     &               vsliquid/km2m,ro(i)/km2m,qp(i),qs(i)
        endif
      enddo
      if(n0rs.gt.0)then
        write(*,*)' The layered model of receiver site:'
        write(*,'(7a)')'    no ','  z(km)  ',
     &               '  vp(km/s) ','  vs(km/s) ',' ro(g/cm^3)',
     &               '    qp   ','    qs'
        depth=0.d0
        do i=1,n0rs
          if(vsrs(i).le.vspmin*vprs(i))then
            vsliquid=0.d0
          else
            vsliquid=vsrs(i)
          endif
          write(*,1000)i,depth/km2m,vprs(i)/km2m,
     &               vsliquid/km2m,rors(i)/km2m,qprs(i),qsrs(i)
          depth=depth+hrs(i)
          if(i.lt.n0rs)then
            write(*,1000)i,depth/km2m,vprs(i)/km2m,
     &               vsliquid/km2m,rors(i)/km2m,qprs(i),qsrs(i)
          endif
        enddo
      endif
c
      call qslayer(ierr)
      n=nno(ls)
      ros=ro(n)
      vps=vp(n)
      vss=vs(n)
      if(iflat.eq.1)then
        rr=rr0*dexp(-zs/rr0)
        ros=ros*(rr0/rr)**ndens
        vps=vps*rr/rr0
        vss=vss*rr/rr0
      endif
      call qssource(ros,vps,vss)
c
      if(v0.gt.0.d0)then
        v00=1.d0/v0
      else
        v00=0.d0
      endif
      do ir=1,nr
        t0(ir)=tstart+r(ir)*v00
      enddo
c
      hpmin=hp(min0(ls,lzr))
      smin=1.d0/vp(nno(min0(ls,lzr)))
      do l=min0(ls,lzr),max0(ls,lzr)-1
        if(smin.gt.1.d0/vp(nno(l)))then
          smin=1.d0/vp(nno(l))
          hpmin=hp(l)
        endif
      enddo
c
c     compare direct p, reflected p and head wave phase
c
      lcut=0
      do l=max0(ls,lzr),lp
c
        twinmin=twindow
c
c       1. direct or reflected p wave
c
        irlast=0
        rlast=0.d0
        tlast=0.d0
        do l1=min0(ls,lzr),max0(ls,lzr)-1
          tlast=tlast+hp(l1)/vp(nno(l1))
        enddo
        do l1=max0(ls,lzr),l-1
          tlast=tlast+2.d0*hp(l1)/vp(nno(l1))
        enddo
        s1=0.d0
        s2=smin/dsqrt(1.d0+0.5d0*(hpmin/r(nr))**2)
        ns=2*nr+10
        ds=s2/dble(ns)
        do is=1,ns
          s=s1+dble(is)*ds
          rnow=0.d0
          tnow=0.d0
          do l1=min0(ls,lzr),max0(ls,lzr)-1
            rdis=hp(l1)*s/dsqrt(1.d0/vp(nno(l1))**2-s**2)
            rnow=rnow+rdis
            tnow=tnow+dsqrt(rdis**2+hp(l1)**2)/vp(nno(l1))
          enddo
          do l1=max0(ls,lzr),l-1
            rdis=hp(l1)*s/dsqrt(1.d0/vp(nno(l1))**2-s**2)
            rnow=rnow+2.d0*rdis
            tnow=tnow+2.d0*dsqrt(rdis**2+hp(l1)**2)/vp(nno(l1))
          enddo
          irnow=irlast
          do ir=irlast+1,nr
            if(r(ir).gt.rnow)then
              goto 400
            else if(r(ir).ge.rlast)then
              t1=(tlast*(rnow-r(ir))+tnow*(r(ir)-rlast))
     &          /(rnow-rlast)
              irnow=ir
              twinmin=dmin1(twinmin,t1-t0(ir))
            endif
          enddo
400       rlast=rnow
          tlast=tnow
          irlast=irnow
        enddo
c
c       2. head wave
c
        shead=1.d0/vp(nno(l))
        if(smin.gt.shead)then
          rlast=0.d0
          tlast=0.d0
          do l1=min0(ls,lzr),max0(ls,lzr)-1
            rdis=hp(l1)*shead/dsqrt(1.d0/vp(nno(l1))**2-shead**2)
            rlast=rlast+rdis
            tlast=tlast+dsqrt(rdis**2+hp(l1)**2)/vp(nno(l1))
          enddo
          do l1=max0(ls,lzr),l-1
            rdis=hp(l1)*shead/dsqrt(1.d0/vp(nno(l1))**2-shead**2)
            rlast=rlast+2.d0*rdis
            tlast=tlast+2.d0*dsqrt(rdis**2+hp(l1)**2)/vp(nno(l1))
          enddo
          irlast=1
          do ir=1,nr
            if(r(ir).lt.rlast)irlast=ir+1
          enddo
          do ir=irlast,nr
            t1=tlast+shead*(r(ir)-rlast)
            twinmin=dmin1(twinmin,t1-t0(ir))
          enddo
          smin=shead
          hpmin=hp(l)
        endif
c
        if(twinmin.lt.twindow)lcut=l
      enddo
      if(lcut.lt.1)then
        stop ' time window too small!'
      else if(lcut.lt.lp)then
        lp=lcut
        hp(lp)=0.d0
        n0=nno(lp)
        write(*,'(a,i3)')' actually used number of layers: ',n0
      endif
c
c     for partial solution only
c
      do i=1,lp
        n=nno(i)
        pup(i)=.true.
        pdw(i)=.true.
        if(vs(n).gt.vspmin*vp(n))then
          svup(i)=.true.
          svdw(i)=.true.
          sh(i)=.true.
        else
          svup(i)=.false.
          svdw(i)=.false.
          sh(i)=.false.
        endif
      enddo
      if(ipartial.eq.1)then
        z=zr
        do i=1,lp-1
          z=z+0.5d0*hp(i)
          do j=1,npar
            if(z.ge.zup(j).and.z.le.zlow(j))then
              if(ipsv(j).eq.1)pup(i)=.false.
              if(ipsv(j).eq.2)pdw(i)=.false.
              if(ipsv(j).eq.3)svup(i)=.false.
              if(ipsv(j).eq.4)svdw(i)=.false.
            endif
          enddo
          z=z+0.5d0*hp(i)
        enddo
      endif
      if(ipath.eq.1)then
        z=dmax1(zs,zr)
        lpath=max0(ls,lzr)
        do i=max0(ls,lzr)+1,lp
          z=z+hp(i-1)
          if(pathdepth.ge.z)lpath=i
        enddo
        if(lpath.eq.lp)then
          print *,'the depth limit for path filter is too large!'
          print *,'=> no signals in the given time window!'
          stop
        endif
      else
        lpath=0
      endif
c
      write(*,'(a)')' The receiver distance profile:'
      do j=1,nr
        if(mod(j,8).eq.0)then
          write(*,'(f10.3)')r(j)/km2m
        else
          write(*,'(f10.3,$)')r(j)/km2m
        endif
      enddo
      write(*,'(a)')' km'
c
      do istp=1,7
        do i=1,4
          if(ssel(istp).ge.1)then
            fsel(i,istp)=1
          else
            fsel(i,istp)=0
          endif
        enddo
      enddo
c
c     no toroidal component if ms = 0
c
      do istp=1,6
        if(ms(istp).eq.0)fsel(3,istp)=0
      enddo
c
      calsh=ssel(2).eq.1.or.ssel(3).eq.1.or.ssel(6).eq.1
c
c     for marine seismic
c
      if(lzr.lt.lp)then
        nup=nno(lzr)
        nlw=nno(lzr+1)
        if(vs(nup).lt.vspmin*vp(nup).and.
     +     vs(nlw).lt.vspmin*vp(nlw))then
          do istp=1,7
            do i=1,3
              fsel(i,istp)=0
            enddo
          enddo
        endif
      endif
1000  format(i5,f12.2,3f11.4,2f8.1)
c
      return
      end
