      subroutine qswvint(srate)
      implicit none
c
      double precision srate
c
      include 'qsglobal.h'
c
      integer istp,n,nrec,l,lf,lf1,i,ir,nrr,nrs,ns
      integer ik,ik1,ik2,nk1,nk2,nbsj,idtrans
      double precision f,fcut,r0,k,kmax,kc,dk,slwn
      double precision pi,pi2,rr,rs,delta,cmax,ymax,yabs
      double precision fac,zdis,rsdis,thickness
      double precision kcut(4),kcut1(4),kcut2(4)
      double complex carec,cbrec,ck,ck2,cdk,c2dk,cdk2,cfac,czdis2
      double complex swap(nbsjmax+ndtransmax),y(4,6,nbsjmax+ndtransmax)
      double complex cy(6,6),yb(4),cics(6),cm2(4,6)
      double precision taper
      character*80 kernels
c
      double precision eps
      double complex c2
      data eps/1.0d-03/
      data c2/(2.d0,0.d0)/
c
c     ics = 1  when the azmuth-factor is cos(ms*theta) for poloidal mode
c             (psv) and sin(ms*theta) for the toroidal mode (sh);
c     ics = -1 otherwise.
c
      pi=4.d0*datan(1.d0)
      pi2=2.d0*pi
c
      do istp=1,6
        cics(istp)=dcmplx(dble(ics(istp)),0.d0)
        cm2(1,istp)=dcmplx(dble(ms(istp)**2),0.d0)
        cm2(2,istp)=dcmplx(dble((ms(istp)-1)**2),0.d0)
        cm2(3,istp)=dcmplx(dble((ms(istp)+1)**2),0.d0)
        cm2(4,istp)=cm2(1,istp)
      enddo
c
      nrec=nno(lzr)
      carec=dcmplx(1.d0/(ro(nrec)*vp(nrec)**2),0.d0)
      cbrec=dcmplx(2.d0*(vs(nrec)/vp(nrec))**2,0.d0)
c
      ns=nno(ls)
      l=min0(ls,lzr)
      n=nno(l)
      if(vs(n).le.vspmin*vp(n).or.
     &   .not.(svup(l).or.svdw(l)))then
        cmax=vp(n)
      else
        cmax=vs(n)
      endif
      do l=min0(ls,lzr)+1,max0(ls,lzr,lpath)-1
        n=nno(l)
        if(vs(n).le.vspmin*vp(n).or.
     &     .not.(svup(l).or.svdw(l)))then
          cmax=dmax1(cmax,vp(n))
        else
          cmax=dmax1(cmax,vs(n))
        endif
      enddo
c
      rsdis=0.d0
      do l=1,max0(ls,lzr,lpath)-1
        rsdis=rsdis+hp(l)
      enddo
c
      fcut=dble(nf)*df
c
      zdis=dabs(zr-zs)
      czdis2=dcmplx(zdis*zdis,0.d0)
c
      thickness=0.d0
      do l=1,lp-1
        thickness=thickness+hp(l)
      enddo
c
      dk=pi/dmax1(srate*r(nr),3.d0*thickness)
      cdk=dcmplx(dk,0.d0)
      c2dk=dcmplx(2.d0*dk,0.d0)
      cdk2=dcmplx(dk*dk,0.d0)
c
      r0=dmax1(cmax/(0.1d0*fcut),1.d0/(200.d0*dk))
c
      if(fullwave)then
        f=0.d0
        call qsqmodel(f)
        call qssource(ro(ns),cvp(ns),cvs(ns))
        kc=0.d0
        ik=1
        ymax=0.d0
101     k=dble(ik)*dk
        call qskern(cy,f,k)
        yabs=0.d0
        do istp=1,4
          do i=1,3
            yabs=yabs+cdabs(cy(2*i-1,istp))**2
          enddo
        enddo
        yabs=dsqrt(yabs)*k
        if(k*r(1).gt.1.d0)yabs=yabs/dsqrt(k*r(1))
        yabs=yabs*dexp(-0.25d0*(k*r0)**2)
        ymax=dmax1(ymax,yabs)
        if(yabs.gt.eps*ymax.and.yabs.gt.0.d0)then
         ik=ik+1
          goto 101
        endif
        kcut1(1)=0.d0
        kcut1(2)=0.d0
        kcut1(3)=0.d0
        kcut1(4)=k
c
        f=fcut
        call qsqmodel(f)
        call qssource(ro(ns),cvp(ns),cvs(ns))
        kc=1.15d0*pi2*f/cmax
        ik=1
        ymax=0.d0
102     k=dble(ik)*dk
        call qskern(cy,f,k)
        yabs=0.d0
        do istp=1,4
          do i=1,3
            yabs=yabs+cdabs(cy(2*i-1,istp))**2
          enddo
        enddo
        yabs=dsqrt(yabs)*k
        if(k*r(1).gt.1.d0)yabs=yabs/dsqrt(k*r(1))
        if(k.gt.kc)yabs=yabs*dexp(-0.25d0*((k-kc)*r0)**2)
        ymax=dmax1(ymax,yabs)
        if(yabs.gt.eps*ymax.and.yabs.gt.0.d0)then
          ik=ik+1
          goto 102
        endif
c
        kcut2(1)=0.d0
        kcut2(2)=0.d0
        kcut2(3)=dmin1(k,1.15d0*pi2*fcut/cmax)
        kcut2(4)=k
        lf1=1
      else
        kcut1(1)=0.d0
        kcut1(2)=0.d0
        kcut1(3)=0.d0
        kcut1(4)=0.d0
        do i=1,4
          kcut2(i)=pi2*fcut*slw(i)
        enddo
        lf1=2
      endif
c
      nbsj=2+idint(dmax1(kcut1(4),kcut2(4))/dk)
c
c     for tests
c
      i=0
      if(i.eq.1)then
        write(*,'(a,$)')' frequency (hz) for kernel functions: '
        read(*,*)f
        write(*,'(a,$)')' file name for kernel functions: '
        read(*,'(a)')kernels
        call qsqmodel(f)
        call qssource(ro(ns),cvp(ns),cvs(ns))
        open(21,file=kernels,status='unknown')
        write(21,'(a)')'    slowness   '
     &   //'    ex1r        ex1i        ex2r        ex2i    '
     &   //'    ex3r        ex3i        ex4r        ex4i    '
     &   //'    ss1r        ss1i        ss2r        ss2i    '
     &   //'    ss3r        ss3i        ss4r        ss4i    '
     &   //'    ds1r        ds1i        ds2r        ds2i    '
     &   //'    ds3r        ds3i        ds4r        ds4i    '
     &   //'    cl1r        cl1i        cl2r        cl2i    '
     &   //'    cl3r        cl3i        cl4r        cl4i    '
     &   //'    fz1r        fz1i        fz2r        fz2i    '
     &   //'    fz3r        fz3i        fz4r        fz4i    '
     &   //'    fh1r        fh1i        fh2r        fh2i    '
     &   //'    fh3r        fh3i        fh4r        fh4i    '
        do i=1,4
          kcut(i)=kcut1(i)+(kcut2(i)-kcut1(i))*f/fcut
        enddo
        nk2=2+idint(kcut(4)/dk)
        ik1=1
        ik2=nk2+ndtrans
        do ik=ik1,ik2
          k=dble(ik)*dk
          ck=dcmplx(k,0.d0)
          call qskern(cy,f,k)
          do istp=1,6
            y(1,istp,ik)=cy(1,istp)
            y(2,istp,ik)=( cy(3,istp)+cics(istp)*cy(5,istp))/c2
            y(3,istp,ik)=(-cy(3,istp)+cics(istp)*cy(5,istp))/c2
            y(4,istp,ik)=carec*cy(2,istp)-cbrec*ck*cy(3,istp)
          enddo
        enddo
c
        do idtrans=1,ndtrans
          ik1=1+idtrans
          ik2=nk2+ndtrans-idtrans
          do istp=1,6
            do i=1,4
              do ik=ik1-1,ik2+1
                swap(ik)=y(i,istp,ik)
              enddo
              do ik=ik1,ik2
                ck=dcmplx(dble(ik)*dk,0.d0)
                ck2=ck*ck
                y(i,istp,ik)=swap(ik)*(czdis2+cm2(i,istp)/ck2)
     &            -(swap(ik+1)-swap(ik-1))/c2dk/ck
     &            -(swap(ik+1)-c2*swap(ik)+swap(ik-1))/cdk2
              enddo
            enddo
          enddo
        enddo
        do i=1,4
          kcut(i)=kcut1(i)+(kcut2(i)-kcut1(i))*f/fcut
        enddo
        ik1=1+ndtrans
        ik2=nk2
        do ik=ik1,ik2
          k=dble(ik)*dk
          cfac=dcmplx(taper(k,kcut(1),kcut(2),kcut(3),kcut(4)),0.d0)
          do istp=1,6
            do i=1,4
              y(i,istp,ik)=y(i,istp,ik)*cfac
            enddo
          enddo
          slwn=1000.d0*dble(ik)*dk/(pi2*f)
          write(21,'(E15.7,48E12.4)')slwn,
     &                 ((y(i,istp,ik),i=1,4),istp=1,6)
        enddo
        close(21)
        stop
      endif
c
c     end tests
c
      if(nbsj.gt.nbsjmax)then
        stop ' parameter nbsjmax defined too small'
      else
        print *,' Calculate Bessel functions for x up to ',
     &       dble(nbsj)*dk*r(nr)
        do ir=1,nr
          geospr(ir)=1.d0/(zdis*zdis+r(ir)*r(ir))**ndtrans
        enddo
        call qsbsj(dk,nbsj)
      endif
c
      do lf=1,nf
        do istp=1,6
          do i=1,4
            do ir=1,nr
              grns(lf,i,ir,istp)=(0.d0,0.d0)
             enddo
          enddo
        enddo
      enddo
c
      write(*,'(a,2(f10.7,a))')' Min./max. slowness at f_cut: ',
     &     1000.d0*kcut2(1)/(pi2*fcut),' / ',
     &     1000.d0*kcut2(4)/(pi2*fcut),' s/km'
c
      do lf=lf1,nf
        f=dble(lf-1)*df
        call qsqmodel(f)
        call qssource(ro(ns),cvp(ns),cvs(ns))
c
        do i=1,4
c          kcut(i)=kcut1(i)+dble(lf-1)*(kcut2(i)-kcut1(i))/dble(nf)
          kcut(i)=kcut1(i)
     &           +(kcut2(i)-kcut1(i))*dsqrt(f**2+(pi*fi)**2)/fcut
        enddo
c
        nk2=min0(2+idint(kcut(4)/dk),nbsj)
        nk1=min0(1+idint(kcut(1)/dk),nk2)
c
        ik1=max0(1,nk1-ndtrans)
        ik2=nk2+ndtrans
        do ik=ik1,ik2
          k=dble(ik)*dk
          ck=dcmplx(k,0.d0)
          call qskern(cy,f,k)
          do istp=1,6
            y(1,istp,ik)=cy(1,istp)
            y(2,istp,ik)=( cy(3,istp)+cics(istp)*cy(5,istp))/c2
            y(3,istp,ik)=(-cy(3,istp)+cics(istp)*cy(5,istp))/c2
            y(4,istp,ik)=carec*cy(2,istp)-cbrec*ck*cy(3,istp)
          enddo
        enddo
        do idtrans=1,ndtrans
          ik1=max0(1,nk1-ndtrans)+idtrans
          ik2=nk2+ndtrans-idtrans
          do istp=1,6
            do i=1,4
              do ik=ik1-1,ik2+1
                swap(ik)=y(i,istp,ik)
              enddo
              do ik=ik1,ik2
                ck=dcmplx(dble(ik)*dk,0.d0)
                ck2=ck*ck
                y(i,istp,ik)=swap(ik)*(czdis2+cm2(i,istp)/ck2)
     &            -(swap(ik+1)-swap(ik-1))/c2dk/ck
     &            -(swap(ik+1)-c2*swap(ik)+swap(ik-1))/cdk2
              enddo
            enddo
          enddo
        enddo
c
        ik1=max0(1,nk1+ndtrans)
        ik2=nk2
        do ik=ik1,ik2
          k=dble(ik)*dk
          fac=k*dk
          if(.not.fullwave)then
            fac=fac*taper(k,kcut(1),kcut(2),kcut(3),kcut(4))
          else if(k.gt.kcut(3))then
            fac=fac*dexp(-0.25d0*((k-kcut(3))*r0)**2)
          endif
          cfac=dcmplx(fac,0.d0)
c
          do istp=1,6
            do i=1,4
              y(i,istp,ik)=y(i,istp,ik)*cfac
            enddo
            do ir=1,nr
              yb(1)=y(1,istp,ik)*dcmplx(bsj(ik,ms(istp),ir),0.d0)
              yb(2)=y(2,istp,ik)*dcmplx(bsj(ik,ms(istp)-1,ir),0.d0)
              yb(3)=y(3,istp,ik)*dcmplx(bsj(ik,ms(istp)+1,ir),0.d0)
              yb(4)=y(4,istp,ik)*dcmplx(bsj(ik,ms(istp),ir),0.d0)
c
              grns(lf,1,ir,istp)=grns(lf,1,ir,istp)+yb(1)
              grns(lf,2,ir,istp)=grns(lf,2,ir,istp)+yb(2)+yb(3)
              grns(lf,3,ir,istp)=grns(lf,3,ir,istp)
     &                          -cics(istp)*(yb(2)-yb(3))
              grns(lf,4,ir,istp)=grns(lf,4,ir,istp)+yb(4)
            enddo
          enddo
        enddo
c
        write(*,'(i6,a,E13.6,a,i7)')lf,'.',f,
     &      'Hz: slowness samples = ',1+nk2-nk1-ndtrans
      enddo
c
      if(iflat.eq.1)then
c
c       amplitude correction when using the flat-earth transform
c       see Mueller (1977) for n = -2
c
        rs=rr0*dexp(-zs/rr0)
        rr=rr0*dexp(-zrrs/rr0)
        nrs=5-ndens
        nrr=3-ndens
        do ir=1,nr
          if(r(ir).gt.0.d0)then
            delta=r(ir)/rr0
            fac=delta/dsin(delta)
          else
            fac=1.d0
          endif
          cfac=dcmplx(dsqrt((rr0/rr)**nrr*(rr0/rs)**nrs*fac),0.d0)
          do istp=1,6
            do i=1,4
              do lf=lf1,nf
                grns(lf,i,ir,istp)=grns(lf,i,ir,istp)*cfac
              enddo
            enddo
          enddo
        enddo
      endif
c
      return
      end