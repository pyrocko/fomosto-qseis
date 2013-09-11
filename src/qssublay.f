      subroutine qssublay(resolut,fcut)
      implicit none
c
      double precision resolut(3),fcut
c
      include 'qsglobal.h'
c
c     work space
c
      integer i,i0,l,ivp,ivs,iro
      double precision dh,dro,dvp,dvs,dqp,dqs,z,dz,halfwvlen
c
      n0=0
      do l=1,l0-1
        dz=z2(l)-z1(l)
        dvp=2.d0*dabs(vp2(l)-vp1(l))/(vp2(l)+vp1(l))
        if(vs1(l).gt.vspmin*vp1(l))then
          dvs=2.d0*dabs(vs2(l)-vs1(l))/(vs2(l)+vs1(l))
        else
          dvs=0.d0
        endif
        dro=2.d0*dabs(ro2(l)-ro1(l))/(ro2(l)+ro1(l))
c
        if(dvp.le.0.d0)then
          ivp=1
        else
          halfwvlen=0.25d0*(vp2(l)+vp1(l))/fcut
          ivp=1+idint(dz/halfwvlen)
          if(resolut(1).gt.0.d0)ivp=min0(ivp,1+idint(dvp/resolut(1)))
        endif
        if(dvs.le.0.d0)then
          ivs=1
        else
          halfwvlen=0.25d0*(vs2(l)+vs1(l))/fcut
          ivs=1+idint(dz/halfwvlen)
          if(resolut(2).gt.0.d0)ivs=min0(ivs,1+idint(dvs/resolut(2)))
        endif
        iro=1
        if(dro.gt.0.d0.and.resolut(3).gt.0.d0)then
          iro=1+idint(dro/resolut(3))
        endif
c
        i0=max0(ivp,ivs,iro)
        dro=(ro2(l)-ro1(l))/dz
        dvp=(vp2(l)-vp1(l))/dz
        dvs=(vs2(l)-vs1(l))/dz
        dqp=(qp2(l)-qp1(l))/dz
        dqs=(qs2(l)-qs1(l))/dz
        dh=dz/dble(i0)
        do i=1,i0
          n0=n0+1
          if(n0.ge.lmax)then
            stop ' Max. number of layers (lmax) too small defined!'
          endif
          h(n0)=dh
          z=(dble(i)-0.5d0)*dh
          ro(n0)=ro1(l)+dro*z
          vp(n0)=vp1(l)+dvp*z
          vs(n0)=vs1(l)+dvs*z
          qp(n0)=qp1(l)+dqp*z
          qs(n0)=qs1(l)+dqs*z
        enddo
      enddo
c
c     last layer is halfspace
c
      n0=n0+1
      h(n0)=0.d0
      ro(n0)=ro1(l0)
      vp(n0)=vp1(l0)
      vs(n0)=vs1(l0)
      qp(n0)=qp1(l0)
      qs(n0)=qs1(l0)
c
c     for receiver site
c
      n0rs=0
      if(l0rs.gt.0)then
        do l=1,l0rs-1
          dz=z2rs(l)-z1rs(l)
          dvp=2.d0*dabs(vp2rs(l)-vp1rs(l))/(vp2rs(l)+vp1rs(l))
          if(vs1rs(l).gt.vspmin*vp1rs(l))then
            dvs=2.d0*dabs(vs2rs(l)-vs1rs(l))/(vs2rs(l)+vs1rs(l))
          else
            dvs=0.d0
          endif
          dro=2.d0*dabs(ro2rs(l)-ro1rs(l))/(ro2rs(l)+ro1rs(l))
c
          if(dvp.le.0.d0)then
            ivp=1
          else
            halfwvlen=0.25d0*(vp2rs(l)+vp1rs(l))/fcut
            ivp=1+idint(dz/halfwvlen)
            if(resolut(1).gt.0.d0)ivp=min0(ivp,1+idint(dvp/resolut(1)))
          endif
          if(dvs.le.0.d0)then
            ivs=1
          else
            halfwvlen=0.25d0*(vs2rs(l)+vs1rs(l))/fcut
            ivs=1+idint(dz/halfwvlen)
            if(resolut(2).gt.0.d0)ivs=min0(ivs,1+idint(dvs/resolut(2)))
          endif
          iro=1
          if(dro.gt.0.d0.and.resolut(3).gt.0.d0)then
            iro=1+idint(dro/resolut(3))
          endif
c
          i0=max0(ivp,ivs,iro)
          dro=(ro2rs(l)-ro1rs(l))/dz
          dvp=(vp2rs(l)-vp1rs(l))/dz
          dvs=(vs2rs(l)-vs1rs(l))/dz
          dqp=(qp2rs(l)-qp1rs(l))/dz
          dqs=(qs2rs(l)-qs1rs(l))/dz
          dh=dz/dble(i0)
          do i=1,i0
            n0rs=n0rs+1
            if(n0rs.ge.lmax)then
              stop ' Max. number of layers (lmax) too small defined!'
            endif
            hrs(n0rs)=dh
            z=(dble(i)-0.5d0)*dh
            rors(n0rs)=ro1rs(l)+dro*z
            vprs(n0rs)=vp1rs(l)+dvp*z
            vsrs(n0rs)=vs1rs(l)+dvs*z
            qprs(n0rs)=qp1rs(l)+dqp*z
            qsrs(n0rs)=qs1rs(l)+dqs*z
          enddo
        enddo
c
c       last layer is halfspace
c
        n0rs=n0rs+1
        hrs(n0rs)=0.d0
        rors(n0rs)=ro1rs(l0rs)
        vprs(n0rs)=vp1rs(l0rs)
        vsrs(n0rs)=vs1rs(l0rs)
        qprs(n0rs)=qp1rs(l0rs)
        qsrs(n0rs)=qs1rs(l0rs)
      endif
c
      return
      end
