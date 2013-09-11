      double precision function taper(k,kcut1,kcut2,kcut3,kcut4)
      implicit none
c
      double precision k,kcut1,kcut2,kcut3,kcut4
      double precision pid2
c
      pid2=2.d0*datan(1.d0)
      if(k.gt.kcut1.and.k.lt.kcut2)then
c
c       within the lower end taper
c
        taper=dcos(pid2*(kcut2-k)/(kcut2-kcut1))
      else if(k.gt.kcut3.and.k.lt.kcut4)then
c
c       within the higher end taper
c
        taper=dcos(pid2*(k-kcut3)/(kcut4-kcut3))
      else if(k.ge.kcut2.and.k.le.kcut3)then
        taper=1.d0
      else
        taper=0.d0
      endif
c
      return
      end