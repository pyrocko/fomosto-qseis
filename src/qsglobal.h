c     GLOBAL INDEX PARAMETERS FOR DEFINING ARRAYS
c     ===========================================
c     nzmax: max. interface index;
c     lmax: max. number of total homogeneous layers (lmax <= nzmax-2);
c     nrmax: max. number of traces;
c     nfmax: max. number of frequency samples
c
      integer nzmax,lmax,nrmax,nfmax,ndtransmax
      parameter(lmax=1500)
      parameter(nzmax=lmax+2)
      parameter(nrmax=101,nfmax=2048)
      parameter(ndtransmax=4)
c
c     INDEX PARAMETERS FOR BESSEL FUNCTION TABLES
c     ===========================================
c
      integer nbsjmax
      parameter(nbsjmax=30000)
c
c     INDEX PARAMETERS FOR SEISMOMETER CHARACTERISTICS
c     ================================================
c     (max. number of roots and poles)
c
      integer nrootmax,npolemax
      parameter(nrootmax=10,npolemax=10)
c
c     EARTH RADIUS IN METER
c     =====================
c
      double precision rr0,km2m
      parameter(rr0=6.371d+06,km2m=1.0d+03)
c
c     ATMOSPHERIC PARAMETERS
c     ======================
c     double precision roair,vpair
c     parameter(roair=0.1300d+01,vpair=0.3318d+03)
      double precision roair,vpair
      parameter(roair=0.d0,vpair=0.d0)
c
c     THE MININUM VS/VP RATIO: VSPMIN
c     ===============================
c     (if vs/vp < vspmin, then fluid medium is assumed)
c
      double precision vspmin
      parameter(vspmin=0.05d0)
c
c     FOR FLAT-EARTH-TRANSFORMATION
c     =============================
c
      integer ndens
      parameter(ndens=1)
c
      double complex accair,cvpair,kpair,comega
      common /airpara/ accair,cvpair,kpair,comega
c
c     zr: receiver depth
c     lzr: sublayer no of receiver
c
      integer lzr,lzrrs
      double precision zr,zrrs
      common /receiver/ zr,zrrs,lzr,lzrrs
c
      integer nr
      double precision r(nrmax)
      common /distance/ r,nr
c
      integer lp,nno(nzmax)
      double precision hp(nzmax)
      common /sublayer/ hp,lp,nno
c
      integer lprs,nnors(nzmax)
      double precision hprs(nzmax)
      common /rssublayer/ hprs,lprs,nnors
c
c     original model parameters
c
      integer l0
      double precision z1(lmax),z2(lmax),ro1(lmax),ro2(lmax)
      double precision vp1(lmax),vp2(lmax),vs1(lmax),vs2(lmax)
      double precision qp1(lmax),qp2(lmax),qs1(lmax),qs2(lmax)
      common /model0/z1,z2,ro1,ro2,vp1,vp2,vs1,vs2,qp1,qp2,qs1,qs2,l0
c
c     original model parameters (receiver site)
c
      integer l0rs
      double precision z1rs(lmax),z2rs(lmax),ro1rs(lmax),ro2rs(lmax)
      double precision vp1rs(lmax),vp2rs(lmax),vs1rs(lmax),vs2rs(lmax)
      double precision qp1rs(lmax),qp2rs(lmax),qs1rs(lmax),qs2rs(lmax)
      common /rsmodel0/z1rs,z2rs,ro1rs,ro2rs,vp1rs,vp2rs,vs1rs,vs2rs,
     +                 qp1rs,qp2rs,qs1rs,qs2rs,l0rs
c       
c     layered model parameter:
c     n0: number of homogeneous layers
c
      integer n0
      double precision h(lmax),ro(lmax),vp(lmax),vs(lmax)
      double precision qp(lmax),qs(lmax)
      common /model/ h,ro,vp,vs,qp,qs,n0
c
      double complex acc(lmax),kp(lmax),ks(lmax),cla(lmax),cmu(lmax)
      double complex cvp(lmax),cvs(lmax),wa(lmax),wb(lmax)
      common /cpara/ acc,kp,ks,cla,cmu,cvp,cvs,wa,wb
c       
c     layered model parameter:
c     n0rs: number of homogeneous layers (receiver site)
c
      integer n0rs
      double precision hrs(lmax),rors(lmax),vprs(lmax),vsrs(lmax)
      double precision qprs(lmax),qsrs(lmax)
      common /rsmodel/ hrs,rors,vprs,vsrs,qprs,qsrs,n0rs
c
      double complex accrs(lmax),kprs(lmax),ksrs(lmax),
     +               clars(lmax),cmurs(lmax)
      double complex cvprs(lmax),cvsrs(lmax),wars(lmax),wbrs(lmax)
      common /rscpara/ accrs,kprs,ksrs,clars,cmurs,cvprs,cvsrs,wars,wbrs
c
c     for partial solution
c
      integer ipartial,npar,ipsv(nzmax)
      double precision zup(nzmax),zlow(nzmax)
      common /partial/ zup,zlow,ipartial,npar,ipsv
      logical pup(nzmax),pdw(nzmax),svup(nzmax),svdw(nzmax),sh(nzmax)
      common /psvfilter/ pup,pdw,svup,svdw,sh
c
c     slowness cut-offs
c
      double precision slw(4)
      logical fullwave,calsh
      common /slwcutoffs/ slw,fullwave,calsh
c
c     source parameters
c
      integer ls,ms(6),ics(6)
      double precision zs
      double precision sfct0(6,6),sfct1(6,6)
      common /source/ zs,sfct0,sfct1,ls,ms,ics
c
c     path filtering
c
      integer iflat,ipath,lpath,isurf,ndtrans
      double precision pathdepth
      common /pathfilter/ pathdepth,iflat,ipath,lpath,isurf,ndtrans
c
      integer nt,nf
      double precision dt,df,fi
      common /sampling/ dt,df,fi,nt,nf
c
      double precision tstart,twindow,tau,v0
      common /tparas/ tstart,twindow,tau,v0
c
      double precision mtensor(6),azimuth(nrmax)
      common /eqparas/ mtensor,azimuth
c
      integer nnmax,nn0,iexist,wdeg
      parameter(nnmax=1024)
      double precision wv0(nnmax)
      common /wavelets/ wv0,nn0,iexist,wdeg
c
c     seismometer filtering
c
      integer nroot,npole
      double precision asm
      double complex root(nrootmax),pole(npolemax) 
      common /seismometer/ root,pole,asm,nroot,npole
c
c     table of J_n(x), n = -1, 0, 1, 2, 3
c
      double precision bsj(nbsjmax,-1:3,nrmax),geospr(nrmax)
      common /bessels/ bsj,geospr
c
c     green's functions
c
      double complex grns(nfmax,4,nrmax,7)
      common /grnfcts/ grns
c
c     title text
c
      character*1 comptxt(4),varbtxt
      character*4 rcvtxt(nrmax)
      common /title/ comptxt,varbtxt,rcvtxt
c
c     input and output data files
c
      character*80 inputfile
      common /inputdata/ inputfile
      integer ssel(7),fsel(4,7),flen(4,7)
      character*83 outfile(4,7)
      common /outsel/ ssel,fsel,flen
      common /outdata/ outfile
