      program qseis
      implicit none
c
      include 'qsglobal.h'
c
c     work space
c
      integer i,istp,nssel,runtime
      double precision pi,srate
      logical grnexist
      integer time
c
c     read input file file
c
      print *,'######################################################'
      print *,'#                                                    #'
      print *,'#               Welcome to the program               #'
      print *,'#                                                    #'
      print *,'#                                                    #'
      print *,'#        QQQ     SSSS    EEEEE    III     SSSS       #'
      print *,'#       Q   Q   S        E         I     S           #'
      print *,'#       Q Q Q    SSS     EEEE      I      SSS        #'
      print *,'#       Q  QQ       S    E         I         S       #'
      print *,'#        QQQQ   SSSS     EEEEE    III    SSSS        #'
      print *,'#                                                    #'
      print *,'#                  (Version 2006)                    #'
      print *,'#                                                    #'
      print *,'#                                                    #'
      print *,'#                      by                            #'
      print *,'#                 Rongjiang Wang                     #'
      print *,'#              (wang@gfz-potsdam.de)                 #'
      print *,'#                                                    #'
      print *,'#           GeoForschungsZentrum Potsdam             #'
      print *,'#           Last modified: October, 2006             #'
      print *,'######################################################'
      print *,'                          '
      write(*,'(a,$)')' the input data file is '
      read(*,'(a)')inputfile
      runtime=time()
c
      pi=4.d0*datan(1.d0)
c
      open(10,file=inputfile,status='old')
      call qsgetinp(10,srate,nssel)
      close(10)
c
      if(nssel-ssel(7).gt.0)then
        grnexist=.false.
        call qswvint(srate)
        call qsmultis(grnexist)
        iexist=0
        do istp=1,7
          do i=1,4
            if(fsel(i,istp).eq.1)then
              call qsfftinv(i,istp)
            endif
          enddo
        enddo
      else
        grnexist=.true.
        call qsmultis(grnexist)
      endif
c
      runtime=time()-runtime
      write(*,'(a)')' #############################################'
      write(*,'(a)')' #                                           #'
      write(*,'(a)')' #      End of computations with qseis06     #'
      write(*,'(a)')' #                                           #'
      write(*,'(a,i10,a)')' #       Run time: ',runtime,
     +                                           ' sec            #'
      write(*,'(a)')' #############################################'
1001  format(2i7,E12.4,a)
1002  format(i4,a,E12.4,a,$)
1003  format(E12.5,$)
1004  format(2E12.4,$)
1005  format(2E12.4)
 500  stop
      end
