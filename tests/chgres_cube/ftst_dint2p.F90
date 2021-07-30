 program ftst_dint2p

 use input_data

 implicit none

 integer, parameter :: npin = 3
 integer, parameter :: npout = 3

 integer :: ier, linlog

 real*8 :: ppin(npin), xxin(npin)
 real*8 :: ppout(npout), xmsg
 real*8 :: xxout(npout)

 data ppin /1000., 700., 300./
 data xxin /1., 2., 3. /
 data ppout /850., 600., 400./

 xmsg = -999.0

 linlog = 1
 
 call dint2p(ppin,xxin,npin,ppout,xxout,npout, &
             linlog, xmsg, ier)

 if (ier /= 0) stop 2
 if (xxout(1) /= 1.5) stop 4
 if (xxout(2) /= 2.25) stop 6
 if (xxout(3) /= 2.75) stop 8

 print*, "OK"
 print*, "SUCCESS!"

 end program ftst_dint2p
