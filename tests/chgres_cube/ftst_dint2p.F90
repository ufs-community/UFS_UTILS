 program ftst_dint2p

 use utilities 

 implicit none

 integer, parameter :: npin = 5
 integer, parameter :: npout1 = 4,npout2 = 6
 real, parameter    :: eps=1E-4

 integer :: ier, linlog

 real*8 :: ppin(npin), xxin(npin)
 real*8 :: ppout1(npout1), xmsg
 real*8 :: xxout1(npout1)
 real*8 :: ppout2(npout2)
 real*8 :: xxout2(npout2)

! input pressure level and variable
 data ppin /1000., 900., 700., 550., 300./
 data xxin /1., -999.0, 2., -999.0, 3. /
! output pressure levels and variable
 data ppout1 /850., 700., 600., 400./
 data ppout2 /850., 700., 600., 400., 200., 100./


 xmsg = -999.0

! linear interpolation; no extrapolation
 linlog = 1
 
 call dint2p(ppin,xxin,npin,ppout1,xxout1,npout1, &
             linlog, xmsg, ier)

 if (ier /= 0) stop 1 
 if (xxout1(1) /= 1.5) stop 11 
 if (xxout1(2) /= 2.0) stop 12
 if (xxout1(3) /= 2.25) stop 13
 if (xxout1(4) /= 2.75) stop 14

! linear interpolation; linear extrapolation
 linlog = -1

 call dint2p(ppin,xxin,npin,ppout2,xxout2,npout2, &
             linlog, xmsg, ier)

 if (ier /= 0) stop 2
 if (xxout2(1) /= 1.5) stop 21 
 if (xxout2(2) /= 2.0) stop 22
 if (xxout2(3) /= 2.25) stop 23
 if (xxout2(4) /= 2.75) stop 24
 if (xxout2(5) /= 3.25) stop 25
 if (xxout2(6) /= 3.5) stop 26



! lnP interpolation; no extrapolation
 linlog = 2

 call dint2p(ppin,xxin,npin,ppout1,xxout1,npout1, &
             linlog, xmsg, ier)

 if (ier /= 0) stop 3 
 if (abs(xxout1(1)-1.45565) .gt. eps) stop 31 
 if (abs(xxout1(2)-2.00000) .gt. eps) stop 32
 if (abs(xxout1(3)-2.18193) .gt. eps) stop 33
 if (abs(xxout1(4)-2.66047) .gt. eps) stop 34

! lnP interpolation; lnP extrapolation
 linlog = -2

 call dint2p(ppin,xxin,npin,ppout2,xxout2,npout2, &
             linlog, xmsg, ier)

 if (ier /= 0) stop 4 
 if (abs(xxout2(1)-1.45565) .gt. eps) stop 41
 if (abs(xxout2(2)-2.00000) .gt. eps) stop 42
 if (abs(xxout2(3)-2.18193) .gt. eps) stop 43
 if (abs(xxout2(4)-2.66047) .gt. eps) stop 44
 if (abs(xxout2(5)-3.47854) .gt. eps) stop 45
 if (abs(xxout2(6)-4.29661) .gt. eps) stop 46




 print*, "OK"
 print*, "SUCCESS!"

 end program ftst_dint2p
