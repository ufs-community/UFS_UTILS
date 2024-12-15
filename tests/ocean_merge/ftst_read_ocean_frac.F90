! Unit test for the read_ocean_frac routine.
!
! Reads a 6x5 version of the MOM6 ocean mask file and
! checks values from the ocean fraction record.
! If differences exceed a threshold, then the test fails.
!
 program read_ocean_info

 implicit none

 integer, parameter     :: lon = 6
 integer, parameter     :: lat = 5
 integer, parameter     :: tile = 1

 character(len=3)       :: pth1="./"
 character(len=3)       :: atmres="C48"
 character(len=5)       :: ocnres="mx500"

 real, parameter        :: thresh = 0.001
 real                   :: ocn_frac(lon,lat)

 print*,"Call routine read_ocean_frac."

 call read_ocean_frac(pth1,atmres,ocnres,tile,lon,lat,ocn_frac)

 print*,"Check records."

 if (abs(ocn_frac(1,1)-1.0) > thresh) stop 2
 if (abs(ocn_frac(6,5)-0.0) > thresh) stop 4
 if (abs(ocn_frac(4,5)-0.0917) > thresh) stop 6
 if (abs(ocn_frac(3,1)-0.162347) > thresh) stop 8

 print*,"OK"

 print*,"SUCCESS"

 end program read_ocean_info
