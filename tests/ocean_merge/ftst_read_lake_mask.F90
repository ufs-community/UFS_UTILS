! Unit test for the read_lake_mask routine.
!
! Reads a 6x4 version of the lake mask NetCDF file and
! checks values from the lake fraction, lake depth
! and latitude records. If differences exceed a
! threshold, then the test fails.
!
 program read_lake_info

 implicit none

 integer, parameter     :: lon = 6
 integer, parameter     :: lat = 4
 integer, parameter     :: tile = 1

 character(len=3)       :: pth2="./"
 character(len=3)       :: atmres="C48"

 real, parameter        :: thresh = 0.001
 real                   :: lake_frac(lon,lat)
 real                   :: lake_depth(lon,lat)
 real                   :: lat2d(lon,lat)

 print*,"Call routine read_lake_mask."

 call read_lake_mask(pth2,atmres,tile,lon,lat,lake_frac, &
                           lake_depth,lat2d)

 print*,"Check records."

 if (abs(lake_depth(1,1)-0.0) > thresh) stop 2
 if (abs(lake_depth(5,2)-68.3817) > thresh) stop 4
 if (abs(lake_frac(4,3)-1.0) > thresh) stop 6
 if (abs(lake_frac(4,4)-0.0) > thresh) stop 8
 if (abs(lat2d(1,4)-0.399218) > thresh) stop 10
 if (abs(lat2d(6,1)-(-1.87402)) > thresh) stop 12

 print*,"OK"

 print*,"SUCCESS"

 end program read_lake_info
