 program quick_sort

! Test the quicksort function of the chgres_cube program.
! The routine sorts numbers in ascending order. Sort
! an array of pressure levels in millibars.
!
! @author George Gayno NOAA/EMC

 use input_data

 implicit none

 real*8 :: pressure_in_mb(5), sorted_pressure_in_mb(5)
 integer :: first, last

 data pressure_in_mb /300.0, 600.0, 200.0, 50.0, 1000.0/
 data sorted_pressure_in_mb /50.0, 200.0, 300.0, 600.0, 1000.0/

 print*,'Starting test of quicksort function.'

 first = 1
 last = 5
 call quicksort(pressure_in_mb,first,last)

 if (pressure_in_mb(1) /= sorted_pressure_in_mb(1)) stop 2
 if (pressure_in_mb(2) /= sorted_pressure_in_mb(2)) stop 3
 if (pressure_in_mb(3) /= sorted_pressure_in_mb(3)) stop 4
 if (pressure_in_mb(4) /= sorted_pressure_in_mb(4)) stop 5
 if (pressure_in_mb(5) /= sorted_pressure_in_mb(5)) stop 6

 print*,"OK"

 print*,"SUCCESS!"

 end program quick_sort
