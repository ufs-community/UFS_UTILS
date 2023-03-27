program msistest

! Test the MSIS library. This test is based on the test that
! is part of the library. It was modified slightly so it could
! be run under Github actions.
!
! It reads a text file (that is also part of the library) of 200 test
! cases. The file contains variables input to routine gtd8d and
! the check values of the variables output from that routine.
! If the computed values differ from the check values by a small
! threshold, the test fails.
!
! @author George.Gayno@noaa.gov

!#######################################################################
! MSIS (NRL-SOF-014-1) SOFTWARE
! NRLMSIS empirical atmospheric model software. Use is governed by the
! Open Source Academic Research License Agreement contained in the file
! nrlmsis2.1_license.txt, which is part of this software package. BY
! USING OR MODIFYING THIS SOFTWARE, YOU ARE AGREEING TO THE TERMS AND
! CONDITIONS OF THE LICENSE.  
!#######################################################################

!!! ===========================================================================
!!! NRLMSIS 2.1:
!!! Neutral atmosphere empirical model from the surface to lower exosphere
!!! ===========================================================================

!==================================================================================================
! MSISTEST: Test program for NRLMSIS 2.1
!==================================================================================================

  use msis_init, only          : msisinit

  implicit none

  integer                     :: i, iyd, isec
  integer                     :: mass ! Not use by v2.

  real(4)                     :: sec, alt, glat, glong, stl, f107a, f107, ap(7), apd
  real(4)                     :: d(10), t(2), d_check(10), t_check

  print*,'Starting test of msis library.'

  !Initialize model
  call msisinit(parmpath='./data/',parmfile='msis21.parm')

! Open file that contains the variables input to routine
! gtd8d and the check values for the variables output from
! that routine.

  open(78,file='./data/msis2.1_test_ref_dp.txt',status='old')
  read(78,*)  ! Ignore first line of file.

! Loop through all 200 records.

  do i = 1,200

    read(78,'(2i7,3f7.1,f7.2,3f7.1,10e13.4,1e13.4,f8.2)')  &
         iyd,isec,alt,glat,glong,stl,f107a,f107,apd,d_check,t_check

    ap(1) = apd
    sec = float(isec)

! Input variables are:
!
!   iyd   - year and day
!   sec   - seconds
!   alt   - altitude in km
!   glat  - latitude in deg
!   glong - longitude in deg
!   stl   - local solar time
!   f107a - 81 day average of solar activity index
!   f107  - daily solar activity indiex
!   ap    - daily geomagnetic activity index
 
    call gtd8d(iyd,sec,alt,glat,glong,stl,f107a,f107,ap,mass,d,t)

    print*,'Check case ',i

! Check He number density
     call checkit(d(1), d_check(1), 'He')

! Check O number density.
     call checkit(d(2), d_check(2), 'O')

! Check N2 number density.
     call checkit(d(3), d_check(3), 'N2')

! Check O2 number density.
     call checkit(d(4), d_check(4), 'O2')

! Check Ar number density.
     call checkit(d(5), d_check(5), 'Ar')

! Check Total mass density.
     call checkit(d(6), d_check(6), 'TMD')

! Check H number density.
     call checkit(d(7), d_check(7), 'H')

! Check N number density.
     call checkit(d(8), d_check(8), 'N')

! Check Anomalous oxygen number density.
     call checkit(d(9), d_check(9), 'AO')

! Check NO number density.
     call checkit(d(10), d_check(10), 'NO')

! Check temperature at altitude.
     if ( abs(t(2)-t_check) > 0.01) stop 28

  enddo

  close(78)

  print*,"OK"
  print*,"SUCCESS!"

end program msistest

! Routine to check the calculated values against the
! check values. A percentage threshold is used. Some
! fields contain missing values (flag value 9.99e-38). 

subroutine checkit(calc_val, check_val, field)

 implicit none

 character(len=*), intent(in) :: field

 real(4), intent(in) :: calc_val, check_val

 real                :: epsilon

 if (check_val > 1) then  ! Check value is not missing.
   epsilon =  abs(calc_val-check_val) / check_val
   if (epsilon > 0.001) then
     print*,'Bad value of ', field
     stop 8
   endif
 else ! Check value is missing. Is computed value also missing?
   if (calc_val > 1) then
     print*,'Value not missing for field ', field
     stop 9
   endif
 endif

end subroutine checkit
