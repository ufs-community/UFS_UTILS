!#######################################################################
! MSIS® (NRL-SOF-014-1) SOFTWARE
! NRLMSIS® empirical atmospheric model software. Use is governed by the
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
program msistest

  use msis_init, only          : msisinit

  implicit none

  integer, parameter          :: nrec = 200

  integer                     :: iyd, mass, iyd_check, sec_check
  real(4)                     :: sec, alt, glat, glong, stl, f107a, f107, ap(7), apd
  real(4)                     :: alt_check, glat_check, glong_check, &
                                 stl_check, f107a_check, f107_check, ap_check
  real(4)                     :: d(10),t(2)
  real(4)                     :: d_check(10),t_check, epsilon
  
  integer                     :: i
  character(128)              :: dummy

  !Initialize model
  call msisinit(parmpath='./data/',parmfile='msis21.parm')

  !Open input and output files, loop through records, and call model
  open(77,file='./data/msis2.1_test_in.txt',status='old')
  open(78,file='/scratch1/NCEPDEV/da/George.Gayno/ufs_utils.git/UFS_UTILS/sorc/chgres_cube.fd/msis2.1.fd/msis2.1_test_ref_dp.txt',status='old')
  read(77,*) dummy
  read(78,*)
  do i = 1,200
    read(77,*) iyd,sec,alt,glat,glong,stl,f107a,f107,apd
    ap(1) = apd
    call gtd8d(iyd,sec,alt,glat,glong,stl,f107a,f107,ap,mass,d,t)
    print*,'after gdt8d ',i,d,t
    read(78,'(2i7,3f7.1,f7.2,3f7.1,10e13.4,1e13.4,f8.2)')  &
      iyd_check,sec_check,alt_check,glat_check,glong_check,stl_check, &
      f107a_check,f107_check,ap_check,d_check(1:10),t_check
      print*,'after read ',i,iyd_check,sec_check,alt_check,glat_check,glong_check, &
               stl_check,f107a_check,f107_check,ap_check,d_check,t_check

! Check He number density
     call checkit(d(1), d_check(1), 2)
! Check O number density.
     call checkit(d(2), d_check(2), 4)
! Check N2 number density.
     call checkit(d(3), d_check(3), 6)
! Check O2 number density.
     call checkit(d(4), d_check(4), 8)
! Check Ar number density.
     call checkit(d(5), d_check(5), 10)
! Check Total mass density.
     call checkit(d(6), d_check(6), 12)
! Check H number density.
     call checkit(d(7), d_check(7), 14)
! Check N number density.
     call checkit(d(8), d_check(8), 16)
! Check Anomalous oxygen number density.
     call checkit(d(9), d_check(9), 18)
! Check NO number density.
     call checkit(d(10), d_check(10), 20)
! Check temperature at altitude.
     if ( abs(t(2)-t_check) > 0.01) stop 28

  enddo

  close(77)
  close(78)

  print*,"SUCCESS!"

end program msistest

subroutine checkit(calc_val, check_val, istat)

 implicit none

 integer, intent(in) :: istat
 real(4), intent(in) :: calc_val, check_val

 real :: epsilon

 if (check_val > 1) then
   epsilon =  abs(calc_val-check_val) / check_val
   if (epsilon > 0.001) stop istat
 else ! Check value is missing. Is computed value also missing?
   if (calc_val > 1) stop istat
 endif

end subroutine checkit
