program ftst_land_increments

! Test "apply_land_da_adjustments" using sample points.
!
! author: George Gayno (george.gayno@noaa.gov)
! author: Yuan Xue: add a total of four testing points to test out
! the newly added frozen soil ice fraction calculations
! under different scenarios after ingesting
! soil temp increments (08/2023)

 use mpi
 use land_increments

 implicit none

 integer :: my_rank, ierr, lsm, isot, ivegsrc
 integer :: lensfc, lsoil, l, lsoil_incr

 real, parameter :: EPSILON=0.001

 real(kind=4), allocatable :: zsoil(:)
 real, allocatable :: rsoiltype(:)
 integer, allocatable :: mask(:)
 real, allocatable :: stc_bck(:,:)
 real, allocatable :: smc_anl(:,:)
 real, allocatable :: slc_anl(:,:)
 real, allocatable :: stc_anl(:,:)
 integer, allocatable :: stc_updated(:), slc_updated(:)

 call mpi_init(ierr)

 call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)

 lsm=2             ! For noah-mp land surface model.
 isot = 1          ! STATSGO soil type.
 ivegsrc = 1       ! IGBP vegetation type.
 lsoil = 4         ! Number of soil layers.
 lsoil_incr = 3    ! Number of soil layers (from top) to apply increments to.
 
 lensfc= 4  ! Number of test points.
 
 allocate(zsoil(lsoil))
 allocate(rsoiltype(lensfc))      ! Soil type.
 allocate(mask(lensfc))           ! Land mask
 allocate(stc_bck(lensfc,lsoil))  ! Background soil temperature (K).
 allocate(smc_anl(lensfc,lsoil))  ! Analyzed total soil moisture.
 allocate(slc_anl(lensfc,lsoil))  ! Analyzed liquid soil moisture.
 allocate(stc_anl(lensfc,lsoil))  ! Analyzed soil temperature (K).
 allocate(stc_updated(lensfc))
 allocate(slc_updated(lensfc))

 zsoil(1) = -0.1
 zsoil(2) = -0.4
 zsoil(3) = -1.0
 zsoil(4) = -2.0
 slc_updated=1
 stc_updated=1

! Point 1 is above freezing before the adjustment
! and above freezing after the adjustment. Therefore,
! the increments to STC will be retained.
! temp: unfrozen ==> unfrozen

 rsoiltype(1) = 5.
 mask(1) = 1
 stc_bck(1,:) =  280.0

 smc_anl(1,:) = .25
 slc_anl(1,:) = .25
 stc_anl(1,1:lsoil_incr) =  281.0

! Point 2 is below freezing before the adjustment
! and above freezing after the adjustment. Therefore,
! all soil ice will be melted
! temp: frozen ==> unfrozen

 rsoiltype(2) = 5.
 mask(2) = 1
 stc_bck(2,:) =  271.0

 smc_anl(2,:) = .25
 slc_anl(2,:) = .20
 stc_anl(2,1:lsoil_incr) =  275.0

! Point 3 freezes before and after the adjustment. Therefore,
! SLC will be recomputed, SMC remained unchanged
! temp: frozen ==> frozen

 rsoiltype(3) = 5.
 mask(3) = 1
 stc_bck(3,:) =  272.0

 smc_anl(3,:) = .25
 slc_anl(3,:) = .20
 stc_anl(3,1:lsoil_incr) =  271.0

! Point 4 is above freezing before and below freezing after the adjustment
! Therfore, SLC will be recomputed, SMC remained unchanged
! temp: unfrozen ==> frozen

 rsoiltype(4) = 5.
 mask(4) = 1
 stc_bck(4,:) =  280.0

 smc_anl(4,:) = .25
 slc_anl(4,:) = .25
 stc_anl(4,1:lsoil_incr) =  271.0

 call apply_land_da_adjustments_soil(lsoil_incr, lsm, isot, ivegsrc,lensfc, &
           lsoil, rsoiltype, mask, stc_bck, stc_anl, smc_anl, slc_anl, &
           stc_updated, slc_updated,zsoil)

 do l = 1,lsoil_incr
  if (abs(smc_anl(1,l) - 0.25) > EPSILON) stop 2
  if (abs(slc_anl(1,l) - 0.25) > EPSILON) stop 3
  if (abs(stc_anl(1,l) - 281.) > EPSILON) stop 4
 enddo

 do l = 1,lsoil_incr
  if (abs(smc_anl(2,l) - 0.25) > EPSILON) stop 5
  if (abs(slc_anl(2,l) - 0.25) > EPSILON) stop 6
  if (abs(stc_anl(2,l) - 275.) > EPSILON) stop 7
 enddo

 do l = 1,lsoil_incr
  if (abs(smc_anl(3,l) - 0.25) > EPSILON) stop 8
  if (abs(slc_anl(3,l) - 0.112) > EPSILON) stop 9
  if (abs(stc_anl(3,l) - 271.) > EPSILON) stop 10
 enddo

 do l = 1,lsoil_incr
  if (abs(smc_anl(4,l) - 0.25) > EPSILON) stop 11
  if (abs(slc_anl(4,l) - 0.112) > EPSILON) stop 12
  if (abs(stc_anl(4,l) - 271.) > EPSILON) stop 13
 enddo

 call mpi_finalize(ierr)

 deallocate(rsoiltype,stc_bck,smc_anl,slc_anl,stc_anl,mask)

 if (my_rank .eq. 0) print*, "OK"
 if (my_rank .eq. 0) print*, "SUCCESS!"

end program ftst_land_increments
