program ftst_land_increments

! Test "apply_land_da_adjustments" using sample points.
!
! author: George Gayno (george.gayno@noaa.gov)

 use mpi
 use land_increments

 implicit none

 integer :: my_rank, ierr, lsm, isot, ivegsrc
 integer :: lensfc, lsoil, l

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
 
 lensfc= 3  ! Number of test points.
 
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
! the increments to STC and SLC will be retained.

 rsoiltype(1) = 5.
 mask(1) = 1
 stc_bck(1,:) =  280.0

 smc_anl(1,:) = .25
 slc_anl(1,:) = .25
 stc_anl(1,1:3) =  281.0 ! DA only updates 3 layers

! Point 2 is below freezing before the adjustment
! and above freezing after the adjustment. Therefore,
! the increment to STC will be removed, and SMC / SLC 
! are unchanged.

 rsoiltype(2) = 5.
 mask(2) = 1
 stc_bck(2,:) =  270.0

 smc_anl(2,:) = .25
 slc_anl(2,:) = .20
 stc_anl(2,1:3) =  274.0 

! Point 3 freezes. Therefore,
! the increment to STC will be removed, and SMC / SLC
! are unchanged.

 rsoiltype(3) = 5.
 mask(3) = 1
 stc_bck(3,:) =  274.0

 smc_anl(3,:) = .25
 slc_anl(3,:) = .25
 stc_anl(3,1:3) =  271.0

 call apply_land_da_adjustments_soil(lsm, isot, ivegsrc,lensfc, &
           lsoil, rsoiltype, mask, stc_bck, stc_anl, smc_anl, slc_anl, &
           stc_updated, slc_updated,zsoil)

 do l = 1, 3
   if (abs(smc_anl(1,l) - 0.25) > EPSILON) stop 2
   if (abs(slc_anl(1,l) - 0.25) > EPSILON) stop 4
   if (abs(stc_anl(1,l) - 281.) > EPSILON) stop 3
 enddo

 do l = 1, 3
   if (abs(smc_anl(2,l) - 0.25) > EPSILON) stop 6
   if (abs(slc_anl(2,l) - 0.20) > EPSILON) stop 8
   if (abs(stc_anl(2,l) - 270.) > EPSILON) stop 5
 enddo

 do l = 1, 3
   if (abs(smc_anl(3,l) - 0.25) > EPSILON) stop 10
   if (abs(slc_anl(3,l) - 0.25) > EPSILON) stop 12
   if (abs(stc_anl(3,l) - 274.) > EPSILON) stop 11
 enddo

 call mpi_finalize(ierr)

 deallocate(rsoiltype,stc_bck,smc_anl,slc_anl,stc_anl,mask)

 if (my_rank .eq. 0) print*, "OK"
 if (my_rank .eq. 0) print*, "SUCCESS!"

end program ftst_land_increments
