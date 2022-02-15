!> @file
!! @brief Define atmospheric target data variables.
!! @author George Gayno NCEP/EMC

!> Module to hold variables and ESMF fields associated 
!! with the target grid atmospheric data.
!!
!! @author George Gayno NCEP/EMC
 module atmosphere_target_data

 use esmf

 implicit none

 private

 integer, public                    :: lev_target       !< Number of vertical levels.
 integer, public                    :: levp1_target     !< Number of vertical levels plus 1.
 integer, public                    :: nvcoord_target   !< Number of vertical coordinate variables.

 real(esmf_kind_r8), allocatable, public :: vcoord_target(:,:)  !< Vertical coordinate.

 type(esmf_field), public               :: delp_target_grid !< Pressure thickness.
 type(esmf_field), public               :: dzdt_target_grid !< Vertical velocity.
 type(esmf_field), public               :: ps_target_grid !< Surface pressure.
 type(esmf_field), public               :: temp_target_grid !< Temperautre.
 type(esmf_field), allocatable, public  :: tracers_target_grid(:) !< Tracers.
 type(esmf_field), public               :: u_s_target_grid !< U-wind, 'south' edge of grid cell.
 type(esmf_field), public               :: v_s_target_grid !< V-wind, 'south' edge of grid cell.
 type(esmf_field), public               :: u_w_target_grid !< U-wind, 'west' edge of grid cell.
 type(esmf_field), public               :: v_w_target_grid !< V-wind, 'west' edge of grid cell.
 type(esmf_field), public               :: zh_target_grid !< 3-d height.
 type(esmf_field), public               :: qnifa_climo_target_grid !< Number concentration of ice
                                           !! friendly aerosols.
 type(esmf_field), public               :: qnwfa_climo_target_grid !< Number concentration of water
                                           !! friendly aerosols.

 public :: cleanup_atmosphere_target_data

 contains

!> Free up memory for fields and variables in this module.
!!
!! @author George.Gayno NOAA/EMC
 subroutine cleanup_atmosphere_target_data

 use program_setup, only      : num_tracers

 implicit none

 integer                     :: i, rc

 print*,"- DESTROY TARGET GRID ATMOSPHERIC FIELDS."

 if (ESMF_FieldIsCreated(delp_target_grid)) call ESMF_FieldDestroy(delp_target_grid, rc=rc)
 if (ESMF_FieldIsCreated(dzdt_target_grid)) call ESMF_FieldDestroy(dzdt_target_grid, rc=rc)
 if (ESMF_FieldIsCreated(ps_target_grid))   call ESMF_FieldDestroy(ps_target_grid, rc=rc)
 if (ESMF_FieldIsCreated(temp_target_grid)) call ESMF_FieldDestroy(temp_target_grid, rc=rc)
 if (ESMF_FieldIsCreated(u_s_target_grid))  call ESMF_FieldDestroy(u_s_target_grid, rc=rc)
 if (ESMF_FieldIsCreated(v_s_target_grid)) call ESMF_FieldDestroy(v_s_target_grid, rc=rc)
 if (ESMF_FieldIsCreated(u_w_target_grid)) call ESMF_FieldDestroy(u_w_target_grid, rc=rc)
 if (ESMF_FieldIsCreated(v_w_target_grid)) call ESMF_FieldDestroy(v_w_target_grid, rc=rc)
 if (ESMF_FieldIsCreated(zh_target_grid))  call ESMF_FieldDestroy(zh_target_grid, rc=rc)

 do i = 1, num_tracers
   if (ESMF_FieldIsCreated(tracers_target_grid(i))) call ESMF_FieldDestroy(tracers_target_grid(i), rc=rc)
 enddo

 if (allocated (tracers_target_grid)) deallocate(tracers_target_grid)

 if (ESMF_FieldIsCreated(qnifa_climo_target_grid)) then
   call ESMF_FieldDestroy(qnifa_climo_target_grid, rc=rc)
 endif

 if (ESMF_FieldIsCreated(qnwfa_climo_target_grid)) then
   call ESMF_FieldDestroy(qnwfa_climo_target_grid, rc=rc)
 endif

 if (allocated (vcoord_target)) deallocate(vcoord_target)

 end subroutine cleanup_atmosphere_target_data

 end module atmosphere_target_data
