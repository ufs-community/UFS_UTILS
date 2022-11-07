!> @file
!! @brief Define target grid surface data variables.
!! @author George Gayno NCEP/EMC

!> Module to hold ESMF fields associated
!! with the target grid surface data.
!!
!! @author George Gayno NCEP/EMC
 module surface_target_data

 use esmf

 implicit none

 private

! surface fields (not including nst)
 type(esmf_field), public   :: canopy_mc_target_grid
                                       !< Canopy moisture content.
 type(esmf_field), public   :: f10m_target_grid
                                       !< log((z0+10)*1/z0)
                                       !< See sfc_diff.f for details.
 type(esmf_field), public   :: ffmm_target_grid
                                       !< log((z0+z1)*1/z0)
                                       !< See sfc_diff.f for details.
 type(esmf_field), public   :: q2m_target_grid
                                       !< 2-m specific humidity.
 type(esmf_field), public   :: seaice_depth_target_grid
                                       !< Sea ice depth.
 type(esmf_field), public   :: seaice_fract_target_grid
                                       !< Sea ice fraction.
 type(esmf_field), public   :: seaice_skin_temp_target_grid
                                       !< Sea ice skin temperature.
 type(esmf_field), public   :: skin_temp_target_grid
                                       !< Skin temperature/sst.
 type(esmf_field), public   :: srflag_target_grid
                                       !< Snow/rain flag.
 type(esmf_field), public   :: snow_liq_equiv_target_grid
                                       !< Liquid equivalent snow depth.
 type(esmf_field), public   :: snow_depth_target_grid
                                       !< Physical snow depth.
 type(esmf_field), public   :: soil_temp_target_grid
                                       !< 3-d soil temperature.
 type(esmf_field), public   :: soilm_liq_target_grid
                                       !< 3-d liquid soil moisture.
 type(esmf_field), public   :: soilm_tot_target_grid
                                       !< 3-d total soil moisture.
 type(esmf_field), public   :: t2m_target_grid
                                       !< 2-m temperatrure.
 type(esmf_field), public   :: tprcp_target_grid
                                       !< Precipitation.
 type(esmf_field), public   :: ustar_target_grid
                                       !< Friction velocity.
 type(esmf_field), public   :: z0_target_grid
                                       !< Roughness length.
 type(esmf_field), public   :: lai_target_grid
                                       !< Leaf area index.

! nst fields
 type(esmf_field), public   :: c_d_target_grid
                                       !< Coefficient 2 to calculate d(tz)/d(ts).
 type(esmf_field), public   :: c_0_target_grid
                                       !< Coefficient 1 to calculate d(tz)/d(ts).
 type(esmf_field), public   :: d_conv_target_grid
                                       !< Thickness of free convection layer.
 type(esmf_field), public   :: dt_cool_target_grid
                                       !< Sub-layer cooling amount.
 type(esmf_field), public   :: ifd_target_grid
                                       !< Model mode index. 0-diurnal model not
                                       !< started; 1-diurnal model started.
 type(esmf_field), public   :: qrain_target_grid
                                       !< Sensible heat flux due to rainfall.
 type(esmf_field), public   :: tref_target_grid
                                       !< Reference temperature.
 type(esmf_field), public   :: w_d_target_grid
                                       !< Coefficient 4 to calculate d(tz)/d(ts).
 type(esmf_field), public   :: w_0_target_grid
                                       !< Coefficient 3 to calculate d(tz)/d(ts).
 type(esmf_field), public   :: xs_target_grid
                                       !< Salinity content in diurnal
                                       !< thermocline layer.
 type(esmf_field), public   :: xt_target_grid
                                       !< Heat content in diurnal thermocline
                                       !< layer.
 type(esmf_field), public   :: xu_target_grid
                                       !< u-current content in diurnal
                                       !< thermocline layer.
 type(esmf_field), public   :: xv_target_grid
                                       !< v-current content in diurnal
                                       !< thermocline layer.
 type(esmf_field), public   :: xz_target_grid
                                       !< Diurnal thermocline layer thickness.
 type(esmf_field), public   :: xtts_target_grid
                                       !< d(xt)/d(ts).
 type(esmf_field), public   :: xzts_target_grid
                                       !< d(xz)/d(ts).
 type(esmf_field), public   :: z_c_target_grid
                                       !< Sub-layer cooling thickness.
 type(esmf_field), public   :: zm_target_grid
                                       !< Oceanic mixed layer depth.

 public :: cleanup_target_nst_data
 public :: cleanup_target_sfc_data

 contains

!> Free up memory once the target grid surface fields are
!! no longer needed.
!!
!! @author George Gayno NOAA/EMC
 subroutine cleanup_target_sfc_data

 implicit none

 integer                     :: rc

 print*,"- DESTROY TARGET GRID SURFACE FIELDS."
 call ESMF_FieldDestroy(t2m_target_grid, rc=rc)
 call ESMF_FieldDestroy(q2m_target_grid, rc=rc)
 call ESMF_FieldDestroy(tprcp_target_grid, rc=rc)
 call ESMF_FieldDestroy(f10m_target_grid, rc=rc)
 call ESMF_FieldDestroy(ffmm_target_grid, rc=rc)
 call ESMF_FieldDestroy(ustar_target_grid, rc=rc)
 call ESMF_FieldDestroy(snow_liq_equiv_target_grid, rc=rc)
 call ESMF_FieldDestroy(snow_depth_target_grid, rc=rc)
 call ESMF_FieldDestroy(seaice_fract_target_grid, rc=rc)
 call ESMF_FieldDestroy(seaice_depth_target_grid, rc=rc)
 call ESMF_FieldDestroy(seaice_skin_temp_target_grid, rc=rc)
 call ESMF_FieldDestroy(srflag_target_grid, rc=rc)
 call ESMF_FieldDestroy(skin_temp_target_grid, rc=rc)
 call ESMF_FieldDestroy(canopy_mc_target_grid, rc=rc)
 call ESMF_FieldDestroy(lai_target_grid,rc=rc)
 call ESMF_FieldDestroy(z0_target_grid, rc=rc)
 call ESMF_FieldDestroy(soil_temp_target_grid, rc=rc)
 call ESMF_FieldDestroy(soilm_tot_target_grid, rc=rc)
 call ESMF_FieldDestroy(soilm_liq_target_grid, rc=rc)

 end subroutine cleanup_target_sfc_data

!> Free up memory once the target grid nst fields are
!! no longer needed.
!!
!! @author George Gayno NOAA/EMC
 subroutine cleanup_target_nst_data

 implicit none

 integer                            :: rc

 print*,"- DESTROY TARGET GRID NST DATA."

 call ESMF_FieldDestroy(c_d_target_grid, rc=rc)
 call ESMF_FieldDestroy(c_0_target_grid, rc=rc)
 call ESMF_FieldDestroy(d_conv_target_grid, rc=rc)
 call ESMF_FieldDestroy(dt_cool_target_grid, rc=rc)
 call ESMF_FieldDestroy(ifd_target_grid, rc=rc)
 call ESMF_FieldDestroy(qrain_target_grid, rc=rc)
 call ESMF_FieldDestroy(tref_target_grid, rc=rc)
 call ESMF_FieldDestroy(w_d_target_grid, rc=rc)
 call ESMF_FieldDestroy(w_0_target_grid, rc=rc)
 call ESMF_FieldDestroy(xs_target_grid, rc=rc)
 call ESMF_FieldDestroy(xt_target_grid, rc=rc)
 call ESMF_FieldDestroy(xu_target_grid, rc=rc)
 call ESMF_FieldDestroy(xv_target_grid, rc=rc)
 call ESMF_FieldDestroy(xz_target_grid, rc=rc)
 call ESMF_FieldDestroy(xtts_target_grid, rc=rc)
 call ESMF_FieldDestroy(xzts_target_grid, rc=rc)
 call ESMF_FieldDestroy(z_c_target_grid, rc=rc)
 call ESMF_FieldDestroy(zm_target_grid, rc=rc)

 end subroutine cleanup_target_nst_data

 end module surface_target_data
