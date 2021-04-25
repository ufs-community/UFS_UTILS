module sfc_target_data_mod
  use esmf
  implicit none

  private

  ! surface fields (not including nst)
  type(esmf_field), public           :: canopy_mc_target_grid
  !< canopy moisture content
  type(esmf_field), public           :: f10m_target_grid
  !< log((z0+10)*1/z0)
  !< See sfc_diff.f for details
  type(esmf_field), public           :: ffmm_target_grid
  !< log((z0+z1)*1/z0)
  !< See sfc_diff.f for details
  type(esmf_field), public           :: q2m_target_grid
  !< 2-m specific humidity
  type(esmf_field), public           :: seaice_depth_target_grid
  !< sea ice depth
  type(esmf_field), public           :: seaice_fract_target_grid
  !< sea ice fraction
  type(esmf_field), public           :: seaice_skin_temp_target_grid
  !< sea ice skin temperature
  type(esmf_field), public           :: skin_temp_target_grid
  !< skin temperature/sst
  type(esmf_field), public           :: srflag_target_grid
  !< snow/rain flag
  type(esmf_field), public           :: snow_liq_equiv_target_grid
  !< liquid equiv snow depth
  type(esmf_field), public           :: snow_depth_target_grid
  !< physical snow depth
  type(esmf_field), public           :: soil_temp_target_grid
  !< 3-d soil temperature
  type(esmf_field), public           :: soilm_liq_target_grid
  !< 3-d liquid soil moisture
  type(esmf_field), public           :: soilm_tot_target_grid
  !< 3-d total soil moisture
  type(esmf_field), public           :: t2m_target_grid
  !< 2-m temperatrure
  type(esmf_field), public           :: tprcp_target_grid
  !< precip
  type(esmf_field), public           :: ustar_target_grid
  !< friction velocity
  type(esmf_field), public           :: z0_target_grid
  !< roughness length
  type(esmf_field), public           :: lai_target_grid
  !< leaf area index

  ! nst fields
  type(esmf_field), public           :: c_d_target_grid
  !< Coefficient 2 to calculate d(tz)/d(ts)
  type(esmf_field), public           :: c_0_target_grid
  !< Coefficient 1 to calculate d(tz)/d(ts)
  type(esmf_field), public           :: d_conv_target_grid
  !< Thickness of free convection layer
  type(esmf_field), public           :: dt_cool_target_grid
  !< Sub-layer cooling amount
  type(esmf_field), public           :: ifd_target_grid
  !< Model mode index. 0-diurnal model not
  !< started; 1-diurnal model started.
  type(esmf_field), public           :: qrain_target_grid
  !< Sensible heat flux due to rainfall
  type(esmf_field), public           :: tref_target_grid
  !< reference temperature
  type(esmf_field), public           :: w_d_target_grid
  !< Coefficient 4 to calculate d(tz)/d(ts)
  type(esmf_field), public           :: w_0_target_grid
  !< Coefficient 3 to calculate d(tz)/d(ts)
  type(esmf_field), public           :: xs_target_grid
  !< Salinity content in diurnal
  !< thermocline layer
  type(esmf_field), public           :: xt_target_grid
  !< Heat content in diurnal thermocline
  !< layer
  type(esmf_field), public           :: xu_target_grid
  !< u-current content in diurnal
  !< thermocline layer
  type(esmf_field), public           :: xv_target_grid
  !< v-current content in diurnal
  !< thermocline layer
  type(esmf_field), public           :: xz_target_grid
  !< Diurnal thermocline layer thickness
  type(esmf_field), public           :: xtts_target_grid
  !< d(xt)/d(ts)
  type(esmf_field), public           :: xzts_target_grid
  !< d(xz)/d(ts)
  type(esmf_field), public           :: z_c_target_grid
  !< Sub-layer cooling thickness
  type(esmf_field), public           :: zm_target_grid
  !< Oceanic mixed layer depth

contains
end module sfc_target_data_mod



