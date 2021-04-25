module atm_target_data_mod
  use esmf
  implicit none

  private
  
  integer, public                    :: lev_target       !< num vertical levels
  integer, public                    :: levp1_target     !< num levels plus 1
  integer, public                    :: nvcoord_target   !< num vertical coordinate variables

  real(esmf_kind_r8), allocatable, public :: vcoord_target(:,:)  !< vertical coordinate

  type(esmf_field), public               :: delp_target_grid !< pressure thickness
  type(esmf_field), public               :: dzdt_target_grid !< vertical velocity

  type(esmf_field), allocatable, public  :: tracers_target_grid(:) !< tracers
  type(esmf_field), public               :: ps_target_grid !< surface pressure
  type(esmf_field), public               :: temp_target_grid !< temperautre
  type(esmf_field), public               :: u_s_target_grid !< u-wind, 'south' edge
  type(esmf_field), public               :: v_s_target_grid !< v-wind, 'south' edge
  type(esmf_field), public               :: u_w_target_grid !< u-wind, 'west' edge
  type(esmf_field), public               :: v_w_target_grid !< v-wind, 'west' edge

  type(esmf_field), public               :: zh_target_grid !< 3-d height

  type(esmf_field), public               :: qnifa_climo_target_grid !< number concentration of ice
                                                                    !! friendly aerosols on target 
                                                                    !! horiz/vert grid.
  
  type(esmf_field), public               :: qnwfa_climo_target_grid !< number concentration of water
                                                                    !! friendly aerosols on target 
                                                                    !! horiz/vert grid.
contains
end module atm_target_data_mod
