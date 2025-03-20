!> Subroutine to read in namelists, and convert
!! values into setupgrid.
!! Also fills in some values, and tests have all
!! needed vals, according to the selected grid type.
!! 
!! @param[in] unt          file unit
!! @param[in] namel        options: input or output
!! @param[out] grid_setup  data structure with grid details

 subroutine readin_setup(unt,namel,grid_setup)

 use grids_IO, only     : grid_setup_type
 use utilities, only    : error_handler

 implicit none

 ! INPUTS
 integer, intent(in) :: unt
 character(*), intent(in) :: namel
 ! OUTPUTS
 type(grid_setup_type), intent(out) :: grid_setup

 character(len=7)   :: gridtype
 character(len=100) :: fname, fname_mask, fname_coord
 character(len=100) :: dir, dir_mask, dir_coord
 character(len=4)   :: default_str="NULL"
 integer            :: ires, jres
 integer            :: ierr

 namelist /input/  fname, dir, &
                   gridtype, &
                   fname_mask, dir_mask, &
                   fname_coord, dir_coord, &
                   ires, jres

 namelist /output/  fname, dir, &
                    gridtype, &
                    fname_mask, dir_mask, &
                    fname_coord, dir_coord,&
                    ires, jres

 ! set defaults
 fname = default_str
 dir = default_str
 fname_mask = default_str
 dir_mask = default_str
 fname_coord = default_str
 dir_coord = default_str
 ires = 0
 jres = 0

 select case (namel)
 case ("input")
     read(unt, nml=input, iostat=ierr)
     if (ierr /= 0) call error_handler("READING input NAMELIST.", ierr)
 case ("output")
     read(unt, nml=output, iostat=ierr)
     if (ierr /= 0) call error_handler("READING output NAMELIST.", ierr)
 case default
     call error_handler("unknown namel in readin_setup", 1)
 end select

 grid_setup%descriptor = gridtype

 grid_setup%dir = dir
 grid_setup%fname = fname

 ! set-up mask details, based on file type
 select case (gridtype)
 case ("fv3_rst") ! for history file and restarts, use veg type
     grid_setup%dir_mask = dir_mask ! get from a fix file
     grid_setup%fname_mask = fname_mask
     grid_setup%mask_variable(1) =  "vegetation_type" ! if getting from fix file
 case ("gau_inc") ! gsi-output incr files only, use calculated mask
     if (trim(fname_mask) == default_str) then ! if not specified, use input file
         grid_setup%dir_mask = dir
         grid_setup%fname_mask = fname
     else
         grid_setup%dir_mask = dir_mask
         grid_setup%fname_mask = fname_mask
     endif
     grid_setup%mask_variable(1) =  "soilsnow_mask  "
 case default
     call error_handler("unknown gridtype in readin_setup", 1)
 end select

 grid_setup%dir_coord = dir_coord
 grid_setup%fname_coord = fname_coord
 grid_setup%ires = ires
 grid_setup%jres = jres

 end subroutine readin_setup
