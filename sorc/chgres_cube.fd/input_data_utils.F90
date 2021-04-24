!> @file
!! @brief Shared routines for input data modules
!! @author George Gayno NCEP/EMC

!> Read atmospheric, surface and nst data on the input grid.
!! Supported formats include fv3 tiled 'restart' files, fv3 tiled 
!! 'history' files, fv3 gaussian history files, spectral gfs
!! gaussian nemsio files, and spectral gfs sigio/sfcio files.
!!
!! Public variables are defined below: "input" indicates field
!! associated with the input grid.
!!
!! @author George Gayno NCEP/EMC
module input_data_utils_mod
  use esmf
  use netcdf
  use chgres_cube_utils_mod

  use program_setup, only: data_dir_input_grid, &
       sfc_files_input_grid, read_from_input

  private
  
  public :: read_fv3_grid_data_netcdf
  public :: handle_grib_error

  integer, public      :: lsoil_input=4  !< number of soil layers, no longer hardwired to allow
  !! for 7 layers of soil for the RUC LSM

  type(ESMF_Field), public :: terrain_input_grid

contains

  !> Read a record from a netcdf file
  !!
  !! @param [in] field  name of field to be read 
  !! @param [in] tile_num  grid tile number
  !! @param [in] imo i-dimension of field
  !! @param [in] jmo j-dimension of field
  !! @param [in] lmo number of vertical levels of field
  !! @param [out] sfcdata 1-d array containing field data
  !! @param [out] sfcdata_3d  3-d array containing field data
  !! @author George Gayno NCEP/EMC   
  SUBROUTINE READ_FV3_GRID_DATA_NETCDF(FIELD,TILE_NUM,IMO,JMO,LMO, &
       SFCDATA, SFCDATA_3D)

    IMPLICIT NONE

    CHARACTER(LEN=*),INTENT(IN)      :: FIELD

    INTEGER, INTENT(IN)   :: IMO, JMO, LMO, TILE_NUM

    REAL(ESMF_KIND_R8), INTENT(OUT), OPTIONAL     :: SFCDATA(IMO,JMO)
    REAL(ESMF_KIND_R8), INTENT(OUT), OPTIONAL     :: SFCDATA_3D(IMO,JMO,LMO)

    CHARACTER(LEN=256)    :: TILEFILE

    INTEGER               :: ERROR, NCID, ID_VAR

    TILEFILE = TRIM(DATA_DIR_INPUT_GRID) // "/" // TRIM(SFC_FILES_INPUT_GRID(TILE_NUM))

    PRINT*,'WILL READ ',TRIM(FIELD), ' FROM: ', TRIM(TILEFILE)

    ERROR=NF90_OPEN(TRIM(TILEFILE),NF90_NOWRITE,NCID)
    CALL NETCDF_ERR(ERROR, 'OPENING: '//TRIM(TILEFILE) )

    ERROR=NF90_INQ_VARID(NCID, FIELD, ID_VAR)
    CALL NETCDF_ERR(ERROR, 'READING FIELD ID' )

    IF (PRESENT(SFCDATA_3D)) THEN
       ERROR=NF90_GET_VAR(NCID, ID_VAR, SFCDATA_3D)
       CALL NETCDF_ERR(ERROR, 'READING FIELD' )
    ELSE
       ERROR=NF90_GET_VAR(NCID, ID_VAR, SFCDATA)
       CALL NETCDF_ERR(ERROR, 'READING FIELD' )
    ENDIF

    ERROR = NF90_CLOSE(NCID)

  END SUBROUTINE READ_FV3_GRID_DATA_NETCDF

  !> Handle GRIB2 read error based on the user selected
  !! method in the varmap file.
  !!
  !! @param [in] vname  grib2 variable name
  !! @param [in] lev    grib2 variable level
  !! @param [in] method  how missing data is handled
  !! @param [in] value   fill value for missing data
  !! @param [in] varnum  grib2 variable number
  !! @param [inout] iret  return status code
  !! @param [inout] var   4-byte array of corrected data
  !! @param [inout] var8  8-byte array of corrected data
  !! @param [inout] var3d 3-d array of corrected data
  !! @author Larissa Reames
  subroutine handle_grib_error(vname,lev,method,value,varnum, iret,var,var8,var3d)

    use, intrinsic :: ieee_arithmetic

    implicit none

    real(esmf_kind_r4), intent(in)    :: value
    real(esmf_kind_r4), intent(inout), optional :: var(:,:)
    real(esmf_kind_r8), intent(inout), optional :: var8(:,:)
    real(esmf_kind_r8), intent(inout), optional  :: var3d(:,:,:)

    character(len=20), intent(in)     :: vname, lev, method

    integer, intent(in)               :: varnum 
    integer, intent(inout)            :: iret

    iret = 0
    if (varnum == 9999) then
       print*, "WARNING: ", trim(vname), " NOT FOUND AT LEVEL ", lev, " IN EXTERNAL FILE ", &
            "AND NO ENTRY EXISTS IN VARMAP TABLE. VARIABLE WILL NOT BE USED."
       iret = 1

       return
    endif

    if (trim(method) == "skip" ) then
       print*, "WARNING: SKIPPING ", trim(vname), " IN FILE"
       read_from_input(varnum) = .false.
       iret = 1
    elseif (trim(method) == "set_to_fill") then
       print*, "WARNING: ,", trim(vname), " NOT AVILABLE AT LEVEL ", trim(lev), &
            ". SETTING EQUAL TO FILL VALUE OF ", value
       if(present(var)) var(:,:) = value
       if(present(var8)) var8(:,:) = value
       if(present(var3d)) var3d(:,:,:) = value
    elseif (trim(method) == "set_to_NaN") then
       print*, "WARNING: ,", trim(vname), " NOT AVILABLE AT LEVEL ", trim(lev), &
            ". SETTING EQUAL TO NaNs"
       if(present(var)) var(:,:) = ieee_value(var,IEEE_QUIET_NAN)
       if(present(var8)) var8(:,:) = ieee_value(var8,IEEE_QUIET_NAN)
       if(present(var3d)) var3d(:,:,:) = ieee_value(var3d,IEEE_QUIET_NAN)
    elseif (trim(method) == "stop") then
       call error_handler("READING "//trim(vname)// " at level "//lev//". TO MAKE THIS NON- &
            FATAL, CHANGE STOP TO SKIP FOR THIS VARIABLE IN YOUR VARMAP &
            FILE.", iret)
    else
       call error_handler("ERROR USING MISSING_VAR_METHOD. PLEASE SET VALUES IN" // &
            " VARMAP TABLE TO ONE OF: set_to_fill, set_to_NaN,"// &
            " , skip, or stop.", 1)
    endif

  end subroutine handle_grib_error

end module input_data_utils_mod
