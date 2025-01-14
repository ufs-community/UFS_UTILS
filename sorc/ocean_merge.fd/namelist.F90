!> Read program namelist.
!!
!! @param[out] ocean_mask_dir Directory containing MOM6 ocean mask file.
!! @param[out] lake_mask_dir Directory containing the lake mask file.
!! @param[out] out_dir Directory where output file will be written.
!! @param[out] atmres Atmosphere grid resolution.
!! @param[out] ocnres Ocean grid resolution.
!! @param[out] binary_lake When '1', treat lake fraction as either 0 or 1. Otherwise, 
!! it is a fraction.
!!
!! @author Rahul Mahajan
!! @author Sanath Kumar
subroutine read_nml(ocean_mask_dir, lake_mask_dir, atmres,ocnres,out_dir,binary_lake)

  implicit none

  integer :: unit=7, io_status

  character(len=200), intent(out) :: ocean_mask_dir
  character(len=200), intent(out) :: lake_mask_dir
  character(len=200), intent(out) :: out_dir
  character(len=10),  intent(out) :: atmres,ocnres
  integer, intent(out):: binary_lake

  namelist/mask_nml/ocean_mask_dir, lake_mask_dir, atmres, ocnres,out_dir,binary_lake
  open(unit=unit, file='input.nml', iostat=io_status )
  read(unit,mask_nml, iostat=io_status )
  close(unit)
  if (io_status > 0) then
        print *,'FATAL ERROR reading input.nml' 
        call handle_err(-1)
  end if      
end subroutine read_nml
