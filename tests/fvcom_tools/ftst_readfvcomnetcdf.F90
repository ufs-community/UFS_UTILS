 program readfvcomnetcdf

! Unit test for the fvcom_tools routines. 
!
! Reads in a 5x5 user generated grid to see if  
! file is read in correctly. Grid generated via
! netCDF4 python library. Expected values come
! from those used in the python generation
! script.
!
! Author David Wright


 use module_ncio, only: ncio
 use module_nwp, only: fcst_nwp

 implicit none

 integer, parameter :: NUM_VALUES=2 !number of values to check

 real, parameter :: EPSILON=0.0001 !error difference to check against

 integer :: nlat, nlon, t1, t2
 integer :: fv3lon, fv3lat, fv3times
 integer :: lbclon, lbclat, lbctimes
 integer :: t2_expected
 logical :: fv3_exists, fvcom_exists

 integer :: lat_lon_expected_values(NUM_VALUES) !expected number of lat/lons
 integer :: fv3mask_expected(NUM_VALUES) !expected fv3 mask values
 integer :: fv3sst_expected(NUM_VALUES) !expected fv3 sst values
 real :: fv3ice_expected(NUM_VALUES) !expected fv3 ice values
 real :: fv3iceT_expected(NUM_VALUES) !expected fv3 ice temp values
 real :: fv3zorl_expected(NUM_VALUES) !expected fv3 surface roughness values
 real :: fv3hice_expected(NUM_VALUES) !exepcted fv3 ice thickness values
 integer :: lbcmask_expected(NUM_VALUES) !expected fvcom mask values
 real :: lbcsst_expected(NUM_VALUES) !expected fvcom sst values
 real :: lbcice_expected(NUM_VALUES) !expected fvcom ice values
 real :: lbciceT_expected(NUM_VALUES) !expected fvcom ice temp values
 real :: lbcvice_expected(NUM_VALUES) !expected fvcom ice thickness values

! Create allocabable arrays to read from .nc files
 real, allocatable :: fv3ice(:,:), fv3sst(:,:)
 real, allocatable :: fv3sfcT(:,:), fv3mask(:,:)
 real, allocatable :: fv3iceT(:,:), fv3sfcTl(:,:)
 real, allocatable :: fv3zorl(:,:), fv3hice(:,:)
 real, allocatable :: lbcice(:,:), lbcsst(:,:)
 real, allocatable :: lbcsfcT(:,:), lbcmask(:,:)
 real, allocatable :: lbciceT(:,:), lbchice(:,:)
 real, allocatable :: lbczorl(:,:)
! Expected values from the dummy files
 data lat_lon_expected_values /5, 5/
 data fv3mask_expected /1, 0/
 data fv3sst_expected /1, 0/
 data fv3ice_expected /.1, 0/
 data fv3iceT_expected /.1, 0/
 data fv3zorl_expected /1.1, 0/
 data fv3hice_expected /.1, 0/
 data lbcmask_expected /1, 0/
 data lbcsst_expected /1, -99.99999/
 data lbcice_expected /1, -99.99999/
 data lbciceT_expected /1, -99.99999/
 data lbcvice_expected /1, -99.99999/
 data t2_expected /2 / !expect second time index from fvcom file

 type(ncio) :: geo !grid data object
 type(fcst_nwp) :: fcst !object to read data into

 character(len=180) :: fv3file !fv3 file name
 character(len=180) :: fvcomfile !fvcom file name
 character(len=180) :: wcstart !warm or cold start
 character(len=180) :: inputFVCOMselStr !time str in fvcom file


 print*,"Starting test of fvcom_tools."
!Set default file names, cold start, and time str
 fv3file = './data/sfcdata_unittest.nc'
 fvcomfile = './data/fvcom_unittest.nc'
 wcstart = 'cold'
 inputFVCOMselStr = '3333-44-55T66:77:88.000000'
 t1 = 1

!If files do not exist, stop
 INQUIRE(FILE=trim(fv3file), EXIST=fv3_exists)
 if(.not.fv3_exists) stop 1
 INQUIRE(FILE=trim(fvcomfile), EXIST=fvcom_exists)
 if(.not.fvcom_exists) stop 2
!Open grid file and read in number of lat/lon points
 call geo%open(trim(fv3file), 'r', 200)
 call geo%get_dim("xaxis_1",nlon)
 call geo%get_dim("yaxis_1",nlat)
 call geo%close
!If file does not have expected lat/lon points, stop
 if (abs(nlon - lat_lon_expected_values(2)) > EPSILON) stop 3
 if (abs(nlat - lat_lon_expected_values(1)) > EPSILON) stop 4

 allocate(fv3ice(nlon,nlat))
 allocate(fv3sfcT(nlon,nlat))
 allocate(fv3sst(nlon,nlat))
 allocate(fv3mask(nlon,nlat))
 allocate(fv3iceT(nlon,nlat))
 allocate(fv3sfcTl(nlon,nlat))
 allocate(fv3zorl(nlon,nlat))
 allocate(fv3hice(nlon,nlat))

 allocate(lbcice(nlon,nlat))
 allocate(lbcsfcT(nlon,nlat))
 allocate(lbcsst(nlon,nlat))
 allocate(lbcmask(nlon,nlat))
 allocate(lbciceT(nlon,nlat))
 allocate(lbczorl(nlon,nlat))
 allocate(lbchice(nlon,nlat))

!Initialize and read in fv3 sfc data
 call fcst%initial('FV3LAM',wcstart)
 call fcst%read_n(trim(fv3file),'FV3LAM',wcstart,fv3lon,fv3lat,fv3times,t1,fv3mask,fv3sst,fv3ice,fv3sfcT,fv3iceT,fv3sfcTl,fv3zorl,fv3hice,1,nlat)
 call fcst%finish('FV3LAM',wcstart)
!If variables in fv3 sfc file do not match expected, stop
 if (abs(fv3mask(1,1) - fv3mask_expected(1)) > EPSILON) stop 5
 if (abs(fv3mask(5,5) - fv3mask_expected(2)) > EPSILON) stop 6

 if (abs(fv3sst(1,1) - fv3sst_expected(1)) > EPSILON) stop 7
 if (abs(fv3sst(5,5) - fv3sst_expected(2)) > EPSILON) stop 8
 
 if (abs(fv3ice(1,1) - fv3ice_expected(1)) > EPSILON) stop 7
 if (abs(fv3ice(5,5) - fv3ice_expected(2)) > EPSILON) stop 8

 if (abs(fv3iceT(1,1) - fv3iceT_expected(1)) > EPSILON) stop 9
 if (abs(fv3iceT(5,5) - fv3iceT_expected(2)) > EPSILON) stop 10

 if (abs(fv3zorl(1,1) - fv3zorl_expected(1)) > EPSILON) stop 11
 if (abs(fv3zorl(5,5) - fv3zorl_expected(2)) > EPSILON) stop 12

 if (abs(fv3hice(1,1) - fv3hice_expected(1)) > EPSILON) stop 13
 if (abs(fv3hice(5,5) - fv3hice_expected(2)) > EPSILON) stop 14

!Initialize and read in fvcom data
 call fcst%initial(' FVCOM',wcstart)
 call fcst%get_time_ind(trim(fvcomfile),inputFVCOMselStr,t2)
!If second time index is not returned, stop
 if  (abs(t2 - t2_expected) > EPSILON) stop 15
 call fcst%read_n(trim(fvcomfile),' FVCOM',wcstart,lbclon,lbclat,lbctimes,t2,lbcmask,lbcsst,lbcice,lbcsfcT,lbciceT,fv3sfcTl,lbczorl,lbchice,1,nlat)
 call fcst%finish(' FVCOM',wcstart)
!If variables in fvcom file do not match expected, stop
 if (abs(lbcmask(1,1) - lbcmask_expected(1)) > EPSILON) stop 16
 if (abs(lbcmask(5,5) - lbcmask_expected(2)) > EPSILON) stop 17

 if (abs(lbcsst(1,1) - lbcsst_expected(1)) > EPSILON) stop 18
 if (abs(lbcsst(5,5) - lbcsst_expected(2)) > EPSILON) stop 19
 
 if (abs(lbcice(1,1) - lbcice_expected(1)) > EPSILON) stop 20
 if (abs(lbcice(5,5) - lbcice_expected(2)) > EPSILON) stop 21

 if (abs(lbciceT(1,1) - lbciceT_expected(1)) > EPSILON) stop 22
 if (abs(lbciceT(5,5) - lbciceT_expected(2)) > EPSILON) stop 23

 if (abs(lbchice(1,1) - lbcvice_expected(1)) > EPSILON) stop 24
 if (abs(lbchice(5,5) - lbcvice_expected(2)) > EPSILON) stop 25

 print*,"OK"

 print*,"SUCCESS!"

 end program readfvcomnetcdf

