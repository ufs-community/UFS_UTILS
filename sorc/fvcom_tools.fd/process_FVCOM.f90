! process_FVCOM.f90
! David Wright
! University of Michigan and GLERL
! 17 Aug 2020
!
! This is the code to put lake surface temperature and aerial ice
! concentration from GLERL-provided FVCOM forecast files (which have
! already been mapped to the FV3-LAM grid) into sfc_data.nc.  
!
! This script will take two variables from the command line:
! 1. Name of FV3 sfc data file (e.g. sfc_data.tile7.halo0.nc)
! 2. Name of FVCOM data file (e.g. fvcom.nc)
!
! To run the script, use the following example, modifying file
! names as needed:
!   ./fvcom_to_FV3.exe sfc_data.tile7.halo0.nc fvcom.nc
!
! Code is strongly based upon Eric James' (ESRL/GSL) work 
!  to update HRRR/WRF Great Lakes' temperature data with FVCOM.
!  Code also relies heavily on Ming Hu's ncio module.

program process_FVCOM

   use mpi
   use kinds, only: r_kind, i_kind, r_single
   use module_ncio, only: ncio
   use module_nwp, only: fcst_nwp

   implicit none

! MPI variables
  integer :: npe, mype, mypeLocal,ierror
!

!  New object-oriented declarations

   type(ncio) :: geo
   type(fcst_nwp) :: fcst

!  Grid variables

   character*180 :: geofile
   character*2 :: workPath
   character*1 :: char1

   integer :: MAP_PROJ, NLON, NLAT
   integer :: fv3lon, fv3lat, fv3times
   integer :: lbclon, lbclat, lbctimes
   integer :: i, j, t1, t2
   integer :: num_args, ix

   real :: rad2deg = 180.0/3.1415926
   real :: userDX, userDY, CEN_LAT, CEN_LON
   real :: userTRUELAT1, userTRUELAT2, MOAD_CEN_LAT, STAND_LON
   real :: truelat1, truelat2, stdlon, lat1, lon1, r_earth
   real :: knowni, knownj, dx
   real :: one, pi, deg2rad

   character(len=180) :: fv3file
   character(len=180) :: fvcomfile
   character(len=180), dimension(:), allocatable :: args

   real(r_kind), allocatable :: fv3ice(:,:), fv3sst(:,:)
   real(r_kind), allocatable :: fv3sfcT(:,:), fv3mask(:,:)
   real(r_kind), allocatable :: fv3iceT(:,:)
   real(r_kind), allocatable :: lbcice(:,:), lbcsst(:,:)
   real(r_kind), allocatable :: lbcsfcT(:,:), lbcmask(:,:)
   real(r_kind), allocatable :: lbciceT(:,:)   

!  Declare namelists
!  SETUP (general control namelist) :

   integer :: update_type

!   namelist/setup/update_type, t2

! MPI setup
  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(mpi_comm_world,npe,ierror)
  call MPI_COMM_RANK(mpi_comm_world,mype,ierror)
!
! NCEP LSF has to use all cores allocated to run this application 
! but this if check can make sure only one core run through the real code.
!
if(mype==0) then

!  Get file names from command line arguements
   num_args = command_argument_count()
   allocate(args(num_args))

   do ix = 1, num_args
       call get_command_argument(ix,args(ix))
   end do

   fv3file=trim(args(1))
   write(*,*) trim(fv3file)
   fvcomfile=trim(args(2))    
   write(*,*) trim(fvcomfile)

   t2 = 1
!  Obtain grid parameters

   workPath='./'
   write(geofile,'(a,a)') trim(workPath), trim(fv3file)
   write(*,*) 'sfc data file', trim(geofile)
   call geo%open(trim(geofile),'r',200)
   call geo%get_dim("xaxis_1",NLON)
   call geo%get_dim("yaxis_1",NLAT)
   write(*,*) 'NLON,NLAT:', NLON, NLAT
   call geo%close

   write(*,*) 'Finished reading sfc_data grid information.'
   write(*,*) ' '

!  Allocate variables for I/O

   allocate(fv3ice(nlon,nlat))
   allocate(fv3sfcT(nlon,nlat))
   allocate(fv3sst(nlon,nlat))
   allocate(fv3mask(nlon,nlat))
   allocate(fv3iceT(nlon,nlat))

   allocate(lbcice(nlon,nlat))
   allocate(lbcsfcT(nlon,nlat))
   allocate(lbcsst(nlon,nlat))
   allocate(lbcmask(nlon,nlat))
   allocate(lbciceT(nlon,nlat))

!  Read fv3 sfc_data.nc before update

!   fv3file='sfc_data.nc'
   t1=1

   call fcst%initial('FV3LAM')
   call fcst%list_initial
   call fcst%read_n(trim(fv3file),'FV3LAM',fv3lon,fv3lat,fv3times,t1,fv3mask,fv3sst,fv3ice,fv3sfcT,fv3iceT)
   call fcst%finish


   write(*,*) 'fv3times: ', fv3times
   write(*,*) 'time to use: ', t1

!  Read FVCOM input datasets

!   fvcomfile='fvcom.nc'
! Space infront of ' FVCOM' below is important!!
   call fcst%initial(' FVCOM')
   call fcst%list_initial
   call fcst%read_n(trim(fvcomfile),' FVCOM',lbclon,lbclat,lbctimes,t2,lbcmask,lbcsst,lbcice,lbcsfcT,lbciceT)
   call fcst%finish

!  Check that the dimensions match

   if (lbclon .ne. nlon .or. lbclat .ne. nlat) then
      write(*,*) 'ERROR: FVCOM/FV3 dimensions do not match:'
     write(*,*) 'lbclon: ', lbclon
     write(*,*) 'nlon: ', nlon
     write(*,*) 'lbclat: ', lbclat
     write(*,*) 'nlat: ', nlat
     stop 135
   endif

   write(*,*) 'lbctimes: ', lbctimes
   write(*,*) 'time to use: ', t2

!  Update with FVCOM fields and process
!   ice cover data. Ice fraction is set
!   to a minimum of 15% due to FV3-LAM
!   raising any value below 15% to 15%.

   do j=1,nlat
      do i=1,nlon
         if (lbcmask(i,j) > 0. .and. lbcsst(i,j) .ge. -90.0) then
            !If ice fraction below 15%, set to 0
            if (lbcice(i,j) < 0.15) then
              lbcice(i,j) = 0.0
            endif
            fv3ice(i,j) = lbcice(i,j)
            !If ice in FVCOM, but not in FV3-LAM, change to ice
            if (lbcice(i,j) > 0. .and. fv3mask(i,j) == 0.) then
              fv3mask(i,j) = 2.
            endif
            !If ice in FV3-LAM and not FVCOM, remove it from FV3-LAM
            if (fv3mask(i,j) == 2. .and. lbcice(i,j) == 0.) then
              fv3mask(i,j) = 0.
            endif
            fv3sst(i,j) = lbcsst(i,j) + 273.15
            fv3sfcT(i,j) = lbcsst(i,j) + 273.15
            fv3iceT(i,j) = lbcsst(i,j) + 273.15
             !If ice exists in FVCOM, change ice surface temp
            if (lbcice(i,j) > 0.) then
              fv3iceT(i,j) = lbciceT(i,j) + 273.15
            end if
         endif
      enddo
   enddo

! Write out sfc file again

   call geo%open(trim(fv3file),'w',300)
   call geo%replace_var("tsea",NLON,NLAT,fv3sst)
   call geo%replace_var("fice",NLON,NLAT,fv3ice)
   call geo%replace_var("slmsk",NLON,NLAT,fv3mask)
   call geo%replace_var("tisfc",NLON,NLAT,fv3iceT)
! Add_New_Var takes names of (Variable,Dim1,Dim2,Dim3,Long_Name,Units)
   call geo%add_new_var('glmsk','xaxis_1','yaxis_1','Time','glmsk','none')
   call geo%replace_var('glmsk',NLON,NLAT,lbcmask)
   call geo%close

   write(6,*) "=== LOWBC UPDATE SUCCESS ==="

endif ! mype==0

call MPI_FINALIZE(ierror)


end program process_FVCOM 
