!> @file
!!
!! @brief This is the code to put lake surface temperature and aerial
!! ice concentration from GLERL-provided FVCOM forecast files (which
!! have already been mapped to the FV3-LAM grid) into sfc_data.nc.

!> Put lake surface temperature and aerial ice concentration from
!! GLERL-provided FVCOM forecast files (which have already been mapped
!! to the FV3-LAM grid) into sfc_data.nc.
!!
!! This script will take four variables from the command line:
!! 1. Name of FV3 sfc data file (e.g. sfc_data.tile7.halo0.nc)
!! 2. Name of FVCOM data file (e.g. fvcom.nc)
!! 3. "warm" or "cold" start. "warm" start will read in
!!    sfc_data.nc files generated from a restart of UFS-SRW.
!!    "cold" start will read in sfc_data.nc files generated
!!    from chgres_cube.
!! 4. String of time slice to use in the fvcom.nc file. This string
!!    should match exactly what is in the Times variable of the .nc file.
!! To run the script, use the following example, modifying file
!! names as needed:
!!   ./fvcom_to_FV3 sfc_data.tile7.halo0.nc fvcom.nc cold \
!!     2020-01-31T18:00:00.000000
!! Code is strongly based upon Eric James' (ESRL/GSL) work to update
!! HRRR/WRF Great Lakes' temperature data with FVCOM. Code also
!! relies heavily on Ming Hu's ncio module.
!!
!! @author David Wright, University of Michigan and GLERL
!! @date 17 Aug 2020
!! @return 0 for success, error code otherwise
!!
program process_FVCOM

   use mpi
   use kinds, only: r_kind, i_kind, r_single
   use module_ncio, only: ncio
   use module_nwp, only: fcst_nwp

   implicit none

! MPI variables
   integer :: npe, mype, mypeLocal,ierror
   integer,allocatable :: mpi_layout_begin(:),mpi_layout_end(:)
!

!  New object-oriented declarations

   type(ncio) :: geo
   type(fcst_nwp) :: fcst

!  Grid variables

   character*180 :: thisfv3file
   character*2 :: workPath
   character*1 :: char1

   integer :: MAP_PROJ, NLON, NLAT
   integer :: fv3lon, fv3lat, fv3times
   integer :: lbclon, lbclat, lbctimes
   integer :: i, j, t1, t2
   integer :: num_args, ix
   integer :: indexFVCOMsel

   real :: rad2deg = 180.0/3.1415926
   real :: userDX, userDY, CEN_LAT, CEN_LON
   real :: userTRUELAT1, userTRUELAT2, MOAD_CEN_LAT, STAND_LON
   real :: truelat1, truelat2, stdlon, lat1, lon1, r_earth
   real :: knowni, knownj, dx
   real :: one, pi, deg2rad
   real :: zero

   character(len=180) :: fv3file
   character(len=180) :: fvcomfile
   character(len=180) :: wcstart
   character(len=180) :: inputFVCOMselStr
   character(len=180), dimension(:), allocatable :: args
   integer            :: fv3_io_layout_y
   integer,allocatable :: fv3_layout_begin(:),fv3_layout_end(:)

   real(r_kind), allocatable :: fv3ice(:,:), fv3sst(:,:)
   real(r_kind), allocatable :: fv3sfcT(:,:), fv3mask(:,:)
   real(r_kind), allocatable :: fv3iceT(:,:), fv3sfcTl(:,:)
   real(r_kind), allocatable :: fv3zorl(:,:), fv3hice(:,:)
   real(r_kind), allocatable :: lbcice(:,:), lbcsst(:,:)
   real(r_kind), allocatable :: lbcsfcT(:,:), lbcmask(:,:)
   real(r_kind), allocatable :: lbciceT(:,:), lbczorl(:,:)
   real(r_kind), allocatable :: lbchice(:,:)   

!  Declare namelists
!  SETUP (general control namelist) :

   integer :: update_type
   character(len=4) :: clayout_y
   integer :: iy,iblock

!   namelist/setup/update_type, t2

! MPI setup
   call MPI_INIT(ierror)
   call MPI_COMM_SIZE(mpi_comm_world,npe,ierror)
   call MPI_COMM_RANK(mpi_comm_world,mype,ierror)
!
! NCEP LSF has to use all cores allocated to run this application 
! but this if check can make sure only one core run through the real code.
!

   write(*,*) 'number of cores=',npe

   zero = 0.0
!  Get file names from command line arguements
   num_args = command_argument_count()
   allocate(args(num_args))

   do ix = 1, num_args
       call get_command_argument(ix,args(ix))
   end do
! fv3file: location of UFS grid
! fvcomfile: location of FVCOM data file
! wcstart: warm (restart) or cold start
! inputFVCOMtimes: string of time to use
   fv3file=trim(args(1))
   write(*,*) 'surface file=',trim(fv3file)
   fvcomfile=trim(args(2))    
   write(*,*) 'fvcom file=',trim(fvcomfile)
   wcstart=trim(args(3))
   write(*,*) 'warm or cold start = ', trim(wcstart)
   inputFVCOMselStr=trim(args(4))
!   read(inputFVCOMselStr,*) inputFVCOMsel
   write(*,*) 'select time = ', trim(inputFVCOMselStr)
   clayout_y=trim(args(5))
   read(clayout_y,'(I2)') fv3_io_layout_y
   if (wcstart == 'cold') then
     fv3_io_layout_y=1
   endif
   write(*,*) 'subdomain number=',fv3_io_layout_y
   if(mype < fv3_io_layout_y) then
      write(thisfv3file,'(a,I3.3)') 'stdout_fvcom.',mype
      open(6, file=trim(thisfv3file),form='formatted',status='unknown')
      write(6,*) '===> start process for core = ', mype
   endif
! 
! figure out subdomain grid begin and end grid index in Y for each core  
! The match the subdomains to each core
!
   allocate(fv3_layout_begin(fv3_io_layout_y))
   allocate(fv3_layout_end(fv3_io_layout_y))
   allocate(mpi_layout_begin(npe))
   allocate(mpi_layout_end(npe))
   fv3_layout_begin=0
   fv3_layout_end=0

   workPath='./'
   iy=0
   do ix=1,fv3_io_layout_y
      if(fv3_io_layout_y > 1) then
          write(thisfv3file,'(a,a,a1,I4.4)') trim(workPath), trim(fv3file),'.',ix-1
      else
          write(thisfv3file,'(a,a)') trim(workPath), trim(fv3file)
      endif
      call geo%open(trim(thisfv3file),'r',200)
      call geo%get_dim("yaxis_1",NLAT)
      call geo%close
      fv3_layout_begin(ix)=iy+1
      fv3_layout_end(ix)=iy+nlat
      iy=fv3_layout_end(ix)
   enddo  ! find dimension
   write(6,'(a,20I5)') 'begin index for each subdomain',fv3_layout_begin
   write(6,'(a,20I5)') '  end index for each subdomain',fv3_layout_end

   mypeLocal=mype+1
   if(npe==1) then
      mpi_layout_begin(mypeLocal)=1
      mpi_layout_end(mypeLocal)=fv3_io_layout_y
   else
      mpi_layout_begin=0
      mpi_layout_end=0
      do ix=1,fv3_io_layout_y
         iy=mod(ix-1,npe)+1
         mpi_layout_end(iy)=mpi_layout_end(iy)+1
      enddo
      iy=0
      do ix=1,npe
         if(mpi_layout_end(ix) > 0) then
           mpi_layout_begin(ix)=iy+1  
           mpi_layout_end(ix)=iy+mpi_layout_end(ix)
           iy=mpi_layout_end(ix)
         else
           mpi_layout_begin(ix)=0
           mpi_layout_end(ix)=0
         endif
      enddo 
   endif
   write(6,'(a)') 'begin and end domain index for each core:'
   write(6,'(20I5)') mpi_layout_begin
   write(6,'(20I5)') mpi_layout_end

if(mypeLocal <= fv3_io_layout_y) then

   do ix=mpi_layout_begin(mypeLocal),mpi_layout_end(mypeLocal)

      write(6,*) 'process subdomain ',ix,' with core ',mype
      if(fv3_io_layout_y > 1) then
          write(thisfv3file,'(a,a,a1,I4.4)') trim(workPath), trim(fv3file),'.',ix-1
      else
          write(thisfv3file,'(a,a)') trim(workPath), trim(fv3file)
      endif
      write(6,*) 'sfc data file', trim(thisfv3file)

      call geo%open(trim(thisfv3file),'r',200)
      call geo%get_dim("xaxis_1",NLON)
      call geo%get_dim("yaxis_1",NLAT)
      write(6,*) 'NLON,NLAT:', NLON, NLAT
      call geo%close

      write(6,*) 'Finished reading sfc_data grid information.'
      write(6,*) ' '
!
!  figure out subdomain grid dimension in Y
!
      write(6,*) 'the Y dimensions=',fv3_layout_begin(ix),fv3_layout_end(ix)

!  Allocate variables for I/O

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
!  Read fv3 sfc_data.nc before update

!   fv3file='sfc_data.nc'
!   fv3times: length of time dimension of UFS atmospheric grid (should be 1)
!   t1: index of time dimension to pull (should be 1)
      fv3times=1  
      t1=1

      call fcst%initial('FV3LAM',wcstart)
      call fcst%list_initial
      call fcst%read_n(trim(thisfv3file),'FV3LAM',wcstart,fv3lon,fv3lat,fv3times,t1,&
                   fv3mask,fv3sst,fv3ice,fv3sfcT,fv3iceT,fv3sfcTl,fv3zorl,fv3hice,  &
                   1,nlat)
      call fcst%finish('FV3LAM',wcstart)


      write(6,*) 'fv3times: ', fv3times
      write(6,*) 'time to use: ', t1

!  Read FVCOM input datasets

!   fvcomfile='fvcom.nc'
!   lbctimes: length of time dimension of FVCOM input data (command line input)
! Space infront of ' FVCOM' below is important!!
      call fcst%initial(' FVCOM',wcstart)
      call fcst%list_initial
      call fcst%get_time_ind(trim(fvcomfile),inputFVCOMselStr,indexFVCOMsel)
!   t2: index of time dimension to pull from FVCOM
      t2=indexFVCOMsel
      write(6,*) 'time asked for =', trim(inputFVCOMselStr)
      write(6,*) 'time index selected = ', t2
      call fcst%read_n(trim(fvcomfile),' FVCOM',wcstart,lbclon,lbclat,lbctimes,t2, &
                   lbcmask,lbcsst,lbcice,lbcsfcT,lbciceT,fv3sfcTl,lbczorl,lbchice, &
                   fv3_layout_begin(ix),fv3_layout_end(ix))
      call fcst%finish(' FVCOM',wcstart)

!  Check that the dimensions match

      if (lbclon .ne. nlon .or. lbclat .ne. nlat) then
        write(6,*) 'ERROR: FVCOM/FV3 dimensions do not match:'
        write(6,*) 'lbclon: ', lbclon
        write(6,*) 'nlon: ', nlon
        write(6,*) 'lbclat: ', lbclat
        write(6,*) 'nlat: ', nlat
        stop 'error'
      endif

      write(6,*) 'lbctimes: ', lbctimes
      write(6,*) 'time to use: ', t2

      if (t2 .gt. lbctimes) then
        write(6,*) 'ERROR: Requested time dimension out of range'
        write(6,*) 'Length of time dimension: ', lbctimes
        write(6,*) 'Time index to use: ', t2
        stop 'error'
      endif

!  Update with FVCOM fields and process
!   ice cover data. Ice fraction is set
!   to a minimum of 15% due to FV3-LAM
!   raising any value below 15% to 15%.
   if (wcstart == 'warm') then
     do j=1,nlat
        do i=1,nlon
           if (lbcmask(i,j) > 0. .and. lbcsst(i,j) .ge. -90.0) then !GL Points
              !If ice fraction below 15%, set to 0
              if (lbcice(i,j) < 0.15) then
                lbcice(i,j) = 0.0
                lbchice(i,j) = 0.0 !remove ice thickness
              endif
              fv3ice(i,j) = lbcice(i,j)
              fv3hice(i,j) = lbchice(i,j)

              !If ice in FVCOM, but not in FV3-LAM, change to ice
              if (lbcice(i,j) > 0. .and. fv3mask(i,j) == 0.) then
                fv3mask(i,j) = 2.
                fv3zorl(i,j) = 1.1
              endif
              !If ice in FV3-LAM and not FVCOM, remove it from FV3-LAM
              if (fv3mask(i,j) == 2. .and. lbcice(i,j) == 0.) then
                fv3mask(i,j) = 0.
                fv3zorl(i,j) = zero / zero !Use Fill_Value
              endif
              fv3sst(i,j) = lbcsst(i,j) + 273.15
              fv3sfcT(i,j) = lbcsst(i,j) + 273.15
              fv3iceT(i,j) = lbcsst(i,j) + 273.15
              fv3sfcTl(i,j)= lbcsst(i,j) + 273.15
               !If ice exists in FVCOM, change ice surface temp
              if (lbcice(i,j) > 0.) then
                fv3iceT(i,j) = lbciceT(i,j) + 273.15
              end if
           end if
        enddo
     enddo
   else if (wcstart == 'cold') then
     do j=1,nlat
        do i=1,nlon
           if (lbcmask(i,j) > 0. .and. lbcsst(i,j) .ge. -90.0) then
              !If ice fraction below 15%, set to 0
              if (lbcice(i,j) < 0.15) then
                lbcice(i,j) = 0.0
                lbchice(i,j) = 0.0 !remove ice thickness
              endif
              fv3ice(i,j) = lbcice(i,j)
              fv3hice(i,j) = lbchice(i,j)
              !If ice in FVCOM, but not in FV3-LAM, change to ice
              if (lbcice(i,j) > 0. .and. fv3mask(i,j) == 0.) then
                fv3mask(i,j) = 2.
                fv3zorl(i,j) = 1.1
              endif
              !If ice in FV3-LAM and not FVCOM, remove it from FV3-LAM
              if (fv3mask(i,j) == 2. .and. lbcice(i,j) == 0.) then
                fv3mask(i,j) = 0.
                fv3zorl(i,j) = zero / zero !Use Fill_Value
              endif
              fv3sst(i,j) = lbcsst(i,j) + 273.15
              fv3sfcT(i,j) = lbcsst(i,j) + 273.15
              fv3iceT(i,j) = lbcsst(i,j) + 273.15
               !If ice exists in FVCOM, change ice surface temp
              if (lbcice(i,j) > 0.) then
                fv3iceT(i,j) = lbciceT(i,j) + 273.15
              end if
           end if
        enddo
     enddo
   else
     write(6,*) 'Variable wcstart is not set to either warm or cold'
   end if

! Write out sfc file again

      write(6,*) "udpate file=",trim(thisfv3file)
      call geo%open(trim(thisfv3file),'w',300)
      call geo%replace_var("tsea",NLON,NLAT,fv3sst)
      call geo%replace_var("fice",NLON,NLAT,fv3ice)
      call geo%replace_var("slmsk",NLON,NLAT,fv3mask)
      call geo%replace_var("tisfc",NLON,NLAT,fv3iceT)
      call geo%replace_var("hice",NLON,NLAT,fv3hice)
      if (wcstart == 'cold') then
! Add_New_Var takes names of (Variable,Dim1,Dim2,Dim3,Long_Name,Units)
         call geo%replace_var("zorl",NLON,NLAT,fv3zorl)
         call geo%add_new_var('glmsk','xaxis_1','yaxis_1','Time','glmsk','none')
         call geo%replace_var('glmsk',NLON,NLAT,lbcmask)
      end if
      if (wcstart == 'warm') then
         call geo%replace_var("zorli",NLON,NLAT,fv3zorl)
         call geo%replace_var("tsfc",NLON,NLAT,fv3sfcT)
         call geo%replace_var("tsfcl",NLON,NLAT,fv3sfcTl)
         call geo%add_new_var('glmsk','xaxis_1','yaxis_1','Time','glmsk','none')
         call geo%replace_var('glmsk',NLON,NLAT,lbcmask)
      end if
      call geo%close

      deallocate(fv3ice)
      deallocate(fv3sfcT)
      deallocate(fv3sst)
      deallocate(fv3mask)
      deallocate(fv3iceT)
      deallocate(fv3sfcTl)
      deallocate(fv3zorl)
      deallocate(fv3hice)

      deallocate(lbcice)
      deallocate(lbcsfcT)
      deallocate(lbcsst)
      deallocate(lbcmask)
      deallocate(lbciceT)
      deallocate(lbczorl)
      deallocate(lbchice)

      write(6,*) "=== LOWBC UPDATE SUCCESS ===",ix

   enddo ! loop through subdomain

endif ! mype==0

call MPI_FINALIZE(ierror)


end program process_FVCOM 
