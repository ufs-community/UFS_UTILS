!> @file
!! @brief Defines FV3LAM and FVCOM forecast data structure.
!! @author David Wright, University of Michigan and GLERL, @date 17 Aug 2020

!> This module defines FV3LAM and FVCOM forecast data structure and
!! the method to read and write observations from and to those data
!! structures. It is used by ingest_FVCOM.f90.
!!
!! This script is strongly based upon Eric James' (ESRL/GSL) work with
!! HRRR/WRF to get FVCOM data into the model.
!!
!! @author David Wright, University of Michigan and GLERL,
!! @date 17 Aug 2020
!!
module module_nwp

   use kinds, only: r_kind, r_single, i_short, rmissing
   use module_nwp_base, only: nwpbase
!   use module_map_utils, only: map_util
   use module_ncio, only: ncio

   implicit none

   public :: fcst_nwp
   public :: nwp_type

   private
   type :: nwp_type
      character(len=6) :: datatype !< Data type.
      integer :: numvar !< Number of variabls.
      integer :: xlat !< Number of latitudes.
      integer :: xlon !< Number of longitudes.
      integer :: xtime !< Number of times.
      integer :: i_mask !< Is var visible (always 1).
      integer :: i_sst !< Index of sst var.
      integer :: i_ice !< Index of ice var.
      integer :: i_sfcT !< Index of sst temp var.
      integer :: i_iceT !< Index of ice temp var.
      character(len=20), allocatable :: varnames(:) !< Variable names.
      character(len=20), allocatable :: latname !< Latitude name.
      character(len=20), allocatable :: lonname !< Longitude name.
      character(len=20), allocatable :: dimnameEW !< East/West dimension name.
      character(len=20), allocatable :: dimnameNS !< North/South dimension name.
      character(len=20), allocatable :: dimnameTIME !< Time dimension name.

      real(r_kind), allocatable :: nwp_mask(:,:,:) !< Land/water mask 3D array
      real(r_kind), allocatable :: nwp_sst(:,:,:) !< SST 3D array
      real(r_kind), allocatable :: nwp_ice(:,:,:) !< Over water ice concentration 3D array
      real(r_kind), allocatable :: nwp_sfcT(:,:,:) !< Skin temperature 3D array
      real(r_kind), allocatable :: nwp_iceT(:,:,:)  !< Ice skin temperature 3D array
   end type nwp_type

   type, extends(nwp_type) :: fcst_nwp
      ! The pointers are carryover from when I inherited the code from
      ! GSL's work with HRRR for a similar use. I am not sure with
      ! object based coding in Fortran if it needs to have parts
      ! initialized to gain access to the procedures within it. - D. Wright.
      type(nwpbase), pointer :: head => NULL() !< Pointer to head of list.
      type(nwpbase), pointer :: tail => NULL() !< Pointer to tail of list.
      contains
         procedure :: initial => initial_nwp !< Defines vars and names. @return
         procedure :: list_initial => list_initial_nwp !< List the setup. @return
         procedure :: read_n => read_nwp !< Initialize arrays, get data. @return
         procedure :: finish => finish_nwp !< Finish and deallocate. @return
   end type fcst_nwp

   type(ncio) :: ncdata !< Wrapper object for netCDF data file.
!   type(map_util) :: map

   contains

     !> This subroutine defines the number of variables and their
     !! names for each NWP data type. The indices of the variables are
     !! also defined for later reference.
     !!
     !! @param this fcst_nwp object
     !! @param[in] itype either ' FVCOM' or 'FV3LAM'.
     !! @author David Wright, University of Michigan and GLERL
      subroutine initial_nwp(this,itype)
         class(fcst_nwp) :: this

         character(len=6), intent(in) :: itype

!        FVCOM grid
         if (itype==' FVCOM') then
            this%datatype = itype
            this%numvar = 4

            this%i_mask = 1
            this%i_sst = 2
            this%i_ice = 3
            this%i_iceT = 4
            this%i_sfcT = 0

            allocate(this%varnames(this%numvar))
            this%varnames(1) = 'glmask'
            this%varnames(2) = 'tsfc'
            this%varnames(3) = 'aice'
            this%varnames(4) = 'tisfc'

            allocate(this%latname)
            allocate(this%lonname)
            this%latname = 'lat'
            this%lonname = 'lon'

            allocate(this%dimnameEW)
            allocate(this%dimnameNS)
            allocate(this%dimnameTIME)
            this%dimnameEW = 'lon'
            this%dimnameNS = 'lat'
            this%dimnameTIME = 'Time'

!        FV3LAM grid

         else if (trim(itype)=='FV3LAM') then
            this%datatype = itype
            this%numvar = 4

            this%i_mask = 1
            this%i_sst = 2
            this%i_ice = 3
            this%i_iceT = 4
            this%i_sfcT = 0

            allocate(this%varnames(this%numvar))
            this%varnames(1) = 'slmsk'
            this%varnames(2) = 'tsea'
            this%varnames(3) = 'fice'
            this%varnames(4) = 'tisfc'

            allocate(this%latname)
            allocate(this%lonname)
            this%latname = 'yaxis_1'
            this%lonname = 'xaxis_1'

            allocate(this%dimnameEW)
            allocate(this%dimnameNS)
            allocate(this%dimnameTIME)
            this%dimnameEW = 'xaxis_1'
            this%dimnameNS = 'yaxis_1'
            this%dimnameTIME = 'Time'

!        If the data type does not match one of the known types, exit.

         else
            write(*,*) 'Unknown data type:', itype
            stop 1234
         end if

         this%head => NULL()
         this%tail => NULL()

         write(*,*) 'Finished initial_nwp'
         write(*,*) ' '

      end subroutine initial_nwp

      !> This subroutine lists the setup for NWP data that was done by
      !! the initial_nwp subroutine.
      !!
      !! @param this fcst_nwp object
      !! @author David Wright, University of Michigan and GLERL
      subroutine list_initial_nwp(this)

         class(fcst_nwp) :: this

         integer :: k

         write(*,*) 'List initial setup for ', this%datatype
         write(*,*) 'number of variables ', this%numvar
         write(*,*) 'variable index: mask, sst, ice, sfcT'
         write(*,'(15x,10I3)') this%i_mask, this%i_sst, this%i_ice, &
      &      this%i_sfcT
         write(*,*) 'variable name:'
         do k=1,this%numvar
            write(*,*) k,trim(this%varnames(k))
         enddo

         write(*,*) 'Finished list_initial_nwp'
         write(*,*) ' '

      end subroutine list_initial_nwp

      !> This subroutine initializes arrays to receive the NWP data,
      !! and opens the file and gets the data.
      !!
      !! @param this fcst_nwp ojbect
      !! @param[in] filename netcdf file name
      !! @param[in] itype either ' FVCOM' or 'FV3LAM'
      !! @param[inout] numlon number of grid points in x-direction
      !! @param[inout] numlat number of grid poinst in y-direction
      !! @param[inout] numtimes length of time dimension
      !! @param[in] time_to_get integer of time dimension to read in
      !! @param[inout] mask Water points mask
      !! @param[inout] sst Water surface temperature
      !! @param[inout] ice Ice concentration (%)
      !! @param[inout] sfcT Skin Temperature
      !! @param[inout] iceT Ice Skin Temperature
      !!
      !! @author David Wright, University of Michigan and GLERL
      subroutine read_nwp(this,filename,itype,numlon,numlat,numtimes,time_to_get,mask,sst,ice,sfcT,iceT)

         class(fcst_nwp) :: this

         character(len=5), intent(in) :: itype
         character(len=*), intent(in) :: filename

         integer, intent(in) :: time_to_get
         integer, intent(inout) :: numlon, numlat, numtimes
!         real(r_single), intent(inout) :: mask(:,:), sst(:,:), ice(:,:), sfcT(:,:)
         real(r_kind), intent(inout) :: mask(:,:),sst(:,:),ice(:,:),sfcT(:,:) &
                                        ,iceT(:,:)

!        Open the file using module_ncio.f90 code, and find the number of
!        lat/lon points

         call ncdata%open(trim(filename),'r',200)
         call ncdata%get_dim(this%dimnameEW,this%xlon)
         call ncdata%get_dim(this%dimnameNS,this%xlat)
         call ncdata%get_dim(this%dimnameTIME,this%xtime)

         write(*,*) 'number of longitudes for file ', filename, this%xlon
         numlon = this%xlon
         write(*,*) 'number of latitudes for file ', filename, this%xlat
         numlat = this%xlat
         write(*,*) 'number of times for file ', filename, this%xtime
         numtimes = this%xtime

!        Allocate all the arrays to receive data

         allocate(this%nwp_mask(this%xlon,this%xlat,this%xtime))
         allocate(this%nwp_sst(this%xlon,this%xlat,this%xtime))
         allocate(this%nwp_ice(this%xlon,this%xlat,this%xtime))
         allocate(this%nwp_sfcT(this%xlon,this%xlat,this%xtime))
         allocate(this%nwp_iceT(this%xlon,this%xlat,this%xtime))

!        Get variables from the data file, but only if the variable is
!        defined for that data type.

         if (this%i_mask .gt. 0) then
            call ncdata%get_var(this%varnames(this%i_mask),this%xlon,  &
                                this%xlat,this%xtime,this%nwp_mask)
            mask = this%nwp_mask(:,:,1)
         end if
         if (this%i_sst .gt. 0) then
            call ncdata%get_var(this%varnames(this%i_sst),this%xlon,  &
                                this%xlat,this%xtime,this%nwp_sst)
            sst = this%nwp_sst(:,:,time_to_get)
         end if
         if (this%i_ice .gt. 0) then
            call ncdata%get_var(this%varnames(this%i_ice),this%xlon,  &
                                this%xlat,this%xtime,this%nwp_ice)
            ice = this%nwp_ice(:,:,time_to_get)
         end if
         if (this%i_sfcT .gt. 0) then
            call ncdata%get_var(this%varnames(this%i_sfcT),this%xlon,  &
                                this%xlat,this%xtime,this%nwp_sfcT)
            sfcT = this%nwp_sfcT(:,:,time_to_get)
         end if
         if (this%i_iceT .gt. 0) then
             call ncdata%get_var(this%varnames(this%i_iceT),this%xlon,  &
                                this%xlat,this%xtime,this%nwp_iceT)
            iceT = this%nwp_iceT(:,:,time_to_get)
         end if 

!        Close the netCDF file.

         call ncdata%close

         write(*,*) 'Finished read_nwp'
         write(*,*) ' '

      end subroutine read_nwp

      !> Finish and deallocate.
      !!
      !! @param this fcst_nwp object
      !! @author David Wright, University of Michigan and GLERL
      subroutine finish_nwp(this)

         class(fcst_nwp) :: this

         type(nwpbase), pointer :: thisobs,thisobsnext

         deallocate(this%varnames)
         deallocate(this%latname)
         deallocate(this%lonname)
         deallocate(this%dimnameEW)
         deallocate(this%dimnameNS)
         deallocate(this%dimnameTIME)
         deallocate(this%nwp_mask)
         deallocate(this%nwp_sst)
         deallocate(this%nwp_ice)
         deallocate(this%nwp_sfcT)
         deallocate(this%nwp_iceT)

         thisobs => this%head
         if(.NOT.associated(thisobs)) then
            write(*,*) 'No memory to release'
            return
         endif
         do while(associated(thisobs))
!            write(*,*) 'destroy ==',thisobs%name

            thisobsnext => thisobs%next
            call thisobs%destroy()
            thisobs => thisobsnext
         enddo

         write(*,*) 'Finished finish_nwp'
         write(*,*) ' '

      end subroutine finish_nwp

end module module_nwp
