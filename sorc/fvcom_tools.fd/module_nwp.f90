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
!   use module_map_utils, only: map_util
   use module_ncio, only: ncio

   implicit none

   public :: fcst_nwp

   private
   type :: fcst_nwp
      character(len=6) :: datatype !< Data type.
      integer :: numvar !< Number of variabls.
      integer :: xlat !< Number of latitudes.
      integer :: xlon !< Number of longitudes.
      integer :: xtime !< Number of times.
      integer :: datelen !< Length of date string.
      integer :: i_mask !< Is var visible (always 1).
      integer :: i_sst !< Index of sst var.
      integer :: i_ice !< Index of ice var.
      integer :: i_sfcT !< Index of sst temp var.
      integer :: i_iceT !< Index of ice temp var.
      integer :: i_sfcTl !< Index of sfcTl
      integer :: i_zorl !< Index of surface roughness
      integer :: i_hice !< Index of ice thickness
      character(len=20), allocatable :: varnames(:) !< Variable names.
      character(len=20), allocatable :: latname !< Latitude name.
      character(len=20), allocatable :: lonname !< Longitude name.
      character(len=20), allocatable :: dimnameEW !< East/West dimension name.
      character(len=20), allocatable :: dimnameNS !< North/South dimension name.
      character(len=20), allocatable :: dimnameTIME !< Time dimension name.
      character(len=20), allocatable :: dimnameDATE !< String dimension name.
      character(len=1), allocatable :: times(:,:) !< Array of times in FVCOM.

      real(r_kind), allocatable :: nwp_mask_c(:,:) !< cold start land/water mask 3d array
      real(r_kind), allocatable :: nwp_sst_c(:,:,:) !< cold start sst 3d array
      real(r_kind), allocatable :: nwp_ice_c(:,:,:) !< cold start over water ice concentration 3d array
      real(r_kind), allocatable :: nwp_sfct_c(:,:,:) !< cold start skin temperature 3d array
      real(r_kind), allocatable :: nwp_icet_c(:,:,:)  !< cold start ice skin temperature 3d array
      real(r_kind), allocatable :: nwp_zorl_c(:,:,:) !< cold start surface roughness
      real(r_kind), allocatable :: nwp_hice_c(:,:,:) !< cold start ice thickness

      real(r_kind), allocatable :: nwp_mask_w(:,:) !< warm start land/water mask 3d array
      real(r_kind), allocatable :: nwp_sst_w(:,:) !< warm start sst 3d array
      real(r_kind), allocatable :: nwp_ice_w(:,:) !< warm start over water ice concentration 3d array
      real(r_kind), allocatable :: nwp_sfct_w(:,:) !< warm start skin temperature 3d array
      real(r_kind), allocatable :: nwp_icet_w(:,:)  !< warm start ice skin temperature 3d array
      real(r_kind), allocatable :: nwp_sfctl_w(:,:) !< warm start skin temperature 3d array
      real(r_kind), allocatable :: nwp_zorl_w(:,:) !< warm start surface roughness
      real(r_kind), allocatable :: nwp_hice_w(:,:) !< warm start ice thickness

      contains
         procedure :: initial => initial_nwp !< Defines vars and names. @return
         procedure :: list_initial => list_initial_nwp !< List the setup. @return
         procedure :: read_n => read_nwp !< Initialize arrays, get data. @return
         procedure :: get_time_ind => get_time_ind_nwp !< Get time ind. @return
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
     !! @param[in] wcstart either 'warm' or 'cold'.
     !! @author David Wright, University of Michigan and GLERL
      subroutine initial_nwp(this,itype,wcstart)
         class(fcst_nwp) :: this

         character(len=6), intent(in) :: itype
         character(len=4), intent(in) :: wcstart

!        FVCOM grid
         if (itype==' FVCOM') then
            this%datatype = itype
            this%numvar = 5

            this%i_mask = 1
            this%i_sst = 2
            this%i_ice = 3
            this%i_iceT = 4
            this%i_hice = 5
            this%i_sfcT = 0
            this%i_zorl = 0

            allocate(this%varnames(this%numvar))
            this%varnames(1) = 'glmask'
            this%varnames(2) = 'tsfc'
            this%varnames(3) = 'aice'
            this%varnames(4) = 'tisfc'
            this%varnames(5) = 'vice'

            allocate(this%latname)
            allocate(this%lonname)
            this%latname = 'lat'
            this%lonname = 'lon'

            allocate(this%dimnameEW)
            allocate(this%dimnameNS)
            allocate(this%dimnameTIME)
            allocate(this%dimnameDATE)
            this%dimnameEW = 'lon'
            this%dimnameNS = 'lat'
            this%dimnameTIME = 'Time'
            this%dimnameDATE = 'DateStrLen'

!        FV3LAM grid

         else if (trim(itype)=='FV3LAM' .AND. wcstart=='warm') then
            this%datatype = itype
            this%numvar = 8

            this%i_mask = 1
            this%i_sst = 2
            this%i_ice = 3
            this%i_iceT = 4
            this%i_sfcT = 5
            this%i_sfcTl= 6
            this%i_zorl = 7
            this%i_hice = 8

            allocate(this%varnames(this%numvar))
            this%varnames(1) = 'slmsk'
            this%varnames(2) = 'tsea'
            this%varnames(3) = 'fice'
            this%varnames(4) = 'tisfc'
            this%varnames(5) = 'tsfc'
            this%varnames(6) = 'tsfcl'
            this%varnames(7) = 'zorli'
            this%varnames(8) = 'hice'

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

         else if (trim(itype)=='FV3LAM' .AND. wcstart=='cold') then
            this%datatype = itype
            this%numvar = 6

            this%i_mask = 1
            this%i_sst = 2
            this%i_ice = 3
            this%i_iceT = 4
            this%i_zorl = 5
            this%i_hice = 6
            this%i_sfcT = 0

            allocate(this%varnames(this%numvar))
            this%varnames(1) = 'slmsk'
            this%varnames(2) = 'tsea'
            this%varnames(3) = 'fice'
            this%varnames(4) = 'tisfc'
            this%varnames(5) = 'zorl'
            this%varnames(6) = 'hice'

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
            write(6,*) 'Unknown data type:', itype
            stop 1234
         end if

         write(6,*) 'Finished initial_nwp'
         write(6,*) ' '

      end subroutine initial_nwp

      !> This subroutine lists the setup for NWP data that was done by
      !! the initial_nwp subroutine.
      !!
      !! @param this fcst_nwp object
      !! @author David Wright, University of Michigan and GLERL
      subroutine list_initial_nwp(this)

         class(fcst_nwp) :: this

         integer :: k

         write(6,*) 'List initial setup for ', this%datatype
         write(6,*) 'number of variables ', this%numvar
         write(6,*) 'variable index: mask, sst, ice, sfcT, sfcTl'
         write(6,'(15x,10I3)') this%i_mask, this%i_sst, this%i_ice, &
      &      this%i_sfcT, this%i_sfcTl
         write(6,*) 'variable name:'
         do k=1,this%numvar
            write(6,*) k,trim(this%varnames(k))
         enddo

         write(6,*) 'Finished list_initial_nwp'
         write(6,*) ' '

      end subroutine list_initial_nwp

      !> This subroutine initializes arrays to receive the NWP data,
      !! and opens the file and gets the data.
      !!
      !! @param this fcst_nwp ojbect
      !! @param[in] filename netcdf file name
      !! @param[in] itype either ' FVCOM' or 'FV3LAM'
      !! @param[in] wcstart either 'warm' or 'cold'.
      !! @param[inout] numlon number of grid points in x-direction
      !! @param[inout] numlat number of grid poinst in y-direction
      !! @param[inout] numtimes length of time dimension
      !! @param[in] time_to_get integer of time dimension to read in
      !! @param[inout] mask Water points mask
      !! @param[inout] sst Water surface temperature
      !! @param[inout] ice Ice concentration (%)
      !! @param[inout] sfcT Skin Temperature
      !! @param[inout] iceT Ice Skin Temperature
      !! @param[inout] sfcTl Skin Temperature in restart file
      !! @param[inout] zorl Surface roughness length
      !! @param[inout] hice Ice thickness
      !! @param[in]    ybegin Start grid point in Y direction for the domain
      !! @param[in]    yend   End grid point in Y direction for the domain
      !!
      !! @author David Wright, University of Michigan and GLERL
      subroutine read_nwp(this,filename,itype,wcstart,numlon,numlat,numtimes,time_to_get,mask,sst,ice,sfcT,iceT,sfcTl,zorl,hice,ybegin,yend)

         class(fcst_nwp) :: this

         character(len=6), intent(in) :: itype
         character(len=*), intent(in) :: filename
         character(len=4), intent(in) :: wcstart

         integer, intent(in) :: time_to_get
         integer, intent(in) :: ybegin,yend
         integer, intent(inout) :: numlon, numlat, numtimes
!         real(r_single), intent(inout) :: mask(:,:), sst(:,:), ice(:,:), sfcT(:,:)
         real(r_kind), intent(inout) :: mask(:,:),sst(:,:),ice(:,:),sfcT(:,:) &
                                        ,iceT(:,:),sfcTl(:,:),zorl(:,:),hice(:,:)

!
!        Open the file using module_ncio.f90 code, and find the number of
!        lat/lon points

         call ncdata%open(trim(filename),'r',200)
         call ncdata%get_dim(this%dimnameEW,this%xlon)
         call ncdata%get_dim(this%dimnameNS,this%xlat)
         call ncdata%get_dim(this%dimnameTIME,this%xtime)

         write(6,*) 'number of longitudes for file ', filename, this%xlon
         numlon = this%xlon
         write(6,*) 'number of latitudes for file ', filename, this%xlat
         !numlat = this%xlat
         numlat = yend-ybegin+1
         write(6,*) 'number of times for file ', filename, this%xtime
         numtimes = this%xtime
         write(6,*) 'the range of Y for this domain is=',ybegin,yend

!        Allocate all the arrays to receive data
         if (wcstart == 'cold' .OR. itype == ' FVCOM') then
            allocate(this%nwp_mask_c(this%xlon,this%xlat))
            allocate(this%nwp_sst_c(this%xlon,this%xlat,this%xtime))
            allocate(this%nwp_ice_c(this%xlon,this%xlat,this%xtime))
            allocate(this%nwp_sfcT_c(this%xlon,this%xlat,this%xtime))
            allocate(this%nwp_iceT_c(this%xlon,this%xlat,this%xtime))
            allocate(this%nwp_zorl_c(this%xlon,this%xlat,this%xtime))
            allocate(this%nwp_hice_c(this%xlon,this%xlat,this%xtime))

!        Get variables from the data file, but only if the variable is
!        defined for that data type.

            write(6,*) 'itype = ', itype
            write(6,*) 'wcstart = ', wcstart
            write(6,*) 'xlat = ', this%xlat
            write(6,*) 'xlon = ', this%xlon
            write(6,*) 'xtime = ', this%xtime

            if (this%i_mask .gt. 0) then
               call ncdata%get_var(this%varnames(this%i_mask),this%xlon,  &
                                   this%xlat,this%nwp_mask_c)
               mask = this%nwp_mask_c(:,ybegin:yend)
            end if
            if (this%i_sst .gt. 0) then
               write(6,*) 'get sst for cold or FVCOM'
               call ncdata%get_var(this%varnames(this%i_sst),this%xlon,  &
                                   this%xlat,this%xtime,this%nwp_sst_c)
               sst = this%nwp_sst_c(:,ybegin:yend,time_to_get)
            end if
            if (this%i_ice .gt. 0) then
               call ncdata%get_var(this%varnames(this%i_ice),this%xlon,  &
                                   this%xlat,this%xtime,this%nwp_ice_c)
               ice = this%nwp_ice_c(:,ybegin:yend,time_to_get)
            end if
            if (this%i_sfcT .gt. 0) then
               call ncdata%get_var(this%varnames(this%i_sfcT),this%xlon,  &
                                   this%xlat,this%xtime,this%nwp_sfcT_c)
               sfcT = this%nwp_sfcT_c(:,ybegin:yend,time_to_get)
            end if
            if (this%i_iceT .gt. 0) then
                call ncdata%get_var(this%varnames(this%i_iceT),this%xlon,  &
                                    this%xlat,this%xtime,this%nwp_iceT_c)
                iceT = this%nwp_iceT_c(:,ybegin:yend,time_to_get)
            end if
            if (this%i_zorl .gt. 0) then
                call ncdata%get_var(this%varnames(this%i_zorl),this%xlon,  &
                                    this%xlat,this%xtime,this%nwp_zorl_c)
                zorl = this%nwp_zorl_c(:,ybegin:yend,time_to_get)
            end if 
            if (this%i_hice .gt. 0) then
                call ncdata%get_var(this%varnames(this%i_hice),this%xlon,  &
                                    this%xlat,this%xtime,this%nwp_hice_c)
                hice = this%nwp_hice_c(:,ybegin:yend,time_to_get)
            end if 
 
         else if (wcstart == 'warm') then
            allocate(this%nwp_mask_w(this%xlon,this%xlat))
            allocate(this%nwp_sst_w(this%xlon,this%xlat))
            allocate(this%nwp_ice_w(this%xlon,this%xlat))
            allocate(this%nwp_sfcT_w(this%xlon,this%xlat))
            allocate(this%nwp_iceT_w(this%xlon,this%xlat))
            allocate(this%nwp_sfcTl_w(this%xlon,this%xlat))
            allocate(this%nwp_zorl_w(this%xlon,this%xlat))
            allocate(this%nwp_hice_w(this%xlon,this%xlat))
!        Get variables from the data file, but only if the variable is
!        defined for that data type.
      
            write(6,*) 'itype = ', itype
            write(6,*) 'wcstart =', wcstart
            write(6,*) 'xlat = ', this%xlat
            write(6,*) 'xlon = ', this%xlon
            write(6,*) 'xtime = ', this%xtime

            if (this%i_mask .gt. 0) then
               call ncdata%get_var(this%varnames(this%i_mask),this%xlon,  &
                                   this%xlat,this%nwp_mask_w)
               mask = this%nwp_mask_w(:,ybegin:yend)
            end if
            if (this%i_sst .gt. 0) then
               call ncdata%get_var(this%varnames(this%i_sst),this%xlon,  &
                                   this%xlat,this%nwp_sst_w)
               sst = this%nwp_sst_w(:,ybegin:yend)
            end if
            if (this%i_ice .gt. 0) then
               call ncdata%get_var(this%varnames(this%i_ice),this%xlon,  &
                                   this%xlat,this%nwp_ice_w)
               ice = this%nwp_ice_w(:,ybegin:yend)
            end if
            if (this%i_sfcT .gt. 0) then
               call ncdata%get_var(this%varnames(this%i_sfcT),this%xlon,  &
                                   this%xlat,this%nwp_sfcT_w)
               sfcT = this%nwp_sfcT_w(:,ybegin:yend)
            end if
            if (this%i_iceT .gt. 0) then
                call ncdata%get_var(this%varnames(this%i_iceT),this%xlon,  &
                                    this%xlat,this%nwp_iceT_w)
                iceT = this%nwp_iceT_w(:,ybegin:yend)
            end if 
            if (this%i_sfcTl .gt. 0) then
               call ncdata%get_var(this%varnames(this%i_sfcTl),this%xlon,  &
                                   this%xlat,this%nwp_sfcTl_w)
               sfcTl = this%nwp_sfcTl_w(:,ybegin:yend)
            end if
            if (this%i_zorl .gt. 0) then
                call ncdata%get_var(this%varnames(this%i_zorl),this%xlon,  &
                                    this%xlat,this%nwp_zorl_w)
                zorl = this%nwp_zorl_w(:,ybegin:yend)
            end if 
            if (this%i_hice .gt. 0) then
                call ncdata%get_var(this%varnames(this%i_hice),this%xlon,  &
                                    this%xlat,this%nwp_hice_w)
                hice = this%nwp_hice_w(:,ybegin:yend)
            end if 

         else
            write(6,*) 'Choose either "warm" or "cold" for file'
            stop 'Error in wcstart. Check spelling or if variable was assigned'
         end if
!        Close the netCDF file.

         call ncdata%close

         write(6,*) 'Finished read_nwp'
         write(6,*) ' '

      end subroutine read_nwp

      !> Finish and deallocate.
      !!
      !! @param this fcst_nwp object
      !! @param[in] itype either ' FVCOM' or 'FV3LAM'
      !! @param[in] wcstart either 'warm' or 'cold'.
      !! @author David Wright, University of Michigan and GLERL
      subroutine finish_nwp(this,itype,wcstart)

         class(fcst_nwp) :: this
         character(len=6), intent(in) :: itype
         character(len=4), intent(in) :: wcstart

         deallocate(this%varnames)
         deallocate(this%latname)
         deallocate(this%lonname)
         deallocate(this%dimnameEW)
         deallocate(this%dimnameNS)
         deallocate(this%dimnameTIME)
         if (wcstart == 'cold' .OR. itype==' FVCOM') then
            deallocate(this%nwp_mask_c)
            deallocate(this%nwp_sst_c)
            deallocate(this%nwp_ice_c)
            deallocate(this%nwp_sfcT_c)
            deallocate(this%nwp_iceT_c)
            deallocate(this%nwp_zorl_c)
            deallocate(this%nwp_hice_c)
            if (itype==' FVCOM') deallocate(this%dimnameDATE)
         else if (wcstart == 'warm') then
            deallocate(this%nwp_mask_w)
            deallocate(this%nwp_sst_w)
            deallocate(this%nwp_ice_w)
            deallocate(this%nwp_sfcT_w)
            deallocate(this%nwp_iceT_w)
            deallocate(this%nwp_sfcTl_w)
            deallocate(this%nwp_zorl_w)
            deallocate(this%nwp_hice_w)
         else
            write(6,*) 'no deallocation'
         end if

         write(6,*) 'Finished finish_nwp'
         write(6,*) ' '

      end subroutine finish_nwp

      !> This subroutine searches the FVCOM 'Times' variable
      !!  and returns the matching index
      !!
      !! @param this fcst_nwp ojbect
      !! @param[in] filename netcdf file name
      !! @param[in] instr string of requested time
      !! @param[out] outindex int index that matches instr
      !!
      !! @author David Wright, University of Michigan and GLERL
      subroutine get_time_ind_nwp(this,filename,instr,outindex)

         class(fcst_nwp) :: this

         character(len=*), intent(in) :: filename
         character(len=*), intent(in) :: instr
         integer, intent(out) :: outindex

         character(len=26) :: temp
         integer :: foundind
         integer :: k,i
 
!        Open the file using module_ncio.f90 code, and find the length of
!        time in the file
         call ncdata%open(trim(filename),'r',200)
         call ncdata%get_dim(this%dimnameTIME,this%xtime)
         call ncdata%get_dim(this%dimnameDATE,this%datelen)
         write(6,*) 'xtime = ', this%xtime
         write(6,*) 'datelen = ', this%datelen
         allocate(this%times(this%datelen,this%xtime))
         call ncdata%get_var('Times',this%datelen,this%xtime,this%times)

         foundind = 0

         do k=1,this%xtime,1
            do i = 1,len(temp),1
               temp(i:i) = this%times(i,k)
            end do
            if (trim(temp) == trim(instr)) then !If times are equal return k
              outindex = k
              foundind = 1
            end if
         end do
         if (foundind == 0) then
            outindex = -999
            deallocate(this%times)
            call ncdata%close
            write(6,*) 'WARNING: Supplied time not found in file: ', trim(instr)
            write(6,*) 'Stoppping fvcom_to_FV3 and proceeding without using FVCOM data'
            stop
         end if

         deallocate(this%times)
         call ncdata%close

      end subroutine get_time_ind_nwp

end module module_nwp
