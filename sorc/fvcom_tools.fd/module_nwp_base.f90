!> @file
!! @brief Defines nwp observation data structure. 
!! @author David Wright, University of Michigan and GLERL @date 17 Aug 2020

!> This module defines nwp observation data structure and the method
!! to read and write observations from and to those data
!! structures. It is used by ingest_FVCOM.f90.
!!
!! This script is strongly based upon Eric James' (ESRL/GSL) work with
!! HRRR/WRF.
!!
!! @author David Wright, University of Michigan and GLERL @date 17 Aug 2020
!!
module module_nwp_base

   use kinds, only: r_kind, r_single, rmissing

   implicit none

   public :: nwpbase
   public :: nwplocation

   private

!  Define a nwp observation type.

   type nwplocation
      real(r_single) :: lon  !< stroke longitude
      real(r_single) :: lat  !< stroke latitiude
   end type nwplocation

!  Define a nwp observation type to contain actual data.

   type, extends(nwplocation) :: nwpbase
!  HOW DOES THIS POINTER THING WORK?
      type(nwpbase), pointer :: next => NULL() !< Pointer.
      real(r_single) :: time                 !< observation time.
      integer :: numvar                      !< number of variables in this obs type.
!      real(r_single), allocatable :: obs(:)  !< observation value (# numvar).
      real(r_kind), allocatable :: obs(:) !< Observations.
      logical :: ifquality                   !< do these obs include quality info? GLM has flash_quality_flag.
      integer, allocatable :: quality(:)     !< if so, quality flags.
      contains
         procedure :: list => list_obsbase !< List contents of obs. @return
         procedure :: alloc => alloc_obsbase !< Allocate memory for observations. @return
         procedure :: destroy => destroy_obsbase  !< Release memory. @return
   end type nwpbase

   contains

     !> This subroutine lists the contents of a base nwp observation.
     !!
     !! @param this the base nwp obervation
     !! @author David Wright, University of Michigan and GLERL
     !! @date 17 Aug 2020
      subroutine list_obsbase(this)
         class(nwpbase) :: this

         integer :: i, numvar

!        Write out the lon, lat, and time of the ob

         write(*,'(a,3f10.3)') 'LIGHTNING OB: longitude, latitude, time =', &
            this%lon, this%lat, this%time

!        Loop through all variables and print out obs and quality

         numvar = this%numvar
         if (numvar >= 1) then
! MULTI-DIMENSIONAL EXAMPLE IN  module_obs_base.f90
            write(*,'(a4,10F12.2)') 'obs=', (this%obs(i),i=1,numvar)
            if(this%ifquality) &
            write(*,'(a4,10I12)') 'qul=', (this%quality(i),i=1,numvar)
         else
            write(*,*) 'No obs for this location'
         endif

      end subroutine list_obsbase

      !> This subroutine allocates memory for base nwp observation
      !! variables.
      !!
      !! @param this the base nwp obervation
      !! @param[in] numvar number of variables in this ob type
      !! @param[in] ifquality does this observation include quality
      !! information?
      !! @author David Wright, University of Michigan and GLERL @date 17 Aug 2020
      subroutine alloc_obsbase(this,numvar,ifquality)


         class(nwpbase) :: this

         integer, intent(in) :: numvar

         logical, intent(in), optional :: ifquality

         if (numvar >= 1) then
            this%numvar = numvar

            if(allocated(this%obs)) deallocate(this%obs)
            allocate(this%obs(numvar))

            this%ifquality=.false.
            if(present(ifquality)) this%ifquality = ifquality
            if(this%ifquality) allocate(this%quality(numvar))

         else
            write(*,*) 'alloc_obsbase Error: dimension must be larger than 0:', numvar
            stop 1234
         endif

      end subroutine alloc_obsbase

      !> This subroutine releases memory associated with nwp
      !! observations.
      !!
      !! @param this the base nwp obervation
      !! @author David Wright, University of Michigan and GLERL @date 17 Aug 2020
      subroutine destroy_obsbase(this)

         class(nwpbase) :: this

         this%numvar = 0
         this%time = 0

         if(allocated(this%obs)) deallocate(this%obs)

         this%ifquality=.false.
         if(allocated(this%quality)) deallocate(this%quality)

         this%next => NULL()

      end subroutine destroy_obsbase

end module module_nwp_base
