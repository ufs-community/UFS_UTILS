!> @file
!! @brief Routines for applyng land DA increments
!! @author Clara Draper ESRL/PSL

module land_increments

    private

    public add_increment_soil
    public add_increment_snow
    public calculate_landinc_mask
    public apply_land_da_adjustments_stc
    public apply_land_da_adjustments_snd

    integer, parameter            :: lsm_noah=1      !< flag for NOAH land surface model
                                                     !! copied from GFS_typedefs.F90
contains

 !> Read in gsi file with soil state  increments (on the gaussian
 !! grid), interpolate increments to the cubed-sphere tile, and
 !! add to the soil states. Adapted from adjust_nsst.
 !! Currently only coded for soil temperature. Soil moisture will
 !! need the model soil moisture paramaters for regridding.
 !!
 !!  Does not make a temperature update if snow differ
 !!  between fg and anal (allow correction of snow to
 !!  address temperature error first), or if snow is present 
 !!  (will eventually updating of snow temperature in this case)
 !!
 !! @param[inout] RLA Latitude on the cubed-sphere tile
 !! @param[inout] RLO Longitude on the cubed-sphere tile
 !! @param[inout] STC_STATE
 !! @param[in] SOILSNOW_TILE Land mask for increments on the cubed-sphere tile
 !! @param[in] SOILSNOW_FG_TILE First guess land mask for increments on the cubed-sphere tile
 !! @param[in] LENSFC Number of land points on a tile
 !! @param[in] LSOIL Number of soil layers
 !! @param[in] IDIM 'I' dimension of a tile
 !! @param[in] JDIM 'J' dimension of a tile
 !! @param[in] MYRANK MPI rank number
 !!
 !! @author Clara Draper. @date March 2021

subroutine add_increment_soil(rla,rlo,stc_state,soilsnow_tile, soilsnow_fg_tile, &
                        lensfc,lsoil,idim,jdim, myrank)

    use utils
    use gdswzd_mod
    use read_write_data, only : idim_gaus, jdim_gaus, &
                             stc_inc_gaus, soilsnow_gaus
    use mpi

    implicit none

    integer, intent(in)      :: lensfc, lsoil, idim, jdim, myrank

    integer, intent(in)         :: soilsnow_tile(lensfc), soilsnow_fg_tile(lensfc)
    real, intent(inout)      :: rla(lensfc), rlo(lensfc)
    real, intent(inout)      :: stc_state(lensfc, lsoil)

    integer                  :: iopt, nret, kgds_gaus(200)
    integer                  :: igaus, jgaus, ij
    integer                  :: mask_tile, mask_fg_tile
    integer                  :: itile, jtile
    integer                  :: j, ierr
    integer                  :: igausp1, jgausp1
    real                     :: fill

    integer, allocatable     :: id1(:,:), id2(:,:), jdc(:,:)

    real                     :: wsum
    real                     :: stc_inc(lsoil)
    real, allocatable        :: xpts(:), ypts(:), lats(:), lons(:)
    real, allocatable        :: dum2d(:,:), lats_rad(:), lons_rad(:)
    real, allocatable        :: agrid(:,:,:), s2c(:,:,:)

    integer                  :: k, nother, nsnowupd, nnosoilnear, nsoilupd, nsnowchange
    logical                  :: gaus_has_soil
        
    integer, parameter       :: lsoil_incr=3 ! number of layers to add incrments to
                                                
    ! this produces the same lat/lon as can be read in from the file

    kgds_gaus     = 0
    kgds_gaus(1)  = 4          ! oct 6 - type of grid (gaussian)
    kgds_gaus(2)  = idim_gaus  ! oct 7-8 - # pts on latitude circle
    kgds_gaus(3)  = jdim_gaus
    kgds_gaus(4)  = 90000      ! oct 11-13 - lat of origin
    kgds_gaus(5)  = 0          ! oct 14-16 - lon of origin
    kgds_gaus(6)  = 128        ! oct 17 - resolution flag
    kgds_gaus(7)  = -90000     ! oct 18-20 - lat of extreme point
    kgds_gaus(8)  = nint(-360000./float(idim_gaus))  ! oct 21-23 - lon of extreme point
    kgds_gaus(9)  = nint((360.0 / float(idim_gaus))*1000.0)
                            ! oct 24-25 - longitude direction incr.
    kgds_gaus(10) = jdim_gaus/2     ! oct 26-27 - number of circles pole to equator
    kgds_gaus(12) = 255        ! oct 29 - reserved
    kgds_gaus(20) = 255        ! oct 5  - not used, set to 255

    print*
    print*,'adjust soil temperature using gsi increments on gaussian grid'
    print*,'adjusting first ', lsoil_incr, ' surface layers only'

    !----------------------------------------------------------------------
    ! call gdswzd to compute the lat/lon of each gsi gaussian grid point.
    !----------------------------------------------------------------------

    iopt = 0
    fill = -9999.
    allocate(xpts(idim_gaus*jdim_gaus))
    allocate(ypts(idim_gaus*jdim_gaus))
    allocate(lats(idim_gaus*jdim_gaus))
    allocate(lons(idim_gaus*jdim_gaus))
    xpts = fill
    ypts = fill
    lats = fill
    lons = fill

    call gdswzd(kgds_gaus,iopt,(idim_gaus*jdim_gaus),fill,xpts,ypts,lons,lats,nret)

    if (nret /= (idim_gaus*jdim_gaus)) then
    print*,'fatal error: problem in gdswzd. stop.'
    call mpi_abort(mpi_comm_world, 12, ierr)
    endif

    deallocate (xpts, ypts)

    allocate(dum2d(idim_gaus,jdim_gaus))
    dum2d = reshape(lats, (/idim_gaus,jdim_gaus/) )
    deallocate(lats)

    allocate(lats_rad(jdim_gaus))

    do j = 1, jdim_gaus
    lats_rad(j) = dum2d(1,jdim_gaus-j+1) * 3.1415926 / 180.0
    enddo

    dum2d = reshape(lons, (/idim_gaus,jdim_gaus/) )
    deallocate(lons)
    allocate(lons_rad(idim_gaus))
    lons_rad = dum2d(:,1) * 3.1415926 / 180.0

    deallocate(dum2d)

    allocate(agrid(idim,jdim,2))
    agrid(:,:,1) = reshape (rlo, (/idim,jdim/) )
    agrid(:,:,2) = reshape (rla, (/idim,jdim/) )
    agrid        = agrid * 3.1415926 / 180.0

    allocate(id1(idim,jdim))
    allocate(id2(idim,jdim))
    allocate(jdc(idim,jdim))
    allocate(s2c(idim,jdim,4))

    !----------------------------------------------------------------------
    ! compute bilinear weights for each model point from the nearest
    ! four gsi/gaussian points.  does not account for mask.  that
    ! happens later.
    !----------------------------------------------------------------------

    call remap_coef( 1, idim, 1, jdim, idim_gaus, jdim_gaus, &
                  lons_rad, lats_rad, id1, id2, jdc, s2c, agrid )

    deallocate(lons_rad, lats_rad, agrid)
    !
    ! initialize variables for counts statitics to be zeros
    !

    !
    nother = 0 ! grid cells not land
    nsnowupd = 0  ! grid cells with snow (temperature not yet updated)
    nsnowchange = 0  ! grid cells where no temp upd made, because snow occurence changed
    nnosoilnear = 0 ! grid cells where model has soil, but 4 closest gaus grids don't
                 ! (no update made here)
    nsoilupd = 0


    ij_loop : do ij = 1, lensfc

        mask_tile    = soilsnow_tile(ij)
        mask_fg_tile = soilsnow_fg_tile(ij)

        !----------------------------------------------------------------------
        ! mask: 1  - soil, 2 - snow, 0 - land-ice, -1 - not land
        !----------------------------------------------------------------------

        if (mask_tile <= 0) then ! skip if neither soil nor snow
         nother = nother + 1
         cycle ij_loop
        endif


        !  get i,j index on array of (idim,jdim) from known ij

        jtile = (ij-1) / idim + 1
        itile = mod(ij,idim)
        if (itile==0) itile = idim

        !----------------------------------------------------------------------
        ! if the snow analysis has chnaged to occurence of snow, skip the
        ! temperature analysis
        !----------------------------------------------------------------------

        if ((mask_fg_tile == 2 .and. mask_tile == 1) .or. &
            (mask_fg_tile == 1 .and. mask_tile == 2) ) then
         nsnowchange = nsnowchange + 1
         cycle ij_loop
        endif

        !----------------------------------------------------------------------
        !  do update to soil temperature grid cells, using bilinear interp
        !----------------------------------------------------------------------

        if (mask_tile == 1) then
           ! these are the four nearest grid cells on the gaus grid
           igaus   = id1(itile,jtile)
           jgaus   = jdc(itile,jtile)
           igausp1 = id2(itile,jtile)
           jgausp1 = jdc(itile,jtile)+1

        ! make sure gaus grid has soil nearby
           gaus_has_soil = .false.
           if (soilsnow_gaus(igaus,jgaus)     == 1 .or. &
               soilsnow_gaus(igausp1,jgaus)   == 1 .or. &
               soilsnow_gaus(igausp1,jgausp1) == 1 .or. &
               soilsnow_gaus(igaus,jgausp1)   == 1)  gaus_has_soil = .true.
           
           if (.not. gaus_has_soil) then
             nnosoilnear = nnosoilnear + 1
             cycle ij_loop
           endif

        ! calcualate weighted increment over nearby grid cells that have soil

        ! Draper: to-do, code adding increments to soil moisture.
        !         will require converting to soil wetness index first
        !         (need to add soil properties to the increment file)

           nsoilupd = nsoilupd + 1

           stc_inc = 0.0
           wsum  = 0.0

           if (soilsnow_gaus(igaus,jgaus) == 1) then
             do k = 1, lsoil_incr
                 stc_inc(k)  = stc_inc(k) + (s2c(itile,jtile,1) * stc_inc_gaus(k,igaus,jgaus))
             enddo
             wsum  = wsum + s2c(itile,jtile,1)
           endif

           if (soilsnow_gaus(igausp1,jgaus) == 1) then
             do k = 1, lsoil_incr
                 stc_inc(k) = stc_inc(k) + (s2c(itile,jtile,2) * stc_inc_gaus(k,igausp1,jgaus))
             enddo
             wsum  = wsum + s2c(itile,jtile,2)
           endif

           if (soilsnow_gaus(igausp1,jgausp1) == 1) then
             do k = 1, lsoil_incr
                 stc_inc(k) = stc_inc(k) + (s2c(itile,jtile,3) * stc_inc_gaus(k,igausp1,jgausp1))
             enddo
             wsum  = wsum + s2c(itile,jtile,3)
           endif

           if (soilsnow_gaus(igaus,jgausp1) == 1) then
             do k = 1, lsoil_incr
                 stc_inc(k) = stc_inc(k) + (s2c(itile,jtile,4) * stc_inc_gaus(k,igaus,jgausp1))
             enddo
             wsum  = wsum + s2c(itile,jtile,4)
           endif

        ! add increment
           do k = 1, lsoil_incr
             stc_inc(k) = stc_inc(k) / wsum
             stc_state(ij,k) = stc_state(ij,k) + stc_inc(k)
        ! todo, apply some bounds?
           enddo

        ! don't update soil states if snow present.
        elseif(mask_tile==2) then
           nsnowupd = nsnowupd + 1

        endif ! if soil/snow point

    enddo ij_loop

    write(*,'(a,i2)') 'statistics of grids number processed for rank : ', myrank
    write(*,'(a,i8)') ' soil grid cells updated = ',nsoilupd
    write(*,'(a,i8)') ' (not updated) soil grid cells, no soil nearby on gsi grid = ',nnosoilnear
    write(*,'(a,i8)') ' (not updated) soil grid cells, change in presence of snow = ', nsnowchange
    write(*,'(a,i8)') ' (not updated yet) snow grid cells = ', nsnowupd
    write(*,'(a,i8)') ' grid cells, without soil or snow = ', nother

    nother = 0 ! grid cells not land
    nsnowupd = 0  ! grid cells where no temp upd made, because snow occurence changed
    nnosoilnear = 0 ! grid cells where model has soil, but 4 closest gaus grids don't
                 ! (no update made here)
    nsoilupd = 0

    deallocate(id1, id2, jdc, s2c)

end subroutine add_increment_soil

 !> Add snow depth increment to model snow depth state,
 !! and limit output to be non-negative. JEDI increments are
 !! calculated globally, so must be screened to land-only locations 
 !! here.
 !!
 !! @param[in] lensfc Number of land points on this tile
 !! @param[in] snd_inc Soil depth increments
 !! @param[in] mask Land mask for increments
 !! @param[inout] snd Soil depth background (in), and analysis (out) 
 !! 
 !! @author Clara Draper. @date August 2021

subroutine add_increment_snow(snd_inc,mask,lensfc,snd)

    implicit none

    integer, intent(in)      :: lensfc
    real, intent(in)         :: snd_inc(lensfc)
    integer, intent(in)      :: mask(lensfc)
    real, intent(inout)      :: snd(lensfc)

    integer                  :: i


    do i =1, lensfc
        if (mask(i) > 0) then ! if land
                snd(i) = max( snd(i) + snd_inc(i), 0.)
        endif
    enddo

end subroutine add_increment_snow

!> Calculate soil mask for land on model grid.
!! Output is 1  - soil, 2 - snow-covered, 0 - land ice, -1  not land.
!!
!! @param[in] lensfc  Number of land points for this tile 
!! @param[in] veg_type_landice Value of vegetion class that indicates land-ice
!! @param[in] smc Model soil moisture.
!! @param[in] swe Model snow water equivalent
!! @param[in] vtype Model vegetation type
!! @param[out] mask Land mask for increments
!! @author Clara Draper @date March 2021
subroutine calculate_landinc_mask(smc,swe,vtype,lensfc,veg_type_landice,mask)
 
    implicit none

    integer, intent(in)           :: lensfc, veg_type_landice
    real, intent(in)              :: smc(lensfc), swe(lensfc)
    real, intent(in)              :: vtype(lensfc)
    integer, intent(out)          :: mask(lensfc)

    integer :: i

    mask = -1 ! not land

    ! land (but not land-ice)
    do i=1,lensfc
        if (smc(i) .LT. 1.0) then
          if (swe(i) .GT. 0.001) then ! snow covered land
                mask(i) = 2
          else                        ! non-snow covered land
                mask(i) = 1
          endif
        end if ! else should work here too
        if ( nint(vtype(i)) ==  veg_type_landice  ) then ! land-ice
                mask(i) = 0
        endif
    end do

end subroutine calculate_landinc_mask

!> Make adjustments to dependent variables after applying land increments.
!! These adjustments are model-dependent, and are currently only coded
!! for Noah LSM.
!! For Noah LSM, copy relevent code blocks from model code (same as has
!! been done in sfc_sub).
!! Here: adjust (frozen) soil moisture to be consistent with changes in
!! soil temperature from DA

!> @param[in] lsm Integer code for the LSM
!! @param[in] isot Integer code for the soil type data set
!! @param[in] ivegsrc Integer code for the vegetation type data set
!! @param[in] lensfc Number of land points for this tile
!! @param[in] lsoil Number of soil layers
!! @param[in] rsoiltype Array of input soil types
!! @param[in] mask Mask indicating surface type
!! @param[in] stc_bck Background soil temperature states
!! @param[in] stc_anl Analysis soil temperature states
!! @param[inout] smc_adj Soil moisture state to be adjusted
!! @param[inout] slc_adj Liquid soil moisture states to be adjusted
!! @author Clara Draper @date April 2021

subroutine apply_land_da_adjustments_stc(lsm, isot, ivegsrc,lensfc, &
                 lsoil, rsoiltype, mask, stc_bck, stc_anl, smc_adj, slc_adj )

    use mpi
    use set_soilveg_snippet_mod, only: set_soilveg
    use sflx_snippet,    only: frh2o

    implicit none
 
    integer, intent(in)           :: lsm, lensfc, lsoil, isot, ivegsrc
    real, intent(in)              :: rsoiltype(lensfc) ! soil types, as real
    integer, intent(in)              :: mask(lensfc)
    real, intent(in)              :: stc_bck(lensfc, lsoil) , stc_anl(lensfc, lsoil)
    real, intent(inout)           :: smc_adj(lensfc,lsoil), slc_adj(lensfc,lsoil)
    

    logical                       :: frzn_bck, frzn_anl

    integer                       :: i, l, n_freeze, n_thaw, ierr
    integer                       :: myrank, soiltype, iret

    real                          :: slc_new

    real, parameter               :: tfreez=273.16 !< con_t0c  in physcons
    real, dimension(30)           :: maxsmc, bb, satpsi

    call mpi_comm_rank(mpi_comm_world, myrank, ierr)

    if (lsm .NE. lsm_noah) then
        print *, 'FATAL ERROR: apply_land_da_adjustments not coded for models other than noah', lsm
        call mpi_abort(mpi_comm_world, 10, ierr)
    endif
       
    ! initialise soil properties
    call set_soilveg(isot, ivegsrc, maxsmc, bb, satpsi, iret)
    if (iret < 0) then
        print *, 'FATAL ERROR: problem in set_soilveg'
        call mpi_abort(mpi_comm_world, 10, ierr)
    endif

    print *, 'Adjusting smc after stc DA update'

    n_freeze = 0
    n_thaw = 0
    
    do i=1,lensfc
      if (mask(i) > 0) then ! if soil location
        do l = 1, lsoil
           frzn_bck = (stc_bck(i,l) .LT. tfreez )
           frzn_anl = (stc_anl(i,l) .LT. tfreez )

           if (frzn_bck .eqv. frzn_anl) then
                cycle
           elseif (frzn_bck .and. .not. frzn_anl) then
                n_thaw = n_thaw + 1
           else
                n_freeze = n_freeze + 1
           endif

           ! make adjustment (same routine for both)
           soiltype = nint(rsoiltype(i))
           ! bb and maxsmc are in the namelist_soilveg, need soiltype index
           call frh2o(stc_anl(i,l), smc_adj(i,l),slc_adj(i,l), maxsmc(soiltype), &
                      bb(soiltype), satpsi(soiltype),slc_new)

           slc_adj(i,l) = max( min( slc_new, smc_adj(i,l)), 0.0 )
        enddo
      endif
    enddo
    print *, 'adjusted: ', n_thaw,' thawed,', n_freeze, ' frozen'

end subroutine apply_land_da_adjustments_stc

!> Make adjustments to dependent variables after applying land increments.
!! These adjustments are model-dependent, and are currently only coded
!! for Noah LSM.
!! Here: adjust SWE to be consistent with updated SND, using snow density 
!! from the forecast.

!> @param[in] lsm Integer code for the LSM
!! @param[in] lensfc  Number of land points for this tile
!! @param[in] mask Land mask for increments
!! @param[in] swe_bck Background SWE 
!! @param[in] snd_bck Background snow depth
!! @param[in] snd_anl Analysis snow depth 
!! @param[inout] swe_adj SWE to be adjusted 
!! @author Clara Draper @date August 2021

subroutine apply_land_da_adjustments_snd(lsm, lensfc, mask, swe_bck, snd_bck, snd_anl, swe_adj)

    use mpi
    use bulk_snow_module, only: calc_density

    implicit none

    integer, intent(in) :: lsm, lensfc
    integer, intent(in) :: mask(lensfc)
    real, intent(in)    :: swe_bck(lensfc), snd_bck(lensfc)
    real, intent(in)    :: snd_anl(lensfc)
    real, intent(inout)    :: swe_adj(lensfc)

    integer  :: ierr, myrank, i

    real     :: density(lensfc)

    call mpi_comm_rank(mpi_comm_world, myrank, ierr)

    if (lsm .NE. lsm_noah) then
        print *, 'FATAL ERROR: apply_land_da_adjustments not coded for models other than noah', lsm
        call mpi_abort(mpi_comm_world, 10, ierr)
    endif

   ! calculate snow density from forecasts
   call calc_density(lensfc, mask, swe_bck, snd_bck, myrank, density)

   do i =1, lensfc
        if ( mask(i)>0 ) then
                swe_adj(i) = snd_anl(i)*density(i)
        endif
   enddo
   

end subroutine apply_land_da_adjustments_snd

end module land_increments
