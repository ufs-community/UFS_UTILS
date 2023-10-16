!> @file
!! @brief Compute lake fraction and depth.
!! @author Ning Wang

!> This program computes lake fraction and depth numbers for FV3 cubed sphere 
!! grid cells, from a high resolution lat/lon data set.
!! 
!! @author Ning Wang @date July 2018
!! 
!!  - Shan Sun, Aug. 2018: Added Caspian Sea and Aral Sea to the lake fraction 
!!                        and lake depth fields.
!!  - Shan Sun, Dec. 2018: Added round up and round down with respect to a 
!!                        numerical minimum value and a cut-off value, for lake 
!!                        fraction number. 
!!  - Ning Wang, Apr. 2019: Extended the program to process the same lake data 
!!                         for FV3 stand-alone regional (SAR) model.
!!     
!! @return 0 for success.
!#define DIAG_N_VERBOSE
#define ADD_ATT_FOR_NEW_VAR
PROGRAM lake_frac
    USE netcdf
    IMPLICIT NONE

    CHARACTER(len=256) :: sfcdata_path
    INTEGER :: cs_res, ncsmp, ncscp, i
    INTEGER :: res_x, res_y

    INTEGER*1, ALLOCATABLE :: lakestatus(:)
    INTEGER*2, ALLOCATABLE :: lakedepth(:)
    REAL, ALLOCATABLE :: cs_grid(:,:)
    REAL, ALLOCATABLE :: cs_lakestatus(:), cs_lakedepth(:)
    REAL, ALLOCATABLE :: src_grid_lon(:), src_grid_lat(:)

    INTEGER :: tile_req, tile_beg, tile_end
    REAL :: lake_cutoff

    INTEGER, PARAMETER :: nlat = 21600, nlon = 43200
    REAL, PARAMETER :: d2r = acos(-1.0) / 180.0
    REAL, PARAMETER :: r2d = 180.0 /acos(-1.0) 
    REAL, PARAMETER :: pi = acos(-1.0) 
    REAL*8, PARAMETER :: gppdeg = 120.0
    REAL*8, PARAMETER :: delta = 1.0 / 120.0

    CHARACTER(len=32) :: arg, lakestatus_srce, lakedepth_srce
    CHARACTER(len=256) :: lakedata_path
    INTEGER :: stat
    
    CALL getarg(0, arg) ! get the program name
    IF (iargc() /= 5 .AND. iargc() /= 6) THEN
      PRINT*, 'Usage: ', trim(arg), & 
       ' [tile_num (0:all tiles, 7:regional)] [resolution (48,96, ...)] &
         [lake data path] [lake status source] [lake depth source]'  
      PRINT*, 'Or: ', trim(arg), & 
       ' [tile_num (0:all tiles, 7:regional)] [resolution (48,96, ...)] &
         [lake data path] [lake status source] [lake depth source] [lake_cutoff]'
      STOP -1
    ENDIF
    CALL getarg(1, arg)
    READ(arg,*,iostat=stat) tile_req
    CALL getarg(2, arg)
    READ(arg,*,iostat=stat) cs_res
    CALL getarg(3, lakedata_path)
    CALL getarg(4, lakestatus_srce)
    CALL getarg(5, lakedepth_srce)

    IF (iargc() == 5) THEN
      lake_cutoff = 0.20
    ELSE
      CALL getarg(6, arg)
      READ(arg,*,iostat=stat) lake_cutoff
    ENDIF

    PRINT*, 'lake status source:', trim(lakestatus_srce) 
    PRINT*, 'lake depth source:', trim(lakedepth_srce) 
    PRINT*, 'lake cutoff:', lake_cutoff 

    IF (tile_req == 0) THEN
      tile_beg = 1; tile_end = 6
      PRINT*, 'Process tile 1 - 6 at resolution C',cs_res
    ELSE IF (tile_req /= 7) THEN
      tile_beg = tile_req; tile_end = tile_req
      PRINT*, 'Process tile',tile_req, 'at resolution C',cs_res
    ELSE
      tile_beg = 1; tile_end = 1
      PRINT*, 'Process regional tile (tile', tile_req, ') at resolution C',cs_res
    ENDIF

    ! read in grid spec data for each tile and concatenate them together

    ncsmp = (2*cs_res+1)*(2*cs_res+1)*6   
    IF (tile_req /= 7) THEN
      PRINT*, 'Read in cubed sphere grid information ... ',ncsmp,'pairs of lat/lons'
    ENDIF

    IF (tile_req /= 7) THEN
      ALLOCATE(cs_grid(ncsmp, 2))
      CALL  read_cubed_sphere_grid(cs_res, cs_grid) 
    ELSE
      CALL  read_cubed_sphere_reg_grid(cs_res, cs_grid, 3, res_x, res_y) 
    ENDIF

    ! allocate and compute source grid 
    ALLOCATE(src_grid_lon(nlon), src_grid_lat(nlat))

    IF (lakestatus_srce == "GLDBV2" .OR. lakestatus_srce == "GLDBV3") THEN 
      !  GLDB data points are at the lower right corners of the grid cells
      DO i = 1, nlon
        src_grid_lon(i) = -180.0 + delta * i
      ENDDO
      DO i = 1, nlat
        src_grid_lat(i) = 90.0 - delta * i 
      ENDDO
    ENDIF 

    IF (lakestatus_srce == "MODISP" .OR. lakestatus_srce == "VIIRS") THEN 
      !  GLDB data points are at the upprt left corners of the grid cells
      DO i = 1, nlon
        src_grid_lon(i) = -180.0 + delta * (i-1)
      ENDDO
      DO i = 1, nlat
        src_grid_lat(i) = 90.0 - delta * (i-1) 
      ENDDO
    ENDIF

    ! read in lake data file
    lakedata_path = trim(lakedata_path) // "/"
    ALLOCATE(lakestatus(nlon*nlat),lakedepth(nlon*nlat)) 
    PRINT*, 'Read in lake data file ...'
    CALL read_lakedata(lakedata_path,lakestatus,lakedepth,nlat,nlon)

    ! calculate fraction numbers for all cs-cells
    ncscp = cs_res*cs_res*6
    ALLOCATE(cs_lakestatus(ncscp))
    ALLOCATE(cs_lakedepth(ncscp))

    PRINT*, 'Calculate lake fraction and average depth for cubed-sphere cells ...'
    CALL cal_lake_frac_depth(lakestatus,cs_lakestatus,lakedepth,cs_lakedepth)

    ! write lake status (in terms of fraction) and lake depth(average, in meters)
    ! to a netcdf file
    IF (tile_req /= 7) THEN
      PRINT*, 'Write lake fraction/depth on cubed sphere grid cells to netCDF files ...'
      CALL write_lakedata_to_orodata(cs_res, cs_lakestatus, cs_lakedepth) 
    ELSE
      PRINT*, 'Write lake fraction/depth on regional FV3 grid cells to a netCDF file ...'
      CALL write_reg_lakedata_to_orodata(cs_res, res_x, res_y, cs_lakestatus, cs_lakedepth) 
    ENDIF

    DEALLOCATE(cs_lakestatus,cs_lakedepth)
    DEALLOCATE(cs_grid)
    DEALLOCATE(lakestatus,lakedepth)
    DEALLOCATE(src_grid_lat, src_grid_lon)

    STOP 0
CONTAINS

!> Calculate lake fraction and depth on the model grid from
!! high-resolution data.
!!
!! @param[in] lakestat High-resolution lake status code.
!! @param[in] lakedpth High-resolution lake depth.
!! @param[out] cs_lakestat Lake fraction on the model grid.
!! @param[out] cs_lakedpth Lake depth on the model grid.
!! @author Ning Wang
SUBROUTINE cal_lake_frac_depth(lakestat,cs_lakestat,lakedpth,cs_lakedpth)
    INTEGER*1, INTENT(IN) :: lakestat(:)
    INTEGER*2, INTENT(IN) :: lakedpth(:)
    REAL, INTENT(OUT) :: cs_lakestat(:), cs_lakedpth(:)

    REAL*8 lolf(2), lort(2), uplf(2), uprt(2), sd_ltmn(4), sd_ltmx(4)
    REAL*8 :: v(2,4), p(2)
    REAL :: latmin1, latmax1
    REAL :: latmin, latmax, lonmin, lonmax, lontmp, lat_sz_max, lon_sz_max
    INTEGER :: tile_num, i, j, gp, row, col, cs_grid_idx, cs_data_idx
    INTEGER :: sidex_res, sidey_res, sidex_sz, sidey_sz 
    INTEGER :: stride_lat, stride_lon
    INTEGER :: src_grid_lat_beg,src_grid_lat_end,src_grid_lon_beg,src_grid_lon_end
    INTEGER :: src_grid_lon_beg1,src_grid_lon_end1,src_grid_lon_beg2,src_grid_lon_end2
    INTEGER :: grid_ct, lake_ct, co_gc, tmp

    INTEGER*1 :: lkst
    INTEGER*2 :: lkdp
    REAL*8 :: lake_dpth_sum, lake_avg_frac
    LOGICAL :: two_section, enclosure_cnvx 

    IF (tile_req /= 7) THEN
      sidex_res = cs_res; sidey_res = cs_res
    ELSE
      sidex_res = res_x; sidey_res = res_y
    ENDIF

    sidex_sz = 2*sidex_res+1; sidey_sz = 2*sidey_res+1

    stride_lat = 1

    lat_sz_max = 0.0
    lon_sz_max = 0.0

    cs_lakestat = 0.0

    DO tile_num = tile_beg, tile_end
      row = 2 + sidex_sz*(tile_num-1); col = 2
      DO gp = 1, sidex_res*sidey_res
        two_section = .false.
        cs_grid_idx = (row-1)*sidex_sz+col
        cs_data_idx = (tile_num-1)*sidex_res*sidey_res+gp
        IF (abs(cs_grid(cs_grid_idx,1)) > 80.0 ) THEN !ignore lakes in very high latitude
          cs_lakestat(cs_data_idx) = 0.0
          cs_lakedpth(cs_data_idx) = 0.0
          ! move to next cs cell
          col = col + 2
          IF (col > sidex_sz) THEN
            col = 2
            row = row + 2
          ENDIF
          CYCLE
        ENDIF
        ! get the four corners of the cs cell 
        lolf(1) = cs_grid(cs_grid_idx-sidex_sz-1, 1)
        lolf(2) = cs_grid(cs_grid_idx-sidex_sz-1, 2)
        IF (lolf(2) > 180.0) lolf(2) = lolf(2) - 360.0
        lort(1) = cs_grid(cs_grid_idx-sidex_sz+1, 1)
        lort(2) = cs_grid(cs_grid_idx-sidex_sz+1, 2)
        IF (lort(2) > 180.0) lort(2) = lort(2) - 360.0
        uplf(1) = cs_grid(cs_grid_idx+sidex_sz-1,1)
        uplf(2) = cs_grid(cs_grid_idx+sidex_sz-1,2)
        IF (uplf(2) > 180.0) uplf(2) = uplf(2) - 360.0
        uprt(1) = cs_grid(cs_grid_idx+sidex_sz+1,1)
        uprt(2) = cs_grid(cs_grid_idx+sidex_sz+1,2)

        v(1,1) = lolf(1); v(2,1) = lolf(2) 
        v(1,2) = lort(1); v(2,2) = lort(2) 
        v(1,3) = uprt(1); v(2,3) = uprt(2) 
        v(1,4) = uplf(1); v(2,4) = uplf(2) 
        v(:,:) = v(:,:) * d2r

        IF (uprt(2) > 180.0) uprt(2) = uprt(2) - 360.0
       ! gather the candidate indices in lakestat
#ifdef LIMIT_CAL
        CALL find_limit (lolf, lort, sd_ltmn(1), sd_ltmx(1))
        CALL find_limit (lort, uprt, sd_ltmn(2), sd_ltmx(2))
        CALL find_limit (uprt, uplf, sd_ltmn(3), sd_ltmx(3))
        CALL find_limit (uplf, lolf, sd_ltmn(4), sd_ltmx(4))
        latmin = min(sd_ltmn(1),min(sd_ltmn(2),min(sd_ltmn(3),sd_ltmn(4))))    
        latmax = max(sd_ltmx(1),max(sd_ltmx(2),max(sd_ltmx(3),sd_ltmx(4))))    
#endif
        latmin = min(lolf(1),min(lort(1),min(uplf(1),uprt(1))))    
        latmax = max(lolf(1),max(lort(1),max(uplf(1),uprt(1))))    
        lonmin = min(lolf(2),min(lort(2),min(uplf(2),uprt(2))))    
        lonmax = max(lolf(2),max(lort(2),max(uplf(2),uprt(2))))    
!        lat_sz_max = max(lat_sz_max, (latmax-latmin))
!        lon_sz_max = max(lon_sz_max, (lonmax-lonmin))

        src_grid_lat_beg = nint((90.0-latmax)*gppdeg+0.5)
        src_grid_lat_end = nint((90.0-latmin)*gppdeg+0.5)
        src_grid_lon_beg = nint((180.0+lonmin)*gppdeg+0.5) 
        src_grid_lon_end = nint((180.0+lonmax)*gppdeg+0.5) 
       
        IF (src_grid_lat_beg > src_grid_lat_end) THEN
          tmp = src_grid_lat_beg
          src_grid_lat_beg = src_grid_lat_end
          src_grid_lat_end = tmp
        ENDIF
        IF (src_grid_lon_beg > src_grid_lon_end) THEN
          tmp = src_grid_lon_beg
          src_grid_lon_beg = src_grid_lon_end
          src_grid_lon_end = tmp
        ENDIF
        IF ((src_grid_lon_end - src_grid_lon_beg) > nlon*0.75) THEN
          two_section = .true.
          src_grid_lon_beg1 = src_grid_lon_end
          src_grid_lon_end1 = nlon
          src_grid_lon_beg2 = 1
          src_grid_lon_end2 = src_grid_lon_beg
        ENDIF

#ifdef DIAG_N_VERBOSE
        PRINT*, 'cell centre lat/lon =', &
          gp, cs_res*cs_res, cs_grid(cs_grid_idx,1),cs_grid(cs_grid_idx,2)
        PRINT*, 'lat index range and stride', &
          src_grid_lat_beg,src_grid_lat_end,stride_lat
        PRINT*, 'lat range ',  &
          src_grid_lat(src_grid_lat_beg),src_grid_lat(src_grid_lat_end)
#endif
        lake_ct = 0; grid_ct = 0
        lake_dpth_sum = 0.0
        lake_avg_frac = 0.0
        DO j = src_grid_lat_beg, src_grid_lat_end, stride_lat
          stride_lon = int(1.0/cos(src_grid_lat(j)*d2r)*REAL(stride_lat))
#ifdef DIAG_N_VERBOSE
          IF (j == src_grid_lat_beg) THEN
            PRINT*, 'lon index range and stride', &
              src_grid_lon_beg,src_grid_lon_end,stride_lon
            PRINT*, 'lon range ', &
              src_grid_lon(src_grid_lon_beg),src_grid_lon(src_grid_lon_end)
            IF (two_section .eqv. .true.) THEN
              PRINT*, 'section1 index lon range and stride', &
                src_grid_lon_beg1,src_grid_lon_end1,stride_lon
              PRINT*, 'section1 lon range ', &
                src_grid_lon(src_grid_lon_beg1),src_grid_lon(src_grid_lon_end1)
              PRINT*, 'section2 index lon range and stride', &
                src_grid_lon_beg2,src_grid_lon_end2,stride_lon
              PRINT*, 'section2 lon range ', &
                src_grid_lon(src_grid_lon_beg2),src_grid_lon(src_grid_lon_end2)
            ENDIF
          ENDIF
#endif
          IF (two_section .eqv. .false.) THEN
            DO i =  src_grid_lon_beg, src_grid_lon_end, stride_lon
              p(1) = src_grid_lat(j); p(2) = src_grid_lon(i) 
              p(:) = p(:)*d2r
              IF(enclosure_cnvx(v, 4, p, co_gc) .eqv. .true.) THEN
                grid_ct = grid_ct+1
                lkst = lakestat((j-1)*nlon+i); lkdp = lakedpth((j-1)*nlon+i)
                CALL lake_cell_comp(lkst, lkdp, lake_ct, lake_avg_frac, lake_dpth_sum)
              ENDIF   
            ENDDO
          ELSE 
            DO i =  src_grid_lon_beg1, src_grid_lon_end1, stride_lon
              p(1) = src_grid_lat(j); p(2) = src_grid_lon(i) 
              p(:) = p(:)*d2r
              IF(enclosure_cnvx(v, 4, p, co_gc) .eqv. .true.) THEN
                grid_ct = grid_ct+1
                lkst = lakestat((j-1)*nlon+i); lkdp = lakedpth((j-1)*nlon+i)
                CALL lake_cell_comp(lkst, lkdp, lake_ct, lake_avg_frac, lake_dpth_sum)
              ENDIF   
            ENDDO
            DO i =  src_grid_lon_beg2, src_grid_lon_end2, stride_lon
              p(1) = src_grid_lat(j); p(2) = src_grid_lon(i) 
              p(:) = p(:)*d2r
              IF(enclosure_cnvx(v, 4, p, co_gc) .eqv. .true.) THEN
                grid_ct = grid_ct+1
                lkst = lakestat((j-1)*nlon+i); lkdp = lakedpth((j-1)*nlon+i)
                CALL lake_cell_comp(lkst, lkdp, lake_ct, lake_avg_frac, lake_dpth_sum)
              ENDIF   
            ENDDO
          ENDIF
        ENDDO
        IF (lakestatus_srce == "GLDBV3" .OR. lakestatus_srce == "GLDBV2" .OR. &
            lakestatus_srce == "VIIRS" ) THEN 
          cs_lakestat(cs_data_idx)=REAL(lake_ct)/REAL(grid_ct)
        ENDIF
        IF (lakestatus_srce == "MODISP") THEN 
          cs_lakestat(cs_data_idx)=lake_avg_frac/REAL(grid_ct)
        ENDIF
        IF (lake_ct /= 0) THEN
          cs_lakedpth(cs_data_idx)=lake_dpth_sum/REAL(lake_ct)/10.0 !convert to meter
        ELSE
          cs_lakedpth(cs_data_idx)=0.0
        ENDIF
#ifdef DIAG_N_VERBOSE
        PRINT*, 'tile_num, row, col:', tile_num, row, col
        PRINT*, 'grid_ct, lake_ct = ', grid_ct, lake_ct
        PRINT*, 'lake_frac= ', cs_lakestat(cs_data_idx) 
        PRINT*, 'lake_depth (avg) = ', cs_lakedpth(cs_data_idx) 
#endif

       ! move to the next control volume
        col = col + 2
        IF (col > sidex_sz) THEN
          col = 2
          row = row + 2
        ENDIF
      ENDDO
      PRINT "('*'$)"  ! progress '*'
    ENDDO
    PRINT*, ''

END SUBROUTINE cal_lake_frac_depth 

!> Compute cumulatively the lake fraction and lake depth for a cell 
!!
!! @param[in] lkst lake status value from a grid point in the source data.
!! @param[in] lkdp lake depth value from a grid point in the source data.
!! @param[out] lake_ct lake points number accumulated for the cell 
!! @param[out] lake_avg_frac lake fraction value accumulated for the cell. 
!! @param[out] lake_dpth_sum is the lake depth value accumulated for the cell. 
!! @author Ning Wang
SUBROUTINE lake_cell_comp(lkst, lkdp, lake_ct, lake_avg_frac, lake_dpth_sum)
    INTEGER*1, INTENT(IN) :: lkst
    INTEGER*2, INTENT(IN) :: lkdp
    INTEGER, INTENT(OUT) :: lake_ct
    REAL*8, INTENT(OUT) :: lake_avg_frac, lake_dpth_sum
    
    IF (lkst /= 0) THEN ! lake point 
      lake_ct = lake_ct+1
      IF (lakestatus_srce == "GLDBV3" .OR. lakestatus_srce == "GLDBV2") THEN 
        IF (lkdp <= 0) THEN
          IF (lkst == 4) THEN
            lake_dpth_sum = lake_dpth_sum+30.0
          ELSE
            lake_dpth_sum = lake_dpth_sum+100.0
          ENDIF
        ELSE
          lake_dpth_sum = lake_dpth_sum+REAL(lkdp)
        ENDIF
      ENDIF
      IF (lakestatus_srce == "MODISP" .OR. lakestatus_srce == "VIIRS") THEN 
        lake_avg_frac = lake_avg_frac + REAL(lkst) / 100.0 
        IF (lkdp <= 0) THEN
          lake_dpth_sum = lake_dpth_sum+100.0
        ELSE
          lake_dpth_sum = lake_dpth_sum+REAL(lkdp)
        ENDIF
      ENDIF 
    ENDIF
END SUBROUTINE lake_cell_comp

!> Read the latitude and longitude for a cubed-sphere
!! grid from the 'grid' files. For global grids, all
!! six sides are returned.
!!
!! @param[in] res The resolution. Example: '96' for C96.
!! @param[out] grid Array containing the latitude and
!! longitude on the 'supergrid'.  Multiple tiles
!! are concatenated.
!! @author Ning Wang
SUBROUTINE read_cubed_sphere_grid(res, grid) 
    INTEGER, INTENT(IN) :: res
    REAL, INTENT(OUT) :: grid(:,:)

    REAL*8, ALLOCATABLE :: lat(:), lon(:)
    INTEGER :: side_sz, tile_sz, ncid, varid 
    INTEGER :: i, tile_beg, tile_end, stat
    CHARACTER(len=256) :: gridfile_path,gridfile
    CHARACTER(len=1) ich
    CHARACTER(len=4) res_ch

    side_sz = 2*res+1
    tile_sz = side_sz*side_sz
    ALLOCATE(lat(tile_sz), lon(tile_sz))

    IF (tile_req == 0) THEN
      tile_beg = 1; tile_end = 6
    ELSE 
      tile_beg = tile_req; tile_end = tile_req
    ENDIF
    WRITE(res_ch,'(I4)') res
    gridfile_path = "./"
    DO i = tile_beg, tile_end
      WRITE(ich, '(I1)') i
      gridfile = trim(gridfile_path)//"C"//trim(adjustl(res_ch))//"_grid.tile"//ich//".nc"

      PRINT*, 'Open cubed sphere grid spec file ', trim(gridfile)

      stat = nf90_open(trim(gridfile), NF90_NOWRITE, ncid)
      CALL nc_opchk(stat,'nf90_open')

      stat = nf90_inq_varid(ncid,'x',varid)
      CALL nc_opchk(stat,'nf90_inq_lon')
      stat = nf90_get_var(ncid,varid,lon,start=(/1,1/),count=(/side_sz,side_sz/))
      CALL nc_opchk(stat,'nf90_get_var_lon')
      stat = nf90_inq_varid(ncid,'y',varid)
      CALL nc_opchk(stat,'nf90_inq_lat')
      stat = nf90_get_var(ncid,varid,lat,start=(/1,1/),count=(/side_sz,side_sz/))
      CALL nc_opchk(stat,'nf90_get_var_lat')
      stat = nf90_close(ncid)
      CALL nc_opchk(stat,'nf90_close')

      tile_beg = (i-1)*tile_sz+1
      tile_end = i*tile_sz
      grid(tile_beg:tile_end,1) = lat(1:tile_sz)
      grid(tile_beg:tile_end,2) = lon(1:tile_sz)
    END DO
    DEALLOCATE(lat,lon)

END SUBROUTINE read_cubed_sphere_grid

!> Read the latitude and longitude for a regional grid
!! from the 'grid' file. 
!!
!! @param[in] res Resolution of grid.  Example: '96' for C96.
!! @param[out] grid Latitude and longitude on the supergrid.
!! @param[in] halo_depth Lateral halo. Not used.
!! @param[out] res_x Number of grid points in the 'x' direction.
!! @param[out] res_y Number of grid points in the 'y' direction.
!! @author Ning Wang
SUBROUTINE read_cubed_sphere_reg_grid(res, grid, halo_depth, res_x, res_y) 
    INTEGER, INTENT(IN) :: res, halo_depth
    INTEGER, INTENT(OUT) :: res_x, res_y
    REAL, ALLOCATABLE, INTENT(OUT) :: grid(:,:)

    REAL*8, ALLOCATABLE :: lat(:), lon(:)
    INTEGER :: sidex_sz, sidey_sz, tile_sz, ncid, varid, dimid 
    INTEGER :: x_start, y_start
    INTEGER :: nxp, nyp, stat 
    CHARACTER(len=256) :: gridfile_path,gridfile
    CHARACTER(len=1) ich
    CHARACTER(len=4) res_ch
    CHARACTER(len=8) dimname

    WRITE(res_ch,'(I4)') res
    gridfile_path = './' 
    gridfile = trim(gridfile_path)//"C"//trim(adjustl(res_ch))//"_grid.tile7.nc"

    PRINT*, 'Open cubed sphere grid spec file ', trim(gridfile)

    stat = nf90_open(trim(gridfile), NF90_NOWRITE, ncid)
    CALL nc_opchk(stat,'nf90_open')

    stat = nf90_inq_dimid(ncid,'nxp',dimid)
    CALL nc_opchk(stat,'nf90_inq_dimid: nxp')
    stat = nf90_inquire_dimension(ncid,dimid,dimname,len=nxp)
    CALL nc_opchk(stat,'nf90_inquire_dimension: nxp')

    stat = nf90_inq_dimid(ncid,'nyp',dimid)
    CALL nc_opchk(stat,'nf90_inq_dimid: nyp')
    stat = nf90_inquire_dimension(ncid,dimid,dimname,len=nyp)
    CALL nc_opchk(stat,'nf90_inquire_dimension: nyp')

    sidex_sz = nxp
    sidey_sz = nyp
    tile_sz = sidex_sz*sidey_sz
    ALLOCATE(lat(tile_sz), lon(tile_sz))

!    x_start = halo_depth+1; y_start = halo_depth+1
    x_start = 1; y_start = 1
    stat = nf90_inq_varid(ncid,'x',varid)
    CALL nc_opchk(stat,'nf90_inq_lon')
    stat = nf90_get_var(ncid,varid,lon,start=(/x_start,y_start/),count=(/sidex_sz,sidey_sz/))
    CALL nc_opchk(stat,'nf90_get_var_lon')
    stat = nf90_inq_varid(ncid,'y',varid)
    CALL nc_opchk(stat,'nf90_inq_lat')
    stat = nf90_get_var(ncid,varid,lat,start=(/x_start,y_start/),count=(/sidex_sz,sidey_sz/))
    CALL nc_opchk(stat,'nf90_get_var_lat')
    stat = nf90_close(ncid)
    CALL nc_opchk(stat,'nf90_close')

    ALLOCATE(grid(tile_sz,2))
    grid(1:tile_sz,1) = lat(1:tile_sz)
    grid(1:tile_sz,2) = lon(1:tile_sz)

    res_x = sidex_sz/2; res_y = sidey_sz/2
    DEALLOCATE(lat,lon)

END SUBROUTINE read_cubed_sphere_reg_grid

!> Read a high-resolution lake depth dataset, and a corresponding 
!! lake status dataset which provides a status code on the
!! reliability of each lake depth point.
!!
!! @param[in] lakedata_path Path to the lake depth and lake status
!! dataset.
!! @param[out] lake_stat Status code.
!! @param[out] lake_dpth Lake depth.
!! @param[in] nlat 'j' dimension of both datasets.
!! @param[in] nlon 'i' dimension of both datasets.
SUBROUTINE read_lakedata(lakedata_path,lake_stat,lake_dpth,nlat,nlon)
    CHARACTER(len=256), INTENT(IN) :: lakedata_path
    INTEGER*1, INTENT(OUT) :: lake_stat(:)
    INTEGER*2, INTENT(OUT) :: lake_dpth(:)
    INTEGER, INTENT(IN) :: nlat, nlon

    CHARACTER(len=256) lakefile
    INTEGER :: data_sz, i

    data_sz = nlon*nlat

! read in 30 arc seconds lake data 
   ! lake fraction data
    lakefile = trim(lakedata_path) // "GlobalLakeStatus_GLDBv3release.dat"
    IF (lakestatus_srce == "GLDBV2") THEN 
      lakefile = trim(lakedata_path) // "GlobalLakeStatus.dat"
    ENDIF
    IF (lakestatus_srce == "GLDBV3") THEN 
      lakefile = trim(lakedata_path) // "GlobalLakeStatus_GLDBv3release.dat"
    ENDIF
    IF (lakestatus_srce == "MODISP") THEN 
!     lakefile = trim(lakedata_path) // "GlobalLakeStatus_MODIS15s.dat"
      lakefile = trim(lakedata_path) // "GlobalLakeStatus_MODISp.dat"
    ENDIF
    IF (lakestatus_srce == "VIIRS") THEN 
      lakefile = trim(lakedata_path) // "GlobalLakeStatus_VIIRS.dat"
    ENDIF
    OPEN(10, file=lakefile, form='unformatted', access='direct', status='old', recl=data_sz*1)
    READ(10,rec=1) lake_stat
    CLOSE(10)

   ! lake depth data
    lakefile = trim(lakedata_path) // "GlobalLakeDepth_GLDBv3release.dat"
    IF (lakedepth_srce == "GLDBV2") THEN 
      lakefile = trim(lakedata_path) // "GlobalLakeDepth.dat"
    ENDIF
    IF (lakedepth_srce == "GLDBV3") THEN 
      lakefile = trim(lakedata_path) // "GlobalLakeDepth_GLDBv3release.dat"
    ENDIF
    IF (lakedepth_srce == "GLOBATHY") THEN 
      lakefile = trim(lakedata_path) // "GlobalLakeDepth_GLOBathy.dat"
    ENDIF
    OPEN(10, file=lakefile, form='unformatted', access='direct', status='old', recl=data_sz*2)
    READ(10,rec=1) lake_dpth
    CLOSE(10)

END SUBROUTINE read_lakedata
   
!> Write lake depth and fraction to an existing model orography file.
!! Also, perform some quality control checks on the lake data.
!! This routine is used for non-regional grids.
!!
!! @param[in] cs_res Resolution. Example: '96' for C96.
!! @param[in] cs_lakestat Lake fraction.
!! @param[in] cs_lakedpth Lake depth.
!! @author Ning Wang
SUBROUTINE write_lakedata_to_orodata(cs_res, cs_lakestat, cs_lakedpth) 
    USE netcdf 
    INTEGER, INTENT(IN) :: cs_res
    REAL, INTENT(IN) :: cs_lakestat(:)
    REAL, INTENT(IN) :: cs_lakedpth(:)
   
    INTEGER :: tile_sz, tile_num
    INTEGER :: stat, ncid, x_dimid, y_dimid, varid, dimids(2)
    INTEGER :: lake_frac_id, lake_depth_id
    INTEGER :: land_frac_id, slmsk_id, inland_id, geolon_id, geolat_id
    CHARACTER(len=256) :: filename,string
    CHARACTER(len=1) :: ich
    CHARACTER(len=4) res_ch
    REAL :: lake_frac(cs_res*cs_res),lake_depth(cs_res*cs_res)
    REAL :: geolon(cs_res*cs_res),geolat(cs_res*cs_res)
    REAL :: land_frac(cs_res*cs_res),slmsk(cs_res*cs_res),inland(cs_res*cs_res)
    real, parameter :: epsil=1.e-6   ! numerical min for lake_frac/land_frac
    real            :: land_cutoff=1.e-4 ! land_frac=0 if it is < land_cutoff

    INTEGER :: i, j

    tile_sz = cs_res*cs_res

    WRITE(res_ch,'(I4)') cs_res
    DO tile_num = tile_beg, tile_end
      WRITE(ich, '(I1)') tile_num
!      filename = "C" // trim(adjustl(res_ch)) // "_oro_data.tile" // ich // ".nc" 
!      filename = "oro_data.tile" // ich // ".nc" 
      filename = "oro.C" // trim(adjustl(res_ch)) // ".tile" // ich // ".nc" 
      print *,'Read, update, and write ',trim(filename)
      stat = nf90_open(filename, NF90_WRITE, ncid)
      CALL nc_opchk(stat, "nf90_open oro_data.nc")
      stat = nf90_inq_dimid(ncid, "lon", x_dimid)
      CALL nc_opchk(stat, "nf90_inq_dim: x")
      stat = nf90_inq_dimid(ncid, "lat", y_dimid)
      CALL nc_opchk(stat, "nf90_inq_dim: y")
!      dimids = (/ y_dimid, x_dimid /)
! original orodata netcdf file uses (y, x) order, so we made change to match it. 
      dimids = (/ x_dimid, y_dimid /)
      stat = nf90_inq_varid(ncid, "land_frac", land_frac_id)
      CALL nc_opchk(stat, "nf90_inq_varid: land_frac")
      stat = nf90_inq_varid(ncid, "slmsk", slmsk_id)
      CALL nc_opchk(stat, "nf90_inq_varid: slmsk")
! define 2 new variables 
      stat = nf90_redef(ncid)
      CALL nc_opchk(stat, "nf90_redef")
      stat = nf90_def_var(ncid,"lake_frac",NF90_FLOAT,dimids,lake_frac_id)
      CALL nc_opchk(stat, "nf90_def_var: lake_frac")
#ifdef ADD_ATT_FOR_NEW_VAR
      stat = nf90_put_att(ncid, lake_frac_id,'coordinates','geolon geolat')
      CALL nc_opchk(stat, "nf90_put_att: lake_frac:coordinates") 
      stat = nf90_put_att(ncid, lake_frac_id,'long_name','lake fraction')
      CALL nc_opchk(stat, "nf90_put_att: lake_frac:long_name") 
      stat = nf90_put_att(ncid, lake_frac_id,'unit','fraction')
      CALL nc_opchk(stat, "nf90_put_att: lake_frac:unit") 
      IF (lakestatus_srce == "GLDBV3") THEN 
        write(string,'(a,es8.1)') 'based on GLDBv3 (Choulga et al. 2019); missing lakes &
        added based on land_frac in this dataset; lake_frac cutoff:',lake_cutoff
      ELSE IF (lakestatus_srce == "GLDBV2") THEN 
        write(string,'(a,es8.1)') 'based on GLDBv2 (Choulga et al. 2014); missing lakes &
        added based on land_frac in this dataset; lake_frac cutoff:',lake_cutoff
      ELSE IF (lakestatus_srce == "MODISP") THEN
        write(string,'(a,es8.1)') 'based on MODIS (2011-2015) product updated with two &
        Landsat products: the JRC water product (2016-2020) and the GLC-FCS30 (2020); &
        the source data set was created by Chengquan Huang of UMD; &  
        lake_frac cutoff:',lake_cutoff
      ELSE IF (lakestatus_srce == "VIIRS") THEN
        write(string,'(a,es8.1)') 'based on multi-year VIIRS global surface type & 
        classification map (2012-2019); the source data set was created by &
        Chengquan Huang of UMD and Michael Barlage of NOAA; &
        lake_frac cutoff:',lake_cutoff
      ENDIF
      stat = nf90_put_att(ncid, lake_frac_id,'description',trim(string))
      CALL nc_opchk(stat, "nf90_put_att: lake_frac:description") 
#endif
      stat = nf90_def_var(ncid,"lake_depth",NF90_FLOAT,dimids,lake_depth_id)
      CALL nc_opchk(stat, "nf90_def_var: lake_depth")
#ifdef ADD_ATT_FOR_NEW_VAR
      stat = nf90_put_att(ncid, lake_depth_id,'coordinates','geolon geolat')
      CALL nc_opchk(stat, "nf90_put_att: lake_depth:coordinates") 
      stat = nf90_put_att(ncid, lake_depth_id,'long_name','lake depth')
      CALL nc_opchk(stat, "nf90_put_att: lake_depth:long_name") 
      stat = nf90_put_att(ncid, lake_depth_id,'unit','meter')
      CALL nc_opchk(stat, "nf90_put_att: lake_depth:long_name") 
      IF (lakedepth_srce == "GLDBV3") THEN 
        stat = nf90_put_att(ncid, lake_depth_id,'description', &
        'based on GLDBv3 (Choulga et al. 2019); missing depth set to 10m &
        (except to 211m in Caspian Sea); spurious large pos. depths are left unchanged.')
      ELSE IF (lakedepth_srce == "GLDBV2") THEN 
        stat = nf90_put_att(ncid, lake_depth_id,'description', &
        'based on GLDBv2 (Choulga et al. 2014); missing depth set to 10m &
        (except to 211m in Caspian Sea); spurious large pos. depths are left unchanged.')
      ELSE IF (lakedepth_srce == "GLOBATHY") THEN
        stat = nf90_put_att(ncid, lake_depth_id,'description', &
        'based on GLOBathy data resampled and projected to the MODIS domain.')
      ENDIF
      CALL nc_opchk(stat, "nf90_put_att: lake_depth:description") 
#endif

      write(string,'(a,es8.1)') 'land_frac and lake_frac are adjusted such that &
         their sum is 1 at points where inland=1; land_frac cutoff is',land_cutoff
      stat = nf90_put_att(ncid, land_frac_id,'description',trim(string))
      CALL nc_opchk(stat, "nf90_put_att: land_frac:description") 

      write(string,'(a)') 'slmsk = nint(land_frac)'
      stat = nf90_put_att(ncid, slmsk_id,'description',trim(string))
      CALL nc_opchk(stat, "nf90_put_att: slmsk:description") 

      stat = nf90_enddef(ncid) 
      CALL nc_opchk(stat, "nf90_enddef")

! read in geolon and geolat and 2 variables from orog data file
      stat = nf90_inq_varid(ncid, "geolon", geolon_id)
      CALL nc_opchk(stat, "nf90_inq_varid: geolon")
      stat = nf90_inq_varid(ncid, "geolat", geolat_id)
      CALL nc_opchk(stat, "nf90_inq_varid: geolat")
      stat = nf90_inq_varid(ncid, "land_frac", land_frac_id)
      CALL nc_opchk(stat, "nf90_inq_varid: land_frac")
      stat = nf90_inq_varid(ncid, "slmsk", slmsk_id)
      CALL nc_opchk(stat, "nf90_inq_varid: slmsk")
      stat = nf90_inq_varid(ncid, "inland", inland_id)
      CALL nc_opchk(stat, "nf90_inq_varid: inland")
      stat = nf90_get_var(ncid, geolon_id, geolon, &
             start = (/ 1, 1 /), count = (/ cs_res, cs_res /) )
      CALL nc_opchk(stat, "nf90_get_var: geolon")
      stat = nf90_get_var(ncid, geolat_id, geolat, &
             start = (/ 1, 1 /), count = (/ cs_res, cs_res /) )
      CALL nc_opchk(stat, "nf90_get_var: geolat")
      stat = nf90_get_var(ncid, land_frac_id, land_frac, &
             start = (/ 1, 1 /), count = (/ cs_res, cs_res /) )
      CALL nc_opchk(stat, "nf90_get_var: land_frac")
      stat = nf90_get_var(ncid, slmsk_id, slmsk, &
             start = (/ 1, 1 /), count = (/ cs_res, cs_res /) )
      CALL nc_opchk(stat, "nf90_get_var: slmsk")
      stat = nf90_get_var(ncid, inland_id, inland, &
             start = (/ 1, 1 /), count = (/ cs_res, cs_res /) )
      CALL nc_opchk(stat, "nf90_get_var: inland")

      lake_frac (:) = cs_lakestat ((tile_num-1)*tile_sz+1:tile_num*tile_sz)
      lake_depth(:) = cs_lakedepth((tile_num-1)*tile_sz+1:tile_num*tile_sz)

! include Caspian Sea and Aral Sea if GLDB data set is used, and 
! exclude lakes in the coastal areas of Antarctica if MODIS data set is used  
      CALL include_exclude_lakes(lake_frac,land_frac,lake_depth,geolat,geolon,tile_num)

! epsil is "numerical" nonzero min for lake_frac/land_frac
! lake_cutoff/land_cutoff is practical min for lake_frac/land_frac
      IF (min(lake_cutoff,land_cutoff) < epsil) then
        PRINT *,'lake_cutoff/land_cutoff cannot be smaller than epsil, reset...'
        lake_cutoff=max(epsil,lake_cutoff)
        land_cutoff=max(epsil,land_cutoff)
      ENDIF

! adjust land_frac and lake_frac, and make sure land_frac+lake_frac=1 at inland points
      DO i = 1, tile_sz
        if (land_frac(i)<   land_cutoff) land_frac(i)=0.
        if (land_frac(i)>1.-land_cutoff) land_frac(i)=1.

        if (inland(i) /= 1.) then 
          lake_frac(i) = 0.
        endif

        if (lake_frac(i) < lake_cutoff) lake_frac(i)=0.
        if (lake_frac(i) > 1.-epsil) lake_frac(i)=1.
      ENDDO

! finalize land_frac/slmsk based on modified lake_frac
      DO i = 1, tile_sz
        if (inland(i) == 1.) then
          land_frac(i) = 1. - lake_frac(i)
        end if

        if (lake_frac(i) < lake_cutoff) then
          lake_depth(i)=0.
        elseif (lake_frac(i) >= lake_cutoff .and. lake_depth(i)==0.) then
          lake_depth(i)=10.
        end if
        slmsk(i) = nint(land_frac(i))
      ENDDO

! write 2 new variables      
      stat = nf90_put_var(ncid, lake_frac_id, lake_frac, &
             start = (/ 1, 1 /), count = (/ cs_res, cs_res /) )
      CALL nc_opchk(stat, "nf90_put_var: lake_frac") 

      stat = nf90_put_var(ncid, lake_depth_id, lake_depth, &
             start = (/ 1, 1 /), count = (/ cs_res, cs_res /) )
      CALL nc_opchk(stat, "nf90_put_var: lake_depth") 

! write back 2 modified variables: land_frac and slmsk
      stat = nf90_put_var(ncid, land_frac_id, land_frac, &
             start = (/ 1, 1 /), count = (/ cs_res, cs_res /) )
      CALL nc_opchk(stat, "nf90_put_var: land_frac") 

      stat = nf90_put_var(ncid, slmsk_id, slmsk, &
             start = (/ 1, 1 /), count = (/ cs_res, cs_res /) )
      CALL nc_opchk(stat, "nf90_put_var: slmsk") 
    
      stat = nf90_close(ncid)
      CALL nc_opchk(stat, "nf90_close") 
    ENDDO
  
END SUBROUTINE write_lakedata_to_orodata

!> Write lake depth and fraction to an existing model orography file.
!! Also, perform some quality control checks on the lake data.
!! This routine is used for regional grids.
!!
!! @param[in] cs_res Resolution. Example: '96' for C96.
!! @param[in] cs_lakestat Lake fraction.
!! @param[in] cs_lakedpth Lake depth.
!! @param[in] tile_x_dim 'x' dimension of the model grid.
!! @param[in] tile_y_dim 'y' dimension of the model grid.
!! @author Ning Wang
SUBROUTINE write_reg_lakedata_to_orodata(cs_res, tile_x_dim, tile_y_dim, cs_lakestat, cs_lakedpth) 
    USE netcdf 
    INTEGER, INTENT(IN) :: cs_res, tile_x_dim, tile_y_dim
    REAL, INTENT(IN) :: cs_lakestat(:)
    REAL, INTENT(IN) :: cs_lakedpth(:)
   
    INTEGER :: tile_sz, tile_num
    INTEGER :: stat, ncid, x_dimid, y_dimid, varid, dimids(2)
    INTEGER :: lake_frac_id, lake_depth_id
    INTEGER :: land_frac_id, slmsk_id, geolon_id, geolat_id, inland_id
    CHARACTER(len=256) :: filename,string
    CHARACTER(len=1) :: ich
    CHARACTER(len=4) res_ch

    REAL, ALLOCATABLE :: lake_frac(:), lake_depth(:), geolon(:), geolat(:)
    REAL, ALLOCATABLE :: land_frac(:), slmsk(:), inland(:)

    real, parameter :: epsil=1.e-6   ! numerical min for lake_frac/land_frac
    real            :: land_cutoff=1.e-6 ! land_frac=0 if it is < land_cutoff

    INTEGER :: i, j, var_id

!    include "netcdf.inc"
    tile_sz = tile_x_dim*tile_y_dim

    ALLOCATE(lake_frac(tile_sz), lake_depth(tile_sz))
    ALLOCATE(geolon(tile_sz), geolat(tile_sz))
    ALLOCATE(land_frac(tile_sz), slmsk(tile_sz), inland(tile_sz))

    WRITE(res_ch,'(I4)') cs_res
    tile_num  = 7
    WRITE(ich, '(I1)') tile_num
!    filename = "C" // trim(adjustl(res_ch)) // "_oro_data.tile" // ich // ".halo4.nc" 
    filename = "oro.C" // trim(adjustl(res_ch)) // ".tile" // ich // ".nc"
    PRINT*, 'Open and write regional orography data netCDF file ', trim(filename)
    stat = nf90_open(filename, NF90_WRITE, ncid)
    CALL nc_opchk(stat, "nf90_open oro_data.nc")
    stat = nf90_inq_dimid(ncid, "lon", x_dimid)
    CALL nc_opchk(stat, "nf90_inq_dim: x")
    stat = nf90_inq_dimid(ncid, "lat", y_dimid)
    CALL nc_opchk(stat, "nf90_inq_dim: y")
    dimids = (/ x_dimid, y_dimid /)

    stat = nf90_inq_varid(ncid, "land_frac", land_frac_id)
    CALL nc_opchk(stat, "nf90_inq_varid: land_frac")
    stat = nf90_inq_varid(ncid, "slmsk", slmsk_id)
    CALL nc_opchk(stat, "nf90_inq_varid: slmsk")

! define 2 new variables, lake_frac and lake_depth 
    stat = nf90_redef(ncid)
    CALL nc_opchk(stat, "nf90_redef")

    IF (nf90_inq_varid(ncid, "lake_frac",lake_frac_id) /= 0) THEN
      stat = nf90_def_var(ncid,"lake_frac",NF90_FLOAT,dimids,lake_frac_id)
      CALL nc_opchk(stat, "nf90_def_var: lake_frac")
#ifdef ADD_ATT_FOR_NEW_VAR
      stat = nf90_put_att(ncid, lake_frac_id,'coordinates','geolon geolat')
      CALL nc_opchk(stat, "nf90_put_att: lake_frac:coordinates") 
      stat = nf90_put_att(ncid, lake_frac_id,'long_name','lake fraction')
      CALL nc_opchk(stat, "nf90_put_att: lake_frac:long_name") 
      stat = nf90_put_att(ncid, lake_frac_id,'unit','fraction')
      CALL nc_opchk(stat, "nf90_put_att: lake_frac:unit") 
      IF (lakestatus_srce == "GLDBV3") THEN 
        write(string,'(a,es8.1)') 'based on GLDBv3 (Choulga et al. 2019); missing lakes &
        added based on land_frac in this dataset; lake_frac cutoff:',lake_cutoff
      ELSE IF (lakestatus_srce == "GLDBV2") THEN 
        write(string,'(a,es8.1)') 'based on GLDBv2 (Choulga et al. 2014); missing lakes &
        added based on land_frac in this dataset; lake_frac cutoff:',lake_cutoff
      ELSE IF (lakestatus_srce == "MODISP") THEN
        write(string,'(a,es8.1)') 'based on MODIS (2011-2015) product updated with two &
        Landsat products: the JRC water product (2016-2020) and the GLC-FCS30 (2020); &
        the source data set was created by Chengquan Huang of UMD; &  
        lake_frac cutoff:',lake_cutoff
      ELSE IF (lakestatus_srce == "VIIRS") THEN
        write(string,'(a,es8.1)') 'based on multi-year VIIRS global surface type & 
        classification map (2012-2019); the source data set was created by &
        Chengquan Huang of UMD and Michael Barlage of NOAA; &
        lake_frac cutoff:',lake_cutoff
      ENDIF
      stat = nf90_put_att(ncid, lake_frac_id,'description',trim(string))
      CALL nc_opchk(stat, "nf90_put_att: lake_frac:description") 
#endif
    ENDIF
    IF (nf90_inq_varid(ncid, "lake_depth",lake_depth_id) /= 0) THEN
      stat = nf90_def_var(ncid,"lake_depth",NF90_FLOAT,dimids,lake_depth_id)
      CALL nc_opchk(stat, "nf90_def_var: lake_depth")
#ifdef ADD_ATT_FOR_NEW_VAR
      stat = nf90_put_att(ncid, lake_depth_id,'coordinates','geolon geolat')
      CALL nc_opchk(stat, "nf90_put_att: lake_depth:coordinates") 
      stat = nf90_put_att(ncid, lake_depth_id,'long_name','lake depth')
      CALL nc_opchk(stat, "nf90_put_att: lake_depth:long_name") 
      stat = nf90_put_att(ncid, lake_depth_id,'unit','meter')
      CALL nc_opchk(stat, "nf90_put_att: lake_depth:long_name") 
      IF (lakedepth_srce == "GLDBV3") THEN 
        stat = nf90_put_att(ncid, lake_depth_id,'description', &
        'based on GLDBv3 (Choulga et al. 2019); missing depth set to 10m &
        (except to 211m in Caspian Sea); spurious large pos. depths are left unchanged.')
      ELSE IF (lakedepth_srce == "GLDBV2") THEN 
        stat = nf90_put_att(ncid, lake_depth_id,'description', &
        'based on GLDBv2 (Choulga et al. 2014); missing depth set to 10m &
        (except to 211m in Caspian Sea); spurious large pos. depths are left unchanged.')
      ELSE IF (lakedepth_srce == "GLOBATHY") THEN
        stat = nf90_put_att(ncid, lake_depth_id,'description', &
        'based on GLOBathy data resampled and projected to the MODIS domain.')
      ENDIF
      CALL nc_opchk(stat, "nf90_put_att: lake_depth:description") 
    ENDIF
#endif
    write(string,'(a,es8.1)') 'land_frac and lake_frac are adjusted such that &
      their sum is 1 at points where inland=1; land_frac cutoff is',land_cutoff
    stat = nf90_put_att(ncid, land_frac_id,'description',trim(string))
    CALL nc_opchk(stat, "nf90_put_att: land_frac:description") 
    write(string,'(a)') 'slmsk = nint(land_frac)'
    stat = nf90_put_att(ncid, slmsk_id,'description',trim(string))
    CALL nc_opchk(stat, "nf90_put_att: slmsk:description") 

    stat = nf90_enddef(ncid) 
    CALL nc_opchk(stat, "nf90_enddef")

! read in geolon and geolat and 2 variables from orog data file
    stat = nf90_inq_varid(ncid, "geolon", geolon_id)
    CALL nc_opchk(stat, "nf90_inq_varid: geolon")
    stat = nf90_inq_varid(ncid, "geolat", geolat_id)
    CALL nc_opchk(stat, "nf90_inq_varid: geolat")
    stat = nf90_inq_varid(ncid, "land_frac", land_frac_id)
    CALL nc_opchk(stat, "nf90_inq_varid: land_frac")
    stat = nf90_inq_varid(ncid, "slmsk", slmsk_id)
    CALL nc_opchk(stat, "nf90_inq_varid: slmsk")
    stat = nf90_inq_varid(ncid, "inland", inland_id)
    CALL nc_opchk(stat, "nf90_inq_varid: inland")

    stat = nf90_get_var(ncid, geolon_id, geolon, &
           start = (/ 1, 1 /), count = (/ tile_x_dim, tile_y_dim /) )
    CALL nc_opchk(stat, "nf90_get_var: geolon")
    stat = nf90_get_var(ncid, geolat_id, geolat, &
           start = (/ 1, 1 /), count = (/ tile_x_dim, tile_y_dim /) )
    CALL nc_opchk(stat, "nf90_get_var: geolat")
    stat = nf90_get_var(ncid, land_frac_id, land_frac, &
           start = (/ 1, 1 /), count = (/ tile_x_dim, tile_y_dim /) )
    CALL nc_opchk(stat, "nf90_get_var: land_frac")
    stat = nf90_get_var(ncid, slmsk_id, slmsk, &
           start = (/ 1, 1 /), count = (/ tile_x_dim, tile_y_dim /) )
    CALL nc_opchk(stat, "nf90_get_var: slmsk")
    stat = nf90_get_var(ncid, inland_id, inland, &
           start = (/ 1, 1 /), count = (/ tile_x_dim, tile_y_dim /) )
    CALL nc_opchk(stat, "nf90_get_var: inland")

    tile_num = 1
    lake_frac(:)  = cs_lakestat((tile_num-1)*tile_sz+1:tile_num*tile_sz)
    lake_depth(:) = cs_lakedepth((tile_num-1)*tile_sz+1:tile_num*tile_sz)

! include Caspian Sea and Aral Sea if GLDB data set is used, and 
! exclude lakes in the coastal areas of Antarctica if MODIS data set is used  
    CALL include_exclude_lakes(lake_frac,land_frac,lake_depth,geolat,geolon,tile_num)

! epsil is "numerical" nonzero min for lake_frac/land_frac
! lake_cutoff/land_cutoff is practical min for lake_frac/land_frac
    IF (min(lake_cutoff,land_cutoff) < epsil) then
        print *,'lake_cutoff/land_cutoff cannot be smaller than epsil, reset...'
        lake_cutoff=max(epsil,lake_cutoff)
        land_cutoff=max(epsil,land_cutoff)
    ENDIF 

! adjust land_frac and lake_frac, and make sure land_frac+lake_frac=1 at inland points
    DO i = 1, tile_sz
      if (land_frac(i)<   land_cutoff) land_frac(i)=0.
      if (land_frac(i)>1.-land_cutoff) land_frac(i)=1.

      if (inland(i) /= 1.) then
        lake_frac(i) = 0.
      endif

      if (lake_frac(i) < lake_cutoff) lake_frac(i)=0.
      if (lake_frac(i) > 1.-epsil) lake_frac(i)=1.
    ENDDO

    DO i = 1, tile_sz
      if (inland(i) == 1.) then
        land_frac(i) = 1. - lake_frac(i)
      end if

      if (lake_frac(i) < lake_cutoff) then
        lake_depth(i)=0.
      elseif (lake_frac(i) >= lake_cutoff .and. lake_depth(i)==0.) then
        lake_depth(i)=10.
      end if
      slmsk(i) = nint(land_frac(i))
    ENDDO

! write 2 new variables      
    stat = nf90_put_var(ncid, lake_frac_id, lake_frac, &
           start = (/ 1, 1 /), count = (/ tile_x_dim, tile_y_dim /) )
    CALL nc_opchk(stat, "nf90_put_var: lake_frac") 

    stat = nf90_put_var(ncid, lake_depth_id, lake_depth, &
           start = (/ 1, 1 /), count = (/ tile_x_dim, tile_y_dim /) )
    CALL nc_opchk(stat, "nf90_put_var: lake_depth") 

! write back 2 modified variables: land_frac and slmsk
    stat = nf90_put_var(ncid, land_frac_id, land_frac, &
           start = (/ 1, 1 /), count = (/ tile_x_dim, tile_y_dim /) )
    CALL nc_opchk(stat, "nf90_put_var: land_frac") 

    stat = nf90_put_var(ncid, slmsk_id, slmsk, &
           start = (/ 1, 1 /), count = (/ tile_x_dim, tile_y_dim /) )
    CALL nc_opchk(stat, "nf90_put_var: slmsk") 
    
    stat = nf90_close(ncid)
    CALL nc_opchk(stat, "nf90_close") 
  
END SUBROUTINE write_reg_lakedata_to_orodata

!> Include Caspian Sea and Aral Sea if GLDB dataset is used, and 
!! exclude lakes in the coastal areas of Antarctica if MODIS dataset is used.
!!
!! @param[inout] lake_frac lake fraction array of the given tile. 
!! @param[inout] lake_depth lake depth array of the given tile. 
!! @param[in] land_frac land fraction array of the given tile.
!! @param[in] geolat latitude array of the given tile.
!! @param[in] geolon longitude array of the given tile.
!! @param[in] tile_num tile number of the given tile.
!! @author Ning Wang
SUBROUTINE include_exclude_lakes(lake_frac,land_frac,lake_depth,geolat,geolon,tile_num)
    REAL, INTENT(INOUT) :: lake_frac(cs_res*cs_res), lake_depth(cs_res*cs_res)
    REAL, INTENT(IN) :: land_frac(cs_res*cs_res)
    REAL, INTENT(IN) :: geolat(cs_res*cs_res), geolon(cs_res*cs_res)
    INTEGER, INTENT(IN) :: tile_num
    
    INTEGER :: tile_sz

    tile_sz = cs_res*cs_res
! add Caspian Sea and Aral Sea 
    IF (tile_num == 2 .OR. tile_num == 3 .OR. tile_num == 7) THEN
      IF (lakedepth_srce == "GLDBV3" .OR. lakedepth_srce == "GLDBV2") THEN 
        IF (lakestatus_srce == "GLDBV3" .OR. lakestatus_srce == "GLDBV2") THEN
          DO i = 1, tile_sz
            IF (land_frac(i) < 0.9 .AND. lake_frac(i) < 0.1) THEN
              IF (geolat(i) > 35.0 .AND. geolat(i) <= 50.0 .AND. &
                  geolon(i) > 45.0 .AND. geolon(i) <= 55.0) THEN
                lake_frac(i) = 1.-land_frac(i)
                lake_depth(i) = 211.0
              ENDIF
              IF (geolat(i) > 35.0 .AND. geolat(i) <= 50.0 .AND. &
                  geolon(i) > 57.0 .AND. geolon(i) <= 63.0) THEN
                lake_frac(i) = 1.-land_frac(i)
                lake_depth(i) = 10.0
              ENDIF
            ENDIF
          ENDDO
        ENDIF
        IF (lakestatus_srce == "MODISP" .OR. lakestatus_srce == "VIIRS") THEN
          DO i = 1, tile_sz
            IF (land_frac(i) < 0.9) THEN
              IF (geolat(i) > 35.0 .AND. geolat(i) <= 50.0 .AND. &
                  geolon(i) > 45.0 .AND. geolon(i) <= 55.0) THEN
                lake_depth(i) = 211.0
              ENDIF
              IF (geolat(i) > 35.0 .AND. geolat(i) <= 50.0 .AND. &
                  geolon(i) > 57.0 .AND. geolon(i) <= 63.0) THEN
                lake_depth(i) = 10.0
              ENDIF
            ENDIF
          ENDDO
        ENDIF
      ENDIF
    ENDIF
! remove lakes below 60 deg south 
    IF (tile_num == 6) THEN
      IF (lakestatus_srce == "MODISP" .OR. lakestatus_srce == "VIIRS") THEN 
        DO i = 1, tile_sz
          IF (geolat(i) < -60.0) THEN
            lake_frac(i) = 0.0
            lake_depth(i) = 0.0
          ENDIF
        ENDDO
      ENDIF
    ENDIF
     
END SUBROUTINE include_exclude_lakes

!> Check NetCDF error code
!!
!! @param[in] stat Error code.
!! @param[in] opname NetCDF operation that failed.
!! @author Ning Wang
SUBROUTINE nc_opchk(stat,opname)
   USE netcdf
   IMPLICIT NONE
   INTEGER stat
   CHARACTER(len=*) opname
   CHARACTER(64) msg

   IF (stat .NE.0)  THEN
     msg = trim(opname) // ' Error, status code and message:'
     PRINT*,trim(msg), stat, nf90_strerror(stat)
     STOP -1
   END IF

END SUBROUTINE nc_opchk

END PROGRAM lake_frac
