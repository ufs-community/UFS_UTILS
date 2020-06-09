!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
! This program computes lake fraction and depth numbers for FV3 cubed sphere 
! grid cells, from a high resolution lat/lon data set.
! 
! Ning Wang, July 2018: Original version
! 
!   Shan Sun, Aug. 2018: Added Caspian Sea and Aral Sea to the lake fraction 
!                        and lake depth fields.
!   Shan Sun, Dec. 2018: Added round up and round down with respect to a 
!                        numerical minimum value and a cut-off value, for lake 
!                        fraction number. 
!   Ning Wang, Apr. 2019: Extended the program to process the same lake data 
!                         for FV3 stand-alone regional (SAR) model.
!     
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!#define DIAG_N_VERBOSE
#define WRITE_TO_ORODATA
#define ADD_ATT_FOR_NEW_VAR
PROGRAM lake_frac
    USE netcdf
    IMPLICIT NONE

    CHARACTER(len=256) :: sfcdata_path
    INTEGER :: cs_res, ncsmp, ncscp, i

    INTEGER*1, ALLOCATABLE :: lakestatus(:)
    INTEGER*2, ALLOCATABLE :: lakedepth(:)
    REAL, ALLOCATABLE :: cs_grid(:,:)
    REAL, ALLOCATABLE :: cs_lakestatus(:), cs_lakedepth(:)
    REAL, ALLOCATABLE :: src_grid_lon(:), src_grid_lat(:)

    INTEGER :: tile_req, tile_beg, tile_end

    INTEGER, PARAMETER :: nlat = 21600, nlon = 43200
    REAL, PARAMETER :: d2r = acos(-1.0) / 180.0
    REAL, PARAMETER :: r2d = 180.0 /acos(-1.0) 
    REAL, PARAMETER :: pi = acos(-1.0) 
    REAL*8, PARAMETER :: gppdeg = 119.99444445
    REAL*8, PARAMETER :: delta = 1.0 / 119.99444445

    CHARACTER(len=32) :: arg
    INTEGER :: stat
    
    CALL getarg(0, arg) ! get the program name
    IF (iargc() /= 2) THEN
      PRINT*, 'Usage: ', trim(arg), & 
       ' [tile_num (0 - all tiles)] [Resolution (48,96, ...)]'  
      STOP
    ENDIF
    CALL getarg(1, arg)
    READ(arg,*,iostat=stat) tile_req
    CALL getarg(2, arg)
    READ(arg,*,iostat=stat) cs_res

    IF (tile_req == 0) THEN
      tile_beg = 1; tile_end = 6
      PRINT*, 'Process tile 1 - 6 at resolution C',cs_res
    ELSE
      tile_beg = tile_req; tile_end = tile_req
      PRINT*, 'Process tile',tile_req, 'at resolution C',cs_res
    ENDIF

!    cs_res = 96 !    tile_beg = 3; tile_end = 3

    ! read in grid spec data for each tile and concatenate them together
    ncsmp = (2*cs_res+1)*(2*cs_res+1)*6   

    ALLOCATE(cs_grid(ncsmp, 2))
    PRINT*, 'Read in cubed sphere grid information ... ',ncsmp,'pairs of lat/lons'
    CALL  read_cubed_sphere_grid(cs_res, cs_grid) 
!    cs_grid = cs_grid * d2r

    ! compute source grid 
    ALLOCATE(src_grid_lon(nlon), src_grid_lat(nlat))
    DO i = 1, nlon
      src_grid_lon(i) = -180.0 + delta * (i-1)
    ENDDO
    DO i = 1, nlat
      src_grid_lat(i) = 90.0 - delta * (i-1) 
    ENDDO

    ! read in lake data file
    sfcdata_path = '/scratch1/BMC/gsd-fv3-dev/data_others/lake_rawdata/'
    ALLOCATE(lakestatus(nlon*nlat),lakedepth(nlon*nlat)) 
    PRINT*, 'Read in lake data file ...'
    CALL read_lakedata(sfcdata_path,lakestatus,lakedepth,nlat,nlon)

    ! calculate fraction numbers for all cs-cells
    ncscp = cs_res*cs_res*6
    ALLOCATE(cs_lakestatus(ncscp))
    ALLOCATE(cs_lakedepth(ncscp))

    PRINT*, 'Calculate lake fraction and average depth for cubed-sphere cells ...'
    CALL cal_lake_frac_depth(cs_res,lakestatus,cs_lakestatus,lakedepth,cs_lakedepth)

    ! write lake status (in terms of fraction) and lake depth(average)
    ! to a netcdf file
    PRINT*, 'Write netcdf file for lake fraction on cubed sphere grid cells ...'
#ifdef WRITE_TO_ORODATA
    CALL write_lakedata_to_orodata(cs_res, cs_lakestatus, cs_lakedepth) 
#else
    CALL write_lakefrac_netcdf(cs_res, cs_lakestatus, cs_lakedepth) 
#endif

    DEALLOCATE(cs_lakestatus,cs_lakedepth)
    DEALLOCATE(cs_grid)
    DEALLOCATE(lakestatus,lakedepth)
    STOP
CONTAINS

SUBROUTINE cal_lake_frac_depth(res,lakestat,cs_lakestat,lakedpth,cs_lakedpth)
    INTEGER, INTENT(IN) :: res
    INTEGER*1, INTENT(IN) :: lakestat(:)
    INTEGER*2, INTENT(IN) :: lakedpth(:)
    REAL, INTENT(OUT) :: cs_lakestat(:), cs_lakedpth(:)

    REAL*8 lolf(2), lort(2), uplf(2), uprt(2), sd_ltmn(4), sd_ltmx(4)
    REAL*8 :: v(2,4), p(2)
    REAL :: latmin1, latmax1
    REAL :: latmin, latmax, lonmin, lonmax, lontmp, lat_sz_max, lon_sz_max
    REAL :: lonmin_180pm, lonmax_180pm
    INTEGER :: tile_num, i, j, gp, side_sz, row, col, cs_grid_idx, cs_data_idx
    INTEGER :: stride_lat, stride_lon
    INTEGER :: src_grid_lat_beg,src_grid_lat_end,src_grid_lon_beg,src_grid_lon_end
    INTEGER :: src_grid_lon_beg1,src_grid_lon_end1,src_grid_lon_beg2,src_grid_lon_end2
    INTEGER :: grid_ct, lake_ct, co_gc, tmp

    REAL*8 :: lake_dpth_sum
    LOGICAL :: two_section, enclosure_cnvx 

    stride_lat = 1
    side_sz = 2*res+1

    lat_sz_max = 0.0
    lon_sz_max = 0.0

    cs_lakestat = 0.0

    DO tile_num = tile_beg, tile_end
      row = 2 + side_sz*(tile_num-1); col = 2
      DO gp = 1, res*res
        two_section = .false.
        cs_grid_idx = (row-1)*side_sz+col
        cs_data_idx = (tile_num-1)*res*res+gp
        IF (abs(cs_grid(cs_grid_idx,1)) > 80.0 ) THEN !ignore lakes in very high latitude
          cs_lakestat(cs_data_idx) = 0.0
          cs_lakedpth(cs_data_idx) = 0.0
          ! move to next cs cell
          col = col + 2
          IF (col > side_sz) THEN
            col = 2
            row = row + 2
          ENDIF
          CYCLE
        ENDIF
        ! get the four corners of the cs cell 
        lolf(1) = cs_grid(cs_grid_idx-side_sz-1, 1)
        lolf(2) = cs_grid(cs_grid_idx-side_sz-1, 2)
        IF (lolf(2) > 180.0) lolf(2) = lolf(2) - 360.0
        lort(1) = cs_grid(cs_grid_idx-side_sz+1, 1)
        lort(2) = cs_grid(cs_grid_idx-side_sz+1, 2)
        IF (lort(2) > 180.0) lort(2) = lort(2) - 360.0
        uplf(1) = cs_grid(cs_grid_idx+side_sz-1,1)
        uplf(2) = cs_grid(cs_grid_idx+side_sz-1,2)
        IF (uplf(2) > 180.0) uplf(2) = uplf(2) - 360.0
        uprt(1) = cs_grid(cs_grid_idx+side_sz+1,1)
        uprt(2) = cs_grid(cs_grid_idx+side_sz+1,2)

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
          print*,'switch!!!'
          tmp = src_grid_lat_beg
          src_grid_lat_beg = src_grid_lat_end
          src_grid_lat_end = tmp
        ENDIF
        IF (src_grid_lon_beg > src_grid_lon_end) THEN
          print*,'switch!!!'
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
          gp, res*res, cs_grid(cs_grid_idx,1),cs_grid(cs_grid_idx,2)
        PRINT*, 'lat index range and stride', &
          src_grid_lat_beg,src_grid_lat_end,stride_lat
        PRINT*, 'lat range ',  &
          src_grid_lat(src_grid_lat_beg),src_grid_lat(src_grid_lat_end)
#endif
        lake_ct = 0; grid_ct = 0
        lake_dpth_sum = 0.0
        DO j = src_grid_lat_beg, src_grid_lat_end, stride_lat
          stride_lon = int(1.0/cos(src_grid_lat(j)*d2r)*REAL(stride_lat))
#ifdef DIAG_N_VERBOSE
        IF (j == src_grid_lat_beg) THEN
          PRINT*, 'lon index range and stride', &
            src_grid_lon_beg,src_grid_lon_end,stride_lon
          PRINT*, 'lon range ', &
            src_grid_lon(src_grid_lon_beg),src_grid_lon(src_grid_lon_end)
          IF (two_section == .true.) THEN
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
          IF (two_section == .false.) THEN
            DO i =  src_grid_lon_beg, src_grid_lon_end, stride_lon
              p(1) = src_grid_lat(j); p(2) = src_grid_lon(i) 
              p(:) = p(:)*d2r
              IF(enclosure_cnvx(v, 4, p, co_gc) == .true.) THEN
                grid_ct = grid_ct+1
                IF (lakestat((j-1)*nlon+i) /= 0) THEN
                  lake_ct = lake_ct+1
                  IF (lakedpth((j-1)*nlon+i) < 0) THEN
                    IF (lakestat((j-1)*nlon+i) == 4) THEN
                      lake_dpth_sum = lake_dpth_sum+30.0
                    ELSE
                      lake_dpth_sum = lake_dpth_sum+100.0
                    ENDIF
                  ELSE
                    lake_dpth_sum = lake_dpth_sum+REAL(lakedpth((j-1)*nlon+i))
                  ENDIF
                ENDIF
              ENDIF   
            ENDDO
          ELSE 
            DO i =  src_grid_lon_beg1, src_grid_lon_end1, stride_lon
              p(1) = src_grid_lat(j); p(2) = src_grid_lon(i) 
              p(:) = p(:)*d2r
              IF(enclosure_cnvx(v, 4, p, co_gc) == .true.) THEN
                grid_ct = grid_ct+1
                IF (lakestat((j-1)*nlon+i) /= 0) THEN
                  lake_ct = lake_ct+1
                  IF (lakedpth((j-1)*nlon+i) < 0) THEN
                    IF (lakestat((j-1)*nlon+i) == 4) THEN
                      lake_dpth_sum = lake_dpth_sum+30.0
                    ELSE
                      lake_dpth_sum = lake_dpth_sum+100.0
                    ENDIF
                  ELSE
                    lake_dpth_sum = lake_dpth_sum+REAL(lakedpth((j-1)*nlon+i))
                  ENDIF
                ENDIF
              ENDIF   
            ENDDO
            DO i =  src_grid_lon_beg2, src_grid_lon_end2, stride_lon
              p(1) = src_grid_lat(j); p(2) = src_grid_lon(i) 
              p(:) = p(:)*d2r
              IF(enclosure_cnvx(v, 4, p, co_gc) == .true.) THEN
                grid_ct = grid_ct+1
                IF (lakestat((j-1)*nlon+i) /= 0) THEN
                  lake_ct = lake_ct+1
                  IF (lakedpth((j-1)*nlon+i) < 0) THEN
                    IF (lakestat((j-1)*nlon+i) == 4) THEN
                      lake_dpth_sum = lake_dpth_sum+30.0
                    ELSE
                      lake_dpth_sum = lake_dpth_sum+100.0
                    ENDIF
                  ELSE
                    lake_dpth_sum = lake_dpth_sum+REAL(lakedpth((j-1)*nlon+i))
                  ENDIF
                ENDIF
              ENDIF   
            ENDDO
          ENDIF
        ENDDO
        cs_lakestat(cs_data_idx)=REAL(lake_ct)/REAL(grid_ct)
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
        IF (col > side_sz) THEN
          col = 2
          row = row + 2
        ENDIF
      ENDDO
      PRINT "('*'$)"  ! progress '*'
    ENDDO
    PRINT*, ''

END SUBROUTINE cal_lake_frac_depth 

 
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
    WRITE(res_ch,'(I4)') res
    if (res == 96) then
      gridfile_path = '/scratch1/NCEPDEV/nems/emc.nemspara/RT/NEMSfv3gfs/develop-20200603/INTEL/FV3_input_data/INPUT/'
    else
      gridfile_path = '/scratch1/NCEPDEV/nems/emc.nemspara/RT/NEMSfv3gfs/develop-20200603/INTEL/FV3_input_data_c'//trim(adjustl(res_ch))//'/INPUT/'
    end if 
    DO i = 1, 6
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


SUBROUTINE read_lakedata(sfcdata_path,lake_stat,lake_dpth,nlat,nlon)
    CHARACTER(len=256), INTENT(IN) :: sfcdata_path
    INTEGER*1, INTENT(OUT) :: lake_stat(:)
    INTEGER*2, INTENT(OUT) :: lake_dpth(:)
    INTEGER, INTENT(IN) :: nlat, nlon

    CHARACTER(len=256) lakefile
    INTEGER :: data_sz, i

    data_sz = nlon*nlat

   ! read in 30 arc seconds lake data 
    lakefile = trim(sfcdata_path) // "GlobalLakeStatus.dat"
    OPEN(10, file=lakefile, form='unformatted', access='direct', recl=data_sz*1)
    READ(10,rec=1) lake_stat
    CLOSE(10)

    lakefile = trim(sfcdata_path) // "GlobalLakeDepth.dat"
    OPEN(10, file=lakefile, form='unformatted', access='direct', recl=data_sz*2)
    READ(10,rec=1) lake_dpth
    CLOSE(10)

END SUBROUTINE read_lakedata
   

SUBROUTINE write_lakedata_to_orodata(cs_res, cs_lakestat, cs_lakedpth) 
    USE netcdf 
    INTEGER, INTENT(IN) :: cs_res
    REAL, INTENT(IN) :: cs_lakestat(:)
    REAL, INTENT(IN) :: cs_lakedpth(:)
   
    INTEGER :: tile_sz, tile_num
    INTEGER :: stat, ncid, x_dimid, y_dimid, varid, dimids(2)
    INTEGER :: lake_frac_id, lake_depth_id
    INTEGER :: land_frac_id, slmsk_id, geolon_id, geolat_id
    CHARACTER(len=256) :: filename,string
    CHARACTER(len=1) :: ich
    CHARACTER(len=4) res_ch
    REAL :: lake_frac(cs_res*cs_res), lake_depth(cs_res*cs_res)
    REAL :: geolon(cs_res*cs_res), geolat(cs_res*cs_res)
    REAL :: land_frac(cs_res*cs_res), slmsk(cs_res*cs_res)
    real, parameter :: epsil=1.e-6   ! numerical min for lake_frac/land_frac
    real            :: cutoff_lake=0.01, cutoff_land=1.e-6 ! if frac<cutoff, frac=0

    INTEGER :: i, j

!    include "netcdf.inc"
    tile_sz = cs_res*cs_res

    WRITE(res_ch,'(I4)') cs_res
    DO tile_num = tile_beg, tile_end
      WRITE(ich, '(I1)') tile_num
      filename = "C" // trim(adjustl(res_ch)) // "_oro_data.tile" // ich // ".nc" 
      filename = "oro_data.tile" // ich // ".nc" 
!      filename = "oro.C" // trim(adjustl(res_ch)) // ".tile" // ich // ".nc" 
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
      stat = nf90_put_att(ncid, lake_frac_id,'long_name','lake fraction')
      CALL nc_opchk(stat, "nf90_put_att: lake_frac:long_name") 
      stat = nf90_put_att(ncid, lake_frac_id,'unit','fraction')
      CALL nc_opchk(stat, "nf90_put_att: lake_frac:unit") 
      write(string,'(a,f5.2)') 'based on GLDBv2 (Choulga et al. 2014); missing Caspian/Aral Sea added based on land_frac in this dataset. lake_frac cutoff is',cutoff_lake
      stat = nf90_put_att(ncid, lake_frac_id,'description',trim(string))
      CALL nc_opchk(stat, "nf90_put_att: lake_frac:description") 
#endif
      stat = nf90_def_var(ncid,"lake_depth",NF90_FLOAT,dimids,lake_depth_id)
      CALL nc_opchk(stat, "nf90_def_var: lake_depth")
#ifdef ADD_ATT_FOR_NEW_VAR
      stat = nf90_put_att(ncid, lake_depth_id,'long_name','lake depth')
      CALL nc_opchk(stat, "nf90_put_att: lake_depth:long_name") 
      stat = nf90_put_att(ncid, lake_depth_id,'unit','meter')
      CALL nc_opchk(stat, "nf90_put_att: lake_depth:long_name") 
      stat = nf90_put_att(ncid, lake_depth_id,'description', &
        'based on GLDBv2 (Choulga et al. 2014); missing depth set to 211m &
         for Caspian Sea and 10m for Aral Sea; neg. depths are replaced by &
         default values of 10m for lake and 3m for river, spurious large pos. &
         depths are left unchanged.')
      CALL nc_opchk(stat, "nf90_put_att: lake_depth:description") 
#endif

      write(string,'(a,es8.1)') 'land_frac with lake; land_frac cutoff is',cutoff_land
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


! add Caspian Sea and Aral Sea to lake_frac and lake_depth
      lake_frac(:) = cs_lakestat((tile_num-1)*tile_sz+1:tile_num*tile_sz)
      lake_depth(:) = cs_lakedepth((tile_num-1)*tile_sz+1:tile_num*tile_sz)
      IF (tile_num == 2 .or. tile_num == 3) THEN
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

      if (min(cutoff_lake,cutoff_land) < epsil) then
        print *,'cutoff_lake/cutoff_land cannot be smaller than epsil, reset...'
        cutoff_lake=max(epsil,cutoff_lake)
        cutoff_land=max(epsil,cutoff_land)
      end if
      DO i = 1, tile_sz
!      IF (tile_num == 5 .and. abs(land_frac(i)-0.318) < 0.01 ) print *,'qq1 i=',i,' land/lake=',land_frac(i),lake_frac(i)
! epsil is "numerical" nonzero min for lake_frac/land_frac
        if (lake_frac(i)>0. .and. lake_frac(i)<   epsil) lake_frac(i)=0.
        if (lake_frac(i)<1. .and. lake_frac(i)>1.-epsil) lake_frac(i)=1.
        if (land_frac(i)>0. .and. land_frac(i)<   epsil) land_frac(i)=0.
        if (land_frac(i)<1. .and. land_frac(i)>1.-epsil) land_frac(i)=1.
! cutoff_lake/cutoff_land is practical min for lake_frac/land_frac
        if (land_frac(i)>0. .and. land_frac(i)<cutoff_land) then
          land_frac(i)=0.
        end if
! recalculate lake_frac, land_frac/slmsk, assuming lake/ocean don't coexist, so that max(lake,ocn)+land=1
!        if (land_frac(i) > 0.) then
!          if (lake_frac(i) >= cutoff_lake) then ! lake_frac dominates over land_frac
!            land_frac(i) = min(1., land_frac(i)+lake_frac(i)) - lake_frac(i)
!            land_frac(i) = max(1., land_frac(i)+lake_frac(i)) - lake_frac(i)
!          else ! land_frac dominiates over lake_frac
!            lake_frac(i) = min(1., land_frac(i)+lake_frac(i)) - land_frac(i)
!            lake_frac(i) = max(1., land_frac(i)+lake_frac(i)) - land_frac(i)
!          end if
!        end if

        if (lake_frac(i)>0. .and. lake_frac(i)<cutoff_lake) then
          lake_frac(i)=0.
          lake_depth(i)=0.
          land_frac(i)=1.
        elseif (lake_frac(i) >= cutoff_lake .and. lake_depth(i)==0.) then
!          print *,'lake without depth at i=',i
          lake_depth(i)=10.
        end if

        slmsk(i) = nint(land_frac(i))
!      IF (tile_num == 5 .and. (i==113375.or.i==98441.or.i==147232.or.i==108080)) print *,'qq2 i=',i,' land/lake=',land_frac(i),lake_frac(i)
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

SUBROUTINE write_lakefrac_netcdf(cs_res, cs_lakestat, cs_lakedpth) 
    USE netcdf 
    INTEGER, INTENT(IN) :: cs_res
    REAL, INTENT(IN) :: cs_lakestat(:)
    REAL, INTENT(IN) :: cs_lakedpth(:)
   
    INTEGER :: tile_sz, tile_num
    INTEGER :: stat, ncid, x_dimid, y_dimid, varid, dimids(2)
    INTEGER :: lake_frac_id, lake_depth_id
    CHARACTER(len=256) :: filename
    CHARACTER(len=1) :: ich
    CHARACTER(len=4) res_ch
    REAL :: data_out(cs_res*cs_res)

!    include "netcdf.inc"
    tile_sz = cs_res*cs_res

    WRITE(res_ch,'(I4)') cs_res
    DO tile_num = 1, 6
      WRITE(ich, '(I1)') tile_num
      filename = "C" // trim(adjustl(res_ch)) // "lake_frac.tile" // ich // ".nc" 
      stat = nf90_create(filename, NF90_CLOBBER, ncid)
      CALL nc_opchk(stat, "nf90_create lake_frac.nc")
      stat = nf90_def_dim(ncid, "x", cs_res, x_dimid)
      CALL nc_opchk(stat, "nf90_def_dim: x")
      stat = nf90_def_dim(ncid, "y", cs_res, y_dimid)
      CALL nc_opchk(stat, "nf90_def_dim: y")
      dimids = (/ y_dimid, x_dimid /)
     
      stat = nf90_def_var(ncid,"lake_frac",NF90_FLOAT,dimids,lake_frac_id)
      CALL nc_opchk(stat, "nf90_def_var: lake_frac")
      stat = nf90_def_var(ncid,"lake_depth",NF90_FLOAT,dimids,lake_depth_id)
      CALL nc_opchk(stat, "nf90_def_var: lake_depth")
      stat = nf90_enddef(ncid) 
      CALL nc_opchk(stat, "nf90_enddef")

      data_out(:) = cs_lakestat((tile_num-1)*tile_sz+1:tile_num*tile_sz)
      stat = nf90_put_var(ncid, lake_frac_id, data_out, &
             start = (/ 1, 1 /), count = (/ cs_res, cs_res /) )
      CALL nc_opchk(stat, "nf90_put_var: lake_frac") 
      data_out(:) = cs_lakedepth((tile_num-1)*tile_sz+1:tile_num*tile_sz)
      stat = nf90_put_var(ncid, lake_depth_id, data_out, &
             start = (/ 1, 1 /), count = (/ cs_res, cs_res /) )
      CALL nc_opchk(stat, "nf90_put_var: lake_depth") 
    
      stat = nf90_close(ncid)
      CALL nc_opchk(stat, "nf90_close") 
    ENDDO
  
END SUBROUTINE write_lakefrac_netcdf

SUBROUTINE nc_opchk(stat,opname)
   USE netcdf
   IMPLICIT NONE
   INTEGER stat
   CHARACTER(len=*) opname
   CHARACTER(64) msg

   IF (stat .NE.0)  THEN
     msg = trim(opname) // ' Error, status code and message:'
     PRINT*,trim(msg), stat, nf90_strerror(stat)
     STOP
   END IF

END SUBROUTINE nc_opchk

END PROGRAM lake_frac
