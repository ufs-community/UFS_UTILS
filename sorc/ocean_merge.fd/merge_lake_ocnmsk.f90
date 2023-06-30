!!! input are from pth1 and pth2, output will be in the local dir
!!! To compile standalone do
!!! ifort merge_lake_ocnmsk.F90 -I$NETCDF/include -L$NETCDF/lib -lnetcdff -lnetcdf
!!! created by Shan.Sun  modified by Rahul Mahajan and Sanath Kumar
!!!
!!!

!!!
!
!!! "Determine the water mask by merging the lake mask with the mapped ocean
!!! mask from MOM6, both are on the FV3 grid. During merge, the ocean mask
!!! dominates the lake mask if there is a conflict. After the merge, the remaining
!!! non-water fraction is the land fraction.m"
!!!
!!!




program merge_lake_ocnmsk
  use netcdf

  implicit none

  character(len=120) :: pth1
  character(len=120) :: pth2,pth3
  character(len=10)  :: atmres,ocnres
  real, parameter    :: min_land=1.e-4, def_lakedp=10.
  logical, parameter :: int_lake=.true.

  character(len=120) :: flnm
  integer :: ncid,ndims,nvars,natts,lat,lon,v1id,v2id,v3id,v4id,start(2),count(2),i,j,latid,lonid,ncid4, dims(2),tile,nodp_pt,lake_pt,vlat
  real, allocatable :: lake_frac(:,:),lake_depth(:,:),land_frac(:,:),ocn_frac(:,:),slmsk(:,:),lat2d(:,:)

  call read_nml(pth1, pth2, atmres, ocnres, pth3)

  nodp_pt=0
  lake_pt=0
  do tile=1,6
    write(flnm,'(5a,i1,a)') trim(pth1),trim(atmres),'.',trim(ocnres),'.tile',tile,'.nc'
    call handle_err (nf90_open (flnm, NF90_NOWRITE, ncid))
    call handle_err (nf90_inquire (ncid, ndimensions=ndims, nvariables=nvars, nattributes=natts))
    write(6,*) 'flnm_ocn=',flnm,' ncid=',ncid, ' ndims=',ndims, 'nvars=',nvars,' natts=',natts
    call handle_err (nf90_inq_dimid (ncid, 'grid_xt', latid))  ! RM: lon is no longer in this file; try grid_xt
    call handle_err (nf90_inq_dimid (ncid, 'grid_yt', lonid))  ! RM: lat is no longer in this file; try grid_tyt
    call handle_err (nf90_inquire_dimension (ncid, latid, len=lat))
    call handle_err (nf90_inquire_dimension (ncid, lonid, len=lon))
    if (tile==1) then
      write(6,*) 'lat=',lat,'lon=',lon
      allocate (lake_frac(lon,lat),lake_depth(lon,lat),land_frac(lon,lat),ocn_frac(lon,lat),slmsk(lon,lat),lat2d(lon,lat))
      start(1:2) = (/1,1/)
      count(1:2) = (/lon,lat/)
    end if
    call handle_err (nf90_inq_varid(ncid, 'land_frac', v1id))
    call handle_err (nf90_get_var (ncid, v1id, ocn_frac, start=start, count=count))

    write(flnm,'(3a,i1,a)') trim(pth2),trim(atmres),'_oro_data.tile',tile,'.nc'
    print *,' flnm2=',trim(flnm)
    call handle_err (nf90_open (flnm, NF90_NOWRITE, ncid))
    call handle_err (nf90_inquire (ncid, ndimensions=ndims, nvariables=nvars, nattributes=natts))
    write(6,*) 'flnm_lake=',flnm,' ncid=',ncid, ' ndims=',ndims, 'nvars=',nvars,' natts=',natts
   !if (tile==1) then
   !  call handle_err (nf90_inq_dimid (ncid, 'lat', latid))
   !  call handle_err (nf90_inq_dimid (ncid, 'lon', lonid))
   !  call handle_err (nf90_inquire_dimension (ncid, latid, len=lat))
   !  call handle_err (nf90_inquire_dimension (ncid, lonid, len=lon))
   !  write(6,*) 'lat=',lat,'lon=',lon
   !  start(1:2) = (/1,1/)
   !  count(1:2) = (/lon,lat/)
   !end if
    call handle_err (nf90_inq_varid(ncid, 'lake_frac', v2id))
    call handle_err (nf90_inq_varid(ncid, 'lake_depth',v3id))
    call handle_err (nf90_inq_varid(ncid, 'geolat'    ,vlat))
    call handle_err (nf90_get_var (ncid, v2id, lake_frac, start=start, count=count))
    call handle_err (nf90_get_var (ncid, v3id, lake_depth,start=start, count=count))
    call handle_err (nf90_get_var (ncid, vlat, lat2d,     start=start, count=count))

    do i=1,lon
    do j=1,lat
      if (int_lake) lake_frac(i,j)=nint(lake_frac(i,j))     ! using integer lake_frac
      if (lat2d(i,j).le.-60.) lake_frac(i,j)=0.             ! ignore lakes on Antarctica
      land_frac(i,j)=1.-ocn_frac(i,j)
      if (land_frac(i,j) <    min_land) land_frac(i,j)=0.   ! ignore land  < min_land
      if (land_frac(i,j) > 1.-min_land) land_frac(i,j)=1.   ! ignore water < min_land
      if (1.-land_frac(i,j) > 0.) lake_frac(i,j)=0.         ! ocn dominates

      if (lake_frac(i,j) > 0.) then
        lake_pt=lake_pt+1            ! calculating total lake points
        if (int_lake) then
          land_frac(i,j)=0.
        else
          land_frac(i,j)=1.-lake_frac(i,j)
        end if
        if (lake_depth(i,j) <= 0.) then
          lake_depth(i,j)=def_lakedp ! set missing lake depth to default value
          nodp_pt=nodp_pt+1          ! calculating total lake points without depth
        end if
      else
        lake_depth(i,j)=0.
      end if
      slmsk(i,j) = ceiling(land_frac(i,j)) ! ceiling is used for orog smoothing
    end do
    end do

    write(flnm,'(4a,i1,a)') trim(atmres),'.',trim(ocnres),'.tile',tile,'.nc'
    print *,'output=',trim(flnm)
    call handle_err (nf90_create (path=trim(pth3)//trim(flnm), &
    cmode=or(NF90_CLOBBER, NF90_64BIT_OFFSET), ncid=ncid4))   ! netcdf3

    call handle_err (nf90_def_dim (ncid4,'lon', lon, dims(1)))
    call handle_err (nf90_def_dim (ncid4,'lat', lat, dims(2)))
    call handle_err (nf90_def_var (ncid4,'land_frac', nf90_float, dims(1:2), v1id))
    call handle_err (nf90_def_var (ncid4,'lake_frac', nf90_float, dims(1:2), v2id))
    call handle_err (nf90_def_var (ncid4,'lake_depth',nf90_float, dims(1:2), v3id))
    call handle_err (nf90_def_var (ncid4,'slmsk',     nf90_float, dims(1:2), v4id))

    call handle_err (nf90_enddef  (ncid4))

    call handle_err (nf90_put_var (ncid4, v1id,land_frac))
    call handle_err (nf90_put_var (ncid4, v2id,lake_frac))
    call handle_err (nf90_put_var (ncid4, v3id,lake_depth))
    call handle_err (nf90_put_var (ncid4, v4id,slmsk))
    call handle_err (nf90_close(ncid4))

!   do i=1,48
!   do j=1,48
!     write(*,'(2i4,4f6.1)') i,j,land_frac(i,j),lake_frac(i,j),lake_depth(i,j),slmsk(i,j)
!   end do
!   end do

  end do ! tile
  write(*,'(a,i8,a,i8,a)') 'total lake point ',lake_pt,' where ',nodp_pt,' has no depth'

end program merge_lake_ocnmsk

subroutine handle_err (ret)
  use netcdf
  integer, intent(in) :: ret

  if (ret /= NF90_NOERR) then
    write(6,*) nf90_strerror (ret)
    stop 999
  end if
end subroutine handle_err

subroutine read_nml(ocean_mask_dir, lake_mask_dir, atmres, ocnres,out_dir)

  integer :: unit=7, io_status

  character(len=120), intent(out) :: ocean_mask_dir
  character(len=120), intent(out) :: lake_mask_dir
  character(len=120), intent(out) :: out_dir
  character(len=10),  intent(out) :: atmres,ocnres

  !ocean_mask_dir='/scratch1/NCEPDEV/stmp4/Sanath.Kumar/CPLD_GRIDGEN/100/' ! MOM6 mask dir
  !lake_mask_dir='/scratch1/NCEPDEV/stmp4/Sanath.Kumar/CPLD_GRIDGEN/100/' ! this is for P8
  !atmres='C96'
  !ocnres='mx100'
  !out_dir='/scratch1/BMC/gsd-fv3-dev/fv3data/edit_landfrac_mom6_nov2021/' 

  namelist/mask_nml/ocean_mask_dir, lake_mask_dir, atmres, ocnres, out_dir
  open(unit=unit, file='input.nml', iostat=io_status )
  read(unit,mask_nml, iostat=io_status )
  close(unit)
  if (io_status > 0) call handle_err(-1, 'Error reading input.nml')

end subroutine read_nml
