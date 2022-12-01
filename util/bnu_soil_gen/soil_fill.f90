use netcdf
implicit none

character*3, parameter  :: level    = "top"  ! top or bot
integer, parameter      :: numlat   = 21600, numlon   = 43200
integer, parameter      :: numlatt  = 3600 , numlont  = 3600

integer(kind=1), dimension(:,:), allocatable :: soil_texture, veg_class
integer,         dimension(numlont,numlatt)  :: insoil

character*3 :: latnames(0:6)  = (/"-90","-60","-30","+00","+30","+60","+90"/)
character*4 :: lonnames(0:12) = (/"-180","-150","-120","-090","-060","-030", &
                           "+000","+030","+060","+090","+120","+150","+180"/)

integer :: countmax,countwater,countsnow,countcur,isoil,maxtype,num_missing
integer :: iret, ncid, varid(2), dimid(3)  ! Variables for NetCDF access

integer(kind=1), parameter :: missing = -9
integer :: itile,jtile,i,j,imax,imin,jmax,jmin
character (len=128)  :: filename

allocate(soil_texture(numlon,numlat))
allocate(veg_class   (numlon,numlat))

do jtile = 0,5     ! 30deg x 30deg tiles
do itile = 0,11

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! create name for soil files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

filename = "lat_"//latnames(jtile+1)//"_"//latnames(jtile)//"_lon_"//lonnames(itile)//"_"//lonnames(itile+1)//".nc"
print*, "Reading: ",trim(filename)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read in netcdf soil data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(jtile > 0) then

  iret = nf90_open("tiles/"//trim(filename),NF90_NOWRITE,ncid)
    if(iret /= nf90_noerr) stop "error opening file"
  if(level == "top") iret = nf90_inq_varid(ncid,"soil_texture_top",varid(1)) 
  if(level == "bot") iret = nf90_inq_varid(ncid,"soil_texture_bot",varid(1)) 
    if(iret /= nf90_noerr) stop "error accessing soil variable"

  iret = nf90_get_var(ncid,varid(1),insoil) 
    if(iret /= nf90_noerr) stop "error reading soil"

  iret = nf90_close(ncid)

else

  insoil = -9   ! set texture to missing for 60S-90S

end if

where(insoil < 0) insoil = -9

soil_texture(itile*3600+1:(itile+1)*3600,jtile*3600+1:(jtile+1)*3600) = insoil

end do  ! itile loop
end do  ! jtile loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read in netcdf vegclass data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

filename = "/scratch2/NCEPDEV/land/data/input_data/vegetation_type/vegetation_type.viirs.igbp.30s.nc"
print*, "Reading: ",filename

  iret = nf90_open(filename,NF90_NOWRITE,ncid)
    if(iret /= nf90_noerr) stop "error opening file"
  iret = nf90_inq_varid(ncid,"vegetation_type",varid(1)) 
    if(iret /= nf90_noerr) stop "error accessing veg variable"

  iret = nf90_get_var(ncid,varid(1),veg_class,start=(/1,1,1/),count=(/43200,21600,1/)) 
    if(iret /= nf90_noerr) stop "error reading veg"

  iret = nf90_close(ncid)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! do some checking
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  where(soil_texture < 0 .and. veg_class == 15) soil_texture = 16  ! set ice to 16
  where(soil_texture < 0 .and. veg_class == 17) soil_texture = 14  ! set water to 14
  where(soil_texture < 0 .and. veg_class == -9) soil_texture = 14  ! set water to 14

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! first check 10x10
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  num_missing = count(soil_texture < 0)
  print *, "Num missing before fill: ", num_missing
  if(num_missing > 0) then
    do i = 1,43200
    do j = 1,21600
      if(soil_texture(i,j) < 0) then
       imax = min(i+5,43200)
       imin = max(i-5,1)
       jmax = min(j+5,21600)
       jmin = max(j-5,1)
       countmax = 0
       maxtype = -1
       do isoil = 1,12
         countcur = count(soil_texture(imin:imax,jmin:jmax) == isoil)
	 if(countcur > countmax) maxtype = isoil
	 if(countcur > countmax) countmax = countcur
       end do
       if(countmax > 0) soil_texture(i,j) = maxtype
      end if
    end do
    end do
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! next check 50x50
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  num_missing = count(soil_texture < 0)
  print *, "Num missing after fill #1: ", num_missing
  if(num_missing > 0) then
    do i = 1,43200
    do j = 1,21600
      if(soil_texture(i,j) < 0) then
       imax = min(i+50,43200)
       imin = max(i-50,1)
       jmax = min(j+50,21600)
       jmin = max(j-50,1)
       countmax = 0
       maxtype = -1
       do isoil = 1,12
         countcur = count(soil_texture(imin:imax,jmin:jmax) == isoil)
	 if(countcur > countmax) maxtype = isoil
	 if(countcur > countmax) countmax = countcur
       end do
       if(countmax > 0) soil_texture(i,j) = maxtype
      end if
    end do
    end do
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! next check 200x200
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  num_missing = count(soil_texture < 0)
  print *, "Num missing after fill #2: ", num_missing
  if(num_missing > 0) then
    do i = 1,43200
    do j = 1,21600
      if(soil_texture(i,j) < 0) then
       imax = min(i+200,43200)
       imin = max(i-200,1)
       jmax = min(j+200,21600)
       jmin = max(j-200,1)
       countmax = 0
       maxtype = -1
       do isoil = 1,12
         countcur = count(soil_texture(imin:imax,jmin:jmax) == isoil)
	 if(countcur > countmax) maxtype = isoil
	 if(countcur > countmax) countmax = countcur
       end do
       if(countmax > 0) soil_texture(i,j) = maxtype
      end if
    end do
    end do
    print *, "Num missing after fill #3: ",count(soil_texture < 0)
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! next just fill with loam type
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  where(soil_texture < 0) soil_texture = 6

  print *, "Num missing after filling: ",count(soil_texture < 0)
  if(count(soil_texture < 0) > 0) stop "too many missing"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! for ufs_utils use, set water to missing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  where(soil_texture == 14) soil_texture = missing

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! write out in netcdf format
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  iret = nf90_create("bnu_soil.nc", NF90_CLOBBER, ncid)
    if(iret /= nf90_noerr) stop "error creating file"

  iret = nf90_def_dim(ncid, "idim", 43200, dimid(1))
  iret = nf90_def_dim(ncid, "jdim", 21600, dimid(2))
  iret = nf90_def_dim(ncid, "time", 1    , dimid(3))

!  iret = nf90_def_var(ncid,  "veg"    , NF90_BYTE, dimid, varid(1))
  iret = nf90_def_var(ncid,  "soil_type"   , NF90_BYTE, dimid, varid(2))
      iret = nf90_put_att(ncid, varid(2), "missing_value", missing)
      iret = nf90_put_att(ncid, varid(2), "landice_category", 16)

  iret = nf90_enddef(ncid)

!  iret = nf90_put_var(ncid, varid(1), veg_class)
  iret = nf90_put_var(ncid, varid(2), soil_texture, start = (/1,1,1/), count = (/43200,21600,1/))
  
  iret = nf90_close(ncid)

end
