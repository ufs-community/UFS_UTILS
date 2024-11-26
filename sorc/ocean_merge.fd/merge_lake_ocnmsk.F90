!> @file
!! @brief Determines the water mask by merging the lake mask with
!! the mapped ocean mask from MOM6.
!! @author Shan Sun
!! @author Rahul Mahajan
!! @author Sanath Kumar

!> Determine the water mask by merging the lake mask with the mapped ocean
!! mask from MOM6, both are on the FV3 grid. During merge, the ocean mask
!! dominates the lake mask if there is a conflict. After the merge, the remaining
!! non-water fraction is the land fraction.
!! 
!! @return 0 for success, error code otherwise.
!! @author Shan Sun
!! @author Rahul Mahajan
!! @author Sanath Kumar
program merge_lake_ocnmsk

  implicit none

  character(len=120) :: pth1
  character(len=120) :: pth2,pth3
  character(len=10)  :: atmres,ocnres

  integer :: binary_lake
  integer :: lat,lon,tile

  real, allocatable :: lake_frac(:,:),lake_depth(:,:),land_frac(:,:),ocn_frac(:,:),slmsk(:,:),lat2d(:,:)

  print*,"- BEGIN OCEAN MERGE PROGRAM."

  call read_nml(pth1, pth2, atmres, ocnres, pth3,binary_lake)

  do tile=1,6

    call read_grid_dims(pth1, atmres, ocnres, tile, lon, lat)

    if (tile==1) then
      write(6,*) 'lat=',lat,'lon=',lon
      allocate (lake_frac(lon,lat),lake_depth(lon,lat),land_frac(lon,lat), &
                ocn_frac(lon,lat),slmsk(lon,lat),lat2d(lon,lat))
    endif

    call read_ocean_frac(pth1,atmres,ocnres,tile,lon,lat,ocn_frac)

    call read_lake_mask(pth2,atmres,tile,lon,lat,lake_frac, &
                        lake_depth,lat2d)

    call merge(lon, lat, binary_lake, lat2d, ocn_frac, lake_frac, lake_depth, land_frac, slmsk)

    call write_data(atmres,ocnres,pth3,tile,lon,lat,land_frac, &
                    lake_frac,lake_depth,slmsk)

  end do ! tile

  deallocate (lake_frac,lake_depth,land_frac,ocn_frac,slmsk,lat2d)

  print*,"- NORMAL TERMINATION."

end program merge_lake_ocnmsk
