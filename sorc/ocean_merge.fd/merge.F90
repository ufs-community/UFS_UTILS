!> Determine the water mask by merging the lake mask with the mapped ocean
!! mask from MOM6. Both are on the FV3 grid. During merge, the ocean mask
!! dominates the lake mask if there is a conflict. After the merge, the remaining
!! non-water fraction is the land fraction.
!! 
!! @param[in] lon The "east/west" dimension of the model grid.
!! @param[in] lat The "north/south" dimension of the model grid.
!! @param[in] binary_lake When '1', lake fraction is either 0 or 1. Otherwise, it is a fraction.
!! @param[in] lat2d Latitude of the model grid points.
!! @param[in] ocn_frac Fraction of the grid point that is ocean.
!! @param[inout] lake_frac Fraction of the grid point that is lake.
!! @param[inout] lake_depth Lake depth in meters.
!! @param[out] land_frac Fraction of the grid point that is land.
!! @param[out] slmsk Land/sea mask. '0' if less than 50% land. Otherwise, '1'.
!!
!! @author Shan Sun
!! @author Rahul Mahajan
!! @author Sanath Kumar
 subroutine merge(lon, lat, binary_lake, lat2d, ocn_frac, &
                  lake_frac, lake_depth, land_frac, slmsk)

 implicit none

 integer, intent(in)        :: lon, lat, binary_lake

 real, intent(in)           :: lat2d(lon,lat)
 real, intent(in)           :: ocn_frac(lon,lat)
 real, intent(inout)        :: lake_frac(lon,lat)
 real, intent(inout)        :: lake_depth(lon,lat)
 real, intent(out)          :: land_frac(lon,lat)
 real, intent(out)          :: slmsk(lon,lat)

 real, parameter            :: min_land=1.e-4, def_lakedp=10.

 integer                    :: i, j, nodp_pt, lake_pt

 nodp_pt=0
 lake_pt=0

 do i=1,lon
 do j=1,lat
   if (binary_lake.eq.1) lake_frac(i,j)=nint(lake_frac(i,j))     ! using integer lake_frac
   if (lat2d(i,j).le.-60.) lake_frac(i,j)=0.             ! ignore lakes on Antarctica
   land_frac(i,j)=1.-ocn_frac(i,j)
   if (land_frac(i,j) <    min_land) land_frac(i,j)=0.   ! ignore land  < min_land
   if (land_frac(i,j) > 1.-min_land) land_frac(i,j)=1.   ! ignore water < min_land
   if (1.-land_frac(i,j) > 0.) lake_frac(i,j)=0.         ! ocn dominates

   if (lake_frac(i,j) > 0.) then
     lake_pt=lake_pt+1            ! calculating total lake points
     if (binary_lake.eq.1) then
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
   slmsk(i,j) = nint(land_frac(i,j)) ! nint got the land pts correct
 end do
 end do

 write(*,'(a,i8,a,i8,a)') 'Total lake point ',lake_pt,' where ',nodp_pt,' has no depth'

 end subroutine merge
