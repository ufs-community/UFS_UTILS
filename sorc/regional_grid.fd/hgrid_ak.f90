!=============================================================================
subroutine hgrid_ak(lx,ly,nx,ny,a,k,plat,plon,pazi, &
                    re,delx,dely,  glat,glon,garea, ff)
!=============================================================================
! Use a and k as the parameters of a generalized Schmidt-transformed
! gnomonic mapping centered at (plat,plon) and twisted about this center
! by an azimuth angle of pazi counterclockwise (these angles in radians).
!
! Assuming the radius of the earth is re, and using the central mapping
! point as the coordinate origin, set up the grid with central x-spacing delx
! and y-spacing dely in physical units, and with the location of the left-lower
! corner of the grid at center-relative grid index pair, (lx,ly) and with
! the number of the grid spaces in x and y directions given by nx and ny.
! (Note that, for a centered rectangular grid lx and ly are negative and, in
! magnitude, half the values of nx and ny respectively.)
! Return the latitude and longitude, in radians again, of the grid points thus
! defined in the arrays, glat and glon, and return a rectangular array, garea,
! of dimensions nx-1 by ny-1, that contains the areas of each of the grid cells
! in the SQUARE of the same physical length unit that was employed to define
! the radius of the earth, re (and the central grid steps, delx and dely).
!
! if all goes well, return a .FALSE. failure flag, ff. If, for some
! reason, it is not possible to complete this task, return the failure flag
! as .TRUE.
!=============================================================================
use pkind, only: dp
use pietc, only: u0,u1,dtor
use pmat4, only: sarea
use pmat5, only: ctog
implicit none
integer,                                  intent(in ):: lx,ly,nx,ny
real(dp),                                 intent(in ):: a,k,plat,plon,pazi, &
                                                        re,delx,dely
real(dp),dimension(lx:lx+nx  ,ly:ly+ny  ),intent(out):: glat,glon
real(dp),dimension(lx:lx+nx-1,ly:ly+ny-1),intent(out):: garea
logical,                                  intent(out):: ff
!-----------------------------------------------------------------------------
real(dp),dimension(3,3):: prot,azirot
real(dp),dimension(3,2):: xcd
real(dp),dimension(3)  :: xc
real(dp),dimension(2)  :: xm
real(dp)               :: clat,slat,clon,slon,cazi,sazi,&
                          rlat,drlata,drlatb,drlatc,    &
                          rlon,drlona,drlonb,drlonc,   rre
integer                :: ix,iy,mx,my
!=============================================================================
clat=cos(plat); slat=sin(plat)
clon=cos(plon); slon=sin(plon)
cazi=cos(pazi); sazi=sin(pazi)

azirot(:,1)=(/ cazi, sazi, u0/)
azirot(:,2)=(/-sazi, cazi, u0/)
azirot(:,3)=(/   u0,   u0, u1/)

prot(:,1)=(/     -slon,       clon,    u0/)
prot(:,2)=(/-slat*clon, -slat*slon,  clat/)
prot(:,3)=(/ clat*clon,  clat*slon,  slat/)
prot=matmul(prot,azirot)
mx=lx+nx ! Index of the 'right' edge of the rectangular grid
my=ly+ny ! Index of the 'top' edge of the rectangular grid
!This code assumes symmetry about the grid center
do iy=ly,my
   !xm(2)=iy*dely/re
   xm(2)=-iy*dely/re
   do ix=lx,mx
      !xm(1)=ix*delx/re
      xm(1)=-ix*delx/re
      call xmtoxc_ak(a,k,xm,xc,xcd,ff)
      if(ff)return
      xcd=matmul(prot,xcd)
      xc =matmul(prot,xc )
      call ctog(xc,glat(ix,iy),glon(ix,iy))
   enddo
enddo

! Convert degrees to radians in the glat and glon arrays:
glat=glat*dtor
glon=glon*dtor

! Compute the areas of the quadrilateral grid cells:
do iy=ly,my-1
   do ix=lx,mx-1
      rlat  =glat(ix  ,iy  )
      drlata=glat(ix+1,iy  )-rlat
      drlatb=glat(ix+1,iy+1)-rlat
      drlatc=glat(ix  ,iy+1)-rlat
      rlon  =glon(ix  ,iy  )
      drlona=glon(ix+1,iy  )-rlon
      drlonb=glon(ix+1,iy+1)-rlon
      drlonc=glon(ix  ,iy+1)-rlon
! If 'I' is the grid point (ix,iy), 'A' is (ix+1,iy); 'B' is (ix+1,iy+1)
! and 'C' is (ix,iy+1), then the area of the grid cell IABC is the sum of
! the areas of the traingles, IAB and IBC (the latter being the negative
! of the signed area of triangle, ICB):
      garea(ix,iy)=sarea(rlat, drlata,drlona, drlatb,drlonb) &
                  -sarea(rlat, drlatc,drlonc, drlatb,drlonb)
   enddo
enddo
! Convert the areas to area units consistent with the length units used for
! the radius, re, of the sphere:
rre=re*re
garea=garea*rre

end subroutine hgrid_ak
