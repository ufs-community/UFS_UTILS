!=============================================================================
subroutine hgrid_ak(lx,ly,nx,ny,a,k,plat,plon,pazi, &
                    re,delx,dely,  glat,glon,garea,dx,dy,angle_dx,angle_dy, ff)
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
! the radius of the earth, re, and the central grid steps, delx and dely.
! In this version, these grid cell areas are computed by 2D integrating the
! scalar jacobian of the transformation, using a 4th-order centered scheme.
! The estimated grid steps, dx snd dy, are returned at the grid cell edges,
! using the same 4th-order scheme to integrate the 1D projected jacobian.
! The angles, relative to local east and north, are returned respectively
! as angle_dx and angle_dy at grid cell corners, in degrees counterclockwise.
!
! if all goes well, return a .FALSE. failure flag, ff. If, for some
! reason, it is not possible to complete this task, return the failure flag
! as .TRUE.
!=============================================================================
use pkind, only: dp
use pietc, only: u0,u1,dtor,rtod
use pmat4, only: cross_product,triple_product
use pmat5, only: ctog
implicit none
integer,                                  intent(in ):: lx,ly,nx,ny
real(dp),                                 intent(in ):: a,k,plat,plon,pazi, &
                                                        re,delx,dely
real(dp),dimension(lx:lx+nx  ,ly:ly+ny  ),intent(out):: glat,glon
real(dp),dimension(lx:lx+nx-1,ly:ly+ny-1),intent(out):: garea
real(dp),dimension(lx:lx+nx-1,ly:ly+ny  ),intent(out):: dx
real(dp),dimension(lx:lx+nx  ,ly:ly+ny-1),intent(out):: dy
real(dp),dimension(lx:lx+nx  ,ly:ly+ny  ),intent(out):: angle_dx,angle_dy
logical,                                  intent(out):: ff
!-----------------------------------------------------------------------------
real(dp),dimension(lx-1:lx+nx+1,ly-1:ly+ny+1):: gat ! Temporary area array
real(dp),dimension(lx-1:lx+nx+1,ly  :ly+ny  ):: dxt ! Temporary dx array
real(dp),dimension(lx  :lx+nx  ,ly-1:ly+ny+1):: dyt ! Temporary dy array
real(dp),dimension(3,3):: prot,azirot
real(dp),dimension(3,2):: xcd,eano
real(dp),dimension(2,2):: xcd2
real(dp),dimension(3)  :: xc,east,north
real(dp),dimension(2)  :: xm
real(dp)               :: clat,slat,clon,slon,cazi,sazi,&
                          rlat,drlata,drlatb,drlatc,    &
                          rlon,drlona,drlonb,drlonc,   delxy,delxore,delyore
integer                :: ix,iy,mx,my,lxm,lym,mxp,myp
!=============================================================================
delxore=delx/re
delyore=dely/re
delxy=delx*dely
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
lxm=lx-1; mxp=mx+1 ! Indices of extra left and right edges
lym=ly-1; myp=my+1 ! Indices of extra bottom and top edges

!-- main body of horizontal grid:
do iy=ly,my
   xm(2)=iy*delyore
   do ix=lx,mx
      xm(1)=ix*delxore
      call xmtoxc_ak(a,k,xm,xc,xcd,ff); if(ff)return
      xcd=matmul(prot,xcd)
      xc =matmul(prot,xc )
      call ctog(xc,glat(ix,iy),glon(ix,iy))
      east=(/-xc(2),xc(1),u0/); east=east/sqrt(dot_product(east,east))
      north=cross_product(xc,east)
      eano(:,1)=east; eano(:,2)=north
      xcd2=matmul(transpose(eano),xcd)
      angle_dx(ix,iy)=atan2( xcd2(2,1),xcd2(1,1))*rtod
      angle_dy(ix,iy)=atan2(-xcd2(1,2),xcd2(2,2))*rtod
      dxt(ix,iy)=sqrt(dot_product(xcd2(:,1),xcd2(:,1)))*delx
      dyt(ix,iy)=sqrt(dot_product(xcd2(:,2),xcd2(:,2)))*dely
      gat(ix,iy)=triple_product(xc,xcd(:,1),xcd(:,2))*delxy
   enddo
enddo

!-- extra left edge, gat, dxt only:
xm(1)=lxm*delxore
do iy=ly,my
   xm(2)=iy*delyore
   call xmtoxc_ak(a,k,xm,xc,xcd,ff); if(ff)return
   xcd=matmul(prot,xcd)
   xc =matmul(prot,xc )
   east=(/-xc(2),xc(1),u0/); east=east/sqrt(dot_product(east,east))
   north=cross_product(xc,east)
   eano(:,1)=east; eano(:,2)=north
   xcd2=matmul(transpose(eano),xcd)
   dxt(lxm,iy)=sqrt(dot_product(xcd2(:,1),xcd2(:,1)))*delx
   gat(lxm,iy)=triple_product(xc,xcd(:,1),xcd(:,2))*delxy
enddo

!-- extra right edge, gat, dxt only:
xm(1)=mxp*delxore
do iy=ly,my
   xm(2)=iy*delyore
   call xmtoxc_ak(a,k,xm,xc,xcd,ff); if(ff)return
   xcd=matmul(prot,xcd)
   xc =matmul(prot,xc )
   east=(/-xc(2),xc(1),u0/); east=east/sqrt(dot_product(east,east))
   north=cross_product(xc,east)
   eano(:,1)=east; eano(:,2)=north
   xcd2=matmul(transpose(eano),xcd)
   dxt(mxp,iy)=sqrt(dot_product(xcd2(:,1),xcd2(:,1)))*delx
   gat(mxp,iy)=triple_product(xc,xcd(:,1),xcd(:,2))*delxy
enddo

!-- extra bottom edge, gat, dyt only:
xm(2)=lym*delyore
do ix=lx,mx
   xm(1)=ix*delxore
   call xmtoxc_ak(a,k,xm,xc,xcd,ff); if(ff)return
   xcd=matmul(prot,xcd)
   xc =matmul(prot,xc )
   east=(/-xc(2),xc(1),u0/); east=east/sqrt(dot_product(east,east))
   north=cross_product(xc,east)
   eano(:,1)=east; eano(:,2)=north
   xcd2=matmul(transpose(eano),xcd)
   dyt(ix,lym)=sqrt(dot_product(xcd2(:,2),xcd2(:,2)))*dely
   gat(ix,lym)=triple_product(xc,xcd(:,1),xcd(:,2))*delxy
enddo

!-- extra top edge, gat, dyt only:
xm(2)=myp*delyore
do ix=lx,mx
   xm(1)=ix*delxore
   call xmtoxc_ak(a,k,xm,xc,xcd,ff); if(ff)return
   xcd=matmul(prot,xcd)
   xc =matmul(prot,xc )
   east=(/-xc(2),xc(1),u0/); east=east/sqrt(dot_product(east,east))
   north=cross_product(xc,east)
   eano(:,1)=east; eano(:,2)=north
   xcd2=matmul(transpose(eano),xcd)
   dyt(ix,myp)=sqrt(dot_product(xcd2(:,2),xcd2(:,2)))*dely
   gat(ix,myp)=triple_product(xc,xcd(:,1),xcd(:,2))*delxy
enddo

! Extra four corners, gat only:
xm(2)=lym*delyore
!-- extra bottom left corner:
xm(1)=lxm*delxore
call xmtoxc_ak(a,k,xm,xc,xcd,ff); if(ff)return
xcd=matmul(prot,xcd)
xc =matmul(prot,xc )
gat(lxm,lym)=triple_product(xc,xcd(:,1),xcd(:,2))*delxy

!-- extra bottom right corner:
xm(1)=mxp*delxore
call xmtoxc_ak(a,k,xm,xc,xcd,ff); if(ff)return
xcd=matmul(prot,xcd)
xc =matmul(prot,xc )
gat(mxp,lym)=triple_product(xc,xcd(:,1),xcd(:,2))*delxy

xm(2)=myp*delyore
!-- extra top left corner:
xm(1)=lxm*delxore
call xmtoxc_ak(a,k,xm,xc,xcd,ff); if(ff)return
xcd=matmul(prot,xcd)
xc =matmul(prot,xc )
gat(lxm,myp)=triple_product(xc,xcd(:,1),xcd(:,2))*delxy

!-- extra top right corner:
xm(1)=mxp*delxore
call xmtoxc_ak(a,k,xm,xc,xcd,ff); if(ff)return
xcd=matmul(prot,xcd)
xc =matmul(prot,xc )
gat(mxp,myp)=triple_product(xc,xcd(:,1),xcd(:,2))*delxy

!-- 4th-order averaging over each central interval using 4-pt. stencils:
dx            =(13*(dxt(lx :mx-1,:)+dxt(lx+1:mx ,:)) &
                  -(dxt(lxm:mx-2,:)+dxt(lx+2:mxp,:)))/24
dy            =(13*(dyt(:,ly :my-1)+dyt(:,ly+1:my )) &
                  -(dyt(:,lym:my-2)+dyt(:,ly+2:myp)))/24
gat(lx:mx-1,:)=(13*(gat(lx :mx-1,:)+gat(lx+1:mx ,:)) &
                  -(gat(lxm:mx-2,:)+gat(lx+2:mxp,:)))/24
garea         =(13*(gat(lx:mx-1,ly :my-1)+gat(lx:mx-1,ly+1:my )) &
                  -(gat(lx:mx-1,lym:my-2)+gat(lx:mx-1,ly+2:myp)))/24
! Convert degrees to radians in the glat and glon arrays:
glat=glat*dtor
glon=glon*dtor

end subroutine hgrid_ak
