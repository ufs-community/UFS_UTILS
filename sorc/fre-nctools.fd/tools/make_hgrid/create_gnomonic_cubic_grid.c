/*******************************************************************************
                             create_gnomonic_cubic_grid.c
  This file creates the 6 tiles for a gnomonic projection of a cubed sphere.
    It also creates nest grids if they are defined.
  Modifications:
  05/10/2020  -- Added multiple nest capability.  Bill Ramstrom, AOML/HRD
                 Nests can be specified on any parent tile, and can each be different sizes.
                 Nests with different refinement ratios have NOT been tested and should be considered unsupported.

*******************************************************************************/



#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mpp.h"
#include "constant.h"
#include "mosaic_util.h"
#include "tool_util.h"
#include "create_hgrid.h"
#define  D2R (M_PI/180.)
#define  R2D (180./M_PI)
#define  EPSLN10 (1.e-10)
#define  EPSLN4 (1.e-4)
#define  EPSLN5 (1.e-5)
#define  EPSLN7 (1.e-7)
#define  EPSLN8  (1.0e-8)

/* private subroutines */
void gnomonic_ed  (int ni, double* lamda, double* theta);
void gnomonic_angl(int ni, double* lamda, double* theta);
void gnomonic_dist(int ni, double* lamda, double* theta);
void cartesian_to_spherical(double x, double y, double z, double *lon, double *lat, double *r);
void spherical_to_cartesian(double lon, double lat, double r, double *x, double *y, double *z);
void symm_ed(int ni, double *lamda, double *theta);
void mirror_grid(int ni, int ntiles, double *x, double *y );
void mirror_latlon(double lon1, double lat1, double lon2, double lat2, double lon0,
		   double lat0, double *lon, double *lat);
void rot_3d(int axis, double x1in, double y1in, double z1in, double angle, double *x2out,
	    double *y2out, double *z2out, int degrees, int convert);
double excess_of_quad2(const double *vec1, const double *vec2, const double *vec3, const double *vec4 );
double angle_between_vectors2(const double *vec1, const double *vec2);
void plane_normal2(const double *P1, const double *P2, double *plane);
void calc_rotation_angle2(int nxp, double *x, double *y, double *angle_dx, double *angle_dy);
void cell_center(int ni, int nj, const double *lonc, const double *latc, double *lont, double *latt);
void cell_east(int ni, int nj, const double *lonc, const double *latc, double *lone, double *late);
void cell_north(int ni, int nj, const double *lonc, const double *latc, double *lonn, double *latn);
void calc_cell_area(int nx, int ny, const double *x, const double *y, double *area);
void direct_transform(double stretch_factor, int i1, int i2, int j1, int j2, double lon_p, double lat_p,
		      int n, double *lon, double *lat);
void setup_aligned_nest(int parent_ni, int parent_nj, const double *parent_xc, const double *parent_yc,
                        int halo, int refine_ratio, int istart, int iend, int jstart, int jend,
                        double *xc, double *yc);

void spherical_linear_interpolation(double beta, const double *p1, const double *p2, double *pb);

/*******************************************************************************
  void create_gnomonic_cubic_grid( int *npoints, int *nratio, char *method, char *orientation, double *x,
                          double *y, double *dx, double *dy, double *area, double *angle_dx,
                          double *angle_dy )
  create nomomic cubic grid. All six tiles grid will be generated.
*******************************************************************************/
void create_gnomonic_cubic_grid( char* grid_type, int *nlon, int *nlat, double *x, double *y,
				 double *dx, double *dy, double *area, double *angle_dx,
				 double *angle_dy, double shift_fac, int do_schmidt, double stretch_factor,
				 double target_lon, double target_lat, int num_nest_grids,
				 int parent_tile[MAX_NESTS], int refine_ratio[MAX_NESTS], int istart_nest[MAX_NESTS],
				 int iend_nest[MAX_NESTS], int jstart_nest[MAX_NESTS], int jend_nest[MAX_NESTS], 
				 int halo)
{
  const int ntiles = 6;
  int verbose = 1;
  int ntiles2, global_nest=0;

  int nx, ny, nxp, nyp, ni, nj, nip, njp;

  int nx_nest[MAX_NESTS], ny_nest[MAX_NESTS];
  int ni_nest[MAX_NESTS], nj_nest[MAX_NESTS];
  int istart[MAX_NESTS], iend[MAX_NESTS], jstart[MAX_NESTS], jend[MAX_NESTS];

  int *nx_nest_arr=NULL;
  int *ny_nest_arr=NULL;
  int *ni_nest_arr=NULL;
  int *nj_nest_arr=NULL;


  int ni2, nj2, ni2p, nj2p, n1, n2;
  int *nxl=NULL, *nyl=NULL, *nil=NULL, *njl=NULL;
  int *tile_offset=NULL;
  int *tile_offset_supergrid=NULL;
  int *tile_offset_supergrid_m=NULL;
  int *tile_offset_area=NULL;

  int i, j, n, npts, nn;
  int npts_supergrid, npts_supergrid_m, npts_area;

  double p1[2], p2[2];
  double *lon=NULL, *lat=NULL;
  double *xc=NULL, *yc=NULL, *xtmp=NULL, *ytmp=NULL;
  double *xc2=NULL, *yc2=NULL;
  int    stretched_grid=0;

  /* make sure the first 6 tiles have the same grid size and 
     the size in x and y-direction are the same
  */


  /* ntiles is a constant always equal to 6.  ntiles2 is variable, and includes the 6 global tiles plus any nests */

  if (verbose) printf("[INFO] A.  Starting create_gnomonic_grid with num_nest_grids=%d\n", num_nest_grids);

  for(n=0; n<ntiles; n++) {
    if(nlon[n] != nlat[n] ) mpp_error("create_gnomonic_cubic_grid: the grid size in x and y-direction "
				      "should be the same for the 6 tiles of cubic sphere grid");
    if( nlon[n]%2 ) mpp_error("create_gnomonic_cubic_grid: supergrid size in x-direction should be divided by 2");
    if( nlat[n]%2 ) mpp_error("create_gnomonic_cubic_grid: supergrid size in y-direction should be divided by 2");
  }
  for(n=1; n<ntiles; n++) {
    if(nlon[n] != nlon[0]) mpp_error("create_gnomonic_cubic_grid: all six tiles should have same size");
  }

  // nx, ny variables correspond to the supergrid, which has twice as many points
  // ni, nj variables correspond to the number of cell centers
  // nip, njp variables correspond to the number of cell edges (i.e. ni+1, nj+1)


  nx  = nlon[0];
  ny  = nx;
  nxp = nx+1;
  nyp = ny+1;
  ni  = nx/2;
  nj  = ni;
  nip = ni+1;
  njp = nip;

  for (nn=0; nn < MAX_NESTS; nn++) {
    ni_nest[nn] = 0;
    nj_nest[nn] = 0;  
  }

  ntiles2=ntiles;
  global_nest=0;

  if (verbose) printf("[INFO] B\n");


  if(num_nest_grids && parent_tile[0]== 0)
    global_nest = 1;
  else {
    for (nn=0; nn < num_nest_grids; nn++) {
      ntiles2 = ntiles+num_nest_grids;
      if( (istart_nest[nn]+1)%2 ) mpp_error("create_gnomonic_cubic_grid: istart_nest+1 is not divisbile by 2");
      if( iend_nest[nn]%2 ) mpp_error("create_gnomonic_cubic_grid: iend_nest is not divisbile by 2");
      if( (jstart_nest[nn]+1)%2 ) mpp_error("create_gnomonic_cubic_grid: jstart_nest+1 is not divisbile by 2");
      if( jend_nest[nn]%2 ) mpp_error("create_gnomonic_cubic_grid: jend_nest is not divisbile by 2");  
      istart[nn] = (istart_nest[nn]+1)/2;
      iend[nn]   = iend_nest[nn]/2;
      jstart[nn] = (jstart_nest[nn]+1)/2;
      jend[nn]   = jend_nest[nn]/2;
      ni_nest[nn] = (iend[nn]-istart[nn]+1)*refine_ratio[nn];
      nj_nest[nn] = (jend[nn]-jstart[nn]+1)*refine_ratio[nn];  
      
      nx_nest[nn] = ni_nest[nn]*2;
      ny_nest[nn] = nj_nest[nn]*2;
    }
  }
  if (verbose) printf("[INFO] C\n");



  if (verbose) printf("[INFO] C1 ntiles: %d num_nest_grids: %d ntiles2: %d\n", ntiles, num_nest_grids, ntiles2);

  /*  Since many of the variables for the global and nest tiles are stored in 1D arrays, 
      we generate a set of indices that navigate to the start of each global and nest tile.
      ntiles is set to 1 for regional or 6 for global
      ntiles2 is ntiles + the number of nests

      nxl, nyl indicate supergrid size for each tile
      nil, njl indicate model grid size for each tile

   */

    

  nxl = (int *)malloc(ntiles2*sizeof(int));
  nyl = (int *)malloc(ntiles2*sizeof(int));
  nil = (int *)malloc(ntiles2*sizeof(int));
  njl = (int *)malloc(ntiles2*sizeof(int));

  
  nx_nest_arr = (int *)malloc(num_nest_grids*sizeof(int));
  ny_nest_arr = (int *)malloc(num_nest_grids*sizeof(int));
  ni_nest_arr = (int *)malloc(num_nest_grids*sizeof(int));
  nj_nest_arr = (int *)malloc(num_nest_grids*sizeof(int));

for (nn=0; nn < num_nest_grids; nn++) {
    nx_nest_arr[nn] = nx_nest[nn];
    ny_nest_arr[nn] = ny_nest[nn];
    ni_nest_arr[nn] = ni_nest[nn];
    nj_nest_arr[nn] = nj_nest[nn];
  }

  tile_offset = (int *)malloc(ntiles2*sizeof(int));
  tile_offset_supergrid = (int *)malloc(ntiles2*sizeof(int));
  tile_offset_supergrid_m = (int *)malloc(ntiles2*sizeof(int));
  tile_offset_area = (int *)malloc(ntiles2*sizeof(int));

  for(n=0; n<ntiles; n++) {
    nxl[n] = nx;
    nyl[n] = ny;
    nil[n] = ni;
    njl[n] = nj;
  }

  /* Use of these arrays permits different sized nests */
  for (nn=0; nn < num_nest_grids; nn++) {
    nxl[nn+ntiles] = nx_nest_arr[nn];
    nyl[nn+ntiles] = ny_nest_arr[nn];
    nil[nn+ntiles] = ni_nest_arr[nn];
    njl[nn+ntiles] = nj_nest_arr[nn];
  }


  if (verbose) {
    printf("[INFO] INDEX ntiles: %d ntiles2: %d\n", ntiles, ntiles2);
    for(n=0; n<ntiles2; n++) {
      printf("[INFO] INDEX n: %d nxl[n]: %d nyl[n]: %d nil[n]: %d njl[n]: %d\n",n, nxl[n], nyl[n], nil[n], njl[n]);
    }
  }
  
  if (verbose) printf("[INFO] D\n");

  /* for global nest grid, set ni to the coarse grid size */
  /* TODO -- can this code handle multiple different refinement ratios for global nests? */
  if(global_nest) {
    ni /= refine_ratio[0];
    nj /= refine_ratio[0];
  }
  nip=ni+1;
  njp=nj+1;
  
  if ( do_schmidt && fabs(stretch_factor-1.) > EPSLN5 ) stretched_grid = 1;
  
  lon = (double *)malloc(nip*nip*sizeof(double));
  lat = (double *)malloc(nip*nip*sizeof(double));
  
  if (verbose) printf("[INFO] E\n");
  
  if(strcmp(grid_type, "gnomonic_ed")==0 )
    gnomonic_ed(  ni, lon, lat);
  else if(strcmp(grid_type,"gnomonic_dist")==0)
    gnomonic_dist(ni, lon, lat);
  else if(strcmp(grid_type,"gnomonic_angl")==0)
    gnomonic_angl(ni, lon, lat);
  else mpp_error("create_gnomonic_cubic_grid: grid type should be 'gnomonic_ed', "
		 "'gnomonic_dist' or 'gnomonic_angl'");

  symm_ed(ni, lon, lat);

  
  if (verbose) printf("[INFO] F\n");
  
  // Cycle through all of the tiles; global and nests, adding enough points based on the dimensions
  // The 6 cubed-sphere tiles are square thus, nil=njl, but the nests can be rectangular
  npts = 0;
  npts_supergrid = 0;
  npts_supergrid_m = 0;
  npts_area = 0;

  for (n=0; n<ntiles2; n++) {
    tile_offset[n] = npts;
    tile_offset_supergrid[n] = npts_supergrid;
    tile_offset_supergrid_m[n] = npts_supergrid_m;
    tile_offset_area[n] = npts_area;

    if (verbose) printf("[INFO] INDEX OFFSET n: %d tile_offset[n]: %d tile_offset_supergrid[n]: %d tile_offset_supergrid_m[n]: %d tile_offset_area[n]: %d \n",
			n, tile_offset[n],  tile_offset_supergrid[n],  tile_offset_supergrid_m[n],  tile_offset_area[n]);

    npts += (nil[n] + 1) * (njl[n] + 1);
    npts_supergrid += (nxl[n] + 1) * (nyl[n] + 1);
    npts_supergrid_m += nxl[n] * (nyl[n] + 1);    // needed for grids of dx, dy
    npts_area += nxl[n] * nyl[n];    // needed for area
  }
  
  if (verbose) printf("[INFO] INDEX OFFSET npts: %d\n", npts);

  xc = (double *)malloc(npts*sizeof(double));
  yc = (double *)malloc(npts*sizeof(double));
  
  for(j=0; j<nip; j++) {
    for(i=0; i<nip; i++) {
      xc[j*nip+i] = lon[j*nip+i] - M_PI;
      yc[j*nip+i] = lat[j*nip+i];
    }
  }
  
  /* mirror_grid assumes that the tile=1 is centered on equator
     and greenwich meridian Lon[-pi,pi]  */
  mirror_grid(ni, ntiles, xc, yc);

  if (verbose) printf("[INFO] G\n");

  // Operate only on the 6 parent tiles
  for(n=0; n<ntiles*nip*nip; n++) {
    /* This will result in the corner close to east coast of china */
    if( do_schmidt == 0 && shift_fac > EPSLN4) xc[n] -= M_PI/18.;
    if(xc[n] < 0.) xc[n] += 2.*M_PI;
    if(fabs(xc[n]) < EPSLN10) xc[n] = 0;
    if(fabs(yc[n]) < EPSLN10) yc[n] = 0;
  }

  if (verbose) printf("[INFO] H\n");
      
  /* ensure consistency on the boundary between tiles */
  for(j=0; j<nip; j++) {
    xc[  nip*nip+j*nip] = xc[j*nip+ni];                 /* 1E -> 2W */
    yc[  nip*nip+j*nip] = yc[j*nip+ni];                 /* 1E -> 2W */
    xc[2*nip*nip+j*nip] = xc[ni*nip+ni-j];              /* 1N -> 3W */
    yc[2*nip*nip+j*nip] = yc[ni*nip+ni-j];              /* 1N -> 3W */      
  }
  for(i=0; i<nip; i++) {
    xc[4*nip*nip+ni*nip+i] = xc[(ni-i)*nip];            /* 1W -> 5N */
    yc[4*nip*nip+ni*nip+i] = yc[(ni-i)*nip];            /* 1W -> 2N */
    xc[5*nip*nip+ni*nip+i] = xc[i];                     /* 1S -> 6N */
    yc[5*nip*nip+ni*nip+i] = yc[i];                     /* 1S -> 6N */      
    xc[2*nip*nip+i]        = xc[nip*nip+ni*nip+i];      /* 2N -> 3S */
    yc[2*nip*nip+i]        = yc[nip*nip+ni*nip+i];      /* 2N -> 3S */
    xc[3*nip*nip+i]        = xc[nip*nip+(ni-i)*nip+ni];  /* 2E -> 4S */
    yc[3*nip*nip+i]        = yc[nip*nip+(ni-i)*nip+ni];  /* 2E -> 4S */      
  }
  for(j=0; j<nip; j++) {
    xc[5*nip*nip+j*nip+ni] = xc[nip*nip+ni-j];          /* 2S -> 6E */
    yc[5*nip*nip+j*nip+ni] = yc[nip*nip+ni-j];          /* 2S -> 6E */  
    xc[3*nip*nip+j*nip]    = xc[2*nip*nip+j*nip+ni];    /* 3E -> 4W */
    yc[3*nip*nip+j*nip]    = yc[2*nip*nip+j*nip+ni];    /* 3E -> 4W */
    xc[4*nip*nip+j*nip]    = xc[2*nip*nip+ni*nip+ni-j]; /* 3N -> 5W */
    yc[4*nip*nip+j*nip]    = yc[2*nip*nip+ni*nip+ni-j]; /* 3N -> 5W */
  }
  for(i=0; i<nip; i++) {
    xc[4*nip*nip+i] = xc[3*nip*nip+ni*nip+i];           /* 4N -> 5S */
    yc[4*nip*nip+i] = yc[3*nip*nip+ni*nip+i];           /* 4N -> 5S */
    xc[5*nip*nip+i] = xc[3*nip*nip+(ni-i)*nip+ni];      /* 4E -> 6S */
    yc[5*nip*nip+i] = yc[3*nip*nip+(ni-i)*nip+ni];      /* 4E -> 6S */
  }
  for(j=0; j<nip; j++) {
    xc[5*nip*nip+j*nip] = xc[4*nip*nip+j*nip+ni];    /* 5E -> 6W */
    yc[5*nip*nip+j*nip] = yc[4*nip*nip+j*nip+ni];    /* 5E -> 6W */  
  }

  if (verbose) printf("[INFO] I\n");


  /* Schmidt transformation */
  if ( do_schmidt ) {
    for(n=0; n<ntiles; n++) {
      
      if (verbose) {
	printf("[INFO] Calling direct_transform for tile %d\n", n);
      }

      direct_transform(stretch_factor, 0, ni, 0, ni, target_lon*D2R, target_lat*D2R,
		       n, xc+n*nip*nip, yc+n*nip*nip);
    }
  }

  if (verbose) printf("[INFO] J\n");

  
  /* get nest grid */
  if(global_nest) {
    npts = ntiles*nip*nip;
    xc2 = (double *)malloc(npts*sizeof(double));
    yc2 = (double *)malloc(npts*sizeof(double));  
    for(n=0; n<npts; n++) {
      xc2[n] = xc[n];
      yc2[n] = yc[n];
    }
    free(xc);
    free(yc);
    ni2  = ni;
    ni2p = nip;
    ni   = nx/2;
    nip  = ni + 1;
    npts = ntiles*nip*nip;
    xc = (double *)malloc(npts*sizeof(double));
    yc = (double *)malloc(npts*sizeof(double));  
    for(n=0; n<ntiles; n++) {
      printf("[INFO] calling setup_aligned_nest, n=%d\n",n);
      setup_aligned_nest(ni2, ni2, xc2+ni2p*ni2p*n, yc2+ni2p*ni2p*n, 0, refine_ratio[n],
      			 1, ni2, 1, ni2, xc+n*nip*nip, yc+n*nip*nip );      



    }
  }
  else if( num_nest_grids > 0 ) {
    for (nn=0; nn < num_nest_grids; nn++) {
      if (verbose) {
	printf("[INFO] Processing setup_aligned_nest for nest %d .  ntiles=%d parent_tile: %d\n", nn, ntiles, parent_tile);
      }
      
      /* Setup aligned nest -- final two arguments are memory locations for data to be returned */
      /* The pointer arithmetic is complicated */
      /* ni = number of points on supergrid */
      /* nip = ni + 1 */
      
      setup_aligned_nest(ni, ni, xc+tile_offset[parent_tile[nn]-1],
			 yc+tile_offset[parent_tile[nn]-1], halo, refine_ratio[nn],
			 istart[nn], iend[nn], jstart[nn], jend[nn],
			 xc+tile_offset[ntiles+nn], yc+tile_offset[ntiles+nn]);
    }
    
    
    if (verbose) {
      printf("[INFO] Completed processing setup_aligned_nest for nest(s)\n");
    }


  }

  if (verbose) printf("[INFO] K\n");

  
  /* calculate grid box center location */

  ni2 = 0;
  nj2 = 0;
  for(n=0; n<ntiles2; n++) {
    if(nil[n]>ni2) ni2 = nil[n];
    if(njl[n]>nj2) nj2 = njl[n];
  }
  ni2p = ni2+1;
  nj2p = nj2+1;
  xtmp = (double *)malloc(ni2p*nj2p*sizeof(double));
  ytmp = (double *)malloc(ni2p*nj2p*sizeof(double));

  if (verbose) printf("[INFO] L\n");

  /* Setting the x, y values for each tile */
  /* Not clear that data is handled correctly for nested tiles, though. */

  /* Iterate over all of the tiles                                                                   */
  /* Copy the lat/lons into the x, y array                                                           */
  /*     C-cell                                                                                      */
  /*     Center                                                                                      */
  /*     East                                                                                        */
  /*     North                                                                                       */
      

  for(n=0; n<ntiles2; n++) {
    int min_n1 = -1;
    int max_n1 = -1;
    

    /* copy C-cell to supergrid */
    if (verbose) printf("[INFO] INDEX fill x and y from C-cell.  n=%d n*nxp*nxp=%d tile_offset[n]: %d tile_offset_supergrid[n]: %d njl[n]: %d nil[n]: %d\n",n,  n*nxp*nxp, tile_offset[n], tile_offset_supergrid[n], njl[n], nil[n]);
    if (verbose) printf("[INFO] START fill x and y from C-cell.  n=%d tile_offset_supergrid[n]: %d \n",n, tile_offset_supergrid[n]);

    for(j=0; j<=njl[n]; j++) for(i=0; i<=nil[n]; i++) {
	/* n1 = n*nxp*nxp+j*2*(2*nil[n]+1)+i*2; */
	//n1 = n*nxp*nxp+j*2*(2*nil[n]+1)+i*2;
	//n2 = n*nip*nip+j*(nil[n]+1)+i;
	n1 = tile_offset_supergrid[n] + j*2*(2*nil[n]+1) + i*2;
	n2 = tile_offset[n] + j*(nil[n]+1)+i;
	//if (verbose) printf("[INFO] CCELL n: %d i: %d j: %d n1: %d npts_supergrid: %d n2: %d npts: %d\n", n, i, j, n1, npts_supergrid, n2, npts);
	x[n1]=xc[n2];
	y[n1]=yc[n2];

	if (verbose){ 
	  if (n1 < min_n1 || min_n1 == -1) min_n1 = n1;
	  if (n1 > max_n1) max_n1 = n1;
	}



      }
    
    /* cell center and copy to super grid */
    //cell_center(nil[n], njl[n], xc+n*nip*nip, yc+n*nip*nip, xtmp, ytmp);
    cell_center(nil[n], njl[n], xc + tile_offset[n], yc + tile_offset[n], xtmp, ytmp);
    if (verbose) printf("[INFO] CENTER n: %d n*nip*nip: %d tile_offset[n]: %d\n", n, n*nip*nip, tile_offset[n]);
    for(j=0; j<njl[n]; j++) for(i=0; i<nil[n]; i++) {
	// Offset of 2 for i=0, j=0
	//n1 = n*nxp*nxp+(j*2+1)*(2*nil[n]+1)+i*2+1;
	//n2 = j*nil[n]+i;
	n1 = tile_offset_supergrid[n] + (j*2+1)*(2*nil[n]+1)+i*2+1;
	n2 = j*nil[n]+i;   // WDR why does this not have a tile_offset??
	x[n1]=xtmp[n2];
	y[n1]=ytmp[n2];

	if (verbose){ 
	  if (n1 < min_n1) min_n1 = n1;
	  if (n1 > max_n1) max_n1 = n1;
	}

      }

    /* cell east and copy to super grid */
    //cell_east(nil[n], njl[n], xc+n*nip*nip, yc+n*nip*nip, xtmp, ytmp);
    cell_east(nil[n], njl[n], xc + tile_offset[n], yc + tile_offset[n], xtmp, ytmp);
    for(j=0; j<njl[n]; j++) for(i=0; i<=nil[n]; i++) {
	// Offset of 2*nil[n] + 1 for i=0, j=0
	//n1 = n*nxp*nxp+(j*2+1)*(2*nil[n]+1)+i*2;
	n1 = tile_offset_supergrid[n] + (j*2+1)*(2*nil[n]+1)+i*2;
	n2 = j*(nil[n]+1)+i;   // WDR why does this not have a tile_offset??
	x[n1]=xtmp[n2];
	y[n1]=ytmp[n2];

	if (verbose){ 
	  if (n1 < min_n1) min_n1 = n1;
	  if (n1 > max_n1) max_n1 = n1;
	}


      }

    /* cell north and copy to super grid */
    //cell_north(nil[n], njl[n], xc+n*nip*nip, yc+n*nip*nip, xtmp, ytmp);
    cell_north(nil[n], njl[n], xc + tile_offset[n], yc + tile_offset[n], xtmp, ytmp);
    for(j=0; j<=njl[n]; j++) for(i=0; i<nil[n]; i++) {
	// Offset of 1 for i=0, j=0
	//n1 = n*nxp*nxp+(j*2)*(2*nil[n]+1)+i*2+1;
	n1 = tile_offset_supergrid[n] + (j*2)*(2*nil[n]+1)+i*2+1;
	n2 = j*nil[n]+i;   // WDR why does this not have a tile_offset??
	x[n1]=xtmp[n2];
	y[n1]=ytmp[n2];

	if (verbose){ 
	  if (n1 < min_n1) min_n1 = n1;
	  if (n1 > max_n1) max_n1 = n1;
	}

      }

    if (verbose) printf("[INFO] INDEX tile: %d min_n1: %d max_n1: %d max_n1 - min_n1: %d sqrt(max_n1 - min_n1 + 1): %f\n", n, min_n1, max_n1, max_n1 - min_n1, sqrt(max_n1 - min_n1 + 1));

  } 

  free(xtmp);
  free(ytmp);

  if (verbose) printf("[INFO] M\n");
  
  /* calculate grid cell length */

  /* Calculate dx */
  for(n=0; n<ntiles2; n++) {
    if (verbose) printf("[INFO] M1 Calculating dx for tile n: %d ntiles2: %d\n", n, ntiles2);
    
    for(j=0; j<=nyl[n]; j++) {
      for(i=0; i<nxl[n]; i++) {

	p1[0] = x[tile_offset_supergrid[n] + j*(nxl[n]+1)+i];
	p1[1] = y[tile_offset_supergrid[n] + j*(nxl[n]+1)+i];
	p2[0] = x[tile_offset_supergrid[n] + j*(nxl[n]+1)+i+1];
	p2[1] = y[tile_offset_supergrid[n] + j*(nxl[n]+1)+i+1];

	dx[tile_offset_supergrid_m[n] + j*nxl[n]+i] = great_circle_distance(p1, p2);

      }
    }
  }

  /* Calculate dy */
  for(n=0; n<ntiles2; n++) {
    if (verbose) printf("[INFO] M2a Calculating dy for tile n: %d ntiles: %d ntiles2: %d\n", n, ntiles, ntiles2);

    //if( stretched_grid || n==ntiles ) { 
    if( stretched_grid || (n >= 6) ) { 
      if (verbose) printf("[INFO] M2b Calculating dy for tile n: %d ntiles: %d ntiles2: %d\n", n, ntiles, ntiles2);
      for(j=0; j<nyl[n]; j++) {
	for(i=0; i<=nxl[n]; i++) {
	  p1[0] = x[tile_offset_supergrid[n] + j*(nxl[n]+1)+i];
	  p1[1] = y[tile_offset_supergrid[n] + j*(nxl[n]+1)+i];
	  p2[0] = x[tile_offset_supergrid[n] + (j+1)*(nxl[n]+1)+i];
	  p2[1] = y[tile_offset_supergrid[n] + (j+1)*(nxl[n]+1)+i];

	  dy[tile_offset_supergrid_m[n] + j*(nxl[n]+1)+i] = great_circle_distance(p1, p2);
	}
      }
    }
    else {
      if (verbose) printf("[INFO] M2c Calculating dy for tile n: %d ntiles: %d ntiles2: %d\n", n, ntiles, ntiles2);
      for(j=0; j<nyp; j++) {
	for(i=0; i<nx; i++) dy[tile_offset_supergrid_m[n] + i*nxp+j] = dx[tile_offset_supergrid_m[n] + j*nx+i];
      }
    }
  }

  if (verbose) printf("[INFO] N\n");

  
  /* ensure consistency on the boundaries between tiles */
  for(j=0; j<nx; j++) {
    dy[j*nxp]             = dx[4*nx*nxp+nx*nx+nx-j-1]; /* 5N -> 1W */
    dy[j*nxp+nx]          = dy[nxp*nx+j*nxp];          /* 2W -> 1E */
    dy[nxp*nx+j*nxp+nx]   = dx[3*nx*nxp+(nx-j-1)];     /* 4S -> 2E */
    dy[2*nxp*nx+j*nxp]    = dx[nx*nx+nx-j-1];          /* 1N -> 3W */
    dy[2*nxp*nx+j*nxp+nx] = dy[3*nxp*nx+j*nxp];        /* 4W -> 3E */
    dy[3*nxp*nx+j*nxp+nx] = dx[5*nx*nxp+(nx-j-1)];     /* 4S -> 2E */
    dy[4*nxp*nx+j*nxp]    = dx[2*nx*nxp+nx*nx+nx-j-1]; /* 3N -> 5W */
    dy[4*nxp*nx+j*nxp+nx] = dy[5*nxp*nx+j*nxp];        /* 6W -> 5E */
    dy[5*nxp*nx+j*nxp+nx] = dx[nx*nxp+(nx-j-1)];       /* 2S -> 6E */    
  }

  if (verbose) printf("[INFO] O\n");


  if(do_schmidt) { /* calculate area for each tile */
    for(n=0; n<ntiles; n++) {
      if (verbose) printf("[INFO] call calc_cell_area do_schmidt for tile n=%d\n", n);
      calc_cell_area(nx, ny, x + tile_offset_supergrid[n], y + tile_offset_supergrid[n], area + tile_offset_area[n]);
    }
  }
  else  {
    if (verbose) printf("[INFO] call calc_cell_area for first tile.\n");
    calc_cell_area(nx, ny, x, y, area);
    for(j=0; j<nx; j++) {
      for(i=0; i<nx; i++) {
	double ar;
	/* all the faces have the same area */
	ar = area[j*nx+i];
	area[nx*nx+j*nx+i] = ar;
	area[2*nx*nx+j*nx+i] = ar;
	area[3*nx*nx+j*nx+i] = ar;
	area[4*nx*nx+j*nx+i] = ar;
	area[5*nx*nx+j*nx+i] = ar;        
      }
    }
  
  }
  if (verbose) printf("[INFO] P\n");


  for (nn=0; nn < num_nest_grids; nn++) {
    if (verbose) printf("[INFO] call calc_cell_area  for nest nn=%d tile n=%d\n", nn, ntiles+nn);
    if (verbose) printf("[INFO] Q1\n");

      calc_cell_area(nx_nest_arr[nn], ny_nest_arr[nn], x + tile_offset_supergrid[ntiles+nn], y + tile_offset_supergrid[ntiles+nn], area + tile_offset_area[ntiles+nn]);
  }

  if (verbose) printf("[INFO] Q2\n");

  
  /*calculate rotation angle, just some workaround, will modify this in the future. */
  calc_rotation_angle2(nxp, x, y, angle_dx, angle_dy );

  if (verbose) printf("[INFO] R\n");


  /* since angle is used in the model, set angle to 0 for nested region */
  for(nn=0; nn < num_nest_grids; nn++) {
    for(i=0; i<=(nx_nest_arr[nn]+1)*(ny_nest_arr[nn]+1); i++) {
      angle_dx[tile_offset_supergrid[ntiles+nn] + i]=0;
      angle_dy[tile_offset_supergrid[ntiles+nn] + i]=0;
    }
  }
      
    
  /* convert grid location from radians to degree */
  if (verbose) printf("[INFO] Convert radians to degrees: npts = %d npts_supergrid: %d\n", npts, npts_supergrid);
  for(i=0; i<npts_supergrid; i++) {
    x[i] = x[i]*R2D;
    y[i] = y[i]*R2D;
  }

  free(xc);
  free(yc);  
  
  if (verbose) printf("[INFO] Z\n");

  
}; /* void create_gnomonic_cubic_grid */

void calc_cell_area(int nx, int ny, const double *x, const double *y, double *area)
{
  int i,j, nxp;
  double p_ll[2], p_ul[2], p_lr[2], p_ur[2];

  nxp = nx+1;
  for(j=0; j<ny; j++) {
    for(i=0; i<nx; i++) {
      //printf("[INFO] calc_cell_area: i: %d j: %d\n", i, j);
      p_ll[0] = x[j*nxp+i];       p_ll[1] = y[j*nxp+i];
      p_ul[0] = x[(j+1)*nxp+i];   p_ul[1] = y[(j+1)*nxp+i];
      p_lr[0] = x[j*nxp+i+1];     p_lr[1] = y[j*nxp+i+1];
      p_ur[0] = x[(j+1)*nxp+i+1]; p_ur[1] = y[(j+1)*nxp+i+1];
      /* all the face have the same area */
      area[j*nx+i] = spherical_excess_area(p_ll, p_ul, p_lr, p_ur, RADIUS);
    }
  }

}



/*-------------------------------------------------------------------------
  void direct_transform(double c, int i1, int i2, int j1, int j2, double lon_p, double lat_p, int n,
  double *lon, double *lat)

  This is a direct transformation of the standard (symmetrical) cubic grid
  to a locally enhanced high-res grid on the sphere; it is an application
  of the Schmidt transformation at the south pole followed by a 
  pole_shift_to_target (rotation) operation

  arguments:
  c            : Stretching factor
  lon_p, lat_p : center location of the target face, radian
  n            : grid face number
  i1,i2,j1,j2  : starting and ending index in i- and j-direction
  lon          : longitude. 0 <= lon <= 2*pi
  lat          : latitude. -pi/2 <= lat <= pi/2
  ------------------------------------------------------------------------*/

void direct_transform(double stretch_factor, int i1, int i2, int j1, int j2, double lon_p, double lat_p,
		      int n, double *lon, double *lat)
{
#ifdef NO_QUAD_PRECISION
  double lat_t, sin_p, cos_p, sin_lat, cos_lat, sin_o, p2, two_pi;
  double c2p1, c2m1;
#else
  long double lat_t, sin_p, cos_p, sin_lat, cos_lat, sin_o, p2, two_pi;
  long double c2p1, c2m1;
#endif
  int i, j, l, nxp;

  nxp = i2-i1+1;
  p2 = 0.5*M_PI;
  two_pi = 2.*M_PI;
  if(n==0) printf("create_gnomonic_cubic_grid: Schmidt transformation: stretching factor=%g, center=(%g,%g)\n",
		  stretch_factor, lon_p, lat_p);

  c2p1 = 1. + stretch_factor*stretch_factor;
  c2m1 = 1. - stretch_factor*stretch_factor;

  sin_p = sin(lat_p);
  cos_p = cos(lat_p);

  for(j=j1; j<=j2; j++) for(i=i1; i<=i2; i++) {
      l = j*nxp+i;
      if ( fabs(c2m1) > EPSLN7 ) {
	sin_lat = sin(lat[l]); 
	lat_t   = asin( (c2m1+c2p1*sin_lat)/(c2p1+c2m1*sin_lat) );
      }
      else {
	lat_t = lat[l];
      }
      sin_lat = sin(lat_t);
      cos_lat = cos(lat_t); 
      sin_o = -(sin_p*sin_lat + cos_p*cos_lat*cos(lon[l]));
      if ( (1.-fabs(sin_o)) < EPSLN7 ) {    /* poles */
	lon[l] = 0.;
	lat[l] = (sin_o < 0) ? -p2:p2;
      }
      else {
	lat[l] = asin( sin_o );
	lon[l] = lon_p + atan2(-cos_lat*sin(lon[l]), -sin_lat*cos_p+cos_lat*sin_p*cos(lon[l]));
	if ( lon[l] < 0. )
	  lon[l] +=two_pi;
	else if( lon[l] >= two_pi )
	  lon[l] -=two_pi;
      }
    }
  
}; /* direct_transform */



/*-----------------------------------------------------
  void gnomonic_ed
  Equal distance along the 4 edges of the cubed sphere
  -----------------------------------------------------
  Properties: 
  * defined by intersections of great circles
  * max(dx,dy; global) / min(dx,dy; global) = sqrt(2) = 1.4142
  * Max(aspect ratio) = 1.06089
  * the N-S coordinate curves are const longitude on the 4 faces with equator 
  For C2000: (dx_min, dx_max) = (3.921, 5.545)    in km unit
  ! Ranges:
  ! lamda = [0.75*pi, 1.25*pi]
  ! theta = [-alpha, alpha]
  --------------------------------------------------------*/
void gnomonic_ed(int ni, double* lamda, double* theta)
{

  int i, j, n, nip;
  double dely;
  double *pp, *x, *y, *z;
  double rsq3, alpha;


  nip = ni + 1;
  rsq3 = 1./sqrt(3.);
  alpha = asin( rsq3 );

  dely = 2.*alpha/ni;

  /* Define East-West edges: */
  for(j=0; j<nip; j++) {
    lamda[j*nip]    = 0.75*M_PI;               /* West edge */
    lamda[j*nip+ni] = 1.25*M_PI;               /* East edge */
    theta[j*nip]    = -alpha + dely*j;       /* West edge */
    theta[j*nip+ni] = theta[j*nip];          /* East edge */
  }

  /* Get North-South edges by symmetry: */

  for(i=1; i<ni; i++) {
    mirror_latlon(lamda[0], theta[0], lamda[ni*nip+ni], theta[ni*nip+ni],  
		  lamda[i*nip], theta[i*nip], &lamda[i], &theta[i] );
    lamda[ni*nip+i] = lamda[i];
    theta[ni*nip+i] = -theta[i];
  }

  x = (double *)malloc(nip*nip*sizeof(double));
  y = (double *)malloc(nip*nip*sizeof(double));
  z = (double *)malloc(nip*nip*sizeof(double));
  /* Set 4 corners: */
  latlon2xyz(1, &lamda[0], &theta[0], &x[0], &y[0], &z[0]);
  latlon2xyz(1, &lamda[ni], &theta[ni], &x[ni], &y[ni], &z[ni]);
  latlon2xyz(1, &lamda[ni*nip], &theta[ni*nip], &x[ni*nip], &y[ni*nip], &z[ni*nip]);
  latlon2xyz(1, &lamda[ni*nip+ni], &theta[ni*nip+ni], &x[ni*nip+ni], &y[ni*nip+ni], &z[ni*nip+ni]);

  /* Map edges on the sphere back to cube:
     Intersections at x=-rsq3   */

  for(j=1; j<ni; j++) {
    n = j*nip;
    latlon2xyz(1, &lamda[n], &theta[n], &x[n], &y[n], &z[n]);
    y[n] = -y[n]*rsq3/x[n];
    z[n] = -z[n]*rsq3/x[n];
  }

  for(i=1; i<ni; i++) {
    latlon2xyz(1, &lamda[i], &theta[i], &x[i], &y[i], &z[i]);
    y[i] = -y[i]*rsq3/x[i];
    z[i] = -z[i]*rsq3/x[i];
  }    

  for(j=0; j<nip; j++)
    for(i=0; i<nip; i++) x[j*nip+i] = -rsq3;

  for(j=1;j<nip; j++) {
    for(i=1; i<nip; i++) {
      y[j*nip+i] = y[i];
      z[j*nip+i] = z[j*nip];
    }
  }

  xyz2latlon(nip*nip, x, y, z, lamda, theta);

}; /* gnomonic_ed */

/*----------------------------------------------------------
  void gnomonic_angl()
  This is the commonly known equi-angular grid
**************************************************************/

void gnomonic_angl(int ni, double* lamda, double* theta)
{



}; /* gnomonic_angl */

/*----------------------------------------------------------
  void gnomonic_dist()
  This is the commonly known equi-distance grid
**************************************************************/

void gnomonic_dist(int ni, double* lamda, double* theta)
{



}; /* gnomonic_dist */

/*------------------------------------------------------------------
  void mirror_latlon
  Given the "mirror" as defined by (lon1, lat1), (lon2, lat2), and center 
  of the sphere, compute the mirror image of (lon0, lat0) as  (lon, lat)
  ---------------------------------------------------------------*/
   
void mirror_latlon(double lon1, double lat1, double lon2, double lat2, double lon0,
		   double lat0, double *lon, double *lat)
{
  double p0[3], p1[3], p2[3], pp[3], nb[3];
  double pdot;
  int k;

  latlon2xyz(1, &lon0, &lat0, &p0[0], &p0[1], &p0[2]);
  latlon2xyz(1, &lon1, &lat1, &p1[0], &p1[1], &p1[2]);
  latlon2xyz(1, &lon2, &lat2, &p2[0], &p2[1], &p2[2]);
  vect_cross(p1, p2, nb);
     
  pdot = sqrt(nb[0]*nb[0]+nb[1]*nb[1]+nb[2]*nb[2]);
  for(k=0; k<3; k++) nb[k] = nb[k]/pdot;

  pdot = p0[0]*nb[0] + p0[1]*nb[1] + p0[2]*nb[2];
  for(k=0; k<3; k++) pp[k] = p0[k] - 2*pdot*nb[k];
  xyz2latlon(1, &pp[0], &pp[1], &pp[2], lon, lat);
    
}; /* mirror_latlon */

/*-------------------------------------------------------------------------
  void symm_ed(int ni, double *lamda, double *theta)
  ! Make grid symmetrical to i=ni/2+1 and j=nj/2+1
  ------------------------------------------------------------------------*/
void symm_ed(int ni, double *lamda, double *theta)
{
  
  int nip, i, j, ip, jp;
  double avg;
  
  nip = ni+1;

  for(j=1; j<nip; j++)
    for(i=1; i<ni; i++) lamda[j*nip+i] = lamda[i];
  
  for(j=0; j<nip; j++) {
    for(i=0; i<ni/2; i++) {
      ip = ni - i;
      avg = 0.5*(lamda[j*nip+i]-lamda[j*nip+ip]);
      lamda[j*nip+i] = avg + M_PI;
      lamda[j*nip+ip] = M_PI - avg;
      avg = 0.5*(theta[j*nip+i]+theta[j*nip+ip]);
      theta[j*nip+i] = avg;
      theta[j*nip+ip] = avg;      
    }
  }

  /* Make grid symmetrical to j=im/2+1 */
  for(j = 0; j<ni/2; j++) {
    jp = ni - j;
    for(i=1; i<ni; i++) {
      avg = 0.5*(lamda[j*nip+i]+lamda[jp*nip+i]);
      lamda[j*nip+i] = avg;
      lamda[jp*nip+i] = avg;
      avg = 0.5*(theta[j*nip+i]-theta[jp*nip+i]);
      theta[j*nip+i] = avg;
      theta[jp*nip+i] = -avg;
    }
  }
}/* symm_ed */

/*------------------------------------------------------------------------------
  void mirror_grid( )
  Mirror Across the 0-longitude
  ----------------------------------------------------------------------------*/
void mirror_grid(int ni, int ntiles, double *x, double *y )
{
  int nip, i, j, ip, jp, nt;
  double x1, y1, z1, x2, y2, z2, ang;

  nip = ni+1;
  
  for(j=0; j<ceil(nip/2.); j++) {
    jp = ni - j;
    for(i=0; i<ceil(nip/2.); i++) {
      ip = ni - i;
      x1 = 0.25 * (fabs(x[j*nip+i]) + fabs(x[j*nip+ip]) + fabs(x[jp*nip+i]) + fabs(x[jp*nip+ip]) );
      x[j*nip+i]   = x1 * (x[j*nip+i]   >=0 ? 1:-1);
      x[j*nip+ip]  = x1 * (x[j*nip+ip]  >=0 ? 1:-1);
      x[jp*nip+i]  = x1 * (x[jp*nip+i]  >=0 ? 1:-1);
      x[jp*nip+ip] = x1 * (x[jp*nip+ip] >=0 ? 1:-1);      

      y1 = 0.25 * (fabs(y[j*nip+i]) + fabs(y[j*nip+ip]) + fabs(y[jp*nip+i]) + fabs(y[jp*nip+ip]) );
      y[j*nip+i]   = y1 * (y[j*nip+i]   >=0 ? 1:-1);
      y[j*nip+ip]  = y1 * (y[j*nip+ip]  >=0 ? 1:-1);
      y[jp*nip+i]  = y1 * (y[jp*nip+i]  >=0 ? 1:-1);
      y[jp*nip+ip] = y1 * (y[jp*nip+ip] >=0 ? 1:-1);      
      
      /* force dateline/greenwich-meridion consitency */
      if( nip%2 ) {
	if( i == (nip-1)/2 ) {
	  x[j*nip+i] = 0.0;
	  x[jp*nip+i] = 0.0;
	}
      }
    }
  }

  /* define the the other five tiles. */
  for(nt=1; nt<ntiles; nt++) {
    for(j=0; j<nip; j++) {
      for(i=0; i<nip; i++) {
	x1 = x[j*nip+i];
	y1 = y[j*nip+i];
	z1 = RADIUS;
	switch (nt) {
	case 1: /* tile 2 */
	  ang = -90.;
	  rot_3d( 3, x1, y1, z1, ang, &x2, &y2, &z2, 1, 1);  /* rotate about the z-axis */
	  break;
	case 2: /* tile 3 */
	  ang = -90.;
	  rot_3d( 3, x1, y1, z1, ang, &x2, &y2, &z2, 1, 1);  /* rotate about the z-axis */
	  ang = 90.;
	  rot_3d( 1, x2, y2, z2, ang, &x1, &y1, &z1, 1, 1); /* rotate about the z-axis */
	  x2=x1;
	  y2=y1;
	  z2=z1;

	  /* force North Pole and dateline/greenwich-meridion consitency */
	  if(nip%2) {
	    if( (i==(nip-1)/2) && (i==j) ) {
	      x2 = 0;
	      y2 = M_PI*0.5;
	    }

	    if( (j==(nip-1)/2) && (i<(nip-1)/2) ) x2 = 0;
	    if( (j==(nip-1)/2) && (i>(nip-1)/2) ) x2 = M_PI;
	  }
	  break;
	case 3: /* tile 4 */
	  ang = -180.;
	  rot_3d( 3, x1, y1, z1, ang, &x2, &y2, &z2, 1, 1); /* rotate about the z-axis */
	  ang = 90.;
	  rot_3d( 1, x2, y2, z2, ang, &x1, &y1, &z1, 1, 1); /* rotate about the z-axis */
	  x2=x1;
	  y2=y1;
	  z2=z1;

	  /* force dateline/greenwich-meridion consitency */
	  if( nip%2 ) {
	    if( j == (nip-1)/2 ) x2 = M_PI;
	  }
	  break;
	case 4: /* tile 5 */
	  ang = 90.;
	  rot_3d( 3, x1, y1, z1, ang, &x2, &y2, &z2, 1, 1); /* rotate about the z-axis */
	  ang = 90.;
	  rot_3d( 2, x2, y2, z2, ang, &x1, &y1, &z1, 1, 1); /* rotate about the z-axis */
	  x2=x1;
	  y2=y1;
	  z2=z1;
	  break;
	case 5: /* tile 6 */
	  ang = 90.;
	  rot_3d( 2, x1, y1, z1, ang, &x2, &y2, &z2, 1, 1); /* rotate about the z-axis */
	  ang = 0.;
	  rot_3d( 3, x2, y2, z2, ang, &x1, &y1, &z1, 1, 1); /* rotate about the z-axis */
	  x2=x1;
	  y2=y1;
	  z2=z1;
	  
	  /* force South Pole and dateline/greenwich-meridion consitency */
	  if(nip%2) {
	    if( (i==(nip-1)/2) && (i==j) ) {
	      x2 = 0;
	      y2 = -M_PI*0.5;
	    }

	    if( (i==(nip-1)/2) && (j>(nip-1)/2) ) x2 = 0;
	    if( (i==(nip-1)/2) && (j<(nip-1)/2) ) x2 = M_PI;
	  }
	  break;
	}
	x[nt*nip*nip+j*nip+i] = x2;
	y[nt*nip*nip+j*nip+i] = y2;
      }
    }
  }
}; /* mirror_grid */


/*-------------------------------------------------------------------------------
  void rot_3d()
  rotate points on a sphere in xyz coords (convert angle from
  degrees to radians if necessary)
  -----------------------------------------------------------------------------*/
void rot_3d(int axis, double x1in, double y1in, double z1in, double angle, double *x2out,
	    double *y2out, double *z2out, int degrees, int convert)
{

  double x1, y1, z1, x2, y2, z2, c, s;
  
  if(convert)
    spherical_to_cartesian(x1in, y1in, z1in, &x1, &y1, &z1);
  else {
    x1=x1in;
    y1=y1in;
    z1=z1in;
  }

  if(degrees) angle = angle*D2R;

  c = cos(angle);
  s = sin(angle);

  switch (axis) {
  case 1:
    x2 =  x1;
    y2 =  c*y1 + s*z1;
    z2 = -s*y1 + c*z1;
    break;
  case 2:
    x2 = c*x1 - s*z1;
    y2 = y1;
    z2 = s*x1 + c*z1;
    break;
  case 3:
    x2 =  c*x1 + s*y1;
    y2 = -s*x1 + c*y1;
    z2 = z1;
    break;
  default:
    mpp_error("Invalid axis: must be 1 for X, 2 for Y, 3 for Z.");
  }
  
  if(convert)
    cartesian_to_spherical(x2, y2, z2, x2out, y2out, z2out);
  else {
    *x2out=x2;;
    *y2out=y2;
    *z2out=z2;
  }
} /* rot_3d */

/*-------------------------------------------------------------
  void cartesian_to_spherical(x, y, z, lon, lat, r)
  may merge with xyz2latlon in the future
  ------------------------------------------------------------*/
void cartesian_to_spherical(double x, double y, double z, double *lon, double *lat, double *r)
{

  *r = sqrt(x*x + y*y + z*z);
  if ( (fabs(x) + fabs(y)) < EPSLN10 )       /* poles */
    *lon = 0.;
  else
    *lon = atan2(y,x);    /* range: [-pi,pi] */


  *lat = acos(z/(*r)) - M_PI/2.;
};/* cartesian_to_spherical */

/*-------------------------------------------------------------------------------
  void spherical_to_cartesian
  convert from spheircal coordinates to xyz coords
  may merge with latlon2xyz in the future
  -----------------------------------------------------------------------------*/
void spherical_to_cartesian(double lon, double lat, double r, double *x, double *y, double *z)
{
  *x = r * cos(lon) * cos(lat);
  *y = r * sin(lon) * cos(lat);

  *z = -r * sin(lat);
} /* spherical_to_cartesian */


/*****************************************************************
   double* excess_of_quad(int ni, int nj, double *vec1, double *vec2, 
                          double *vec3, double *vec4 )
*******************************************************************/
double excess_of_quad2(const double *vec1, const double *vec2, const double *vec3, const double *vec4 )
{
  double plane1[3], plane2[3], plane3[3], plane4[3];
  double angle12, angle23, angle34, angle41, excess;
  double ang12, ang23, ang34, ang41;
  
  plane_normal2(vec1, vec2, plane1);
  plane_normal2(vec2, vec3, plane2);
  plane_normal2(vec3, vec4, plane3);
  plane_normal2(vec4, vec1, plane4);
  angle12 = angle_between_vectors2(plane2,plane1);
  angle23 = angle_between_vectors2(plane3,plane2);
  angle34 = angle_between_vectors2(plane4,plane3);
  angle41 = angle_between_vectors2(plane1,plane4);
  ang12 = M_PI-angle12;
  ang23 = M_PI-angle23;
  ang34 = M_PI-angle34;
  ang41 = M_PI-angle41;
  excess = ang12+ang23+ang34+ang41-2*M_PI;
  /*  excess = 2*M_PI - angle12 - angle23 - angle34 - angle41; */

  return excess;

}; /* excess_of_quad */

/******************************************************************************* 
double angle_between_vectors(const double *vec1, const double *vec2)
*******************************************************************************/

double angle_between_vectors2(const double *vec1, const double *vec2) {
  int n;
  double vector_prod, nrm1, nrm2;
  double angle;
  
  vector_prod=vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2];
  nrm1=pow(vec1[0],2)+pow(vec1[1],2)+pow(vec1[2],2);
  nrm2=pow(vec2[0],2)+pow(vec2[1],2)+pow(vec2[2],2);
  if(nrm1*nrm2>0)
    angle = acos( vector_prod/sqrt(nrm1*nrm2) );
  else
    angle = 0;
  
  return angle;
}; /* angle_between_vectors */


/***********************************************************************
   double* plane_normal(int ni, int nj, double *P1, double *P2)
***********************************************************************/

void plane_normal2(const double *P1, const double *P2, double *plane)
{
  double mag;
  
  plane[0] = P1[1] * P2[2] - P1[2] * P2[1];
  plane[1] = P1[2] * P2[0] - P1[0] * P2[2];
  plane[2] = P1[0] * P2[1] - P1[1] * P2[0];
  mag=sqrt(pow(plane[0],2) + pow(plane[1],2) + pow(plane[2],2));
  if(mag>0) {
    plane[0]=plane[0]/mag;
    plane[1]=plane[1]/mag;
    plane[2]=plane[2]/mag;
  }
  
}; /* plane_normal */

/******************************************************************

  void calc_rotation_angle()

******************************************************************/

void calc_rotation_angle2(int nxp, double *x, double *y, double *angle_dx, double *angle_dy)
{
  int ip1, im1, jp1, jm1, tp1, tm1, i, j, n, ntiles, nx;
  double lon_scale;

  nx = nxp-1;
  ntiles = 6;
  for(n=0; n<ntiles; n++) {
    for(j=0; j<nxp; j++) {
      for(i=0; i<nxp; i++) {
	lon_scale = cos(y[n*nxp*nxp+j*nxp+i]*D2R);
	tp1 = n;
	tm1 = n;
	ip1 = i+1;
	im1 = i-1;
	jp1 = j;
	jm1 = j;

        if(ip1 >= nxp) {  /* find the neighbor tile. */
	  if(n % 2 == 0) { /* tile 1, 3, 5 */
	    tp1 = n+1;
	    ip1 = 0;
	  }
	  else { /* tile 2, 4, 6 */
	    tp1 = n+2;
	    if(tp1 >= ntiles) tp1 -= ntiles;
	    ip1 = nx-j-1;
	    jp1 = 0;
	  }
	}        
        if(im1 < 0) {  /* find the neighbor tile. */
	  if(n % 2 == 0) { /* tile 1, 3, 5 */
	    tm1 = n-2;
	    if(tm1 < 0) tm1 += ntiles;
	    jm1 = nx;
	    im1 = nx-j;
	  }
	  else { /* tile 2, 4, 6 */
	    tm1 = n-1;
	    im1 = nx;
	  }
	}

	angle_dx[n*nxp*nxp+j*nxp+i] = atan2(y[tp1*nxp*nxp+jp1*nxp+ip1]-y[tm1*nxp*nxp+jm1*nxp+im1],
					    (x[tp1*nxp*nxp+jp1*nxp+ip1]-x[tm1*nxp*nxp+jm1*nxp+im1])*lon_scale )*R2D;
	tp1 = n;
	tm1 = n;
	ip1 = i;
	im1 = i;
	jp1 = j+1;
	jm1 = j-1;
        if(jp1 >=nxp) {  /* find the neighbor tile. */
	  if(n % 2 == 0) { /* tile 1, 3, 5 */
	    tp1 = n+2;
	    if(tp1 >= ntiles) tp1 -= ntiles;
	    jp1 = nx-i;
	    ip1 = 0;
	  }
	  else { /* tile 2, 4, 6 */
	    tp1 = n+1;
	    if(tp1 >= ntiles) tp1 -= ntiles;
	    jp1 = 0;
	  }
	}        
        if(jm1 < 0) {  /* find the neighbor tile. */
	  if(n % 2 == 0) { /* tile 1, 3, 5 */
	    tm1 = n-1;
	    if(tm1 < 0) tm1 += ntiles;
	    jm1 = nx;
	  }
	  else { /* tile 2, 4, 6 */
	    tm1 = n-2;
	    if(tm1 < 0) tm1 += ntiles;
	    im1 = nx;
	    jm1 = nx-i;
	  }
	}	

	angle_dy[n*nxp*nxp+j*nxp+i] = atan2(y[tp1*nxp*nxp+jp1*nxp+ip1]-y[tm1*nxp*nxp+jm1*nxp+im1],
					    (x[tp1*nxp*nxp+jp1*nxp+ip1]-x[tm1*nxp*nxp+jm1*nxp+im1])*lon_scale )*R2D;
      }
    }
  }

}; /* calc_rotation_angle2 */


/* This routine calculate center location based on the vertices location */
void cell_center(int ni, int nj, const double *lonc, const double *latc, double *lont, double *latt)
{

  int    nip, njp, i, j, p, p1, p2, p3, p4;
  double *xc, *yc, *zc, *xt, *yt, *zt;
  double dd;
  
  nip = ni+1;
  njp = nj+1;
  xc = (double *)malloc(nip*njp*sizeof(double));
  yc = (double *)malloc(nip*njp*sizeof(double));
  zc = (double *)malloc(nip*njp*sizeof(double));
  xt = (double *)malloc(ni *nj *sizeof(double));
  yt = (double *)malloc(ni *nj *sizeof(double));
  zt = (double *)malloc(ni *nj *sizeof(double));  
  latlon2xyz(nip*njp, lonc, latc, xc, yc, zc); 

  for(j=0; j<nj; j++) for(i=0; i<ni; i++) {
      p =  j*ni+i;
      p1 = j*nip+i;
      p2 = j*nip+i+1;
      p3 = (j+1)*nip+i+1;
      p4 = (j+1)*nip+i;
      xt[p] = xc[p1] + xc[p2] + xc[p3] + xc[p4];
      yt[p] = yc[p1] + yc[p2] + yc[p3] + yc[p4];
      zt[p] = zc[p1] + zc[p2] + zc[p3] + zc[p4];

      dd = sqrt(pow(xt[p],2) + pow(yt[p],2) + pow(zt[p],2) );
      xt[p] /= dd;
      yt[p] /= dd;
      zt[p] /= dd;
    }
  xyz2latlon(ni*nj, xt, yt, zt, lont, latt);
  free(zt);
  free(yt);
  free(xt);
  free(zc);
  free(yc);
  free(xc);
  
}; /* cell_center */


/* This routine calculate east location based on the vertices location */
void cell_east(int ni, int nj, const double *lonc, const double *latc, double *lone, double *late)
{

  int    nip, njp, i, j, p, p1, p2;
  double *xc, *yc, *zc, *xe, *ye, *ze;
  double dd;
  
  nip = ni+1;
  njp = nj+1;
  xc = (double *)malloc(nip*njp*sizeof(double));
  yc = (double *)malloc(nip*njp*sizeof(double));
  zc = (double *)malloc(nip*njp*sizeof(double));
  xe = (double *)malloc(nip*nj *sizeof(double));
  ye = (double *)malloc(nip*nj *sizeof(double));
  ze = (double *)malloc(nip*nj *sizeof(double));  
  latlon2xyz(nip*njp, lonc, latc, xc, yc, zc); 

  for(j=0; j<nj; j++) for(i=0; i<nip; i++) {
      p =  j*nip+i;
      p1 = j*nip+i;
      p2 = (j+1)*nip+i;
      xe[p] = xc[p1] + xc[p2];
      ye[p] = yc[p1] + yc[p2];
      ze[p] = zc[p1] + zc[p2];

      dd = sqrt(pow(xe[p],2) + pow(ye[p],2) + pow(ze[p],2) );
      xe[p] /= dd;
      ye[p] /= dd;
      ze[p] /= dd;
    }
  xyz2latlon(nip*nj, xe, ye, ze, lone, late);
  free(ze);
  free(ye);
  free(xe);
  free(zc);
  free(yc);
  free(xc);
  
}; /* cell_east */


/* This routine calculate center location based on the vertices location */
void cell_north(int ni, int nj, const double *lonc, const double *latc, double *lonn, double *latn)
{

  int    nip, njp, i, j, p, p1, p2;
  double *xc, *yc, *zc, *xn, *yn, *zn;
  double dd;
  
  nip = ni+1;
  njp = nj+1;
  xc = (double *)malloc(nip*njp*sizeof(double));
  yc = (double *)malloc(nip*njp*sizeof(double));
  zc = (double *)malloc(nip*njp*sizeof(double));
  xn = (double *)malloc(ni *njp*sizeof(double));
  yn = (double *)malloc(ni *njp*sizeof(double));
  zn = (double *)malloc(ni *njp*sizeof(double));  
  latlon2xyz(nip*njp, lonc, latc, xc, yc, zc); 

  for(j=0; j<njp; j++) for(i=0; i<ni; i++) {
      p =  j*ni+i;
      p1 = j*nip+i;
      p2 = j*nip+i+1;
      xn[p] = xc[p1] + xc[p2];
      yn[p] = yc[p1] + yc[p2];
      zn[p] = zc[p1] + zc[p2];

      dd = sqrt(pow(xn[p],2) + pow(yn[p],2) + pow(zn[p],2) );
      xn[p] /= dd;
      yn[p] /= dd;
      zn[p] /= dd;
    }
  xyz2latlon(ni*njp, xn, yn, zn, lonn, latn);
  free(zn);
  free(yn);
  free(xn);
  free(zc);
  free(yc);
  free(xc);
  
}; /* cell_north */

/*-------------------------------------------------------------------------------------------
  void spherical_linear_interpolation
  This formula interpolates along the great circle connecting points p1 and p2.
  This formula is taken from http://en.wikipedia.org/wiki/Slerp and is
  attributed to Glenn Davis based on a concept by Ken Shoemake.
  -------------------------------------------------------------------------------------------*/
void spherical_linear_interpolation(double beta, const double *p1, const double *p2, double *pb)
{

  double pm[2];
  double e1[3], e2[3], eb[3];
  double dd, alpha, omega;
 
  if ( fabs(p1[0] - p2[0]) < EPSLN8 && fabs(p1[1] - p2[1]) < EPSLN8 ) {
    printf("WARNING from create_gnomonic_cubic_grid: spherical_linear_interpolation was passed two colocated points.\n");
    pb[0] = p1[0];
    pb[1] = p1[1];
    return ;
  }

  latlon2xyz(1, p1, p1+1, e1, e1+1, e1+2);
  latlon2xyz(1, p2, p2+1, e2, e2+1, e2+2);

  dd = sqrt( e1[0]*e1[0] + e1[1]*e1[1] + e1[2]*e1[2]);

  e1[0] /= dd;
  e1[1] /= dd;
  e1[2] /= dd;

  dd = sqrt( e2[0]*e2[0] + e2[1]*e2[1] + e2[2]*e2[2]);
  e2[0] /= dd;
  e2[1] /= dd;
  e2[2] /= dd;

  alpha = 1. - beta;

  omega = acos( e1[0]*e2[0] + e1[1]*e2[1] + e1[2]*e2[2] );

  if ( fabs(omega) < EPSLN5 ) {
    printf("spherical_linear_interpolation: omega=%g, p1 = %g,%g, p2 = %g,%g\n",
	   omega, p1[0], p1[1], p2[0], p2[1]);
    mpp_error("spherical_linear_interpolation: interpolation not well defined between antipodal points");
  }

  eb[0] = sin( beta*omega )*e2[0] + sin(alpha*omega)*e1[0];
  eb[1] = sin( beta*omega )*e2[1] + sin(alpha*omega)*e1[1];
  eb[2] = sin( beta*omega )*e2[2] + sin(alpha*omega)*e1[2];

  eb[0] /= sin(omega);
  eb[1] /= sin(omega);
  eb[2] /= sin(omega);
  
  xyz2latlon(1, eb, eb+1, eb+2, pb, pb+1);
}



/* void setup_aligned_nest

   /*

   parent_ni    : (input) parent grid size in x-direction.
   parent_nj    : (input) parent grid size in y-direction.
   parent_xc    : (input) parent array in x-direction
   parent_yc    : (input) parent array in y-direction
   halo         : (input) halo size
   refine_ratio : (input) refinement ratio
   istart       : (input) start of nest in x direction
   iend         : (input) end of nest in x direction
   jstart       : (input) start of nest in y direction
   jend         : (input) end of nest in y direction

   xc           : (output) nest array in x-direction
   yc           : (output) nest array in y-direction

*/


void setup_aligned_nest(int parent_ni, int parent_nj, const double *parent_xc, const double *parent_yc,
			int halo, int refine_ratio, int istart, int iend, int jstart, int jend,
			double *xc, double *yc)	
{
  double q1[2], q2[2], t1[2], t2[2], p1[0], p2[0];
  double two_pi;
  int    ni, nj, npi, npj;
  int    parent_npi, i, j, ic, jc, imod, jmod;
  
  int verbose = 1;

  two_pi = 2.*M_PI;
    
  /* Check that the grid does not lie outside its parent */
  if( (jstart - halo) < 1 || (istart - halo) < 1 ||
      (jend + halo) > parent_nj || (iend + halo) > parent_ni )
    mpp_error("create_gnomonic_cubic_grid(setup_aligned_nest): nested grid lies outside its parent");

  if (verbose) {
    printf("[INFO] setup_aligned nest: parent_ni: %d parent_nj: %d refine_ratio: %d parent_xc: %p parent_yc: %p\n", 
	   parent_ni, parent_nj, refine_ratio, parent_xc, parent_yc);

  }

  ni = (iend-istart+1)*refine_ratio;
  nj = (jend-jstart+1)*refine_ratio;
  npi = ni+1;
  npj = nj+1;
  parent_npi = parent_ni+1;
  
  for(j=0; j<npj; j++) {
    jc = jstart - 1 + j/refine_ratio;
    jmod = j%refine_ratio;
    for(i=0; i<npi; i++) {
      ic = istart - 1 + i/refine_ratio;
      imod = i%refine_ratio;

      if(jmod == 0) {
	q1[0] = parent_xc[jc*parent_npi+ic];
	q1[1] = parent_yc[jc*parent_npi+ic];
	q2[0] = parent_xc[jc*parent_npi+ic+1];
	q2[1] = parent_yc[jc*parent_npi+ic+1];
      }
      else {
	t1[0] = parent_xc[jc*parent_npi+ic];
	t1[1] = parent_yc[jc*parent_npi+ic];
	t2[0] = parent_xc[(jc+1)*parent_npi+ic];
	t2[1] = parent_yc[(jc+1)*parent_npi+ic];
	spherical_linear_interpolation( (double)jmod/refine_ratio, t1, t2, q1);
	t1[0] = parent_xc[jc*parent_npi+ic+1];
	t1[1] = parent_yc[jc*parent_npi+ic+1];
	t2[0] = parent_xc[(jc+1)*parent_npi+ic+1];
	t2[1] = parent_yc[(jc+1)*parent_npi+ic+1];	
	spherical_linear_interpolation( (double)jmod/refine_ratio, t1, t2, q2);
      }

      if (imod == 0) {
	xc[j*npi+i] = q1[0];
	yc[j*npi+i] = q1[1];
      }
      else {
	spherical_linear_interpolation( (double)imod/refine_ratio, q1, q2, t1 );
	xc[j*npi+i] = t1[0];
	yc[j*npi+i] = t1[1];
      }

      if( xc[j*npi+i] > two_pi ) xc[j*npi+i] -= two_pi;
      if( xc[j*npi+i] < 0. ) xc[j*npi+i] += two_pi;


      if ((i==0) && (j==0)) {
	printf("setup_aligned_nest xc[0]: %f yc[0]: %f\n", xc[0], yc[0]);
	}
      if ((i==1) && (j==0)) {
	printf("setup_aligned_nest xc[1]: %f yc[1]: %f\n", xc[1], yc[1]);
	}


    }
  }

}

