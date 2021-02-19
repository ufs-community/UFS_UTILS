/* The following is a test program to test subroutines in create_xgrid.c */
#include "create_xgrid.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define D2R (M_PI/180)
#define R2D (180/M_PI)
#define MAXPOINT 1000

int main(int argc, char* argv[])
{

    double lon1_in[MAXPOINT], lat1_in[MAXPOINT];
    double lon2_in[MAXPOINT], lat2_in[MAXPOINT];
    double x1_in[MAXPOINT], y1_in[MAXPOINT], z1_in[MAXPOINT];
    double x2_in[MAXPOINT], y2_in[MAXPOINT], z2_in[MAXPOINT];
    double lon_out[20], lat_out[20];
    double x_out[20], y_out[20], z_out[20];
    int    n1_in, n2_in, n_out, i, j;
    int    nlon1=0, nlat1=0, nlon2=0, nlat2=0;
    int    n;
    int    ntest = 11;

    printf("Testing create_xgrid.\n");

    for(n=11; n<=ntest; n++) {

        switch (n) {
        case 1:
            /****************************************************************

       test clip_2dx2d_great_cirle case 1:
       box 1: (20,10), (20,12), (22,12), (22,10)
       box 2: (21,11), (21,14), (24,14), (24,11)
       out  : (21, 12.0018), (22, 12), (22, 11.0033), (21, 11)

            ****************************************************************/
            n1_in = 4; n2_in = 4;
            /* first a simple lat-lon grid box to clip another lat-lon grid box */
            lon1_in[0] = 20; lat1_in[0] = 10;
            lon1_in[1] = 20; lat1_in[1] = 12;
            lon1_in[2] = 22; lat1_in[2] = 12;
            lon1_in[3] = 22; lat1_in[3] = 10;
            lon2_in[0] = 21; lat2_in[0] = 11;
            lon2_in[1] = 21; lat2_in[1] = 14;
            lon2_in[2] = 24; lat2_in[2] = 14;
            lon2_in[3] = 24; lat2_in[3] = 11;
            break;

        case 2:
            /****************************************************************

        test clip_2dx2d_great_cirle case 2: two identical box
        box 1: (20,10), (20,12), (22,12), (22,10)
        box 2: (20,10), (20,12), (22,12), (22,10)
        out  : (20,10), (20,12), (22,12), (22,10)

            ****************************************************************/
            lon1_in[0] = 20; lat1_in[0] = 10;
            lon1_in[1] = 20; lat1_in[1] = 12;
            lon1_in[2] = 22; lat1_in[2] = 12;
            lon1_in[3] = 22; lat1_in[3] = 10;

            for(i=0; i<n2_in; i++) {
                lon2_in[i] = lon1_in[i];
                lat2_in[i] = lat1_in[i];
            }
            break;

        case 3:
            /****************************************************************

       test clip_2dx2d_great_cirle case 3: one cubic sphere grid close to the pole with lat-lon grid.
       box 1: (251.7, 88.98), (148.3, 88.98), (57.81, 88.72), (342.2, 88.72)
       box 2: (150, 88), (150, 90), (152.5, 90), (152.5, 88)
       out  : (152.5, 89.0642), (150, 89.0165), (0, 90)

            ****************************************************************/
            n1_in = 4; n2_in = 4;
            /* first a simple lat-lon grid box to clip another lat-lon grid box */
            lon1_in[0] = 251.7; lat1_in[0] = 88.98;
            lon1_in[1] = 148.3; lat1_in[1] = 88.98;
            lon1_in[2] = 57.81; lat1_in[2] = 88.72;
            lon1_in[3] = 342.2; lat1_in[3] = 88.72;

            lon2_in[0] = 150; lat2_in[0] = 88;
            lon2_in[1] = 150; lat2_in[1] = 90;
            lon2_in[2] = 152.5; lat2_in[2] = 90;
            lon2_in[3] = 152.5; lat2_in[3] = 88;
            /*
              for(i=0; i<4; i++) {
              lon2_in[i] = lon1_in[i];
              lat2_in[i] = lat1_in[i];
              }
            */
            break;

        case 4:
            /****************************************************************

       test clip_2dx2d_great_cirle case 4: One box contains the pole
       box 1: (-160, 88.5354), (152.011, 87.8123) , (102.985, 88.4008), (20, 89.8047)
       box 2: (145,88), (145,90), (150,90), (150,88)
       out  : (145.916, 88.0011), (145, 88.0249), (0, 90), (150, 88)

            ****************************************************************/
            n1_in = 4; n2_in = 4;
            /* first a simple lat-lon grid box to clip another lat-lon grid box */

            lon1_in[0] = -160;  lat1_in[0] = 88.5354;
            lon1_in[1] = 152.011; lat1_in[1] = 87.8123;
            lon1_in[2] = 102.985; lat1_in[2] = 88.4008;
            lon1_in[3] = 20; lat1_in[3] = 89.8047;

            lon2_in[0] = 145; lat2_in[0] = 88;
            lon2_in[1] = 145; lat2_in[1] = 90;
            lon2_in[2] = 150; lat2_in[2] = 90;
            lon2_in[3] = 150; lat2_in[3] = 88;
            break;

        case 5:
            /****************************************************************

       test clip_2dx2d_great_cirle case 5: One tripolar grid around the pole with lat-lon grid.
       box 1: (-202.6, 87.95), (-280, 89.56), (-100, 90), (-190, 88)
       box 2: (21,11), (21,14), (24,14), (24,11)
       out  : (150, 88.7006), (145,  88.9507), (0, 90)

            ****************************************************************/
            n1_in = 4; n2_in = 4;
            /* first a simple lat-lon grid box to clip another lat-lon grid box */

            lon1_in[0] = -202.6;  lat1_in[0] = 87.95;
            lon1_in[1] = -280.;   lat1_in[1] = 89.56;
            lon1_in[2] = -100.0; lat1_in[2] = 90;
            lon1_in[3] = -190.; lat1_in[3] = 88;

            lon2_in[0] = 145; lat2_in[0] = 88;
            lon2_in[1] = 145; lat2_in[1] = 90;
            lon2_in[2] = 150; lat2_in[2] = 90;
            lon2_in[3] = 150; lat2_in[3] = 88;
            break;

        case 6:
            /****************************************************************

       test clip_2dx2d_great_cirle case 6: One cubic sphere grid arounc the pole with one tripolar grid box
                                       around the pole.
       box 1: (-160, 88.5354), (152.011, 87.8123) , (102.985, 88.4008), (20, 89.8047)
       box 2: (-202.6, 87.95), (-280, 89.56), (-100, 90), (-190, 88)
       out  : (170, 88.309), (157.082, 88.0005), (83.714, 89.559), (80, 89.6094), (0, 90), (200, 88.5354)


            ****************************************************************/
            n1_in = 4; n2_in = 4;
            /* first a simple lat-lon grid box to clip another lat-lon grid box */

            lon1_in[0] = -160;  lat1_in[0] = 88.5354;
            lon1_in[1] = 152.011; lat1_in[1] = 87.8123;
            lon1_in[2] = 102.985; lat1_in[2] = 88.4008;
            lon1_in[3] = 20; lat1_in[3] = 89.8047;

            lon2_in[0] = -202.6;  lat2_in[0] = 87.95;
            lon2_in[1] = -280.;   lat2_in[1] = 89.56;
            lon2_in[2] = -100.0;  lat2_in[2] = 90;
            lon2_in[3] = -190.;   lat2_in[3] = 88;
            break;

        case 7:
            /****************************************************************

       test clip_2dx2d_great_cirle case 7: One small grid box inside a big grid box.
       box 1: (20,10), (20,12), (22,12), (22,10)
       box 2: (18,8), (18,14), (24,14), (24,8)
       out  : (20,10), (20,12), (22,12), (22,10)

            ****************************************************************/
            n1_in = 4; n2_in = 4;
            /* first a simple lat-lon grid box to clip another lat-lon grid box */
            lon1_in[0] = 20; lat1_in[0] = 10;
            lon1_in[1] = 20; lat1_in[1] = 12;
            lon1_in[2] = 22; lat1_in[2] = 12;
            lon1_in[3] = 22; lat1_in[3] = 10;
            lon2_in[0] = 18; lat2_in[0] = 8;
            lon2_in[1] = 18; lat2_in[1] = 14;
            lon2_in[2] = 24; lat2_in[2] = 14;
            lon2_in[3] = 24; lat2_in[3] = 8;
            break;

        case 8:
            /****************************************************************

       test clip_2dx2d_great_cirle case 8: Cubic sphere grid at tile = 1, point (i=25,j=1)
          with N45 at (i=141,j=23)
       box 1:
       box 2:
       out  : None

            ****************************************************************/
            n1_in = 4; n2_in = 4;
            /* first a simple lat-lo
               n grid box to clip another lat-lon grid box */
            lon1_in[0] = 350.0; lat1_in[0] = -45;
            lon1_in[1] = 350.0; lat1_in[1] = -43.43;
            lon1_in[2] = 352.1; lat1_in[2] = -43.41;
            lon1_in[3] = 352.1; lat1_in[3] = -44.98;
            lon2_in[0] = 350.0;   lat2_in[0] = -46;
            lon2_in[1] = 350.0;   lat2_in[1] = -44;
            lon2_in[2] = 352.5; lat2_in[2] = -44;
            lon2_in[3] = 352.5; lat2_in[3] = -46;
            break;

        case 9:
            /****************************************************************

       test clip_2dx2d_great_cirle case 9: Cubic sphere grid at tile = 1, point (i=1,j=1)
          with N45 at (i=51,j=61)
       box 1:
       box 2:
       out  : None

            ****************************************************************/
            n1_in = 4; n2_in = 4;

            lon1_in[0] = 305.0; lat1_in[0] = -35.26;
            lon1_in[1] = 305.0; lat1_in[1] = -33.80;
            lon1_in[2] = 306.6; lat1_in[2] = -34.51;
            lon1_in[3] = 306.6; lat1_in[3] = -35.99;
            lon2_in[0] = 125;   lat2_in[0] = 32;
            lon2_in[1] = 125;   lat2_in[1] = 34;
            lon2_in[2] = 127.5; lat2_in[2] = 34;
            lon2_in[3] = 127.5; lat2_in[3] = 32;
            break;

        case 10:
            /****************************************************************

       test clip_2dx2d_great_cirle case 10: Cubic sphere grid at tile = 3, point (i=24,j=1)
          with N45 at (i=51,j=46)
       box 1:
       box 2:
       out  : None

            ****************************************************************/
            n1_in = 4; n2_in = 4;

            lon1_in[0] = 125.0; lat1_in[0] = 1.46935;
            lon1_in[1] = 126.573; lat1_in[1] = 1.5091;
            lon1_in[2] = 126.573; lat1_in[2] = 0;
            lon1_in[3] = 125.0; lat1_in[3] = 0;
            lon2_in[0] = 125;   lat2_in[0] = 0;
            lon2_in[1] = 125;   lat2_in[1] = 2;
            lon2_in[2] = 127.5; lat2_in[2] = 2;
            lon2_in[3] = 127.5; lat2_in[3] = 0;
            break;

        case 11:
            /****************************************************************

       test clip_2dx2d_great_cirle case 10: Cubic sphere grid at tile = 3, point (i=24,j=1)
          with N45 at (i=51,j=46)
       box 1:
       box 2:
       out  :

            ****************************************************************/
            nlon1 = 1;
            nlat1 = 1;
            nlon2 = 1;
            nlat2 = 1;
            n1_in = (nlon1+1)*(nlat1+1);
            n2_in = (nlon2+1)*(nlat2+1);

            lon1_in[0] = 350.0; lat1_in[0] = 90.00;
            lon1_in[1] = 170.0; lat1_in[1] = 87.92;
            lon1_in[2] = 260.0; lat1_in[2] = 87.92;
            lon1_in[3] = 215.0;  lat1_in[3] = 87.06;

/*       lon1_in[0] = 35.0; lat1_in[0] = 87.06; */
/*       lon1_in[1] = 80.0; lat1_in[1] = 87.92; */
/*       lon1_in[2] = 125.0; lat1_in[2] = 87.06; */
/*       lon1_in[3] = 350.0; lat1_in[3] = 87.92; */
/*       lon1_in[4] = 350.0; lat1_in[4] = 90.00; */
/*       lon1_in[5] = 170.0; lat1_in[5] = 87.92; */
/*       lon1_in[6] = 305.0; lat1_in[6] = 87.06; */
/*       lon1_in[7] = 260.0; lat1_in[7] = 87.92; */
/*       lon1_in[8] = 215.0;  lat1_in[8] = 87.06; */

            lon2_in[0] = 167.5; lat2_in[0] = 88;
            lon2_in[1] = 170;   lat2_in[1] = 88;
            lon2_in[2] = 167.5; lat2_in[2] = 90;
            lon2_in[3] = 170;   lat2_in[3] = 90;

/*       nlon1 = 3; */
/*       nlat1 = 2; */
/*       nlon2 = 1; */
/*       nlat2 = 1; */
/*       n1_in = (nlon1+1)*(nlat1+1); */
/*       n2_in = (nlon2+1)*(nlat2+1); */

/*       lon1_in[0] = 35.00;     lat1_in[0] = -59.90; */
/*       lon1_in[1] = 37.64;     lat1_in[1] = -58.69; */
/*       lon1_in[2] = 40.07;     lat1_in[2] = -57.44; */
/*       lon1_in[3] = 42.32;     lat1_in[3] = -56.15; */
/*       lon1_in[4] = 32.36;     lat1_in[4] = -58.69; */
/*       lon1_in[5] = 35.00;     lat1_in[5] = -57.56; */
/*       lon1_in[6] = 37.45;     lat1_in[6] = -56.39; */
/*       lon1_in[7] = 39.74;     lat1_in[7] = -55.18; */
/*       lon1_in[8] = 29.93;     lat1_in[8] = -57.44; */
/*       lon1_in[9] = 32.55;     lat1_in[9] = -56.39; */
/*       lon1_in[10] = 35.00;     lat1_in[10] = -55.29; */
/*       lon1_in[11] = 37.30;     lat1_in[11] = -54.16; */
/*       lon2_in[0] = 35;   lat2_in[0] = -58; */
/*       lon2_in[1] = 37.5; lat2_in[1] = -58; */
/*       lon2_in[2] = 35;   lat2_in[2] = -56; */
/*       lon2_in[3] = 37.5; lat2_in[3] = -56; */

/*       nlon1 = 1; */
/*       nlat1 = 1; */
/*       nlon2 = 1; */
/*       nlat2 = 1; */
/*       n1_in = (nlon1+1)*(nlat1+1); */
/*       n2_in = (nlon2+1)*(nlat2+1); */

/*       lon1_in[0] = 305;     lat1_in[0] = -35.26; */
/*       lon1_in[1] = 306;     lat1_in[1] = -35.99; */
/*       lon1_in[2] = 305;     lat1_in[2] = -33.80; */
/*       lon1_in[3] = 306;     lat1_in[3] = -34.51; */
/*       lon2_in[0] = 305;   lat2_in[0] = -34; */
/*       lon2_in[1] = 307.5; lat2_in[1] = -34; */
/*       lon2_in[2] = 305;   lat2_in[2] = -32; */
/*       lon2_in[3] = 307.5; lat2_in[3] = -32; */

            nlon1 = 2;
            nlat1 = 2;
            nlon2 = 1;
            nlat2 = 1;
            n1_in = (nlon1+1)*(nlat1+1);
            n2_in = (nlon2+1)*(nlat2+1);

            lon1_in[0] = 111.3; lat1_in[0] = 1.591;
            lon1_in[1] = 109.7; lat1_in[1] = 2.926;
            lon1_in[2] = 108.2; lat1_in[2] = 4.256;
            lon1_in[3] = 110.0; lat1_in[3] = 0.000;
            lon1_in[4] = 108.4; lat1_in[4] = 1.335;
            lon1_in[5] = 106.8; lat1_in[5] = 2.668;
            lon1_in[6] = 108.7; lat1_in[6] = -1.591;
            lon1_in[7] = 107.1; lat1_in[7] = -0.256;
            lon1_in[8] = 105.5;  lat1_in[8] = 1.078;

            lon2_in[0] = 107.5; lat2_in[0] = 0;
            lon2_in[1] = 110;   lat2_in[1] = 0;
            lon2_in[2] = 107.5; lat2_in[2] = 2;
            lon2_in[3] = 110;   lat2_in[3] = 2;

            break;

        case 12:
            /****************************************************************

       test : create_xgrid_great_circle
       box 1: (20,10), (20,12), (22,12), (22,10)
       box 2: (21,11), (21,14), (24,14), (24,11)
       out  : (21, 12.0018), (22, 12), (22, 11.0033), (21, 11)

            ****************************************************************/
            nlon1 = 2;
            nlat1 = 2;
            nlon2 = 3;
            nlat2 = 3;
            n1_in = (nlon1+1)*(nlat1+1);
            n2_in = (nlon2+1)*(nlat2+1);

            /* first a simple lat-lon grid box to clip another lat-lon grid box */
            for(j=0; j<=nlat1; j++) for(i=0; i<=nlon1; i++){
                    lon1_in[j*(nlon1+1)+i] = 20.0 + (i-1)*2.0;
                    lat1_in[j*(nlon1+1)+i] = 10.0 + (j-1)*2.0;
                }
            for(j=0; j<=nlat2; j++) for(i=0; i<=nlon2; i++){
                    lon2_in[j*(nlon2+1)+i] = 19.0 + (i-1)*2.0;
                    lat2_in[j*(nlon2+1)+i] = 9.0 + (j-1)*2.0;
                }

            break;

        case 13:

            nlon1 = 1;
            nlat1 = 1;
            nlon2 = 1;
            nlat2 = 1;
            n1_in = (nlon1+1)*(nlat1+1);
            n2_in = (nlon2+1)*(nlat2+1);

/*       lon1_in[0] = ; lat1_in[0] = ; */
/*       lon1_in[1] = ; lat1_in[1] = ; */
/*       lon1_in[2] = ; lat1_in[2] = ; */
/*       lon1_in[3] = ; lat1_in[3] = ; */
/*       lon2_in[0] = ; lat2_in[0] = ; */
/*       lon2_in[1] = ; lat2_in[1] = ; */
/*       lon2_in[2] = ; lat2_in[2] = ; */
/*       lon2_in[3] = ; lat2_in[3] = ;     */

/*       lon1_in[0] = 1.35536; lat1_in[0] = 1.16251; */
/*       lon1_in[1] = 1.36805; lat1_in[1] = 1.15369; */
/*       lon1_in[2] = 1.37843; lat1_in[2] = 1.16729; */
/*       lon1_in[3] = 1.39048; lat1_in[3] = 1.15826; */
/*       lon2_in[0] = 1.34611; lat2_in[0] = 1.16372; */
/*       lon2_in[1] = 1.35616; lat2_in[1] = 1.15802;    */
/*       lon2_in[2] = 1.35143; lat2_in[2] = 1.16509; */
/*       lon2_in[3] = 1.36042; lat2_in[3] = 1.15913; */

/*       lon1_in[0] = 12.508065121288551; lat1_in[0] = -87.445883646793547; */
/*       lon1_in[1] = 325.425637772; lat1_in[1] = -86.481216821859505; */
/*       lon1_in[2] = 97.5; lat1_in[2] = -89.802136057677174; */
/*       lon1_in[3] = 277.5; lat1_in[3] = -87.615232005344637; */

/*       for(j=0; j<=nlat2; j++) for(i=0; i<=nlon2; i++) { */
/*      lon2_in[j*(nlon2+1)+i] = -280.0 + i*1.0; */
/*      lat2_in[j*(nlon2+1)+i] = -90.0 + j*8.0; */
/*       } */
            lon1_in[0] = 120.369397984526174; lat1_in[0] = 16.751543427495864;
            lon1_in[1] = 119.999999999999986; lat1_in[1] = 16.751871929590038;
            lon1_in[2] = 120.369397846883501; lat1_in[2] = 16.397797979598028;
            lon1_in[3] = 119.999999999999986; lat1_in[3] = 16.398120477217255;
            lon2_in[0] = 120.369415056522087; lat2_in[0] = 16.752176828509153;
            lon2_in[1] = 119.999999999999986; lat2_in[1] = 16.752505523196167;
            lon2_in[2] = 120.369415056522087; lat2_in[2] = 16.397797949548146;
            lon2_in[3] = 119.999999999999986; lat2_in[3] = 16.398120477217255;

            break;


        default:
            error_handler("test_create_xgrid: incorrect case number");
        }

        /* convert to radian */

        for(i=0; i<n1_in; i++) {
            lon1_in[i] *= D2R; lat1_in[i] *=D2R;
        }
        for(i=0; i<n2_in; i++) {
            lon2_in[i] *= D2R; lat2_in[i] *=D2R;
        }


        printf("\n*********************************************************\n");
        printf("\n               Case %d                                    \n", n);
        printf("\n*********************************************************\n");


        if( n > 10 ) {
            int nxgrid;
            int *i1, *j1, *i2, *j2;
            double *xarea, *xclon, *xclat, *mask1;

            mask1 = (double *)malloc(nlon1*nlat1*sizeof(double));
            i1    = (int    *)malloc(MAXXGRID*sizeof(int));
            j1    = (int    *)malloc(MAXXGRID*sizeof(int));
            i2    = (int    *)malloc(MAXXGRID*sizeof(int));
            j2    = (int    *)malloc(MAXXGRID*sizeof(int));
            xarea = (double *)malloc(MAXXGRID*sizeof(double));
            xclon = (double *)malloc(MAXXGRID*sizeof(double));
            xclat = (double *)malloc(MAXXGRID*sizeof(double));

            for(i=0; i<nlon1*nlat1; i++) mask1[i] = 1.0;

            nxgrid = create_xgrid_great_circle(&nlon1, &nlat1, &nlon2, &nlat2, lon1_in, lat1_in,
                                               lon2_in, lat2_in, mask1, i1, j1, i2, j2,
                                               xarea, xclon, xclat);
            printf("\n*********************************************************\n");
            printf("\n     First input grid box longitude, latitude   \n \n");
            for(i=0; i<n1_in; i++) printf(" %g,  %g \n", lon1_in[i]*R2D, lat1_in[i]*R2D);

            printf("\n     Second input grid box longitude, latitude \n \n");
            for(i=0; i<n2_in; i++) printf(" %g,  %g \n", lon2_in[i]*R2D, lat2_in[i]*R2D);

            printf("\n  Number of exchange grid is %d\n", nxgrid);
            for(i=0; i<nxgrid; i++) {
                printf("(i1,j1)=(%d,%d), (i2,j2)=(%d, %d), xgrid_area=%g, xgrid_clon=%g, xgrid_clat=%g\n",
                       i1[i], j1[i], i2[i], j2[i], xarea[i], xclon[i], xclat[i]);
            }

            /* comparing the area sum of exchange grid and grid1 area */
            {
                double *x1, *y1, *z1, *area1;
                double area_sum;
                int    i;
                area_sum = 0.0;

                for(i=0; i<nxgrid; i++) {
                    area_sum+= xarea[i];
                }

                area1 = (double *)malloc((nlon1)*(nlat1)*sizeof(double));
                get_grid_great_circle_area_(&nlon1, &nlat1, lon1_in, lat1_in, area1);

                printf("xgrid area sum is %g, grid 1 area is %g\n", area_sum, area1[0]);
            }

            printf("\n");
            free(i1);
            free(i2);
            free(j1);
            free(j2);
            free(xarea);
            free(xclon);
            free(xclat);
            free(mask1);
        }
        else {
            latlon2xyz(n1_in, lon1_in, lat1_in, x1_in, y1_in, z1_in);
            latlon2xyz(n2_in, lon2_in, lat2_in, x2_in, y2_in, z2_in);

            n_out = clip_2dx2d_great_circle(x1_in, y1_in, z1_in, 4, x2_in, y2_in, z2_in, n2_in,
                                            x_out, y_out,  z_out);
            xyz2latlon(n_out, x_out, y_out, z_out, lon_out, lat_out);

            printf("\n*********************************************************\n");
            printf("\n     First input grid box longitude, latitude   \n \n");
            for(i=0; i<n1_in; i++) printf(" %g,  %g \n", lon1_in[i]*R2D, lat1_in[i]*R2D);

            printf("\n     Second input grid box longitude, latitude \n \n");
            for(i=0; i<n2_in; i++) printf(" %g,  %g \n", lon2_in[i]*R2D, lat2_in[i]*R2D);

            printf("\n     output clip grid box longitude, latitude for case 1 \n \n");
            for(i=0; i<n_out; i++) printf(" %g,  %g \n", lon_out[i]*R2D, lat_out[i]*R2D);
            printf("\n");
        }
    }

    printf("SUCCESS!\n");
}
