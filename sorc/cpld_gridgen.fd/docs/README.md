# GridGen
Generate fixed grid files from MOM6 supergrid

 Denise.Worthen@noaa.gov

 This code generates files the fixed grid and land mask file for the CICE
 model on the MOM6 tripole grid. It also generates a fixed grid file which
 contains all the vertice locations on the tripole grid. This fixed grid
 file is used to create the interpolation weights for regridding between
 various combinations of tripole and rectilinear grids.

 Information on MOM6 supergrid can be found at
 https://gist.github.com/adcroft/c1e207024fe1189b43dddc5f1fe7dd6c

 also: https://mom6.readthedocs.io/en/latest/api/generated/modules/mom_grid.html

 and also from:

 MOM_grid_initialize.F90 :
  MOM6 variable geoLonBu <==> CICE variable ulon
  MOM6 variable geoLatBu <==> CICE variable ulat
  MOM6 variable     dxCv <==> CICE variable htn
  MOM6 variable     dyCu <==> CICE variable hte

 MOM6 code snippets follow:

 from MOM_grid_initialize.F90  (tmpZ = x)
```
  do J=G%JsdB,G%JedB ; do I=G%IsdB,G%IedB ; i2 = 2*I ; j2 = 2*J
    G%geoLonBu(I,J) = tmpZ(i2,j2)
```
 so....
```
          ulon(I,J) = x(i2,j2)
```
 from MOM_grid_initialize.F90  (tmpZ = y)
```
  do J=G%JsdB,G%JedB ; do I=G%IsdB,G%IedB ; i2 = 2*I ; j2 = 2*J
    G%geoLatBu(I,J) = tmpZ(i2,j2)
```
 so....
```
          ulat(I,J) = y(i2,j2)
```
 from MOM_grid_initialize.F90  (tmpV = dx)
```
  do J=G%JsdB,G%JedB ; do i=G%isd,G%ied ; i2 = 2*i ; j2 = 2*j
    dxCv(i,J) = tmpV(i2-1,j2) + tmpV(i2,j2)
```
 so....
```
     htn(i,J) =   dx(i2-1,j2) +   dx(i2,j2)
```

 from MOM_grid_initialize.F90  (tmpU = dy)
```
  do J=G%JsdB,G%JedB ; do i=G%isd,G%ied ; i2 = 2*i ; j2 = 2*j
    dyCu(I,j) = tmpU(i2,j2-1) + tmpU(i2,j2)
```
 so....
```
     hte(I,j) =   dy(i2,j2-1) +   dy(i2,j2)
```

 rotation angle on supergrid vertices can be found
 using the formula in MOM_shared_initialization.F90, accounting
 for indexing difference between reduced grid and super grid

```
         SuperGrid                  Reduced grid

  i-1,j+1         i+1,j+1
     X-------X-------X             I-1,J      I,J
     |       |       |                X-------X
     |       |       |                |       |
     |       | i,j   |                |   T   |
     X-------X-------X                |       |
     |       |       |                X-------X
     |       |       |             I-1,J-1   I,J-1
     |       |       |
     X-------X-------X
  i-1,j-1         i+1,j-1

```
 so that in angle formulae
```
         I==>i+1,I-1==>i-1
         J==>j+1,J-1==>j-1
```

 CICE expects angle to be XY -> LatLon so change the sign from MOM6 value
 This has been determined from the HYCOM code: ALL/cice/src/grid2cice.f

            anglet(i,j) =    -pang(i+i0,  j+j0)   !radians
c           pang is from lon-lat to x-y, but anglet is the reverse

 where anglet is the angle variable being written to the CICE grid file
 and pang is HYCOM's own rotation angle.

 Area of the T-grid cell is obtained as in MOM_grid_initialize where
 tmpV = dx on SG and tmpU is dy on SG

```
    dxT(i,j) = tmpV(i2-1,j2-1) + tmpV(i2,j2-1)
    dyT(i,j) = tmpU(i2-1,j2-1) + tmpU(i2-1,j2)
```

 This code utilizes a "seam flip" to obtain the required values across
 the tripole seam. If ipL an ipR are the i-indices of the pole along the
 last j-row of the reduced grid, then:

```
 ipL-1     ipL    ipL+1       ipR-1     ipR    ipR+1
    x-------x-------x     |||    x-------x-------x
```

 Fold over; ipL must align with ipR
```

  ipR+1     ipR    ipR-1
     x-------x-------x
  ipL-1     ipL    ipL+1
     x-------x-------x
```


 SCRIP requires that the vertices be ordered counter-clockwise so that
 the center grid point is always to the left of the vertex. Here,
 Vertices are defined counter-clockwise from upper right. Ct-grid vertices
 are located on the Bu grid; Cu vertices on the Cv grid, Cv vertices on the Cu
 grid and Bu vertices on the Ct grid. For example, for the Ct-grid, the vertices
 are:
```
             Vertex #2             Vertex #1
             Bu(i-1,j)             Bu(i,j)
                         Ct(i,j)
           Bu(i-1,j-1)             Bu(i,j-1)
             Vertex #3             Vertex #4
```

 so that the vertices of any Ct(i,j) are found as off-sets of the i,j index on the
 Bu grid

```
     iVertCt(4) = (/0, -1, -1, 0/)
     jVertCt(4) = (/0, 0, -1, -1/)
```

 Careful examination of the Cu,Cv and Bu grids lead to similar definitions for the
 i,j offsets required to extract the other grid stragger vertices locations, all of
 which can be defined in terms of the iVertCt and jVertCt values

 Special treatment is require at the bottom of the grid, where the verticies of the
 Ctand Cu grid must be set manually (note, these points are on land.) The top of
 the grid also requires special treatment because the required verticies are located
 across the tripole seam. This is accomplished by creating 1-d arrays which hold
 the Ct and Cu grid point locations across the matched seam.
