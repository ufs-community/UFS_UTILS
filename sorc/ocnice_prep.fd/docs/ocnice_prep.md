# ocnice_prep

# Introduction

The cpld_gengrid program and associated script related functions create the files required for Fix and IC files for the coupled model.

This document is part of the <a href="../index.html">UFS_UTILS
documentation</a>.

The cpld_gengrid program is part of the
[UFS_UTILS](https://github.com/ufs-community/UFS_UTILS) project.

## Creating Fix and IC files required for the Coupled Model

For the UFS coupled model application S2S or S2SW, the following fix files are required:

- The CICE6 grid and mask file

- The mesh file for the desired OCN/ICE resolution, which is identical for MOM6 and CICE6.

- The mapped ocean mask on the FV3 tiles

- The ESMF regridding weights required to create the CICE6 IC from CPC (SIS2) reanalysis or to map a 1/4 deg MOM6 or CICE6 tripole restart file
to a lower tripole resolution.

- The latitude,longitude,depth and mask arrays required by WW3 to create a mod_def file.

- The ESMF regridding weights required to remap the CICE6 or MOM6 output from tripole grid to a rectilinear grid (optional).

Since MOM6 creates the model grid at runtime (including adjusting the land mask, if required), the required files for CICE and UFSAtm must be created in a pre-processing step using only the MOM6 supergrid, topography and land mask files as input. This allows the mapped ocean mask (used for creating of the ATM ICs) and the CICE6 grid and mask files to be consistent with the run-time configuration of MOM6.

## Background:

### MOM6 grids

The MOM6 supergrid contains a MOM6 grid at twice the desired resolution. The indexing of the supergrid vs the reduced grid is:


            Super Grid               Reduced Grid


        I-1,J+1     I+1,J+1
           X─────X─────X            I-1,J   i,j
           │     │     │               X─────X
           │     │     │               │     │
           │    i│j    │               │  T  │
           X─────X─────X               │     │
           │     │     │               X─────X
           │     │     │          I-1,J-1    I,J-1
           │     │     │
           X─────X─────X
        I-1,J-1     I+1,J-1


MOM6 uses an Arakawa C grid. Within cpld_gridgen, these are referred to as "stagger" locations, and named as follows:

                 Bu────Cv─────Bu
                 │            │
                 │            │
                 Cu    Ct     Cu
                 │            │
                 │            │
                 Bu────Cv─────Bu


### Rotation angles

For the tripole grid, a rotation angle is defined to translate vectors to/from the grid (i-j) orientation from/to true E-W. The rotation angle on ``Ct`` grid points is calculated at run-time in MOM6 (src/initialization/MOM_shared_initialization.F90). However, CICE6 requires a rotation at the corner (``Bu``) grid points. To find these angles, the rotation angle on ``Ct`` points on the opposite side of the tripole fold are used. In cpld_gridgen, these values are found by "flipping over" and changing the sign of the values on the last row of the MOM6 grid. If ``ipL`` and ``ipR`` are the i-indices of the poles along the last j-row:


                ipL-1     ipL    ipL+1            ipR-1     ipR    ipR+1
                   x-------x-------x     |||        x-------x-------x

then after folding along the tripole seam, ``ipL`` and ``ipR`` must align:


                               ipR+1     ipR    ipR-1
                                  x-------x-------x
                               ipL-1     ipL    ipL+1
                                  x-------x-------x


Using the folded seam, the values of the rotation on ``Ct`` points across the seam are known. The same procedure that CICE uses internally to calculate the ``Ct`` angles from the ``Bu`` angles can be used to instead calculate the ``Bu`` angles knowing the ``Ct`` angles.

### SCRIP format files

For calculating interpolation weights using ESMF, a SCRIP file needs to be provided. A SCIP file contains the both the grid locations of any stagger grid location (e.g. ``Ct``) and the associated grid vertices for that point. As seen from the above diagram, for the ``Ct`` points, those grid vertices are given by the ``Bu`` grid locations.

SCRIP requires that the vertices be ordered counter-clockwise so that the center grid point is always to the left of the vertex. In cpld_gridgen, vertices are defined counter-clockwise from upper right. ``Ct`` vertices are located on the ``Bu`` grid (as shown above), ``Cu`` vertices on the ``Cv`` grid, ``Cv`` vertices on the ``Cu`` grid and ``Bu`` vertices on the ``Ct`` grid. For example, for the ``Ct`` grid, the vertices are:

             Vertex #2             Vertex #1
             Bu(i-1,j)             Bu(i,j)
                         Ct(i,j)
           Bu(i-1,j-1)             Bu(i,j-1)
             Vertex #3             Vertex #4


so that the vertices for the ``Ct`` grid are found as off-sets of the i,j index of the ``Bu`` grid

     iVertCt(4) = (/0, -1, -1,  0/)
     jVertCt(4) = (/0,  0, -1, -1/)

Careful examination of the remaining stagger locations lead to similar definitions for the i,j offsets required to extract the vertices, all of which can be defined in terms of the ``iVertCt`` and ``jVertCt`` values.

Special treatment is require at the bottom of the grid, where the vertices of the ``Ct`` and ``Cu`` grid must be set manually (note, these points are on land.) The top of the grid also requires special treatment because the required vertices are located across the tripole seam. This is accomplished by creating 1-d arrays which hold the ``Ct`` and ``Cu`` grid point locations across the matched seam.


## Generating the grid files

The cpld_gridgen program and associated script related functions perform the following tasks:

1. read the MOM6 supergrid and ocean mask file and optionally creates the required *topo_edits* file if the land mask for MOM6 is to be changed at runtime.
2. create a master grid file containing all stagger locations of the grid fully defined
3. create the CICE6 grid variables and writes the required CICE6 grid file
4. create a SCRIP file for the center stagger (``Ct``) grid points and a second SCRIP file also containing the land mask
5. create the ESMF conservative regridding weights to map the ocean mask to the FV3 tiles and write the mapped mask to 6 tile files
6. create the ESMF regridding weights to map a 1/4 deg ice or ocean restart file to a lower resolution tripole grid
7. optionally call a routine to generate ESMF regridding weights to map the tripole grid to a set of rectilinear grids
8. use the command line command *ESMF_Scrip2Unstruct* to generate the ocean mesh from the SCRIP file containing the land mask (item 4)
9. use an NCO command line command to generate the CICE6 land mask file from the CICE6 grid file

## The generated files

The exact list of files produced by the *cpld_gridgen.sh* script will vary depending on several factors. For example, if the *DO_POSTWGHTS* flag is true, then a SCRIP format file will be produced for each rectilinear destination grid desired and a file containing the regridding weights to map from the center ``Ct`` stagger point to the rectilinear grid will also be written. Because both MOM6 and CICE6 velocity variables are located at ``Cu``,``Cv`` or ``Bu`` locations, additional files will also be created to regrid from the velocity stagger locations to the center ``Ct`` location. If an OCN/ICE grid resolution less than 1/4 degree is chosen, then a file containing regridding weights from the 1/4 degree grid to a lower resolution grid will also be written. Note also that multiple intermediate SCRIP format files may be produced depending on the options chosen.

<br>

* Executing the script for the 1/4 deg OCN/ICE resolution will result in the following files being produced in the output location:


<table>
<caption id="foutmx025">Output files for 1/4 deg</caption>
<tr><th>File name                      <th>Description                              <th>Function
<tr><td row=1>tripole.mx025.nc         <td>master grid file                         <td>Creating all subsequent grid or mapping files
<tr><td row=2>grid_cice_NEMS_mx025.nc  <td>the CICE grid file                       <td>used at runtime by CICE6
<tr><td row=3>kmtu_cice_NEMS_mx025.nc  <td>the CICE mask file                       <td>used at runtime by CICE6
<tr><td row=4>mesh.mx025.nc            <td>the ocean and ice mesh file              <td>used at runtime by CICE6, MOM6, and CMEPS
<tr><td row=5>C384.mx025.tile[1-6].nc  <td>the mapped ocean mask on the ATM tiles   <td>used to create ATM ICs consistent with the <br>                                                                                           fractional grid
</table>

<br>

* If the optional post-weights are generated, the following files will be produced in the output location:


<table>
<caption id="foutpost">Optional post-weights files for 1/4deg</caption>
<tr><th>File name                                                                       <th>Function
<tr><td row=1>tripole.mx025.[Cu][Cv][Bu].to.Ct.bilinear.nc                              <td>the ESMF weights for mapping OCN or ICE <br>                                                                                              output fields from the various stagger <br>                                                                                               locations on the tripole grid to the <br>                                                                                                 center (Ct) stagger location on the <br>
                                                                                            same tripole grid using bilinear mapping
<tr><td row=2>tripole.mx025.Ct.to.rect.[destination resolution].[bilinear][conserve].nc <td>the ESMF weights for mapping variables <br>                                                                                               on the center (Ct) stagger location on <br>                                                                                               the tripole grid to a rectilinear grid <br>                                                                                               with [destination resolution] using <br>                                                                                                  either bilinear or conservative mapping
</table>

<br>

* The following file will be produced in the output location to map 1/4 degree tripole values to a tripole grid of lower resolution.

<table>
<caption id="fouttripole">Output files for down-sampled IC creation at tripole destination resolution</caption>
<tr><th>File name                                                               <th>Function
<tr><td row=1>tripole.mx025.Ct.to.mx[destination resolution].Ct.neareststod.nc  <td>the ESMF weights for mapping 1/4 deg tripole ICs to <br>
                                                                                    a lower tripole destination resolution using nearest <br>                                                                                                  source-to-destination mapping
<tr><td row=1>tripole.mx025.Ct.to.mx[destination resolution].Ct.bilinear.nc  <td>the ESMF weights for mapping 1/4 deg tripole ICs to <br>
                                                                                 a lower tripole destination resolution using <br>                                                                                                          bilinear mapping
<tr><td row=1>tripole.[destination resolution].Ct.to.[Cu][Cv][Bu].bilinear.nc <td> the ESMF weights for mapping downscaled IC values on a <br>
				   					      tripole grid from Ct locations to the native stagger locations
</table>

<br>

* If run-time land mask changes for MOM6 are requested, the following file will be produced in the output location:


<table>
<caption id="fouttopo">Output files for run-time modification of MOM6 land mask</caption>
<tr><th>File name                       <th>Function
<tr><td row=1>ufs.[Default filename].nc <td>Topo-edits required for UFS application. These are appended to the existing default topo <br>                                             edits file and implemented at run time with the parameter flag <br>                                                                       ``ALLOW_LANDMASK_CHANGES=true``. All files produced by the *cpld_gridgen.sh* will be <br>                                                 consistent with this run-time land mask.
</table>
