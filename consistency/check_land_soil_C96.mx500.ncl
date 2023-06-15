begin

;fixpath = "/scratch1/NCEPDEV/stmp2/Michael.Barlage/FV3_RT/rt_147523/control_p8/INPUT/"
;sfcpath = "/scratch1/NCEPDEV/stmp2/Michael.Barlage/FV3_RT/rt_147523/control_p8/"
;fixpath = "/scratch1/NCEPDEV/nems/emc.nemspara/RT/NEMSfv3gfs/input-data-20221101/FV3_fix_tiled/C9"
;sfcpath = "/scratch1/NCEPDEV/stmp2/Michael.Barlage/FV3_RT/rt_147523/control_p8/INPUT/"
;sfcpath = "/scratch1/NCEPDEV/nems/emc.nemspara/RT/NEMSfv3gfs/input-data-20221101/FV3_input_frac/C96_L127.mx100_frac/"

fixpath="/scratch2/NCEPDEV/stmp1/Sanath.Kumar/step3_test/fix_sfc/C96.mx100/"
oropath="/scratch2/NCEPDEV/stmp1/Sanath.Kumar/step3_test/C96.mx100/"
#sfcpath="/scratch2/NCEPDEV/stmp1/Sanath.Kumar/data/OUTPUT/C96.mx500/merged/gfs.20210323/12/atmos/INPUT/"


lake_present = False

 num_oro_land_frac_gt0 = 0
 num_oro_lake_frac = 0
 num_oro_land_frac_100 = 0
 num_oro_land_frac_ge50 = 0
 num_oro_slmsk = 0
 num_oro_land_frac_gt0_and_no_fix_veg = 0
 num_oro_land_frac_gt0_and_no_fix_soil = 0
 num_oro_land_frac_gt0_and_no_sfc_veg = 0
 num_oro_land_frac_gt0_and_no_sfc_soil = 0
 num_fix_valid_veg = 0
 num_fix_valid_soil = 0
 num_fix_veg_17 = 0
 num_sfc_slmsk = 0
 num_sfc_valid_veg = 0
 num_sfc_valid_soil = 0
 num_sfc_veg_0 = 0
 num_sfc_soil_0 = 0
 num_sfc_veg_17 = 0
 num_sfc_soil_14 = 0
 fix_veg_sfc_veg_check = 0
 fix_soil_sfc_soil_check = 0
 fix_veg_sfc_soil_check = 0
 fix_soil_sfc_veg_cheeck = 0
 oro_slmsk_sfc_veg_check = 0
 oro_slmsk_sfc_soil_check = 0
 
do itile = 1,6

print("Starting tile: "+itile)

 fix_oro_file       = addfile(oropath+"oro_C96.mx100.tile"+itile+".nc","r")
 fix_vegtype_file   = addfile(fixpath+"C96.vegetation_type.tile"+itile+".nc","r")
 fix_soil_file      = addfile(fixpath+"C96.soil_type.tile"+itile+".nc","r")
 #sfc_tile_files     = addfile(sfcpath+"merged_sfc_data.tile"+itile+".nc","r")
 
 fix_oro_slmsk      = fix_oro_file->slmsk 
 fix_oro_land_frac  = fix_oro_file->land_frac
if(lake_present) then
 fix_oro_lake_frac  = fix_oro_file->lake_frac
end if
 fix_veg_type       = fix_vegtype_file->vegetation_type(0,:,:)
 fix_soil_type      = fix_soil_file->soil_type(0,:,:)
 
 #sfc_slmsk          = sfc_tile_files->slmsk(0,:,:)
 #sfc_veg_type       = sfc_tile_files->vtype(0,:,:)
 #sfc_soil_type      = sfc_tile_files->stype(0,:,:)

 fix_valid_veg  = fix_veg_type
 fix_valid_veg  = where(fix_veg_type.gt.0,1,0)
 fix_valid_soil = fix_veg_type
 fix_valid_soil = where((fix_soil_type.gt.0.and.fix_soil_type.le.13).or.fix_soil_type.ge.15,1,0)
 
 sfc_valid_veg  = sfc_veg_type
 sfc_valid_veg  = where(sfc_veg_type.gt.0,1,0)
 sfc_valid_soil = sfc_veg_type
 sfc_valid_soil = where((sfc_soil_type.gt.0.and.sfc_soil_type.le.13).or.sfc_soil_type.ge.15,1,0)
 
 num_oro_slmsk          = num_oro_slmsk + num(fix_oro_slmsk.eq.1) 
 num_oro_land_frac_gt0  = num_oro_land_frac_gt0 + num(fix_oro_land_frac.le.1.0.and.fix_oro_land_frac.gt.0) 
 num_oro_land_frac_100  = num_oro_land_frac_100 + num(fix_oro_land_frac.eq.1.0) 
 num_oro_land_frac_ge50 = num_oro_land_frac_ge50 + num(fix_oro_land_frac.le.1.0.and.fix_oro_land_frac.ge.0.5) 
if(lake_present) then
 num_oro_lake_frac      = num_oro_lake_frac + num(fix_oro_lake_frac.le.1.0.and.fix_oro_lake_frac.gt.0) 
end if
 num_fix_valid_veg      = num_fix_valid_veg  + num(.not.ismissing(fix_veg_type)) 
 num_fix_valid_soil     = num_fix_valid_soil + num(.not.ismissing(fix_soil_type)) 
 num_oro_land_frac_gt0_and_no_fix_veg = num_oro_land_frac_gt0_and_no_fix_veg + num(fix_oro_land_frac.le.1.0.and.fix_oro_land_frac.gt.0.and.ismissing(fix_veg_type)) 
 num_oro_land_frac_gt0_and_no_fix_soil = num_oro_land_frac_gt0_and_no_fix_soil + num(fix_oro_land_frac.le.1.0.and.fix_oro_land_frac.gt.0.and.ismissing(fix_soil_type)) 
 num_oro_land_frac_gt0_and_no_sfc_veg = num_oro_land_frac_gt0_and_no_sfc_veg + num(fix_oro_land_frac.le.1.0.and.fix_oro_land_frac.gt.0.and.sfc_valid_veg.eq.0) 
 num_oro_land_frac_gt0_and_no_sfc_soil = num_oro_land_frac_gt0_and_no_sfc_soil + num(fix_oro_land_frac.le.1.0.and.fix_oro_land_frac.gt.0.and.sfc_valid_soil.eq.0) 
 num_fix_veg_17         = num_fix_veg_17  + num(fix_veg_type.eq.17) 

 num_sfc_slmsk = num_sfc_slmsk + num(sfc_slmsk.eq.1)
 num_sfc_valid_veg = num_sfc_valid_veg + num(sfc_valid_veg.eq.1)
 num_sfc_valid_soil = num_sfc_valid_soil + num(sfc_valid_soil.eq.1)
 num_sfc_veg_0 = num_sfc_veg_0 + num(sfc_veg_type.eq.0)
 num_sfc_soil_0 = num_sfc_soil_0 + num(sfc_soil_type.eq.0)
 num_sfc_veg_17 = num_sfc_veg_17 + num(sfc_veg_type.eq.17)
 num_sfc_soil_14 = num_sfc_soil_14 + num(sfc_soil_type.eq.14)

 fix_veg_sfc_veg_check = fix_veg_sfc_veg_check + num(fix_valid_veg .ne. sfc_valid_veg)
 fix_soil_sfc_soil_check = fix_soil_sfc_soil_check + num(fix_valid_soil .ne. sfc_valid_soil)
 fix_veg_sfc_soil_check = fix_veg_sfc_soil_check + num(fix_valid_veg .ne. sfc_valid_soil)
 fix_soil_sfc_veg_cheeck = fix_soil_sfc_veg_cheeck + num(fix_valid_soil .ne. sfc_valid_veg)
 oro_slmsk_sfc_veg_check = oro_slmsk_sfc_veg_check + num(fix_oro_slmsk .ne. sfc_valid_veg)
 oro_slmsk_sfc_soil_check = oro_slmsk_sfc_soil_check + num(fix_oro_slmsk .ne. sfc_valid_soil)

 print("number of cumulative fix_oro_slmsk = 1        : " + num_oro_slmsk)
 print("number of cumulative fix_oro_land_frac (0,1]  : " + num_oro_land_frac_gt0)
 print("number of cumulative fix_oro_land_frac = 1    : " + num_oro_land_frac_100)
 print("number of cumulative fix_oro_land_frac [0.5,1]: " + num_oro_land_frac_ge50)
if(lake_present) then
 print("number of cumulative fix_oro_lake_frac (0,1]:   " + num_oro_lake_frac)
end if
 print("number of cumulative valid fix veg type                    : " + num_fix_valid_veg)
 print("number of cumulative valid fix soil type                   : " + num_fix_valid_soil)
 print("number of cumulative missing fix veg type with land_frac>0 : "+num_oro_land_frac_gt0_and_no_fix_veg)
 print("number of cumulative missing fix soil type with land_frac>0: "+num_oro_land_frac_gt0_and_no_fix_soil)
 print("number of cumulative missing sfc veg type with land_frac>0 : "+num_oro_land_frac_gt0_and_no_sfc_veg)
 print("number of cumulative missing sfc soil type with land_frac>0: "+num_oro_land_frac_gt0_and_no_sfc_soil)
 print("number of cumulative fix land type = 17                    : "+num_fix_veg_17)

 print("number of cumulative sfc_slmsk = 1      : "+num_sfc_slmsk)
 print("number of cumulative valid sfc veg type : "+num_sfc_valid_veg)
 print("number of cumulative valid sfc soil type: "+num_sfc_valid_soil)
 
 print("number of cumulative sfc veg  type >0  : "+num_sfc_valid_veg)
 print("number of cumulative sfc soil type /=14: "+num_sfc_valid_soil)
 print("number of cumulative sfc veg  type =0  : "+num_sfc_veg_0)
 print("number of cumulative sfc soil type =0  : "+num_sfc_soil_0)
 print("number of cumulative sfc veg  type =17 : "+num_sfc_veg_17)
 print("number of cumulative sfc soil type =14 : "+num_sfc_soil_14)

 print("number of cumulative mismatch veg between fix and sfc  : "+fix_veg_sfc_veg_check)
 print("number of cumulative mismatch soil between fix and sfc : "+fix_soil_sfc_soil_check)
 print("number of cumulative mismatch fix_veg and sfc_soil     : "+fix_veg_sfc_soil_check)
 print("number of cumulative mismatch fix_soil and sfc_veg     : "+fix_soil_sfc_veg_cheeck)
 print("number of cumulative mismatch oro_slmsk and sfc_veg    : "+oro_slmsk_sfc_veg_check)
 print("number of cumulative mismatch oro_slmsk and sfc_soil   : "+oro_slmsk_sfc_soil_check)

end do

end
