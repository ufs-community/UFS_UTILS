load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  

begin

latnames = (/"+90","+60"/)
lonnames = (/"-180","-150","-120","-090","-060","-030","+000","+030","+060","+090","+120","+150","+180"/)
wtshal = (/0.045,0.046,0.075,0.123,0.204,0.336,0.171/)
wtdeep = (/0.383,0.617/)

do ilat = 0,0
do ilon = 0,11

  latbeg = 0
  latend = 24*120 - 1   ; original data available to 84N
  lonbeg = ilon*30*120
  lonend = ilon*30*120 + 30*120 - 1

  filename = "lat_"+latnames(ilat)+"_"+latnames(ilat+1)+"_lon_"+lonnames(ilon)+"_"+lonnames(ilon+1)
  print("Starting: "+filename)

  infile = addfile("original/CLAY1.nc","r")
    inclay1 = infile->CLAY(0,latbeg:latend,lonbeg:lonend)
    inclay2 = infile->CLAY(1,latbeg:latend,lonbeg:lonend)
    inclay3 = infile->CLAY(2,latbeg:latend,lonbeg:lonend)
    inclay4 = infile->CLAY(3,latbeg:latend,lonbeg:lonend)
  infile = addfile("original/CLAY2.nc","r")
    inclay5 = infile->CLAY(0,latbeg:latend,lonbeg:lonend)
    inclay6 = infile->CLAY(1,latbeg:latend,lonbeg:lonend)
    inclay7 = infile->CLAY(2,latbeg:latend,lonbeg:lonend)
    inclay8 = infile->CLAY(3,latbeg:latend,lonbeg:lonend)
  infile = addfile("original/SAND1.nc","r")
    insand1 = infile->SAND(0,latbeg:latend,lonbeg:lonend)
    insand2 = infile->SAND(1,latbeg:latend,lonbeg:lonend)
    insand3 = infile->SAND(2,latbeg:latend,lonbeg:lonend)
    insand4 = infile->SAND(3,latbeg:latend,lonbeg:lonend)
  infile = addfile("original/SAND2.nc","r")
    insand5 = infile->SAND(0,latbeg:latend,lonbeg:lonend)
    insand6 = infile->SAND(1,latbeg:latend,lonbeg:lonend)
    insand7 = infile->SAND(2,latbeg:latend,lonbeg:lonend)
    insand8 = infile->SAND(3,latbeg:latend,lonbeg:lonend)
  infile = addfile("original/SILT1.nc","r")
    insilt1 = infile->SILT(0,latbeg:latend,lonbeg:lonend)
    insilt2 = infile->SILT(1,latbeg:latend,lonbeg:lonend)
    insilt3 = infile->SILT(2,latbeg:latend,lonbeg:lonend)
    insilt4 = infile->SILT(3,latbeg:latend,lonbeg:lonend)
  infile = addfile("original/SILT2.nc","r")
    insilt5 = infile->SILT(0,latbeg:latend,lonbeg:lonend)
    insilt6 = infile->SILT(1,latbeg:latend,lonbeg:lonend)
    insilt7 = infile->SILT(2,latbeg:latend,lonbeg:lonend)
    insilt8 = infile->SILT(3,latbeg:latend,lonbeg:lonend)

sand_top = new((/3600,3600/),float)
silt_top = new((/3600,3600/),float)
clay_top = new((/3600,3600/),float)
sand_bot = new((/3600,3600/),float)
silt_bot = new((/3600,3600/),float)
clay_bot = new((/3600,3600/),float)

; 720 here because 6 (degrees) * 120 (pixels per degree) for 84N - 90N

sand_top(latbeg+720:latend+720,:) = wtshal(0) * insand1 + wtshal(1) * insand2 + wtshal(2) * insand3 + wtshal(3) * insand4 + wtshal(4) * insand5 + wtshal(5) * insand6 + wtshal(6) * insand7
silt_top(latbeg+720:latend+720,:) = wtshal(0) * insilt1 + wtshal(1) * insilt2 + wtshal(2) * insilt3 + wtshal(3) * insilt4 + wtshal(4) * insilt5 + wtshal(5) * insilt6 + wtshal(6) * insilt7
clay_top(latbeg+720:latend+720,:) = wtshal(0) * inclay1 + wtshal(1) * inclay2 + wtshal(2) * inclay3 + wtshal(3) * inclay4 + wtshal(4) * inclay5 + wtshal(5) * inclay6 + wtshal(6) * inclay7

sand_top = mask(sand_top,sand_top.ge.0.and.sand_top.le.100,True)
silt_top = mask(silt_top,silt_top.ge.0.and.silt_top.le.100,True)
clay_top = mask(clay_top,clay_top.ge.0.and.clay_top.le.100,True)

sand_bot(latbeg+720:latend+720,:) = wtdeep(0) * insand7 + wtdeep(1) * insand8
silt_bot(latbeg+720:latend+720,:) = wtdeep(0) * insilt7 + wtdeep(1) * insilt8
clay_bot(latbeg+720:latend+720,:) = wtdeep(0) * inclay7 + wtdeep(1) * inclay8

sand_bot = mask(sand_bot,sand_bot.ge.0.and.sand_bot.le.100,True)
silt_bot = mask(silt_bot,silt_bot.ge.0.and.silt_bot.le.100,True)
clay_bot = mask(clay_bot,clay_bot.ge.0.and.clay_bot.le.100,True)

soil_texture_top = new((/3600,3600/),integer)
  soil_texture_top@_FillValue = -999
soil_texture_bot = new((/3600,3600/),integer)
  soil_texture_bot@_FillValue = -999

; SOILPARM:1  sand

soil_texture_top = where((silt_top + 1.5*clay_top) .lt. 15, 1, soil_texture_top) 
soil_texture_bot = where((silt_bot + 1.5*clay_bot) .lt. 15, 1, soil_texture_bot) 

; SOILPARM:2  loamy sand

soil_texture_top = where((silt_top + 1.5*clay_top .ge. 15) .and. (silt_top + 2*clay_top .lt. 30), 2, soil_texture_top)
soil_texture_bot = where((silt_bot + 1.5*clay_bot .ge. 15) .and. (silt_bot + 2*clay_bot .lt. 30), 2, soil_texture_bot)

; SOILPARM:3  sandy loam

soil_texture_top = where((clay_top .ge. 7 .and. clay_top .lt. 20) .and. (sand_top .gt. 52) .and. ((silt_top + 2*clay_top) .ge. 30) .or. (clay_top .lt. 7 .and. silt_top .lt. 50 .and. (silt_top+2*clay_top).ge.30), 3, soil_texture_top)
soil_texture_bot = where((clay_bot .ge. 7 .and. clay_bot .lt. 20) .and. (sand_bot .gt. 52) .and. ((silt_bot + 2*clay_bot) .ge. 30) .or. (clay_bot .lt. 7 .and. silt_bot .lt. 50 .and. (silt_bot+2*clay_bot).ge.30), 3, soil_texture_bot)

; SOILPARM:6  loam

soil_texture_top = where((clay_top .ge. 7 .and. clay_top .lt. 27) .and. (silt_top .ge. 28 .and. silt_top .lt. 50) .and. (sand_top .le. 52), 6, soil_texture_top)
soil_texture_bot = where((clay_bot .ge. 7 .and. clay_bot .lt. 27) .and. (silt_bot .ge. 28 .and. silt_bot .lt. 50) .and. (sand_bot .le. 52), 6, soil_texture_bot)

; SOILPARM:4  silt loam

soil_texture_top = where((silt_top .ge. 50 .and. (clay_top .ge. 12 .and. clay_top .lt. 27)) .or. ((silt_top .ge. 50 .and. silt_top .lt. 80) .and. clay_top .lt. 12), 4, soil_texture_top)
soil_texture_bot = where((silt_bot .ge. 50 .and. (clay_bot .ge. 12 .and. clay_bot .lt. 27)) .or. ((silt_bot .ge. 50 .and. silt_bot .lt. 80) .and. clay_bot .lt. 12), 4, soil_texture_bot)

; SOILPARM:5  silt

soil_texture_top = where(silt_top .ge. 80 .and. clay_top .lt. 12, 5, soil_texture_top)
soil_texture_bot = where(silt_bot .ge. 80 .and. clay_bot .lt. 12, 5, soil_texture_bot)

; SOILPARM:7  sandy clay loam

soil_texture_top = where((clay_top .ge. 20 .and. clay_top .lt. 35) .and. (silt_top .lt. 28) .and. (sand_top .gt. 45), 7, soil_texture_top) 
soil_texture_bot = where((clay_bot .ge. 20 .and. clay_bot .lt. 35) .and. (silt_bot .lt. 28) .and. (sand_bot .gt. 45), 7, soil_texture_bot) 

; SOILPARM:9  clay loam

soil_texture_top = where((clay_top .ge. 27 .and. clay_top .lt. 40) .and. (sand_top .gt. 20 .and. sand_top .le. 45), 9, soil_texture_top)
soil_texture_bot = where((clay_bot .ge. 27 .and. clay_bot .lt. 40) .and. (sand_bot .gt. 20 .and. sand_bot .le. 45), 9, soil_texture_bot)

; SOILPARM:8  silty clay loam

soil_texture_top = where((clay_top .ge. 27 .and. clay_top .lt. 40) .and. (sand_top  .le. 20), 8, soil_texture_top)
soil_texture_bot = where((clay_bot .ge. 27 .and. clay_bot .lt. 40) .and. (sand_bot  .le. 20), 8, soil_texture_bot)

; SOILPARM:10 sandy clay

soil_texture_top = where(clay_top .ge. 35 .and. sand_top .gt. 45, 10, soil_texture_top)
soil_texture_bot = where(clay_bot .ge. 35 .and. sand_bot .gt. 45, 10, soil_texture_bot)

; SOILPARM:11 silty clay

soil_texture_top = where(clay_top .ge. 40 .and. silt_top .ge. 40, 11, soil_texture_top)
soil_texture_bot = where(clay_bot .ge. 40 .and. silt_bot .ge. 40, 11, soil_texture_bot)

; SOILPARM:12 clay

soil_texture_top = where(clay_top .ge. 40 .and. sand_top .le. 45 .and. silt_top .lt. 40, 12, soil_texture_top)
soil_texture_bot = where(clay_bot .ge. 40 .and. sand_bot .le. 45 .and. silt_bot .lt. 40, 12, soil_texture_bot)

soil_texture_top = mask(soil_texture_top,soil_texture_top.ge.0,True)
soil_texture_bot = mask(soil_texture_bot,soil_texture_bot.ge.0,True)

outname = "tiles/"+filename+".nc"
 system("if [ -e "+outname+" ]; then rm -f "+outname+ ";fi")
outfile = addfile(outname,"c")
outfile->soil_texture_top = soil_texture_top(::-1,:)
outfile->soil_texture_bot = soil_texture_bot(::-1,:)

delete(sand_top)
delete(silt_top)
delete(clay_top)
delete(sand_bot)
delete(silt_bot)
delete(clay_bot)
delete(soil_texture_top)
delete(soil_texture_bot)

end do
end do

end

