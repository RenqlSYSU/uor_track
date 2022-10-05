;*************************************************
; shapefiles_4.ncl
;
; Concepts illustrated:
;   - Drawing the Mississippi River Basin using data from a shapefile
;   - Masking a data array based on a geographical area obtained from a shapefile
;   - Attaching markers to a map
;   - Attaching polylines to a map plot
;
;*************************************************
; This script shows the "new" way (post NCL V6.0.0) of  masking
; data and adding shapefile outlines to an existing NCL map.
;*************************************************
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;
; This file still has to be loaded manually
load "/disk/hq247/yhuangci/resource/ncl-resource/shapefile_utils.ncl"

begin
    ;grid file
    lat1d = ispan(1800,5400 ,25)/100.0
    lon1d = ispan(7300,13500,25)/100.0
    lat1d@units = "degrees_north"
    lon1d@units = "degrees_east"
    nlat = dimsizes(lat1d)
    nlon = dimsizes(lon1d)
    
    data = new((/nlat,nlon/),float) 
    data = 1
    data@_FillValue=-9999
    data!0      = "lat"
    data!1      = "lon"
    data&lat    = lat1d
    data&lon    = lon1d

    lat2d = ndtooned(conform_dims((/nlat,nlon/),lat1d,0))
    lon2d = ndtooned(conform_dims((/nlat,nlon/),lon1d,1))
    lat2d@_FillValue = -9999
    lon2d@_FillValue = -9999

;---Open shapefile and read lat/lon values.
  outfile  = "/home/lzhenn/cooperate/script/input_china_east_lv2.csv"
  shp_file = "/home/lzhenn/array74/tracacode/UTILITY-2016/shp/cnmap/gadm36_CHN_2.shp"
  f = addfile(shp_file,"r")
  city_name = f->NAME_2
  prov_name = f->NAME_1

  f_sp = addfile("/home/lzhenn/cooperate/data/ERA5_annual_sp_1979-2020.nc","r")
  sp = f_sp->sp ;Pa
  sp = sp/100.0-50.0 ;hPa

;---Set all hgt values to missing except for those over Ohio.
  opt             = True
  opt@debug       = True
  opt@shape_var   = "NAME_2"

numb = 0
write_table(outfile,"w",[/"index,lat,lon,pressure",",city,prov"/],"%s%s")
do i = 0,dimsizes(city_name)-1
  if (prov_name(i).eq."Qinghai".or.prov_name(i).eq."Xizang" \
      .or.prov_name(i).eq."Xinjiang Uygur") then 
    print("skip "+prov_name(i))
    continue
  end if
  fmft = "%d,%5.2f,%6.2f,%6.2f,"+city_name(i)+","+prov_name(i)
  opt@shape_names = city_name(i) 
  data_mask  = shapefile_mask_data(data,shp_file,opt)
  indx := ind(.not.ismissing(ndtooned(data_mask)))
  npot = dimsizes(indx)
  do np = 0,npot-1,1
    write_table(outfile,"a",[/numb,lat2d(indx(np)),lon2d(indx(np)),sp({lat2d(indx(np))},{lon2d(indx(np))})/],fmft)
    numb = numb + 1 
  end do
  ;write_table(outfile,"a",[/ispan(numb,numb+npot-1,1),lat2d(indx),lon2d(indx)/],fmft)
end do

end

