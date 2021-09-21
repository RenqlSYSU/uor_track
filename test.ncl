begin
lonl=0  
lonr=150
lats=15 
latn=70 
months   = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)

path = "/home/users/qd201969/data/"
f=addfile(path+"ERA5_mon_u_1979-2020.nc","r")
var = short2flt(f->u(:,0,{200},{lats:latn},{lonl:lonr}))
printVarSummary(var)

plot   = new((/12/),graphic)
;define the thickness and height of the title,tick and label
   tick_font      = 0.03
   subtitle_font  = 0.03
   title_font     = 0.03
   vcref_font     = 0.015
   tick_thick     = 8.0
   subtitle_thick = 8.0  ;not have this attribute
   title_thick    = 5.0
   vc_thick       = 3.0
   mp_thick       = 2.5;3.0;
   font_index     = 22

wks = gsn_open_wks("png","test")

    res                  = True
    res@gsnFrame         = False
    res@gsnDraw          = False
    res@gsnMaximize      = True
    res@gsnAddCyclic     = True;False
    res@gsnPaperOrientation = "landscape" ;"portrait"
    ;res@vpHeightF        = 0.25
    ;res@vpWidthF         = 0.8

;-------set the attribute of the map-------------------
    res@mpCenterLonF     = 180.
    res@mpMaxLatF        = latn
    res@mpMinLatF        = lats
    res@mpMaxLonF        = lonr
    res@mpMinLonF        = lonl
    res@gsnMajorLatSpacing = 20 
    res@gsnMajorLonSpacing = 30
    res@gsnAddCyclic     = False
    res@mpFillOn         = False
    res@mpGeophysicalLineThicknessF = mp_thick
    res@mpGeophysicalLineColor      = "black";"white"

;--------set the form of the coordinate-----------------
    res@tmBorderThicknessF  = mp_thick
    res@tmXBLabelFont   = font_index
    res@tmYLLabelFont   = font_index
    res@tmXBLabelFontHeightF     = tick_font
    res@tmXBLabelFontThicknessF  = tick_thick
    res@tmYLLabelFontHeightF     = tick_font
    res@tmYLLabelFontThicknessF  = tick_thick
    
    res@tmXBMajorThicknessF     = mp_thick
    res@tmYLMajorThicknessF     = mp_thick
    res@tmXBMinorThicknessF     = mp_thick
    res@tmYLMinorThicknessF     = mp_thick

;---------set the attribute of the contour-----------------
    res@cnFillOn         = True
    res@cnLinesOn        = False ;True
    res@cnLineLabelsOn   = False
    res@cnInfoLabelOn    = False
    res@cnLevelSelectionMode  = "ManualLevels" ; set manual contour levels
    res@cnMinLevelValF        = 20          ; set min contour level
    res@cnMaxLevelValF        = 52 ;32        ; set max contour level
    res@cnLevelSpacingF       = 2        ; set contour interval
    res@lbLabelBarOn         = False   ; the label bar

    res@gsnSpreadColors      = False  ;True              ; turn of the full color map
    res@gsnSpreadColorStart  = 11                ; start at color 17
    res@gsnSpreadColorEnd    = 2                ; end at color 200
    res@gsnCenterString      = ""
    res@gsnStringFont         = font_index
    res@gsnStringFontHeightF = subtitle_font
    res@cnFillPalette = "precip2_17lev"

time = cd_calendar(f->time,0)
do nm = 0,11,1
    select_time := ind(time(:,1).eq.nm+1)
    var1 = dim_avg_n_Wrap(var(select_time,:,:),0) ;get (nyear,nmonth,nlat,nlon)

    res@gsnLeftString   = months(nm)
    plot(nm)    = gsn_csm_contour_map_ce(wks,var1,res)
end do
printVarSummary(var1)
    
resp  = True    
resp@gsnMaximize   = True
resp@gsnPaperOrientation = "landscaoe";"portrait"
resp@gsnPanelLabelBar = True
resp@lbOrientation     = "Horizontal" ;"Vertical"
resp@pmLabelBarWidthF  = 0.8 ;0.05
resp@pmLabelBarHeightF = 0.05;0.8
resp@lbLabelFont        = font_index              ; make labels smaller
resp@lbLabelFontHeightF = 0.015 
resp@lbLabelFontThicknessF = tick_thick
resp@txFontHeightF = title_font 
resp@txFontThicknessF = title_thick
gsn_panel(wks,plot,(/4,3/),resp)
end 
