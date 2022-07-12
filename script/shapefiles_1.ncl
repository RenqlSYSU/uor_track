;*************************************************
; shapefiles_1.ncl
;
; Concepts illustrated:
;   - Reading shapefiles
;   - Plotting data from shapefiles
;   - Using shapefile data to plot unemployment percentages in the U.S.
;   - Drawing a custom labelbar on a map
;   - Drawing filled polygons over a Lambert Conformal plot
;   - Drawing the US with a Lambert Conformal projection
;   - Zooming in on a particular area on a Lambert Conformal map
;   - Centering the labels under the labelbar boxes
;
;*************************************************
; 
; Simple example of how to draw selected geometry from a shapefile,
; based upon properties of an associated non-spatial variable.
;
; This example color-fills the states based upon "percent unemployment",
; which is calculated from several of the non-spatial variables in the
; file.
; 
;   "states.shp" is from the National Atlas (http://www.nationalatlas.gov/)
;  
; You must also have the files "states.dbf" and "states.shx" for this
; example to run.
;*************************************************
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
  shp_file = "/home/lzhenn/array74/tracacode/UTILITY-2016/shp/cnmap/gadm36_CHN_2.shp"
  f = addfile(shp_file, "r")   ; Open shapefile
  
  segments = f->segments
  geometry = f->geometry
  segsDims = dimsizes(segments) ;2D
  geomDims = dimsizes(geometry) ;2D

  geom_segIndex = f@geom_segIndex ;=0
  geom_numSegs  = f@geom_numSegs  ;=1
  segs_xyzIndex = f@segs_xyzIndex ;=0
  segs_numPnts  = f@segs_numPnts  ;=1

  lines       = new(segsDims(0),graphic)   ; Array to hold polygons
  numFeatures = geomDims(0) ;=dimsizes(NAME_0)

  lon   = f->x
  lat   = f->y
  name  = f->NAME_2

  segNum = 0
  do i=0, numFeatures-1  
     startSegment = geometry(i, geom_segIndex)
     numSegments  = geometry(i, geom_numSegs)
     print(name(i)+" startSegment:"+startSegment+" numSegments:"+numSegments)
     do seg=startSegment, startSegment+numSegments-1
        startPT = segments(seg, segs_xyzIndex)
        endPT = startPT + segments(seg, segs_numPnts) - 1
        print("startPT:"+startPT+" endPT:"+endPT) 
        ;+" lon(startPT:endPT):"+lon(startPT:endPT)+"lat(startPT:endPT):"+lat(startPT:endPT))
        segNum = segNum + 1
     end do
  end do
end

