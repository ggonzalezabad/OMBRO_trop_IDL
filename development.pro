SPAWN, 'ls ../L3_he5_files_OMI_BrO_v2/OMBRO_L3_2005*', files, count=nfiles

; Define Area of study:
latmin = 35.0 & latmax = 55.0
lonmin = 55.0 & lonmax = 70.0

; Variable to hold final data over Aral Sea
final_data = FLTARR(7,nfiles)

; Aral Sea box
aral_lon_min = 58.0 & aral_lon_max = 62.0
aral_lat_min = 43.0 & aral_lat_max = 47.0

; Sigma threshold value
factor = 0.0

; Smooth width
width = [15,15]

; Initial margins
x_margin = !x.margin
y_margin = !y.margin

OPEN_DEVICE, /PS, /COLOR, filename='Development_v2.ps'
myct, 22

For ifile = 0, nfiles-1 do begin

print, files[ifile]
omi_read_gridded_l2_he5, files[ifile], swathname, nlon, nlat, dlon, dlat, $
                         colcor=colcor, colerr=colerr, rms=rms,       $
                         numsam=numsam, lats=lats, lons=lons


;Obtain mask for area of study
mask1 = WHERE(lats GE latmin AND lats LE latmax)
mask2 = WHERE(lons GE lonmin AND lons LE lonmax)

dummy_lat = FLTARR(n_elements(mask2),n_elements(mask1))
For i = 0, n_elements(mask2)-1 do dummy_lat[i,*] = lats[mask1]
dummy_lon = FLTARR(n_elements(mask2),n_elements(mask1))
For i = 0, n_elements(mask1)-1 do dummy_lon[*,i] = lons[mask2]
dummy_col = colcor[min(mask2):max(mask2),   $
                   min(mask1):max(mask1)]

multipanel, cols=2, rows=1
TVMAP, dummy_col, lons[mask2], lats[mask1], /cbar, /sample, /noadvance, $
       /continents, /coasts, mindata = 0.0, title = strmid(files[ifile],33,9), $
       maxdata = max(dummy_col), /countries

OPLOT, [aral_lon_min,aral_lon_max,aral_lon_max,aral_lon_min,aral_lon_min], $
       [aral_lat_min,aral_lat_min,aral_lat_max,aral_lat_max,aral_lat_min], $
       color = 1

!x.margin = x_margin
!y.margin = y_margin

; Work out mean, median and 1 sigma for area of study
mask7 = WHERE(dummy_col GT -1e29)
print, mean(dummy_col[mask7]), median(dummy_col[mask7]), stddev(dummy_col[mask7])
test = stddev(dummy_col[mask7])

; Define latitude grids
lat_resolutions = [0.5]
latgrid         = FINDGEN(( (latmax-latmin)/lat_resolutions)+1.0) * lat_resolutions[0] + latmin
nlatgrid        = n_elements(latgrid)

; Create variable to save, lat_mean, median and 1 sigma, and median+1sigma
data_lat_mean = FLTARR(5,nlatgrid,n_elements(lat_resolutions))

For iresolution = 0, n_elements(lat_resolutions)-1 do begin
   
   For ilat = 0, nlatgrid-2 do begin
      ; Define new masks
      mask3 = WHERE(dummy_lat GE latgrid[ilat] AND dummy_lat LE latgrid[ilat+1] AND $
                    dummy_col GT -1e29, count)

      IF (COUNT LT 3) THEN CONTINUE
      ; Work out mean, median and 1 sigma for area of study
      data_lat_mean[0,ilat,iresolution] = mean(dummy_lat[mask3])
      data_lat_mean[1,ilat,iresolution] = median(dummy_col[mask3])
      data_lat_mean[2,ilat,iresolution] = STDDEV(dummy_col[mask3])
      data_lat_mean[3,ilat,iresolution] = mean(dummy_col[mask3])

;      print, histogram(colcor[min(mask2):max(mask2),min(mask3):max(mask3)],binsize = 0.05e13)

      ; For each lat grid find threshold value
      data_lat_mean[3,ilat,iresolution] = data_lat_mean[3,ilat,iresolution] + $
                                          factor*data_lat_mean[2,ilat,iresolution]

      ; Find pixels that fall with in that threshold value
      mask3 = WHERE(dummy_lat GE latgrid[ilat] AND dummy_lat LE latgrid[ilat+1] AND $
                    dummy_col GT data_lat_mean[3,ilat,iresolution])

      ;; print, 'Latitudes ='+ $
      ;;        STRING(latgrid[ilat],   format = '(F5.1)') + $
      ;;        STRING(latgrid[ilat+1], format = '(F5.1)') +' ; '+ $
      ;;        STRING(data_lat_mean[1,ilat,iresolution], format ='(E10.2)') + $
      ;;        STRING(data_lat_mean[2,ilat,iresolution], format ='(E10.2)') + $
      ;;        STRING(n_elements(mask3), FORMAT = '(I4)')+ $
      ;;        STRING(count, FORMAT = '(I6)')

   Endfor

Endfor

; Start plot with empty lat, column
good_plot = WHERE(data_lat_mean[1,*,0] GT 0)
cgplot, data_lat_mean[1,good_plot,0], $
        data_lat_mean[0,good_plot,0], $
        err_xlow = data_lat_mean[2,good_plot,0], $
        err_xhigh = data_lat_mean[2,good_plot,0], $
        color = 1, ytitle = 'Latitude', xtitle = 'BrO [molecules cm!u-2!n]', $
        xstyle = 1, ystyle = 1, xrange = [1e12,10e13], yrange = [latmin,latmax]

cgplot, data_lat_mean[3,good_plot,0], $
        data_lat_mean[0,good_plot,0], $
        err_xlow = data_lat_mean[2,good_plot,0], $
        err_xhigh = data_lat_mean[2,good_plot,0], $
        color = 2, /Overplot

multipanel, cols=3, rows=2
; For each resolution use mask and plot results on a pixel plot
For iresolution = 0, n_elements(lat_resolutions)-1 do begin


   TVMAP, dummy_col, lons[mask2], lats[mask1], $
          /sample, /cbar, $
          /continents, /coasts, /countries, mindata=0, maxdata=max(dummy_col)
   
   dummy_mask = dummy_col*0.0
   dummy_back = dummy_col*!Values.F_nan

   For ilat = 0, nlatgrid-2 do begin
      mask4 = WHERE(dummy_lat GE latgrid[ilat] AND dummy_lat LE latgrid[ilat+1] AND $
                    dummy_col GT data_lat_mean[3,ilat,iresolution])
      dummy_mask[mask4] = 1.0
      mask4 = WHERE(dummy_col LT -1e29)
      dummy_mask[mask4] = 1.0
   Endfor
   
   TVMAP, dummy_mask, lons[mask2], lats[mask1], $
          /sample, /cbar, $
          /continents, /coasts, /countries

   For ilat = 0, nlatgrid-2 do begin
      mask5 = WHERE(dummy_lat GE latgrid[ilat] AND dummy_lat LE latgrid[ilat+1] AND $
                    dummy_mask EQ 0 AND dummy_col GT -1e29, count)
      mask6 = WHERE(dummy_lat GE latgrid[ilat] AND dummy_lat LE latgrid[ilat+1])
      IF (COUNT LT 3) THEN CONTINUE
      dummy_back[mask6] = MEAN(dummy_col[mask5])
   Endfor

   TVMAP, dummy_back, lons[mask2], lats[mask1], $
          /sample, /cbar, $
          /continents, /coasts, /countries, $
          mindata = 0, maxdata = max(dummy_col)

   mask7 = WHERE(dummy_col GT -1e29)
   dummy_final = dummy_col & dummy_final[mask7] = dummy_col[mask7]-dummy_back[mask7]

   mask_smooth = WHERE(dummy_mask EQ 1) 
   dummy_smooth = dummy_col & dummy_smooth[mask_smooth] = !Values.F_nan
   mask_smooth = WHERE(dummy_smooth EQ -1e30) & dummy_smooth[mask_smooth] = !Values.F_Nan
   dummy_back_smooth = smooth(dummy_smooth,width,/nan,/Edge_mirror)
   dummy_final_smooth = dummy_col
   dummy_final_smooth[mask7] = dummy_col[mask7]-dummy_back_smooth[mask7]

   maxdata = MAX(dummy_final)
   IF (max(dummy_final_smooth) GT maxdata) THEN maxdata = max(dummy_final_smooth)

   TVMAP, dummy_final, lons[mask2], lats[mask1], $
          /sample, /cbar, /continents, /coasts, /countries, mindata = 0, $
          maxdata = 0.8e13, $
          limit = [41,56,49,64]

   TVMAP, dummy_back_smooth, lons[mask2], lats[mask1], $
          /sample, /cbar, $
          /continents, /coasts, /countries, $
          mindata = 0, maxdata = max(dummy_col)

   TVMAP, dummy_final_smooth, lons[mask2], lats[mask1], $
          /sample, /cbar, /continents, /coasts, /countries, mindata = 0, $
          maxdata = 0.8e13, $
          limit = [41,56,49,64]

   ;; myct, 63, /reverse
   ;; TVMAP, (dummy_back-dummy_back_smooth)/dummy_back*100.0, lons[mask2], lats[mask1], $
   ;;        /sample, /cbar, /continents, /coasts, /countries, mindata = -20.0, maxdata = 20.0, $
   ;;        limit = [41,56,49,64]
   ;; TVMAP, dummy_back-dummy_back_smooth, lons[mask2], lats[mask1], $
   ;;        /sample, /cbar, /continents, /coasts, /countries, mindata=-1e13, $
   ;;        maxdata=1e13, $
   ;;        limit = [41,56,49,64]

   ;; TVMAP, dummy_final-dummy_final_smooth, lons[mask2], lats[mask1], $
   ;;        /sample, /cbar, /continents, /coasts, /countries, $
   ;;        limit = [41,56,49,64]


Endfor

final_mask = WHERE(dummy_lat GE aral_lat_min AND dummy_lat LE aral_lat_max AND $
                   dummy_lon GE aral_lon_min AND dummy_lon LE aral_lon_max AND $
                   dummy_final GE -1e29 AND dummy_final_smooth GE -1e29, count)

mask_data = WHERE(dummy_mask EQ 0, count_mask)
final_data[6,ifile] = count_mask

IF (COUNT LT 1) THEN CONTINUE
final_data[0,ifile] = count
final_data[1,ifile] = MEAN(dummy_final[final_mask])
final_data[2,ifile] = MEDIAN(dummy_final[final_mask])
final_data[3,ifile] = MEAN(dummy_final_smooth[final_mask])
final_data[4,ifile] = MEDIAN(dummy_final_smooth[final_mask])
final_data[5,ifile] = STDDEV(dummy_final[final_mask])


multipanel, /off
!x.margin = x_margin
!y.margin = y_margin
Endfor
multipanel, cols=1, rows=1
mask = WHERE(final_data[0,*] NE 0, count)
cgplot, findgen(count)+1, final_data[1,mask], color = 1, xtitle ='Day #', $
        ytitle = 'BrO [molecules cm!u-2!n]'
;gplot, findgen(count)+1, final_data[2,mask], color = 2, /Overplot
cgplot, findgen(count)+1, final_data[3,mask], color = 3, /Overplot
;gplot, findgen(count)+1, final_data[4,mask], color = 4, /Overplot
;print, final_data[6,mask]
multipanel, /off
CLOSE_DEVICE

End

