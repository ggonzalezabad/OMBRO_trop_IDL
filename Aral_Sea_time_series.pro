; This script intends to create a long term time series of Aral Sea
; tropospheric BrO columns inferred from OMBRO data. OMBRO data has
; been gridded to 0.2 by 0.2 degrees. Cloud fraction is below 0.30,
; xtrack pixels considered go from 4 to 56, pixels affected by the row
; anomally are excluded, and contributions are weighted by are and fitting
; uncertainty. The time series is build using daily monthly running
; means extending from 2005 till 2015.

;The idea is for each file compute a background for an extended region
;that is going to be representative of the stratospheric column and
;then we substract. I'm applying a box-car smoothing filter to
;obtain the background. Hopefully comments along the code will help to
;understand each step.

; Gonzalo Gonzalez Abad, February 2016, SAO
; ggonzalezabad@cfa.harvard.edu

PRO Aral_Sea_time_series

; Obtain list of L3 files to be processed
SPAWN, 'ls ../L3_he5_files_OMI_BrO/OMBRO_L3_20??m*_0p2x0p2_SZA85_XT4-56_CF30.he5', $
       files, count=nfiles

; Index of date string start location
date_idx = 33

; Define geographic region to study. I'm giving some room
; around the Aral Sea to allow for a more effective use of the
; smoothing kernel
lonmin = -170.0 & lonmax = 170.0
latmin =  -85.0 & latmax = 85.0

; Define geographic region of the Aral Sea
aral_lon_min = 58.5 & aral_lon_max = 61.5
aral_lat_min = 43.5 & aral_lat_max = 46.5

; Select other random region
rand_lon_min = -120.0 & rand_lon_max = -115.0
rand_lat_min =   -2.5 & rand_lat_max =    2.5

; Sigma threshold factor. This factor applies to the following
; expression which defines the masking of the points that will not be
; used to compute the background (stratospheric column of the region
; of interest) for a given latitude bin:
; mask = mean(lat) + sigma_factor * stardard_deviation(lat)
sigma_factor = 0.0

; How many pixels to be considered in the box car smoothing
; function. Each pixel has a resolution of 0.2x0.2 degrees. Trying to
; find a balance between the expected change of the stratospheric
; column in the horizontal plane and the efficiency of the smoothing I
; have decided to use +-0.75 degrees or 15 pixels in each
; dimenssion. This assumption can be revisited
smooth_width = [400,20]

; Define latitude grid to work out background correction
lat_resolution = 0.5
latgrid        = FINDGEN(( (latmax-latmin)/lat_resolution)+1.0) * $
                 lat_resolution[0] + latmin
nlatgrid       = n_elements(latgrid)

; L3 missing value
miss_val = -1.0e30

; Create variables to hold time series data for the Aral Sea
; Number of points, mean, variance, skewness, kurtosis and standard
; deviation of tropospheric columns and mean and standard deviation of
; stratospheric columns
Aral_final_data = FLTARR(8,nfiles) * !Values.F_nan
Rand_final_data = FLTARR(8,nfiles) * !Values.F_nan
date_int = LONARR(nfiles)
date_str = STRARR(nfiles)

; Create variable to save BrO mean, standard deviation, number of
; points included in background, threshold value, and number of pixels
; contributing to background for each latitude bin and for each file
Background_lat = FLTARR(5,nlatgrid-1,nfiles) * !Values.F_nan

; Create variable to store smooth stratospheric background and tropospheric
; VCDs. First I open a L3 file to get the dimensions
omi_read_gridded_l2_he5, files[0], swathname, nlon, nlat, dlon, dlat, $
                         colcor=colcor
stratospheric_vcds = FLTARR(nlon,nlat,nfiles) * !Values.F_nan
tropospheric_vcds  = FLTARR(nlon,nlat,nfiles) * !Values.F_nan

; Keep original plotting margins to prevent extrange behaviour
x_margin = !x.margin & y_margin = !y.margin

; Open postcript document to save plots
OPEN_DEVICE, /PS, /COLOR, filename = 'Aral_Sea_time_series_plots_0p2x0p2_SZA70_XT1-23_CF30_global.ps'

; Load color table, Hue Sat Value 2
MYCT, 22

; Loop over L3 files
FOR ifile = 0, nfiles-1 do begin
   ; Print control line to track program running
   print, nfiles-ifile, files[ifile]
   
   ; Fill in values of date arrays, STRMID positions will 
   ; need to be updated depending on the file name and
   ; path to folder containing L3 files. Defined above
   date_int[ifile] = LONG( STRMID(files[ifile], date_idx, 4)+ $
                           STRMID(files[ifile], date_idx+5,4))
   date_str[ifile] = STRMID(files[ifile], date_idx, 9)

   ; Read L3 file: colcor = Destriped column
   ;               colerr = Fitting uncertainty
   ;               rms    = Root mean Square of the fitting
   ;               numsam = Number of samples used in each pixel
   ;               alt    = Surface altitude
   omi_read_gridded_l2_he5, files[ifile], swathname, nlon, nlat,      $
                            dlon, dlat, colcor=colcor, colerr=colerr, $
                            rms=rms, numsam=numsam, alt=alt,          $
                            lats=lats, lons=lons

   ; Create 2D arrays for longitudes and latitudes.
   ; It simplifies things later on. A bit slower thouhg
   lon_2D = FLTARR(nlon,nlat)
   lat_2D = FLTARR(nlon,nlat)
   For i = 0, nlat-1 do lon_2D[*,i] = lons
   For i = 0, nlon-1 do lat_2D[i,*] = lats

   ; Create 2D array to keep mask for background values
   mask_2D = INTARR(nlon,nlat)

   ; Create 2D array to keep stratospheric background values
   strat_back = FLTARR(nlon,nlat) * !Values.F_nan

   ; Work out mean, standard deviation and threshold level for each
   ; latitude grid
   For ilat = 0, nlatgrid-2 do begin
      
      ; Define mask
      mask = WHERE(lon_2D GE lonmin        AND lon_2D LE lonmax AND        $
                   lat_2D GE latgrid[ilat] AND lat_2D LE latgrid[ilat+1] AND $
                   colcor NE miss_val, count)

      ; Be sure we have at least 3 ponts in each latitudinal bin,
      ; otherwise skip it
      IF (COUNT LT 3) THEN CONTINUE
      ; Work out mean, median and 1 sigma for area of study
      Background_lat[0,ilat,ifile] = mean(colcor[mask])
      Background_lat[1,ilat,ifile] = STDDEV(colcor[mask])
      Background_lat[2,ilat,ifile] = count
      Background_lat[3,ilat,ifile] = Background_lat[0,ilat,ifile] + $
                                     sigma_factor * Background_lat[1,ilat,ifile]

      ; Using the values computed above to work out mask_2D
      mask = WHERE(lon_2D GE lonmin        AND lon_2D LE lonmax          AND $
                   lat_2D GE latgrid[ilat] AND lat_2D LE latgrid[ilat+1] AND $
                   colcor LT Background_lat[3,ilat,ifile]                AND $
                   colcor NE miss_val, count)
      IF (COUNT GE 1) THEN BEGIN
         mask_2D[mask] = 1
         Background_lat[4,ilat,ifile] = count
      ENDIF
   Endfor

   ; Work out values of stratospheric VCDs
   mask = WHERE (mask_2D EQ 1, count)
   strat_back[mask] = colcor[mask]

   ; Smooth stratospheric background
   stratospheric_vcds[*,*,ifile] = SMOOTH(strat_back,smooth_width,/nan,/Edge_truncate)

   ; Compute tropospheric column only for pixel where we have valid data
   mask = WHERE (colcor EQ miss_val)
   colcor[mask] = !Values.F_nan
   tropospheric_vcds[*,*,ifile] = colcor - stratospheric_vcds[*,*,ifile]

   ; Work out mean, variance, skewness, kurtosis, 
   ; and standard deviation for the are defined by
   ; aral_lon_min, aral_lon_max, aral_lat_min, and
   ; aral_lat_max of the tropospheric column.
   mask_lon = WHERE(lons GE aral_lon_min AND lons LE aral_lon_max)
   mask_lat = WHERE(lats GE aral_lat_min AND lats LE aral_lat_max)

   dummy = MOMENT(tropospheric_vcds[min(mask_lon):max(mask_lon), $
                                    min(mask_lat):max(mask_lat), $
                                    ifile], /nan)

   Aral_final_data[1,ifile] = dummy[0] ;Save mean trophosphere
   Aral_final_data[2,ifile] = dummy[1] ;Save variance trophosphere
   Aral_final_data[3,ifile] = dummy[2] ;Save skewness trophosphere
   Aral_final_data[4,ifile] = dummy[3] ;Save kurtosis trophosphere
   Aral_final_data[5,ifile] = SQRT(dummy[1]) ;Save standard deviation trophosphere

   dummy = WHERE(FINITE(tropospheric_vcds[min(mask_lon):max(mask_lon), $
                                          min(mask_lat):max(mask_lat), $
                                          ifile]),count)
   Aral_final_data[0,ifile] = count ; Save number of data points averaged over Aral Sea

   dummy = MOMENT(stratospheric_vcds[min(mask_lon):max(mask_lon), $
                                     min(mask_lat):max(mask_lat), $
                                     ifile], /nan)
   Aral_final_data[6,ifile] = dummy[0] ;Save mean stratosphere
   Aral_final_data[7,ifile] = SQRT(dummy[1]) ;Save standard deviation stratosphere

   ; Work out mean, variance, skewness, kurtosis, 
   ; and standard deviation for the are defined by
   ; Rand_lon_min, Rand_lon_max, rand_lat_min, and
   ; rand_lat_max of the tropospheric column.
   mask_lon = WHERE(lons GE rand_lon_min AND lons LE rand_lon_max)
   mask_lat = WHERE(lats GE rand_lat_min AND lats LE rand_lat_max)

   dummy = MOMENT(tropospheric_vcds[min(mask_lon):max(mask_lon), $
                                    min(mask_lat):max(mask_lat), $
                                    ifile], /nan)

   Rand_final_data[1,ifile] = dummy[0] ;Save mean trophosphere
   Rand_final_data[2,ifile] = dummy[1] ;Save variance trophosphere
   Rand_final_data[3,ifile] = dummy[2] ;Save skewness trophosphere
   Rand_final_data[4,ifile] = dummy[3] ;Save kurtosis trophosphere
   Rand_final_data[5,ifile] = SQRT(dummy[1]) ;Save standard deviation trophosphere

   dummy = WHERE(FINITE(tropospheric_vcds[min(mask_lon):max(mask_lon), $
                                          min(mask_lat):max(mask_lat), $
                                          ifile]),count)
   Rand_final_data[0,ifile] = count ; Save number of data points averaged over Aral Sea

   dummy = MOMENT(stratospheric_vcds[min(mask_lon):max(mask_lon), $
                                     min(mask_lat):max(mask_lat), $
                                     ifile], /nan)
   Rand_final_data[6,ifile] = dummy[0] ;Save mean stratosphere
   Rand_final_data[7,ifile] = SQRT(dummy[1]) ;Save standard deviation stratosphere

   ; Plot only once each 5 files
   IF ( (ifile MOD 5.0) EQ 0.0) THEN BEGIN
      ; Plot only region of interest
      mask_lon = WHERE(lons GE lonmin AND lons LE lonmax)
      mask_lat = WHERE(lats GE latmin AND lats LE latmax)
      
      ; Several panels on each page
      MULTIPANEL, cols=3, rows=3
      TVMAP, colcor[min(mask_lon):max(mask_lon),min(mask_lat):max(mask_lat)],    $
             lons[mask_lon], lats[mask_lat], /Sample, /Cbar, /coast, /countries, $
             /continents, title = 'BrO VCDs '+date_str[ifile], mindata = 0
      
      TVMAP, mask_2d[min(mask_lon):max(mask_lon),min(mask_lat):max(mask_lat)],   $
             lons[mask_lon], lats[mask_lat], /Sample, /Cbar, /coast, /countries, $
             /continents, title = 'Background mask'+                             $
             STRING(TOTAL(Background_lat[4,*,ifile],/Nan),FORMAT = '(I8)')
      
      TVMAP, colerr[min(mask_lon):max(mask_lon),min(mask_lat):max(mask_lat)],   $
             lons[mask_lon], lats[mask_lat], /Sample, /Cbar, /coast, /countries, $
             /continents, title = 'BrO VCDs error', mindata = 0
      
      TVMAP, numsam[min(mask_lon):max(mask_lon),min(mask_lat):max(mask_lat)],   $
             lons[mask_lon], lats[mask_lat], /Sample, /Cbar, /coast, /countries, $
             /continents, title = '# of samples'
      
      TVMAP, alt[min(mask_lon):max(mask_lon),min(mask_lat):max(mask_lat)],   $
             lons[mask_lon], lats[mask_lat], /Sample, /Cbar, /coast, /countries, $
             /continents, title = 'Altitude', mindata = 0
      
      TVMAP, strat_back[min(mask_lon):max(mask_lon),min(mask_lat):max(mask_lat)],   $
             lons[mask_lon], lats[mask_lat], /Sample, /Cbar, /coast, /countries, $
             /continents, title = 'Raw stratospheric background', $
             mindata = 0, nan_color=13, $
             maxdata = MAX(colcor[min(mask_lon):max(mask_lon), $
                                  min(mask_lat):max(mask_lat)], /nan)
      
      TVMAP, stratospheric_vcds[min(mask_lon):max(mask_lon),   $
                                min(mask_lat):max(mask_lat),   $
                                ifile],                        $
          lons[mask_lon], lats[mask_lat], /Sample, /Cbar,   $
             /coast, /countries, /continents,                  $
             title = 'Smoothed stratospheric background',      $
             mindata = 0, nan_color=13,                        $
             maxdata = MAX(stratospheric_vcds[min(mask_lon):max(mask_lon), $
                                              min(mask_lat):max(mask_lat)],/nan)
      
      ;; TVMAP, tropospheric_vcds[min(mask_lon):max(mask_lon),  $
      ;;                          min(mask_lat):max(mask_lat),   $
      ;;                          ifile],                        $
      ;;        lons[mask_lon], lats[mask_lat], /Sample, /Cbar, $
      ;;        /coast, /countries, /continents, /noadvance,    $
      ;;        title = 'Tropospheric VCDs',                    $
      ;;        mindata = 0, nan_color=13, $
      ;;        maxdata = MAX(tropospheric_vcds[min(mask_lon):max(mask_lon), $
      ;;                                        min(mask_lat):max(mask_lat), $
      ;;                                        ifile])
      
      TVMAP, tropospheric_vcds[min(mask_lon):max(mask_lon),   $
                               min(mask_lat):max(mask_lat),   $
                               ifile],                        $
             lons[mask_lon], lats[mask_lat], /Sample, /Cbar,  $
             /coast, /countries, /continents,                 $
             title = 'Tropospheric VCDs',                     $
             mindata = 0, nan_color=13,                       $
             maxdata = 2e13
             ;; maxdata = MAX(tropospheric_vcds[min(mask_lon):max(mask_lon), $
             ;;                                 min(mask_lat):max(mask_lat), $
             ;;                                 ifile],/nan)

      mask_lon = WHERE(lons GE 56 AND lons LE 64)
      mask_lat = WHERE(lats GE 41 AND lats LE 49)
      TVMAP, tropospheric_vcds[min(mask_lon):max(mask_lon),    $
                               min(mask_lat):max(mask_lat),    $
                               ifile],                         $
             lons[mask_lon], lats[mask_lat], /Sample, /Cbar,   $
             /coast, /countries, /continents,                  $
             title = 'Tropospheric VCDs',                      $
             mindata = 0, nan_color=13, limit = [41,56,49,64], $
             maxdata = MAX(tropospheric_vcds[min(mask_lon):max(mask_lon), $
                                             min(mask_lat):max(mask_lat), $
                                             ifile],/nan)
      
   ENDIF

; End loop over files
ENDFOR
MULTIPANEL, /OFF
!x.margin = x_margin
!y.margin = y_margin
MULTIPANEL, cols = 1, rows = 3

mask = WHERE(date_int NE 0, count)
cgplot, INDGEN(count)+1, Aral_final_data[0,mask], color = 1, xtitle ='Day #', $
        ytitle = 'Aral # of data points in average', xrange = [0,count+2], xstyle = 1
cgplot, INDGEN(count)+1, Aral_final_data[1,mask], color = 1, xtitle ='Day #',  $
        ytitle = 'BrO [molecules cm!u-2!n]', xrange = [0,count+2], xstyle = 1, $
        err_yhigh=Aral_final_data[5,mask], err_ylow=Aral_final_data[5,mask],   $
        yrange=[0.0e13,5e13], err_color=9, err_width=0,                     $
        title = 'Aral Tropospheric BrO columns'
cgplot, INDGEN(count)+1, Aral_final_data[6,mask], color = 1, xtitle ='Day #',  $
        ytitle = 'BrO [molecules cm!u-2!n]', xrange = [0,count+2], xstyle = 1, $
        err_yhigh=Aral_final_data[7,mask], err_ylow=Aral_final_data[7,mask],   $
        yrange=[1.0e13,5e13], err_color=9, err_width=0,                     $
        title = 'Aral Stratospheric BrO columns'

cgplot, INDGEN(count)+1, Rand_final_data[0,mask], color = 1, xtitle ='Day #', $
        ytitle = 'Pacific # of data points in average', xrange = [0,count+2], xstyle = 1
cgplot, INDGEN(count)+1, Rand_final_data[1,mask], color = 1, xtitle ='Day #',  $
        ytitle = 'BrO [molecules cm!u-2!n]', xrange = [0,count+2], xstyle = 1, $
        err_yhigh=Rand_final_data[5,mask], err_ylow=Rand_final_data[5,mask],   $
        yrange=[0.0e13,5e13], err_color=9, err_width=0,                     $
        title = 'Pacific Tropospheric BrO columns'
cgplot, INDGEN(count)+1, Rand_final_data[6,mask], color = 1, xtitle ='Day #',  $
        ytitle = 'BrO [molecules cm!u-2!n]', xrange = [0,count+2], xstyle = 1, $
        err_yhigh=Rand_final_data[7,mask], err_ylow=Rand_final_data[7,mask],   $
        yrange=[1.0e13,5e13], err_color=9, err_width=0,                     $
        title = 'Pacific Stratospheric BrO columns'

MULTIPANEL, /OFF
MULTIPANEL, /RESET
; Close postcript document
!x.margin = x_margin
!y.margin = y_margin
CLOSE_DEVICE

; Save data to file
SAVE, stratospheric_vcds, tropospheric_vcds, Background_lat, $
      date_int, date_str, $
      filename = 'Aral_strat_trop_vcds_0p2x0p2_SZA70_XT1-23_CF30_global.sav'
SAVE, aral_final_data, rand_final_data, date_int, date_str, $
      filename = 'Aral_time_series_0p2x0p2_SZA70_XT1-23_CF30_global.sav'

END
