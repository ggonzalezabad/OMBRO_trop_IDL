;restore, 'Aral_strat_trop_vcds_0p2x0p2_SZA70_XT1-23_CF30.sav.sav'
lonmin = 57 & lonmax = 73
latmin = 37 & latmax = 53

mask_lon = WHERE(findgen(1800)*0.2-180.0 GT lonmin AND findgen(1800)*0.2-180.0 LT lonmax)
mask_lat = WHERE(findgen( 900)*0.2- 90.0 GT latmin AND findgen( 900)*0.2- 90.0 LT latmax)

; Convert date_int to Julian dates
ndates   = n_elements(date_int)
date_jul = LONARR(ndates) & year = INTARR(ndates)
month    = INTARR(ndates) & day  = INTARR(ndates)

FOR idate = 0, ndates-1 do begin
   year[idate]     = FIX(STRMID(STRING(date_int[idate],FORMAT='(I8)'),0,4))
   month[idate]    = FIX(STRMID(STRING(date_int[idate],FORMAT='(I8)'),4,2))
   day[idate]      = FIX(STRMID(STRING(date_int[idate],FORMAT='(I8)'),6,2))
   date_jul[idate] = JULDAY(month[idate],day[idate],year[idate])
ENDFOR

; Create array to save mean and stddev
data = FLTARR(4,n_elements(STRATOSPHERIC_VCDS[0,0,*]))

FOR iday = 0, n_elements(STRATOSPHERIC_VCDS[0,0,*])-1 DO BEGIN

   data[0,iday] = MEAN(STRATOSPHERIC_VCDS[min(mask_lon):max(mask_lon), $
                                          min(mask_lat):max(mask_lat),iday], /nan)
   data[1,iday] = STDDEV(STRATOSPHERIC_VCDS[min(mask_lon):max(mask_lon), $
                                            min(mask_lat):max(mask_lat),iday], /nan)

   data[2,iday] = MEAN(TROPOSPHERIC_VCDS[min(mask_lon):max(mask_lon), $
                                         min(mask_lat):max(mask_lat),iday], /nan)
   data[3,iday] = STDDEV(TROPOSPHERIC_VCDS[min(mask_lon):max(mask_lon), $
                                           min(mask_lat):max(mask_lat),iday], /nan)
ENDFOR
data = data/1e13

; Don't plot December, January and February data
mask = WHERE(MONTH NE 12 AND MONTH NE  1 AND MONTH NE 2 AND $
             MONTH NE 3)

OPEN_DEVICE, /PS, /Color, filename ='Other_XT1-23.ps'
!P.Multi = [0,1,3]

cgplot, date_jul[mask], data[2,mask], color = 1,                  $
        ytitle = 'BrO [molecules cm!u-2!n]', xstyle = 1,                     $
        err_yhigh=data[3,mask], err_ylow=data[3,mask], $
        yrange=[min(data[2,mask]-data[3,mask]),        $
                max(data[2,mask]+data[3,mask])],       $
        err_color=9, err_width=0, title = 'Tropospheric BrO columns',        $
        xtickunits = 'Years', CHARSIZE = 1.7

cgplot, date_jul[mask], data[0,mask], color = 1,                  $
        ytitle = 'BrO [molecules cm!u-2!n]', xstyle = 1,                     $
        err_yhigh=data[1,mask], err_ylow=data[1,mask], $
        yrange=[min(data[0,mask]-data[1,mask]),        $
                max(data[0,mask]+data[1,mask])],       $
        err_color=9, err_width=0, title = 'Stratospheric BrO columns',        $
        xtickunits = 'Years', CHARSIZE = 1.7

cgplot, date_jul[mask], data[0,mask]+data[2,mask], color = 1,                 $
        ytitle = 'BrO [molecules cm!u-2!n]', xstyle = 1,                      $
        err_yhigh=data[1,mask], err_ylow=data[1,mask], $
        yrange=[min(data[0,mask]-data[1,mask]),        $
                max(data[0,mask]+data[1,mask])],       $
        err_color=9, err_width=0, title = 'Total BrO columns',        $
        xtickunits = 'Years', CHARSIZE = 1.7

!P.Multi = 0
CLOSE_DEVICE
END
