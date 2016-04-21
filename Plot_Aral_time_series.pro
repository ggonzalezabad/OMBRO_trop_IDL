restore, 'Aral_time_series_0p2x0p2_SZA70_XT1-23_CF30_global.sav'
restore, 'aral_avg_reflectance_05.sav'

; Divide by 1e13
aral_final_data[1:7,*] = aral_final_data[1:7,*]/1e13

; Convert date_int to Julian dates (OMBRO)
ndates   = n_elements(date_int)
date_jul = LONARR(ndates) & year = INTARR(ndates)
month    = INTARR(ndates) & day  = INTARR(ndates)

FOR idate = 0, ndates-1 do begin
   year[idate]     = FIX(STRMID(STRING(date_int[idate],FORMAT='(I8)'),0,4))
   month[idate]    = FIX(STRMID(STRING(date_int[idate],FORMAT='(I8)'),4,2))
   day[idate]      = FIX(STRMID(STRING(date_int[idate],FORMAT='(I8)'),6,2))
   date_jul[idate] = JULDAY(month[idate],day[idate],year[idate])
ENDFOR

; Convert dayofyear (MODIS REFLECTANCE)
dayofyear      = dayofyear + 8
modis_date_jul = LONARR(n_elements(dayofyear))

FOR idate = 0, n_elements(dayofyear)-1 DO BEGIN

   dummy    = FLTARR(5)-1
   dummy[0] = 2005 & dummy[1] = dayofyear[idate]
   modis_date_jul[idate] = DATE_CONV(dummy,'J')
   
ENDFOR

; Don't plot December, January and February data
mask = WHERE(MONTH NE 12 AND MONTH NE  1 AND MONTH NE 2 AND $
             MONTH NE 3  AND YEAR  EQ 2005)
mask = WHERE(YEAR GE 2005 AND MONTH NE 12 AND MONTH NE 1 AND MONTH NE 2)
mask_modis = WHERE(modis_date_jul GT 0)

OPEN_DEVICE, /PS, /Color, $
             filename='Aral_time_series_0p2x0p2_SZA70_XT1-23_CF30.ps'
!P.Multi = [0,1,3]
cgplot, date_jul[mask], Aral_final_data[0,mask], color = 1, $
        ytitle = '# of data points in average', xstyle = 1, $
        xtickunits = 'Months', CHARSIZE = 1.7

cgplot, date_jul[mask], Aral_final_data[1,mask], color = 1,                  $
        ytitle = 'BrO [molecules cm!u-2!n]', xstyle = 1,                     $
        err_yhigh=Aral_final_data[5,mask], err_ylow=Aral_final_data[5,mask], $
        yrange=[min(aral_final_data[1,mask]-aral_final_data[5,mask]),        $
                max(aral_final_data[1,mask]+aral_final_data[5,mask])],       $
        err_color=9, err_width=0, title = 'Tropospheric BrO columns',        $
        xtickunits = 'Months', CHARSIZE = 1.7

cgaxis, YAxis = 1, yrange = [0.05,0.15], color = 1, /save
cgplot, modis_date_jul[mask_modis], aral_avg[mask_modis], color = 3, /Overplot

cgplot, date_jul[mask], Aral_final_data[6,mask], $
        color = 1, xtitle ='Time',    $
        ytitle = 'BrO [molecules cm!u-2!n]', xstyle = 1,                     $
        err_yhigh=Aral_final_data[7,mask], err_ylow=Aral_final_data[7,mask], $
        yrange=[min(aral_final_data[6,mask]-aral_final_data[7,mask]),        $
                max(aral_final_data[6,mask]+aral_final_data[7,mask])],       $
        err_color=9, err_width=0, title = 'Stratospheric BrO columns',  $
        xtickunits = 'Months', CHARSIZE = 1.7

cgaxis, YAxis = 1, yrange = [0.05,0.15], color = 1, /save
cgplot, modis_date_jul[mask_modis], aral_avg[mask_modis], color = 3, /Overplot

!P.Multi = 0
CLOSE_DEVICE
END
