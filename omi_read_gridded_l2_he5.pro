; ===================================================
; Example: Access an OMI HCHO data file to read
;          - AMF
;          - Total column amount
;          - Total column corrected
;          - Total column uncertainty
;          - Grid Area
;          - Number of Samples
;          - Main data quality flag
;          - Slant column amount
;          - Slant corrected column amount
;          - Slant colunm uncertainty
;          - Surface albedo
;          - Pixel center latitude
;          - Pixel center longitude
;          - Pixel corner latitudes
;          - Pixel corner longitudes
; ===================================================
;

PRO omi_read_gridded_l2_he5, he5file,                                       $
                             swathname, nlon, nlat, dlon, dlat,             $
                             AMF=amf, COLFIT=colfit, COLCOR=colcor,         $
                             COLERR=colerr, GRIDAR=gridar, NUMSAM = numsam, $
                             QFLAG=qflag, RMS=rms, SLTFIT=sltfit,           $
                             SLTCOR=SLTCOR, SLTERR=slterr, ALBEDO=albedo,   $
                             ALT=alt, CFR=cfr, CTP=ctp, LATS=lats,          $
                             LONS=lons, HELP=help

   ; -------------------------------------
   ; Initialize mandatory output variables
   ; -------------------------------------
   nlon = -1   & nlat = -1
   dlon = -1.0 & dlat = -1.0

   ; -------------------------------------
   ; Print HELP screen and return
   ; -------------------------------------
   IF ( KEYWORD_SET(HELP)) THEN BEGIN
      PRINT, 'Use of this routine:'
      PRINT, ''
      PRINT, '  OMI_READ_GRIDDED_L2_HE5, he5file, swathname, nlon, nlat [, KEY1=key1, KEY2=key2, ...]'
      PRINT, ''
      PRINT, 'Input:  he5file'
      PRINT, 'Output: swathname, nlon, nlat, dlon, dlat [, key1, ...]'
      PRINT, ''
      PRINT, '----------------'
      PRINT, 'The Data Fields:'
      PRINT, '----------------'
      PRINT, 'ALT       Surface altitude'
      PRINT, 'AMF       AMF'
      PRINT, 'CFR       Cloud Fraction'
      PRINT, 'CTP       Cloud Pressure'
      PRINT, 'COLFIT    ColumnFitted'
      PRINT, 'COLCOR    ColumnFittedDestriped'
      PRINT, 'COLERR    ColumnUncertainty'
      PRINT, 'GRIDAR    GridArea'
      PRITN, 'NUMSAM    NumberOfSamples'
      PRINT, 'QFLAG     QualityFlag'
      PRINT, 'RMS       RMS'
      PRINT, 'SLTFIT    SlantFitted'
      PRINT, 'SLTCOR    SlantFittedDestriped'
      PRINT, 'SLTERR    SlantUncertainty'
      PRINT, 'ALBEDO    SurfaceAlbedo'
      PRINT, ''
      PRINT, '-----------------------'
      PRINT, 'The Geolocation Fields:'
      PRINT, '-----------------------'
      PRINT, 'LATS       Latitudes' 
      PRINT, 'LONS       Longitudes'

      RETURN
   ENDIF


   ; ========================================================
   ; Read data fields from the OMI L2 gridded HE5 file
   ; ========================================================

   ; -------------------
   ; Open the HDF5 file.
   ; -------------------
   file_id = H5F_OPEN( he5file ) 
   IF ( file_id LE 0L ) THEN BEGIN
      PRINT, 'ERROR opening file >>'+he5file+'<<'
      RETURN
   ENDIF
   
   ; ----------------------------
   ; Compose the total swath name
   ; ----------------------------
   swathname = ''
   sidx = 0
   swathname = H5G_GET_MEMBER_NAME(file_id, 'HDFEOS/SWATHS', sidx)
   swname    = '/HDFEOS/SWATHS/'+swathname+'/'
   sfname    = 'ScaleFactor'
   mvname    = 'MissingValue'
   gsname    = 'GridSpacing'

   ; =======================
   ; >>>>> Data Fields <<<<<
   ; =======================
   ; -----------------------
   ; AMF
   ; -----------------------
   IF ( ARG_PRESENT(AMF) ) THEN BEGIN
      dfname = 'Data Fields/AMF'
      dsid   = H5D_OPEN(file_id, swname + dfname)
      amf    = H5D_READ(dsid)

      atid = H5A_OPEN_NAME(dsid, sfname)
      sfac = H5A_READ(atid)
      H5A_CLOSE, atid

      atid = H5A_OPEN_NAME(dsid, mvname)
      mval = H5A_READ(atid)
      H5A_CLOSE, atid

      idx = WHERE ( amf GT mval, cnt )
      IF ( cnt GT 0 ) THEN amf(idx)  = amf(idx) * sfac(0)

      nlon = N_ELEMENTS(amf(*,0))
      nlat = N_ELEMENTS(amf(0,*))

      H5D_CLOSE, dsid
   ENDIF

   ; -----------------------
   ; Fitted Columns, Regular
   ; -----------------------
   IF ( ARG_PRESENT(COLFIT) ) THEN BEGIN
      dfname = 'Data Fields/ColumnFitted'
      dsid   = H5D_OPEN(file_id, swname + dfname)
      colfit = H5D_READ(dsid)

      atid = H5A_OPEN_NAME(dsid, sfname)
      sfac = H5A_READ(atid)
      H5A_CLOSE, atid

      atid = H5A_OPEN_NAME(dsid, mvname)
      mval = H5A_READ(atid)
      H5A_CLOSE, atid

      idx = WHERE ( colfit GT mval, cnt )
      IF ( cnt GT 0 ) THEN colfit(idx)  = colfit(idx) * sfac(0)

      nlon = N_ELEMENTS(colfit(*,0))
      nlat = N_ELEMENTS(colfit(0,*))

      H5D_CLOSE, dsid
   ENDIF

   ; -------------------------
   ; Fitted Columns, corrected
   ; -------------------------
   IF ( ARG_PRESENT(COLCOR) ) THEN BEGIN
      dfname = 'Data Fields/ColumnFittedDestriped'
      dsid   = H5D_OPEN(file_id, swname + dfname)
      colcor = H5D_READ(dsid)

      atid = H5A_OPEN_NAME(dsid, sfname)
      sfac = H5A_READ(atid)
      H5A_CLOSE, atid

      atid = H5A_OPEN_NAME(dsid, mvname)
      mval = H5A_READ(atid)
      H5A_CLOSE, atid

      idx = WHERE ( colcor GT mval, cnt )
      IF ( cnt GT 0 ) THEN colcor(idx)  = colcor(idx) * sfac(0)

      nlon = N_ELEMENTS(colcor(*,0))
      nlat = N_ELEMENTS(colcor(0,*))

      H5D_CLOSE, dsid
   ENDIF

   ; --------------------
   ; Column Uncertainties
   ; --------------------
   IF ( ARG_PRESENT(COLERR) ) THEN BEGIN
      dfname = 'Data Fields/ColumnUncertainty'
      dsid   = H5D_OPEN(file_id, swname + dfname)
      colerr = H5D_READ(dsid)

      atid = H5A_OPEN_NAME(dsid, sfname)
      sfac = H5A_READ(atid)
      H5A_CLOSE, atid

      atid = H5A_OPEN_NAME(dsid, mvname)
      mval = H5A_READ(atid)
      H5A_CLOSE, atid

      idx = WHERE ( colerr GT mval, cnt )
      IF ( cnt GT 0 ) THEN colerr(idx)  = colerr(idx) * sfac(0)

      nlon   = N_ELEMENTS(colerr(*,0))
      nlat   = N_ELEMENTS(colerr(0,*))

      H5D_CLOSE, dsid
   ENDIF

   ; -------------------------
   ; Grid Areas
   ; -------------------------
   IF ( ARG_PRESENT(GRIDAR) ) THEN BEGIN
      dfname = 'Data Fields/GridArea'
      dsid   = H5D_OPEN(file_id, swname + dfname)
      gridar = H5D_READ(dsid)

      atid = H5A_OPEN_NAME(dsid, sfname)
      sfac = H5A_READ(atid)
      H5A_CLOSE, atid

      atid = H5A_OPEN_NAME(dsid, mvname)
      mval = H5A_READ(atid)
      H5A_CLOSE, atid

      idx = WHERE ( gridar GT mval, cnt )
      IF ( cnt GT 0 ) THEN gridar(idx)  = gridar(idx) * sfac(0)

      nlon   = N_ELEMENTS(gridar(*,0))
      nlat   = N_ELEMENTS(gridar(0,*))

      H5D_CLOSE, dsid
   ENDIF

   ; -----------------
   ; Number of samples
   ; -----------------
   IF ( ARG_PRESENT(NUMSAM) ) THEN BEGIN
      dfname = 'Data Fields/NumberOfSamples'
      dsid   = H5D_OPEN(file_id, swname + dfname)
      numsam = H5D_READ(dsid)

      atid = H5A_OPEN_NAME(dsid, sfname)
      sfac = H5A_READ(atid)
      H5A_CLOSE, atid

      numsam(idx)  = numsam(idx) * sfac(0)

      nlon   = N_ELEMENTS(numsam(*,0))
      nlat   = N_ELEMENTS(numsam(0,*))

      H5D_CLOSE, dsid
   ENDIF

   ; ------------------------------
   ; Quality Flag
   ; ------------------------------
   IF ( ARG_PRESENT(QFLAG) ) THEN BEGIN
      dfname = 'Data Fields/QualityFlag'
      dsid  = H5D_OPEN(file_id, swname + dfname)
      qflag = H5D_READ(dsid)

      nlon  = N_ELEMENTS(qflag(*,0))
      nlat  = N_ELEMENTS(qflag(0,*))

      H5D_CLOSE, dsid
   ENDIF

   ; -----------------
   ; RMS
   ; -----------------
   IF ( ARG_PRESENT(RMS) ) THEN BEGIN
      dfname = 'Data Fields/RMS'
      dsid   = H5D_OPEN(file_id, swname + dfname)
      rms    = H5D_READ(dsid)

      atid = H5A_OPEN_NAME(dsid, sfname)
      sfac = H5A_READ(atid)
      H5A_CLOSE, atid

      atid = H5A_OPEN_NAME(dsid, mvname)
      mval = H5A_READ(atid)
      H5A_CLOSE, atid

      idx = WHERE ( rms GT mval, cnt )
      IF ( cnt GT 0 ) THEN rms(idx)  = rms(idx) * sfac(0)

      nlon   = N_ELEMENTS(rms(*,0))
      nlat   = N_ELEMENTS(rms(0,*))

      H5D_CLOSE, dsid
   ENDIF

   ; ----------------------
   ; Slant Columns, Regular
   ; ----------------------
   IF ( ARG_PRESENT(SLTFIT) ) THEN BEGIN
      dfname = 'Data Fields/SlantFitted'
      dsid   = H5D_OPEN(file_id, swname + dfname)
      sltfit = H5D_READ(dsid)

      atid = H5A_OPEN_NAME(dsid, sfname)
      sfac = H5A_READ(atid)
      H5A_CLOSE, atid

      atid = H5A_OPEN_NAME(dsid, mvname)
      mval = H5A_READ(atid)
      H5A_CLOSE, atid

      idx = WHERE ( sltfit GT mval, cnt )
      IF ( cnt GT 0 ) THEN sltfit(idx)  = sltfit(idx) * sfac(0)

      nlon = N_ELEMENTS(sltfit(*,0))
      nlat = N_ELEMENTS(sltfit(0,*))

      H5D_CLOSE, dsid
   ENDIF

   ; -------------------------
   ; Slant Columns, corrected
   ; -------------------------
   IF ( ARG_PRESENT(SLTCOR) ) THEN BEGIN
      dfname = 'Data Fields/SlantFittedDestriped'
      dsid   = H5D_OPEN(file_id, swname + dfname)
      sltcor = H5D_READ(dsid)

      atid = H5A_OPEN_NAME(dsid, sfname)
      sfac = H5A_READ(atid)
      H5A_CLOSE, atid

      atid = H5A_OPEN_NAME(dsid, mvname)
      mval = H5A_READ(atid)
      H5A_CLOSE, atid

      idx = WHERE ( sltcor GT mval, cnt )
      IF ( cnt GT 0 ) THEN sltcor(idx)  = sltcor(idx) * sfac(0)

      nlon = N_ELEMENTS(sltcor(*,0))
      nlat = N_ELEMENTS(sltcor(0,*))

      H5D_CLOSE, dsid
   ENDIF

   ; -------------------
   ; Slant Uncertainties
   ; -------------------
   IF ( ARG_PRESENT(COLERR) ) THEN BEGIN
      dfname = 'Data Fields/SlantUncertainty'
      dsid   = H5D_OPEN(file_id, swname + dfname)
      slterr = H5D_READ(dsid)

      atid = H5A_OPEN_NAME(dsid, sfname)
      sfac = H5A_READ(atid)
      H5A_CLOSE, atid

      atid = H5A_OPEN_NAME(dsid, mvname)
      mval = H5A_READ(atid)
      H5A_CLOSE, atid

      idx = WHERE ( slterr GT mval, cnt )
      IF ( cnt GT 0 ) THEN slterr(idx)  = slterr(idx) * sfac(0)

      nlon   = N_ELEMENTS(slterr(*,0))
      nlat   = N_ELEMENTS(slterr(0,*))

      H5D_CLOSE, dsid
   ENDIF

   ; -----------------------
   ; Surface Albedo
   ; -----------------------
   IF ( ARG_PRESENT(ALBEDO) ) THEN BEGIN
      dfname = 'Data Fields/SurfaceAlbedo'
      dsid   = H5D_OPEN(file_id, swname + dfname)
      albedo = H5D_READ(dsid)

      atid = H5A_OPEN_NAME(dsid, sfname)
      sfac = H5A_READ(atid)
      H5A_CLOSE, atid

      atid = H5A_OPEN_NAME(dsid, mvname)
      mval = H5A_READ(atid)
      H5A_CLOSE, atid

      idx = WHERE ( albedo GT mval, cnt )
      IF ( cnt GT 0 ) THEN albedo(idx)  = albedo(idx) * sfac(0)

      nlon = N_ELEMENTS(albedo(*,0))
      nlat = N_ELEMENTS(albedo(0,*))

      H5D_CLOSE, dsid
   ENDIF

   ; ----------------
   ; Surface altitude
   ; ----------------
   IF ( ARG_PRESENT(ALT) ) THEN BEGIN
      dfname = 'Data Fields/SurfaceAltitude'
      dsid   = H5D_OPEN(file_id, swname + dfname)
      alt    = H5D_READ(dsid)

      atid = H5A_OPEN_NAME(dsid, sfname)
      sfac = H5A_READ(atid)
      H5A_CLOSE, atid

      atid = H5A_OPEN_NAME(dsid, mvname)
      mval = H5A_READ(atid)
      H5A_CLOSE, atid

      idx = WHERE ( alt GT mval, cnt )
      IF ( cnt GT 0 ) THEN alt(idx)  = alt(idx) * sfac(0)

      nlon = N_ELEMENTS(alt(*,0))
      nlat = N_ELEMENTS(alt(0,*))

      H5D_CLOSE, dsid
   ENDIF

   ; --------------
   ; Cloud fraction
   ; --------------
   IF ( ARG_PRESENT(CFR) ) THEN BEGIN
      dfname = 'Data Fields/CloudFraction'
      dsid   = H5D_OPEN(file_id, swname + dfname)
      cfr    = H5D_READ(dsid)

      atid = H5A_OPEN_NAME(dsid, sfname)
      sfac = H5A_READ(atid)
      H5A_CLOSE, atid

      atid = H5A_OPEN_NAME(dsid, mvname)
      mval = H5A_READ(atid)
      H5A_CLOSE, atid

      idx = WHERE ( cfr GT mval, cnt )
      IF ( cnt GT 0 ) THEN cfr(idx)  = cfr(idx) * sfac(0)

      nlon = N_ELEMENTS(cfr(*,0))
      nlat = N_ELEMENTS(cfr(0,*))

      H5D_CLOSE, dsid
   ENDIF

   ; --------------
   ; Cloud pressure
   ; --------------
   IF ( ARG_PRESENT(CTP) ) THEN BEGIN
      dfname = 'Data Fields/CloudPressure'
      dsid   = H5D_OPEN(file_id, swname + dfname)
      ctp    = H5D_READ(dsid)

      atid = H5A_OPEN_NAME(dsid, sfname)
      sfac = H5A_READ(atid)
      H5A_CLOSE, atid

      atid = H5A_OPEN_NAME(dsid, mvname)
      mval = H5A_READ(atid)
      H5A_CLOSE, atid

      idx = WHERE ( ctp GT mval, cnt )
      IF ( cnt GT 0 ) THEN ctp(idx)  = ctp(idx) * sfac(0)

      nlon = N_ELEMENTS(ctp(*,0))
      nlat = N_ELEMENTS(ctp(0,*))

      H5D_CLOSE, dsid
   ENDIF

   ; ==============================
   ; >>>>> Geolocation Fields <<<<<
   ; ==============================

   ; ------------
   ; Latitudes
   ; ------------
   IF ( ARG_PRESENT(LATS) ) THEN BEGIN
      gfname = 'Geolocation Fields/Latitudes'
      dsid   = H5D_OPEN(file_id, swname + gfname) 

      lats   = H5D_READ(dsid)
      nlat   = N_ELEMENTS(lats(*))

      atid = H5A_OPEN_NAME(dsid, gsname)
      dlat = H5A_READ(atid)
      H5A_CLOSE, atid

      H5D_CLOSE, dsid
   ENDIF

   ; -------------
   ; Longitudes
   ; -------------
   IF ( ARG_PRESENT(LONS) ) THEN BEGIN
      gfname = 'Geolocation Fields/Longitudes'
      dsid   = H5D_OPEN(file_id, swname + 'Geolocation Fields/Longitudes')

      lons   = H5D_READ(dsid)
      nlon   = N_ELEMENTS(lons(*))

      atid = H5A_OPEN_NAME(dsid, gsname)
      dlon = H5A_READ(atid)
      H5A_CLOSE, atid

      H5D_CLOSE, dsid
   ENDIF

END
