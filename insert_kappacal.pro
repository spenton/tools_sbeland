;+
; NAME:   INSERT_KAPPACAL
;
; PURPOSE: 
;    This routine populates the SORCE/SimPrismDegKappaCal database table with the 
;    new Kappa values obtained from comparing the degradation of ESRA and B and the 
;    UVA and UVB.
;
; CALLING SEQUENCE:
;    insert_kappacal, 54, /dbinsert
;
; OPTIONAL INPUT PARAMETERS:
;    none
;
; OPTIONAL INPUT KEYWORDS:
;    dbinsert - if set will insert the re-formatted data in the database
;
; OUTPUT PARAMETERS:
;    none
;
; OPTIONAL OUTPUT PARAMETERS:
;    NONE
;
; EXAMPLE:  
;
; COMMON BLOCKS:
;    NONE
;
; NOTES:
;    The program expects as input parameter the result from get_esrkappa and get_uvkappa
;    if the flag USEUV is set.  This will combine the results from the UV and ESR data.
;    For the ESR, an exponential fit is used for wavelengths graeter than 400nm. For 
;    wavelengths less than 400, the results from a bspline fit with both the ESR and UV 
;    data is used.
;
; REVISION HISTORY:
;   Revision: $Id: insert_kappacal.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;-
;*****************************************************************************
FUNCTION INSERT_KAPPACAL, instrumentModeId, dbinsert=dbinsert, verbose=verbose

    ;INSERT INTO dbo.SimPrismDegKappaCal(calibrationSetId, x, y, yUncertainty) VALUES(0, 0, 0, 0)

    modes = [31,32,41,43,44,45,47,48,54, 55]
    pos = where(modes eq instrumentModeId, count)
    if count eq 0 then begin
        print,'Error: invalid InstrumentModeId ',instrumentModeId,' [',modes,']'
        return,-1
    endif

    ;version=19
    version=20
    if instrumentModeId eq 31 or instrumentModeId eq 32 or instrumentModeId eq 54 or instrumentModeId eq 55 then begin
        ;if instrumentModeId eq 31 then newCalId=1051 
        ;if instrumentModeId eq 32 then newCalId=1052
        ;if instrumentModeId eq 54 then newCalId=1031
        ;if instrumentModeId eq 55 then newCalId=1032
        if instrumentModeId eq 31 then newCalId=1165 
        if instrumentModeId eq 32 then newCalId=1166
        wave= dindgen(2800)+200d
        ; exponential coefficients for F-Function from ESR data
        ; coeffs=[0.0095408430d, -0.0070394634d, 1.7073941d-05]
        ; latest fit below returns "better" match to data points from 310nm and larger (2013/04/19)
        ;coeffs=[0.011592055d, -0.0075241204d, 1.8037498d-05]

        ; the following coefficients come from iterating the Kappa and F-function over 20 times
        ; converging on new values for Kappa and F-function
        ;coeffs=[0.023980000d,   -0.0070074030d,   3.7123102d-05]

        ; read the actual data and do a bspline to smooth out the curve
        ;readcol,'~/SORCE/data/ESRAB_degradation_453-1570_kappa_nooutliers.txt',wave,kvalue,format='(d,d)'
        ;readcol,'ESRAB_degradation_453-1570_kappa_19.txt',wave,kvalue,format='(d,d)'
        ;readcol,'~/SORCE/data/ESRAB_degradation_453-1570_kappa_19_with_neg.txt',wave,kvalue,format='(d,d)'

        ; LATEST fit (May 10, 2013) with only ESR and FFUNCT=[1.4907, -6.962e-4, 1.6403e-7] and ESR raypath from raytrace
        ;readcol,'~/SORCE/data/ESRAB_degradation_453-1570_kappa_ffunct_raypath.txt',wave,kvalue,format='(d,d)'
        ;coeffs=[0.010193445d, -0.0070877304d, 1.6031084d-05]
        ; new coeffs to match IR end better (forcing last coeff to 1d-5 in get_esrkappa.pro)
        ;coeffs=[0.010014673d, -0.0070126508d, 1.0000000d-05]

        ; new coeffs by letting all 3 coeffs float and then moving the fit down by
        ; simply changing the 3rd coeff to 1d-5 (since we know that value to be good in the IR)
        ;coeffs=[0.009843545d, -0.0070887388d, 1.0000000d-05]

        ; revert to best fit ESR Kappa without adjusting last parameter and letting match data (Version 68)
        ; we'll rely on the raypath to correct for differences
        ;coeffs=[0.010507225d, -0.0071907915d, 1.8650240d-05]

        ; this new kappa was determined with a FFUNCT += SolarExpB/SolarExpA (version 69)
        ;coeffs=[0.0085282734d,   -0.0071666449d,   1.5236341d-05]

        ; this new Kappa was obtained after scaling the solarExp55 by (1 - sum55/sum54) (Version 70)
        ; coeffs=[0.0096235817d,   -0.0071923836d,   1.7092305d-05]

        ; this new kappa was determined with a FFUNCT *= (1 - SolarExpB/SolarExpA) (version 71)
        ;coeffs=[0.013541562d,   -0.0072183885d,   2.3887407d-05]

        ;kappa = coeffs[0] * exp(coeffs[1] * wave) + coeffs[2]

        ; this new kappa is obtained from version 1003 SimUncorrectedIrradiance using the raytrace raypath
        ;readcol,'~/SORCE/data/ESRAB_453-1570_kappa_raypath_raytrace_1003.txt',wave,kappa,format='(d,d)'
        ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2002.txt',wave,kappa,format='(d,d)'
        ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2011_mod_407nm.txt',wave,kappa,format='(d,d)'
        readcol,'~/SORCE/data/ESR_kappa_453_1570_2011_uv_aligned.txt',wave,kappa,format='(d,d)'
        version=21
        if instrumentModeId eq 31 then newCalId=2009
        if instrumentModeId eq 32 then newCalId=2010
        ww=[200d,220d,250d,260d,wave,3000d]
        rr=interpol(kappa,wave,ww)
        kappa=rr
        wave=ww
        wavelog10 = alog10(wave)

    endif else if instrumentModeId eq 41 or instrumentModeId eq 45 then begin
        ;if instrumentModeId eq 41 then newCalId=1053 else newCalId=1056
        if instrumentModeId eq 41 then newCalId=1167 else newCalId=1168
        wave= dindgen(1540)/2d +300d
        wavelog10 = alog10(wave)
        ; for now use the fit from ESR 
        ;coeffs=[0.010193445d, -0.0070877304d, 1.6031084d-05]
        ;kappa = coeffs[0] * exp(coeffs[1] * wave) + coeffs[2]

        ; use the 4th order polynomial fit found by running get_pdkappa with a=0
        ;coeffs=[0.0045833331d, -2.1434220d-05,  3.8946063d-08, -3.2086643d-11,  1.0050902d-14]

        ; use the 4th order polynomial fit found by running get_pdkappa with afact=0.4
        ;coeffs=[0.0068613348d, -3.2330968d-05,  5.9185976d-08, -4.9130765d-11,  1.5505889d-14]
        ;kappa=poly(wave,coeffs)

        ;readcol,'~/SORCE/data/VISAB_kappa_453_1570_ffunct.txt',ww,kvalue,format='(d,d)'
        ;kappa=interpol(kvalue,ww,wave)

        ; as a test, try out the Kappa from the ESR over the VIS wavelength range (we'll adjust the raypath accordingly)
        ;coeffs=[0.0098435447, -0.0070887388d, 1.0000000d-05]

        ; this new Kappa was obtained after scaling the solarExp55 by (1 - sum55/sum54) (Version 70)
        ;coeffs=[0.0096235817d,   -0.0071923836d,   1.7092305d-05]
        ;kappa = coeffs[0] * exp(coeffs[1] * wave) + coeffs[2]

        ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2002.txt',wave0,kappa0,format='(d,d)'
        readcol,'~/SORCE/data/ESR_kappa_453_1570_2011_mod_407nm.txt',wave0,kappa0,format='(d,d)'
        kappa = interpol(kappa0,wave0, wave, /spline)

    endif else if instrumentModeId eq 43 or instrumentModeId eq 47 then begin
        if instrumentModeId eq 43 then newCalId=2074 else newCalId=2075
        ;readcol,'~/SORCE/data/UVAB_degradation_453-1570_kappa_ffunct_noraypath.txt',wave,kappa,format='(d,d)'
        ;readcol,'~/SORCE/data/UVAB_degradation_453-1570_kappa_ffunct_noraypath_filtered.txt',wave,kappa,format='(d,d)'
        ; the latest kappa with afact=0 and which was used to align the kappa for VIS with afact=0.4 20130523 SBeland
        ;readcol,'~/SORCE/data/UVAB_degradation_453-1570_kappa_ffunct_noraypath_cleaned.txt',wave,kappa,format='(d,d)'
        ;readcol,'~/SORCE/data/UVAB_degradation_453-1570_kappa_ffunct_noraypath.txt',wave,kappa,format='(d,d)'
        ;readcol,'~/SORCE/data/UVAB_degradation_453-1570_kappa_ffunct_V70.txt',wave,kappa,format='(d,d)'
        ;readcol,'~/SORCE/data/UVAB_degradation_453-1570_kappa_ffunct_V70.txt',wave,kappa,format='(d,d)'
        ;readcol,'~/SORCE/data/UVAB_453-1570_kappa_2011_160_s5455mod2.txt',wave,kappa,format='(d,d)'
        restore,'~/SORCE/data/UVAB_453-1570_kappa_2011_160_s5455mod4.sav'
        wave=kappa_uv_mod4.x
        kappa=kappa_uv_mod4.y
        version=23
        wavelog10 = alog10(wave)

    endif else if instrumentModeId eq 44 or instrumentModeId eq 48 then begin
        ;if instrumentModeId eq 44 then newCalId=1055 else newCalId=1058
        if instrumentModeId eq 44 then newCalId=1169 else newCalId=1170
        wave= dindgen(900)+800d
        wavelog10 = alog10(wave)
        ; for now use the fit from ESR 
        ;coeffs=[0.010193445d, -0.0070877304d, 1.6031084d-05]
        ;coeffs=[0.0098435447, -0.0070887388d, 1.0000000d-05]

        ; this new Kappa was obtained after scaling the solarExp55 by (1 - sum55/sum54) (Version 70)
        ;coeffs=[0.0096235817d,   -0.0071923836d,   1.7092305d-05]
        ;kappa = coeffs[0] * exp(coeffs[1] * wave) + coeffs[2]

        ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2002.txt',wave0,kappa0,format='(d,d)'
        readcol,'~/SORCE/data/ESR_kappa_453_1570_2011_uv_aligned.txt',wave0,kappa0,format='(d,d)'
        kappa = interpol(kappa0,wave0, wave, /spline)
    endif


    user='sorce'
    password='sorcedb'
    server='sorce-db'
    database='SORCE'
    if keyword_set(dbinsert) then insert=1 else insert=0

    ; make sure we have an entry in the CalibrationMetadata table for this calibrationSetId
    q1="SELECT * from CalibrationMetadata where calibrationSetId="+strtrim(string(newCalId),2)
    if insert then begin
        query_database,/reset
        date_str=string(systime(/jul),format='(C(CYI,"-",CMOI2.2,"-",CDI2.2," ",CHI2.2,":",CMI2.2,":",CSI2.2))')
        query_database,q1,data,nrows,user=user,password=password,server=server
        if nrows eq 0 then begin
            ; insert a new row in the CalibrationMetadata table
            q1="INSERT INTO dbo.CalibrationMetadata(calibrationSetId, calibrationTableName, instrumentModeId, "
            q1+="version, releaseDate, effectiveDate, rationaleForNewVersion, responsibleIndividual, "
            q1+="dataDescription, sourceOfData, applicableReferences, additionalInformation) VALUES("
            q1=q1+strtrim(string(newCalId),2)+", 'SimPrismDegKappaCal', "
            q1=q1+strtrim(string(instrumentModeId),2)+", "+strtrim(string(version),2)
            q1=q1+", '"+date_str+"', '2003-01-01 00:00:00.0'"
            q1=q1+", 'New Kappa for version "+strtrim(string(version),2)+" from updated wavelength knowledge for each separate instrumentModeId instead of SIMA,B'"
            q1=q1+", 'Stephane Beland'"
            q1=q1+", 'Prism degradation kappa as a function of wavelength'"
            q1=q1+", 'Stephane processing and A-B comparison'"
            if keyword_set(esrexp) then begin
                q1=q1+", 'Generated from bspline fit'"
            endif else begin
                q1=q1+", 'Generated from actual measurements with smoothing and bspline fit'"
            endelse
            q1=q1+", 'none')"
            if keyword_set(verbose) then print,q1
            if insert then query_database,q1,result
        endif else begin
            ; update the releaseDate only
            q1="UPDATE CalibrationMetadata SET releaseDate='"+date_str+"'"
            q1=q1+" WHERE calibrationSetId="+strtrim(string(newCalId),2)
            if keyword_set(verbose) then print,q1
            if insert then query_database,q1,result
        endelse
    endif

    ; first delete all entries in the table for this particular version (version 1 only for now)
    q1 = "DELETE FROM SimPrismDegKappaCal WHERE calibrationSetId="+strtrim(string(newcalId),2)
    if keyword_set(verbose) then print,q1
    if insert then query_database,q1,result

    ; insert one row at a time
    q1="INSERT INTO dbo.SimPrismDegKappaCal(calibrationSetId, x, y, yUncertainty) VALUES("
    nwave=n_elements(wave)
    for i=0L,n_elements(wave)-1 do begin
        print,i+1,double(i+1)/double(nwave+1)*100d,format='("   ",i0,"  (",i0," % )",$,%"\r")'
        str = strtrim(string(newCalId,format='(I)'),2)
        str = str+', '+strtrim(string(wavelog10[i]),2)
        str = str+', '+strtrim(string(kappa[i]),2)
        str = str+', 0.0'
        q2 = q1+str+')'
        if keyword_set(verbose) then print,q2
        if insert then query_database,q2,result
    endfor
    print,' '

    return, {wavelength:wave, kappa:kappa}

END

