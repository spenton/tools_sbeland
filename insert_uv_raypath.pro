;+
; NAME:   INSERT_UV_RAYPATH
;
; PURPOSE: 
;    This routine populates the SORCE/SimRayPathDegradationParams database table with the 
;    new RayPath degradation factor for the UV diodes obtained when comparing the irradiance 
;    of the SimA and SimB data (with adjust_pd_raypath.pro).
;
; CALLING SEQUENCE:
;    res=insert_raypath(/dbinsert)
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
;   The raypath values are obtained from the code adjust_pd_raypath.pro
;
; REVISION HISTORY:
;   Revision: $Id: insert_uv_raypath.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;-
;*****************************************************************************
PRO update_calibMetadata, instMode, calId, effDate, version=version, conn17=conn17, $
    verbose=verbose, dbinsert=insert

    ; make sure we have an entry in the CalibrationMetadata table for this calibrationSetId
    q1="SELECT * from CalibrationMetadata where calibrationSetId="+strtrim(string(CalId),2)
    if insert then begin
        query_database,/reset
        query_database,q1,data,nrows
        date_str=string(systime(/jul),format='(C(CYI,"-",CMOI2.2,"-",CDI2.2," ",CHI2.2,":",CMI2.2,":",CSI2.2))')
        if nrows eq 0 then begin
            ; insert a new row in the CalibrationMetadata table
            q1="INSERT INTO dbo.CalibrationMetadata(calibrationSetId, calibrationTableName, instrumentModeId, "
            q1+="version, releaseDate, effectiveDate, rationaleForNewVersion, responsibleIndividual, "
            q1+="dataDescription, sourceOfData, applicableReferences, additionalInformation) VALUES("
            q1=q1+strtrim(string(CalId),2)+", 'SimRayPathDegradationParams', "
            q1=q1+strtrim(string(instMode),2)+", "+strtrim(string(version),2)
            q1=q1+", '"+date_str+"', '"+effDate+"'"
            q1=q1+", 'New RayPath for version 20 from updated wavelength knowledge and new KappaDegradation'"
            q1=q1+", 'Stephane Beland'"
            q1=q1+", 'RayPath degradation as (1-A)*exp(K*degCol)+A*exp(K*degCol*B) '"
            q1=q1+", 'Stephane processing and ESR to diode comparison'"
            q1=q1+", 'Generated by straightening the difference in diode and ESR signal time series as a function of wavelength'"
            q1=q1+", 'none')"
            if keyword_set(verbose) then print,q1
            if insert then result=conn17->execute(q1)
        endif else begin
            q1="UPDATE CalibrationMetadata SET releaseDate='"+date_str+"'"
            q1=q1+", effectiveDate='"+effDate+"'"
            q1=q1+", version="+strtrim(string(version),2)
            q1=q1+" WHERE calibrationSetId="+strtrim(string(CalId),2)
            if keyword_set(verbose) then print,q1
            if insert then result=conn17->execute(q1)
        endelse
    endif

END
;-----------------------------------------------------------------------------
PRO my_insertRay, CalId, wave, singlePass, firstSurf, conn17=conn17, verbose=verbose, dbinsert=insert

    ; insert one row at a time
    q1="INSERT INTO dbo.SimRayPathDegradationParams(calibrationSetId, wavelength, singlePassAreaFraction, "+$
       "singlePassAreaFractionUnc, firstSurfaceDegradation, firstSurfaceDegradationUnc) VALUES("+$
       strtrim(string(CalId,format='(I)'),2)+', '
    nwave=n_elements(wave)
    for k=0L,n_elements(wave)-1 do begin
        print,k+1,double(k+1)/double(nwave+1)*100d,format='("   ",i0,"  (",i0," % )",$,%"\r")'
        str = strtrim(string(wave[k]),2)
        str = str+', '+strtrim(string(singlePass[k]),2)
        str = str+', 0.0'
        str = str+', '+strtrim(string(firstSurf[k]),2)
        str = str+', 0.0'
        q2 = q1+str+')'
        if keyword_set(verbose) then print,q2
        if insert then result=conn17->execute(q2)
    endfor
    print,''

END
;-----------------------------------------------------------------------------
FUNCTION INSERT_UV_RAYPATH, dbinsert=dbinsert, verbose=verbose

    ;INSERT INTO dbo.SimRayPathDegradationParams(calibrationSetId, wavelength, singlePassAreaFraction, 
    ;    singlePassAreaFractionUnc, firstSurfaceDegradation, firstSurfaceDegradationUnc) VALUES(0, 0, 0, 0, 0, 0)


    ; for version 20 we settled on a single raypath for the whole mission
    version=6
    newCalId=[1198,1199]
    instrumentModeId=[43,47]
    effectiveDate='2003-01-01 00:00:00.0'

    readcol,'~/SORCE/data/overlap_uv_smooth_pos_160.txt',wave0,singlePass,format='(d,d)'
    firstSurf = singlePass*0d +0.5d

    stmt17 = fjava_get_jdbc_statement(user='sorce', password='sorcedb', server='sorce-db', database="SORCE")
    conn17 = OBJ_NEW("oJava_DbExchange", stmt17)
    if keyword_set(dbinsert) then insert=1 else insert=0

    for i=0, n_elements(instrumentModeId)-1 do begin
        update_calibMetadata, instrumentModeId[i], newCalId[i], effectiveDate, $
            version=version, conn17=conn17, dbinsert=insert, verbose=verbose

        ; first delete all entries in the table for this particular version (version 1 only for now)
        q1 = "DELETE FROM SimRayPathDegradationParams WHERE calibrationSetId="+strtrim(string(newCalId[i]),2)
        if keyword_set(verbose) then print,q1
        if insert then result=conn17->execute(q1)

        my_insertRay, newCalId[i], wave0, singlePass, firstSurf, conn17=conn17, verbose=verbose, dbinsert=insert
    endfor

    return,1





    version=4
    newCalId=[1107,1108]
    ;readcol,'~/SORCE/data/UVAB_raypath_50-178_v16.txt',wave0,singlePass,format='(d,d)'
    ;readcol,'~/SORCE/data/UVAB_raypath_50-453_v16_5.txt',wave0,singlePass,format='(d,d)'
    ;readcol,'~/SORCE/data/UVAB_raypath_50-178_v70.txt',wave0,singlePass,format='(d,d)'
    readcol,'~/SORCE/data/UVAB_raypath_50-178_v93.txt',wave0,singlePass,format='(d,d)'
    singlePass -= 0.030d
    firstSurf = singlePass*0d +0.5d
    effectiveDate='2003-01-01 00:00:00.0'
    all_calId=intarr(n_elements(wave0))+newCalId[0]
    all_wave=wave0
    all_single=singlePass
    all_first=firstSurf

    stmt17 = fjava_get_jdbc_statement(user='sorce', password='sorcedb', server='sorce-db', database="SORCE")
    conn17 = OBJ_NEW("oJava_DbExchange", stmt17)
    if keyword_set(dbinsert) then insert=1 else insert=0

    for i=0, n_elements(instrumentModeId)-1 do begin
        update_calibMetadata, instrumentModeId[i], newCalId[i], effectiveDate, $
            version=version, conn17=conn17, dbinsert=insert, verbose=verbose

        ; first delete all entries in the table for this particular version (version 1 only for now)
        q1 = "DELETE FROM SimRayPathDegradationParams WHERE calibrationSetId="+strtrim(string(newCalId[i]),2)
        if keyword_set(verbose) then print,q1
        if insert then result=conn17->execute(q1)

        my_insertRay, newCalId[i], wave0, singlePass, firstSurf, conn17=conn17, verbose=verbose, dbinsert=insert
    endfor

    ;readcol,'~/SORCE/data/UVAB_raypath_0-453.txt',wave0,singlePass,format='(d,d)'
    ; adjust raypath to match UVA and UVB TSI
    ;readcol,'~/SORCE/data/UVAB_raypath_178-453_v16.txt',wave0,singlePass,format='(d,d)'
    ;readcol,'~/SORCE/data/UVAB_raypath_453-1570_v16.txt',wave0,singlePass,format='(d,d)'
    ;readcol,'~/SORCE/data/UVAB_raypath_178-453_v70.txt',wave0,singlePass,format='(d,d)'
    readcol,'~/SORCE/data/UVAB_raypath_178-453_v93.txt',wave0,singlePass,format='(d,d)'
    singlePass -= 0.030d
    firstSurf = singlePass*0d +0.5d
    newCalId=[1059,1063]
    effectiveDate='2003-07-21 00:00:00.0'
    all_calId=intarr(n_elements(wave0))+newCalId[0]
    all_wave=wave0
    all_single=singlePass
    all_first=firstSurf

    stmt17 = fjava_get_jdbc_statement(user='sorce', password='sorcedb', server='sorce-db', database="SORCE")
    conn17 = OBJ_NEW("oJava_DbExchange", stmt17)
    if keyword_set(dbinsert) then insert=1 else insert=0

    for i=0, n_elements(instrumentModeId)-1 do begin
        update_calibMetadata, instrumentModeId[i], newCalId[i], effectiveDate, $
            version=version, conn17=conn17, dbinsert=insert, verbose=verbose

        ; first delete all entries in the table for this particular version (version 1 only for now)
        q1 = "DELETE FROM SimRayPathDegradationParams WHERE calibrationSetId="+strtrim(string(newCalId[i]),2)
        if keyword_set(verbose) then print,q1
        if insert then result=conn17->execute(q1)

        my_insertRay, newCalId[i], wave0, singlePass, firstSurf, conn17=conn17, verbose=verbose, dbinsert=insert
    endfor

    ;readcol,'~/SORCE/data/UVAB_raypath_453-1570.txt',wave0,singlePass,format='(d,d)'
    ; adjust raypath to match UVA and UVB TSI
    ;readcol,'~/SORCE/data/UVAB_raypath_453-1570_v16.txt',wave0,singlePass,format='(d,d)'
    ;readcol,'~/SORCE/data/UVAB_raypath_453-1570_v70.txt',wave0,singlePass,format='(d,d)'
    readcol,'~/SORCE/data/UVAB_raypath_453-1570_v93.txt',wave0,singlePass,format='(d,d)'
    singlePass -= 0.030d
    firstSurf = singlePass*0d +0.5d
    newCalId=[1060,1064]
    effectiveDate='2004-04-21 00:00:00.0'
    all_calId=[all_calId, intarr(n_elements(wave0))+newCalId[0]]
    all_wave=[all_wave, wave0]
    all_single=[all_single, singlePass]
    all_first=[all_first, firstSurf]

    stmt17 = fjava_get_jdbc_statement(user='sorce', password='sorcedb', server='sorce-db', database="SORCE")
    conn17 = OBJ_NEW("oJava_DbExchange", stmt17)
    if keyword_set(dbinsert) then insert=1 else insert=0

    for i=0, n_elements(instrumentModeId)-1 do begin
        update_calibMetadata, instrumentModeId[i], newCalId[i], effectiveDate, $
            version=version, conn17=conn17, dbinsert=insert, verbose=verbose

        ; first delete all entries in the table for this particular version (version 1 only for now)
        q1 = "DELETE FROM SimRayPathDegradationParams WHERE calibrationSetId="+strtrim(string(newCalId[i]),2)
        if keyword_set(verbose) then print,q1
        if insert then result=conn17->execute(q1)

        my_insertRay, newCalId[i], wave0, singlePass, firstSurf, conn17=conn17, verbose=verbose, dbinsert=insert
    endfor


    ;readcol,'~/SORCE/data/UVAB_raypath_1570-2173.txt',wave0,singlePass,format='(d,d)'
    ; adjust raypath to match UVA and UVB TSI
    ;readcol,'~/SORCE/data/UVAB_raypath_1570-2173_v16.txt',wave0,singlePass,format='(d,d)'
    ;readcol,'~/SORCE/data/UVAB_raypath_1571-2173_v16_5.txt',wave0,singlePass,format='(d,d)'
    ;readcol,'~/SORCE/data/UVAB_raypath_1570-2173_v70.txt',wave0,singlePass,format='(d,d)'
    readcol,'~/SORCE/data/UVAB_raypath_1570-2173_v93.txt',wave0,singlePass,format='(d,d)'
    singlePass -= 0.030d
    firstSurf = singlePass*0d +0.5d
    newCalId=[1061,1065]
    effectiveDate='2007-05-13 00:00:00.0'
    all_calId=[all_calId, intarr(n_elements(wave0))+newCalId[0]]
    all_wave=[all_wave, wave0]
    all_single=[all_single, singlePass]
    all_first=[all_first, firstSurf]

    stmt17 = fjava_get_jdbc_statement(user='sorce', password='sorcedb', server='sorce-db', database="SORCE")
    conn17 = OBJ_NEW("oJava_DbExchange", stmt17)
    if keyword_set(dbinsert) then insert=1 else insert=0

    for i=0, n_elements(instrumentModeId)-1 do begin
        update_calibMetadata, instrumentModeId[i], newCalId[i], effectiveDate, $
            version=version, conn17=conn17, dbinsert=insert, verbose=verbose

        ; first delete all entries in the table for this particular version (version 1 only for now)
        q1 = "DELETE FROM SimRayPathDegradationParams WHERE calibrationSetId="+strtrim(string(newCalId[i]),2)
        if keyword_set(verbose) then print,q1
        if insert then result=conn17->execute(q1)

        my_insertRay, newCalId[i], wave0, singlePass, firstSurf, conn17=conn17, verbose=verbose, dbinsert=insert
    endfor

    ;readcol,'~/SORCE/data/UVAB_raypath_2173-2455.txt',wave0,singlePass,format='(d,d)'
    ;readcol,'~/SORCE/data/UVAB_raypath_2173-3200.txt',wave0,singlePass,format='(d,d)'
    ; adjust raypath to match UVA and UVB TSI
    ;readcol,'~/SORCE/data/UVAB_raypath_2173-2455_v16.txt',wave0,singlePass,format='(d,d)'
    ;readcol,'~/SORCE/data/UVAB_raypath_2173-3200_v16_5.txt',wave0,singlePass,format='(d,d)'
    ;readcol,'~/SORCE/data/UVAB_raypath_2173-3147_v70.txt',wave0,singlePass,format='(d,d)'
    readcol,'~/SORCE/data/UVAB_raypath_2173-3147_v93.txt',wave0,singlePass,format='(d,d)'
    singlePass -= 0.030d
    firstSurf = singlePass*0d +0.5d
    newCalId=[1062,1066]
    effectiveDate='2009-01-05 00:00:00.0'
    all_calId=[all_calId, intarr(n_elements(wave0))+newCalId[0]]
    all_wave=[all_wave, wave0]
    all_single=[all_single, singlePass]
    all_first=[all_first, firstSurf]

    stmt17 = fjava_get_jdbc_statement(user='sorce', password='sorcedb', server='sorce-db', database="SORCE")
    conn17 = OBJ_NEW("oJava_DbExchange", stmt17)
    if keyword_set(dbinsert) then insert=1 else insert=0

    for i=0, n_elements(instrumentModeId)-1 do begin
        update_calibMetadata, instrumentModeId[i], newCalId[i], effectiveDate, $
            version=version, conn17=conn17, dbinsert=insert, verbose=verbose

        ; first delete all entries in the table for this particular version (version 1 only for now)
        q1 = "DELETE FROM SimRayPathDegradationParams WHERE calibrationSetId="+strtrim(string(newCalId[i]),2)
        if keyword_set(verbose) then print,q1
        if insert then result=conn17->execute(q1)

        my_insertRay, newCalId[i], wave0, singlePass, firstSurf, conn17=conn17, verbose=verbose, dbinsert=insert
    endfor

    ;readcol,'~/SORCE/data/UVAB_raypath_2455-2803.txt',wave0,singlePass,format='(d,d)'
    ;readcol,'~/SORCE/data/UVAB_raypath_2173-3200.txt',wave0,singlePass,format='(d,d)'
    ; adjust raypath to match UVA and UVB TSI
    ;readcol,'~/SORCE/data/UVAB_raypath_2455-2803_v16.txt',wave0,singlePass,format='(d,d)'
    ;readcol,'~/SORCE/data/UVAB_raypath_2173-3200_v16_5.txt',wave0,singlePass,format='(d,d)'
    ;readcol,'~/SORCE/data/UVAB_raypath_2173-3147_v70.txt',wave0,singlePass,format='(d,d)'
    readcol,'~/SORCE/data/UVAB_raypath_2173-3147_v93.txt',wave0,singlePass,format='(d,d)'
    singlePass -= 0.030d
    firstSurf = singlePass*0d +0.5d
    newCalId=[1067,1068]
    effectiveDate='2009-10-14 00:00:00.0'
    all_calId=[all_calId, intarr(n_elements(wave0))+newCalId[0]]
    all_wave=[all_wave, wave0]
    all_single=[all_single, singlePass]
    all_first=[all_first, firstSurf]

    stmt17 = fjava_get_jdbc_statement(user='sorce', password='sorcedb', server='sorce-db', database="SORCE")
    conn17 = OBJ_NEW("oJava_DbExchange", stmt17)
    if keyword_set(dbinsert) then insert=1 else insert=0

    for i=0, n_elements(instrumentModeId)-1 do begin
        update_calibMetadata, instrumentModeId[i], newCalId[i], effectiveDate, $
            version=version, conn17=conn17, dbinsert=insert, verbose=verbose

        ; first delete all entries in the table for this particular version (version 1 only for now)
        q1 = "DELETE FROM SimRayPathDegradationParams WHERE calibrationSetId="+strtrim(string(newCalId[i]),2)
        if keyword_set(verbose) then print,q1
        if insert then result=conn17->execute(q1)

        my_insertRay, newCalId[i], wave0, singlePass, firstSurf, conn17=conn17, verbose=verbose, dbinsert=insert
    endfor

    ;readcol,'~/SORCE/data/UVAB_raypath_2803-3147.txt',wave0,singlePass,format='(d,d)'
    ;readcol,'~/SORCE/data/UVAB_raypath_2173-3200.txt',wave0,singlePass,format='(d,d)'
    ; adjust raypath to match UVA and UVB TSI
    ;readcol,'~/SORCE/data/UVAB_raypath_2803-2894_v16.txt',wave0,singlePass,format='(d,d)'
    ;readcol,'~/SORCE/data/UVAB_raypath_2173-3200_v16_5.txt',wave0,singlePass,format='(d,d)'
    ;readcol,'~/SORCE/data/UVAB_raypath_2173-3147_v70.txt',wave0,singlePass,format='(d,d)'
    readcol,'~/SORCE/data/UVAB_raypath_2173-3147_v93.txt',wave0,singlePass,format='(d,d)'
    singlePass -= 0.030d
    firstSurf = singlePass*0d +0.5d
    newCalId=[1069,1070]
    effectiveDate='2010-09-27 00:00:00.0'
    all_calId=[all_calId, intarr(n_elements(wave0))+newCalId[0]]
    all_wave=[all_wave, wave0]
    all_single=[all_single, singlePass]
    all_first=[all_first, firstSurf]

    stmt17 = fjava_get_jdbc_statement(user='sorce', password='sorcedb', server='sorce-db', database="SORCE")
    conn17 = OBJ_NEW("oJava_DbExchange", stmt17)
    if keyword_set(dbinsert) then insert=1 else insert=0

    for i=0, n_elements(instrumentModeId)-1 do begin
        update_calibMetadata, instrumentModeId[i], newCalId[i], effectiveDate, $
            version=version, conn17=conn17, dbinsert=insert, verbose=verbose

        ; first delete all entries in the table for this particular version (version 1 only for now)
        q1 = "DELETE FROM SimRayPathDegradationParams WHERE calibrationSetId="+strtrim(string(newCalId[i]),2)
        if keyword_set(verbose) then print,q1
        if insert then result=conn17->execute(q1)

        my_insertRay, newCalId[i], wave0, singlePass, firstSurf, conn17=conn17, verbose=verbose, dbinsert=insert
    endfor



    ;readcol,'~/SORCE/data/UVAB_raypath_2894-3034_v16.txt',wave0,singlePass,format='(d,d)'
    ; the previous raypath is clearly wrong - use the next one instead
    ;readcol,'~/SORCE/data/UVAB_raypath_3034-3147_v16.txt',wave0,singlePass,format='(d,d)'
    ;readcol,'~/SORCE/data/UVAB_raypath_2173-3200_v16_5.txt',wave0,singlePass,format='(d,d)'
    ;readcol,'~/SORCE/data/UVAB_raypath_2173-3147_v70.txt',wave0,singlePass,format='(d,d)'
    readcol,'~/SORCE/data/UVAB_raypath_2173-3200_v93.txt',wave0,singlePass,format='(d,d)'
    singlePass -= 0.030d
    firstSurf = singlePass*0d +0.5d
    newCalId=[1109,1110]
    effectiveDate='2010-12-27 00:00:00.0'
    all_calId=[all_calId, intarr(n_elements(wave0))+newCalId[0]]
    all_wave=[all_wave, wave0]
    all_single=[all_single, singlePass]
    all_first=[all_first, firstSurf]

    stmt17 = fjava_get_jdbc_statement(user='sorce', password='sorcedb', server='sorce-db', database="SORCE")
    conn17 = OBJ_NEW("oJava_DbExchange", stmt17)
    if keyword_set(dbinsert) then insert=1 else insert=0

    for i=0, n_elements(instrumentModeId)-1 do begin
        update_calibMetadata, instrumentModeId[i], newCalId[i], effectiveDate, $
            version=version, conn17=conn17, dbinsert=insert, verbose=verbose

        ; first delete all entries in the table for this particular version (version 1 only for now)
        q1 = "DELETE FROM SimRayPathDegradationParams WHERE calibrationSetId="+strtrim(string(newCalId[i]),2)
        if keyword_set(verbose) then print,q1
        if insert then result=conn17->execute(q1)

        my_insertRay, newCalId[i], wave0, singlePass, firstSurf, conn17=conn17, verbose=verbose, dbinsert=insert
    endfor



    ;readcol,'~/SORCE/data/UVAB_raypath_3034-3147_v16.txt',wave0,singlePass,format='(d,d)'
    ;readcol,'~/SORCE/data/UVAB_raypath_2173-3200_v16_5.txt',wave0,singlePass,format='(d,d)'
    ;readcol,'~/SORCE/data/UVAB_raypath_2173-3147_v70.txt',wave0,singlePass,format='(d,d)'
    readcol,'~/SORCE/data/UVAB_raypath_2173-3200_v93.txt',wave0,singlePass,format='(d,d)'
    singlePass -= 0.030d
    firstSurf = singlePass*0d +0.5d
    newCalId=[1111,1112]
    effectiveDate='2011-05-16 00:00:00.0'
    all_calId=[all_calId, intarr(n_elements(wave0))+newCalId[0]]
    all_wave=[all_wave, wave0]
    all_single=[all_single, singlePass]
    all_first=[all_first, firstSurf]

    stmt17 = fjava_get_jdbc_statement(user='sorce', password='sorcedb', server='sorce-db', database="SORCE")
    conn17 = OBJ_NEW("oJava_DbExchange", stmt17)
    if keyword_set(dbinsert) then insert=1 else insert=0

    for i=0, n_elements(instrumentModeId)-1 do begin
        update_calibMetadata, instrumentModeId[i], newCalId[i], effectiveDate, $
            version=version, conn17=conn17, dbinsert=insert, verbose=verbose

        ; first delete all entries in the table for this particular version (version 1 only for now)
        q1 = "DELETE FROM SimRayPathDegradationParams WHERE calibrationSetId="+strtrim(string(newCalId[i]),2)
        if keyword_set(verbose) then print,q1
        if insert then result=conn17->execute(q1)

        my_insertRay, newCalId[i], wave0, singlePass, firstSurf, conn17=conn17, verbose=verbose, dbinsert=insert
    endfor




    ;readcol,'~/SORCE/data/UVAB_raypath_3147-3568.txt',wave0,singlePass,format='(d,d)'
    ;readcol,'~/SORCE/data/UVAB_raypath_2173-3200.txt',wave0,singlePass,format='(d,d)'
    ; adjust raypath to match UVA and UVB TSI
    ;readcol,'~/SORCE/data/UVAB_raypath_3147-3568_v16.txt',wave0,singlePass,format='(d,d)'
    ;readcol,'~/SORCE/data/UVAB_raypath_2173-3200_v16_5.txt',wave0,singlePass,format='(d,d)'
    ;readcol,'~/SORCE/data/UVAB_raypath_3147-3820_v70.txt',wave0,singlePass,format='(d,d)'
    readcol,'~/SORCE/data/UVAB_raypath_2173-3200_v93.txt',wave0,singlePass,format='(d,d)'
    singlePass -= 0.030d
    firstSurf = singlePass*0d +0.5d
    newCalId=[1071,1072]
    effectiveDate='2011-09-06 00:00:00.0'
    all_calId=[all_calId, intarr(n_elements(wave0))+newCalId[0]]
    all_wave=[all_wave, wave0]
    all_single=[all_single, singlePass]
    all_first=[all_first, firstSurf]

    stmt17 = fjava_get_jdbc_statement(user='sorce', password='sorcedb', server='sorce-db', database="SORCE")
    conn17 = OBJ_NEW("oJava_DbExchange", stmt17)
    if keyword_set(dbinsert) then insert=1 else insert=0

    for i=0, n_elements(instrumentModeId)-1 do begin
        update_calibMetadata, instrumentModeId[i], newCalId[i], effectiveDate, $
            version=version, conn17=conn17, dbinsert=insert, verbose=verbose

        ; first delete all entries in the table for this particular version (version 1 only for now)
        q1 = "DELETE FROM SimRayPathDegradationParams WHERE calibrationSetId="+strtrim(string(newCalId[i]),2)
        if keyword_set(verbose) then print,q1
        if insert then result=conn17->execute(q1)

        my_insertRay, newCalId[i], wave0, singlePass, firstSurf, conn17=conn17, verbose=verbose, dbinsert=insert
    endfor


    return, {wavelength:all_wave, singlePassAreaFraction:all_single, firstSurfaceDegradation:all_first, $
        calibrationSetId:all_calId}

END

