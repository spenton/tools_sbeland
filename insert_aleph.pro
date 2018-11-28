;+
; NAME:   INSERT_ALEPH
;
; PURPOSE: 
;    This routine populates the SORCE/SimAlephFactorCal database table with the 
;    new aleph factor obtained with get_aleph.pro by comparing data from version 17 on
;    our nominal day (453.67) for SimA.  The spectrum from version 20 form both SimA
;    and SimB are compared with the same SimA reference spectra.
;
; CALLING SEQUENCE:
;    insert_aleph, 41, /dbinsert
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
;  
; REVISION HISTORY:
;   Revision: $Id: insert_aleph.pro,v 1.2 2018/11/26 16:16:02 spenton Exp $
;-
;*****************************************************************************
PRO update_calibMetadata_aleph, instMode, calId, effDate, version=version, conn17=conn17, $
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
            q1=q1+strtrim(string(CalId),2)+", 'SimAlephFactorCal', "
            q1=q1+strtrim(string(instMode),2)+", "+strtrim(string(version),2)
            q1=q1+", '"+date_str+"', '"+effDate+"'"
            q1=q1+", 'Aleph factor for version 21 with latest calibrated irradiance and modified solar exposure record'"
            q1=q1+", 'Stephane Beland'"
            q1=q1+", 'Stores the Aleph factor for each detector as a function of wavelength'"
            q1=q1+", 'IDL text file'"
            q1=q1+", 'none'"
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
PRO INSERT_ALEPH, instrumentModeId, dbinsert=dbinsert, verbose=verbose, version=version

    pos = where(instrumentModeId, count)
    if instrumentModeId lt 31 or instrumentModeId gt 48 then begin
        print,'Error: invalid InstrumentModeId ',instrumentModeId,' [',mode,']'
        return
    endif

    if n_elements(version) eq 0 then version=2

    ; read the save file and re-arrange the coefficients in a 1D array
    if instrumentModeid eq 41 then begin
        ;readcol,'~/SORCE/data/sima_vis_aleph_2011.txt',wavelength, aleph,format='(d,d)'
        readcol,'~/SORCE/data/sima_vis_aleph_2021.txt',wavelength, aleph,format='(d,d)'
        corr = aleph * 0d
    endif else if instrumentModeId eq 45 then begin
        ;readcol,'~/SORCE/data/simb_vis_aleph_2011.txt',wavelength, aleph,format='(d,d)'
        readcol,'~/SORCE/data/simb_vis_aleph_2021.txt',wavelength, aleph,format='(d,d)'
        corr = aleph * 0d
    endif else if instrumentModeid eq 43 then begin
        ;readcol,'~/SORCE/data/sima_uv_aleph_2011.txt',wavelength, aleph,format='(d,d)'
        readcol,'~/SORCE/data/sima_uv_aleph_2021.txt',wavelength, aleph,format='(d,d)'
        corr = aleph * 0d
    endif else if instrumentModeid eq 47 then begin
        ;readcol,'~/SORCE/data/simb_uv_aleph_2011.txt',wavelength, aleph,format='(d,d)'
        readcol,'~/SORCE/data/simb_uv_aleph_2021.txt',wavelength, aleph,format='(d,d)'
        corr = aleph * 0d
    endif else if instrumentModeid eq 44 then begin
        ;readcol,'~/SORCE/data/sima_ir_aleph_2011.txt',wavelength, aleph,format='(d,d)'
        readcol,'~/SORCE/data/sima_ir_aleph_2021.txt',wavelength, aleph,format='(d,d)'
        corr = aleph * 0d
    endif else if instrumentModeid eq 48 then begin
        ; no simb_ir calibrated data available for version 17 - use 1.0
        readcol,'~/SORCE/data/sima_ir_aleph_2011.txt',wavelength, aleph,format='(d,d)'
        aleph = aleph * 0d + 1d
        corr = aleph * 0d
    endif else if instrumentModeid eq 31 or instrumentModeId eq 32 then begin
        ; use the default aleph which aligns very well with the calibrated IR photodiode data
        ;wavelength=[200d,400d,1000d,3000d]
        ;aleph = [1d,1d,1d,1d]
        query_database,'select * from SimAlephFactorCal where calibrationSetId=659', data, nrows
        wavelength=data.wavelength
        aleph=data.correctfactor
        corr=data.correctfactorunc
    endif else begin
        print, 'No data available yet for this instrumentModeId'
        return
    endelse

    query="select * from CalibrationMetadata where calibrationTableName='SimAlephFactorCal' and "
    query+=" version="+strtrim(string(version),2)+" and instrumentModeId="+strtrim(string(instrumentModeId),2)

    query_database,query, res, nrows

    if nrows gt 0 then begin
        calId = res[0].CALIBRATIONSETID
    endif else begin
        ; get a new calibrationSetId since the version to update doesn't exist
        query_database,"select max(calibrationSetId) from CalibrationMetadata", maxcal
        calId = maxcal.(0)+1
    endelse
    print,'  inserting data for calibrationSetId=',calId

    stmt17 = fjava_get_jdbc_statement(user='sorce', password='sorcedb', server='sorce-db', database="SORCE")
    conn17 = OBJ_NEW("oJava_DbExchange", stmt17)
    if keyword_set(dbinsert) then insert=1 else insert=0

    ; first delete all entries in the table for this particular version (version 1 only for now)
    q1 = "DELETE FROM SimAlephFactorCal WHERE calibrationSetId="+strtrim(string(calId),2)
    if keyword_set(verbose) then print,q1
    if insert then result=conn17->execute(q1)

    ; make sure we have an entry in the CalibrationMetadata table for this calibrationSetId
    effectiveDate='2003-01-01 00:00:00.0'
    update_calibMetadata_aleph, instrumentModeId, calId, effectiveDate, version=version, conn17=conn17, dbinsert=insert, verbose=verbose

    ; insert one row at a time
    q1="INSERT INTO dbo.SimAlephFactorCal(calibrationSetId, wavelength, correctFactor, correctFactorUnc) VALUES("

    nwaves=n_elements(WAVELENGTH)
    for i=0L,nwaves-1 do begin
        print,i+1,double(i+1)/double(nwaves+1)*100d,format='("   ",i0,"  (",i0," % )",$,%"\r")'
        q2 = strtrim(string(calId,format='(I)'),2)
        q2 = q2+','+strtrim(string(wavelength[i],format='(d0.8)'),2)
        q2 = q2+','+strtrim(string(aleph[i],format='(d0.8)'),2)
        q2 = q2+','+strtrim(string(corr[i],format='(d0.8)'),2)+')'
        if keyword_set(verbose) then print,q1+q2
        if insert then result=conn17->execute(q1+q2)
    endfor
    print,''

    return

END





