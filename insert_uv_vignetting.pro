;+
; NAME:   INSERT_UV_VIGNETTING
;
; PURPOSE: 
;    This routine populates the SORCE/SimVignettingCal database table with the 
;    new vignetting for the UV diode before day 192 (cause by mirror not fully 
;    removed from optical path).
;
; CALLING SEQUENCE:
;    res=insert_uv_vignetting(/dbinsert)
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
;   The vignetting is obtained form the CORRECTED CCD POSITIONS, and must
;  also interpolated in the same reference frame.
;
; REVISION HISTORY:
;   Revision: $Id: insert_uv_vignetting.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;-
;*****************************************************************************
PRO update_vignetMetadata, instMode, calId, effDate, version=version, conn17=conn17, $
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
            q1=q1+strtrim(string(CalId),2)+", 'SimVignettingCal', "
            q1=q1+strtrim(string(instMode),2)+", "+strtrim(string(version),2)
            q1=q1+", '"+date_str+"', '"+effDate+"'"
            q1=q1+", 'need a correction for early mission'"
            q1=q1+", 'Stephane Beland'"
            q1=q1+", 'SIM UV Vignetting correction when cal_in=0'"
            q1=q1+", 'Stephane early mission data'"
            q1=q1+", 'none'"
            q1=q1+", 'Generated from the CORRECTED PRISM POSITIONS (apply to the corrected positions)')"
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
PRO my_insertVignetting, CalId, ccdpos, correction, conn17=conn17, verbose=verbose, dbinsert=insert

    ; insert one row at a time
    q1="INSERT INTO dbo.SimVignettingCal(calibrationSetId, ccdPos, correction) VALUES("+strtrim(string(CalId,format='(I)'),2)+', '
    npts=n_elements(ccdpos)
    for k=0L,npts-1 do begin
        print,k+1,double(k+1)/double(npts+1)*100d,format='("   ",i0,"  (",i0," % )",$,%"\r")'
        str = strtrim(string(ccdpos[k],format='(F)'),2)
        str = str+', '+strtrim(string(correction[k],format='(F)'),2)
        q2 = q1+str+')'
        if keyword_set(verbose) then print,q2
        if insert then result=conn17->execute(q2)
    endfor
    print,''

END
;-----------------------------------------------------------------------------
FUNCTION INSERT_UV_VIGNETTING, instrumentModeId, version=version, dbinsert=dbinsert, verbose=verbose

    ;INSERT INTO dbo.SimVignettingCal(calibrationSetId, ccdpos, correction) VALUES(0, 0, 0)

    if n_elements(version) eq 0 then version=1

    if instrumentModeId eq 43 then begin
        readcol,'~/SORCE/data/sima_uv_vignetting_4.txt',ccdpos,correction,format='(d,d)'
    endif else if instrumentModeId eq 47 then begin
        readcol,'~/SORCE/data/simb_uv_vignetting_4.txt',ccdpos,correction,format='(d,d)'
    endif else begin
        print,'Error: can only process modes 43 and 47'
        return,-1
    endelse

    query="select * from CalibrationMetadata where calibrationTableName='SimVignettingCal' and "
    query+=" version="+strtrim(string(version),2)+" and instrumentModeId="+strtrim(string(instrumentModeId),2)

    query_database,query, res, nrows

    if nrows gt 0 then begin
        calId = res[0].CALIBRATIONSETID
    endif else begin
        ; get a new calibrationSetId since the version to update doesn't exist
        query_database,"select max(calibrationSetId) from CalibrationMetadata", maxcal
        calId = maxcal.(0)+1
    endelse

    effectiveDate='2003-01-25 00:00:00.0'

    stmt17 = fjava_get_jdbc_statement(user='sorce', password='sorcedb', server='sorce-db', database="SORCE")
    conn17 = OBJ_NEW("oJava_DbExchange", stmt17)
    if keyword_set(dbinsert) then insert=1 else insert=0

    update_vignetMetadata, instrumentModeId, CalId, effectiveDate, version=version, $
        conn17=conn17, dbinsert=insert, verbose=verbose

    ; first delete all entries in the table for this particular version 
    q1 = "DELETE FROM SimVignettingCal WHERE calibrationSetId="+strtrim(string(CalId),2)
    if keyword_set(verbose) then print,q1
    if insert then result=conn17->execute(q1)

    my_insertVignetting, CalId, ccdpos, correction, conn17=conn17, verbose=verbose, dbinsert=insert


    return, {calibrationSetId:calId, ccdpos:ccdpos, correction:correction}

END

