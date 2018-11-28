;+
; NAME:   INSERT_PROINT_19
;
; PURPOSE: 
;    This routine populates the SORCE/SimProfileIntegralCal database table with the 
;    data for version 19.
;
; CALLING SEQUENCE:
;    insert_proint_19, instrumentModeId, /dbinsert
;
; OPTIONAL INPUT PARAMETERS:
;    none
;
; OPTIONAL INPUT KEYWORDS:
;    dbinsert - if set will insert the re-formatted data in the database
;
; OUTPUT PARAMETERS:
;   none
;
; OPTIONAL OUTPUT PARAMETERS:
;      NONE
;
; EXAMPLE:  
;
; COMMON BLOCKS:
;      NONE
;
; NOTES:
;
; REVISION HISTORY:
;   Revision: $Id: insert_proint_19.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;-
;*****************************************************************************
PRO INSERT_PROINT_19, instrumentModeId, dbinsert=dbinsert, verbose=verbose
    ; database table SORCE.dbo.SimProfileIntegralCal

    ; get the latest calibrationSetId from the CalibrationMetadata table

    ; prepare the database connection
    jstmt = fjava_get_jdbc_statement(user='sorce', password='sorcedb', server='sorce-db', database='SORCE')
    sorceDbExchange = OBJ_NEW("oJava_DbExchange", jstmt)
    
    query = 'SELECT max(calibrationSetId) as calibId FROM CalibrationMetadata'
    res=sorceDbExchange->getAllvalues(query)
    if size(res,/tname) ne 'STRUCT' then begin
        print,'can not access CalibrationMetadata table'
        exit,-1
    endif
    calibId = fix(re[0].calibId) + 1

    ; get the reference entry for version 19 (corresponding to version 17)
    query = "SELECT calibrationSetId FROM dbo.CalibrationMetadata where calibrationTableName='SimProfileIntegralCal' and "
    query += "instrumentModeId="+strtrim(string(instrumentModeId),2) 
    query += " and version=19 and effectiveDate<='2004-01-01'"
    res=sorceDbExchange->getAllvalues(query)
    ref_calibId = res[0].calibrationSetId
    query = "SELECT * FROM SimProfileIntegralCal where calibrationSetId="+strtrim(string(ref_calibId),2)
    query += " order by y4"
    ref_data=sorceDbExchange->getAllvalues(query)

    ; get the list of rows we'll need to copy from version 18 to version 19
    query = "SELECT * FROM dbo.CalibrationMetadata where calibrationTableName='SimProfileIntegralCal' and "
    query += "instrumentModeId="+strtrim(string(instrumentModeId),2) 
    query += " and version=18 and effectiveDate>='2004-01-01'"
    res18=sorceDbExchange->getAllvalues(query)

    ; loop through every calibrationSetId
    for i=0,n_elements(res)-1 do begin
        ; extract the data from the SimProfileIntegralCal for this calibrationSetId
        query = "SELECT * FROM SimProfileIntegralCal where calibrationSetId="+strtrim(string(res18[i].calibrationSetId),2)
        query += " order by y4"
        data=sorceDbExchange->getAllvalues(query)

        ; the ccd position for a specific wavelength needs to stay the same throughout
        ; since we are now applying a ccd shift based on the date
        ; keep the wavelength the same for the new row but interpolate the corresponding
        ; ccd position from the reference row (from day 453).
        ccdpos = interpol(ref_data.x, ref_data.y4, data.y4,/double)
        data.x = ccdpos

        ; insert the data in the CalibrationMetadata table
        res18[i].calibrationsetid = calibId
        res18[i].version = 19
        res18[i].releaseDate = jd2syb(systime(/jul,/utc))
        insert_str=''
        for j=0,n_tags(res18[i])-1 do begin
            if j eq 0 then $
                insert_str = insert_str + '"' + strtrim(string(res18[i].(j)),2) + '"' $
            else $ 
                insert_str = insert_str + ', "' + strtrim(string(res18[i].(j)),2) + '"' $
        endfor
        stmt = "INSERT INTO CalibrationMetadata VALUES (" + insert_str + ")"

        ; now insert the data in the SimProfileIntegralCal table
        q1 = "INSERT INTO SimProfileIntegralCal VALUES("
        for j=0L,n_elements(data)-1L do begin
            for k=0,n_tags(data)-1 do begin
                if k eq 0 then $
                    insert_str = strtrim(string(data[j].(k),format='(E20.12)'),2) $
                else $ 
                    insert_str = insert_str + ', ' + strtrim(string(data[j].(k),format='(E20.12)'),2) $
            endfor
            q2 = q1+insert_str+')'
            if keyword_set(verbose) then print,q2
            if insert then res=sorceDbExchange->execute(q2)
        endfor

    endfor

END

