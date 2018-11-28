;+
; NAME:   INSERT_DETDEG
;
; PURPOSE: 
;    This routine populates the SORCE/SimDiodeDegCalPolynomialCoeffs database table with the 
;    new detector degradation column values obtained from comparing the degradation of the 
;    specific diode to the corresponding ESR data (using adjust_diodedeg.pro and fitting a
;    surface to the combined SimA and SimB results).
;
; CALLING SEQUENCE:
;    insert_detdeg, 41, /dbinsert
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
;    The coefficients is an array of 16 values and they are applied as:
;  
;  deg = c[0] + c[1]*w + c[2]*w^2 + c[3]*w^3 + c[4]*t + c[5]*t*w + c[6]*t*w^2 + c[7]*t*w^3 +
;        c[8]*t^2 + c[9]*t^2*w + c[10]*t^2*w^2 + c[11]*t^2*w^3 + c[12]*t^3 + c[13]*t^3*w + 
;        c[14]*t^3*w^2 + c[15]*t^3*w^3
;
;    where t is in microSecondsSinceGPS and w is in nm for simab_vis_deg_coeffs.sav
;    where t is in SORCE MISSION DAYS and w is in nm for simab_vis_deg_coeffs_afact48.sav
;
; REVISION HISTORY:
;   Revision: $Id: insert_detdeg.pro,v 1.2 2018/11/26 16:16:02 spenton Exp $
;-
;*****************************************************************************
PRO update_calibMetadata_detdeg, instMode, calId, effDate, version=version, $
    ncoeffs=ncoeffs, verbose=verbose, dbinsert=insert

    ; make sure we have an entry in the CalibrationMetadata table for this calibrationSetId
    q1="SELECT * from CalibrationMetadata where calibrationSetId="+strtrim(string(CalId),2)
    if insert then begin
        ;query_database,/reset
        query_database,q1,data,nrows
        date_str=string(systime(/jul),format='(C(CYI,"-",CMOI2.2,"-",CDI2.2," ",CHI2.2,":",CMI2.2,":",CSI2.2))')
        if nrows eq 0 then begin
            ; insert a new row in the CalibrationMetadata table
            q1="INSERT INTO dbo.CalibrationMetadata(calibrationSetId, calibrationTableName, instrumentModeId, "
            q1+="version, releaseDate, effectiveDate, rationaleForNewVersion, responsibleIndividual, "
            q1+="dataDescription, sourceOfData, applicableReferences, additionalInformation) VALUES("
            q1=q1+strtrim(string(CalId),2)+", 'SimDiodeDegCalPolynomialCoeffs', "
            q1=q1+strtrim(string(instMode),2)+", "+strtrim(string(version),2)
            q1=q1+", '"+date_str+"', '"+effDate+"'"
            q1=q1+", 'diode degradation from combined SimA and SimB compared to ESR data for version 23'"
            q1=q1+", 'Stephane Beland'"
            ;q1=q1+", 'Stores the nxn polynomial coefficients representing the surface over the gpsmicroseconds and wavelength'"
            q1=q1+", 'Stores the NxN polynomial coefficients representing the surface over the MISSION DAYS and WAVELENGTH'"
            q1=q1+", 'IDL save file with NxN coefficients'"
            q1=q1+", 'none'"
            q1=q1+", 'none')"
            if keyword_set(verbose) then print,q1
            if insert then query_database, q1, result

            ; create a new entry in the SimDiodeDegCalModels
            q2="INSERT INTO SimDiodeDegCalModels(calibrationSetId, node0, nodeIncrement, nCoeffs) VALUES("+strtrim(string(calId),2)
            q2=q2+", 0, 0, "+strtrim(string(ncoeffs),2)+")"
            if keyword_set(verbose) then print,q2
            if insert then query_database, q2, result

            ; create a new entry in the SimDiodeDegCalSplineCoeffs
            q3="INSERT INTO SimDiodeDegCalSplineCoeffs(calibrationSetId, splineCoeffIndex, polynomialDegree) VALUES("+strtrim(string(calId),2)
            q3=q3+", 0, "+strtrim(string(ncoeffs),2)+")"
            if keyword_set(verbose) then print,q3
            if insert then query_database, q3, result

        endif else begin
            q1="UPDATE CalibrationMetadata SET releaseDate='"+date_str+"'"
            q1=q1+", effectiveDate='"+effDate+"'"
            q1=q1+", version="+strtrim(string(version),2)
            q1=q1+" WHERE calibrationSetId="+strtrim(string(CalId),2)
            if keyword_set(verbose) then print,q1
            if insert then query_database, q1, result

            q2="UPDATE SimDiodeDegCalSplineCoeffs SET splineCoeffIndex=0, polynomialDegree="+strtrim(string(ncoeffs),2)
            q2=q2+" WHERE calibrationSetId="+strtrim(string(CalId),2)
            if keyword_set(verbose) then print,q2
            if insert then query_database, q2, result

        endelse
    endif

END

;-----------------------------------------------------------------------------
PRO INSERT_DETDEG, instrumentModeId, dbinsert=dbinsert, verbose=verbose, version=version

    pos = where(instrumentModeId, count)
    if instrumentModeId lt 41 or instrumentModeId gt 48 then begin
        print,'Error: invalid InstrumentModeId ',instrumentModeId,' [',mode,']'
        return
    endif

    if n_elements(version) eq 0 then version=2

    ; read the save file and re-arrange the coefficients in a 1D array
    if instrumentModeid eq 41 or instrumentModeId eq 45 then begin
        ;restore,'/Users/sbeland/SORCE/data/simab_vis_deg_coeffs.sav'
        ;restore,'/Users/sbeland/SORCE/data/simab_vis_deg_coeffs_afact48.sav'
        ;restore,'/Users/sbeland/SORCE/data/simab_vis_deg_coeffs_mod_407nm.sav'
        ;coeffs = double(reform(visab_coeffs, n_elements(visab_coeffs)))
        restore,'/Users/sbeland/SORCE/data/simab_vis_deg_coeffs_mod_407nm_v23.sav'
        coeffs = double(reform(sfit_coeffs, n_elements(sfit_coeffs)))
        ncoeffs=n_elements(coeffs)
    endif else if instrumentModeid eq 43 or instrumentModeId eq 47 or $
                  instrumentModeid eq 44 or instrumentModeId eq 48 then begin
        ; use the same size as VIS degradation but zero out all values except coeffs[0]=1.0
        ;restore,'/Users/sbeland/SORCE/data/simab_vis_deg_coeffs_afact48.sav'
        restore,'/Users/sbeland/SORCE/data/simab_vis_deg_coeffs_mod_407nm.sav'
        coeffs = double(reform(visab_coeffs, n_elements(visab_coeffs)))
        ncoeffs=n_elements(coeffs)
        coeffs *= 0d
        coeffs[0] = 1.0d
    endif else begin
        print, 'No data available yet for this instrumentModeId'
        return
    endelse

    query="select * from CalibrationMetadata where calibrationTableName='SimDiodeDegCalPolynomialCoeffs' and "
    query+=" version="+strtrim(string(version),2)+" and instrumentModeId="+strtrim(string(instrumentModeId),2)

    query_database,query, res, nrows, user='sorce', password='sorcedb', server='sorce-db', database="SORCE"

    if nrows gt 0 then begin
        calId = res[0].CALIBRATIONSETID
    endif else begin
        ; get a new calibrationSetId since the version to update doesn't exist
        query_database,"select max(calibrationSetId) from CalibrationMetadata", maxcal
        calId = maxcal.(0)+1
    endelse

    if keyword_set(dbinsert) then insert=1 else insert=0

    ; first delete all entries in the table for this particular version (version 1 only for now)
    q1 = "DELETE FROM SimDiodeDegCalPolynomialCoeffs WHERE calibrationSetId="+strtrim(string(calId),2)
    if keyword_set(verbose) then print,q1
    if insert then query_database,q1,result

    ; make sure we have an entry in the CalibrationMetadata table for this calibrationSetId
    effectiveDate='2003-01-01 00:00:00.0'
    update_calibMetadata_detdeg, instrumentModeId, calId, effectiveDate, ncoeffs=ncoeffs, $
        version=version, dbinsert=insert, verbose=verbose

    ; insert one row at a time
    q1="INSERT INTO dbo.SimDiodeDegCalPolynomialCoeffs(calibrationSetId, splineCoeffIndex, coeffIndex, coeffValue, coeffUncertainty) VALUES("

    for i=0L,ncoeffs-1 do begin
        print,i+1,double(i+1)/double(ncoeffs+1)*100d,format='("   ",i0,"  (",i0," % )",$,%"\r")'
        str = strtrim(string(calId,format='(I)'),2)
        str = str+', 0'
        str = str+', '+strtrim(string(i),2)
        str = str+', '+strtrim(string(coeffs[i],format='(G25.17)'),2)
        str = str+', 0.0'
        q2 = q1+str+')'
        if keyword_set(verbose) then print,q2
        if insert then query_database,q2,result
    endfor
    print,''

    return

END





