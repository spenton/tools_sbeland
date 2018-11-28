;+
; NAME:   INSERT_FDEG_KAPPA
;
; PURPOSE: 
;    This routine populates the SORCE/SimPrismDegColumnCalCoeffs database table with the 
;    new COMBINED F_Factor*Kappa values obtained from get_visdiode_fdeg.pro.
;    The file read contains the coefficients of a 3D surface representing the values of 
;    C=F_Factor * Kappa as a function of wavelength and time (Mission DAYS).
;    These are picked up by the Java code to determine the prism transmission degradation.
;
;    NOTE: This is used starting with version 20.  Version 19 included the F_Factor with the
;          degradation column (SimPrismDegColumnCalTable) and used a constant Kappa throughout
;          the mission.  In version 20, we combine F_Factor * Kappa and let this value
;          change with both mission day and wavelength.
;
; CALLING SEQUENCE:
;    insert_fdeg_kappa, 41, /dbinsert
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
;    where t is in SORCE MISSION DAYS and w is in nm
;
; REVISION HISTORY:
;   Revision: $Id: insert_fdeg_kappa.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;-
;*****************************************************************************
PRO update_calibMetadata_fdeg, instMode, calId, effDate, tablename, version=version, $
    conn17=conn17, ncoeffs=ncoeffs, verbose=verbose, dbinsert=insert

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
            q1=q1+strtrim(string(CalId),2)+", '"+strtrim(tablename,2)+"', "
            q1=q1+strtrim(string(instMode),2)+", "+strtrim(string(version),2)
            q1=q1+", '"+date_str+"', '"+effDate+"'"
            q1=q1+", 'Combined F_Factor and Kappa vs MissionDays and wavelength from SimA and SimB data'"
            q1=q1+", 'Stephane Beland'"
            ;q1=q1+", 'Stores the nxn polynomial coefficients representing the surface over the gpsmicroseconds and wavelength'"
            q1=q1+", 'Stores the nxn polynomial coefficients representing the F*Kappa surface over the MISSION DAYS and WAVELENGTH'"
            q1=q1+", 'IDL save file with nxn coefficients (simab_vis_fdeg_coeffs_afact48.sav)'"
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
PRO INSERT_FDEG_KAPPA, instrumentModeId, dbinsert=dbinsert, verbose=verbose, version=version

    pos = where(instrumentModeId, count)
    if instrumentModeId lt 31 or instrumentModeId gt 48 then begin
        print,'Error: invalid InstrumentModeId ',instrumentModeId,' [',mode,']'
        return
    endif

    if n_elements(version) eq 0 then version=1
    tablename = 'SimFDegKappaPolyCoeffs'

    ; read the save file and re-arrange the coefficients in a 1D array
    if instrumentModeid eq 41 or instrumentModeId eq 45 or instrumentModeId eq 31 or instrumentModeId eq 32 then begin
        restore,'/Users/sbeland/SORCE/data/simab_vis_fdeg_coeffs_afact48.sav'
        coeffs = double(reform(VIS_FDEG_COEFFS, n_elements(VIS_FDEG_COEFFS)))
        ncoeffs=n_elements(coeffs)
    endif else begin
        print, 'No data available yet for this instrumentModeId'
        return
    endelse

    query="select * from CalibrationMetadata where calibrationTableName='"+tablename+"' and "
    query+=" version="+strtrim(string(version),2)+" and instrumentModeId="+strtrim(string(instrumentModeId),2)

    query_database,query, res, nrows

    if nrows gt 0 then begin
        calId = res[0].CALIBRATIONSETID
    endif else begin
        ; get a new calibrationSetId since the version to update doesn't exist
        query_database,"select max(calibrationSetId) from CalibrationMetadata", maxcal
        calId = maxcal.(0)+1
    endelse


    stmt17 = fjava_get_jdbc_statement(user='sorce', password='sorcedb', server='sorce-db', database="SORCE")
    conn17 = OBJ_NEW("oJava_DbExchange", stmt17)
    if keyword_set(dbinsert) then insert=1 else insert=0


    ; first delete all entries in the table for this particular version (version 1 only for now)
    q1 = "DELETE FROM "+tablename+" WHERE calibrationSetId="+strtrim(string(calId),2)
    if keyword_set(verbose) then print,q1
    if insert and nrows gt 0 then result=conn17->execute(q1)

    ; make sure we have an entry in the CalibrationMetadata table for this calibrationSetId
    effectiveDate='2003-01-01 00:00:00.0'
    update_calibMetadata_fdeg, instrumentModeId, calId, effectiveDate, tablename, ncoeffs=ncoeffs, $
        version=version, conn17=conn17, dbinsert=insert, verbose=verbose

    ; insert one row at a time
    q1="INSERT INTO "+tablename+" (calibrationSetId, coeffIndex, coeffValue, coeffUncertainty) VALUES("

    for i=0L,ncoeffs-1 do begin
        print,i+1,double(i+1)/double(ncoeffs+1)*100d,format='("   ",i0,"  (",i0," % )",$,%"\r")'
        str = strtrim(string(calId,format='(I)'),2)
        str = str+', '+strtrim(string(i),2)
        str = str+', '+strtrim(string(coeffs[i],format='(G25.17)'),2)
        str = str+', 0.0'
        q2 = q1+str+')'
        if keyword_set(verbose) then print,q2
        if insert then result=conn17->execute(q2)
    endfor
    print,''

    return

END





