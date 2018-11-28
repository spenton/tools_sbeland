;+
; NAME:   INSERT_RAYPATH
;
; PURPOSE: 
;    This routine populates the SORCE/SimRayPathDegradationParams database table with the 
;    new RayPath degradation factor obtained when comparing the irradiance of the ESR and
;    the photodiode.  The RayPath degradation contains the change in degradation from the 
;    ESR to the photodiode but also the diode degradation with time.
;
; CALLING SEQUENCE:
;    insert_raypath, 54, /dbinsert
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
;   The raypath values are obtained from the code adjust_diodedeg.pro
;
; REVISION HISTORY:
;   Revision: $Id: insert_raypath.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;-
;*****************************************************************************
FUNCTION INSERT_RAYPATH, instrumentModeId, dbinsert=dbinsert, verbose=verbose

    ;INSERT INTO dbo.SimRayPathDegradationParams(calibrationSetId, wavelength, singlePassAreaFraction, 
    ;    singlePassAreaFractionUnc, firstSurfaceDegradation, firstSurfaceDegradationUnc) VALUES(0, 0, 0, 0, 0, 0)

    modes = [31,32,41,43,44,45,47,48]
    pos = where(modes eq instrumentModeId, count)
    if count eq 0 then begin
        print,'Error: invalid InstrumentModeId ',instrumentModeId,' [',mode,']'
        return,-1
    endif

    if instrumentModeId eq 31 or instrumentModeId eq 32 then begin
        w0=240.0d
        w1=3000.0d
        step=10d
        npts=ceil((w1-w0)/step)
        wavelength = dindgen(npts+1)*step + w0
        singlePass = wavelength * 0d
        ; read the file with the ray trace results to start with
        readcol,'~/SORCE/data/overlap_esr.txt',wavelength,singlePass,format='(d,d)'
        singlePass=abs(singlePass)
        firstSurf = singlePass*0d +0.5d
        if instrumentModeId eq 31 then newCalId=1043 else newCalId=1044
    endif else if instrumentModeId eq 41 or instrumentModeId eq 45 then begin
        w0=300.0d
        w1=1020.0d
        ;step=0.25d
        step=10d
        npts=ceil((w1-w0)/step)
        wavelength = dindgen(npts+1)*step + w0
        ; adjusted afact for VIS for the VIS kappa to match the UV kappa (with the uv afact=0.0)
        singlePass = wavelength*0d +0.4d
        ;dd=[0.00053165011d, 0.0090682349d, 0.90287432d]
        ;  use the latest a-factor from the kappa and f-funct after 20 iterations
        ;dd=[0.00039227191d, 0.0094851024d, 0.89668452d]
        ; updated A-Factor with VIS processed without the Y11 in profileIntegral
        ;dd=[0.00039639601d, 0.0094760814d, 0.89548653d]
        ;singlePass = 1d - (dd[0]*exp(dd[1]*wavelength)+dd[2])

        ;dd=[6.3133921d-05, 0.011147337d, 0.99384342d]
        ;singlePass = 1d - (dd[0]*exp(dd[1]*wavelength)+dd[2])
        firstSurf = singlePass*0d +0.5d
        if instrumentModeId eq 41 then newCalId=1045 else newCalId=1046
    endif else if instrumentModeId eq 43 then begin
        w0=195.0d
        w1=310.0d
        ;step=0.1d
        step=1d
        npts=ceil((w1-w0)/step)
        wavelength = dindgen(npts+1)*step + w0
        singlePass = wavelength * 0d
        ; use latest a-factor from iterating 20 times on kappa and ffunct
        ;dd=[-2.2680018d, 0.007046815d]
        ; updated A-Factor processed without the Y11 in profileIntegral
        ;dd=[-1.9258139d, 0.005939250d]
        ;singlePass = poly(wavelength,dd)
        ;dd=[-3.8336794d, 0.0118392867d]
        ;
        ; new linear model of UV raypath for data at SD>1570 (20130524)
        ;dd=[1.633333d, -0.0041666667d]
        ;singlePass = poly(wavelength,dd)
        firstSurf = singlePass*0d +0.5d
        newCalId=1049
    endif else if instrumentModeId eq 47 then begin
        w0=195.0d
        w1=310.0d
        ;step=0.1d
        step=1d
        npts=ceil((w1-w0)/step)
        wavelength = dindgen(npts+1)*step + w0
        singlePass = wavelength * 0d
        ;
        ; new linear model of UV raypath for data at SD>1570 (20130524)
        ;dd=[1.633333d, -0.0041666667d]
        ;singlePass = poly(wavelength,dd)
        firstSurf = singlePass*0d +0.5d
        newCalId=1050
    endif else if instrumentModeId eq 44 or instrumentModeId eq 48 then begin
        w0=850.0d
        w1=1700.0d
        step=1d
        npts=ceil((w1-w0)/step)
        wavelength = dindgen(npts+1)*step + w0
        ;singlePass = wavelength*0d
        ; polynomial fit for the IR photodiode
        ;dd=[2.9015510d, -0.0053566140d, 1.8269480d-06]
        ; use latest a-factor from iterating 20 times on kappa and ffunct
        ;
        ; updated A-Factor processed without the Y11 in profileIntegral
        ;dd=[3.6328230d, -0.0065217392d, 2.2113357d-6]
        ;singlePass = poly(wavelength,dd)
        dd=[0.75693086d, -0.0013653438d, 4.8538028d-7]
        ;dd=[0.025562854d, -0.00015262944d]
        singlePass = poly(wavelength,dd)
        firstSurf = singlePass*0d +0.5d
        if instrumentModeId eq 44 then newCalId=1047 else newCalId=1048
    endif
    version=3

    stmt17 = fjava_get_jdbc_statement(user='sorce', password='sorcedb', server='sorce-db', database="SORCE")
    conn17 = OBJ_NEW("oJava_DbExchange", stmt17)
    if keyword_set(dbinsert) then insert=1 else insert=0

    ; make sure we have an entry in the CalibrationMetadata table for this calibrationSetId
    q1="SELECT * from CalibrationMetadata where calibrationSetId="+strtrim(string(newCalId),2)
    if insert then begin
        query_database,/reset
        query_database,q1,data,nrows
        date_str=string(systime(/jul),format='(C(CYI,"-",CMOI2.2,"-",CDI2.2," ",CHI2.2,":",CMI2.2,":",CSI2.2))')
        if nrows eq 0 then begin
            ; insert a new row in the CalibrationMetadata table
            q1="INSERT INTO dbo.CalibrationMetadata(calibrationSetId, calibrationTableName, instrumentModeId, "
            q1+="version, releaseDate, effectiveDate, rationaleForNewVersion, responsibleIndividual, "
            q1+="dataDescription, sourceOfData, applicableReferences, additionalInformation) VALUES("
            q1=q1+strtrim(string(newCalId),2)+", 'SimRayPathDegradationParams', "
            q1=q1+strtrim(string(instrumentModeId),2)+", "+strtrim(string(version),2)
            q1=q1+", '"+date_str+"', '2003-01-01 00:00:00.0'"
            q1=q1+", 'New RayPath for version 19 from updated wavelength knowledge and new KappaDegradation'"
            q1=q1+", 'Stephane Beland'"
            q1=q1+", 'RayPath degradation as (1-A)*exp(K*degCol)+A*exp(K*degCol*B) '"
            q1=q1+", 'Stephane processing and ESR to diode comparison'"
            q1=q1+", 'Generated by straightening the difference in diode and ESR signal time series as a function of wavelength'"
            q1=q1+", 'none')"
            if keyword_set(verbose) then print,q1
            if insert then result=conn17->execute(q1)
        endif else begin
            q1="UPDATE CalibrationMetadata SET applicableReferences="
            q1=q1+"'Generated by straightening the difference in diode and ESR signal time series as a function of wavelength'"
            q1=q1+", releaseDate='"+date_str+"'"
            q1=q1+" WHERE calibrationSetId="+strtrim(string(newCalId),2)
            if keyword_set(verbose) then print,q1
            if insert then result=conn17->execute(q1)
        endelse
    endif

    ; first delete all entries in the table for this particular version (version 1 only for now)
    q1 = "DELETE FROM SimRayPathDegradationParams WHERE calibrationSetId="+strtrim(string(newcalId),2)
    if keyword_set(verbose) then print,q1
    if insert then result=conn17->execute(q1)

    ; insert one row at a time
    q1="INSERT INTO dbo.SimRayPathDegradationParams(calibrationSetId, wavelength, singlePassAreaFraction, "+$
       "singlePassAreaFractionUnc, firstSurfaceDegradation, firstSurfaceDegradationUnc) VALUES("
    for i=0L,n_elements(wavelength)-1 do begin
        str = strtrim(string(newCalId,format='(I)'),2)
        str = str+', '+strtrim(string(wavelength[i]),2)
        str = str+', '+strtrim(string(singlePass[i]),2)
        str = str+', 0.0'
        str = str+', '+strtrim(string(firstSurf[i]),2)
        str = str+', 0.0'
        q2 = q1+str+')'
        if keyword_set(verbose) then print,q2
        if insert then result=conn17->execute(q2)
    endfor

    return, {wavelength:wavelength, singlePassAreaFraction:singlePass, firstSurfaceDegradation:firstSurf}

END

