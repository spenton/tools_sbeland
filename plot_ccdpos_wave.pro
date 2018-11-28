;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Plot the Wavelength vs CCD spot position for the specified instrument/detector.
;
; CALLING SEQUENCE:
;   PLOT_CCDPOS_WAVE, sima=sima, simb=simb, detector=['uv','vis1','ir','esr']
;
; INPUT PARAMETERS:
;   None -
;      If no parameters are provided, will default to SimA 
;      with uv,vis1, vis2, ,ir and esr.
;
; OPTIONAL INPUT PARAMETERS:
;   sima -
;      If set, will process the data from SimA
;   simb -
;      If set, will process the data from SimB
;   detector - 
;      String array with the list of detectors to plot.
;      Possible values are uv, vis1, vis2, ir, esr.
;   version -
;      Version number of the profileIntegral to use.
;      Defaults to the latest one.
;
; OPTIONAL OUTPUT PARAMETERS:
;      NONE
;
; EXAMPLE:  
;
; COMMON BLOCKS:
;      NONE
;
; REVISION HISTORY:
;   Revision: $Id: plot_ccdpos_wave.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;-
;
pro plot_ccdpos_wave, sima=sima, simb=simb, detector=detector, version=version

    detector_list=["UV","VIS1","VIS2","IR","ESR"]
    if n_elements(detector) ne 0 then begin
        detector=strupcase(detector)
        ; validate the list of detectors
        for i=0,n_elements(detector)-1 do begin
           p=where(detector_list eq detector[i],count)
           if count eq 0 then begin
               print,'Invalid detector: ',detector[i]
               return
           endif
        endfor
    endif else detector=detector_list

    instrument='SIM_A'
    if keyword_set(simb) then instrument='SIM_B'

    ; prepare the database connection 
    jstmt = fjava_get_jdbc_statement(user=name, password=password, dburl=dbURL, $
        dbdriver=dbdriver, server=server, database="SORCE")
    oJavaDbExchange = OBJ_NEW("oJava_DbExchange", jstmt)

    ; process each detector at a time
    for det=0,n_elements(detector)-1 do begin
        ;first get the corresponding instrumentModeId
        query1 = "SELECT instrumentModeId FROM InstrumentModes WHERE instrument='"+instrument+"'"
        query1 = query1+" AND channel='"+detector[det]+"'"
        res=oJavaDbExchange->getAllValues(query1[0])
        if size(res,/tname) ne 'STRUCT' then begin
            print,'Error: no corresponding instrumentModeId'
            return
        endif
        instrumentModeId=res[0].instrumentModeId

        ; get the corresponding calibrationSetId for this instrument mode
        query1 = "SELECT calibrationSetId FROM CalibrationMetadata WHERE "
        query1 = query1+"calibrationTableName='SimProfileIntegralCal' AND instrumentModeId="+strtrim(string(instrumentModeId),2)
        if keyword_set(version) then query1=query1+" AND version<="+strtrim(string(version),2)
        query1 = query1+" ORDER BY version DESC"
        res=oJavaDbExchange->getAllValues(query1[0])
        if size(res,/tname) ne 'STRUCT' then begin
            print,'Error: no CalibrationMetadata found for the requested date/version'
            return
        endif
        calibrationSetId=res[0].calibrationSetId

        ; now get the SimProfileIntegralCal entry for this calibrationSetId
        query2="SELECT * FROM SimProfileIntegralCal WHERE calibrationSetId="+strtrim(string(calibrationSetId),2)
        profile_data=oJavaDbExchange->getAllValues(query2)

        spotRef = profile_data.x
        wavelengthRef = profile_data.y4

        ; plot the data
        lineplot,spotRef, wavelengthRef, ptitle=instrument+' Wavelength vs CCDPOS',$
           xtitle='CCD Spot Position (sub-pixels)',ytitle='Wavelength (nm)',title=detector[det]

    endfor

 
    return

end
