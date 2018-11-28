;+
; Author: Stephane Beland
;
; PURPOSE: 
;  Plots the prism transmission degradation for a set of wavelengths
;  for the specified instrument using the Kappa, raypath
;  and degradation Column in the database.
;
; CALLING SEQUENCE:
;
; INPUT PARAMETERS:
;   startTime -
;      The lower time range for which data will be returned.
;      May be specified as julian days, gps microseconds, or
;      SORCE mission day numbers.
;
;   stopTime -
;      The upper time range for which data will be returned
;      May be specified as julian days, gps microseconds, or
;      SORCE mission day numbers.
;
;   instrumentModeId -
;      Id of the instrument to process.
;
; OPTIONAL INPUT PARAMETERS:
;   version -
;      Requests a specific version of the SimProfileIntegralCal, otherwise
;      pick the latest version valid for the specified time range.
;
; RETURNED PARAMETERS:
;
; OPTIONAL OUTPUT PARAMETERS:
;
; EXAMPLE:  
;
; COMMON BLOCKS:
;      NONE
;
; NOTES:
;
; REVISION HISTORY:
;   Revision: $Id: plot_prismdeg.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;------------------------------------------------------------------

pro plot_prismdeg, instrumentModeId, waves=waves, $
    kappa=kappa, degcol=degcol, raypath=raypath

    ; limit the wavelength range according to the specified instrumentModeId
    if (instrumentModeId ge 41 and instrumentModeId le 44) or instrumentModeId eq 31 then begin
        instrument='SIMA'
        expmode=54
    endif else if (instrumentModeId ge 45 and instrumentModeId le 48) or instrumentModeId eq 32 then begin
        instrument='SIMB'
        expmode=55
    endif else begin
        print,'Error: wrong instrumentModeId (should be between 41 and 48)'
        return
    endelse

    ; get the kappa
    if n_elements(kappa) eq 0 then begin
        query_database, /reset
        q2 = "SELECT calibrationSetId from CalibrationMetadata where calibrationTableName='SimPrismDegKappaCal' "+$
            "and version=19 and instrumentModeId="+strtrim(string(instrumentModeId),2)
        query_database, q2, res, info
        q3='select x,y from SimPrismDegKappaCal where calibrationSetId='+strtrim(string(res.(0)),2)
        query_database, q3, esr_kappa, info
        esr_kappa.x = 10d ^ esr_kappa.x

        if expmode eq 54 then uvmode=43 else uvmode=47
        q2 = "SELECT calibrationSetId from CalibrationMetadata where calibrationTableName='SimPrismDegKappaCal' "+$
            "and version=19 and instrumentModeId="+strtrim(string(uvmode),2)
        query_database, q2, res, info
        q3='select x,y from SimPrismDegKappaCal where calibrationSetId='+strtrim(string(res.(0)),2)
        query_database, q3, uv_kappa, info
        uv_kappa.x = 10d ^ uv_kappa.x
        p0=where(uv_kappa.x le 308d)
        p1 = where(esr_kappa.x gt 308d)
        kappa=[uv_kappa[p0],esr_kappa[p1]]
    endif

    ; get the degradationColumn
    if n_elements(degcol) eq 0 then begin
        query_database, /reset
        q2 = "SELECT calibrationSetId from CalibrationMetadata where calibrationTableName='SimPrismDegColumnCalTable' "+$
            "and version=19 and instrumentModeId="+strtrim(string(expmode),2)
        query_database, q2, res, info
        q3='select x,y from SimPrismDegColumnCalTable where calibrationSetId='+strtrim(string(res.(0)),2)
        query_database, q3, degcol, info
    endif

    ; get the raypath
    if n_elements(raypath) eq 0 then begin
        query_database, /reset
        q2 = "SELECT calibrationSetId from CalibrationMetadata where calibrationTableName='SimRayPathDegradationParams' "+$
            "and version=4 and effectiveDate<='2004-04-21 00:00:00.0' and instrumentModeId="+strtrim(string(instrumentModeId),2)
        query_database, q2, res, info
        q3='select * from SimRayPathDegradationParams where calibrationSetId='+strtrim(string(res.(0)),2)
        query_database, q3, raypath, info
    endif

    if n_elements(waves) eq 0 then $
        waves=[220d, 230, 280, 350, 450, 550, 650, 850, 950, 1200, 1400, 1600, 2200, 2400]

    sd = gps2jd(degCol.x)

    res=label_date(date_format='%M %Y')

    for w=0, n_elements(waves)-1 do begin
       w_kappa = interpol(kappa.y, kappa.x, waves[w], /lsq)
       w_raypath0 = interpol(raypath.singlePassAreaFraction, raypath.wavelength, waves[w], /spline)
       w_raypath1 = interpol(raypath.firstSurfaceDegradation, raypath.wavelength, waves[w])
       prismDeg = (1d - w_raypath0[0]) * exp(-w_kappa[0] * degCol.y) + w_raypath0[0] * exp(-w_kappa[0] * degCol.y * w_raypath1[0])
       ;prismDeg = exp(-w_kappa[0] * degCol.y) 
       lineplot,sd,prismDeg, xtitle='Date', ytitle='Prism Degradation',title='Wavelength='+strtrim(string(waves[w],format='(F0.1)'),2),charsize=1.4,ptitle=instrument+' Prism Degradation Trend',xtickformat='LABEL_DATE'
    endfor

end
