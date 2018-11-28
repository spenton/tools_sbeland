;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Retreive and format the prism degradation calculated and saved during
;   the daily processing of the production code. The data is restored from
;   the database tables Wavelength and SimPrismTransDegradation according
;   to the instrumentModeId and version for the time range requested.
;
; CALLING SEQUENCE:
;   result=get_prismdeg(80,3850,41,version=20)
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
;      The instrument ID (for SIM: 31,32,41,43,44,45,47,48)
;
; OPTIONAL INPUT PARAMETERS:
;   version -
;      Requests a specific version of the SimProfileIntegralCal, otherwise
;      pick the latest version valid for the specified time range.
;   indata -
;      Structure provided as input, used in combination with wavelength.
;   wavelength -
;      Calculates the prism degradation at the specific wavelength.
;
; RETURNED VALUE:
;   result -
;      Returns a structure with an array for the wavelength, an array for
;      starting date/time of the scan for which the degradation was calculated, 
;      and a 2D array of the degradation with one dimension for the date, the
;      other for the wavelength.
;
; OPTIONAL OUTPUT PARAMETERS:
;   prismdeg-
;      If wavelength is specified, the prism degradation is caluculate
;      for that wavelength and returned as an array
;
; EXAMPLE:  
;
; COMMON BLOCKS:
;      NONE
;
; NOTES:
;
; REVISION HISTORY:
;   Revision: $Id: get_prismdeg.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;------------------------------------------------------------------

function get_prismdeg, startTime, stopTime, instrumentModeId, version=version, help=help, $
    mission=mission, julian=julian, gps=gps, verbose=verbose, indata=indata, wavelength=wavelength, $
    prismdeg=prismdeg

    if size(indata,/tname) ne 'STRUCT' then begin
        if n_params() lt 3 or keyword_set(help) then begin
            doc_library, 'get_prismdeg'
            return,-1
        endif

        match, instrumentModeId, [31,32,41,43,44,45,47,48],a
        if a[0] lt 0 then begin
            doc_library, 'get_prismdeg'
            return,-1
        endif

        if n_elements(version) eq 0 then version=20

        if keyword_set(verbose) then print,'Extracting the prismTransDegradation from the database ...'
        dbTables=['Wavelength','SimPrismTransDegradation']
        res=get_science_product(dbTables,startTime, stopTime, instrumentModeId, nrows,$
            mission_days=mission_days, julian_days=julian_days, gps=gps,version=version)

        if size(res,/tname) ne 'STRUCT' then begin
            print,'No prism degradation data found for the requested time range'
            return,-1
        endif

        if keyword_set(verbose) then print,'Got '+strtrim(string(n_elements(res)),2)+' rows  back ...'
        ; sort the data per scan (look at data spacing of more than 1 minute)
        sd=gps2sd(res.MICROSECONDSSINCEGPSEPOCH/1d6)
        diff=sd[1:-1]-sd[0:-2]
        k=where(diff gt 1d/60d/24d)
        startScan=[sd[0],sd[k]+1]
        if keyword_set(verbose) then print,'   from '+strtrim(string(n_elements(startScan)),2)+' scans ...'
        
        ; get the default set of wavelengths
        plot_simspect,outdata=outdata,/noplot
        if instrumentModeId eq 31 or instrumentModeId eq 32 then begin
            allwaves=outdata.esrtable.wavelength
        endif else if instrumentModeId eq 41 or instrumentModeId eq 45 then begin
            allwaves=outdata.visspect.wavelength
        endif else if instrumentModeId eq 43 or instrumentModeId eq 47 then begin
            allwaves=outdata.uvspect.wavelength
        endif else if instrumentModeId eq 44 or instrumentModeId eq 48 then begin
            allwaves=outdata.irspect.wavelength
        endif

        if keyword_set(verbose) then print,'Reformatting the data ...'
        data=dblarr(n_elements(startScan), n_elements(allwaves))
        k0=0L
        for i=0L,n_elements(startScan)-2L do begin 
            k1=k[i] 
            scantimes=sd[k0:k1]
            if k1-k0 lt 3 then splval=0 else splval=1
            scandeg = interpol(res[k0:k1].prismtransdegradation, res[k0:k1].wavelengthref, allwaves, spline=splval)
            data[i,*]=scandeg
            k0=k1+1 
        endfor
        indata={startScanSD:startscan, wavelength:allwaves, prismTransDeg:data}

    endif

    if n_elements(wavelength) ne 0 then begin
        prismdeg = dblarr(n_elements(indata.startScanSD))
        for i=0L, n_elements(prismdeg)-1L do begin
            prismdeg[i] = interpol(indata.prismtransdeg[i,*], indata.wavelength, wavelength, /spl)
        endfor
        prismdeg={wavelength:wavelength, startScanSD:indata.startScanSD, prismtransDeg:prismdeg}
    endif

    if keyword_set(verbose) then print,'Done'
    return, indata
end

