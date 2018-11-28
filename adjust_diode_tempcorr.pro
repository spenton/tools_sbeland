;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Starting with a known diode temperature correction, tweak the correction 
;   to flatten the irradiance over the specified time range at each wavelength.
;   Use a calibrated version with the diode tempearture correction already applied.
;
; CALLING SEQUENCE:
;   result = ADJUST_DIODE_tempcorr(instrumentModeId, starttime=starttime, stoptime=stoptime)
;
; INPUT PARAMETERS:
;   instrumentModeId -
;      41 for VIS1, 43 for UV or 44 for IR (SimA)
;
; OPTIONAL INPUT PARAMETERS:
;   startTime -
;      The lower time range for which data will be returned.
;      May be specified as julian days, gps microseconds, or
;      SORCE mission day numbers.
;
;   stopTime -
;      The upper time range for which data will be returned
;      May be specified as julian days, gps microseconds, or
;      SORCE mission day numbers.
;   gps - 
;      If set, indicates that input times are to be interpreted as
;      microseconds since Jan 6, 1980 midnight UT.  (default)
;   missionDays - 
;      If set, indicates that input times are to be interpreted as
;      days elapsed since launch, i.e. mission days.  If only the startTime
;      argument contains a non-zero value (stopTime is not supplied), then
;      all data for the specified mission day number will be retrieved.
;   julianDays - 
;      If set, indicates that the input times are to be interpreted as
;      julian day numbers.
;
; RETURNED PARAMETERS:
;   
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
;   Revision: $Id: adjust_diode_tempcorr.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;-
;*******************************************************************
function adjust_diode_tempcorr_amoeba, coeffs
    common adjust_diode_common, timestamp, dettemp, irrad, mytempcorr, myrefTemp

    ; apply a multiplier to the diode temperature correction
    new_irrad = irrad / (1d + mytempcorr*coeffs[0]*(myrefTemp - detTemp))

    c=robust_poly_fit(timestamp, new_irrad,1,/double)
    yfit=poly(timestamp,c)
    return,stddev(yfit - new_irrad)*1d6

end
;*******************************************************************
;
function adjust_diode_tempcorr, instrumentModeId, starttime=startTime, stoptime=stopTime, $
         gps=gps, missionDays=missionDays, julianDays=julianDays, version=version, $
         corrfile=corrfile, alldata=alldata, noplot=noplot

    common adjust_diode_common

    ; convert the time range in mission days
    if n_elements(starttime) eq 0 or n_elements(stoptime) eq 0 then begin
        ;starttime = 1577.7970d
        ;stoptime =  1578.4390d
        ;stoptime =  1578.0842d
        ;starttime=3646d
        ;stoptime=3668d
        starttime=4271d
        stoptime=4470d
        missionDays=1
    endif

    if keyword_set(missionDays) then begin
       ; do nothing here since already in mission days
    endif else if keyword_set(julianDays) then begin
       ; user specified time in julian days
       starttime = jd2sd(starttime)
       stoptime = jd2sd(stoptime)
    endif else begin
       ; user specified time in gps micro seconds
       starttime = gps2sd(starttime/1.d6)
       stoptime = gps2sd(starttime/1.d6)
    endelse

    if instrumentModeId lt 41 or instrumentModeId gt 45 then begin
        print,'Error: invalid instrumentModeId (should be between 41 and 45)'
        return,-1
    endif

    if n_elements(user) eq 0 then     user='sbeland'
    if n_elements(password) eq 0 then password='sbeland1'
    if n_elements(dbdriver) eq 0 then dbDriver='com.sybase.jdbc3.jdbc.SybDriver'
    ;if n_elements(version) eq 0 then version=1005
    ;if n_elements(dburl) eq 0 then    dbUrl='jdbc:sybase:Tds:lasp-db-dev:4100/SORCE_SIM_V20'
    ;if n_elements(corrfile) eq 0 then corrfile='/Users/sbeland/SORCE/data/sima_vis_tempcorr_20.txt'
    if n_elements(version) eq 0 then version=23
    if n_elements(dburl) eq 0 then    dburl='jdbc:sybase:Tds:lasp-db1:4100/SORCE_SIM_V22'
    if n_elements(corrfile) eq 0 then corrfile='/Users/sbeland/SORCE/data/sima_vis_tempcorr_6.txt'

    ; get the data
    if n_elements(alldata) eq 0 then begin
        plans=get_sorce_plan(starttime, stoptime,/mission,sima=sima,simb=simb,activity='SolarQuickScan24')
        alldata=ptrarr(n_elements(plans))

        table1=['Wavelength','SimCalibratedIrradiance','SimConvertedDataNumbers']
        table2=['ModeledDetectorTemperature']
        print,strtrim(string(n_elements(plans)),2)+" plans to process"
        for i=0,n_elements(plans)-1 do begin
            print,i
            spect = get_science_product(table1, plans[i].starttime, plans[i].stoptime, instrumentModeId, /mission, version=version, $
                    user=user, password=password, dburl=dburl, dbdriver=dbdriver) 
            temp = get_science_product(table2, plans[i].starttime, plans[i].stoptime, instrumentModeId, /mission, version=version, $
                    user=user, password=password, dburl=dburl, dbdriver=dbdriver) 
            if size(spect,/tname) ne 'STRUCT' or size(temp,/tname) ne 'STRUCT' then continue

            temperature = interpol(temp.temperature, temp.MICROSECONDSSINCEGPSEPOCH, spect.MICROSECONDSSINCEGPSEPOCH)
            alldata[i] = ptr_new({timestamp:spect.microsecondssincegpsepoch, wavelengthref:spect.wavelengthref, $
                                  irradiance:spect.irradiance, prismposition:spect.prismposition, $
                                  temperature:temperature, version:version})
        endfor
        p=where(ptr_valid(alldata) eq 1,count)
        if count eq 0 then begin
            print,'Error: no valid data within specified time range'
            return,-1
        endif else if count lt n_elements(plans) then alldata=alldata[p]
    endif

    ; get the applied diode temperature corrections
    readcol,corrfile,wref,drdt,tref,format='(d,d,d)'
    myrefTemp = median(tref)

    ; find all the unique prismpositions 
    prismpos=[]
    for i=0L,n_elements(alldata)-1L do prismpos=[prismpos,(*alldata[i]).prismposition]
    s=sort(prismpos)
    q=uniq(prismpos[s])
    prismpos=prismpos[s[q]]
    p=where(prismpos gt 0d,count)
    if count eq 0 then begin
        print,'No valid prismpositions to process'
        return,-1
    endif
    prismpos=prismpos[p]
    mult_fact=dblarr(n_elements(prismpos))
    mult_wave=mult_fact

    print,strtrim(string(n_elements(prismpos)),2)+" wavelengths to process"
    for w=0,n_elements(prismpos)-1 do begin
        timestamp=[]
        irrad=[]
        dettemp=[]
        wave=[]
        for s=0,n_elements(alldata)-1 do begin
            p= where((*alldata[s]).prismposition eq prismpos[w],count)
            if count eq 0 then continue
            wave=[wave,(*alldata[s]).wavelengthref[p]]
            timestamp=[timestamp,(*alldata[s]).timestamp[p]]
            irrad=[irrad,(*alldata[s]).irradiance[p]]
            dettemp=[dettemp,(*alldata[s]).temperature[p]]
        endfor
        if n_elements(wave) lt 8 then continue

        ; prepare to call amoeba to minimize the deviation from a linear trend over the time range
        timestamp=gps2sd(timestamp/1d6)
        mytempcorr = interpol(drdt, wref, wave, /spline)
        ; remove the previous diode temperature correction
        irrad *= (1d + mytempcorr*(myrefTemp - detTemp))

        init_coef=[0d]
        coeffs = AMOEBA(1.0d-8, P0=init_coef, scale=[1d-3], FUNCTION_VALUE=fval, FUNCTION_NAME='adjust_diode_tempcorr_amoeba', nmax=maxiter, ncalls=ncalls)
        mult_fact[w]=median(mytempcorr)+coeffs[0]
        mult_wave[w]=median(wave)
        if NOT keyword_set(noplot) then begin
            irrad0 = irrad / (1d + mytempcorr*(myrefTemp - detTemp))
            irrad1 = irrad / (1d + mytempcorr*coeffs[0]*(myrefTemp - detTemp))
            plot_multi,wave,irrad0,wave,irrad1,/xst,/yst,xtitle='Wavelength',ytitle='Irradiance',$
                psym=[-4,-4],title='Irradiance at '+strtrim(string(median(wave),format='(F0.1)'),2)+' nm',$
                label=['Original Fit','Optimized Fit']
        endif

    endfor

    ; plot the results
    k=where(mult_wave gt 0d,count)
    if count eq 0 then begin
        print,'Not enough data at each wavelength for a fit'
        return,-1
    endif
    mult_fact=mult_fact[k]
    mult_wave=mult_wave[k]
    plot_multi,wref,drdt,mult_wave,mult_fact,xtitle='Wavelength (nm)',ytitle='Diode Temperature Correction',$
        psym=[-3,4],charsize=1.4,/xst,/yst, label=['Original drdt','Optimized drdt']


    return, {wavelength:mult_wave, drdt:mult_fact}
   
end
