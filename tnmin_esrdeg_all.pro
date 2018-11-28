;+
; Author: Stephane Beland
;
; PURPOSE: 
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
; OPTIONAL INPUT PARAMETERS:
;   version -
;      Requests a specific version of the SimProfileIntegralCal, otherwise
;      pick the latest version valid for the specified time range.
;   closest -
;      Skips the interpolation to get the irradiance at requested wavelength
;      but simply uses the value from the closest wavelenegth
;
; RETURNED PARAMETERS:
;
; OPTIONAL OUTPUT PARAMETERS:
;   coeffs -
;      The coefficients of the best polynomial fit.
;   status -
;      Returned STATUS of the mpfitfun routine.
;
; EXAMPLE:  
;
; COMMON BLOCKS:
;      NONE
;
; NOTES:

;
; REVISION HISTORY:
;   Revision: $Id: tnmin_esrdeg_all.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;------------------------------------------------------------------
function mytnmin_func_esrall, x, df, solexp31=solexp31, ird31=ird31, solexp32=solexp32, ird32=ird32, $
    time31=time31, time32=time32, new_ird31=new_ird31, new_ird32=new_ird32
    ; measured_irrad = SolarIrrad * ((exp(-abs(x[0]*solexp)) + exp(-abs(x[1]*solexp))) / 2)
    ; SolarIrrad = measured_irrad / ((exp(-abs(x[0]*solexp)) + exp(-abs(x[1]*solexp))) / 2)
    new_ird31 = ird31 / ((exp(-abs(x[0]*solexp31)) + exp(-abs(x[1]*solexp31))) / 2d)
    new_ird32 = ird32 / ((exp(-abs(x[0]*solexp32)) + exp(-abs(x[1]*solexp32))) / 2d)

    ; since the ESR Table scan data is at different times for A and B
    ; do a 2nd order fit to each iMode and then use the diff at same time
    coeffs31 = robust_poly_fit(time31, new_ird31, 2d, /double)
    coeffs32 = robust_poly_fit(time32, new_ird32, 2d, /double)
    if n_elements(new_ird31) ge n_elements(new_ird32) then begin
        fit31 = poly(time31,coeffs31)
        fit32 = poly(time31,coeffs32)
    endif else begin
        fit31 = poly(time32,coeffs31)
        fit32 = poly(time32,coeffs32)
    endelse
    F = fit32 - fit31
    out_value = abs(mean(F))
    out_value = stddev(F)
    ;print,x[0],x[1],out_value
    return, out_value
end

;*******************************************************************
function mytnmin_func_esrall2, x, df, solexp31=solexp31, ird31=ird31, solexp32=solexp32, ird32=ird32, $
    time31=time31, time32=time32, new_ird31=new_ird31, new_ird32=new_ird32
    new_ird31 = ird31 / exp(-abs(x[0]* solexp31))
    new_ird32 = ird32 / exp(-abs(x[0]* solexp32))

    F = 1.0d - ((ird32[0]/ird32)/(ird31[0]/ird31)) * exp(x[0] * (solexp31[0]-solexp31 - solexp32[0]+solexp32))
    out_value = abs(mean(F))
    ;print,x[0],out_value
    return, out_value
end

;*******************************************************************
function mytnmin_func_esrall3, x, df, solexp31=solexp31, ird31=ird31, solexp32=solexp32, ird32=ird32, $
    time31=time31, time32=time32, new_ird31=new_ird31, new_ird32=new_ird32, dark31=dark31, dark32=dark32, $
    recovery=recovery
    ; measured_irrad = SolarIrrad * exp(-abs(x[1]*solexp))
    ; SolarIrrad = measured_irrad / exp(-abs(x[1]*solexp))
    if n_elements(dark31) eq 0 then begin
        new_ird31 = ird31 / exp(-abs(x[0]*solexp31))
        new_ird32 = ird32 / exp(-abs(x[0]*solexp32))
    endif else begin
        ; account for some recovery
        ; if dark is longer than recovery time - then no degradation left
        temp31 = ((solexp31/(1d -recovery) - dark31)>0.0) < solexp31
        temp32 = ((solexp32/(1d -recovery) - dark32)>0.0) < solexp32
        new_ird31 = ird31 / exp(-abs(x[0] * temp31))
        new_ird32 = ird32 / exp(-abs(x[0] * temp32))
    endelse

    ; since the ESR Table scan data is at different times for A and B
    ; do a 2nd order fit to each iMode and then use the diff at same time
    coeffs31 = robust_poly_fit(time31, new_ird31, 2d, /double)
    coeffs32 = robust_poly_fit(time32, new_ird32, 2d, /double)
    if n_elements(new_ird31) ge n_elements(new_ird32) then begin
        fit31 = poly(time31,coeffs31)
        fit32 = poly(time31,coeffs32)
    endif else begin
        fit31 = poly(time32,coeffs31)
        fit32 = poly(time32,coeffs32)
    endelse
    F = fit32 - fit31
    ;out_value = abs(mean(F))
    out_value = stddev(F)
    ;out_value = MEANABSDEV(F, /median)
    ;out_value = robust_sigma(F, /zero)
    ;out_value = total(abs(f))
    ;print,x[0],x[1],out_value
    return, out_value
end

;*******************************************************************
function mytnmin_func_esrall4, x, df, solexp31=solexp31, ird31=ird31, solexp32=solexp32, ird32=ird32, $
    time31=time31, time32=time32, new_ird31=new_ird31, new_ird32=new_ird32
    ; measured_irrad = SolarIrrad * ((1d + x[0])*exp(x[1]*solexp) + x[0]*exp(x[1]*solexp/2d))
    ; SolarIrrad = measured_irrad / ((1d + x[0])*exp((x[1]*solexp) + x[0]*exp(x[1]*solexp/2d))
    new_ird31 = ird31 / ((1d + x[1])*exp(-abs(x[0]*solexp31)) + x[1]*exp(-abs(x[0]*solexp31/2d)))
    new_ird32 = ird32 / ((1d + x[1])*exp(-abs(x[0]*solexp32)) + x[1]*exp(-abs(x[0]*solexp32/2d)))

    ; since the ESR Table scan data is at different times for A and B
    ; do a 2nd order fit to each iMode and then use the diff at same time
    coeffs31 = robust_poly_fit(time31, new_ird31, 2d, /double)
    coeffs32 = robust_poly_fit(time32, new_ird32, 2d, /double)
    if n_elements(new_ird31) ge n_elements(new_ird32) then begin
        fit31 = poly(time31,coeffs31)
        fit32 = poly(time31,coeffs32)
    endif else begin
        fit31 = poly(time32,coeffs31)
        fit32 = poly(time32,coeffs32)
    endelse
    F = fit32 - fit31
    ;out_value = total(abs(F))
    ;out_value = abs(mean(F))
    ;out_value = stddev(F)
    out_value = MEANABSDEV(F, /median)
    ;print,x,out_value
    return, out_value
end

;*******************************************************************
function mytnmin_func_esrall5, x, df, solexp31=solexp31, ird31=ird31, solexp32=solexp32, ird32=ird32, $
    time31=time31, time32=time32, new_ird31=new_ird31, new_ird32=new_ird32, $
    oneau31=oneau31, oneau32=oneau32
    ; measured_irrad = SolarIrrad * exp(-abs(x[1]*solexp))
    ; SolarIrrad = measured_irrad / exp(-abs(x[1]*solexp))
    if n_elements(x) eq 1 then begin
        new_ird31 = ird31 / exp(-abs(x[0]*solexp31))
        new_ird32 = ird32 / exp(-abs(x[0]*solexp32))
    endif else begin
        new_ird31 = ird31 / (exp(-abs(x[0]*solexp31)) * (1d + au31*x[1])/(x[1]+1d))
        new_ird32 = ird32 / (exp(-abs(x[0]*solexp32)) * (1d + au32*x[1])/(x[1]+1d))
    endelse

    ; since the ESR Table scan data is at different times for A and B
    ; do a 2nd order fit to each iMode and then use the diff at same time
    coeffs31 = robust_poly_fit(time31, new_ird31, 2d, /double)
    coeffs32 = robust_poly_fit(time32, new_ird32, 2d, /double)
    if n_elements(new_ird31) ge n_elements(new_ird32) then begin
        fit31 = poly(time31,coeffs31)
        fit32 = poly(time31,coeffs32)
    endif else begin
        fit31 = poly(time32,coeffs31)
        fit32 = poly(time32,coeffs32)
    endelse
    F = fit32 - fit31
    out_value = abs(mean(F))
    out_value = stddev(F)
    ;print,x[0],x[1],out_value
    return, out_value
end

;*******************************************************************
function mytnmin_func_esrall6, x, df, solexp31=solexp31, ird31=ird31, solexp32=solexp32, ird32=ird32, $
    time31=time31, time32=time32, new_ird31=new_ird31, new_ird32=new_ird32, afact31=afact31, afact32=afact32
    ; measured_irrad = SolarIrrad * exp(-abs(x[1]*solexp))
    ; SolarIrrad = measured_irrad / exp(-abs(x[1]*solexp))
    new_ird31 = ird31 / ((1d - afact31)*exp(-abs(x[0]*solexp31)) + afact31*exp(-abs(x[0]*solexp31)*0.5d))
    new_ird32 = ird32 / ((1d - afact32)*exp(-abs(x[0]*solexp32)) + afact32*exp(-abs(x[0]*solexp32)*0.5d))

    ; since the ESR Table scan data is at different times for A and B
    ; do a 2nd order fit to each iMode and then use the diff at same time
    if n_elements(time31) gt 3 then begin
        coeffs31 = robust_poly_fit(time31, new_ird31, 2d, /double)
    endif else begin
        coeffs31 = ladfit(time31, new_ird31, /double)
    endelse
    if n_elements(time32) gt 3 then begin
        coeffs32 = robust_poly_fit(time32, new_ird32, 2d, /double)
    endif else begin
        coeffs32 = ladfit(time32, new_ird32, /double)
    endelse
    if n_elements(new_ird31) ge n_elements(new_ird32) then begin
        fit31 = poly(time31,coeffs31)
        fit32 = poly(time31,coeffs32)
    endif else begin
        fit31 = poly(time32,coeffs31)
        fit32 = poly(time32,coeffs32)
    endelse
    F = fit32 - fit31
    ;out_value = abs(mean(F))
    out_value = stddev(F)
    ;out_value = MEANABSDEV(F, /median)
    ;out_value = robust_sigma(F, /zero)
    ;out_value = total(abs(f))
    ;print,x[0],x[1],out_value
    return, out_value
end
;*******************************************************************



function tnmin_esrdeg_all, starttime, stoptime, wavelength, version=version, coeffs=coeffs, $
    deltaw=deltaw, status=status, nfev=nfev, niter=niter, errmsg=errmsg, oneau=oneau, $
    dburl=dburl, dbdriver=dbdriver, user=user, password=password, $
    gps=gps, missionDays=missionDays, julianDays=julianDays, inspectrum=inspectrum, $
    solar54=solar54, solar55=solar55, submean=submean, meandiff=meandiff, $
    fit_goodness=fit_goodness, recovery=recovery, au_corr=au_corr, usecoeff=usecoeff, ffunct=ffunct, $
    afact=afact, verbose=verbose, waveb=wavelengthb, prismposa=prismposa, prismposb=prismposb

    ; list of wavelengths from ESRA table scan
    wavea = [280.5d, 283.5, 288, 293.5, 304.5, 319, 333, 342, 355, 369, 375, 385, 395, $
             407, 431, 471, 481, 489, 519, 568, 593, 662, 713, 761, 808, 862, 872, 898, $
             906, 970, 1019, 1070, 1122, 1175, 1285, 1362, 1417, 1482, 1502, 1555, 1597, $
             1627, 1697, 1746, 1822, 1887, 1915, 1950, 1993, 2019, 2111, 2175, 2292, 2402, 2506]

    if keyword_set(missionDays) then begin
       ; user specified time in mission (sorce) days
       t0 = sd2gps(startTime)*1.d6
       t1 = sd2gps(stopTime)*1.d6
    endif else if keyword_set(julianDays) then begin
       ; user specified time in julian days
       t0 = jd2gps(startTime)*1.d6
       t1 = jd2gps(stopTime)*1.d6
    endif else begin
       ; user specified timetags in gps microseconds
       t0 = startTime
       t1 = stopTime
    endelse

    if n_elements(version) eq 0 then version=1003
    if n_elements(user) eq 0 then     user='sbeland'
    if n_elements(password) eq 0 then password='sbeland1'
    if n_elements(dburl) eq 0 then    dbUrl='jdbc:sybase:Tds:lasp-db-dev:4100/SORCE_SIM_V20'
    if n_elements(dbdriver) eq 0 then dbDriver='com.sybase.jdbc3.jdbc.SybDriver'

    ; get the SimSolarExposureData for modes 31 and 32
    if n_elements(solar54) eq 0 then begin
        if keyword_set(verbose) then print,'   getting solar54 ...'
        q1 = 'SELECT orbitStartTimeGps t0, orbitEndTimeGps t1, cumulativeSolarExpHrtOut solar_exp, '+$
             'solarExposureHrtOut as solar_exp_orbit from SimSolarExposureData where instrumentModeId=54'
        query_database, q1, solar54, info
    endif

    if n_elements(solar55) eq 0 then begin
        if keyword_set(verbose) then print,'   getting solar55 ...'
        q1 = 'SELECT orbitStartTimeGps t0, orbitEndTimeGps t1, cumulativeSolarExpHrtOut solar_exp, '+$
             'solarExposureHrtOut as solar_exp_orbit from SimSolarExposureData where instrumentModeId=55'
        query_database, q1, solar55, info
    endif

    if (keyword_set(oneau) or n_elements(ffunct) ge 2) and max(strpos(tag_names(solar54),'ONEAU')) eq -1 then begin
        if keyword_set(verbose) then print,'  getting SUNOBSERVERDISTANCECORRECTION from database ...'
        q1="SELECT microsecondsSinceGpsEpoch, sunObserverDistanceCorrection FROM SolarDistAndDopplerFixedStep2 "+$
           "WHERE instrumentModeId=50 and version=7"
        query_database, q1, solardist, info
        ;time54 = (solar54.(0)+solar54.(1))/2d
        ;time55 = (solar55.(0)+solar55.(1))/2d
        ; the solar exposure time is valid at end of orbit
        time54 = solar54.(1)
        time55 = solar55.(1)
        oneau54 = interpol(solardist.SUNOBSERVERDISTANCECORRECTION, solardist.MICROSECONDSSINCEGPSEPOCH, time54)
        append_tag,solar54,'oneau',oneau54,/slim
        oneau55 = interpol(solardist.SUNOBSERVERDISTANCECORRECTION, solardist.MICROSECONDSSINCEGPSEPOCH, time55)
        append_tag,solar55,'oneau',oneau55,/slim
    endif

    if keyword_set(oneau) or n_elements(ffunct) ge 2 then begin
        ; apply the F-Function and the scaled 1AU to the degradation Column as in insert_degcol.pro
        ; the ffunct is expected to contain the polynomial coefficients 
        ;if n_elements(time54) eq 0 then time54 = (solar54.(0)+solar54.(1))/2d
        ;if n_elements(time55) eq 0 then time55 = (solar55.(0)+solar55.(1))/2d
        if keyword_set(verbose) then print,'   applying ffunct ...'
        if n_elements(time54) eq 0 then time54 = solar54.(1)
        if n_elements(time55) eq 0 then time55 = solar55.(1)
        if n_elements(ffunct) lt 2 then begin
            corr54=1d
            corr55=1d
        endif else begin
            corr54 = poly(gps2sd(time54/1d6), ffunct)
            corr55 = poly(gps2sd(time55/1d6), ffunct)

            ; AS a test (SBeland 2013/07/18), we add the ratio of exposure of SimB/SimA to our F_Function (Version 69)
            ; we then get a new Kappa and raypath
            sum54 = total(solar54.solar_exp_orbit,/cum)
            sum55 = total(solar55.solar_exp_orbit,/cum)
            p=where(sum54 gt 0.0)
            ;corr54[p] += (sum55[p] / sum54[p])  ; Version 69
            ;corr55[p] += (sum55[p] / sum54[p])  ; Version 69
            corr54[p] *= (1d - sum55[p] / sum54[p])  ; Version 71
            corr55[p] *= (1d - sum55[p] / sum54[p])  ; Version 71
 
        endelse
        ; the cumulative solar exposure is in days to go with our Kappa function
        solar54.solar_exp=total(solar54.solar_exp_orbit*corr54,/cum) / (1d +(solar54.oneau-1d)/4d)
        solar55.solar_exp=total(solar55.solar_exp_orbit*corr55,/cum) / (1d +(solar55.oneau-1d)/4d) 

        ; As a different test, scale the solar_exp for SimB with the profile of the ratio B/A (Version 70)
        ;sum55 = interpol(solar55.solar_exp, solar55.t0, solar54.t0)
        ;p=where(solar54.solar_exp gt 0d)
        ;sum55[p] *= ((1d - sum55[p] / sum54[p]) > 0d)
        ;solar55.solar_exp = interpol(sum55, solar54.t0, solar55.t0)
    endif

    if n_elements(recovery) gt 0 then begin
        if keyword_set(verbose) then print,'   applying recovery ...'
        ; a fractional recovery was provided 
        ; the dark times are calculated from the amount of time in the dark since the last 
        ; solar exposure
        ; the darks will correspond to the number of fractional days since last exposure
        ;dark54=dblarr(n_elements(solar54))
        ;dark55=dblarr(n_elements(solar55))
        ;pos=where(solar54.solar_exp_orbit gt 0.0d, count)
        ;p0=0L 
        ;for i=1L,n_elements(solar54)-1L do begin
        ;    dark54[i] = (solar54[i].t0 - solar54[i-1].t1)/1d6
        ;    ; we account for th eportion of the orbit that is not exposed
        ;    ; even if we don't really know wher in the orbit the exposure took place
        ;    ; not perfect but probably better than ignoring it
        ;    if p0 gt 0 then dark54[i] = (solar54[i].t0 - solar54[p0].t0)/1d6 - solar54[p0].solar_exp_orbit
        ;    if solar54[i].solar_exp_orbit gt 0.0 then p0=i
        ;endfor
        ;dark54/=86400d

        ;p0=0L 
        ;for i=1L,n_elements(solar55)-1L do begin
        ;    dark55[i] = (solar55[i].t0 - solar55[i-1].t1)/1d6
        ;    ; we account for th eportion of the orbit that is not exposed
        ;    ; even if we don't really know wher in the orbit the exposure took place
        ;    ; not perfect but probably better than ignoring it
        ;    if p0 gt 0 then dark55[i] = (solar55[i].t0 - solar55[p0].t0)/1d6 - solar55[p0].solar_exp_orbit
        ;    if solar55[i].solar_exp_orbit gt 0.0 then p0=i
        ;endfor
        ;dark55/=86400d

        ; New dark test - cumulative dark time (in days) since first solar exposure
        dark54 = [0d, (solar54[1:-1].t0 - solar54[0:-2].t0)/1d6 - solar54[0:-2].solar_exp_orbit]
        p=where(solar54.solar_exp eq 0)
        dark54[p]=0d 
        dark54=total(dark54,/cum)/86400d
        dark55 = [0d, (solar55[1:-1].t0 - solar55[0:-2].t0)/1d6 - solar55[0:-2].solar_exp_orbit]
        p=where(solar55.solar_exp eq 0)
        dark55[p]=0d 
        dark55=total(dark55,/cum)/86400d
    endif


    if size(inspectrum,/tname) ne 'STRUCT' then begin
        ; get all of the ESR data within the specified timerange and only keep around the requested wavelength
        if keyword_set(verbose) then print,'   getting esrA ...'
        ;esrA=get_science_product(['Wavelength','SimCalibratedIrradiance'],t0,t1,31,/gps,$
        ;esrA=get_science_product(['Wavelength','SimCorrectedIrradiance'],t0,t1,31,/gps,$
        esrA=get_science_product(['Wavelength','SimUncorrectedIrradiance','SimConvertedDataNumbers'],t0,t1,31,/gps,$
            version=version,user=user,password=password,dburl=dburl,dbdriver=dbdriver)
        if keyword_set(verbose) then print,'   getting esrB ...'
        ;esrB=get_science_product(['Wavelength','SimCalibratedIrradiance'],t0,t1,32,/gps,$
        esrB=get_science_product(['Wavelength','SimUncorrectedIrradiance','SimConvertedDataNumbers'],t0,t1,32,/gps,$
            version=version,user=user,password=password,dburl=dburl,dbdriver=dbdriver)
        inspectrum={sp31:esrA, sp32:esrB}
    endif

    if n_elements(prismposa) gt 0 and n_elements(prismposb) gt 0 then begin
        wa = where(inspectrum.sp31.prismposition eq prismposa[0] and $
                   inspectrum.sp31.MICROSECONDSSINCEGPSEPOCH ge t0 and $
                   inspectrum.sp31.MICROSECONDSSINCEGPSEPOCH le t1,count)
    endif else begin
        wa = where(abs(inspectrum.sp31.wavelengthRef - wavelength) le 1d and $
                   inspectrum.sp31.MICROSECONDSSINCEGPSEPOCH ge t0 and $
                   inspectrum.sp31.MICROSECONDSSINCEGPSEPOCH le t1,count)
        if count lt 10 then $
            wa = where(abs(inspectrum.sp31.wavelengthRef - wavelength) le 4d and $
                   inspectrum.sp31.MICROSECONDSSINCEGPSEPOCH ge t0 and $
                   inspectrum.sp31.MICROSECONDSSINCEGPSEPOCH le t1,count)
        if n_elements(deltaw) eq 0 then deltaw = abs(median(inspectrum.sp31[wa].wavelengthRef - wavelength))*4.0d > 0.5d
        ;wa = where(abs(inspectrum.sp31.wavelengthRef - wavelength) le deltaw,count)
        if count eq 0 then begin
            print,'Error: no data from SimA matching the wavelength within '+strtrim(string(deltaw),2)+'nm'
            return,-1
        endif

        ; clean up the data set by removing the points with spikes in irradiance 
        if count gt 3 then begin
            coeffs31=robust_poly_fit(inspectrum.sp31[wa].MICROSECONDSSINCEGPSEPOCH,inspectrum.sp31[wa].irradiance,2,yfit,/double)
            resistant_mean,(inspectrum.sp31[wa].irradiance-yfit),2.0,mean,goodvec=keep0
            if n_elements(keep0) lt n_elements(wa)/2 then $
                resistant_mean,(inspectrum.sp31.irradiance-yfit),3.0,mean,goodvec=keep0
            if n_elements(keep0) lt n_elements(wa)/2 then $
                resistant_mean,(inspectrum.sp31.irradiance-yfit),5.0,mean,goodvec=keep0
            wa=wa[keep0]
        endif
    endelse

    irrad31=[]
    time31=[]
    wave31=[]
    solarexp31=[]
    time54 = [solar54.t0, solar54.t1]
    solarexp54=[0.0,solar54[0:-2].solar_exp, solar54.solar_exp]
    s=sort(time54)
    time54=time54[s]
    solarexp54=solarexp54[s]

    if n_elements(recovery) gt 0 then dark31=[]

    ; look for multiple points within 2 hours (from ESRFullScan)
    if keyword_set(verbose) then print,'   removing duplicates wavelengths for esrA ...'
    yhist=histogram(inspectrum.sp31[wa].(0)/1d6,binsize=2d*3600d,reverse_indices=rev_ind,location=xhist)
    p=where(yhist gt 1.0,count)
    ;if count eq 0 then begin
    if count ge 0 then begin
        ; no duplicate entries within 2 hours
        irrad31=inspectrum.sp31[wa].irradiance
        time31=gps2sd(inspectrum.sp31[wa].MICROSECONDSSINCEGPSEPOCH/1d6)
        wave31=inspectrum.sp31[wa].wavelengthRef
        ; get the solar exposure at these times
        solarexp31 = interpol(solarexp54,time54, inspectrum.sp31[wa].MICROSECONDSSINCEGPSEPOCH)
        solarexp31/=86400d
        if n_elements(recovery) gt 0 then dark31=interpol(dark54, solar54.t1, inspectrum.sp31[wa].MICROSECONDSSINCEGPSEPOCH)
    endif else begin
        ; select one from the duplicates
        p=where(yhist gt 0.0,count)
        for i=0,count-1 do begin
            if yhist[p[i]] le 1.0 then begin
                pos = rev_ind[rev_ind[p[i]] : rev_ind[p[i]+1]-1]
            endif else begin
                pos = rev_ind[rev_ind[p[i]] : rev_ind[p[i]+1]-1]
                ; get the data with the wavelength closest to our ref wavelength
                mn=min(abs(inspectrum.sp31[wa[pos]].wavelengthRef - wavelength),q)
                pos = pos[q]
            endelse
            irrad31=[irrad31,inspectrum.sp31[wa[pos]].irradiance]
            time31=[time31,gps2sd(inspectrum.sp31[wa[pos]].MICROSECONDSSINCEGPSEPOCH/1d6)]
            wave31=[wave31,inspectrum.sp31[wa[pos]].wavelengthRef]
            temp = interpol(solarexp54,time54, inspectrum.sp31[wa[pos]].MICROSECONDSSINCEGPSEPOCH)
            solarexp31=[solarexp31, temp/86400d]
            if n_elements(recovery) gt 0 then $
                dark31=[dark31,interpol(dark54, solar54.t1, inspectrum.sp31[wa[pos]].MICROSECONDSSINCEGPSEPOCH)]
        endfor
    endelse

    ; clean up the data set by removing the points with spikes in irradiance 
    if keyword_set(verbose) then print,'   clean up esrA ...'
    if n_elements(irrad31) gt 3 then begin
        coeffs31=robust_poly_fit(time31,irrad31,2,yfit,/double)
        resistant_mean,(irrad31-yfit),3.0,mean,goodvec=keep0
        if n_elements(keep0) lt n_elements(irrad31)/2 then $
            resistant_mean,(irrad31-yfit),5.0,mean,goodvec=keep0
        irrad31=irrad31[keep0]
        time31=time31[keep0]
        wave31=wave31[keep0]
        solarexp31=solarexp31[keep0]
        if n_elements(recovery) gt 0 then dark31=dark31[keep0]
    endif


    ; find the closest wavelength to match with SimA
    if n_elements(prismposa) gt 0 and n_elements(prismposb) gt 0 then begin
        wb = where(inspectrum.sp32.prismposition eq prismposb[0] and $
                   inspectrum.sp32.MICROSECONDSSINCEGPSEPOCH ge t0 and $
                   inspectrum.sp32.MICROSECONDSSINCEGPSEPOCH le t1,count)
    endif else if n_elements(wavelengthb) eq 0 then begin
        if keyword_set(verbose) then print,'   find matching ESRB points ...'
        posb = where(inspectrum.sp32.MICROSECONDSSINCEGPSEPOCH ge t0 and $
                   inspectrum.sp32.MICROSECONDSSINCEGPSEPOCH le t1,count)
        yhist=histogram(inspectrum.sp32[posb].wavelengthRef, binsize=0.5, reverse_indices=rev_ind,location=xhist)
        p=where(abs(xhist-wavelength) lt 5.0 and yhist ge 10, count)
        if count eq 0 then p=where(abs(xhist-wavelength) lt 12.0 and yhist ge 8, count)
        if count eq 0 then p=where(abs(xhist-wavelength) lt 20.0 and yhist ge 5, count)
        if count eq 0 then p=where(abs(xhist-wavelength) lt 20.0 and yhist ge 2, count)
        ;mn=(min(yhist[p]) * 2.0) >12 
        mn=min(yhist[p],max=mx)
        mn=float(mx)/3.0
        ; adjust width with dispersion
        p=where(abs(xhist-wavelength) lt 1.0 and yhist ge mn, count)
        if count eq 0 or total(yhist[p]) lt mn then p=where(abs(xhist-wavelength) lt 2.0 and yhist ge mn, count)
        if count eq 0 or total(yhist[p]) lt mn then p=where(abs(xhist-wavelength) lt 5.0 and yhist ge mn, count)
        if count eq 0 or total(yhist[p]) lt mn then p=where(abs(xhist-wavelength) lt 10.0 and yhist ge mn, count)
        if count eq 0 or total(yhist[p]) lt mn then p=where(abs(xhist-wavelength) lt 15.0 and yhist ge mn, count)
        if count eq 0 or total(yhist[p]) lt mn then p=where(abs(xhist-wavelength) lt 20.0 and yhist ge mn, count)
        ; get the wavelength from the highest peak within this range
        mx=max(yhist[p],pos)
        wavelength32 = xhist[p[pos]]
        wb = where(abs(inspectrum.sp32[posb].wavelengthRef - wavelength32) le deltaw,count)
        if count lt 3 then wb = where(abs(inspectrum.sp32[posb].wavelengthRef - wavelength32) le deltaw*2d,count)
        if count lt 3 then wb = where(abs(inspectrum.sp32[posb].wavelengthRef - wavelength32) le deltaw*4d,count)
        if count lt 2 then begin
            print,'Error: no data from SimB matching the wavelength within '+strtrim(string(deltaw),2)+'nm'
            return,-1
        endif
        wb=posb[wb]
        medw=median(inspectrum.sp32[wb].wavelengthRef)
    endif else begin
        if keyword_set(verbose) then print,'   find matching ESRB points ...'
        wb = where(abs(inspectrum.sp32.wavelengthRef - wavelengthb) le 1.0 and $
                   inspectrum.sp32.MICROSECONDSSINCEGPSEPOCH ge t0 and $
                   inspectrum.sp32.MICROSECONDSSINCEGPSEPOCH le t1,count)
    endelse
    if n_elements(wb) lt 2 then begin
        print,'Error: no matching data from SimB'
        return,-1
    endif

    irrad32=[]
    time32=[]
    wave32=[]
    solarexp32=[]
    time55 = [solar55.t0, solar55.t1]
    solarexp55=[0.0,solar55[0:-2].solar_exp, solar55.solar_exp]
    s=sort(time55)
    time55=time55[s]
    solarexp55=solarexp55[s]
    if n_elements(recovery) gt 0 then dark32=[]

    if keyword_set(verbose) then print,'   removing duplicates wavelengths for esrB ...'
    yhist=histogram(inspectrum.sp32[wb].(0)/1d6,binsize=2d*3600d,reverse_indices=rev_ind,location=xhist)
    p=where(yhist gt 1.0,count)
    ;if count eq 0 then begin
    if count ge 0 then begin
        ; no duplicate entries within 2 hours
        irrad32=inspectrum.sp32[wb].irradiance
        time32=gps2sd(inspectrum.sp32[wb].MICROSECONDSSINCEGPSEPOCH/1d6)
        wave32=inspectrum.sp32[wb].wavelengthRef
        ; get the solar exposure at these times
        solarexp32 = interpol(solarexp55,time55, inspectrum.sp32[wb].MICROSECONDSSINCEGPSEPOCH)
        solarexp32/=86400d
        if n_elements(recovery) gt 0 then dark32=interpol(dark55, solar55.t1, inspectrum.sp32[wb].MICROSECONDSSINCEGPSEPOCH)
    endif else begin
        ; select one from the duplicates
        p=where(yhist gt 0.0,count)
        for i=0,count-1 do begin
            if yhist[p[i]] le 1.0 then begin
                pos = rev_ind[rev_ind[p[i]] : rev_ind[p[i]+1]-1]
            endif else begin
                pos = rev_ind[rev_ind[p[i]] : rev_ind[p[i]+1]-1]
                ; get the data with the wavelength closest to our ref wavelength
                mn=min(abs(inspectrum.sp32[wb[pos]].wavelengthRef - wavelength),q)
                pos = pos[q]
            endelse
            irrad32=[irrad32,inspectrum.sp32[wb[pos]].irradiance]
            time32=[time32,gps2sd(inspectrum.sp32[wb[pos]].MICROSECONDSSINCEGPSEPOCH/1d6)]
            wave32=[wave32,inspectrum.sp32[wb[pos]].wavelengthRef]
            temp = interpol(solarexp55,time55, inspectrum.sp32[wb[pos]].MICROSECONDSSINCEGPSEPOCH)
            solarexp32=[solarexp32, temp/86400d]
            if n_elements(recovery) gt 0 then $
                dark32=[dark32,interpol(dark55, solar55.t1, inspectrum.sp32[wb[pos]].MICROSECONDSSINCEGPSEPOCH)]
        endfor
    endelse

    ; clean up the data set by removing the points with spikes in irradiance 
    if n_elements(irrad32) gt 3 then begin
        coeffs32=robust_poly_fit(time32,irrad32,2,yfit,/double)
        resistant_mean,(irrad32-yfit),3.0,mean,goodvec=keep0
        if n_elements(keep0) lt n_elements(irrad32)/2 then $
            resistant_mean,(irrad32-yfit),4.0,mean,goodvec=keep0
        if n_elements(keep0) lt n_elements(irrad32)/2 then $
            resistant_mean,(irrad32-yfit),5.0,mean,goodvec=keep0
        irrad32=irrad32[keep0]
        time32=time32[keep0]
        wave32=wave32[keep0]
        solarexp32=solarexp32[keep0]
        if n_elements(recovery) gt 0 then dark32=dark32[keep0]
    endif

    ; check if we want to add a 1 AU correction to the irradiance
    if keyword_set(au_corr) and n_elements(oneau) eq 0 then begin
        if keyword_set(verbose) then print,'   getting the 1AU correction ... '
        ; unfortunately we can't use a single sql query to extract all of the 
        ; data since the timestamps for ESR are averages and don't match the
        ; timestamps in the SolarDistance table
        ; the instrumentId for SORCE is 50 and we're at version=6  (querying for the max(version) is VERY SLOW)
        query1="SELECT TOP 1 * FROM SolarDistAndDopplerFixedStep2 WHERE instrumentModeId=50 and version=6 and "
        for i=0,n_elements(time31)-1 do begin
            query2 = query1+"microsecondsSinceGpsEpoch>="+strtrim(string(ulong64(t0)),2)
            query_database, query2, solarDistDop, info
            sunObserverDistance = interpol(solarDistDop.sunObserverDistanceCorrection, $
                solarDistDop.microsecondsSinceGpsEpoch, outdata.timetag,/spline)
            outdata.irradiance /= sunObserverDistance
        endfor
    endif


    if n_elements(recovery) gt 0 then begin
        ; Should we account for the recovery for only the last
        ; exposure or for for all the previously exposed time?
        ; Can only the last layer evaporate during the long darks
        ; leaving a residual that gets cooked on the surface at 
        ; the next exposure?
        ; From some papers, the recovery is also temperature
        ; dependant, so during the cold OBC, there should have 
        ; been very little recovery.
        ; As a test, assume the recovery has an exponential 
        ; profile in the form  exp(abs(dark/constant)) where
        ; the constant is provided by the user (something around 2d6)
        functargs = {solexp31:solarexp31, ird31:irrad31, solexp32:solarexp32, ird32:irrad32, $
             time31:time31, time32:time32, dark31:dark31, dark32:dark32, recovery:recovery}
    endif else begin
        functargs = {solexp31:solarexp31, ird31:irrad31, solexp32:solarexp32, ird32:irrad32, $
             time31:time31, time32:time32}
    endelse

    parinfo=replicate({value:0.0d, fixed:0, limited:[0,0], limits:[0d,0d], step:0d, tnside:2},1)
    parinfo[0].value=1.0d-4
    ;parinfo[0].fixed=1
    ;parinfo[0].step=1.0d-8
    ;parinfo[1].value=0.0d-3
    ;parinfo[1].fixed=1
    ;parinfo[1].step=1.0d-8
    ; a lower limit of 0 will assure the parameter remains positive
    ;parinfo[1].limited=[1,1]
    ;parinfo[1].limits=[-1d-11,1d-11]
    ;parinfo[1].value=1.0d-10
    ;parinfo[1].value=0d
    ;parinfo[1].fixed=1

    if keyword_set(au_corr) then begin
        ; get the fraction of light extending beyond burned spot
        if keyword_set(verbose) then print,'   applying  the 1AU correction ... '
        closest_dist=0.9832898912d 
        farthest_dist=1.0167103335d
        oneau_corr = farthest_dist - closest_dist
        ;oneau_corr = 1.05d
        ;oneau_corr = (1d + (farthest_dist - closest_dist)) ^2d
        ;oneau_corr = (1d + (farthest_dist - closest_dist))
        oneau31=interpol(oneau.SUNOBSERVERDISTANCECORRECTION, oneau.MICROSECONDSSINCEGPSEPOCH, sd2gps(time31)*1d6)
        oneau32=interpol(oneau.SUNOBSERVERDISTANCECORRECTION, oneau.MICROSECONDSSINCEGPSEPOCH, sd2gps(time32)*1d6)
        parinfo=replicate(parinfo[0],2)
        parinfo[1].value=oneau_corr
        parinfo[1].fixed=1
        parinfo[1].limited=[0,0]
        ;parinfo[1].limited=[1,1]
        ;parinfo[1].limits=[0.05d,0.0d]
        functargs = create_struct(functargs, 'oneau31', oneau31, 'oneau32', oneau32)
    endif

    if n_elements(afact) eq 0 then begin
        afact31=0d
        afact32=0d
    endif else begin
        afact31=interpol(afact.a, afact.wavelength, wave31, /spline)
        afact32=interpol(afact.a, afact.wavelength, wave32, /spline)
    endelse
    functargs = create_struct(functargs, 'afact31', afact31, 'afact32', afact32)

    ; now we want to minimize the difference in irradiance between 31 and 32
    if n_elements(usecoeff) eq 0 then begin
        if keyword_set(verbose) then print,'  initial guess: ',parinfo.value
        coeffs = tnmin('mytnmin_func_esrall6', functargs=functargs, bestmin=f0, status=status, $
            nfev=nfev, niter=niter, autoderivative=1, errmsg=errmsg, parinfo=parinfo, /quiet)
    endif else begin
        coeffs=usecoeff
    endelse

    ;fit_goodness=mytnmin_func_esrall3(coeffs, df, solexp31=solarexp31, ird31=irrad31, solexp32=solarexp32, $
    ;    ird32=irrad32, time31=time31, time32=time32, new_ird31=new_ird31, new_ird32=new_ird32, $
    ;    dark31=dark31, dark32=dark32, recovery=recovery); ,oneau31=oneau31, oneau32=oneau32)
    if keyword_set(verbose) then print,' minimizing ',median(wave31),' ...'
    fit_goodness=mytnmin_func_esrall6(coeffs, df, solexp31=solarexp31, ird31=irrad31, solexp32=solarexp32, $
        ird32=irrad32, time31=time31, time32=time32, new_ird31=new_ird31, new_ird32=new_ird32, $
        afact31=afact31, afact32=afact32)

    if keyword_set(submean) then begin
        ; subtract the difference of the mean between esra and esrb to esrb and re-process
        toto=inspectrum
        meandiff = mean(new_ird32) - mean(new_ird31)
        toto.sp31.irradiance+=meandiff
        result=tnmin_esrdeg_all(starttime, stoptime, wavelength, version=version, coeffs=coeffs, $
            deltaw=deltaw, status=status, nfev=nfev, niter=niter, errmsg=errmsg, $
            dburl=dburl, dbdriver=dbdriver, user=user, password=password, $
            gps=gps, missionDays=missionDays, julianDays=julianDays, inspectrum=toto, oneau=oneau, $
            solar54=solar54, solar55=solar55, meandiff=meandiff, fit_goodness=fitgoodness, $
            au_corr=au_corr, usecoeff=usecoeff, recovery=recovery, afact=afact)
        new_ird31=result.sp31.new_irrad
        new_ird32=result.sp32.new_irrad

        ;mod_irrad31 = irrad31 + meandiff
        ;fit_goodness=mytnmin_func_esrall3(coeffs, df, solexp31=solarexp31, ird31=mod_irrad31, solexp32=solarexp32, $
        ;    ird32=irrad32, time31=time31, time32=time32, new_ird31=new_ird31, new_ird32=new_ird32, $
        ;    dark31=dark31, dark32=dark32, recovery=recovery)
        if keyword_set(verbose) then print,'  difference of means: ',meandiff
    endif
    if keyword_set(verbose) then print,'  final value: ',coeffs
    if keyword_set(verbose) then print,'  fit goodness: ', fit_goodness

    if keyword_set(verbose) then begin
        ; plot the fit in vebose mode
        plot_multi,time31,irrad31, time32,irrad32, time31,new_ird31, time32,new_ird32, /xst,/yst,$
            xtitle='Mission Day',ytitle='Uncorrected Irradiance',psym=[-4,-4,-5,-5], $
            title='ESRA vs ESRB @ '+strtrim(string(median(wave31),format='(F0.1)'),2), $
            label=['ESRA','ESRB','Corrected ESRA','Corrected ESRB']
    endif

    sp31={MICROSECONDSSINCEGPSEPOCH:sd2gps(time31)*1d6, wavelength:wave31, solarexp:solarexp31, $
        irradiance:irrad31, new_irrad:new_ird31}
    sp32={MICROSECONDSSINCEGPSEPOCH:sd2gps(time32)*1d6, wavelength:wave32, solarexp:solarexp32, $
        irradiance:irrad32, new_irrad:new_ird32}
    return, {sp31:sp31, sp32:sp32}

end
