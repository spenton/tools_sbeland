;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Extract the degradation information from the UV diode matching
;   the data for the requested wavelength taken at the "same" time for 
;   SIMA and SIMB.  This is done for two consecutive SolarQuickScan24 of 
;   SIMBESRB to get a running slope of the degradation with mission time.
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
;    The table used for the irradiance is expected to have no degradation
;    correction applied to it (as in version 15 of our development database for version 19).
;
;
; REVISION HISTORY:
;   Revision: $Id: get_diode_fdeg.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;------------------------------------------------------------------
function get_diode_fdeg, starttime, stoptime, instrumentModeId, wavelength, version=version, $
    bin55=bin55, step55=step55, $
    gps=gps, missionDays=missionDays, julianDays=julianDays, $
    dburl=dburl, dbdriver=dbdriver, user=user, password=password, $
    inspectrum=inspectrum, alignobc=alignobc, smooth=smooth, $
    noplot=noplot, solar54=solar54, solar55=solar55, kappa=kappa

    obctimes = sd2gps([0.0d, 441.5, 1570.0d, 2173d, 2455d, 2803d, 5000d])*1d6

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

    if n_elements(version) eq 0 then version=15
    if n_elements(user) eq 0 then     user='sbeland'
    if n_elements(password) eq 0 then password='sbeland1'
    if n_elements(dburl) eq 0 then    dbUrl='jdbc:sybase:Tds:lasp-db-dev:4100/SORCE_SIM_V19'
    if n_elements(dbdriver) eq 0 then dbDriver='com.sybase.jdbc3.jdbc.SybDriver'
    if n_elements(bin55) eq 0 then bin55=1

    ; limit the wavelength range according to the specified instrumentModeId
    if instrumentModeId eq 41 or instrumentModeId eq 45 then begin
        instrument=54
        instModeIdA = 41
        instModeIdB = 45
    endif else if instrumentModeId eq 43 or instrumentModeId eq 47 then begin
        instrument=54
        instModeIdA = 43
        instModeIdB = 47
    endif else begin
        print,'Error: wrong instrumentModeId specified'
        return,-1
    endelse

    nsum=3

    ; get the SimSolarExposureData for modes 31 and 32
    if n_elements(solar54) eq 0 then begin
        q1 = 'SELECT orbitStartTimeGps t0, orbitEndTimeGps t1, cumulativeSolarExpHrtOut solar_exp, '+$
             'solarExposureHrtOut as solar_exp_orbit from SimSolarExposureData where instrumentModeId=54'
        query_database, q1, solar54, info
    endif

    if n_elements(solar55) eq 0 then begin
        q1 = 'SELECT orbitStartTimeGps t0, orbitEndTimeGps t1, cumulativeSolarExpHrtOut solar_exp, '+$
             'solarExposureHrtOut as solar_exp_orbit from SimSolarExposureData where instrumentModeId=55'
        query_database, q1, solar55, info
    endif

    ; the order here is important for the query_database call later
    irradColumns=['Wavelength','SimCorrectedIrradiance']

    query_database,/reset

    ; extract the spectrum for SimA and SimB if they are not provided
    if size(inspectrum,/tname) ne 'STRUCT' then begin
        plansA = get_sorce_plan(t0, t1, /gps,instrument='SIMA',activity='SolarQuickScan24')
        plansB = get_sorce_plan(t0, t1, /gps,instrument='SIMB',activity='SolarQuickScan24')
        match,plansA.starttime, plansB.starttime,suba,subb,count=nspect
        if nspect lt 4 then begin
            print,'Error: Not enough matching SolarQuickScan24 for SimA and SimB'
            return,-1
        endif
        plansA=plansA[suba]
        plansB=plansB[subb]
        diodea_data=ptrarr(nspect)
        diodeb_data=ptrarr(nspect)
        diodea_starttime=dblarr(nspect)
        diodeb_starttime=dblarr(nspect)

        print,'extracting '+strtrim(string(nspect),2)+' x2 spectrum ...'
        for pl=0, nspect-1L do begin
            gps0 = plansA[pl].starttime*1d6
            gps1 = plansA[pl].stoptime*1d6
            q1='SELECT w.microsecondsSinceGpsEpoch,w.instrumentModeId,w.version,w.wavelengthRef,ird.irradiance '+$
               'FROM '+irradColumns[0]+' w, '+ irradColumns[1]+' ird  where '+$
               'w.version='+strtrim(string(version),2)+' and w.version=ird.version and '+$
               'w.instrumentModeId='+strtrim(string(instModeIdA),2)+' and '+$ 
               'w.instrumentModeId=ird.instrumentModeId and '+$
               'w.microsecondsSinceGpsEpoch>='+strtrim(string(ulong64(gps0)),2)+' and '+$
               'w.microsecondsSinceGpsEpoch<='+strtrim(string(ulong64(gps1)),2)+' and '+$
               'w.microsecondsSinceGpsEpoch=ird.microsecondsSinceGpsEpoch'
            query_database, q1, resA, user=user,password=password,dburl=dburl,dbdriver=dbdriver
            if size(resA,/tname) ne 'STRUCT' then begin
                print,'Error: no spectra found for this plan ',plansA[pl]
                continue
            endif

            gps0 = plansB[pl].starttime*1d6
            gps1 = plansB[pl].stoptime*1d6
            q1='SELECT w.microsecondsSinceGpsEpoch,w.instrumentModeId,w.version,w.wavelengthRef,ird.irradiance '+$
               'FROM '+irradColumns[0]+' w, '+ irradColumns[1]+' ird  where '+$
               'w.version='+strtrim(string(version),2)+' and w.version=ird.version and '+$
               'w.instrumentModeId='+strtrim(string(instModeIdB),2)+' and '+$ 
               'w.instrumentModeId=ird.instrumentModeId and '+$
               'w.microsecondsSinceGpsEpoch>='+strtrim(string(ulong64(gps0)),2)+' and '+$
               'w.microsecondsSinceGpsEpoch<='+strtrim(string(ulong64(gps1)),2)+' and '+$
               'w.microsecondsSinceGpsEpoch=ird.microsecondsSinceGpsEpoch'
            query_database, q1, resB, user=user,password=password,dburl=dburl,dbdriver=dbdriver
            if size(resB,/tname) ne 'STRUCT' then begin
                print,'Error: no spectra found for this plan ',plansB[pl]
                continue
            endif

            ; sort in wavelength but keep same starttime (A & B scan in reverse wavelength direction)
            diodea_starttime[pl]= plansA[pl].starttime*1d6
            diodeb_starttime[pl]= plansB[pl].starttime*1d6
            s=sort(resA.wavelengthref)
            resA= resA[s]
            s=sort(resB.wavelengthref)
            resB= resB[s]
            diodea_data[pl]= ptr_new(resA)
            diodeb_data[pl]= ptr_new(resB)
        endfor

        keep = where(ptr_valid(diodea_data) and ptr_valid(diodeb_data),count)
        if count eq 0 then begin
            print,'Error:  no valid diode data found'
            return,-1
        endif
        simA = {starttime:diodea_starttime[keep], diode:diodea_data[keep]}
        simB = {starttime:diodeb_starttime[keep], diode:diodeb_data[keep]}
        inspectrum = {simA:temporary(simA), simB:temporary(simB)}

    endif

    ; get the corresponding irradiance at the specific wavelength for each spectrum
    nspect=n_elements(inspectrum.simA.starttime)
    solarexpA=dblarr(nspect)
    solarexpB=dblarr(nspect)
    irradA=dblarr(nspect)
    irradB=dblarr(nspect)
    timeA=dblarr(nspect)
    timeB=dblarr(nspect)
    for j=0L,nspect-1L do begin
        irradA[j]=interpol((*inspectrum.simA.diode[j]).irradiance,(*inspectrum.simA.diode[j]).wavelengthref,wavelength,/spline)
        val=min(abs((*inspectrum.simA.diode[j]).wavelengthref - wavelength),pos)
        timeA[j]=((*inspectrum.simA.diode[j]).(0))[pos]
        q=where(solar54.t1 le inspectrum.simA.starttime[j],count)
        solarexpA[j]=solar54[q[-1]].solar_exp/86400d
        irradB[j]=interpol((*inspectrum.simB.diode[j]).irradiance,(*inspectrum.simB.diode[j]).wavelengthref,wavelength,/spline)
        val=min(abs((*inspectrum.simB.diode[j]).wavelengthref - wavelength),pos)
        timeB[j]=((*inspectrum.simB.diode[j]).(0))[pos]
        q=where(solar55.t1 le inspectrum.simB.starttime[j],count)
        solarexpB[j]=solar55[q[-1]].solar_exp/86400d
    endfor

    ; if alignobc was specified, we align each segment by moving the latter piece
    ; in irradiance to match the average extrapolated points from both segments.
    if keyword_set(alignobc) then begin
        ; process SIMA 
        for i=0,n_elements(obctimes)-3 do begin
            pleft = where(timeA ge obctimes[i] and timeA lt obctimes[i+1],count)
            if count eq 0 then continue
            pright = where(timeA ge obctimes[i+1] and timeA lt obctimes[i+2],count)
            if count eq 0 then continue
            cleft = robust_poly_fit(timeA[pleft],irradA[pleft],2,/double)
            cright = robust_poly_fit(timeA[pright],irradA[pright],2,/double)
            ;p0=n_elements(pleft)-21 > 0
            ;cleft=ladfit(timeA[pleft[p0:-2]], irradA[pleft[p0:-2]])
            ;p1=n_elements(pright)-1 <20 
            ;cright=ladfit(timeA[pright[1:p1]], irradA[pright[1:p1]])
            v0=poly(timeA[pright[0]],cleft)
            v1=poly(timeA[pright[0]],cright)
            irradA[pright] += (v0-v1)
        endfor
        ; process SIMB
        for i=0,n_elements(obctimes)-3 do begin
            pleft = where(timeB ge obctimes[i] and timeB lt obctimes[i+1],count)
            if count eq 0 then continue
            pright = where(timeB ge obctimes[i+1] and timeB lt obctimes[i+2],count)
            if count eq 0 then continue
            cleft = robust_poly_fit(timeB[pleft],irradB[pleft],2,/double)
            cright = robust_poly_fit(timeB[pright],irradB[pright],2,/double)
            ;p0=n_elements(pleft)-11 > 0
            ;cleft=ladfit(timeB[pleft[p0:-2]], irradB[pleft[p0:-2]])
            ;p1=n_elements(pright)-1 < 8
            ;cright=ladfit(timeB[pright[1:p1]], irradB[pright[1:p1]])
            v0=poly(timeB[pright[0]],cleft)
            v1=poly(timeB[pright[0]],cright)
            irradB[pright] += (v0-v1)
        endfor
    endif


    ; if smooth was flagged, apply a bspline smoothing (from SDSS idlutils/bspline_iterfit.pro)
    timeA=gps2sd(timeA/1d6)
    timeB=gps2sd(timeB/1d6)
    if keyword_set(smooth) then begin
        ; fit a 2nd order curve to data and remove the outliers
        ; UVA
        coeffA=robust_poly_fit(timeA,irradA,2,yfit,/double)
        resistant_mean,(irradA-yfit),3.0,mean,goodvec=keep0
        sset=bspline_iterfit(timeA[keep0],irradA[keep0],maxiter=0,requiren=10,bkspace=5)
        yfitA=bspline_valu(timeA,sset)
        ; UVB
        coeffB=robust_poly_fit(timeB,irradB,2,yfit,/double)
        resistant_mean,(irradB-yfit),3.0,mean,goodvec=keep0
        sset=bspline_iterfit(timeB[keep0],irradB[keep0],maxiter=0,requiren=10,bkspace=5)
        yfitB=bspline_valu(timeB,sset)
    endif else begin
        yfitA=-1
        yfitB=-1
    endelse


    nfact=nspect-bin55
    f_factor=dblarr(nfact) -1d
    time_f=dblarr(nfact) -1d

    for j=0L,nfact-1L do begin
        if n_elements(step55) eq 0 then begin
            pos0=lindgen(bin55)+j
            pos1=pos0+1
        endif else begin
            pos0=lindgen(bin55)+j
            p=where(inspectrum.simA.starttime ge inspectrum.simA.starttime[j]+step55*1d6*86400d,count)
            if count eq 0 then continue
            pos1=lindgen(bin55)+p[0]
        endelse
        ;f_factor1[j-1] = alog(irradB[j]/irradB[j-1]) - alog(irradA[p1]/irradA[p0])
        ;f_factor1[j-1] /= ((solarexpB[j-1]-solarexpB[j]) - (solarexpA[p0] - solarexpA[p1]))
        if keyword_set(smooth) then begin
            f_factor[j] = alog(mean(yfitB[pos1])/mean(yfitB[pos0])) - alog(mean(yfitA[pos1])/mean(yfitA[pos0]))
            f_factor[j] /= ((mean(solarexpB[pos0])-mean(solarexpB[pos1])) - (mean(solarexpA[pos0]) - mean(solarexpA[pos1])))
        endif else begin
            f_factor[j] = alog(mean(irradB[pos1])/mean(irradB[pos0])) - alog(mean(irradA[pos1])/mean(irradA[pos0]))
            f_factor[j] /= ((mean(solarexpB[pos0])-mean(solarexpB[pos1])) - (mean(solarexpA[pos0]) - mean(solarexpA[pos1])))
        endelse
        time_f[j] = mean(timeB[pos0])
    endfor

    p=where(f_factor ne -1d and time_f ne -1d,count)
    if count eq 0 then begin
        print,'No matching UVA and UVB data '
        return,-1
    endif
    f_factor=f_factor[p]
    time_f=time_f[p]

    if size(kappa,/tname) ne 'STRUCT' then begin
        ; get the kappa
        query_database, /reset
        q2 = "SELECT calibrationSetId from CalibrationMetadata where calibrationTableName='SimPrismDegKappaCal' "+$
            "and version=19 and instrumentModeId="+strtrim(string(instrument),2)
        query_database, q2, result, info 
        q3='select x,y from SimPrismDegKappaCal where calibrationSetId='+strtrim(string(result.(0)),2)
        query_database, q3, kappa, info
        kappa.x = 10d ^ kappa.x
    endif
    kappadeg =interpol(kappa.y, kappa.x, wavelength, /spline)

    if not keyword_set(noplot) then begin
        resistant_mean,f_factor,3.0,mean,goodvec=keep0
        ;keep0=lindgen(n_elements(f_factor))
        label='F_Factor @ '+strtrim(string(wavelength,format='(F0.2)'),2)+'nm  '
        label=label+'Kappa='+strtrim(string(mean(kappadeg),format='(E0.4)'),2)
        ;if n_elements(step55) gt 0 then label=label+'   Step='+strtrim(string(step55,format='(F0.1)'),2)+' days'
        ;plot_multi,time_f[keep0],f_factor[keep0]/kappa[keep0],psym=-4,title=label, xtitle='Mission Day',ytitle='F/Kappa',$
        ;    yrange=[0,2],xrange=[0,4000], charsize=1.4
        if n_elements(step55) gt 0 then title='Step='+strtrim(string(step55,format='(F0.1)'),2)+' days' else title=''
        title += '  '+strtrim(string(wavelength,format='(F0.2)'),2)+'nm  InstMode='+strtrim(string(instrumentModeId),2)
        lineplot,time_f[keep0],f_factor[keep0]/kappadeg[0],psym=-3,ptitle='F_Factor', xtitle='Mission Day',ytitle='F/Kappa',$
            yrange=[0,2],xrange=[0,4000], charsize=1.4,title=title
    endif

    return, {timeA:timeA, timeB:timeB, irradA:irradA, irradB:irradB, solarexpA:solarexpA, solarexpB:solarexpB, $
             time_f:time_f, f_factor:f_factor, f_factor_scaled:f_factor/kappadeg[0]}

end

