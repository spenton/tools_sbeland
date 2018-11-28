;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Using an existing version of the diode data, modify the rayPath factor 
;   and apply a new one to the time series for SimA and SimB at some
;   wavelengths.  The new rayPath value at each wavelength is 
;   ajusted so that the differences in the new "calibratedIrradiance" 
;   between SimA and SimB is constant for the specified time range.
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
;      Requests a specific version of the SimCalibratedIrradiance otherwise
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
;    We use the SimUncorrectedIrradiance as a starting point and apply the 1AU,
;    the degradationColumn, and the kappa value.
;
;
; REVISION HISTORY:
;   Revision: $Id: adjust_pd_raypath.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;------------------------------------------------------------------
function func_getpdslope, x, df, sd=sd, kappa_A=kappa_A, kappa_b=kappa_b, degcol_a=degcol_a, degcol_b=degcol_b, $
    irrad_a=irrad_a, irrad_b=irrad_b, new_irrad_a=new_irrad_a, new_irrad_b=new_irrad_b, delta_irrad=delta_irrad

    ; quiet down the warnings in poly_fit
    prismTransDegA = exp(-kappa_a[0] * degcol_a) - x[0] * exp(-kappa_a[0] * degcol_a) + x[0] * exp(-kappa_a[0] * degcol_a / 2d)
    prismTransDegB = exp(-kappa_b[0] * degcol_b) - x[0] * exp(-kappa_b[0] * degcol_b) + x[0] * exp(-kappa_b[0] * degcol_b / 2d)
    detectorDeg = 1.0d ;  - x[1] * exp(x[2] * sd)
    ;detectorDeg = x[1] + x[2]*sd
    new_irrad_a = irrad_a /  (prismTransDegA * detectorDeg)
    new_irrad_b = irrad_b /  (prismTransDegB * detectorDeg)
    delta_irrad = new_irrad_a - new_irrad_b

    ; fit a line 
    ;coeff = robust_poly_fit(sd, delta_irrad, 1, /double)
    coeff = ladfit(sd, delta_irrad, /double)

    ; we want to minimize the slope
    ;print,coeff[1]
    return, abs(coeff[1]*1d20)

end
;*******************************************************************

function adjust_pd_raypath, starttime, stoptime, instrumentModeId, version=version, $
    dburl=dburl, dbdriver=dbdriver, user=user, password=password, $
    gps=gps, missionDays=missionDays, julianDays=julianDays, $
    solardist=solardist, waverange=waverange, plotpd=plotpd, noplot=noplot, $
    inspectrum=inspectpd, degcolA=degcolA, degcolB=degColb, kappa=kappa, wavestep=wavestep, $
    order=order


    modes = [41,43,44]
    match,modes,instrumentModeId,suba,subb,count=count
    if count eq 0 then begin
        print,'Error: wrong instrumentModeId (expecting '+modes,+')'
        return,-1
    endif
    ; define the modes for SimA and SimB
    modes=[modes[suba],modes[suba]+4]

    obctimes = gps2sd([0.0d, 441.5, 1570.0d, 2173d, 2455d, 2803d, 5000d])*1d6

    if keyword_set(gps) then begin
       ; user specified time in gps
       t0 = gps2sd(startTime/1.d6)
       t1 = gps2sd(stopTime/1.d6)
    endif else if keyword_set(julianDays) then begin
       ; user specified time in julian days
       t0 = jd2sd(startTime)
       t1 = jd2gsd(stopTime)
    endif else begin
       ; user specified timetags in mission days (default)
       t0 = startTime
       t1 = stopTime
    endelse

    if n_elements(version) eq 0 then version=16
    if n_elements(user) eq 0 then     user='sbeland'
    if n_elements(password) eq 0 then password='sbeland1'
    if n_elements(dburl) eq 0 then    dbUrl='jdbc:sybase:Tds:lasp-db-dev:4100/SORCE_SIM_V19'
    if n_elements(dbdriver) eq 0 then dbDriver='com.sybase.jdbc3.jdbc.SybDriver'

    if modes[0] eq 41 then begin
        if n_elements(kappa) eq 0 then begin
            query_database, 'SELECT * FROM SimPrismDegKappaCal WHERE calibrationSetId=1053', kappa
            kappa.x=10d ^ kappa.x
        endif
    endif else if modes[0] eq 43 then begin
        if n_elements(kappa) eq 0 then begin
            query_database, 'SELECT * FROM SimPrismDegKappaCal WHERE calibrationSetId=1054', kappa
            kappa.x=10d ^ kappa.x
        endif
    endif else if modes[0] eq 44 then begin
        if n_elements(kappa) eq 0 then begin
            query_database, 'SELECT * FROM SimPrismDegKappaCal WHERE calibrationSetId=1055', kappa
            kappa.x=10d ^ kappa.x
        endif
        wave0=860d
        wave1=1640d
    endif

    if n_elements(degColA) eq 0 then begin
        ;query_database, 'SELECT * FROM SimPrismDegColumnCalTable WHERE calibrationSetId=1035', degColA
        ; testing new degCol for version 20
        query_database, 'SELECT * FROM SimPrismDegColumnCalTable WHERE calibrationSetId=1123', degColA
    endif
    if n_elements(degColB) eq 0 then begin
        ;query_database, 'SELECT * FROM SimPrismDegColumnCalTable WHERE calibrationSetId=1036', degColB
        ; testing new degCol for version 20
        query_database, 'SELECT * FROM SimPrismDegColumnCalTable WHERE calibrationSetId=1124', degColB
    endif

    if n_elements(solardist) eq 0 then begin
        print,'  getting SUNOBSERVERDISTANCECORRECTION from database ...'
        q1="SELECT microsecondsSinceGpsEpoch, sunObserverDistanceCorrection FROM SolarDistAndDopplerFixedStep2 "+$
           "WHERE instrumentModeId=50 and version=7"
        query_database, q1, solardist, info
    endif

    ; get the spectrum and solar exposure data for a single wavelength 
    ; and use the first spectra to determine the wavelengths array
    if modes[0] eq 41 then awave=400d else if modes[0] eq 43 then awave=285d else if modes[0] eq 44 then awave=1000d
    res=tnmin_pddeg(t0,t1,modes,awave,coeffs=c,inspectrum=inspectpd, solar54=solar54, solar55=solar55, $
                    irrad_table='SimUncorrectedIrradiance',version=version)

    waves = (*inspectpd.sima.diode[-1]).wavelengthref
    minwave=min(waves,max=maxwave)
    mn_wb = min((*inspectpd.simb.diode[-1]).wavelengthref, max=mx_wb)

    ; process every other wavelength
    if n_elements(waverange) eq 2 then begin
        p=where(waves ge min(waverange) and waves le max(waverange),count)
        if count eq 0 then begin
            print,'Error: no data within specified waverange'
            return,-1
        endif
        waves=waves[p]
        if count gt 20 then begin
            if n_elements(wavestep) eq 0 then wavestep=2
            p=lindgen(count/wavestep)*wavestep+1
            waves=waves[p]
        endif
    endif else if n_elements(waverange) eq 1 then begin
        waves=waverange[0]
    endif else begin
        if n_elements(wavestep) eq 0 then wavestep=4
        p=lindgen(n_elements(waves)/wavestep)*wavestep+1
        waves=waves[p]
    endelse

    nwaves=n_elements(waves)

    ; process one wavelength at a time
    parinfo=replicate({value:0.0d, fixed:0, limited:[0,0], limits:[0d,0d], step:1d-3, tnside:2},1)
    parinfo[0].value=0.3d

    a_factor=replicate(1000d,nwaves)
    nspect=n_elements(inspectpd.sima.diode)

    ; for each ESR point wavelength
    print,' processing '+strtrim(string(n_elements(waves)),2)+' wavelengths ...'
    for w=0,nwaves-1 do begin
        print,w+1,nwaves,waves[w],format='(I0," of ",I0,"  (",F0.2,")")'
        if waves[w] lt mn_wb or waves[w] gt mx_wb then continue

        ; get the irradiance at this wavelength for each spectra (form our time series)
        irrad_a=dblarr(nspect)
        timetag_a=irrad_a
        irrad_b=dblarr(nspect)
        timetag_b=irrad_b
        for i=0L, nspect-1 do begin
            irrad_a[i] = interpol((*inspectpd.sima.diode[i]).irradiance,(*inspectpd.sima.diode[i]).wavelengthref, $
                waves[w], /spline)
            mn=min(abs((*inspectpd.sima.diode[i]).wavelengthref - waves[w]),pos)
            timetag_a[i] = ((*inspectpd.sima.diode[i])[pos]).(0)

            irrad_b[i] = interpol((*inspectpd.simb.diode[i]).irradiance,(*inspectpd.simb.diode[i]).wavelengthref, $
                waves[w], /spline)
            mn=min(abs((*inspectpd.simb.diode[i]).wavelengthref - waves[w]),pos)
            timetag_b[i] = ((*inspectpd.simb.diode[i])[pos]).(0)
        endfor

        ; recalculate the prismTransmission Degradation that was used
        degCol_pa = interpol(degcolA.y, degcolA.x, timetag_a/1d6)
        degCol_pb = interpol(degcolB.y, degcolB.x, timetag_b/1d6)
        kappa_pa = interpol(kappa.y, kappa.x, waves[w], /spline)
        kappa_pb = interpol(kappa.y, kappa.x, waves[w], /spline)
        soldist_a = interpol(solardist.SUNOBSERVERDISTANCECORRECTION, solardist.(0), timetag_a)
        soldist_b = interpol(solardist.SUNOBSERVERDISTANCECORRECTION, solardist.(0), timetag_b)

        irrad_a /= soldist_a
        irrad_b /= soldist_b

        ; find the best multiplation factor for the rayPath at this wavelength
        functargs = {sd:timetag_a, kappa_a:kappa_pa, kappa_b:kappa_pb, $
                     degcol_a:degcol_pa, degcol_b:degcol_pb, irrad_a:irrad_a, irrad_b:irrad_b}
        coeffs = tnmin('func_getpdslope', functargs=functargs, bestmin=f0, status=status, $
            nfev=nfev, niter=niter, autoderivative=1, errmsg=errmsg, parinfo=parinfo)
        a_factor[w]=coeffs[0]

        ; get the delta_irrad from linear fit after correction
        res=func_getpdslope(coeffs, df, sd=timetag_a, kappa_a=kappa_pa, kappa_b=kappa_pb, degcol_a=degcol_pa, $
            degcol_b=degcol_pb, irrad_a=irrad_a, irrad_b=irrad_b, new_irrad_A=new_irradA, new_irrad_B=new_irradB, $
            delta_irrad=delta_irrad)

        if not keyword_set(noplot) then begin
            resistant_mean,delta_irrad,5.0,mean,goodvec=k
            title='SimA vs SimB'+' @ '+strtrim(string(waves[w],format='(F10.1)'),2)+' raypath='+strtrim(string(a_factor[w]),2)
            lineplot,gps2sd(timetag_a[k]/1d6),delta_irrad[k],title=title,xtitle='Mission Day',ytitle='Delta Irradiance',$
                ptitle='Differences in Irradiance between modes '+strtrim(string(modes[0]),2)+' and '+strtrim(string(modes[1]),2)
            if keyword_set(plotpd) then begin
                ; plot both SIMA and SIMB irradiance
                title='SimA '+' @ '+strtrim(string(waves[w],format='(F10.1)'),2)+' raypath='+strtrim(string(a_factor[w]),2)
                lineplot,gps2sd(timetag_a[k]/1d6), new_irrada, title=title
                title='SimB '+' @ '+strtrim(string(waves[w],format='(F10.1)'),2)+' raypath='+strtrim(string(a_factor[w]),2)
                lineplot,gps2sd(timetag_b[k]/1d6), new_irradb, title=title
            endif
        endif

    endfor  ; loop for each wavelength
    p=where(abs(a_factor) lt 100,count)
    if count gt 0 then begin
        a_factor=a_factor[p] 
        waves=waves[p]
    endif

    ; clean up the data to prepare fot the bspline fit
    ;coeff=ladfit(waves, a_factor,/double)
    if n_elements(wave0) gt 0 then p=where(abs(a_factor) lt 20 and waves ge wave0 and waves le wave1) else p=where(abs(a_factor) lt 20d)
    coeff=robust_poly_fit(waves[p],a_factor[p], 4,/double)
    resistant_mean,a_factor[p]-poly(waves[p],coeff),5.0,mean,good=k0
    sset=bspline_iterfit(waves[p[k0]],a_factor[p[k0]],maxiter=0,requiren=10,bkspace=5,nord=8) 
    afit0=bspline_valu(waves,sset)
    resistant_mean,a_factor[p]-afit0[p],5.0,mean,good=q0
    ww=dindgen((maxwave-minwave+10d)*5d)/5d + minwave -5d
    if n_elements(order) eq 0 then begin
        ; perform a bspline fit by default
        sset=bspline_iterfit(waves[p[q0]],a_factor[p[q0]],maxiter=0,requiren=10,bkspace=5,nord=8)
        afit0=bspline_valu(ww,sset)
    endif else begin
        ; perform a quadratic polynomial fit (for times from 0->453)
        order=order>1
        coeff=robust_poly_fit(waves[p[q0]],a_factor[p[q0]], order,/double)
        afit0 = poly(ww, coeff)
        print,'order ',order,' coeff=[',coeff,']'
    endelse

    if n_elements(waves) ge 2 and not keyword_set(noplot) then begin
        plot_multi,waves,a_factor, ww, afit0, /xst,/yst,xtitle='Wavelength (nm)', $
            ytitle=strtrim(string(modes[0]),2)+' and '+strtrim(string(modes[1]),2)+' RayPath A-Factor',$
            title='SimA vs SimB RayPath from '+strtrim(string(t0,format='(f0.1)'),2)+' to '+$
            strtrim(string(t1,format='(f0.1)'),2),charsize=1.4, psym=[-4,-3],thick=[0,2]
    endif


    return,{wavelength:waves, a_factor:a_factor, good:q0, wavefit:ww, a_fit:afit0}

end
