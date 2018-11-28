;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Estimates the degradation by comparing the data from
;   the specified diode detector from SimA and SimB.
;
; CALLING SEQUENCE:
;
; INPUT PARAMETERS:
;   startTime -
;      The lower time range for which data will be returned.
;      Expects times to be in mission days.
;   stopTime -
;      The upper time range for which data will be returned.
;      Expects times to be in mission days.
;   instrumentModeId -
;      Instrument modes to process. Expecting a 2 element array
;      for Sima and SimB.
;   wavelength - 
;      Wavelength to process.
;
; OPTIONAL INPUT PARAMETERS:
;   version -
;      Requests a specific version of the SimProfileIntegralCal, otherwise
;      pick the latest version valid for the specified time range.
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
;   Revision: $Id: tnmin_pddeg.pro,v 1.2 2018/11/26 16:16:02 spenton Exp $
;------------------------------------------------------------------
;
function mytnmin_pdfunc, x, df, solexpA=solexpA, irdA=irdA, solexpB=solexpB, irdB=irdB, $
    new_irdA=new_irdA, new_irdB=new_irdB
    ; measured_irrad = SolarIrrad * (exp(-abs(x[0]*solexp)) + exp(-abs(x[1]*solexp))) / 2
    ; SolarIrrad = measured_irrad / ((exp(-abs(x[0]*solexp)) + exp(-abs(x[1]*solexp))) / 2)
    new_irdA = irdA / ((exp(-abs(x[0]*solexpA)) + exp(-abs(x[1]*solexpA))) / 2d)
    new_irdB = irdB / ((exp(-abs(x[0]*solexpB)) + exp(-abs(x[1]*solexpB))) / 2d)
    F = new_irdB - new_irdA
    out_value = abs(mean(F))
    ; evaluate the derivative 
    ;tempA = irdA * exp(abs(x[0]*solexpA + x[1])) * (x[0]*solexpA + x[1]) / abs(x[0]*solexpA + x[1])
    ;tempB = irdB * exp(abs(x[0]*solexpB + x[1])) * (x[0]*solexpB + x[1]) / abs(x[0]*solexpB + x[1])
    ;dfd0 = tempA * solexpA - tempB * solexpB
    ;dfd1 = tempA - tempB
    ;df=moment(df0, mdev=out0)
    ;df=moment(df1, mdev=out1)
    ;df = [out0,out1]
    print,x[0],x[1],out_value
    return, out_value
end

;*******************************************************************
function mytnmin_pdfunc2, x, df, solexpA=solexpA, irdA=irdA, solexpB=solexpB, irdB=irdB, $
    new_irdA=new_irdA, new_irdB=new_irdB
    ; measured_irrad = SolarIrrad * (exp(-abs(x[0]*solexp)) + exp(-abs(x[1]*solexp))) / 2
    ; SolarIrrad = measured_irrad / ((exp(-abs(x[0]*solexp)) + exp(-abs(x[1]*solexp))) / 2)
    new_irdA = irdA / ((exp(-abs(x[0]*solexpA)) + exp(-abs(x[1]*solexpA))) / 2d)
    new_irdB = irdB / ((exp(-abs(x[0]*solexpB)) + exp(-abs(x[1]*solexpB))) / 2d)
    F = new_irdA / new_irdB
    out_value = abs(mean(1.0d - F))
    ; evaluate the derivative 
    ;tempA = irdA * exp(abs(x[0]*solexpA + x[1])) * (x[0]*solexpA + x[1]) / abs(x[0]*solexpA + x[1])
    ;tempB = irdB * exp(abs(x[0]*solexpB + x[1])) * (x[0]*solexpB + x[1]) / abs(x[0]*solexpB + x[1])
    ;dfd0 = tempA * solexpA - tempB * solexpB
    ;dfd1 = tempA - tempB
    ;df=moment(df0, mdev=out0)
    ;df=moment(df1, mdev=out1)
    ;df = [out0,out1]
    print,x[0],x[1],out_value
    return, out_value
end

;*******************************************************************
function mytnmin_pdfunc3, x, df, solexpA=solexpA, irdA=irdA, solexpB=solexpB, irdB=irdB, $
    new_irdA=new_irdA, new_irdB=new_irdB, timeA=timeA, timeB=timeB
    ; measured_irrad = SolarIrrad * exp(-abs(x[0]*solexp))
    ; SolarIrrad = measured_irrad / exp(-abs(x[0]*solexp))
    new_irdA = irdA / exp(-abs(x[0]*solexpA))
    new_irdB = irdB / exp(-abs(x[0]*solexpB))

    ; since the same wavelengths data is at different times for A and B
    ; do a 2nd order fit to each iMode and then use the diff at same time
    coeffsA = robust_poly_fit(timeA, new_irdA, 2d, /double)
    coeffsB = robust_poly_fit(timeB, new_irdB, 2d, /double)
    if n_elements(new_irdA) ge n_elements(new_irdB) then begin
        fitA = poly(timeA,coeffsA)
        fitB = poly(timeA,coeffsB)
    endif else begin
        fitA = poly(timeB,coeffsA)
        fitB = poly(timeB,coeffsB)
    endelse
    F = fitB - fitA
    out_value = abs(mean(F))
    out_value = stddev(F)
    out_value = MEANABSDEV(F, /median)
    ;out_value=robust_sigma(F, /zero)
    ;print,x,out_value
    return, out_value
end

;*******************************************************************
function mytnmin_pdfunc4, x, df, solexpA=solexpA, irdA=irdA, solexpB=solexpB, irdB=irdB, $
    new_irdA=new_irdA, new_irdB=new_irdB, timeA=timeA, timeB=timeB, afactA=afactA, afactB=afactB
    ; measured_irrad = SolarIrrad * exp(-abs(x[0]*solexp))
    ; SolarIrrad = measured_irrad / exp(-abs(x[0]*solexp))
    new_irdA = irdA / ((1d - afactA)*exp(-abs(x[0]*solexpA)) + afactA*exp(-abs(x[0]*solexpA)/2d))
    new_irdB = irdB / ((1d - afactB)*exp(-abs(x[0]*solexpB)) + afactB*exp(-abs(x[0]*solexpB)/2d))

    ; since the same wavelengths data is at different times for A and B
    ; do a 2nd order fit to each iMode and then use the diff at same time
    coeffsA = robust_poly_fit(timeA, new_irdA, 2d, /double)
    coeffsB = robust_poly_fit(timeB, new_irdB, 2d, /double)
    if n_elements(new_irdA) ge n_elements(new_irdB) then begin
        fitA = poly(timeA,coeffsA)
        fitB = poly(timeA,coeffsB)
        F = fitB - fitA
        cc=ladfit(timeA,F)
    endif else begin
        fitA = poly(timeB,coeffsA)
        fitB = poly(timeB,coeffsB)
        F = fitB - fitA
        cc=ladfit(timeA,F)
    endelse
    ;out_value = abs(mean(F))
    out_value = stddev(F)
    ;out_value = MEANABSDEV(F, /median)
    ;out_value=robust_sigma(F, /zero)
    ;print,x,out_value
    ;return, out_value
    return,abs(cc[1])*1d6

end

;*******************************************************************
function mytnmin_pdfunc5, x, df, solexpA=solexpA, irdA=irdA, solexpB=solexpB, irdB=irdB, $
    new_irdA=new_irdA, new_irdB=new_irdB, timeA=timeA, timeB=timeB, afactA=afactA, afactB=afactB
    ; measured_irrad = SolarIrrad * exp(-abs(x[0]*solexp))
    ; SolarIrrad = measured_irrad / exp(-abs(x[0]*solexp))
    new_irdA = irdA / ((1d - afactA)*exp(-abs(x[0]*solexpA)) + afactA*exp(-abs(x[0]*solexpA)/2d))
    new_irdB = irdB / ((1d - afactB)*exp(-abs(x[0]*solexpB)) + afactB*exp(-abs(x[0]*solexpB)/2d))

    ; try to remove some of the weighting on days where we have lots of data points
    ; by doing a linear interpolation for every 0.5 day between the start and end
    mn=min([timeA,timeB],max=mx)
    tt= dindgen((ceil(mx)-floor(mn))*2+1)/2d + floor(mn)
    da=interpol(new_irdA, timeA, tt)
    db=interpol(new_irdB, timeB, tt)

    ; since the same wavelengths data is at different times for A and B
    ; do a 2nd order fit to each iMode and then use the diff at same time
    coeffsA = robust_poly_fit(tt, da, 3, /double)
    coeffsB = robust_poly_fit(tt, db, 3, /double)
    fitA = poly(tt,coeffsA)
    fitB = poly(tt,coeffsB)
    F = fitB - fitA
    cc=ladfit(tt,F)

    return,abs(cc[1])*1d6

end

;*******************************************************************

function tnmin_pddeg, t0, t1, instrumentModeId, wavelength, version=version, bversion=bversion, coeffs=coeffs, $
    dburl=dburl, dbdriver=dbdriver, user=user, password=password, server=server, nfev=nfev, niter=niter, $
    status=status, errmsg=errmsg, inspectrum=inspectrum, solar54=solar54, solar55=solar55, $
    fit_goodness=fit_goodness, submean=submean, ffunct=ffunct, afact=afact, irrad_table=irrad_table

    ; get all the data from specified instrumentModeId
    if n_elements(instrumentModeId) ne 2 then begin
        print,'Error:  instrumentModeId requires a value for SimA and a value for SimB'
        return,-1
    endif

    if n_elements(user) eq 0 then     user='sbeland'
    if n_elements(password) eq 0 then password='sbeland1'
    if n_elements(dburl) eq 0 then    dbUrl='jdbc:sybase:Tds:lasp-db-dev:4100/SORCE_SIM_DEV'
    if n_elements(server) eq 0 then   server='lasp-db-dev'
    if n_elements(dbdriver) eq 0 then dbDriver='com.sybase.jdbc3.jdbc.SybDriver'

    ; version 15 of development database contains SimCalibratedIrradiance without degradation
    if n_elements(version) eq 0 then version=2541
    if n_elements(bversion) eq 0 then bversion=version

    ; for now only process one wavelength at a time
    if n_elements(wavelength) gt 1 then wavelength=wavelength[0]

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

    if n_elements(irrad_table) eq 0 then irrad_table='SimUncorrectedIrradiance'

    if n_elements(ffunct) ge 2 and max(strpos(tag_names(solar54),'ONEAU')) eq -1 then begin
        print,'  getting SUNOBSERVERDISTANCECORRECTION from database ...'
        q1="SELECT microsecondsSinceGpsEpoch, sunObserverDistanceCorrection FROM SolarDistAndDopplerFixedStep2 "+$
           "WHERE instrumentModeId=50 and version=7"
        query_database, q1, solardist, info
        ;time54 = (solar54.(0)+solar54.(1))/2d
        ;time55 = (solar55.(0)+solar55.(1))/2d
        time54 = solar54.(1)
        time55 = solar55.(1)
        oneau54 = interpol(solardist.SUNOBSERVERDISTANCECORRECTION, solardist.MICROSECONDSSINCEGPSEPOCH, time54)
        append_tag,solar54,'oneau',oneau54,/slim
        oneau55 = interpol(solardist.SUNOBSERVERDISTANCECORRECTION, solardist.MICROSECONDSSINCEGPSEPOCH, time55)
        append_tag,solar55,'oneau',oneau55,/slim
    endif

    if n_elements(ffunct) ge 2 then begin
        ; apply the F-Function and the scaled 1AU to the degradation Column as in insert_degcol.pro
        ; the ffunct is expected to contain the polynomial coefficients 
        ;if n_elements(time54) eq 0 then time54 = (solar54.(0)+solar54.(1))/2d
        ;if n_elements(time55) eq 0 then time55 = (solar55.(0)+solar55.(1))/2d
        if n_elements(time54) eq 0 then time54 = solar54.(1)
        if n_elements(time55) eq 0 then time55 = solar55.(1)
        corr54 = poly(gps2sd(time54/1d6), ffunct)
        corr55 = poly(gps2sd(time55/1d6), ffunct)
        ; the cumulative solar exposure is in days to go with our Kappa function
        solar54.solar_exp=total(solar54.solar_exp_orbit*corr54,/cum) / (1d +(solar54.oneau-1d)/4d) /86400d
        solar55.solar_exp=total(solar55.solar_exp_orbit*corr55,/cum) / (1d +(solar55.oneau-1d)/4d) /86400d
    endif else begin
        solar54.solar_exp=total(solar54.solar_exp_orbit,/cum) /86400d
        solar55.solar_exp=total(solar55.solar_exp_orbit,/cum) /86400d
    endelse

    ; As a different test, scale the solar_exp for SimB with the profile of the ratio B/A (Version 70)
    ;sum55 = interpol(solar55.solar_exp, solar55.t0, solar54.t0)
    ;p=where(solar54.solar_exp gt 0d)
    ;sum55[p] *= ((1d - sum55[p] / solar54[p].solar_exp) > 0d)
    ;solar55.solar_exp = interpol(sum55, solar54.t0, solar55.t0)

    ; get the data from SimB first and then find SimA data from same orbit
    if size(inspectrum,/tname) ne 'STRUCT' then begin
        ; look for all the SolarQuickScan24 within this time range
        query_database,/reset
        plansb = get_sorce_plan(t0, t1, /mission,/simb,activity='SolarQuickScan24')

        ; get rid of data at 2562.63
        p=where(abs(plansb.starttime - 2562.645) lt 0.02, complement=cp,count)
        if count gt 0 then plansb=plansb[cp]

        print,'   '+strtrim(string(n_elements(plansb)),2)+' spectrum to process ...'
        nspect = n_elements(plansb)
        inspectrum = {sima:{starttime:dblarr(nspect), diode:ptrarr(nspect)}, simb:{starttime:dblarr(nspect), diode:ptrarr(nspect)}}
        query_database,/reset
    endif else nspect=n_elements(inspectrum.simb.(0))

    ; extract the spectrum
    wv0 = wavelength - 10d
    wv1 = wavelength + 10d
    irradA=[]
    solarexpA=[]
    timeA=[]
    irradB=[]
    solarexpB=[]
    timeB=[]
    for pl=0, nspect-1L do begin
        ; extract SIMB data fisrt
        if NOT PTR_VALID(inspectrum.simb.diode[pl]) then begin
            gps0 = ulong64(sd2gps(plansb[pl].starttime)*1d6)
            gps1 = ulong64(sd2gps(plansb[pl].stoptime)*1d6)
            inst = max(instrumentModeId)
            q1='SELECT w.microsecondsSinceGpsEpoch as timestamp,w.instrumentModeId,w.version,w.wavelengthRef as wavelength,ird.irradiance '+$
               'FROM Wavelength w,'+irrad_table[0]+' ird where '+$
               'w.version='+strtrim(string(bversion),2)+' and w.version=ird.version and '+$
               'w.instrumentModeId='+strtrim(string(inst),2)+' and '+$ 
               'w.instrumentModeId=ird.instrumentModeId and '+$
               'w.microsecondsSinceGpsEpoch>='+strtrim(string(gps0),2)+' and '+$
               'w.microsecondsSinceGpsEpoch<='+strtrim(string(gps1),2)+' and '+$
               'w.microsecondsSinceGpsEpoch=ird.microsecondsSinceGpsEpoch'
            query_database,q1,res,user=user, password=password, dburl=dburl
            if size(res,/tname) ne 'STRUCT' then begin
                print,'Error: no spectra found for SimB for this plan ',plansb[pl]
                continue
            endif
            inspectrum.simb.diode[pl]= ptr_new(res)
            inspectrum.simb.starttime[pl]= sd2gps(plansb[pl].starttime)*1d6
        endif

        ; extract SIMA data 
        if NOT PTR_VALID(inspectrum.sima.diode[pl]) then begin
            gps0 = ulong64(sd2gps(plansb[pl].starttime)*1d6)
            gps1 = ulong64(sd2gps(plansb[pl].stoptime)*1d6)
            inst = min(instrumentModeId)
            q1='SELECT w.microsecondsSinceGpsEpoch as timestamp,w.instrumentModeId,w.version,w.wavelengthRef as wavelength,ird.irradiance '+$
               'FROM Wavelength w,'+irrad_table[0]+' ird where '+$
               'w.version='+strtrim(string(version),2)+' and w.version=ird.version and '+$
               'w.instrumentModeId='+strtrim(string(inst),2)+' and '+$ 
               'w.instrumentModeId=ird.instrumentModeId and '+$
               'w.microsecondsSinceGpsEpoch>='+strtrim(string(gps0),2)+' and '+$
               'w.microsecondsSinceGpsEpoch<='+strtrim(string(gps1),2)+' and '+$
               'w.microsecondsSinceGpsEpoch=ird.microsecondsSinceGpsEpoch'
            query_database,q1,res,user=user, password=password, dburl=dburl
            if size(res,/tname) ne 'STRUCT' then begin
                print,'Error: no spectra found for SimA for this plan ',plansb[pl]
                continue
            endif
            inspectrum.sima.diode[pl]= ptr_new(res)
            inspectrum.sima.starttime[pl]= sd2gps(plansb[pl].starttime)*1d6
        endif

        ; if here, we have valid spectrum
        ;s=sort((*inspectrum.simb.diode[pl]).wavelength)
        ;y2=spl_init((*inspectrum.simb.diode[pl])[s].wavelength,(*inspectrum.simb.diode[pl])[s].irradiance)
        ;irrad=spl_interp((*inspectrum.simb.diode[pl])[s].wavelength,(*inspectrum.simb.diode[pl])[s].irradiance, y2, wavelength)
        irrad = interpol((*inspectrum.simb.diode[pl]).irradiance, (*inspectrum.simb.diode[pl]).wavelength, wavelength, /spline) ;,/lsq)
        ;spline_p,(*inspectrum.simb.diode[pl]).wavelength,(*inspectrum.simb.diode[pl]).irradiance,xr,yr,TAN0=[1,0], TAN1=[1,0],interval=0.05
        ;irrad = interpol(yr, xr, wavelength)
        irradB=[irradB,irrad]
        ; use the cumulative solar_exposure at the end of the previous orbit
        pmin=min(abs((*inspectrum.simb.diode[pl]).wavelength-wavelength),pos)
        p=where(solar55.t1 le (*inspectrum.simb.diode[pl])[pos].timestamp,count)
        solarexpB=[solarexpB, solar55[p[-1]].solar_exp]
        timeB=[timeB,(*inspectrum.simb.diode[pl])[pos].timestamp]

        ;s=sort((*inspectrum.simA.diode[pl]).wavelength)
        ;y2=spl_init((*inspectrum.simA.diode[pl])[s].wavelength,(*inspectrum.simA.diode[pl])[s].irradiance)
        ;irrad=spl_interp((*inspectrum.simA.diode[pl])[s].wavelength,(*inspectrum.simA.diode[pl])[s].irradiance, y2, wavelength)
        irrad = interpol((*inspectrum.sima.diode[pl]).irradiance, (*inspectrum.sima.diode[pl]).wavelength, wavelength, /spline) ;,/lsq)
        ;spline_p,(*inspectrum.simA.diode[pl]).wavelength,(*inspectrum.simA.diode[pl]).irradiance,xr,yr,TAN0=[1,0], TAN1=[1,0],interval=0.05
        ;irrad = interpol(yr, xr, wavelength)
        irradA=[irradA,irrad]
        ; use the cumulative solar_exposure at the end of the previous orbit
        pmin=min(abs((*inspectrum.sima.diode[pl]).wavelength-wavelength),pos)
        p=where(solar54.t1 le (*inspectrum.sima.diode[pl])[pos].timestamp,count)
        solarexpA=[solarexpA, solar54[p[-1]].solar_exp]
        timeA=[timeA,(*inspectrum.sima.diode[pl])[pos].timestamp]

    endfor

    ; get rid of missing spectrum from inspectrum array
    pos = where(ptr_valid(inspectrum.sima.diode) eq 1, count)
    if count eq 0 then begin
        print,'Error: no valid spectrum found for SimA and SimB within time range'
        return,-1
    endif else if count gt 0 and count ne n_elements(inspectrum.sima.diode) then begin
        inspectrum = {sima:{starttime:inspectrum.sima.starttime[pos], diode:inspectrum.sima.diode[pos]}, $
                      simb:{starttime:inspectrum.simb.starttime[pos], diode:inspectrum.simb.diode[pos]}}
    endif

    ; now we want to minimize the difference in irradiance between SimA and SimB
    ; when fitting the same exponential degradation model
    ; convert the cumulative exposure in days
    functargs = {solexpA:solarexpA, irdA:irradA, solexpB:solarexpB, irdB:irradB, timeA:timeA, timeB:timeB}
    if n_elements(afact) eq 0 then begin
        afact_valA=0d
        afact_valB=0d
    endif else begin
        ; afactor value is a single element at the requested wavelength
        afact_valA=(interpol(afact.a, afact.wavelength, wavelength, /spline))[0]
        afact_valB = afact_valA
        p=where(strpos(tag_names(solar54),'ONEAU') ge 0,oneau_flag)
        if (oneau_flag) gt 0 then begin
            afact_valA *= (interpol(solar54.oneau,solar54.t1,timeA))^2 
            afact_valB *= (interpol(solar54.oneau,solar54.t1,timeB))^2 
        endif
    endelse
    timeA=gps2sd(timeA/1d6)
    timeB=gps2sd(timeB/1d6)

    ; now we want to minimize the difference in irradiance between SimA and SimB
    ; when fitting the same exponential degradation model
    ; convert the cumulative exposure in days
    functargs = {solexpA:solarexpA, irdA:irradA, solexpB:solarexpB, irdB:irradB, timeA:timeA, timeB:timeB, $
                 afactA:afact_valA, afactB:afact_valB}

    parinfo=replicate({value:0.0d, fixed:0, limited:[0,0], limits:[0d,0d], step:0d, tnside:2},1)
    parinfo[0].value=1.0d-3
    ;parinfo[1].limited=[1,1]
    ;parinfo[1].limits=[-1d-11,1d-11]
    ;parinfo[1].value=1.0d-10

    ;res=ladfit([solarexpB,solarexpA],[irradB,irradA])
    ;print,res
    ;print,'  initial guess: ',parinfo[*].value
    coeffs = tnmin('mytnmin_pdfunc5', coeff0, functargs=functargs, bestmin=f0, status=status, $
             nfev=nfev, niter=niter, autoderivative=1, errmsg=errmsg, parinfo=parinfo,/quiet)
    coeffs=abs(coeffs)
    fit_goodness=mytnmin_pdfunc5(coeffs, df, solexpA=solarexpA, irdA=irradA, solexpB=solarexpB, irdB=irradB, $
        new_irdA=new_irdA, new_irdB=new_irdB, timeA=timeA, timeB=timeB, afactA=afact_valA, afactB=afact_valB)

    if keyword_set(submean) then begin
        ; subtract the difference of the mean between SimA and SimB to SIMA and re-process
        meandiff = mean(new_irdB) - mean(new_irdA)
        mod_irradA = irradA + meandiff
        fit_goodness=mytnmin_pdfunc5(coeffs, df, solexpA=solarexpA, irdA=mod_irradA, solexpB=solarexpB, $
            irdB=irradB, new_irdA=new_irdA, new_irdB=new_irdB, timeA=timeA, timeB=timeB, afactA=afact_valA, afactB=afact_valB)
        print,'  difference of means: ',meandiff
    endif

    print,coeffs, 1d - fit_goodness, format='("  final value: ",D0.8, "    (",D0.8,")")'

    nA=n_elements(new_irdA)
    nB=n_elements(new_irdB)
    spA={timestamp:sd2gps(timeA)*1d6, wavelength:replicate(wavelength,nA), solarexp:solarexpA, $
        irradiance:irradA, new_irrad:new_irdA}
    spB={timestamp:sd2gps(timeB)*1d6, wavelength:replicate(wavelength,nB), solarexp:solarexpB, $
        irradiance:irradB, new_irrad:new_irdB}
    return, {simA:spA, simB:spB}

end
