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
;   Revision: $Id: tnmin_pdtau.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;------------------------------------------------------------------
function mytnmin_pdtau, x, df, solexpA=solexpA, irdA=irdA, solexpB=solexpB, irdB=irdB, $
    new_irdA=new_irdA, new_irdB=new_irdB, timeA=timeA, timeB=timeB, afact=afact
    ; measured_irrad = SolarIrrad * exp(-abs(x[0]*solexp))
    ; SolarIrrad = measured_irrad / exp(-abs(x[0]*solexp))
    new_irdA = irdA / ((1d - afact[0])*exp(-abs(x[0]*solexpA)) + afact[0]*exp(-abs(x[0]*solexpA)/2d))
    new_irdB = irdB / ((1d - afact[0])*exp(-abs(x[0]*solexpB)) + afact[0]*exp(-abs(x[0]*solexpB)/2d))

    ; since the same wavelengths data is at different times for A and B
    ; do a 2nd order fit to each iMode and then use the diff at same time
    coeffsA = robust_poly_fit(timeA, new_irdA, 2d, /double)
    coeffsB = robust_poly_fit(timeB, new_irdB, 2d, /double)
    if coeffsA[0] eq 0d or coeffsB[0] eq 0d then begin
        coeffsA = poly_fit(timeA, new_irdA, 2d, /double)
        coeffsB = poly_fit(timeB, new_irdB, 2d, /double)
    endif
    if n_elements(new_irdA) ge n_elements(new_irdB) then begin
        fitA = poly(timeA,coeffsA)
        fitB = poly(timeA,coeffsB)
    endif else begin
        fitA = poly(timeB,coeffsA)
        fitB = poly(timeB,coeffsB)
    endelse
    F = fitB - fitA
    ;out_value = abs(mean(F))
    out_value = stddev(F)
    ;out_value = MEANABSDEV(F, /median)
    ;out_value=robust_sigma(F, /zero)
    ;print,x,out_value
    return, out_value
end

;*******************************************************************

function tnmin_pdtau, wavelength, inspectrum, coeffs=coeffs, $
    status=status, errmsg=errmsg, solar54=solar54, solar55=solar55, $
    fit_goodness=fit_goodness, submean=submean, ffunct=ffunct, afact=afact, kappa=kappa, noplot=noplot

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

    if max(strpos(tag_names(solar54),'ONEAU')) eq -1 then begin
        print,'  getting SUNOBSERVERDISTANCECORRECTION from database ...'
        q1="SELECT microsecondsSinceGpsEpoch, sunObserverDistanceCorrection FROM SolarDistAndDopplerFixedStep2 "+$
           "WHERE instrumentModeId=50 and version=7"
        query_database, q1, solardist, info
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
        solar54.solar_exp=total(solar54.solar_exp_orbit,/cum) / (1d +(solar54.oneau-1d)/4d) /86400d
        solar55.solar_exp=total(solar55.solar_exp_orbit,/cum) / (1d +(solar54.oneau-1d)/4d) /86400d
    endelse

    ; extract the spectrum
    irradA=[]
    solarexpA=[]
    timeA=[]
    irradB=[]
    solarexpB=[]
    timeB=[]
    nspect = n_elements(inspectrum.simb)

    for pl=0, nspect-1L do begin
        ; if here, we have valid spectrum
        s=sort(inspectrum.simb[pl].wavelength)
        irrad = interpol(inspectrum.simb[pl].irradiance[s], inspectrum.simb[pl].wavelength[s], wavelength, /spline) ;,/lsq)
        irradB=[irradB,irrad]
        ; use the cumulative solar_exposure at the end of the previous orbit
        pmin=min(abs(inspectrum.simb[pl].wavelength-wavelength),pos)
        p=where(solar55.t1 le inspectrum.simb[pl].TIMESTAMP[pos],count)
        solarexpB=[solarexpB, solar55[p[-1]].solar_exp]
        timeB=[timeB,inspectrum.simb[pl].TIMESTAMP[pos]]

        ; find the closest sima spectra in time with simB
        mn=min(abs(inspectrum.simb[pl].TIMESTAMP[0] - inspectrum.sima.TIMESTAMP[0]),pos0)
        s=sort(inspectrum.sima[pos0].wavelength)
        irrad = interpol(inspectrum.sima[pos0].irradiance[s], inspectrum.sima[pos0].wavelength[s], wavelength, /spline) ;,/lsq)
        irradA=[irradA,irrad]
        ; use the cumulative solar_exposure at the end of the previous orbit
        pmin=min(abs(inspectrum.sima[pos0].wavelength-wavelength),pos1)
        p=where(solar54.t1 le inspectrum.sima[pos0].TIMESTAMP[pos1],count)
        solarexpA=[solarexpA, solar54[p[-1]].solar_exp]
        timeA=[timeA,inspectrum.sima[pos0].TIMESTAMP[pos1]]

    endfor

    resistant_mean, irradA, 5.0, mean, good=ka
    if n_elements(ka) ne n_elements(irradA) then begin
        irradA = irradA[ka]
        timeA = timeA[ka]
        solarExpA = solarExpA[ka]
        irradB = irradB[ka]
        timeB = timeB[ka]
        solarExpB = solarExpB[ka]
    endif
    resistant_mean, irradB, 5.0, mean, good=kb
    if n_elements(ka) ne n_elements(irradB) then begin
        irradA = irradA[kb]
        timeA = timeA[kb]
        solarExpA = solarExpA[kb]
        irradB = irradB[kb]
        timeB = timeB[kb]
        solarExpB = solarExpB[kb]
    endif

    timeA=gps2sd(timeA/1d6)
    timeB=gps2sd(timeB/1d6)

    if n_elements(irradA) lt 4 then return,-1

    ; now we want to minimize the difference in irradiance between SimA and SimB
    ; when fitting the same exponential degradation model
    functargs = {solexpA:solarexpA, irdA:irradA, solexpB:solarexpB, irdB:irradB, timeA:timeA, timeB:timeB}
    if n_elements(afact) eq 0 then begin
        afact_val=0d
    endif else begin
        ; afactor value is a single element at the requested wavelength
        afact_val=interpol(afact.a, afact.wavelength, wavelength)
    endelse
    functargs = create_struct(functargs, 'afact', afact_val)

    if n_elements(kappa) gt 0 then kval=interpol(kappa.y, kappa.x, wavelength) else kval=1d-3

    parinfo=replicate({value:0.0d, fixed:0, limited:[0,0], limits:[0d,0d], step:0d, tnside:2},1)
    parinfo[0].value=kval
    ;parinfo[1].limited=[1,1]
    ;parinfo[1].limits=[-1d-11,1d-11]
    ;parinfo[1].value=1.0d-10

    ;print,'  initial guess: ',parinfo[*].value
    coeffs = tnmin('mytnmin_pdtau', coeff0, functargs=functargs, bestmin=f0, status=status, $
             nfev=nfev, niter=niter, autoderivative=1, errmsg=errmsg, parinfo=parinfo, /quiet)
    coeffs=abs(coeffs)
    if finite(coeffs) eq 0 then return,-1   ; in case it doesn't converge

    fit_goodness=mytnmin_pdtau(coeffs, df, solexpA=solarexpA, irdA=irradA, solexpB=solarexpB, irdB=irradB, $
        new_irdA=new_irdA, new_irdB=new_irdB, timeA=timeA, timeB=timeB, afact=afact_val)

    if keyword_set(submean) then begin
        ; subtract the difference of the mean between SimA and SimB to SIMA and re-process
        meandiff = mean(new_irdB) - mean(new_irdA)
        mod_irradA = irradA + meandiff
        fit_goodness=mytnmin_pdtau(coeffs, df, solexpA=solarexpA, irdA=mod_irradA, solexpB=solarexpB, $
            irdB=irradB, new_irdA=new_irdA, new_irdB=new_irdB, timeA=timeA, timeB=timeB, afact=afact_val)
        ;print,'  difference of means: ',meandiff
    endif

    ;print,'  initial guess:',parinfo[*].value,'  final value:',coeffs,'     fit goodness:', 1d -fit_goodness

    if NOT keyword_set(noplot) then begin
        plot_multi,timeA, irradA, timeA, new_irda, timeb, irradB, timeB, new_irdB, $
            /xst,/yst,label=['IrradA','New_IrradA','IrradB','New_IrradB']
    endif

    nA=n_elements(new_irdA)
    nB=n_elements(new_irdB)
    spA={MICROSECONDSSINCEGPSEPOCH:sd2gps(timeA)*1d6, wavelength:replicate(wavelength,nA), solarexp:solarexpA, $
        irradiance:irradA, new_irrad:new_irdA}
    spB={MICROSECONDSSINCEGPSEPOCH:sd2gps(timeB)*1d6, wavelength:replicate(wavelength,nB), solarexp:solarexpB, $
        irradiance:irradB, new_irrad:new_irdB}
    return, {simA:spA, simB:spB}

end
