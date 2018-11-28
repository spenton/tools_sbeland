;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Find the best fit 1AU factor at the specified wavelength for the time
;   range requested using the MPFITFUN routine to minimize the difference 
;   between the F_FACTOR_SCALED and the polynomial fit.
;   This is all doen for the ESR data for now using the get_esrdeg routine.
;
; CALLING SEQUENCE:
;   spectra = MPFIT_ONEAU(t0, t1, /missionDay, wavelength)
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
;
; OPTIONAL INPUT PARAMETERS:
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
;   version -
;      Requests a specific version of the SimProfileIntegralCal, otherwise
;      pick the latest version valid for the specified time range.
;   order - 
;      Order of the polynomial correction.
;
; RETURNED PARAMETERS:
;   A structure with the spectra aligned to the reference spectra in wavelength.
;
; OPTIONAL OUTPUT PARAMETERS:
;
; EXAMPLE:  
;
; COMMON BLOCKS:
;      NONE
;
;
; REVISION HISTORY:
;   Revision: $Id: mpfit_oneau.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;------------------------------------------------------------------
;
;  Function to apply a polynomial to the wavelength scale of the spectra to compare to
function myfitoneau, x, df, t0=t0, t1=t1, solar54=solar54, solar55=solar55, $
        inspectrum=inspectrum, order=order, step55=step55, result=result, wavelength=wavelength

    solar54.solar_exp=total(solar54.solar_exp_orbit / (1d +(solar54.oneau-1d)/x[0]),/cum)
    solar55.solar_exp=total(solar55.solar_exp_orbit / (1d +(solar55.oneau-1d)/x[0]),/cum)
    result=get_esrdeg(t0, t1, wavelength, /mission, solar54=solar54, solar55=solar55, $
        inspectrum=inspectrum, /alignobc, /smooth,step55=step55, /noplot)

    ; calculate the residuals from the polynomial fit up to SD=2805
    p=where(result.time_f lt 2805)
    coeff=robust_poly_fit(result.time_f[p], result.f_factor_scaled[p], order,/double)
    yfit=poly(result.time_f,coeff)
    value = variance((yfit[p]-result.f_factor_scaled[p]))
    print,x,value
    return,value
end

;*******************************************************************

function mpfit_oneau, startTime, stopTime, wavelength, $
         gps=gps, missionDays=missionDays, julianDays=julianDays, $
         version=version, order=order, inspectrum=inspectrum, $
         solar54=solar54, solar55=solar55, step55=step55


    if keyword_set(missionDays) then begin
       ; user specified time in mission (sorce) days
       t0 = startTime
       t1 = stopTime
    endif else if keyword_set(julianDays) then begin
       ; user specified time in julian days
       t0 = jd2sd(startTime)
       t1 = jd2sd(stopTime)
    endif else begin
       ; user specified timetags in gps microseconds
       t0 = gps2sd(startTime/1.d6)
       t1 = gps2sd(stopTime/1.d6)
    endelse

    if n_elements(step55) eq 0 then step55=90
    if n_elements(order) eq 0 then order=2
    result=[]

    if size(solar54,/tname) ne 'STRUCT'  then begin
        print,'  getting solar54 from database ...'
        q1 = 'SELECT orbitStartTimeGps t0, orbitEndTimeGps t1, solarExposureHrtOut as solar_exp_orbit, '+$
             'cumulativeSolarExpHrtOut as solar_exp from SimSolarExposureData where instrumentModeId=54'
        query_database, q1, solar54, info
    endif

    if size(solar55,/tname) ne 'STRUCT'  then begin
        print,'  getting solar55 from database ...'
        q1 = 'SELECT orbitStartTimeGps t0, orbitEndTimeGps t1, solarExposureHrtOut as solar_exp_orbit, '+$
             'cumulativeSolarExpHrtOut as solar_exp from SimSolarExposureData where instrumentModeId=55'
        query_database, q1, solar55, info
    endif

    if max(strpos(tag_names(solar54),'ONEAU')) eq -1 then begin
        print,'  getting SUNOBSERVERDISTANCECORRECTION from database ...'
        q1="SELECT microsecondsSinceGpsEpoch, sunObserverDistanceCorrection FROM SolarDistAndDopplerFixedStep2 "+$
           "WHERE instrumentModeId=50 and version=7"
        query_database, q1, solardist, info
        oneau54 = interpol(solardist.SUNOBSERVERDISTANCECORRECTION, solardist.MICROSECONDSSINCEGPSEPOCH, solar54.t1)
        append_tag,solar54,'oneau',oneau54,/slim
        oneau55 = interpol(solardist.SUNOBSERVERDISTANCECORRECTION, solardist.MICROSECONDSSINCEGPSEPOCH, solar55.t1)
        append_tag,solar55,'oneau',oneau55,/slim
    endif

    solar54.solar_exp=total(solar54.solar_exp_orbit,/cum)
    solar55.solar_exp=total(solar55.solar_exp_orbit,/cum)

    result=get_esrdeg(t0, t1, /mission,  wavelength, solar54=solar54, solar55=solar55, $
        inspectrum=inspectrum, /alignobc, /smooth, step55=step55, /noplot)

    parinfo=replicate({value:5.0d, fixed:0, limited:[0,0], limits:[0d,0d], step:0.1d, tnside:2},1)
    ;parinfo[0].value=1.0d-4
    ;parinfo[0].fixed=1
    ;parinfo[0].step=1.0d-8
    parinfo[0].limited=[1,1]
    parinfo[0].limits=[0.2d,100]

    functargs = {t0:t0, t1:t1, solar54:solar54, solar55:solar55, inspectrum:inspectrum, order:order, $
        step55:step55, result:result, wavelength:wavelength}

    print,'  initial guess: ',parinfo.value
    coeffs = tnmin('myfitoneau', functargs=functargs, bestmin=f0, status=status, $
        nfev=nfev, niter=niter, autoderivative=1, errmsg=errmsg, parinfo=parinfo)

    return, {wavelength:wavelength, oneau_corr:coeffs}

end
