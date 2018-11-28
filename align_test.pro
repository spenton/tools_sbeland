;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Fit a spectra to a reference spectra by using the
;   MPFITFUN routine to minimize the difference with the reference when 
;   adjusting the prismposition scale with a polynomial.  The routine
;   finds the coefficients of the best fit and returns the original spectra 
;   with the new wavelength scale.
;
; CALLING SEQUENCE:
;   spectra = MPFIT_SPECTRA(t0, t1, instrumentModeId, /missionDays, refT0=t0, refT1=t1)
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
;      The instrument mode of interest are currently limited to:
;      41	SIM_A	VIS1
;      43	SIM_A	UV
;      44	SIM_A	IR
;      45	SIM_B	VIS1
;      47	SIM_B	UV
;      31	SIM_A	ESR
;      32	SIM_B	ESR
;
; OPTIONAL INPUT PARAMETERS:
;   spect - 
;      Structure containing the spectra to align. If provided, the
;      startTime, stopTime and instrumentModeId are ignored.
;   refSpect - 
;      Structure containg the reference spectra to align to.
;   refT0 -
;      StartTime of the reference spectra in the same units as T0.
;      If not specified, uses our "standard" reference time range.
;   refT1 -
;      StopTime of the reference spectra in the same units as T1.
;      If not specified, uses our "standard" reference time range.
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
;   quiet -
;      Do not print extra messages
;   noderiv -
;      By default, the program fits the derivative of the spectrum. This flag
;      will let the program use the original spectrum.
;
; RETURNED PARAMETERS:
;   A structure with the spectra aligned to the reference spectra in wavelength.
;
; OPTIONAL OUTPUT PARAMETERS:
;   coeffs -
;      The coefficients of the best polynomial fit.
;   status -
;      Returned STATUS of the mpfitfun routine.
;
; EXAMPLE:  
;   Adjust the wavelength scale of the SimB UV diode to the SimA UV diode for 
;   the SolarQuickScan24 on day 453.67:
;
;   IDL> T0=453.67840d
;   IDL> T1=453.69524d
;   IDL> dbtables=['SimProfileIntegral','SimCalibratedIrradiance']
;   IDL> res43=get_science_product(dbtables, T0, T1, 43, /mission, version=17)
;   IDL> res=align_spectra(T0,T1,47,/mission,refspect=res43,coeffs=coeffs,status=status)
;   IDL> help,status
;   IDL> STATUS          INT       =        2
;   IDL> help,res
;   ** Structure <3877418>, 2 tags, length=18000, data length=18000, refs=1:
;      WAVELENGTH      DOUBLE    Array[1125]
;      IRRADIANCE      DOUBLE    Array[1125]
;
; COMMON BLOCKS:
;      NONE
;
; NOTES:
;  You may want to run get_sorce_plan before hand to make sure an activity
;  with observations in the desired mode was performed during the time
;  span of interest.
;
;
;   STILL A WORK IN PROGRESS:  WORKS FOR SOME DATA BUT NOT WITH OTHER
;   -----------------------------------------------------------------
;
; REVISION HISTORY:
;   Revision: $Id: align_test.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;------------------------------------------------------------------
;
;  Function to apply a polynomial to the wavelength scale of the spectra to compare to
function myfitfunc, x, coeffs, spectx=spectx, specty=specty
    newx=poly(spectx,coeffs)
    newy=interpol(specty,newx,x,/lsq)
    return, newy
end
;*******************************************************************

function align_test, startTime, stopTime, instrumentModeId, $
         gps=gps, missionDays=missionDays, julianDays=julianDays, $
         refT0=refT0, refT1=refT1, refSpect=refSpect, spect=spect, $
         version=version, quiet=quiet, order=order, coeffs=coeffs, $
         noderiv=noderiv, status=status, profile_data=profile_data, $
         wrange=wrange

    if n_elements(version) eq 0 then version=19

    ; if a reference spectra was provided use that instead
    if size(refSpect,/tname) ne 'STRUCT' then begin
        if n_elements(refT0) eq 0 or n_elements(refT1) eq 0 then begin
            ; this is our default REFERENCE SolarQuickScan24 spectra 
            if instrumentModeId eq 31 or instrumentModeId eq 32 then begin
                ; full ESR scan
                refT0 = 453.02065d
                refT1 = 453.98991d
                ;columnNames=['SimProfileIntegral','SimCalibratedIrradiance','SimPhaseDetectedDn']
            endif else begin
                refT0 = 453.67840d
                refT1 = 453.69524d
                ;columnNames=['SimProfileIntegral','SimCalibratedIrradiance','SimConvertedDataNumbers']
            endelse
            if keyword_set(missionDays) then begin
               ; do nothing here since already in mission days
            endif else if keyword_set(julianDays) then begin
               ; user specified time in julian days
               refT0 = sd2jd(refT0)
               refT1 = sd2jd(refT1)
            endif else begin
               ; user specified time in gps days (default)
               refT0 = sd2gps(refT0)*1.d6
               refT1 = sd2gps(refT1)*1.d6
            endelse
        endif

        ; get the reference spectra
        ;refSpect=get_science_product(columnNames, refT0, refT1, instrumentModeId, mission=missionDays, $
        ;    gps=gps, julian=julianDays, version=version)
        refSpect = get_sim_spectra(refT0, refT1, instrumentModeId, gps=gps, missionDays=missionDays, $
            julianDays=julianDays, version=version, /noCcdCorr, /dn_only, profile_data=profile_data)
        if size(refSpect,/tname) ne 'STRUCT' then begin
            print,'No valid reference spectra found'
            return,-1
        endif

    endif

    ; get the spectra to fit to the reference
    if size(spect,/tname) ne 'STRUCT' then begin
        ;spect=get_science_product(columnNames, startTime, stopTime, instrumentModeId, mission=missionDays, $
        ;   gps=gps, julian=julianDays, version=version)
        spect = get_sim_spectra(startTime, stopTime, instrumentModeId, gps=gps, missionDays=missionDays, $
            julianDays=julianDays, version=version, /noCcdCorr, /noTempCorr, /dn_only, profile_data=profile_data)
        if size(spect,/tname) ne 'STRUCT' then begin
            print,'No valid spectra found'
            return,-1
        endif
    endif

    ; re-map the refspect wavelength by adding the dLambda/dTemperature corresponding to the spect temperature
    match, refSpect.ccdpos, spect.ccdpos, suba, subb
    refWave = refSpect[suba].wavelength

    ; get the wavelength correction at each spect temperature and apply to our refSpect
    dndt = interpol(profile_data.y13, profile_data.y4, refWave, /spline)
    ; since our call to get_sim_spectra for refSpect, corrects for the prism temp, the wavelengths are now
    ; for the reference temperature
    refWave = refWave - (spect[subb].prismtemp - profile_data[0].y0) * dndt
    ; now convert our refWave to ccdpos - corresponding to refSpect with same temperature profile as spect
    ccdPos = interpol(profile_data.x, profile_data.y4, refWave, /spline)

    if n_elements(wrange) eq 2 then begin
        wpos=where(refWave ge wrange[0] and refWave le wrange[1],count)
        if count eq 0 then begin
            print,'Error:  no wavelengths found within provided wrange'
            return,-1
        endif
        rsp={x:ccdPos[wpos], y:(refSpect[suba].dn)[wpos]}
        sp={x:(spect[subb].ccdpos)[wpos], y:(spect[subb].dn)[wpos]}
    endif else begin
        rsp={x:ccdPos, y:refSpect[suba].dn}
        sp={x:spect[subb].ccdpos, y: spect[subb].dn}
    endelse

    ; now try to align the spect with the temperature modified refSpect

    if NOT keyword_set(noderiv) then begin
        ; compare the derivatives instead of the actual spectrum
        rsp.y = deriv(rsp.x,deriv(rsp.x,rsp.y))
        sp.y = deriv(sp.x,deriv(sp.x,sp.y))
    endif

    if n_elements(order) eq 0 then order=1 else order=order>1
    ; get the global offset using cross_correlate (in number of steps)
    irrad=interpol(sp.y,sp.x,rsp.x,/lsq)
    cross_correlate,rsp.y,irrad,offset
    ; get the offset in wavelength
    nm_per_step = median(rsp.x[1:-1]-rsp.x[0:-2])
    parinfo=replicate({value:1.0d-4, fixed:0, limited:[0,0], limits:[0d,0d], step:0, mpside:2},order+1)
    parinfo[0].value= nm_per_step * offset
    parinfo[1].value= 1d
    ;parinfo[1].fixed= 1
    print,'Initial offset=',parinfo[0].value

    ;sy = sqrt(ABS(rsp.y))  ; Poisson errors
    sy = (abs(rsp.y) * 0.01d) > 1.d-8
    functargs = {spectx:sp.x, specty:sp.y}
    coeffs = mpfitfun('myfitfunc', rsp.x, rsp.y, sy, parinfo=parinfo, functargs=functargs, $
        quiet=quiet, status=status)

    ; update the spectra's wavelength according the refSpect
    new_x = poly(spect.ccdpos,coeffs)
    outspect=spect
    outspect.ccdpos = new_x

    ; convert ccdpos to wavelength and apply the prism temperature correction
    outspect.wavelength = interpol(profile_data.y4, profile_data.x, outspect.ccdpos,/spline)
    dndt = interpol(profile_data.y13, profile_data.y4, outspect.wavelength, /spline)
    outspect.wavelength = outspect.wavelength + (outspect.prismtemp - profile_data[0].y0) * dndt
 
    return, outspect

end
