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
;   Revision: $Id: align_spectra.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;------------------------------------------------------------------
;
function my_amoeba, coeffs
    common myamoeba_common, sx, sy, rx, ry, no_deriv

    if n_elements(coeffs) gt 1 then $
        newx=poly(sx,coeffs) $
    else $
        newx=poly(sx,[coeffs,1d])
    newy=interpol(sy,newx,rx,/lsq)

    if no_deriv eq 0 then newy = deriv(rx,deriv(rx,newy))

    corr_value = abs(1d - correlate(newy, ry, /double)) * 1d10
    ;print,coeffs,corr_value
    return, corr_value
end

;*******************************************************************

function align_spectra, startTime, stopTime, instrumentModeId, $
         gps=gps, missionDays=missionDays, julianDays=julianDays, $
         refT0=refT0, refT1=refT1, refSpect=refSpect, spect=spect, $
         version=version, quiet=quiet, order=order, coeffs=coeffs, $
         noderiv=noderiv, status=status, profile_data=profile_data, $
		 verbose=verbose, goodness=goodness, wrange=wrange, ccdfit=ccdfit

    common myamoeba_common

    if n_elements(version) eq 0 then version=20
    if keyword_set(noderiv) then no_deriv=1 else no_deriv=0

    ; if a reference spectra was provided use that instead
    if size(refSpect,/tname) ne 'STRUCT' then begin
        if n_elements(refT0) eq 0 or n_elements(refT1) eq 0 then begin
            ; this is our default REFERENCE SolarQuickScan24 spectra 
            if instrumentModeId eq 31 or instrumentModeId eq 32 then begin
                ; full ESR scan
                refT0 = 453.02065d
                refT1 = 453.98991d
            endif else begin
                refT0 = 453.67840d
                refT1 = 453.69524d
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
        refSpect = get_sim_spectra(refT0, refT1, instrumentModeId, gps=gps, missionDays=missionDays, $
            julianDays=julianDays, version=version, /noCcdCorr,/noTempCorr, /dn_only, /no1au, profile_data=profile_data)
        if size(refSpect,/tname) ne 'STRUCT' then begin
            if keyword_set(verbose) then print,'No valid reference spectra found'
            return,-1
        endif

    endif

    ; get the spectra to fit to the reference
    if size(spect,/tname) ne 'STRUCT' then begin
        spect = get_sim_spectra(startTime, stopTime, instrumentModeId, gps=gps, missionDays=missionDays, $
            julianDays=julianDays, version=version, /noCcdCorr, /noTempCorr, /dn_only, /no1au, profile_data=profile_data)
        if size(spect,/tname) ne 'STRUCT' then begin
            if keyword_set(verbose) then print,'No valid spectra found'
            return,-1
        endif
    endif

    ; re-map the refspect wavelength by adding the dLambda/dTemperature corresponding to the spect temperature
    match, refSpect.ccdpos, spect.ccdpos, suba, subb
    refWave = refSpect[suba].wavelength

    if size(profile_data,/tname) ne 'STRUCT' then begin
        rsp={x:refSpect[suba].ccdPos, y:refSpect[suba].dn}
        sp={x:spect[subb].ccdpos, y:spect[subb].dn}
    endif else begin
        ; get the wavelength correction at each spect temperature and apply to our refSpect
        dndt = interpol(profile_data.y18, profile_data.y4, refWave, /lsq);, /spline)
        ; since our call to get_sim_spectra for refSpect, corrects for the prism temp, the wavelengths are now
        ; for the prism reference temperature
        refWave = refWave - (spect[subb].prismtemp - profile_data[0].y0) * dndt
        ; now convert our refWave to ccdpos - corresponding to refSpect with same temperature profile as spect
        ccdPos = interpol(profile_data.x, profile_data.y4, refWave);, /spline)
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
    endelse

    ; now try to align the spect with the temperature modified refSpect

    if n_elements(order) eq 0 then order=1 else order=order>1
    ; get the global offset using cross_correlate (in number of steps)
    irrad=interpol(sp.y,sp.x,rsp.x,/lsq)
    cross_correlate,rsp.y,irrad,offset
    ; get the offset in wavelength
    nm_per_step = median(rsp.x[1:-1]-rsp.x[0:-2])

    if order eq 0 then begin
        coeffs = nm_per_step * offset
        scale = 80d
    endif else begin
        coeffs = replicate(0d, order+1)
        coeffs[0] = nm_per_step * offset
        coeffs[1] = 1d
        scale=coeffs
        scale = scale * 0d + 1d-2
        scale[0]=80d
        scale[1]=1d-2
    endelse
    if keyword_set(verbose) then print,'Initial offset=',coeffs

    sx=sp.x
    sy=sp.y
    rx=rsp.x
    if no_deriv eq 0 then ry=deriv(rsp.x, deriv(rsp.x, rsp.y)) else ry=rsp.y

    coeffs = AMOEBA(1.0e-8, P0=coeffs, scale=scale, FUNCTION_VALUE=fval, FUNCTION_NAME='my_amoeba')

    if keyword_set(verbose) then print,'Final coefficients = ',coeffs
    if keyword_set(verbose) then print,'Fit goodness=',1d - fval / 1d6

    ; update the spectra's wavelength according the refSpect
    new_x = poly(spect.ccdpos,coeffs)
    outspect=spect
    outspect.ccdpos = new_x
    delta_ccdpos = new_x - outspect.ccdpos_uncorr
    ccdfit = robust_poly_fit(outspect.ccdpos_uncorr, delta_ccdpos, order, /double)

    if n_elements(profile_data) eq 0 then begin
        outspect.wavelength = interpol(refSpect.wavelength, refSpect.ccdpos, outspect.ccdpos,/spline)
    endif else begin
        ; convert ccdpos to wavelength and apply the prism temperature correction
        outspect.wavelength = interpol(profile_data.y4, profile_data.x, outspect.ccdpos,/spline)
        dndt = interpol(profile_data.y18, profile_data.y4, outspect.wavelength, /spline)
        outspect.wavelength = outspect.wavelength + (outspect.prismtemp - profile_data[0].y0) * dndt * 1.2d
    endelse

	; define the goodness of the fit by the correlate
	goodness=0d
	if n_elements(outspect) eq n_elements(refspect) then begin
		newy = interpol(outspect.dn, outspect.wavelength, refspect.wavelength, /lsq)
		goodness=correlate(deriv(outspect.wavelength, deriv(outspect.wavelength,newy)), $
            deriv(refspect.wavelength, deriv(refspect.wavelength,refspect.dn)))
	endif
 
    return, outspect

end



