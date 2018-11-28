;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Fit a SOLSTICE spectra to a reference SOLSTICE spectra by using the
;   AMOEBA routine to minimize the difference with the reference when 
;   adjusting the GRATING POSITION with a polynomial.  The routine
;   finds the coefficients of the best fit and returns the original spectra 
;   with the new wavelength scale.
;
; CALLING SEQUENCE:
;   spectra = ALIGN_SOL_SPECTRA(t0, t1, instrumentModeId, /missionDays)
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
;      7  FUV_A
;      9  MUV_A
;      11 FUV_B
;      13 MUV_B
;
; OPTIONAL INPUT PARAMETERS:
;   spect - 
;      Structure containing the spectra to align. If provided, the
;      startTime, stopTime and instrumentModeId are ignored.
;   refSpect - 
;      Structure containg the reference spectra to align to.
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
; REVISION HISTORY:
;   Revision: $Id: align_sol_spectra.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;------------------------------------------------------------------
;
function my_sol_amoeba, coeffs
    common mysol_amoeba_common, sx, sy, rx, ry, no_deriv, use_rms

    if n_elements(coeffs) gt 1 then $
        newx=poly(sx,coeffs) $
    else $
        newx=sx+coeffs[0]

    if no_deriv eq 0 then newy=deriv(newx,sy) else newy=sy
    y2=spl_init(newx, newy)
    newy=spl_interp(newx,newy, y2, rx)

    if use_rms eq 1 then begin
        corr_value = sqrt(total((newy - ry)^2d))
    endif else begin
        corr_value = abs(1d - abs(correlate(newy, ry, /double)))
    endelse

    ;print,coeffs,corr_value
    return, corr_value
end

;*******************************************************************

function align_sol_spectra, startTime, stopTime, instrumentModeId, $
         gps=gps, missionDays=missionDays, julianDays=julianDays, $
         refSpect=refSpect, spect=spect, $
         quiet=quiet, order=order, coeffs=coeffs, $
         noderiv=noderiv, status=status, rms=rms, $
		 verbose=verbose, goodness=goodness, wrange=wrange, ccdfit=ccdfit

    common mysol_amoeba_common

    if keyword_set(noderiv) then no_deriv=1 else no_deriv=0
    if keyword_set(rms) then use_rms=1 else use_rms=0

    if instrumentModeId lt 7 or instrumentModeId gt 13 then begin
        print,'Error: wrong instrumentModeId provided [7,9,11,13]'
        return,-1
    endif

    ; if a reference spectra was provided use that instead
    if size(refSpect,/tname) ne 'STRUCT' then begin
        case instrumentModeId of
            7: refSpectId=106
            9: refSpectId=116
            11: refSpectId=115
            13: refSpectId=117
            else: return,-1
        endcase
        query_database,'select * from ReferenceSpectrumData where referenceSpectrumId='+strtrim(string(refSpectId)),refData
        refSpect=replicate({wavelength:0d, ccdpos:0d, dn:0d},n_elements(refData))
        refSpect.dn = refData.intensity
        refSpect.wavelength=refData.abscissa
        refSpect.ccdpos=lambda2grt(refSpect.wavelength)
    endif

    ; get the spectra to fit to the reference
    if size(spect,/tname) ne 'STRUCT' then begin
        spect=get_solstice_spectra(startTime, stopTime, gps=gps, missionDays=missionDays, julianDays=julianDays, instrumentModeId)
        ; re-scale the data when either of the filter is in place (they are both ND1)
        p=where(spect.filter1_in eq 1,count)
        if count gt 0 then spect[p].dn*=10d
        p=where(spect.filter2_in eq 1,count)
        if count gt 0 then spect[p].dn*=10d
    endif

    ; the reference spectrum 10x the dispersion - smooth it out
    ; and resample it to the spectra to compare to
    y=interpol(smooth(refSpect.dn,10), refSpect.ccdpos, spect.ccdpos, /spline)
    rsp={x:spect.ccdpos, y:y, wavelength:grt2lambda(spect.ccdpos)}
    sp ={x:spect.ccdpos, y:spect.dn, wavelength:grt2lambda(spect.ccdpos)}

    ;if n_elements(order) eq 0 then order=1 else order=order>1

    ; get the global offset using cross_correlate (in number of steps)
    cross_correlate,rsp.y,sp.y,offset
    ; get the offset in wavelength
    nm_per_step = median(sp.x[1:-1]-sp.x[0:-2])

    if order eq 0 then begin
        coeffs = offset
        scale = 5d
    endif else begin
        coeffs = replicate(0d, order+1)
        coeffs[0] = offset
        coeffs[1] = 1d
        scale=coeffs
        scale = scale * 0d + 1d-2
        scale[0]=5d
        scale[1]=1d-2
    endelse
    if keyword_set(verbose) then print,'Initial offset=',coeffs

    s=sort(sp.x)
    q=uniq(sp.x[s])
    sx=sp.x[s[q]] 
    sy=sp.y[s[q]]
    s=sort(rsp.x)
    q=uniq(rsp.x[s])
    rx=rsp.x[s[q]]
    ry=rsp.y[s[q]]
    if no_deriv eq 0 then ry=deriv(rx, ry)

    coeffs = AMOEBA(1.0e-8, P0=coeffs, scale=scale, FUNCTION_VALUE=fval, FUNCTION_NAME='my_sol_amoeba')

    if keyword_set(verbose) then print,'Final coefficients = ',coeffs
    if keyword_set(verbose) then print,'Fit goodness=',1d - fval / 1d6

    ; update the spectra's wavelength according the refSpect
    if n_elements(coeffs) gt 1 then new_x = poly(spect.ccdpos,coeffs) else new_x=spect.ccdpos+coeffs[0]
    outspect=spect
    outspect.wavelength = grt2lambda(new_x)
    outspect.ccdpos=new_x
    ;delta_ccdpos = new_x - outspect.ccdpos_uncorr
    ;ccdfit = robust_poly_fit(outspect.ccdpos_uncorr, delta_ccdpos, order, /double)

	; define the goodness of the fit by the correlate
	goodness=0d
    goodness=correlate(deriv(outspect.wavelength,outspect.dn), deriv(rsp.wavelength,rsp.y))
 
    return, outspect

end



