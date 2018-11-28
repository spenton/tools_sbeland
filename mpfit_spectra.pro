;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Fit a spectra in wavelength space to a reference spectra by using the
;   MPFITFUN routine to minimize the difference with the reference when 
;   adjusting the wavelength scale with a polynomial.  The routine
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
;   IDL> res=mpfit_spectra(T0,T1,47,/mission,refspect=res43,coeffs=coeffs,status=status)
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
; REVISION HISTORY:
;   Revision: $Id: mpfit_spectra.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;------------------------------------------------------------------
;
;  Function to apply a polynomial to the wavelength scale of the spectra to compare to
function myfitfunc, x, coeffs, spectx=spectx, specty=specty
    newx=poly(spectx,coeffs)
    newy=interpol(specty,newx,x,/lsq)
    return, newy
end

;*******************************************************************

function mpfit_spectra, startTime, stopTime, instrumentModeId, $
         gps=gps, missionDays=missionDays, julianDays=julianDays, $
         refT0=refT0, refT1=refT1, refSpect=refSpect, spect=spect, $
         version=version, quiet=quiet, order=order, coeffs=coeffs, $
         noderiv=noderiv, status=status

    ; if a reference spectra was provided use that instead
    if size(refSpect,/tname) ne 'STRUCT' then begin
        if n_elements(refT0) eq 0 or n_elements(refT1) eq 0 then begin
            ; this is our default REFERENCE SolarQuickScan24 spectra 
            if instrumentModeId eq 31 or instrumentModeId eq 32 then begin
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
        refSpect=get_science_product(['SimProfileIntegral','SimCalibratedIrradiance'],refT0, refT1,$
            instrumentModeId, mission=missionDays, gps=gps, julian=julianDays, version=version)
        if size(refSpect,/tname) ne 'STRUCT' then begin
            print,'No valid reference spectra found'
            return,-1
        endif

    endif

    ; get the spectra to fit to the reference
    if size(spect,/tname) ne 'STRUCT' then begin
        spect=get_science_product(['SimProfileIntegral','SimCalibratedIrradiance'],startTime, stopTime, $
            instrumentModeId, /mission, version=version)
        if size(spect,/tname) ne 'STRUCT' then begin
            print,'No valid spectra found'
            return,-1
        endif
    endif
    
    tnames = tag_names(refSpect)
    refx=where(strpos(tnames,'WAVELENGTH') ge 0,count)
    if count eq 0 then begin
        print,'Error: refSpect does not have a WAVELENGTH tag'
        return,-1
    endif else if count gt 1 then begin
        p=where(strpos(tnames,'WAVELENGTHREF') ge 0,count)
        if count eq 1 then refx=p else refx=refx[0]
    endif
    refy=where(strpos(tnames,'DN') ge 0,count)
    if count eq 0 then begin
        refy=where(strpos(tnames,'IRRADIANCE') ge 0,count)
        if count eq 0 then begin
            print,'Error: refSpect does not have a DN or IRRADIANCE tag'
            return,-i
        endif
    endif else if count gt 1 then begin
        p = where(tnames eq 'DN',count)
        if count eq 1 then refy=p else refy=refy[0]
    endif
    s=sort(refSpect.(refx))
    if n_elements(refSpect) eq 1 then begin
        q=uniq(refSpect.(refx)[s])
        rsp={wavelength:refSpect.(refx)[s[q]], irradiance:refSpect.(refy)[s[q]]}
    endif else begin
        q=uniq(refSpect[s].(refx))
        rsp={wavelength:refSpect[s[q]].(refx), irradiance:refSpect[s[q]].(refy)}
    endelse

    tnames = tag_names(spect)
    spx=where(strpos(tnames,'WAVELENGTH') ge 0,count)
    if count eq 0 then begin
        print,'Error: Spect does not have a WAVELENGTH tag'
        return,-1
    endif else if count gt 1 then begin
        p=where(strpos(tnames,'WAVELENGTHREF') ge 0,count)
        if count eq 1 then spx=p else spx=spx[0]
    endif
    spy=where(strpos(tnames,'DN') ge 0,count)
    if count eq 0 then begin
        spy=where(strpos(tnames,'IRRADIANCE') ge 0,count)
        if count eq 0 then begin
            print,'Error: Spect does not have a DN or IRRADIANCE tag'
            return,-1
        endif
    endif else if count gt 1 then begin
        p = where(tnames eq 'DN',count)
        if count eq 1 then spy=p else spy=spy[0]
    endif
    s=sort(spect.(spx))
    if n_elements(spect) eq 1 then begin
        q=uniq(spect.(spx)[s])
        sp={wavelength:spect.(spx)[s[q]], irradiance:spect.(spy)[s[q]]}
    endif else begin
        q=uniq(spect[s].(spx))
        sp={wavelength:spect[s[q]].(spx), irradiance:spect[s[q]].(spy)}
    endelse

    ; limit the 2 spectra to the overlapping area
    rmin=min(rsp.wavelength,max=rmax)
    smin=min(sp.wavelength,max=smax)
    smin=max([smin,rmin])
    smax=min([smax,rmax])
    p=where(rsp.wavelength ge smin and rsp.wavelength le smax,count)
    if count eq 0 then begin
        print,'Error: not overlapping area between the 2 spectrum'
        return,-1
    endif
    rsp={wavelength:rsp.wavelength[p], irradiance:rsp.irradiance[p]}

    p=where(sp.wavelength ge smin and sp.wavelength le smax,count)
    if count eq 0 then begin
        print,'Error: not overlapping area between the 2 spectrum'
        return,-1
    endif
    sp={wavelength:sp.wavelength[p], irradiance:sp.irradiance[p]}
    sp_wave=sp.wavelength
    sp_irrad=sp.irradiance
    rsp_irrad=rsp.irradiance

    if NOT keyword_set(noderiv) then begin
        ; compare the derivatives instead of the actual spectrum
        rsp.irradiance = deriv(rsp.wavelength,deriv(rsp.wavelength,rsp.irradiance))
        sp.irradiance = deriv(sp.wavelength,deriv(sp.wavelength,sp.irradiance))
        ; clean up the derivative with a sliding window and rejecting
        ; points with stdev larger than 9 sigma
        ; break up the wavelength region in two since the first part 
        ; will have a much higher sigma than the longer wavelengths
        keep0=[]
        keep1=[]
        p0=where(rsp.wavelength lt 1000.0,count0,comp=p1)
        if p0[0] ne -1 then resistant_mean,rsp.irradiance[p0],9.0,mean,goodvec=keep0
        if p1[0] ne -1 then resistant_mean,rsp.irradiance[p1],5.0,mean,goodvec=keep1
        keep = [p0[keep0],p1[keep1]]
        if n_elements(keep) gt n_elements(rsp.irradiance)/2 then begin
            temp={wavelength:rsp.wavelength[keep], irradiance:rsp.irradiance[keep]}
            rsp=temporary(temp)
            rsp_irrad=rsp_irrad[keep]
        endif
        rsp.wavelength=rsp.wavelength[keep]
        rsp.irradiance=rsp.irradiance[keep]
        keep0=[]
        keep1=[]
        p0=where(sp.wavelength lt 1000.0,count0,comp=p1)
        if p0[0] ne -1 then resistant_mean,sp.irradiance[p0],9.0,mean,goodvec=keep0
        if p1[0] ne -1 then resistant_mean,sp.irradiance[p1],5.0,mean,goodvec=keep1
        keep = [p0[keep0],p1[keep1]]
        if n_elements(keep) gt n_elements(sp.irradiance)/2 then begin
            temp={wavelength:sp.wavelength[keep], irradiance:sp.irradiance[keep]}
            sp=temporary(temp)
        endif
        sp.wavelength=sp.wavelength[keep]
        sp.irradiance=sp.irradiance[keep]
    endif

    if n_elements(order) eq 0 then order=3
    coeff0=dblarr(order+1)
    ; get the global offset using cross_correlate (in number of steps)
    irrad=interpol(sp.irradiance,sp.wavelength,rsp.wavelength,/lsq)
    cross_correlate,rsp.irradiance,irrad,offset
    ; get the offset in wavelength
    nm_per_step = median(rsp.wavelength[1:-1]-rsp.wavelength[0:-2])/2d
    coeff0[0]=-1d * nm_per_step * offset
    coeff0[1]=1.0d
    ;sy = sqrt(ABS(rsp.irradiance))  ; Poisson errors
    sy = 0.01d*rsp.irradiance
    functargs = {spectx:sp.wavelength, specty:sp.irradiance}
    coeffs = mpfitfun('myfitfunc', rsp.wavelength, rsp.irradiance, sy, coeff0, functargs=functargs, $
        quiet=quiet, status=status, xtol=1d-15, mpprint=0)

    ; update the spectra's wavelength according the refSpect
    new_wavelength = poly(sp_wave,coeffs)

    return, {wavelength:new_wavelength, irradiance:sp_irrad}

end
