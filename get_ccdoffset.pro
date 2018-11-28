;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Return a structure containing the CCD Spot position offset from a
;   reference spectra as a function of CCD spot position for a specified 
;   instrumentModeId and time span.
;
; CALLING SEQUENCE:
;   spectra = GET_CCDOFFSET(t0, t1, instrumentModeId, /missionDays)
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
;
; OPTIONAL INPUT PARAMETERS:
;   notempCorr - 
;      If specified, will skip the prism temperature correction when 
;      assigning a wavelength, returning the interpolated value only.
;   noDark - 
;      If specified, will skip the diode dark value correction
;   noDeriv -
;      If specified, will skip the cross-correlation fit with the 
;      derivative of the spectra. By default, does a 2nd order DERIV.
;   refT0 -
;      StartTime of the reference spectra in the same units as T0.
;      If not specified, uses our previously defined time range.
;   refT1 -
;      StopTime of the reference spectra in the same units as T1.
;      If not specified, uses our previously defined time range.
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
;   width -
;      Width of the portion of the spectra to match in CCDPos index number 
;      (not in subpixels). Defaults to 40 points.
;   dowave -
;      Instead of calculating the offsets in CCDPos, do it in wavelength space,
;      wnd then convert back to pixel space (taking care of prismTempCorr)
;
; RETURNED PARAMETERS:
;   A structure with the wavelength and irradiance (in raw dn).
;
; OPTIONAL OUTPUT PARAMETERS:
;      NONE
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
; REVISION HISTORY:
;   Revision: $Id: get_ccdoffset.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;-
;
function get_ccdoffset, startTime, stopTime, instrumentModeId, $
         gps=gps, missionDays=missionDays, julianDays=julianDays, $
         noTempCorr=noTempCorr, noDark=noDark, version=version, $
         refT0=refT0, refT1=refT1, refSpect=refSpect, spect=spect, $
         noDeriv=noDeriv, width=width, noCcdCorr=noCcdCorr, $
         stepSize=stepSize, fitcoeff=fitcoeff, quiet=quiet, dowave=dowave, $
         profile_data=profile_data

    if n_elements(refT0) eq 0 or n_elements(refT1) eq 0 then begin
        ; this is our default REFERENCE SolarQuickScan24 spectra 
        refT0 = 453.67840d
        refT1 = 453.69524d
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

    if n_elements(width) eq 0 then width=40.0d else width=double(width)
    if n_elements(stepSize) eq 0 then stepSize=10.0d else stepSize=double(stepSize)
    ; force the prism temperature correction if dowave is set
    if keyword_set(dowave) then noTempCorr=0

    ; if a reference spectra was provided use that instead
    if size(refSpect,/tname) ne 'STRUCT' then begin
        ; get the reference spectra
        refSpect = get_sim_spectra(refT0, refT1, instrumentModeId, gps=gps, $
            missionDays=missionDays, julianDays=julianDays, /uncorr, $
            noDark=noDark, version=version, /noCcdCorr, /dn_only, /no1au, profile_data=profile_data)
        if size(refSpect,/tname) ne 'STRUCT' then begin
            print,'No valid reference spectra found'
            return,-1
        endif
        ; some of the early scans have unsorted ccdpos at the beginning
        s=sort(refspect.ccdpos)
        q=uniq(refspect.ccdpos,s)
        refspect=refspect[s[q]]
    endif
    if NOT keyword_set(noDeriv) then begin
        if NOT keyword_set(dowave) then derivRefSpect=deriv(refSpect.ccdpos, refSpect.dn) else $
            derivRefSpect=deriv(refSpect.wavelength, refSpect.dn)
    endif

    if size(spect,/tname) ne 'STRUCT' then begin
        ; get the spectra to fit to the reference
        spect = get_sim_spectra(startTime, stopTime, instrumentModeId, gps=gps, missionDays=missionDays, $
            julianDays=julianDays, noDark=noDark, version=version, fitcoeff=fitcoeff, $
            noCcdCorr=noCcdCorr, /uncorr, /no1au, profile_data=profile_data)
    endif
    if size(spect,/tname) ne 'STRUCT' then begin
        print,'No valid spectra found'
        return,-1
    endif
    ; some of the early scans have unsorted ccdpos at the beginning
    s=sort(spect.ccdpos)
    q=uniq(spect.ccdpos,s)
    spect=spect[s[q]]
    
    ; interpolate the spect dn to the refspect ccdpos
    minRef=min(refSpect.ccdpos,max=maxRef)
    minSpect=min(spect.ccdpos,max=maxSpect)
    if minSpect lt minRef or maxSpect gt maxRef then begin
        ; interpolation will fail if trying to extrapolate
        print,'Warning: trying to extrapolate'
        ;return,-1
    endif
    
    if keyword_set(dowave) then begin
        newdn = interpol(spect.dn , spect.wavelength, refSpect.wavelength,/spline)
        prismTemp = interpol(spect.prismTemp, spect.wavelength, refSpect.wavelength,/spline)
    endif else begin
        ;newdn = spline(spect.ccdpos, spect.dn, refSpect.ccdpos)
        newdn = interpol(spect.dn , spect.ccdpos, refSpect.ccdpos,/spline)
        prismTemp = interpol(spect.prismTemp, spect.ccdpos, refSpect.ccdpos,/spline)
    endelse

    if NOT keyword_set(noDeriv) then begin
        if NOT keyword_set(dowave) then derivSpect=deriv(refSpect.ccdpos, newdn) else $
            derivSpect=DERIV(refSpect.wavelength, newdn)
    endif 

    ; we'll calculate a new offset at every 10 points
    nslices = floor(float(n_elements(refspect.ccdpos)) / stepSize)
    npts = n_elements(refSpect.ccdpos)

    ; make sure the data spacing is constant and the same for both spectrum
    if NOT keyword_set(dowave) then $
        spacing = refSpect[1:-1].ccdpos - refSpect[0:-2].ccdpos $
    else $
        spacing = median(refSpect[1:-1].wavelength - refSpect[0:-2].wavelength)
    minSpc = min(spacing,max=maxSpc)
    if minSpc ne maxSpc then begin
        ; the spacing should be constant for the whole data sample
        print,'Warning: the spacing in the refSpect is not uniform !!'
        return,-1
    endif
    refSpacing=minSpc

    ; start by calculating an initial offset for first 1/4 of spectra
    cross_correlate, refSpect[0:long(npts/4)].dn, newdn[0:long(npts/4)], offset, corr
    ; lag = dindgen(200.0)/10.0d - 10.0d
    ; result=c_correlate(refSpect[0:long(npts/4)].dn, spect[0:long(npts/4)].dn,lag,/double)
    ; perform a gaussfit on the correlation result to find best offset with Lag
    ; fit = gaussfit(lag,result,coeff,nterms=4)
    ; offset = coeff[1]

    if NOT keyword_set(quiet) then print,'Initial OFFSET=',offset*refSpacing
    prismRefTemp=median(profile_data.y0)

    outdata=replicate({meanCCDPos:0.0d, offset:0.0d, deriv_offset:0.0d, corr:0.0d, deriv_corr:0.0d}, nslices)

    p0=0L
    p1=p0+width
    deriv_offset=offset
    for slc=0L, nslices-1L do begin
        if NOT keyword_set(quiet) then print,slc,nslices,format='("slice ",I5," / ",I5)'
        cross_correlate,newdn[p0:p1], refSpect[p0:p1].dn, offset, corr, ishift=offset
        ; result=c_correlate(refSpect[p0:p1].dn, dnspect[p0:p1],lag,/double)
        ; fit = gaussfit(lag,result,coeff,nterms=4)
        ; outdata[slc].offset = coeff[1]*refSpacing
        ; outdata[slc].corr=coeff[0]
        outdata[slc].meanCCDPos = mean(refSpect[p0:p1].ccdPos)
        outdata[slc].corr=max(corr)
        if NOT keyword_set(dowave) then begin
            refSpacing = median(refSpect[p0+1:p1].ccdpos - refSpect[p0:p1-1].ccdpos)
            outdata[slc].offset=offset*refSpacing
        endif else begin
            refSpacing = median(refSpect[p0+1:p1].wavelength - refSpect[p0:p1-1].wavelength)
            ; convert our wavelength offset to pixel offset
            meanRefWave = mean(refSpect[p0:p1].wavelength)
            meanWave = meanRefWave + offset*refSpacing
            ; remove the prismtempCorr for this wavelength
            ;outdata.wavelength = outdata.wavelength + (prismTemp - prismRefTemp) * dndt
            dndt = interpol(profile_data.y13,profile_data.y4,meanWave)
            meanWave -= (median(prismTemp[p0:p1]) - prismRefTemp) * dndt
            ; convert this wavelength to ccd pixel position
            outdata[slc].offset = interpol(profile_data.x,profile_data.y4,meanWave,/spline) - $
                                  interpol(profile_data.x,profile_data.y4,meanRefWave,/spline)
        endelse
        if NOT keyword_set(noDeriv) then begin
            cross_correlate,derivspect[p0:p1], derivRefSpect[p0:p1], deriv_offset, deriv_corr, ishift=deriv_offset
            ; result=c_correlate(derivRefSpect[p0:p1], derivspect[p0:p1], lag,/double)
            ; fit = gaussfit(lag,result,coeff,nterms=4)
            ; deriv_offset = coeff[1]*refSpacing
            ; outdata[slc].deriv_corr=coeff[0]
            outdata[slc].deriv_corr=max(deriv_corr)
            if NOT keyword_set(dowave) then begin
                refSpacing = median(refSpect[p0+1:p1].ccdpos - refSpect[p0:p1-1].ccdpos)
                outdata[slc].deriv_offset=deriv_offset*refSpacing
            endif else begin
                refSpacing = median(refSpect[p0+1:p1].wavelength - refSpect[p0:p1-1].wavelength)
                ; convert our wavelength offset to pixel offset
                meanRefWave = mean(refSpect[p0:p1].wavelength)
                meanWave = meanRefWave + offset*refSpacing
                ; remove the prismtempCorr for this wavelength
                ;outdata.wavelength = outdata.wavelength + (prismTemp - prismRefTemp) * dndt
                dndt = interpol(profile_data.y13,profile_data.y4,meanWave)
                meanWave -= (median(prismTemp[p0:p1]) - prismRefTemp) * dndt
                ; convert this wavelength to ccd pixel position
                outdata[slc].deriv_offset = interpol(profile_data.x,profile_data.y4,meanWave,/spline) - $
                                      interpol(profile_data.x,profile_data.y4,meanRefWave,/spline)
            endelse
        endif
        p0=p0+stepSize > 0L
        p1=p1+stepSize < (n_elements(refSpect)-1)
        if p1 gt npts then break
    endfor

    p=where(outdata.meanccdpos gt 0.0,count)
    if count eq 0 then return,-1
    return, outdata[p]

end
