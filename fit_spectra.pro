;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Fit a spectra in wavelength space to a reference spectra by measuring
;   the change in cross-correlation offset for a window along the spectra
;   and fitting a 3rd order polynomial fit to the spectra to bring it back
;   to the reference spectra wavelength space.
;
; CALLING SEQUENCE:
;   spectra = FIT_SPECTRA(t0, t1, instrumentModeId, /missionDays, refT0=t0, refT1=t1)
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
;   width -
;      Width of the portion of the spectra to match in CCDPos index number 
;      (not in subpixels). Defaults to 40 points.
;
; RETURNED PARAMETERS:
;   A structure with the spectra aligning the reference spectra in wavelength.
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
;   Revision: $Id: fit_spectra.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;-
;
function fit_spectra, startTime, stopTime, instrumentModeId, $
         gps=gps, missionDays=missionDays, julianDays=julianDays, $
         refT0=refT0, refT1=refT1, refSpect=refSpect, $
         version=version, width=width, $
         stepSize=stepSize, quiet=quiet

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

    ; if a reference spectra was provided use that instead
    if size(refSpect,/tname) ne 'STRUCT' then begin
        ; get the reference spectra
        refSpect = get_sim_spectra(refT0, refT1, instrumentModeId, gps=gps, $
            missionDays=missionDays, julianDays=julianDays, version=version, /noCcdCorr, /dn_only)
        if size(refSpect,/tname) ne 'STRUCT' then begin
            print,'No valid reference spectra found'
            return,-1
        endif
    endif
    derivRefSpect = deriv(refSpect.ccdpos, refSpect.dn)

    ; get the spectra to fit to the reference
    spect = get_sim_spectra(startTime, stopTime, instrumentModeId, gps=gps, $
        missionDays=missionDays, julianDays=julianDays, version=version, /noCcdCorr, /dn_only)
    if size(spect,/tname) ne 'STRUCT' then begin
        print,'No valid spectra found'
        return,-1
    endif
    
    ; interpolate the spect dn to the refspect ccdpos
    minRef=min(refSpect.ccdpos,max=maxRef)
    minSpect=min(spect.ccdpos,max=maxSpect)
    if minSpect lt minRef or maxSpect gt maxRef then begin
        ; interpolation will fail if trying to extrapolate
        print,'Warning: trying to extrapolate'
        ;return,-1
    endif
    newdn = interpol(spect.dn , spect.ccdpos, refSpect.ccdpos,/spline)
    ;newdn = spline(spect.ccdpos, spect.dn, refSpect.ccdpos)

    if NOT keyword_set(noDeriv) then derivSpect = deriv(refSpect.ccdpos, newdn)

    ; we'll calculate a new offset at every 10 points
    nslices = floor(float(n_elements(refspect.ccdpos)) / stepSize)
    npts = n_elements(refSpect.ccdpos)

    ; make sure the data spacing is constant and the same for both spectrum
    spacing = refSpect[1:-1].ccdpos - refSpect[0:-2].ccdpos
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

    outdata=replicate({meanCCDPos:0.0d, offset:0.0d, deriv_offset:0.0d, corr:0.0d, deriv_corr:0.0d}, nslices)

    p0=0L
    p1=p0+width
    deriv_offset=offset
    for slc=0L, nslices do begin
        if NOT keyword_set(quiet) then print,slc,nslices,format='("slice ",I5," / ",I5)'
        cross_correlate,refSpect[p0:p1].dn, newdn[p0:p1], offset, corr, ishift=offset
        ; result=c_correlate(refSpect[p0:p1].dn, dnspect[p0:p1],lag,/double)
        ; fit = gaussfit(lag,result,coeff,nterms=4)
        ; outdata[slc].offset = coeff[1]*refSpacing
        ; outdata[slc].corr=coeff[0]
        outdata[slc].meanCCDPos = mean(refSpect[p0:p1].ccdPos)
        outdata[slc].offset=offset*refSpacing
        outdata[slc].corr=max(corr)
        if NOT keyword_set(noDeriv) then begin
            cross_correlate,derivRefSpect[p0:p1], derivspect[p0:p1], deriv_offset, deriv_corr, ishift=deriv_offset
            ; result=c_correlate(derivRefSpect[p0:p1], derivspect[p0:p1], lag,/double)
            ; fit = gaussfit(lag,result,coeff,nterms=4)
            ; deriv_offset = coeff[1]*refSpacing
            ; outdata[slc].deriv_corr=coeff[0]
            outdata[slc].deriv_offset=deriv_offset*refSpacing
            outdata[slc].deriv_corr=max(deriv_corr)
        endif
        p0=p0+stepSize > 0L
        p1=p1+stepSize 
        if p1 gt npts then break
    endfor

    p=where(outdata.meanccdpos gt 0.0,count)
    if count eq 0 then return,-1
    outdata=outdata[p]

    ; get the 3rd order fit for the offset as a function of ccd position
    coeff = poly_fit(outdata.meanccdpos, outdata.deriv_offset, 2)
    ; update the spectra's wavelength according the refSpect
    ccdpos_corr = poly(spect.ccdpos,coeff)
    new_ccdpos = spect.ccdpos + ccdpos_corr
    ; get the corresponding wavelength by interpolating from the refSpect
    new_wavelength = interpol(refSpect.wavelength, refSpect.ccdpos, new_ccdpos,/spline)
    newspect = replicate(create_struct(spect[0],'ref_ccdpos',0.0d, 'ref_wavelength', 0.0d), n_elements(spect))
    newspect.timetag = spect.timetag
    newspect.ccdpos  = spect.ccdpos
    newspect.ccdpos_uncorr = spect.ccdpos_uncorr
    newspect.wavelength = spect.wavelength
    newspect.dn = spect.dn
    newspect.irradiance = spect.irradiance
    newspect.prismtemp = spect.prismtemp
    newspect.detectortemp = spect.detectortemp
    newspect.ref_ccdpos = new_ccdpos
    newspect.ref_wavelength = new_wavelength

    return, newspect

end
