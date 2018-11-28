;+
; Author: Stephane Beland
;
;   STILL A WORK IN PROGRESS:  WORKS FOR SOME DATA BUT NOT WITH OTHER
;   -----------------------------------------------------------------
;
; REVISION HISTORY:
;-
;*******************************************************************
function align_spectra4, x0,y0,x1,y1,order=order,nsegments=nsegments

    s=sort(x0)
    nx0=x0[s]
    ny0=y0[s]

    ; reorder the reference data in linear spacing
    mn=min(nx0,max=mx)
    inc=min(abs(nx0[1:-1] - nx0[0:-2]))
    dim = CEIL((mx-mn)/inc)
    new_x0=make_array(increment=inc,start=mn, /double,dim=dim)
    ; perform a linear interpolation to the new grid
    new_y0=interpol(ny0,nx0,new_x0)

    new_x1=new_x0
    new_y1=interpol(y1, x1, new_x1)

    if n_elements(order) eq 0 then order=3
    if n_elements(nsegments) eq 0 then nsegments=10

    npts = FLOOR(n_elements(x0) / nsegments)

    ;get a global shift first
    cross_correlate, new_y0, new_y1, offset, corr, ishift=0d, width=n_elements(new_y0)/2

    p0=0L
    offset=dblarr(nsegments+1)
    for i=0,nsegments do begin
        p1=(p0L + npts -1L) < (n_elements(x0)-1L)
        cross_correlate, y0[p0:p1], y1[p0:p1], offset[i], corr
    endfor


    if n_elements(version) eq 0 then version=22
	if n_elements(order) eq 0 then order=3
	if n_elements(stepsize) eq 0 then stepsize=10d else stepsize=double(stepsize)
    goodness=0d
    ccdfit=0d

    valid_modes=[31,32,41,43,44,45,47,48]
    match,instrumentModeId,valid_modes,sa,sb
    if sa[0] eq -1 then return,-1 else iMode=instrumentModeId


    result = get_ccdoffset(startTime, stopTime, instrumentModeId, $
         gps=gps, missionDays=missionDays, julianDays=julianDays, $
         refT0=refT0, refT1=refT1, refSpect=refSpect, spect=spect, $
         noTempCorr=noTempCorr, version=version, quiet=quiet, stepsize=stepsize, $
         noDeriv=noDeriv, /noCcdCorr, /noDark, profile_data=profile_data)
         

    ; if not able to fit simply return
    if size(result,/tname) ne 'STRUCT' then return,-1

    ; prepare the output aligned spectra
    outspect = spect
    if keyword_set(noderiv) then begin
        offset = result.offset
        goodness = meadian(result.corr)
    endif else begin
        offset = result.deriv_offset
        goodness = median(result.deriv_corr)
    endelse
    ccdfit0=robust_poly_fit(result.meanccdpos, offset, order, /double)
    outspect.ccdpos = outspect.ccdpos_uncorr - poly(outspect.ccdpos_uncorr, ccdfit0)

    ; get the corresponding wavelength using the instrument model and the reference spectra median prism temperature
    ; since it was align to the reference spectra
    outspect.wavelength = ccd2lambda(instrumentModeId, outspect.ccdpos, median(refspect.prismtemp))

    ; now get the ccd positions for these wavelengths at the corresponding prism temperature
    outspect.ccdpos = lambda2ccd(instrumentModeId, outspect.wavelength, outspect.prismtemp)
    
    ; and get the ccdfit prior to the prism temperature correction so we can use it in the production code
    ccdfit = robust_poly_fit(outspect.ccdpos, outspect.ccdpos_uncorr - outspect.ccdpos, order, /double)

    if ~keyword_set(no_plot) then begin
        title='Mode='+strtrim(string(instrumentModeId),2)+' Align_Spectra4 on '+strtrim(string(startTime,format='(F0.2)'),2)
        plot_multi,result.meanccdpos, offset, result.meanccdpos, poly(result.meanccdpos, ccdfit0), $
            /xst,/yst,xtitle='CCD Position (subpixels)',ytitle='Delta CCDPOS',title=title, psym=[4,-3], $
            label=['Delta CCD Pos','Polynomial Fit'], charsize=1.4, thickness=[1,3]
    endif

    return, outspect

end

