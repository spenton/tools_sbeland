;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Compare a spectra with the reference spectra by matching the 
;   corresponding zero-crossing of the first and seconds derivatives.
;
; CALLING SEQUENCE:
;   list = MATCH_DISPERSION(instrumentModeId, refZero, t0=t0, t1=t1, spect=spect, /missionDays)
;
; INPUT PARAMETERS:
;   refZero - 
;      Structure containing the list of zero-crossing points to match.
;      Comes a previous run of find_ref_peaks.pro
;   instrumentModeId -
;      The instrument mode of interest are currently limited to:
;      41	SIM_A	VIS1
;      43	SIM_A	UV
;      44	SIM_A	IR
;      45	SIM_B	VIS1
;      47	SIM_B	UV
;
; OPTIONAL INPUT PARAMETERS:
;   T0 -
;      StartTime of the spectra to match (if SPECT is not provided).
;   T1 -
;      StopTime of the spectra to match (if SPECT is not provided).
;   spect -
;      Spectrum to match (structure returned from get_sim_spectra)
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
;   d1height -
;      Minimum height the derivative has to increase for a valey or peak
;      to be considered significant.  Any bumps smaller than that will be
;      ignored.
;   d2height -
;      Minimum height the 2nd derivative has to increase for a valey or peak
;      to be considered significant.  Any bumps smaller than that will be
;      ignored.
;
; RETURNED PARAMETERS:
;   A structure with the list of peaks and valeys. 
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
;  The default reference spectra if refT0 and refT1 are not provided are
;  from:  453.67840, 453.69524 which was validated
;
; REVISION HISTORY:
;   Revision: $Id: match_dispersion.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;-
;
function match_dispersion, instrumentModeId, refZero, t0=t0, t1=t1, $
         spect=spect, refSpect=refSpect, refT0=refT0, refT1=refT1, $
         gps=gps, missionDays=missionDays, julianDays=julianDays, $
         version=version, quiet=quiet, d1height=d1height, d2height=d2height, $
         noplot=noplot, coeff=coeff

    ; get the spectrum to match
    if size(spect,/tname) ne 'STRUCT' and n_elements(t0) eq 1 and n_elements(t1) eq 1 then begin
        ; the spectra to match was not provided -> go get it
        spect = get_sim_spectra(T0, T1, instrumentModeId, gps=gps, julianDays=julianDays, $
                missionDays=missionDays, version=version, /noCcdCorr, /dn_only)
    endif

    if size(spect,/tname) ne 'STRUCT' then begin
        print,'No valid reference spectra found'
        return,-1
    endif

    if n_elements(d1height) eq 0 and n_elements(d2height) eq 0 then begin
        case instrumentModeId of
            41: begin
                ; SIMA VIS1
                d1height = 0.002
                d2height = 0.002
                d2range=[488.0d,980.0d]
            end
            42: begin
                ; SIMA VIS2
                d1height = 0.002
                d2height = 0.002
                d2range=[488.0d,980.0d]
            end
            43: begin
                ; SIMA UV
                d1height = 0.02
                d2height = 0.06
                d2range=[0.0d,260.0d]
            end
            44: begin
                ; SIMA IR
                d1height = 0.0002
                d2height = 0.0002
                d2range=[0.0d,1700.0d]
            end
            45: begin
                ; SIMB VIS1
                d1height = 0.002
                d2height = 0.002
                d2range=[488.0d,980.0d]
            end
            46: begin
                ; SIMB VIS2
                d1height = 0.002
                d2height = 0.002
                d2range=[488.0d,980.0d]
            end
            47: begin
                ; SIMB UV
                d1height = 0.02
                d2height = 0.12
                d2range=[0.0d,260.0d]
            end
            48: begin
                ; SIMB IR
                d1height = 00.002
                d2height = 0.0002
                d2range=[0.0d,1700.0d]
            end
        endcase
    endif

    ; get the list of 1st and 2nd order derivatives from our reference spectra
    rz1 = where(refZero.deriv_order eq 1,rz1_count)
    rz2 = where(refZero.deriv_order eq 2,rz2_count)
    ; form the list of distances to 3 neighbors before and after for node matching
    dist_rz1 = replicate({distance:dblarr(4)},rz1_count)
    for i=0,rz1_count-1 do begin
        for j=-2,-1 do if i+j ge 0 then $
            dist_rz1[i].distance[j+2]=refZero[rz1[i]].wavelength - refZero[rz1[(i+j)>0]].wavelength
        for j=1,2 do if i+j lt rz1_count then $
            dist_rz1[i].distance[j+1]=refZero[rz1[i]].wavelength - refZero[rz1[(i+j)<(rz1_count-1)]].wavelength
    endfor

    dist_rz2 = replicate({distance:dblarr(4)},rz2_count)
    for i=0,rz2_count-1 do begin
        for j=-2,-1 do if i+j ge 0 then $
            dist_rz2[i].distance[j+2]=refZero[rz2[i]].wavelength - refZero[rz2[(i+j)>0]].wavelength
        for j=1,2 do if i+j lt rz2_count then $
            dist_rz2[i].distance[j+1]=refZero[rz2[i]].wavelength - refZero[rz2[(i+j)<(rz2_count-1)]].wavelength
    endfor

    zero_cross = find_ref_peaks(instrumentModeId, refSpect=spect, version=version, $
                 d1height=d1height, d2height=d2height)

    if size(zero_cross,/tname) ne 'STRUCT' then return,-1

    zc1 = where(zero_cross.deriv_order eq 1,zc1_count)
    zc2 = where(zero_cross.deriv_order eq 2 and zero_cross.wavelength ge d2range[0] and zero_cross.wavelength le d2range[1],zc2_count)
    dist_zc1 = replicate({distance:dblarr(4)},zc1_count)
    for i=0,zc1_count-1 do begin
        for j=-2,-1 do if i+j ge 0 then $
            dist_zc1[i].distance[j+2]=zero_cross[zc1[i]].wavelength - zero_cross[zc1[(i+j)>0]].wavelength
        for j=1,2 do if i+j lt rz1_count then $
            dist_zc1[i].distance[j+1]=zero_cross[zc1[i]].wavelength - zero_cross[zc1[(i+j)<(zc1_count-1)]].wavelength
    endfor

    ;dist_zc2 = replicate({distance:dblarr(4)},zc2_count)
    ;for i=0,zc2_count-1 do begin
    ;    for j=-2,-1 do if i+j ge 0 then $
    ;        dist_zc2[i].distance[j+2]=zero_cross[zc2[i]].wavelength - zero_cross[zc2[(i+j)>0]].wavelength
    ;    for j=1,2 do if i+j lt rz2_count then $
    ;        dist_zc2[i].distance[j+1]=zero_cross[zc2[i]].wavelength - zero_cross[zc2[(i+j)<(zc2_count-1)]].wavelength
    ;endfor

    ; define the output structure
    outdata = replicate({wavelength:0.0d, dn:0.0d, deriv_order:0, wave_match:0.0d, distance_sum:0.0d},n_elements(refZero))
    outdata.wavelength = refZero.wavelength
    outdata.dn = refZero.dn
    outdata.deriv_order = refZero.deriv_order

    ; first loop through the list of 1st derivative zero-crossings
    for i=0,rz1_count-1 do begin
        ; find the corresponding point from the spacing with adjacent peaks
        ; look for the closest match in wavelength and test the nearest 4 peaks
        delta = ABS(refZero[rz1[i]].wavelength - zero_cross[zc1].wavelength)
        mn=min(delta,pos)
        tst_pos = (([pos-4, pos-3, pos-2, pos-1, pos, pos+1, pos+2, pos+3, pos+4]) > 0) < (zc1_count-1L)
        tst_pos=tst_pos[uniq(tst_pos)]
        distance_sum = dblarr(n_elements(tst_pos))
        for j=0,n_elements(tst_pos)-1 do begin
            ; ignore the points where the distance is 0.0
            p=where(dist_zc1[tst_pos[j]].distance ne 0.0 and dist_rz1[i].distance ne 0.0)
            distance=ABS(dist_zc1[tst_pos[j]].distance - dist_rz1[i].distance)
            distance_sum[j]=TOTAL(distance)
        endfor
        mn = min(distance_sum,pos)
        outdata[rz1[i]].wave_match = zero_cross[zc1[tst_pos[pos]]].wavelength
        outdata[rz1[i]].distance_sum = mn
    endfor

    ; get the initial polynomial fit from the first derivative results
    order = n_elements(rz1)-2L < 3
    coeff=dblarr(order+1)
    coeff = robust_poly_fit(outdata[rz1].wavelength, outdata[rz1].wave_match, order, /double)

    ; we should get an estimate of the new wavelengths but for now
    ; simply use the reference spectra positions
    new_refZero2_wave = refZero[rz2].wavelength
    ;new_refZero2_wave = poly(refZero[rz2].wavelength, coeff)
    for i=0,rz2_count-1 do begin
        ; look for closest point to expected position
        delta = ABS(new_refZero2_wave[i] - zero_cross[zc2].wavelength)
        mn=min(delta,pos)
        ; the corresponding point should be about 1/10 of spacing between position
        p=where(spect.wavelength gt zero_cross[zc2[pos]].wavelength)
        if p[0] eq n_elements(spect)-1 then p[0]=p[0]-1L
        ;if mn gt ABS(spect[p[0]+1].wavelength-spect[p[0]].wavelength)/4.0 then continue
        outdata[rz2[i]].wave_match = zero_cross[zc2[pos]].wavelength
        outdata[rz2[i]].distance_sum = mn
    endfor


    ; loop through the list of 2nd derivative zero-crossings
    ; for i=0,rz2_count-1 do begin
    ;     ; find the corresponding point from the spacing with adjacent peaks
    ;     ; look for the closest match in wavelength and test the nearest 4 peaks
    ;     delta = ABS(refZero[rz2[i]].wavelength - zero_cross[zc2].wavelength)
    ;     mn=min(delta,pos)
    ;     ; if off by more than 5nm skip to next one
    ;     if mn gt 5.0 then continue
    ;     tst_pos = (([pos-4, pos-3, pos-2, pos-1, pos, pos+1, pos+2, pos+3, pos+4]) > 0) < (zc2_count-1L)
    ;     tst_pos=tst_pos[uniq(tst_pos)]
    ;     distance_sum = dblarr(n_elements(tst_pos))
    ;     for j=0,n_elements(tst_pos)-1 do begin
    ;         ; ignore the points where the distance is 0.0
    ;         p=where(dist_zc2[tst_pos[j]].distance ne 0.0 and dist_rz2[i].distance ne 0.0)
    ;         distance=ABS(dist_zc2[tst_pos[j]].distance[p] - dist_rz2[i].distance[p])
    ;         distance_sum[j]=TOTAL(distance)
    ;     endfor
    ;     mn = min(distance_sum,pos)
    ;     outdata[rz2[i]].wave_match = zero_cross[zc2[tst_pos[pos]]].wavelength
    ;     outdata[rz2[i]].distance_sum = mn
    ; endfor

    ; get the fit and plot reference and test spectra with the residuals
    p=where(outdata.wave_match gt 0.0 and outdata.distance_sum lt 1.0,count)
    if count lt 2 then begin
        print,'No matching peaks found - skip plot'
        return, outdata
    endif
    order = n_elements(p)-2L < 3
    coeff=dblarr(order+1)
    coeff = robust_poly_fit(outdata[p].wave_match, outdata[p].wavelength, order, /double)
    yfit = poly(outdata[p].wave_match, coeff)
    residuals = outdata[p].wavelength - poly(outdata[p].wave_match, coeff)
    ;new_spect_wave = poly(spect.wavelength,coeff)

    if keyword_set(noplot) then return,outdata

    title='InstrumentModeId='+strtrim(string(instrumentModeId),2)+' at SD='+strtrim(string(gps2sd(spect.(0)[0]/1d6),format='(F0.3)'),2)
    plot_multi,outdata[p].wavelength, outdata[p].wave_match,outdata[p].wavelength, yfit, $
        outdata[p].wavelength, residuals, ytitle='Matched Wavelengths (nm)',xtitle='Reference Wavelength (nm)',$
        psym=[-4,-3,-5],label=[' ','3rd order fit','Residuals'],title=title,/xst,/yst

    ;lineplot,outdata[p].wavelength, outdata[p].wave_match, ytitle='Matched Wavelengths (nm)',xtitle='Reference Wavelength (nm)',psym=-4
    ;lineplot,outdata[p].wavelength, yfit, title='3rd order fit',psym=-3
    ;lineplot,outdata[p].wavelength, residuals, title='Residuals',psym=-5

    ;lineplot,outdata[p].wave_match, outdata[p].wavelength, xtitle='Matched Wavelengths (nm)',ytitle='Reference Wavelength (nm)',psym=-4
    ;lineplot,outdata[p].wave_match, yfit, title='3rd order fit',psym=-3
    ;lineplot,outdata[p].wave_match, residuals, title='Residuals',psym=-5

    return, outdata

end
