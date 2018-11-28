;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Using a known time range, determine the change in the prism index 
;   of refraction as a function of temperature to map against what is
;   used in George Laurence's equations (LAM_N.pro, REF_INDEX.pro)
;
; CALLING SEQUENCE:
;
; INPUT PARAMETERS:
;
; OPTIONAL INPUT PARAMETERS:
;   version -
;      Requests a specific version of the SimProfileIntegralCal, otherwise
;      pick the latest version valid for the specified time range.
;
; RETURNED PARAMETERS:
;
; OPTIONAL OUTPUT PARAMETERS:
;
; EXAMPLE:  
;
; COMMON BLOCKS:
;      NONE
;
; NOTES:
;
; REVISION HISTORY:
;   Revision: $Id: get_dindex_dt.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;------------------------------------------------------------------
function get_dindex_dt, version=version, outdata=newSpect

    if n_elements(version) eq 0 then version=20

    t0 = 3640.0d
    t1 = 3684.0d
    refT0 = 3654.2623d
    refT1 = 3654.2754d

    ; change in index of refraction with temperature modeled by a 5th order polynomial
    c_tcn = [10.91141546d, 4.130107857d, -7.400374025d, 4.110286051d, -0.937939315d, 0.078332038d]

    TREF_SILICA = 20.0D0   ; ZEMAX AND MALITSON SELLMEIER TEMPERATURE
    ; SELLMEIER COEFS
    KS1 = 0.6961663D
    LS1 = 0.004679148D
    KS2 = 0.4079426D  
    LS2 = 0.01351206D
    KS3 = 0.8974794D
    LS3 = 97.934D

    ; get the refSpect aligned with the day 453 spectra
    refSpect = align_spectra(refT0, refT1, 43, /mission, refspect=spect453, version=version, profile_data=profile, order=2)
    nwaves = n_elements(refSpect)
    ref_index=dblarr(nwaves)
    refTemp = 17.22d  ; refspect was aligned with this reference temperature

    if size(newSpect,/tname) ne 'POINTER' then begin
        pla = get_sorce_plan(t0, t1, /mission, /sima, activity='SolarQuickScan24')

        ; now align every SolarQuickScan24 spectrum wihthin the time period
        nscans = n_elements(pla)
        newSpect = ptrarr(nscans)
        for i=0,nscans-1 do begin 
            print,i 
            tmp=align_spectra(pla[i].starttime,pla[i].stoptime, 43, /mission, refspect=spect453, version=version, profile_data=profile, order=2) 
            if size(tmp,/tname) ne 'STRUCT' then continue
            newSpect[i]=ptr_new(tmp) 
        endfor
    endif else begin
        nscans = n_elements(newSpect)
    endelse

    ; We now have a known wavelength at every ccd position and we assume the only wavelength shift over this time period
    ; is due to prism temperature only
    ; Calculate the index of refraction as a function of temperature
    index = dblarr(nscans, nwaves)
    prism_temp = index
    delta_index=index
    waves =index
    delta_temp=index
    timetag=index

    for j=0,nwaves-1 do begin
        LSQ = (refSpect[j].wavelength/1d3) ^ 2d
        ref_index[j] = SQRT(KS1*LSQ/(LSQ-LS1) + KS2*LSQ/(LSQ-LS2) + KS3*LSQ/(LSQ-LS3) + 1d)

    endfor

    for i=0,nscans-1L do begin
        if not ptr_valid(newSpect[i]) then continue
        for j=0,nwaves-1 do begin
            p = where((*newSpect[i]).ccdpos_uncorr eq refSpect[j].ccdpos_uncorr, count)
            if count eq 0 then continue
            LSQ = ((*newSpect[i])[p].wavelength/1d3) ^ 2d
            index[i,j] = SQRT(KS1*LSQ/(LSQ-LS1) + KS2*LSQ/(LSQ-LS2) + KS3*LSQ/(LSQ-LS3) + 1d)
            delta_index[i,j] = index[i,j] - ref_index[j]
            prism_temp[i,j] = (*newSpect[i])[p].prismtemp
            delta_temp[i,j] = prism_temp[i,j] - refTemp
            waves[i,j] =(*newSpect[i])[p].wavelength
            timetag[i,j] = (*newSpect[i])[p].timetag
        endfor
    endfor

    timetag = gps2sd(timetag/1d6)
    step=20d
    wstep = ceil(nwaves/step)
    for i=0, step-1L do begin
        ;lineplot,prism_temp[*,i*step], index[*,i*step], xtitle='Prism Temperature', ytitle='Prism Index of Refraction',psym=4, $
        ;    title='Wavelength='+strtrim(string(median(waves[*,i*step]),format='(F0.2)'))
        k=where(waves[*,i*wstep] gt 0d)
        lineplot,1d/(waves[k,i*wstep]/1d3), delta_index[k,i*wstep] / delta_temp[k,i*wstep], xtitle='1 / wavelength (um-1)', ytitle='Delta_Index',psym=4, $
            title='Wavelength='+strtrim(string(median(waves[k,i*wstep]),format='(F0.2)'))
    endfor

    return,{index:index, prism_temp:prism_temp, wavelength:waves, delta_index:delta_index, timetag:timetag} 

end
