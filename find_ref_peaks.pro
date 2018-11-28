;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Find the peaks from the reference spectra from the derivative zero
;   crossings to provide a list of peaks and valeys used to adjust the
;   dispersion solution for other spectra.
;
; CALLING SEQUENCE:
;   list = FIND_REF_PEAKS(refT0=t0, refT1=t1, instrumentModeId, /missionDays)
;
; INPUT PARAMETERS:
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
;   Revision: $Id: find_ref_peaks.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;-
;
function find_ref_peaks, instrumentModeId, x=xval, y=yval, refSpect=refSpect, $
         refT0=refT0, refT1=refT1, $
         gps=gps, missionDays=missionDays, julianDays=julianDays, $
         version=version, quiet=quiet, d1height=d1height, d2height=d2height

    if n_elements(xval) eq 0 or n_elements(yval) eq 0 then begin
        if size(refSpect,/tname) ne 'STRUCT' then begin
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

           refSpect = get_sim_spectra(refT0, refT1, instrumentModeId, gps=gps, $
                missionDays=missionDays, julianDays=julianDays, version=version, /noCcdCorr, /dn_only, /noTempCorr)
            if size(refSpect,/tname) ne 'STRUCT' then begin
                print,'No valid reference spectra found'
                return,-1
            endif
        endif
        xval = refSpect.wavelength
        yval = refSpect.dn
    endif

    ; get the derivative of the normalized spectra
    mn=min(yval,max=mx)
    norm_dn = (yval-mn)/(mx-mn)
    derivYval = deriv(xval, norm_dn)
    npos=n_elements(xval)

    if n_elements(d1height) eq 0 and n_elements(d2height) eq 0 then begin
        d1height = 2.0d-4
        d2height=d1height
    endif else if n_elements(d1height) eq 0 then begin
        d1height=d2height
    endif else if n_elements(d2height) eq 0 then begin
        d2height=d1height
    endif

    ; look for zero crossings
    zero_pos = where((derivYval[0:-2] gt 0.0 and derivYval[1:-1] le 0.0) OR $
                     (derivYval[0:-2] le 0.0 and derivYval[1:-1] gt 0.0), count)
    ; get rid of consecutive indexes (too sharp of features - probably just noise)
    if count lt 2 then begin
        print,'No zero-crossing found in spectrum'
        return,-1
    endif

    p=where(zero_pos[1:-1]-zero_pos[0:-2] le 1,cc)
    if cc gt 0 then begin
        ; remove both 
        p=[p,p+1]
        zero_pos[p]=-1
        p=where(zero_pos ge 0,count)
        zero_pos=zero_pos[p]
    endif
    ; start from the beginning
    zero_pos=[zero_pos,npos-1L]
    p=where(zero_pos[1:-1]-zero_pos[0:-2] gt 0,count)
    if count gt 0 then zero_pos=zero_pos[p]
    peakwave=[]
    peakdn=[]
    peakderiv=[]
    for i=0L,count-1L do begin
       ; validate each zero-crossing
       if i eq 0 then p0=0 else p0=zero_pos[i-1]
       if i eq count-1L then p1=npos-1L else p1=zero_pos[i+1]
       if max(ABS(derivYval[p0:zero_pos[i]])) GT d1height AND $
          max(ABS(derivYval[zero_pos[i]:p1])) GT d1height then begin
          ; found a valid peak - find the wavelength at the zero-crossing
          p0=zero_pos[i]-2 > 0
          p1=zero_pos[i]+2 < (n_elements(xval)-1L)
          zero_wave=interpol(xval[p0:p1], derivYval[p0:p1], 0.0)
          peakwave = [peakwave,zero_wave]
          zero_pos_dn=interpol(yval[p0:p1], xval[p0:p1], zero_wave)
          peakdn = [peakdn, zero_pos_dn]
          peakderiv = [peakderiv,1]
       endif
    endfor

    ; look for second derivative zero-crossing (inflection points)
    norm_deriv = (derivYval-min(derivYval))/(max(derivYval) - min(derivYval))
    dderivYval = deriv(xval, norm_deriv)
    zero_pos = where((dderivYval[0:-2] gt 0.0 and dderivYval[1:-1] le 0.0) OR $
                     (dderivYval[0:-2] le 0.0 and dderivYval[1:-1] gt 0.0), count)
    
    if count gt 1 then begin
        ; get rid of consecutive indexes (too sharp of features - probably just noise)
        p=where(zero_pos[1:-1]-zero_pos[0:-2] eq 1,cc)
        if cc gt 0 then begin
            ; remove both 
            p=[p,p+1]
            zero_pos[p]=-1
            p=where(zero_pos ge 0,count)
            zero_pos=zero_pos[p]
        endif
        zero_pos=[zero_pos,npos-1L]
        p=where(zero_pos[1:-1]-zero_pos[0:-2] gt 0,count)
        if count gt 0 then zero_pos=zero_pos[p]
        for i=0L,count-1L do begin
           if i eq 0 then p0=0 else p0=zero_pos[i-1]
           if i eq count-1L then p1=npos-1L else p1=zero_pos[i+1]
           ; validate each zero-crossing
           if max(ABS(dderivYval[p0:zero_pos[i]])) GT d2height AND $
              max(ABS(dderivYval[zero_pos[i]:p1])) GT d2height then begin
              ; found a valid peak - find the wavelength at the zero-crossing
              p0=zero_pos[i]-2 > 0
              p1=zero_pos[i]+2 < (n_elements(xval)-1L)
              zero_wave=interpol(xval[p0:p1], dderivYval[p0:p1], 0.0)
              peakwave = [peakwave,zero_wave]
              zero_pos_dn=interpol(yval[p0:p1], xval[p0:p1], zero_wave)
              peakdn = [peakdn, zero_pos_dn]
              peakderiv = [peakderiv,2]
           endif
        endfor
    endif

    if n_elements(peakdn) eq 0 then begin
        print,'No peaks found in spectrum'
        return,-1
    endif

    outdata=replicate({x:0.0d, y:0.0d, deriv_order:0}, n_elements(peakdn))
    s=sort(peakwave)
    outdata.x = peakwave[s]
    outdata.y = peakdn[s]
    outdata.deriv_order = peakderiv[s]
    return, outdata

end
