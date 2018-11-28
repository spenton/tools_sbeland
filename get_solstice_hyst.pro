;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Fit a Hyteresis function by minimizing the differences between the 
;   UP and DOWN scans.
;
; CALLING SEQUENCE:
;   spectra = GET_SOLSTICE_HYST(t0, t1, instrumentModeId, /missionDays)
;
; INPUT PARAMETERS:
;   startTimes -
;      A 2-element array with the lower time range for which data will 
;      be used.
;      May be specified as julian days, gps microseconds, or
;      SORCE mission day numbers.
;
;   stopTimes -
;      A 2-element array with the upper time range for which data will 
;      be used.
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
;
; RETURNED PARAMETERS:
;   
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
;   Revision: $Id: get_solstice_hyst.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;------------------------------------------------------------------
;
function my_solhyst_amoeba, coeffs
    common mysolhyst_amoeba_common, spup, spdwn, filt1up, filt1dwn, filter_scale

    yup = [spup[0].dn, spup[1:-1].dn + coeffs[0] * spup[0:-2].dn]
    ydwn = [spdwn[0].dn, spdwn[1:-1].dn + coeffs[0] * spdwn[0:-2].dn]

    sup=sort(spup.wavelength)
    sdwn=sort(spdwn.wavelength)
    ydwn = interpol(ydwn[sdwn], spdwn[sdwn].wavelength, spup[sup].wavelength)

    if n_elements(filt1up) gt 1 then yup[filt1up] *= filter_scale
    if n_elements(filt1dwn) gt 1 then ydwn[filt1dwn] *= filter_scale

    ;p=where(yup ne 0d)
    ;delta = (yup[p] - ydwn[p])/yup[p]
    ;corr_value=stdev(gauss_smooth(delta,width=10,/edge_mirror))
    corr_value = sqrt(total(smooth(yup[sup] - ydwn, 20, /edge_mirror)^2d))

    print,coeffs[0],corr_value[0]
    return, corr_value
end

;*******************************************************************

function GET_SOLSTICE_HYST, startTimes, stopTimes, instrumentModeId, $
         gps=gps, missionDays=missionDays, julianDays=julianDays, $
         spectup=myspup, spectdwn=myspdwn, wrange=wrange

    common mysolhyst_amoeba_common

    if keyword_set(noderiv) then no_deriv=1 else no_deriv=0
    if keyword_set(rms) then use_rms=1 else use_rms=0

    if instrumentModeId lt 7 or instrumentModeId gt 13 then begin
        print,'Error: wrong instrumentModeId provided [7,9,11,13]'
        return,-1
    endif

    filter_scale=10d  

    ; get the data
    if n_elements(myspup) eq 0 then $
        myspup =get_solstice_spectra(startTimes[0],stopTimes[0],missionDays=missionDays, gps=gps, julianDays=julianDays, instrumentModeId)
    if n_elements(myspdwn) eq 0 then $
        myspdwn=get_solstice_spectra(startTimes[1],stopTimes[1],missionDays=missionDays, gps=gps, julianDays=julianDays, instrumentModeId)

    ; only consider the overlapping wavelengths
    mnup=min(myspup.wavelength,max=mxup)
    mndwn=min(myspdwn.wavelength,max=mxdwn)
    my_wrange=[max([mnup,mndwn]), min([mxup,mxdwn])]
    if n_elements(wrange) eq 0 then wrange=my_wrange

    posup = where(myspup.wavelength ge wrange[0] and myspup.wavelength le wrange[1],countup)
    if countup eq 0 then begin
        print,'Error: no data in spup covering the requested wavelength range'
        return,-1
    endif
    posdwn = where(myspdwn.wavelength ge wrange[0] and myspdwn.wavelength le wrange[1],countdwn)
    if countdwn eq 0 then begin
        print,'Error: no data in spdwn covering the requested wavelength range'
        return,-1
    endif
    spup=myspup[posup]
    spdwn=myspdwn[posdwn]

    s=sort(spup.timetag)
    spup=spup[s]
    s=sort(spdwn.timetag)
    spdwn=spdwn[s]

    ; get_solstice_spectra applies the filter correction - remove it here
    filt1up=where(spup.filter1_in eq 1,count)
    if count gt 0 then spup[filt1up].dn /= filter_scale
    filt1dwn=where(spdwn.filter1_in eq 1,count)
    if count gt 0 then spdwn[filt1dwn].dn /= filter_scale
 
    coeffs = AMOEBA(1.0e-8, P0=10.0d, scale=0.5, FUNCTION_VALUE=fval, FUNCTION_NAME='my_solhyst_amoeba')

    if keyword_set(verbose) then print,'Final coefficients = ',coeffs
    if keyword_set(verbose) then print,'Fit goodness=',1d - fval / 1d6

    yup = [spup[0].dn, spup[1:-1].dn + coeffs[0] * spup[0:-2].dn]
    ydwn = [spdwn[0].dn, spdwn[1:-1].dn + coeffs[0] * spdwn[0:-2].dn]
    if n_elements(filt1up) gt 1 then begin
        yup[filt1up] *= filter_scale
        spup[filt1up].dn *= filter_scale
    endif
    if n_elements(filt1dwn) gt 1 then begin
        ydwn[filt1dwn] *= filter_scale
        spdwn[filt1dwn].dn *= filter_scale
    endif
    sup=sort(spup.wavelength)
    sdwn=sort(spdwn.wavelength)
    newydwn = interpol(ydwn[sdwn], spdwn[sdwn].wavelength, spup[sup].wavelength)
    newy=interpol(spdwn[sdwn].dn,spdwn[sdwn].wavelength,spup[sup].wavelength)

    delta = (yup[sup] - newydwn)/yup[sup]

    plot_multi,spup.wavelength,smooth((spup[sup].dn-newy)/spup[sup].dn,20),$
        spup.wavelength,smooth(delta,20),/xst,/yst, psym=[-4,-5], $
        xtitle='Wavelength',ytitle='Fractional Diffferences',$
        title='Solstice Up and Down Scans on SD='+strtrim(string(gps2sd(spup[0].(0)/1d6)),2), $
        label=['Original Diff','Minimized Diff Coeff='+strtrim(string(coeffs[0]),2)]

 
    return, coeffs

end



