;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Find the best fit F-Function at the specified wavelength for the time
;   range requested using either get_esrdeg or get_diodedeg. If no Kappa 
;   value is provided, get_esrdeg or get_diodedeg will use the default values.
;
; CALLING SEQUENCE:
;   spectra = get_global_degcol(t0, t1, /missionDay)
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
;   version -
;      Requests a specific version of the SimProfileIntegralCal, otherwise
;      pick the latest version valid for the specified time range.
;   order - 
;      Order of the polynomial correction.
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
;
; REVISION HISTORY:
;   Revision: $Id: get_global_degcol.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;------------------------------------------------------------------
;
function get_global_degcol, startTime, stopTime, wavelength=wavelength, $
         gps=gps, missionDays=missionDays, julianDays=julianDays, $
         version=version, order=order, inspectrum=inspectrum, $
         solar54=solar54, solar55=solar55, step55=step55, lya=lya, smooth=smooth, $
         coeff_corr=coeff_corr, uvdiode=uvdiode, visdiode=visdiode, kappa=kappa, outdata=outdata, multifit=multifit


    if n_elements(wavelength) eq 0 then begin
        if keyword_set(uvdiode) then begin
            ; this array comes from the uv spectra on day 453 keeping only every 8th wavelength (142 positions)
            wavelength = [201.14d, 201.46, 201.78, 202.10, 202.43, 202.76, 203.09, 203.43, 203.76, 204.11, $
                          204.45, 204.79, 205.14, 205.50, 205.85, 206.21, 206.57, 206.94, 207.30, 207.67, $
                          208.05, 208.43, 208.81, 209.19, 209.58, 209.97, 210.37, 210.76, 211.17, 211.57, $
                          211.98, 212.40, 212.81, 213.24, 213.66, 214.09, 214.52, 214.96, 215.40, 215.85, $
                          216.30, 216.76, 217.22, 217.68, 218.15, 218.62, 219.10, 219.58, 220.07, 220.56, $
                          221.06, 221.56, 222.07, 222.58, 223.10, 223.63, 224.16, 224.69, 225.23, 225.78, $
                          226.33, 226.89, 227.46, 228.03, 228.61, 229.19, 229.78, 230.38, 230.98, 231.59, $
                          232.21, 232.84, 233.47, 234.11, 234.76, 235.41, 236.08, 236.75, 237.43, 238.11, $
                          238.81, 239.51, 240.23, 240.95, 241.68, 242.42, 243.17, 243.93, 244.70, 245.48, $
                          246.27, 247.07, 247.88, 248.70, 249.53, 250.37, 251.23, 252.10, 252.98, 253.87, $
                          254.77, 255.69, 256.62, 257.56, 258.52, 259.49, 260.48, 261.48, 262.50, 263.53, $
                          264.58, 265.64, 266.72, 267.82, 268.94, 270.07, 271.22, 272.39, 273.58, 274.79, $
                          276.02, 277.28, 278.55, 279.84, 281.16, 282.50, 283.87, 285.26, 286.67, 288.12, $
                          289.58, 291.08, 292.61, 294.16, 295.75, 297.36, 299.01, 300.69, 302.41, 304.16, $
                          305.95, 307.78]
            ; limit the wavelength range - start at 220nm (everything below is very noisy)
            p=where(wavelength gt 220d)
            wavelength=wavelength[p]
        endif else if keyword_set(visdiode) then begin
            ; this array comes from the visA spectra on day 453 keeping only every 4th wavelength (142 positions)
            wavelength = [310.02, 310.95, 311.89, 312.84, 313.81, 314.78, 315.76, 316.75, 317.76, 318.78, $
                         319.80, 320.84, 321.90, 322.96, 324.04, 325.13, 326.23, 327.34, 328.47, 329.62, $
                         330.77, 331.94, 333.13, 334.33, 335.55, 336.78, 338.03, 339.29, 340.57, 341.86, $
                         343.18, 344.51, 345.86, 347.23, 348.61, 350.02, 351.44, 352.89, 354.35, 355.84, $
                         357.34, 358.87, 360.42, 362.00, 363.59, 365.22, 366.86, 368.53, 370.23, 371.95, $
                         373.70, 375.48, 377.28, 379.12, 380.98, 382.88, 384.81, 386.76, 388.76, 390.78, $
                         392.84, 394.94, 397.07, 399.25, 401.46, 403.71, 406.00, 408.33, 410.71, 413.14, $
                         415.61, 418.12, 420.69, 423.31, 425.98, 428.70, 431.48, 434.32, 437.22, 440.18, $
                         443.20, 446.29, 449.45, 452.67, 455.97, 459.35, 462.80, 466.33, 469.95, 473.65, $
                         477.45, 481.34, 485.32, 489.41, 493.60, 497.90, 502.32, 506.85, 511.50, 516.28, $
                         521.20, 526.25, 531.45, 536.80, 542.31, 547.98, 553.82, 559.84, 566.05, 572.45, $
                         579.06, 585.87, 592.92, 600.19, 607.71, 615.48, 623.52, 631.83, 640.43, 649.34, $
                         658.56, 668.11, 678.00, 688.25, 698.87, 709.88, 721.29, 733.11, 745.36, 758.05, $
                         771.19, 784.79, 798.87, 813.43, 828.47, 844.01, 860.03, 876.55, 893.55, 911.03, $
                         928.98, 947.39]
            ; limit the wavelengths to range where AMEOBA converges 
            p=where(wavelength ge 325.0) ; and wavelength le 350.0)
            wavelength=wavelength[p]
        endif else begin
            ;wavelength=[287.0d,   292.5d,   303.3d,   317.7d,   332.0d,   340.0d,   $
            ;            353.0d,   367.0d,   373.0d,   382.0d,   404.0d,   428.0d,   $
            ;            467.0d,   476.0d,   484.0d,   515.0d,   561.0d,   587.0d] ;,   650.0d]
            wavelength=[317.7d,   332.0d,   340.0d,   $
                        353.0d,   367.0d,   373.0d,   382.0d,   404.0d,   428.0d,   $
                        467.0d,   476.0d,   484.0d,   515.0d,   561.0d,   587.0d,   $
                        650.0d,   699.5d,   746.0d,   792.0d,   844.5d]
        endelse
    endif

    if keyword_set(missionDays) then begin
       ; user specified time in mission (sorce) days
       t0 = startTime
       t1 = stopTime
    endif else if keyword_set(julianDays) then begin
       ; user specified time in julian days
       t0 = jd2sd(startTime)
       t1 = jd2sd(stopTime)
    endif else begin
       ; user specified timetags in gps microseconds
       t0 = gps2sd(startTime/1.d6)
       t1 = gps2sd(stopTime/1.d6)
    endelse

    if n_elements(step55) eq 0 then step55=90d

    ; get the SimSolarExposureData for modes 31 and 32
    query_database,/reset
    if n_elements(solar54) eq 0 then begin
        q1 = 'SELECT orbitStartTimeGps t0, orbitEndTimeGps t1, cumulativeSolarExpHrtOut solar_exp, '+$
             'solarExposureHrtOut as solar_exp_orbit from SimSolarExposureData where instrumentModeId=54'
        query_database, q1, solar54, info
    endif
    if n_elements(solar55) eq 0 then begin
        q1 = 'SELECT orbitStartTimeGps t0, orbitEndTimeGps t1, cumulativeSolarExpHrtOut solar_exp, '+$
             'solarExposureHrtOut as solar_exp_orbit from SimSolarExposureData where instrumentModeId=55'
        query_database, q1, solar55, info
    endif

    if max(strpos(tag_names(solar54),'ONEAU')) lt 0 then begin
        print,'  getting SUNOBSERVERDISTANCECORRECTION from database ...'
        q1="SELECT microsecondsSinceGpsEpoch, sunObserverDistanceCorrection FROM SolarDistAndDopplerFixedStep2 "+$
           "WHERE instrumentModeId=50 and version=7"
        query_database,/reset
        query_database, q1, solardist, info
        oneau54 = interpol(solardist.SUNOBSERVERDISTANCECORRECTION, solardist.MICROSECONDSSINCEGPSEPOCH, solar54.t1)
        append_tag,solar54,'oneau',oneau54,/slim
        oneau55 = interpol(solardist.SUNOBSERVERDISTANCECORRECTION, solardist.MICROSECONDSSINCEGPSEPOCH, solar55.t1)
        append_tag,solar55,'oneau',oneau55,/slim
    endif

    if keyword_set(lya) and max(strpos(tag_names(solar54),'LYALPHA_CORR')) lt 0 then begin
        print,'  getting Lyman Alpha dosage from Level3 data ...'
        ; extract the Lyman alpha irradiance from SOLSTICE data normalizing around the solar minimum
        res=get_level3_spectrum(0d,5000d, 11, 24d, /mission, min_wave=121d, max_wave=122d, /released)
        ; normalized to irradiance from "quiet sun"
        quiet_time=sd2gps([2130d,2230d])*1d6
        p=where(res.(1) ge quiet_time[0] and res.(1) le quiet_time[1])
        avg_irrad = mean(res[p].irradiance)
        res.irradiance /= avg_irrad
        ; add the oneau to Lyman Alpha to have a more accurate "dosage" factor (in reality makes little difference)
        lyalpha_corr = interpol(res.irradiance, res.(1), solar54.t1)
        min_time=min(res.(1),pos)
        p=where(solar54.t1 lt min_time,count)
        if count gt 0 then lyalpha_corr[p]=res[pos].irradiance
        lyalpha_corr *= solar54.oneau
        append_tag,solar54,'LYALPHA_CORR',lyalpha_corr,/slim
        lyalpha_corr = interpol(res.irradiance, res.(1), solar55.t1)
        p=where(solar55.t1 lt min_time,count)
        if count gt 0 then lyalpha_corr[p]=res[pos].irradiance
        lyalpha_corr *= solar55.oneau
        append_tag,solar55,'LYALPHA_CORR',lyalpha_corr,/slim
    endif

    ; with Ly alpha correction, the fit is linear (2nd order without)
    ;coeff_corr=[0.92418784d, -7.5035887d-05]
    ;coeff_corr=[1.4907146d, -0.00069623171d, 1.6402941d-07]

    if n_elements(coeff_corr) gt 0 then begin
        corr54 = poly(gps2sd(solar54.(1)/1d6),coeff_corr)
        corr55 = poly(gps2sd(solar55.(1)/1d6),coeff_corr)
    endif else begin
        corr54=1d
        corr55=1d
    endelse

    temp54=solar54
    temp55=solar55
    if keyword_set(lya) then begin
        temp54.solar_exp=total(temp54.solar_exp_orbit*temp54.LYALPHA_CORR*corr54,/cum) / (1d +(temp54.oneau-1d)/4d)
        temp55.solar_exp=total(temp55.solar_exp_orbit*temp55.LYALPHA_CORR*corr55,/cum) / (1d +(temp55.oneau-1d)/4d)
        def_order=1
    endif else begin
        temp54.solar_exp=total(temp54.solar_exp_orbit*corr54,/cum) / (1d +(temp54.oneau-1d)/4d)
        temp55.solar_exp=total(temp55.solar_exp_orbit*corr55,/cum) / (1d +(temp55.oneau-1d)/4d)
        ;temp54.solar_exp=total(temp54.solar_exp_orbit*corr54,/cum)
        ;temp55.solar_exp=total(temp55.solar_exp_orbit*corr55,/cum)
        def_order=4
    endelse

    if n_elements(order) eq 0 then order=def_order

    ; define a user defined symbol with psym=8
    A = FINDGEN(17) * (!PI*2/16.)
    USERSYM, COS(A), SIN(A), /FILL

    nwaves=n_elements(wavelength)
    time_f=[]
    f_factor=[]
    f_factor_scaled=[]
    npts=dblarr(nwaves+1)
    for i=0,nwaves-1 do begin
        print,' processing ',wavelength[i]
        if keyword_set(uvdiode) then begin
            ; the kappa in this case is expected to be n array containing x and y
            result=get_uvdiode_fdeg(t0, t1, /mission,  wavelength[i], solar54=temp54, solar55=temp55, $
                inspectrum=inspectrum, step55=step55, smooth=smooth, /noplot, version=version, kappa=kappa)
        endif else if keyword_set(visdiode) then begin
            ; the kappa in this case is expected to be n array containing x and y
            result=get_visdiode_fdeg(t0, t1, /mission,  wavelength[i], solar54=temp54, solar55=temp55, $
                inspectrum=inspectrum, step55=step55, smooth=smooth,  kappa=kappa)
        endif else begin
            ; the kappa in this case is expected to be the coefficients of exponent fit
            result=get_esrdeg(t0, t1, /mission,  wavelength[i], solar54=temp54, solar55=temp55, $
                inspectrum=inspectrum, /alignobc, smooth=smooth, step55=step55, /noplot, version=version,kappa_coeff=kappa)
        endelse
        if size(result,/tname) ne 'STRUCT' then continue
        time_f=[time_f,result.time_f]
        f_factor=[f_factor,result.f_factor]
        f_factor_scaled=[f_factor_scaled,result.f_factor_scaled]
        npts[i+1]=n_elements(result.time_f)
    endfor

    if n_elements(time_f) eq 0 then begin
        ; no data was found for all of these wavelengths
        print,'Error: no degradation data was returned for all requested wavelengths'
        return,-1
    endif


    x=time_f+step55/2d
    if keyword_set(uvdiode) or keyword_set(visdiode) then begin
        step=10
        npts=total(npts,/cum)
        for i=0,(nwaves/step-1)>0 do begin
            p0=i*step
            if step gt 1.0 then begin
                p1=(p0+step -1)<(n_elements(wavelength)-1)
                title='F_Factor for '+strtrim(string(wavelength[p0],format='(F0.2)'),2)+', '+$
                    strtrim(string(wavelength[p1],format='(F0.2)'),2)+' nm'
            endif else begin
                p1=(p0+1)<(n_elements(wavelength)-1)
                title='F_Factor for '+strtrim(string(wavelength[p0],format='(F0.2)'),2)
            endelse
            mn=min(wavelength[p0:p1],max=mx)
            title='F_Factor for '+strtrim(string(mn,format='(F0.2)'),2)+' to '+strtrim(string(mx,format='(F0.2)'),2)+' nm'
            xx=x[npts[p0]:npts[p1]-1]
            ;yy=f_factor_scaled[npts[p0]:npts[p1]-1]
            yy=f_factor[npts[p0]:npts[p1]-1]
            k=where(finite(yy),count)
            if count lt 10 then continue
            s=sort(xx[k])
            lineplot,xx[k[s]],yy[k[s]],psym=floor(i/4)+4,charsize=1.5,xtitle='Mission Day',$
                ytitle='F_Factor', title=title, ptitle='Combined F_Factor for Diodes'
            if keyword_set(multifit) then begin
                ; fit this group of data
                p=where(xx[k[s]]-step55/2d lt 4300d)
                coeff=robust_poly_fit(xx[k[s[p]]],yy[k[s[p]]],order,/double)
                lineplot,dindgen(4300),poly(dindgen(4300),coeff), psym=-3, thick=2.0, title=title
            endif
        endfor
    endif else begin
        ; group the data in wavelengths
        step=1.0
        npts=total(npts,/cum)
        for i=0,(nwaves/step-1)>0 do begin
            p0=i*step
            if step gt 1.0 then begin
                p1=p0+step -1 
                title='F_Factor for '+strtrim(string(wavelength[p0],format='(F0.2)'),2)+', '+$
                    strtrim(string(wavelength[p1],format='(F0.2)'),2)+' nm'
            endif else begin
                p1=p0+1
                title='F_Factor for '+strtrim(string(wavelength[p0],format='(F0.2)'),2)
            endelse
            ; make sure the data is plotted with increasing wavelength
            xx=x[npts[p0]:npts[p1]-1]
            yy=f_factor_scaled[npts[p0]:npts[p1]-1]
            s=sort(xx)
            lineplot,xx[s],yy[s],psym=floor(i/4)+4,charsize=1.5,xrange=[0d,4000d],yrange=[0d,2d], xtitle='Mission Day',$
                ytitle='F/Kappa', title=title, ptitle='Combined F_Factor for wavelengths '+$
                strtrim(string(wavelength[0],format='(I0)'),2)+'->'+strtrim(string(wavelength[-1],format='(I0)'),2),font=-1
            if keyword_set(multifit) then begin
                ; fit this group of data
                p=where(xx[s]-step55/2d lt 2805d)
                coeff=robust_poly_fit(xx[s[p]],yy[s[p]],order,/double)
                lineplot,dindgen(4000),poly(dindgen(4000),coeff), psym=-3, thick=2.0, title=title, font=-1
            endif
        endfor
    endelse
    k=where(finite(f_factor_scaled),count)
    if count ne n_elements(f_factor_scaled) then begin
        time_f=time_f[k]
        f_factor=f_factor[k]
        f_factor_scaled=f_factor_scaled[k]
    endif
    ;s=sort(time_f)
    ;time_f=time_f[s]
    ;f_factor=f_factor[s]
    ;f_factor_scaled=f_factor_scaled[s]
    ; only fit the data prior to the OBC on day 2805 
    ;p=where(time_f lt 2805d)
    ;p=where(time_f lt 2600d)
    ;coeff=robust_poly_fit(time_f[p]+step55/2d,f_factor_scaled[p],order,/double)
    p=lindgen(n_elements(time_f))
    coeff=robust_poly_fit(time_f[p]+step55/2d,f_factor[p],order,/double)
    title='fit up to day '+strtrim(string(max(time_f[p]),format='(F6.0)'),2)+', coeff= '
    for i=0,n_elements(coeff)-1 do title=title+strcompress(string(coeff[i]),/rem)+', '
    newx=dindgen(4300/10)*10d
    lineplot,newx,poly(newx,coeff), psym=-3, thick=2.0, title=title, xrange=[0d,4000d],yrange=[0d,2d],font=-1
    k=where(newx ge min(time_f[p]+step55/2d,max=mx) and newx le mx)
    myy=interpol(f_factor,time_f+step55/2d,newx[k])
    sset=bspline_iterfit(newx[k],myy,maxiter=0,requiren=6,bkspace=5)
    newy=bspline_valu(newx[k],sset)
    lineplot,newx[k],newy,psym=-3,thick=4.0, title='BSpline Fit'

    outdata={time_f:time_f, f_factor:f_factor, f_factor_scaled:f_factor_scaled, npts:npts, wavelength:wavelength}

    return, coeff

end
