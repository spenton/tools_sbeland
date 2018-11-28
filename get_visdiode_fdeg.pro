;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Extract the degradation information from the UV diode matching
;   the data for the requested wavelength taken at the "same" time for 
;   SIMA and SIMB.  This is done for two consecutive SolarQuickScan24 of 
;   SIMBESRB to get a running slope of the degradation with mission time.
;
; CALLING SEQUENCE:
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
;    The table used for the irradiance is expected to have no degradation
;    correction applied to it (as in version 15 of our development database for version 19).
;
;
; REVISION HISTORY:
;   Revision: $Id: get_visdiode_fdeg.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;------------------------------------------------------------------
function my_amoeba_fdeg, coeffs
    common myamoeba_fdeg_common, irradRatio, solExpA0, solExpA1, solExpB0, solExpB1, rayPath, kappadeg

    ;coeffs=abs(coeffs)
    temp54=total(solExpA1.solar_exp_orbit * coeffs[0],/cum) / (1d +(solExpA1.oneau-1d)/4d) / 86400d
    temp55=total(solExpB1.solar_exp_orbit * coeffs[0],/cum) / (1d +(solExpB1.oneau-1d)/4d) / 86400d
    solA = solExpA0 + temp54[-1]
    solB = solExpB0 + temp55[-1]
    valueA0 = ( (1d - rayPath)*exp(-kappadeg*solExpA0) + rayPath*exp(-kappadeg*solExpA0/2d) )
    valueA1 = ( (1d - rayPath)*exp(-kappadeg*solA) + rayPath*exp(-kappadeg*solA/2d) )
    valueB0 = ( (1d - rayPath)*exp(-kappadeg*solExpB0) + rayPath*exp(-kappadeg*solExpB0/2d) )
    valueB1 = ( (1d - rayPath)*exp(-kappadeg*solB) + rayPath*exp(-kappadeg*solB/2d) )

    residual = ABS(irradRatio - (valueB1 / ValueB0) / (valueA1 / ValueA0))

    ;print,coeffs, residual,format='(2(G,","),G)'
    return, residual

end
;*******************************************************************

function get_visdiode_fdeg, starttime, stoptime, wavelength, bin55=bin55, step55=step55, $
    gps=gps, missionDays=missionDays, julianDays=julianDays, $
    inspectrum=inspectrum, alignobc=alignobc, smooth=smooth, $
    noplot=noplot, solar54=solar54, solar55=solar55, kappa=kappa

    common myamoeba_fdeg_common

    ; for some obscure reasons, IDL doesn't find the function MEAN used below and
    ; stops execution on that line (even though this is an IDL supplied routine ???)
    FORWARD_FUNCTION mean

    if keyword_set(gps) then begin
       ; user specified time in mission (sorce) days
       t0 = gps2sd(startTime/1.d6)
       t1 = gps2sd(stopTime/1.d6)
    endif else if keyword_set(julianDays) then begin
       ; user specified time in julian days
       t0 = jd2sd(startTime)
       t1 = jd2sd(stopTime)
    endif else begin
       ; user specified timetags in gps microseconds
       t0 = startTime
       t1 = stopTime
    endelse

    ;convert times to microsecondssince gpsecpoch
    t0=sd2gps(t0)*1d6
    t1=sd2gps(t1)*1d6

    ; get the vis diode raypath from the raytrace data
    readcol,'~/SORCE/data/overlap_vis1.txt',pd_aw,pd_apath,format='(d,d)',/silent
    ;just as a test multiply the raytrace raypath by 4.8
    pd_apath *= 4.8d

    if size(inspectrum,/tname) ne 'STRUCT' then begin
        ; use the UNCORRECTED IRRADIANCE from version 2000 (which has the 1AU factor applied)
        ; restore SimA and SimB data
        restore,'~/SORCE/data/sima_vis_uncorr_2011.sav'
        pa=where(visa_uncorr_2011.spect20[*].timestamp[0] ge t0 and visa_uncorr_2011.spect20[*].timestamp[0] le t1,count)
        if count eq 0 then begin
            print,'Error: no VISA data within requested time range'
            return,-1
        endif

        restore,'~/SORCE/data/simb_vis_uncorr_2011.sav'
        pb=where(visb_uncorr_2011.spect20[*].timestamp[0] ge t0 and visb_uncorr_2011.spect20[*].timestamp[0] le t1,count)
        if count eq 0 then begin
            print,'Error: no VISB data within requested time range'
            return,-1
        endif
        inspectrum = {simA:visa_uncorr_2011.spect20[pa], simB:visb_uncorr_2011.spect20[pb]}
    endif

    if size(kappa,/tname) ne 'STRUCT' then begin
        ; for the VIS diode, we use the same Kappa as the ESR
        ;readcol,'~/SORCE/data/ESR_kappa_453_1570_1004.txt',pdkappa_w, pdkappa_k,format='(d,d)'
        ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2001.txt',pdkappa_w, pdkappa_k,format='(d,d)'
        readcol,'~/SORCE/data/ESR_kappa_453_1570_2002.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
        kappa={x:pdkappa_W, y:pdkappa_k}
    endif

    if n_elements(bin55) eq 0 then bin55=1

    ; get the SimSolarExposureData 
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

    ; get the corresponding irradiance at the specific wavelength for each spectrum
    nspect=n_elements(inspectrum.simB)
    solarexpA=dblarr(nspect)
    solarexpB=dblarr(nspect)
    irradA=dblarr(nspect)
    irradB=dblarr(nspect)
    timeA=dblarr(nspect)
    timeB=dblarr(nspect)
    for j=0L,nspect-1L do begin
        s=sort(inspectrum.simB[j].wavelength)
        irradB[j]=interpol(inspectrum.simB[j].irradiance[s],inspectrum.simB[j].wavelength[s],wavelength,/spline)
        val=min(abs(inspectrum.simB[j].wavelength- wavelength),pos)
        timeB[j]=(inspectrum.simB[j].(0))[pos]
        q=where(solar55.t1 le timeB[j],count)
        solarexpB[j]=solar55[q[-1]].solar_exp/86400d

        ; get the closest data to SimB for SimA 
        val=min(abs(timeB[j] - inspectrum.simA[*].timestamp[0]),k)
        s=sort(inspectrum.simA[k].wavelength)
        irradA[j]=interpol(inspectrum.simA[k].irradiance[s],inspectrum.simA[k].wavelength[s],wavelength,/spline)
        val=min(abs(inspectrum.simA[k].wavelength- wavelength),pos)
        timeA[j]=(inspectrum.simA[k].(0))[pos]
        q=where(solar54.t1 le timeA[j],count)
        solarexpA[j]=solar54[q[-1]].solar_exp/86400d
    endfor

    ; if align was provided, re-align the data at the OBC events
    if keyword_set(alignobc) then begin
        ; align SimA and SimB together (making sure they both get the same adjustment)
        posA=replicate('A',n_elements(timeA))
        posB=replicate('B',n_elements(timeB))
        posAB=[posA,posB]
        timeAB=[timeA,timeB]
        irradAB=[irradA,irradB]
        solarexpAB=[solarexpA,solarexpB]
        s=sort(timeAB)
        posAB=posAB[s]
        timeAB=timeAB[s]
        irradAB=irradAB[s]
        solarexpAB=solarexpAB[s]
        align_irrad,timeAB,irradAB,/gps
        posA=where(posAB eq 'A')
        posB=where(posAB eq 'B')
        irradA=irradAB[posA]
        timeA=timeAB[posA]
        solarexpA=solarexpAB[posA]
        irradB=irradAB[posB]
        timeB=timeAB[posB]
        solarexpB=solarexpAB[posB]
    endif

    ; if smooth was flagged, apply a bspline smoothing (from SDSS idlutils/bspline_iterfit.pro)
    timeA=gps2sd(timeA/1d6)
    timeB=gps2sd(timeB/1d6)
    if keyword_set(smooth) then begin
        ; fit a 2nd order curve to data and remove the outliers
        ; SimA
        coeffA=robust_poly_fit(timeA,irradA,2,yfit,/double)
        resistant_mean,(irradA-yfit),3.0,mean,goodvec=keep0
        sset=bspline_iterfit(timeA[keep0],irradA[keep0],maxiter=0,requiren=10,bkspace=5)
        yfitA=bspline_valu(timeA,sset)
        ; SimB
        coeffB=robust_poly_fit(timeB,irradB,2,yfit,/double)
        resistant_mean,(irradB-yfit),3.0,mean,goodvec=keep0
        sset=bspline_iterfit(timeB[keep0],irradB[keep0],maxiter=0,requiren=10,bkspace=5)
        yfitB=bspline_valu(timeB,sset)
    endif else begin
        yfitA=-1
        yfitB=-1
    endelse


    nfact=nspect-bin55
    f_factor=dblarr(nfact) -1d
    time_f=dblarr(nfact) -1d
    rayPath  = (interpol(pd_apath, pd_aw, wavelength, /spline))[0]
    kappadeg = (interpol(kappa.y, kappa.x, wavelength, /spline))[0]

    for j=0L,nfact-1L do begin
        if n_elements(step55) eq 0 then begin
            pos0=lindgen(bin55)+j
            pos1=pos0+1
        endif else begin
            pos0=lindgen(bin55)+j
            p=where(timeA ge (timeA[pos0[0]])[0]+step55,count)
            if count eq 0 then continue
            pos1=lindgen(bin55)+p[0]
        endelse
        if keyword_set(smooth) then begin
            irradRatio = mean(yfitB[pos1]) / mean(yfitB[pos0]) 
            irradRatio /= (mean(yfitA[pos1]) /mean(yfitA[pos0]))
        endif else begin
            irradRatio = mean(irradB[pos1]) / mean(irradB[pos0]) 
            irradRatio /= (mean(irradA[pos1])/mean(irradA[pos0]))
        endelse
        ; we recalculate the solar exposure based on the calculated f_factor in the previous iteration of the for loop
        if j eq 0 then begin
            ; calculate the cummulative solar exposure before the start
            p=where(solar54.t1 lt (sd2gps(timeA[pos0])*1d6)[0],count)
            if count gt 0 then begin
                temp54=total(solar54[p].solar_exp_orbit,/cum) / (1d +(solar54[p].oneau-1d)/4d) / 86400d
                solarexpA[0] = temp54[-1]
            endif else solarexpA[0]=0d
            pA=where(solar54.t1 ge (sd2gps(timeA[pos0])*1d6)[0] and solar54.t1 lt (sd2gps(timeA[pos1])*1d6)[0],count)
            if count gt 0 then begin
                temp54=total(solar54[pA].solar_exp_orbit,/cum) / (1d +(solar54[pA].oneau-1d)/4d) / 86400d
                solarexpA[1] = solarexpA[0] + temp54[-1]
            endif else solarexpA[1]=solarexpA[0]

            p=where(solar55.t1 lt (sd2gps(timeB[pos0])*1d6)[0],count)
            if count gt 0 then begin
                temp55=total(solar55[p].solar_exp_orbit,/cum) / (1d +(solar55[p].oneau-1d)/4d) / 86400d
                solarexpB[0] = temp55[-1]
            endif else solarexpB[0]=0d
            pB=where(solar55.t1 ge (sd2gps(timeB[pos0])*1d6)[0] and solar55.t1 lt (sd2gps(timeB[pos1])*1d6)[0],count)
            if count gt 0 then begin
                temp55=total(solar55[pB].solar_exp_orbit,/cum) / (1d +(solar55[pB].oneau-1d)/4d) / 86400d
                solarexpB[1] = solarexpB[0] + temp55[-1]
            endif else solarexpB[1]=solarexpB[0]
        endif else begin
            pA=where(solar54.t1 ge (sd2gps(timeA[pos0])*1d6)[0] and solar54.t1 lt (sd2gps(timeA[pos1])*1d6)[0],count)
            pB=where(solar55.t1 ge (sd2gps(timeB[pos0])*1d6)[0] and solar55.t1 lt (sd2gps(timeB[pos1])*1d6)[0],count)
        endelse

        solExpA0 = solarexpA[j]
        solExpA1 = solar54[pA]
        solExpB0 = solarexpB[j]
        solExpB1 = solar55[pB]

        ;solExpA0=mean(solarexpA[pos0])
        ;solExpA1=mean(solarexpA[pos1])
        ;solExpB0=mean(solarexpB[pos0])
        ;solExpB1=mean(solarexpB[pos1])
        ; the solarexposures 
        coeffs = AMOEBA(1d-6, P0=1d, scale=1d, FUNCTION_VALUE=fval, FUNCTION_NAME='my_amoeba_fdeg', nmax=maxiter)
        if coeffs[0] eq -1 then begin
            ;print,'AMOEBA failed to converge'
            continue
        endif
        ;print,'F_Factor['+strtrim(string(j),2)+']=',coeffs
        f_factor[j] = coeffs[0]
        time_f[j] = mean(timeB[pos0])
        temp54=total(solar54[pA].solar_exp_orbit * f_factor[j],/cum) / (1d +(solar54[pA].oneau-1d)/4d) / 86400d
        temp55=total(solar55[pB].solar_exp_orbit * f_factor[j],/cum) / (1d +(solar55[pB].oneau-1d)/4d) / 86400d
        solarexpA[j+1] = solarexpA[j] + temp54[-1]
        solarexpB[j+1] = solarexpB[j] + temp55[-1]
    endfor

    p=where(f_factor ne -1d and time_f ne -1d,count)
    if count eq 0 then begin
        print,'No matching VISA and VISB data '
        return,-1
    endif
    f_factor=f_factor[p]
    time_f=time_f[p]

    k=where(abs(f_factor) lt 3.0,count)
    if count gt 10 then begin
        coeff6=robust_poly_fit(time_f[k],f_factor[k],6,/double)
        resistant_mean,f_factor-poly(time_f,coeff6),3.0,mean,goodvec=keep0
    endif else begin
        k=lindgen(n_elements(f_factor))
        keep0=k
    endelse
    coeff2=robust_poly_fit(time_f[k[keep0]],f_factor[k[keep0]],2,/double)
    coeff3=robust_poly_fit(time_f[k[keep0]],f_factor[k[keep0]],3,/double)
    coeff4=robust_poly_fit(time_f[k[keep0]],f_factor[k[keep0]],4,/double)
    coeff6=robust_poly_fit(time_f[k[keep0]],f_factor[k[keep0]],6,/double)
    ; add points to the beginning of mission to smooth out the bspline
    mn=min(time_f[k[keep0]])
    xx=dindgen(400)
    yy=poly(xx,coeff2)
    sset=bspline_iterfit([xx,time_f[k[keep0[1:-1]]]],gauss_smooth([yy,f_factor[k[keep0[1:-1]]]],width=20,/edge_mirror),maxiter=0,requiren=6,bkspace=5)

    if not keyword_set(noplot) then begin
        ;keep0=lindgen(n_elements(f_factor))
        ;label='F_Factor @ '+strtrim(string(wavelength,format='(F0.2)'),2)+'nm  '
        ;label=label+'Kappa='+strtrim(string(mean(kappadeg),format='(E0.4)'),2)
        ;if n_elements(step55) gt 0 then label=label+'   Step='+strtrim(string(step55,format='(F0.1)'),2)+' days'
        ;plot_multi,time_f[keep0],f_factor[keep0]/kappa[keep0],psym=-4,title=label, xtitle='Mission Day',ytitle='F/Kappa',$
        ;    yrange=[0,2],xrange=[0,4000], charsize=1.4
        if n_elements(step55) gt 0 then title='Step='+strtrim(string(step55,format='(F0.1)'),2)+' days' else title=''
        title += '  '+strtrim(string(wavelength,format='(F0.2)'),2)+'nm'
        ;lineplot,time_f[k[keep0]],f_factor[k[keep0]],psym=-3,ptitle='F_Factor', xtitle='Mission Day',ytitle='F_Factor',$
        ;    xrange=gps2sd([t0,t1]/1d6), charsize=1.4,title=title,font=-1
        mx=long(max(time_f[k[keep0]]))
        xx=dindgen(mx)
        ;lineplot,xx,poly(xx,coeff3), psym=-3, thick=2.0, title=title+" 3rd",font=-1
        ;lineplot,xx,poly(xx,coeff4), psym=-3, thick=2.0, title=title+" 4th",font=-1
        ;lineplot,xx,bspline_valu(xx,sset), psym=-3, thick=2.0, title=title+" BSLINE",font=-1
        plot_multi,time_f[k[keep0]],gauss_smooth(f_factor[k[keep0]],width=20,/edge_mirror),xx,poly(xx,coeff3),$
            xx,bspline_valu(xx,sset),/xst,/yst,psym=[-3,-3,-3], yrange=[-1d,3], $
            xtitle='Mission Day',ytitle='F_Factor',title=title,color=[4,5,6],xrange=gps2sd([t0,t1]/1d6), charsize=1.4
    endif

    return, {timeA:timeA[k[keep0]], timeB:timeB[k[keep0]], irradA:irradA[k[keep0]], irradB:irradB[k[keep0]], $
             solarexpA:solarexpA[k[keep0]], solarexpB:solarexpB[k[keep0]], coeffs:coeff4, sset:sset, $
             time_f:time_f[k[keep0]], f_factor:f_factor[k[keep0]], f_factor_scaled:f_factor[k[keep0]]*kappadeg[0]}

end

