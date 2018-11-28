;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Estimates the degradation by comparing the data from
;   the ESR detector from SimA and SimB.
;
; CALLING SEQUENCE:
;
; INPUT PARAMETERS:
;   startTime -
;      The lower time range for which data will be returned.
;      Expects times to be in mission days.
;   stopTime -
;      The upper time range for which data will be returned.
;      Expects times to be in mission days.
;   wavelength - 
;      Wavelength to process.
;
; OPTIONAL INPUT PARAMETERS:
;   version -
;      Requests a specific version of the SimProfileIntegralCal, otherwise
;      pick the latest version valid for the specified time range.
;
; RETURNED PARAMETERS:
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

;
; REVISION HISTORY:
;   Revision: $Id: tnmin_esrtau2.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;------------------------------------------------------------------
function mytnmin_esrtau2_amoeba, x
    common my_esrtau2_common, subSolA, subIrdA, subSolB, subIrdB, afacta_val, afactb_val, kappa_val
    ; measured_irrad = SolarIrrad * exp(-abs(x[0]*solexp))
    ; SolarIrrad = measured_irrad / exp(-abs(x[0]*solexp))
    prismDeg_A0 = (1d - afacta_val[0])*exp(-abs(x[0]*kappa_val*subSolA[0])) + afacta_val[0]*exp(-abs(x[0]*kappa_val*subSolA[0])/2d)
    prismDeg_A1 = (1d - afacta_val[0])*exp(-abs(x[0]*kappa_val*subSolA[1])) + afacta_val[0]*exp(-abs(x[0]*kappa_val*subSolA[1])/2d)

    prismDeg_B0 = (1d - afactb_val[0])*exp(-abs(x[0]*kappa_val*subSolB[0])) + afactb_val[0]*exp(-abs(x[0]*kappa_val*subSolB[0])/2d)
    prismDeg_B1 = (1d - afactb_val[0])*exp(-abs(x[0]*kappa_val*subSolB[1])) + afactb_val[0]*exp(-abs(x[0]*kappa_val*subSolB[1])/2d)

    ratioB = subIrdB[1] / subIrdB[0]
    ratioA = subIrdA[1] / subIrdA[0]

    out_value = ratioB/ratioA - (prismDeg_A0/prismDeg_A1) / (prismDeg_B0/prismDeg_B1)
    ;out_value = abs(mean(F))
    ;out_value = stddev(F)
    ;out_value = MEANABSDEV(F, /median)
    ;out_value=robust_sigma(F, /zero)
    ;print,x,abs(out_value)*1d6
    return, abs(out_value)*1d6
end

;*******************************************************************

function tnmin_esrtau2, wavelength, inspectrum=inspectrum, coeffs=coeffs, $
    solar54=solar54, solar55=solar55, step=step, kappa=kappa, $
    fit_goodness=fit_goodness, afact=afact, noplot=noplot, order=order

    common my_esrtau2_common

    ; get the SimSolarExposureData for modes 31 and 32
    if n_elements(solar54) eq 0 or n_elements(solar55) eq 0 then begin
        ;restore,'~/SORCE/data/everyOrbitExpos_AB.sav'
        restore,'~/SORCE/data/solarexp_54_55_modified_v20.sav'
    endif

    if size(inspectrum,/tname) ne 'STRUCT' then begin
        restore,'~/SORCE/data/sima_vis_uncorr_2011.sav'
        restore,'~/SORCE/data/simb_vis_uncorr_2011.sav'
        restore,'~/SORCE/data/sima_esr_uncorr_2011.sav'
        restore,'~/SORCE/data/simb_esr_uncorr_2011.sav'
        inspectrum={sima:ESRA_UNCORR_2011.spect20, simb:ESRB_UNCORR_2011.spect20, pda:visa_uncorr_2011.spect20, pdb:visb_uncorr_2011.spect20}
    endif

    if n_elements(afact) eq 0 then begin
        ;readcol,'~/SORCE/data/overlap_esr_fit.txt',ray_w,ray_a,format='(d,d)'
        readcol,'~/SORCE/data/overlap_esr_v20.txt',ray_w,ray_a,format='(d,d)'
        afact_a={wavelength:ray_w, a:ray_a}
        afact_b = afact_a
    endif else begin
        afact_a = afact
        afact_b = afact
    endelse

    if size(kappa,/tname) ne 'STRUCT' then begin
        ; for the VIS diode, we use the same Kappa as the ESR
        ;readcol,'~/SORCE/data/ESR_kappa_453_1570_1004.txt',pdkappa_w, pdkappa_k,format='(d,d)'
        ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2000.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
        ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2001.txt',pdkappa_w, pdkappa_k,format='(d,d)'
        ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2002.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
        ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2011.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
        readcol,'~/SORCE/data/ESR_kappa_453_1570_2011_mod_v20.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
        kappa={x:pdkappa_W, y:pdkappa_k}
    endif

    if n_elements(order) eq 0 then order=3

    nspecta=n_elements(inspectrum.sima)
    nspectb=n_elements(inspectrum.simb)

    esra=dblarr(nspecta)
    gpsa=dblarr(nspecta)

    esrb=dblarr(nspectb)
    gpsb=dblarr(nspectb)

    for i=0,nspecta-1 do begin 
        ; find the closest diode scan in time to this ESR scan
        mn=min(abs(inspectrum.pda[*].timestamp[0] - inspectrum.sima[i].timestamp[0]),posa)
        mn=min(abs(inspectrum.sima[i].wavelength - wavelength),posw)
        p=where(inspectrum.pda[posa].wavelength gt 0d)
        if min(inspectrum.pda[posa].wavelength[p],max=mx) gt wavelength or mx lt wavelength then continue
        p=where(inspectrum.sima[i].wavelength gt 0d)
        if min(inspectrum.sima[i].wavelength[p],max=mx) gt wavelength or mx lt wavelength then continue
        s=sort(inspectrum.pda[posa].wavelength) 
        ; get the diode irradiance at requested wavelength and at the closest ESR table wavelength
        tmp0 = interpol(inspectrum.pda[posa].irradiance[s], inspectrum.pda[posa].wavelength[s], wavelength) 
        tmp1 = interpol(inspectrum.pda[posa].irradiance[s], inspectrum.pda[posa].wavelength[s], inspectrum.sima[i].wavelength[posw]) 
        esra[i] = inspectrum.sima[i].irradiance[posw] + (tmp0 - tmp1)
        gpsa[i] = inspectrum.sima[i].timestamp[posw]
    endfor
    p=where(finite(esra) eq 1 and esra gt 0d and gpsa gt 0d)
    cc=robust_poly_fit(gpsa[p], esra[p],3,yfit,/double)
    resistant_mean,abs(esra[p]-yfit),3.0,my_mean,good=keep
    gpsa=gpsa[p[keep]]
    esra=esra[p[keep]]

    for i=0,nspectb-1 do begin 
        ; find the closest diode scan in time to this ESR scan
        mn=min(abs(inspectrum.pdb[*].timestamp[0] - inspectrum.simb[i].timestamp[0]),posb)
        mn=min(abs(inspectrum.simb[i].wavelength - wavelength),posw)
        p=where(inspectrum.pdb[posb].wavelength gt 0d)
        if min(inspectrum.pdb[posb].wavelength[p],max=mx) gt wavelength or mx lt wavelength then continue
        p=where(inspectrum.simb[i].wavelength gt 0d)
        if min(inspectrum.simb[i].wavelength[p],max=mx) gt wavelength or mx lt wavelength then continue
        s=sort(inspectrum.pdb[posb].wavelength) 
        ; get the diode irradiance at requested wavelength and at the closest ESR table wavelength
        tmp0 = interpol(inspectrum.pdb[posb].irradiance[s], inspectrum.pdb[posb].wavelength[s], wavelength) 
        tmp1 = interpol(inspectrum.pdb[posb].irradiance[s], inspectrum.pdb[posb].wavelength[s], inspectrum.simb[i].wavelength[posw]) 
        esrb[i] = inspectrum.simb[i].irradiance[posw] + (tmp0 - tmp1)
        gpsb[i] = inspectrum.simb[i].timestamp[posw]
    endfor
    p=where(finite(esrb) eq 1 and esrb gt 0d and gpsb gt 0d)
    cc=robust_poly_fit(gpsb[p], esrb[p],3,yfit,/double)
    resistant_mean,abs(esrb[p]-yfit),3.0,my_mean,good=keep
    gpsb=gpsb[p[keep]]
    esrb=esrb[p[keep]]

    ; average the data in 3 days increment
    irradA=esra[0]
    sda=gpsa[0]
    irradB=esrb[0]
    sdb=gpsb[0]
    avgdays=3d
    for i=1L,n_elements(gpsb)-1L do begin
        p=where(abs(gpsb[i] - gpsb) le avgdays*86400d6,count)
        if count eq 1 then begin
            tmpird=esrb[p]
            tmpsd=gpsb[p]
        endif else if count gt 1 then begin
            resistant_mean, esrb[p], 3.0, tmpird
            tmpsd=mean(gpsb[p])
        endif
        ; avoid duplicate
        if tmpsd eq sdb[-1] then continue
        sdb=[sdb,tmpsd]
        irradB=[irradB,tmpird]
    endfor

    for i=1L,n_elements(gpsa)-1L do begin
        p=where(abs(gpsa[i] - gpsa) le avgdays*86400d6,count)
        if count eq 1 then begin
            tmpird=esra[p]
            tmpsd=gpsa[p]
        endif else if count gt 1 then begin
            resistant_mean, esra[p], 3.0, tmpird
            tmpsd=mean(gpsa[p])
        endif
        ; avoid duplicate
        if tmpsd eq sda[-1] then continue
        sda=[sda,tmpsd]
        irradA=[irradA,tmpird]
    endfor

    ;interpolate the irradiance of SimA to match times of SimB
    irradA = interpol(irradA, sda, sdb)
    sda=sdb

    solar54.solar_exp=total(solar54.solar_exp_orbit,/cum) / 86400d
    solar55.solar_exp=total(solar55.solar_exp_orbit,/cum) / 86400d
    solexpA = interpol(solar54.solar_exp, solar54.t1, sda)
    solexpB = interpol(solar55.solar_exp, solar55.t1, sdb)

    afacta_wave=(interpol(afact_a.a,afact_a.wavelength, wavelength))[0]
    afactb_wave=(interpol(afact_b.a,afact_b.wavelength, wavelength))[0]
    if where(strpos(tag_names(solar55), 'ONEAU') eq 0) ge 0 then $
        solar55_oneau = interpol(solar55.oneau, solar55.t1, sdb)

    ;afacta_wave = 0.0d
    ;afactb_wave = 0.5d

    kappa_val=(interpol(kappa.y, kappa.x, wavelength, /spline))[0]

    sda=gps2sd(sda/1d6)
    sdb=gps2sd(sdb/1d6)

    ; compare data every step days apart
    if n_elements(step) eq 0 then step=3L
    pos0=0L
    pos1=pos0 + step
    fdeg=[]
    fdeg_sd=[]

    while pos1 lt n_elements(sdb)-2 do begin
        subIrdA=irradA[[pos0,pos1]]
        subIrdB=irradB[[pos0,pos1]]
        subSolA=solexpA[[pos0,pos1]]
        subSolB=solexpB[[pos0,pos1]]

        if n_elements(solar55_oneau) gt 0 then begin
            afacta_val = afacta_wave[0] * (solar55_oneau[[pos0,pos1]])^2
            afactb_val = afactb_wave[0] * (solar55_oneau[[pos0,pos1]])^2
        endif else begin
            afacta_val = afacta_wave[0]
            afactb_val = afactb_wave[0]
        endelse

        cc = AMOEBA(1d-8, P0=1d, scale=1d, FUNCTION_VALUE=fval, FUNCTION_NAME='mytnmin_esrtau2_amoeba', nmax=maxiter)
        if cc[0] ne -1 then begin
            ;print,'AMOEBA failed to converge'

            fdeg=[fdeg,abs(cc[0])]
            fdeg_sd=[fdeg_sd, (sdb[pos1]+sdb[pos0])/2d]
            ;print,fdeg_sd[-1],fdeg[-1]
        endif

        pos0+=1L
        pos1=pos0+step
    endwhile

    fdeg=[fdeg[0],fdeg]

    ;fdeg = fdeg*0d + 1d

    fdeg_sd=[sdb[0],fdeg_sd]
    p=where(fdeg gt 1d-20,count)
    resistant_mean,fdeg[p],3.0,my_mean,good=keep
    fdeg_smooth = gauss_smooth(fdeg[p[keep]],width=10,/edge_mirror)
    fdeg_sd_smooth=fdeg_sd[p[keep]]
    coeffs=robust_poly_fit(fdeg_sd_smooth,fdeg_smooth,order,/double,target_fit)

    ;SIM A
    f_function_a=interpol(fdeg_smooth,fdeg_sd_smooth,gps2sd(solar54.t1/1.D6))
    tau_x1=total(f_function_a*solar54.solar_exp_orbit,/cum) /86400d
    tau_a=interpol(tau_x1, solar54.t1, gpsa)
    if n_elements(solar55_oneau) gt 0 then afacta_val= afacta_wave[0]*(interpol(solar54.oneau,solar54.t1,gpsa))^2

    p=where(gpsa ge sd2gps(2172d)*1d6)
    ;afacta_val[p] *= (1d/afacta_val[p])
    ;afacta_val[p] += 0.25d

    new_esra = esra / ((1d - afacta_val)*exp(-abs(kappa_val*tau_a)) + afacta_val*exp(-abs(kappa_val*tau_a)/2d))
    ; try with polyfit for FFunc
    f_function_a=poly(gps2sd(solar54.t1/1.D6), coeffs)
    tau_x1=total(f_function_a*solar54.solar_exp_orbit,/cum) /86400d
    tau_a=interpol(tau_x1, solar54.t1, gpsa)
    new_esra_fit = esra / ((1d - afacta_val)*exp(-abs(kappa_val*tau_a)) + afacta_val*exp(-abs(kappa_val*tau_a)/2d))

    ;SIM B
    f_function_b=interpol(fdeg_smooth,fdeg_sd_smooth,gps2sd(solar55.t1/1.D6))
    tau_x2=total(f_function_b*solar55.solar_exp_orbit,/cum) /86400d
    tau_b=interpol(tau_x2, solar55.t1, gpsb)
    if n_elements(solar55_oneau) gt 0 then afactb_val= afactb_wave[0]*(interpol(solar55.oneau,solar55.t1,gpsb))^2

    p=where(gpsb ge sd2gps(2172d)*1d6)
    ;afactb_val[p] *= 0.0d
    ;afactb_val[p] += 0.25d

    new_esrb = esrb / ((1d - afactb_val)*exp(-abs(kappa_val*tau_b)) + afactb_val*exp(-abs(kappa_val*tau_b)/2d))
    ; try with polyfit for FFunc
    f_function_b=poly(gps2sd(solar55.t1/1.D6), coeffs)
    tau_x2=total(f_function_b*solar55.solar_exp_orbit,/cum) /86400d
    tau_b=interpol(tau_x2, solar55.t1, gpsb)
    new_esrb_fit = esrb / ((1d - afactb_val)*exp(-abs(kappa_val*tau_b)) + afactb_val*exp(-abs(kappa_val*tau_b)/2d))
    delta=(new_esrb_fit - new_esra_fit)
    p=where(finite(delta) eq 1)
    fit_goodness=1d -stddev(delta[p])/median(new_esra_fit[p])

    if NOT keyword_set(noplot) then begin
        mx=max(gps2sd(gpsa/1d6),min=mn)
        allsd=dindgen(ceil(mx)-floor(mn))+floor(mn)
        wv=strtrim(string(wavelength,format='(F0.2)'),2)+'nm'
        plot_multi, fdeg_sd,fdeg, fdeg_sd[p[keep]],fdeg[p[keep]], fdeg_sd_smooth, fdeg_smooth,$
            allsd,poly(allsd,coeffs),/xst,/yst,label=['FDeg','FDEG Cleaned','FDeg Smooth','3rd order fit'], $
            xtitle='Mission Days',ytitle='F*Kappa',title='ESR F*Kappa @ '+wv

        pa=where(gpsa gt sd2gps(670d)*1d6 and gpsa lt sd2gps(680d)*1d6,counta)
        pb=where(gpsb gt sd2gps(670d)*1d6 and gpsb lt sd2gps(680d)*1d6,countb)
        if counta gt 0 and countb gt 0 then delta1=mean(new_esra[pa])-mean(new_esrb[pb]) else delta1=0d
        if counta gt 0 and countb gt 0 then delta2=mean(new_esra_fit[pa])-mean(new_esrb_fit[pb]) else delta2=0d

        ;lineplot,gps2sd(gpsa/1d6),esra, title='ESR A Uncorrected @ '+wv, psym=-4
        ;lineplot,gps2sd(gpsb/1d6),esrb, title='ESR B Uncorrected @ '+wv, psym=-4
        ;lineplot,gps2sd(gpsa/1d6),new_esra, title='ESR A Corrected @ '+wv, psym=-3, nsum=2
        ;lineplot,gps2sd(gpsb/1d6),new_esrb+delta1, title='ESR B Corrected @ '+wv, psym=-4

        lineplot,gps2sd(gpsa/1d6),new_esra_fit, title='ESR A Corrected FFunc Modified SolarExp poly @ '+wv, psym=-3, nsum=2
        lineplot,gps2sd(gpsb/1d6),new_esrb_fit+delta2, title='ESR B Corrected FFunc Modified SolarExp poly @ '+wv, psym=-4
    endif

    nA=n_elements(esra)
    nB=n_elements(esrb)
    spA={MICROSECONDSSINCEGPSEPOCH:gpsa, wavelength:replicate(wavelength,nA), tau:tau_a, irradiance:esra, new_irrad:new_esra_fit}
    spB={MICROSECONDSSINCEGPSEPOCH:gpsb, wavelength:replicate(wavelength,nB), tau:tau_b, irradiance:esrb, new_irrad:new_esrb_fit}
    return, {simA:spA, simB:spB, coeffs:coeffs, wavelength:wavelength, goodness:fit_goodness}

end
