;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Estimates the degradation by comparing the data from
;   the specified diode detector from SimA and SimB.
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
;   instrumentModeId -
;      Instrument modes to process. Expecting a 2 element array
;      for Sima and SimB.
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
;   Revision: $Id: tnmin_pdtau2.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;------------------------------------------------------------------
function mytnmin_pdtau2_amoeba, x
    common my_pdtau2_common, subSolA, subIrdA, subSolB, subIrdB, afacta_val, afactb_val, kappa_val
    ; measured_irrad = SolarIrrad * exp(-abs(x[0]*solexp))
    ; SolarIrrad = measured_irrad / exp(-abs(x[0]*solexp))
    prismDeg_A0 = (1d - afacta_val[0])*exp(-abs(x[0]*kappa_val*subSolA[0])) + afacta_val[0]*exp(-abs(x[0]*kappa_val*subSolA[0])/2d)
    prismDeg_A1 = (1d - afacta_val[-1])*exp(-abs(x[0]*kappa_val*subSolA[1])) + afacta_val[-1]*exp(-abs(x[0]*kappa_val*subSolA[1])/2d)

    prismDeg_B0 = (1d - afactb_val[0])*exp(-abs(x[0]*kappa_val*subSolB[0])) + afactb_val[0]*exp(-abs(x[0]*kappa_val*subSolB[0])/2d)
    prismDeg_B1 = (1d - afactb_val[-1])*exp(-abs(x[0]*kappa_val*subSolB[1])) + afactb_val[-1]*exp(-abs(x[0]*kappa_val*subSolB[1])/2d)

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

function tnmin_pdtau2, wavelength, inspectrum=inspectrum, coeffs=coeffs, $
    solar54=solar54, solar55=solar55, step=step, kappa=kappa, no_detdeg=no_detdeg, $
    fit_goodness=fit_goodness, afact=afact, noplot=noplot, order=order, alignobc=alignobc, $
    sdrange=sdrange, gettau=gettau

    common my_pdtau2_common

    ; get the SimSolarExposureData for modes 31 and 32
    if n_elements(solar54) eq 0 or n_elements(solar55) eq 0 then begin
        ;restore,'~/SORCE/data/solarexp_54_55.sav'
        ;restore,'~/SORCE/data/solarexp_54_55_modified_v20.sav'
        ;restore,'~/SORCE/data/solarexp_54_55_mod_407nm.sav'
        ;restore,'~/SORCE/data/solarexp_uv_54_55_adjusted_sol260nm.sav'
    endif

    if size(inspectrum,/tname) ne 'STRUCT' then begin
        ;restore,'~/SORCE/data/sima_vis_uncorr_2011.sav'
        ;restore,'~/SORCE/data/simb_vis_uncorr_2011.sav'
        ;inspectrum={sima:visa_uncorr_2011.spect20, simb:visb_uncorr_2011.spect20}
        ;restore,'~/SORCE/data/sima_uv_uncorr_2011.sav'
        ;restore,'~/SORCE/data/simb_uv_uncorr_2011.sav'
        ;inspectrum={sima:uva_uncorr_2011.spect20, simb:uvb_uncorr_2011.spect20}
        ;restore,'~/SORCE/data/sima_uv_uncorr_2011_ccdshift.sav'
        ;restore,'~/SORCE/data/simb_uv_uncorr_2011_ccdshift.sav'
        ;inspectrum={sima:uva_uncorr_2011_ccdshift.spect20, simb:uvb_uncorr_2011_ccdshift.spect20}
    endif

    if size(afact,/tname) ne 'STRUCT' then begin
        ;readcol,'~/SORCE/data/overlap_vis1_smooth.txt',afact_w,afact_a,format='(d,d)'
        ;readcol,'~/SORCE/data/overlap_vis1_v20.txt',ray_w,ray_a,format='(d,d)'
        ;readcol,'~/SORCE/data/overlap_vis1a_v20_2.txt',ray_w,ray_a,format='(d,d)'
        ;readcol,'~/SORCE/data/overlap_vis1a_v20_reduced55.txt',ray_w,ray_a,format='(d,d)'
        ;readcol,'~/SORCE/data/overlap_vis1_smooth_mod_v20.txt',ray_w,ray_a,format='(d,d)'
        ;readcol,'~/SORCE/data/overlap_vis1_mod_407nm.txt',ray_w,ray_a,format='(d,d)'
        ;afact_a={wavelength:ray_w, a:ray_a}
        ;readcol,'~/SORCE/data/overlap_vis1b_v20.txt',ray_w,ray_a,format='(d,d)'
        ;readcol,'~/SORCE/data/overlap_vis1b_v20_2.txt',ray_w,ray_a,format='(d,d)'
        ;readcol,'~/SORCE/data/overlap_vis1a_v20_reduced55.txt',ray_w,ray_a,format='(d,d)'
        ;readcol,'~/SORCE/data/overlap_vis_smooth_mod_v20.txt',ray_w,ray_a,format='(d,d)'
        ;afact_b={wavelength:ray_w, a:ray_a}
        ;afact_b = afact_a
    endif else begin
        ; use the same raypath for SimA and SimB
        afact_a=afact
        afact_b=afact
    endelse

    if size(kappa,/tname) ne 'STRUCT' then begin
        ; for the VIS diode, we use the same Kappa as the ESR
        ;readcol,'~/SORCE/data/ESR_kappa_453_1570_1004.txt',pdkappa_w, pdkappa_k,format='(d,d)'
        ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2000.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
        ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2001.txt',pdkappa_w, pdkappa_k,format='(d,d)'
        ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2002.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
        ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2011.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
        ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2011_reduced55.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
        ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2011_mod_v20.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
        ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2011_mod_407nm.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
        ;kappa={x:pdkappa_W, y:pdkappa_k}
    endif

    if n_elements(order) eq 0 then order=3

    nspecta=n_elements(inspectrum.sima)
    nspectb=n_elements(inspectrum.simb)

    visa=dblarr(nspecta)
    gpsa=dblarr(nspecta)

    visb=dblarr(nspectb)
    gpsb=dblarr(nspectb)

    for i=0,nspecta-1 do begin 
        ; in the hybrid mode we have many segments with partial wavelength coverage
        if ptr_valid(inspectrum.sima[i]) then begin
            if min((*inspectrum.sima[i]).wavelength,max=mx) gt wavelength or mx lt wavelength then continue
            s=sort((*inspectrum.sima[i]).wavelength)
            q=uniq((*inspectrum.sima[i]).wavelength[s])
            ;visa[i] = spline(inspectrum.sima[i].wavelength[good], inspectrum.sima[i].irradiance[good], wavelength, 1d)
            visa[i]=interpol((*inspectrum.sima[i]).irradiance[s[q]], (*inspectrum.sima[i]).wavelength[s[q]], wavelength);, /spline) 
            mn=min(abs((*inspectrum.sima[i]).wavelength[s[q]]-wavelength),pos) 
            gpsa[i]=((*inspectrum.sima[i]).timestamp)[s[q[pos]]]
        endif else begin
            k=where(inspectrum.sima[i].wavelength gt 0d)
            if min(inspectrum.sima[i].wavelength[k],max=mx) gt wavelength or mx lt wavelength then continue
            good = where(inspectrum.sima[i].wavelength gt 0d)
            s=sort(inspectrum.sima[i].wavelength[good])
            q=uniq(inspectrum.sima[i].wavelength[good[s]])
            good=good[s[q]]
            ;visa[i] = spline(inspectrum.sima[i].wavelength[good], inspectrum.sima[i].irradiance[good], wavelength, 1d)
            visa[i]=interpol(inspectrum.sima[i].irradiance[good], inspectrum.sima[i].wavelength[good], wavelength);, /spline) 
            mn=min(abs(inspectrum.sima[i].wavelength[good]-wavelength),pos) 
            gpsa[i]=(inspectrum.sima[i].timestamp)[good[pos]]
        endelse
    endfor
    p=where(finite(visa) eq 1 and visa gt 0.0 and visa lt 3.0,count)
    if count lt 10 then return,-1
    cc=robust_poly_fit(gps2sd(gpsa[p]/1d6), visa[p],6,yfit,/double)
    resistant_mean,abs(visa[p]-yfit),9.0,my_mean,good=keep
    gpsa=gpsa[p[keep]]
    visa=visa[p[keep]]
    sda=gps2sd(gpsa)

    for i=0,nspectb-1 do begin 
        if ptr_valid(inspectrum.simb[i]) then begin
            if min((*inspectrum.simb[i]).wavelength,max=mx) gt wavelength or mx lt wavelength then continue
            s=sort((*inspectrum.simb[i]).wavelength)
            q=uniq((*inspectrum.simb[i]).wavelength[s])
            visb[i]=interpol((*inspectrum.simb[i]).irradiance[s[q]], (*inspectrum.simb[i]).wavelength[s[q]], wavelength);, /spline) 
            mn=min(abs((*inspectrum.simb[i]).wavelength[s[q]]-wavelength),pos) 
            gpsb[i]=((*inspectrum.simb[i]).timestamp)[s[q[pos]]]
        endif else begin
            k=where(inspectrum.simb[i].wavelength gt 0d)
            if min(inspectrum.simb[i].wavelength[k],max=mx) gt wavelength or mx lt wavelength then continue
            good = where(inspectrum.simb[i].wavelength gt 0d)
            s=sort(inspectrum.simb[i].wavelength[good])
            q=uniq(inspectrum.simb[i].wavelength[good[s]])
            good=good[s[q]]
            ;visb[i] = spline(inspectrum.simb[i].wavelength[good], inspectrum.simb[i].irradiance[good], wavelength, 1d)
            visb[i]=interpol(inspectrum.simb[i].irradiance[good], inspectrum.simb[i].wavelength[good], wavelength);, /spline) 
            mn=min(abs(inspectrum.simb[i].wavelength[good]-wavelength),pos) 
            gpsb[i]=(inspectrum.simb[i].timestamp)[good[pos]]
        endelse
    endfor
    p=where(finite(visb) eq 1 and visb gt 0.0 and visb lt 3.0,count)
    if count lt 10 then return,-1
    ; robust_poly_fit doesn't handle very large numbers very well - convert to mission days
    cc=robust_poly_fit(gps2sd(gpsb[p]/1d6), visb[p],6,yfit,/double)
    resistant_mean,abs(visb[p]-yfit),9.0,my_mean,good=keep
    gpsb=gpsb[p[keep]]
    visb=visb[p[keep]]
    sdb=gps2sd(gpsb)

    ;resample sima to same time stamps as simb
    visa_tmp = interpol(gauss_smooth(visa,2.0,/edge_mirror), sda, sdb)
    visb_tmp = gauss_smooth(visb,2.0,/edge_mirror)

    ; average the data in 3 days increment
    irradA=visa[0]
    sda=gpsa[0]
    irradB=visb[0]
    sdb=gpsb[0]
    avgdays=3.0d

    for i=1L,n_elements(gpsb)-1L do begin
        p=where(abs(gpsb[i] - gpsb) le avgdays*86400d6,count)
        if count eq 1 then begin
            tmpirdA=visa_tmp[p]
            tmpirdB=visb_tmp[p]
            tmpsd=gpsb[p]
        endif else if count gt 1 then begin
            resistant_mean, visa_tmp[p], 3.0, tmpirdA
            resistant_mean, visb_tmp[p], 3.0, tmpirdB
            tmpsd=mean(gpsb[p],/double)
        endif
        ; avoid duplicate
        if tmpsd eq sdb[-1] then continue
        sdb=[sdb,tmpsd]
        irradB=[irradB,tmpirdB]
        irradA=[irradA,tmpirdA]
    endfor

    ;for i=1L,n_elements(gpsa)-1L do begin
    ;    p=where(abs(gpsa[i] - gpsa) le avgdays*86400d6,count)
    ;    if count eq 1 then begin
    ;        tmpird=visa[p]
    ;        tmpsd=gpsa[p]
    ;    endif else if count gt 1 then begin
    ;        resistant_mean, visa[p], 3.0, tmpird
    ;        tmpsd=mean(gpsa[p],/double)
    ;    endif
    ;    ; avoid duplicate
    ;    if tmpsd eq sda[-1] then continue
    ;    sda=[sda,tmpsd]
    ;    irradA=[irradA,tmpird]
    ;endfor

    ;interpolate the irradiance of SimA to match times of SimB
    ;irradA = spline(sda, irradA, sdb, 1d)
    ;irradA = interpol(irradA, sda, sdb)
    sda=sdb

    solar54.solar_exp=total(solar54.solar_exp_orbit,/cum) ;/86400d
    solar55.solar_exp=total(solar55.solar_exp_orbit,/cum) ;/86400d
    solexpA = interpol(solar54.solar_exp/86400d, solar54.t1, sda)
    solexpB = interpol(solar55.solar_exp/86400d, solar55.t1, sdb)

    afacta_wave=(interpol(afact_a.a,afact_a.wavelength, wavelength))[0]
    afactb_wave=(interpol(afact_b.a,afact_b.wavelength, wavelength))[0]
    if where(strpos(tag_names(solar55), 'ONEAU') eq 0) ge 0 then $
        solar55_oneau = interpol(solar55.oneau, solar55.t1, sdb)
     
    kappa_val=(interpol(kappa.y, kappa.x, wavelength))[0]
    sda=gps2sd(sda/1d6)
    sdb=gps2sd(sdb/1d6)

    ; compare data every step days apart
    if n_elements(step) eq 0 then step=3L
    pos0=0L
    pos1=pos0 + step
    fdeg=[]
    fdeg_sd=[]

    initial_val=1d
    if keyword_set(gettau) then begin
        ; when set, combine Kappa and FFunct into a single value (set kappa_val=1 and let x float
        initial_val=kappa_val
        kappa_val=1d
    endif

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

        cc = AMOEBA(1d-8, P0=initial_val, scale=1d, FUNCTION_VALUE=fval, FUNCTION_NAME='mytnmin_pdtau2_amoeba', nmax=maxiter)
        if cc[0] ne -1 then begin
            fdeg=[fdeg,abs(cc[0])]
            fdeg_sd=[fdeg_sd, (sdb[pos1]+sdb[pos0])/2d]
            ;print,fdeg_sd[-1],fdeg[-1]
        endif

        pos0+=1L
        pos1=pos0+step
    endwhile

    if keyword_set(no_detdeg) then begin
        detdega=1d
        detdegb=1d
    endif else begin
        ;restore,'/Users/sbeland/SORCE/data/simab_vis_deg_coeffs_afact48.sav'
        ;restore,'/Users/sbeland/SORCE/data/simab_vis_deg_coeffs_2011_new_raypath.sav'
        ;restore,'/Users/sbeland/SORCE/data/simab_vis_deg_coeffs_mod_v20.sav'
        restore,'/Users/sbeland/SORCE/data/simab_vis_deg_coeffs_mod_407nm.sav'
        detdega = (sfit_poly(gps2sd(gpsa/1d6), wavelength[0], visab_coeffs)) < 1.0d
        detdegb = (sfit_poly(gps2sd(gpsb/1d6), wavelength[0], visab_coeffs)) < 1.0d
    endelse

;print,'detector degradation:  ',detdega[-1],detdegb[-1],format='(a,f,f)'

    if n_elements(fdeg) eq 0 then return,-1

    fdeg=[fdeg[0],fdeg]

    ;fdeg = fdeg * 0d + 1d

    fdeg_sd=[sdb[0],fdeg_sd]
    p=where(finite(fdeg) eq 1 and fdeg gt 1d-20,count)
    resistant_mean,fdeg[p],3.0,my_mean,good=keep
    fdeg_smooth = gauss_smooth(fdeg[p[keep]],2.0,/edge_mirror)
    k=where(fdeg_smooth gt 0d)
    fdeg_smooth=[1.8d,fdeg_smooth,fdeg_smooth[k[-1]]]
    fdeg_sd_smooth=fdeg_sd[p[keep]]
    fdeg_sd_smooth=[0d,fdeg_sd_smooth,5200d]
    ;coeffs=robust_poly_fit(fdeg_sd[p],fdeg[p],order,/double,target_fit)
    if n_elements(sdrange) eq 2 then begin
        p=where(fdeg_sd_smooth ge sdrange[0] and fdeg_sd_smooth le sdrange[1],count)
        if count eq 0 then $
            coeffs=robust_poly_fit(fdeg_sd_smooth,fdeg_smooth,order,/double,target_fit) $
        else $
            coeffs=robust_poly_fit(fdeg_sd_smooth[p],fdeg_smooth[p],order,/double,target_fit)
    endif else begin
        coeffs=robust_poly_fit(fdeg_sd_smooth,fdeg_smooth,order,/double,target_fit)
    endelse

    ;SIM A
    f_function_a=interpol(fdeg_smooth,fdeg_sd_smooth,gps2sd(solar54.t1/1.D6))
    tau_x1=total(solar54.solar_exp_orbit*f_function_a,/cum) /86400d; / (1d +(solar54.oneau-1d)/4d)
    tau_a=interpol(tau_x1, solar54.t1, gpsa)
    if n_elements(solar55_oneau) gt 0 then afacta_val= afacta_wave[0]*(interpol(solar54.oneau,solar54.t1,gpsa))^2
    new_visa = visa / ((1d - afacta_val)*exp(-abs(kappa_val*tau_a)) + afacta_val*exp(-abs(kappa_val*tau_a)/2d)) / detdega
    ; try with polyfit for FFunc
    f_function_a=poly(gps2sd(solar54.t1/1.D6), coeffs)
    tau_x1=total(solar54.solar_exp_orbit*f_function_a,/cum) /86400d; / (1d +(solar54.oneau-1d)/4d)
    tau_a=interpol(tau_x1, solar54.t1, gpsa)
    new_visa_fit = visa / ((1d - afacta_val)*exp(-abs(kappa_val*tau_a)) + afacta_val*exp(-abs(kappa_val*tau_a)/2d)) / detdega

    ; generate the degCol for each orbit where we have a solar exposure
    p=where(solar54.solar_exp_orbit gt 0d)
    sd54=gps2sd(solar54[p].t1/1.D6)
    f_function_a=interpol(fdeg_smooth,fdeg_sd_smooth,sd54) > 0d   ; prevent weird polynomial extrapolation results
    degCol54=total(solar54[p].solar_exp_orbit*f_function_a,/cum) /86400d;


    ;SIM B
    f_function_b=interpol(fdeg_smooth,fdeg_sd_smooth,gps2sd(solar55.t1/1.D6))
    tau_x2=total(solar55.solar_exp_orbit*f_function_b,/cum) /86400d; / (1d +(solar55.oneau-1d)/4d)
    tau_b=interpol(tau_x2, solar55.t1, gpsb)
    if n_elements(solar55_oneau) gt 0 then afactb_val= afactb_wave[0]*(interpol(solar55.oneau,solar55.t1,gpsb))^2
    new_visb = visb / ((1d - afactb_val)*exp(-abs(kappa_val*tau_b)) + afactb_val*exp(-abs(kappa_val*tau_b)/2d)) / detdegb
    ; try with polyfit for FFunc
    f_function_b=poly(gps2sd(solar55.t1/1.D6), coeffs)
    tau_x2=total(solar55.solar_exp_orbit*f_function_b,/cum) /86400d; / (1d +(solar55.oneau-1d)/4d)
    tau_b=interpol(tau_x2, solar55.t1, gpsb)
    new_visb_fit = visb / ((1d - afactb_val)*exp(-abs(kappa_val*tau_b)) + afactb_val*exp(-abs(kappa_val*tau_b)/2d)) / detdegb
    delta=(new_visb_fit - new_visa_fit)
    p=where(finite(delta) eq 1)
    fit_goodness=1d -stddev(delta[p])/median(new_visa_fit[p])

    ; generate the degCol for each orbit where we have a solar exposure
    p=where(solar55.solar_exp_orbit gt 0d)
    sd55=gps2sd(solar55[p].t1/1.D6)
    f_function_b=interpol(fdeg_smooth,fdeg_sd_smooth,sd55) > 0d   ; prevent weird polynomial extrapolation results
    degCol55=total(solar55[p].solar_exp_orbit*f_function_b,/cum) /86400d;


    if NOT keyword_set(noplot) then begin
        res=label_date(date_format='%M %Y')
        mx=max(gps2sd(gpsa/1d6),min=mn)
        allsd=dindgen(ceil(mx)-floor(mn))+floor(mn)
        wv=strtrim(string(wavelength,format='(F0.2)'),2)+'nm'

        if keyword_set(gettau) then yrange=[0d,0.01d] else yrange=[0d,2d]
        plot_multi, fdeg_sd,fdeg, fdeg_sd[p[keep]],fdeg[p[keep]], fdeg_sd_smooth, fdeg_smooth,$
            allsd,poly(allsd,coeffs),/xst,/yst,label=['FDeg','FDEG Cleaned','FDeg Smooth',strtrim(string(order),2)+' order fit'], $
            xtitle='Mission Days',ytitle='F*Kappa',title='Diode F*Kappa @ '+wv, charsize=1.4, yrange=yrange
        ;plot_multi, sd2jd(fdeg_sd),fdeg, sd2jd(fdeg_sd[p[keep]]),fdeg[p[keep]], sd2jd(fdeg_sd_smooth), fdeg_smooth,$
        ;    sd2jd(allsd),poly(allsd,coeffs),/xst,/yst,label=['FDeg','FDEG Cleaned','FDeg Smooth','3rd order fit'], xminor=12, $
        ;    ytitle='F*Kappa',title='Diode F*Kappa @ '+wv, charsize=1.4, xtickinterval=730.5, yrange=[0d,2d], xtickformat='LABEL_DATE'

        pa=where(gpsa gt sd2gps(670d)*1d6 and gpsa lt sd2gps(680d)*1d6,counta)
        pb=where(gpsb gt sd2gps(670d)*1d6 and gpsb lt sd2gps(680d)*1d6,countb)
        if counta gt 0 and countb gt 0 then delta1=mean(new_visa[pa])-mean(new_visb[pb]) else delta1=0d
        if counta gt 0 and countb gt 0 then delta2=mean(new_visa_fit[pa])-mean(new_visb_fit[pb]) else delta2=0d

        ;lineplot,gps2sd(gpsa/1d6),visa, title='VIS A Uncorrected @ '+wv, psym=-4
        ;lineplot,gps2sd(gpsb/1d6),visb, title='VIS B Uncorrected @ '+wv, psym=-4
        lineplot,gps2sd(gpsa/1d6),new_visa, title='VIS A Corrected @ '+wv, psym=-3, nsum=2
        lineplot,gps2sd(gpsb/1d6),new_visb+delta1, title='VIS B Corrected @ '+wv, psym=-4

        if keyword_set(alignobc) then begin
            align_irrad, gpsa, new_visa_fit, /gps
            align_irrad, gpsb, new_visb_fit+delta2, /gps
        endif

        ;p=where(finite(new_visa_fit) eq 1)
        ;lineplot,gps2sd(gpsa[p]/1d6),new_visa_fit[p], title='VIS A Corrected FFunc poly @ '+wv, psym=-3, nsum=2
        ;p=where(finite(new_visb_fit) eq 1)
        ;lineplot,gps2sd(gpsb[p]/1d6),new_visb_fit[p]+delta2[p], title='VIS B Corrected FFunc poly @ '+wv, psym=-4
    endif

    nA=n_elements(visa)
    nB=n_elements(visb)
    spA={MICROSECONDSSINCEGPSEPOCH:gpsa, wavelength:replicate(wavelength,nA), tau:tau_a, irradiance:visa, new_irrad:new_visa}
    spB={MICROSECONDSSINCEGPSEPOCH:gpsb, wavelength:replicate(wavelength,nB), tau:tau_b, irradiance:visb, new_irrad:new_visb}
    return, {simA:spA, simB:spB, coeffs:coeffs, wavelength:wavelength, goodness:fit_goodness, $
             fdeg_sd:fdeg_sd_smooth, fdeg:fdeg_smooth, fdeg_sd_all:fdeg_sd, fdeg_all:fdeg, $
             sd54:sd54, degcol54:degcol54, sd55:sd55, degcol55:degcol55, kappa:kappa}

end
