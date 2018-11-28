;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Using an existing version of the uv diode corected irradiance data, 
;   modify the degradation column (solar exposure * FFunct), to match
;   the time series from SOLSTICE at the same day/wavelength.
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
;      Requests a specific version of the SimCalibratedIrradiance otherwise
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
;    We use the SimCorrectedIrradiance at the specified wavelength.
;    The program expects the simspect and solspect to be time series
;    at the same wavelength.
;
; REVISION HISTORY:
;   Revision: $Id: adjust_uv_exp.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;------------------------------------------------------------------
function my_adjust_uvexp_amoeba, x
    common my_adjust_uvexp_common, simIrrad, degcol_val, solirrad_val, afact_val, kappa_val, my_verbose
    ; measured_irrad = SolarIrrad * exp(-abs(x[0]*solexp))
    ; SolarIrrad = measured_irrad / exp(-abs(x[0]*solexp))
    prismDeg = (1d - afact_val)*exp(-abs(x[0]*kappa_val*degcol_val)) + afact_val*exp(-abs(x[0]*kappa_val*degcol_val)/2d)
    corrIrrad = simIrrad / prismDeg 

    out_value = abs(corrIrrad - solirrad_val)
    if keyword_set(my_verbose) then print,x[0],out_value[0]
    return, out_value
end
;*******************************************************************

function my_adjust_uvexp_tnmin, x, df=df, simIrrad=simIrrad, degcol_val=degcol_val, $
    solIrrad_val=solIrrad_val, afact_val=afact_val, kappa_val=kappa_val, $
    corrIrrad=corrIrrad, my_verbose=my_verbose

    ; measured_irrad = SolarIrrad * exp(-abs(x[0]*solexp))
    ; SolarIrrad = measured_irrad / exp(-abs(x[0]*solexp))
    prismDeg = (1d - afact_val)*exp(-abs((1d + x[0])*kappa_val*degcol_val)) + $
        afact_val*exp(-abs((1d + x[0])*kappa_val*degcol_val)/2d)
    corrIrrad = simIrrad / prismDeg 

    out_value = abs(corrIrrad - solirrad_val)
    ;if keyword_set(my_verbose) then print,x[0],out_value[0]
    return, out_value
end
;*******************************************************************

function adjust_uv_exp, instrumentModeId, wavelength, simspect=simspect,  solspect=solspect, degcol=degcol, $
    afact=afact, kappa=kappa, noplot=noplot, verbose=verbose, corrsim=corrsim
    common my_adjust_uvexp_common


    modes = [43,47]
    match,modes,instrumentModeId,suba,subb,count=count
    if count eq 0 then begin
        print,'Error: wrong instrumentModeId (expecting ',modes,')'
        return,-1
    endif

    new_wave=wavelength

    if keyword_set(verbose) then my_verbose=1 else my_verbose=0
    if size(solspect,/tname) ne 'STRUCT' then begin
        ; get the SOLSTICE data
        ;readcol,'~/SORCE/data/solstice_9_v12_level3.txt',solsd,solwave,solirrad,format='(d,d,d)'
        ;readcol,'~/SORCE/data/solstice_9_v13_level3.txt',solsd,solwave,solirrad,format='(d,d,d)'
        ; get the closest wavelength to the requested one
        ;tmp=min(abs(solwave[0:1000] - wavelength), pos)
        ;new_wave=solwave[pos[0]]
        ;p=where(solwave eq new_wave)
        ;solspect={timestamp:solsd[p], wavelength:solwave[p], irradiance:solirrad[p]}
        restore,'~/SORCE/data/solstice_m9_v13.sav'
        solspect = compare_19_20(0,1,43,/skip19,v20=13,alldata=SOLSTICE_M9_V13, plotwave=new_wave,/noplot)
    endif

    new_wave=solspect.plotwave
    print,'Processing closest SOSLTICE wavelength @ '+strtrim(string(new_wave,format='(d0.3)'),2)

    if size(simspect,/tname) ne 'STRUCT' then begin
        if instrumentModeId eq 43 then begin
            restore,'~/SORCE/data/sima_uv_uncorr_2011_ccdshift.sav'
            ;inspectrum = temporary(UVA_UNCORR_2011_CCDSHIFT)
            simspect = compare_19_20(0,1,43,/skip19,v20=2011,alldata=UVA_UNCORR_2011_CCDSHIFT, plotwave=new_wave,/noplot)

        endif else if instrumentModeId eq 47 then begin
            ;restore,'~/SORCE/data/simb_uv_uncorr_2011.sav'
            ;inspectrum = temporary(UVB_UNCORR_2011)
            restore,'~/SORCE/data/simb_uv_uncorr_2011_ccdshift.sav'
            ;inspectrum = temporary(UVB_UNCORR_2011_CCDSHIFT)
            simspect = compare_19_20(0,1,47,/skip19,v20=2011,alldata=UVB_UNCORR_2011_CCDSHIFT, plotwave=new_wave,/noplot)
        endif
    endif

    if n_elements(afact) eq 0 then begin
        ;readcol,'~/SORCE/data/overlap_uv_smooth_pos.txt',ray_w,ray_a,format='(d,d)'
        ;readcol,'~/SORCE/data/overlap_uv_smooth_pos_15.txt',ray_w,ray_a,format='(d,d)'
        ;readcol,'~/SORCE/data/overlap_uv_smooth_pos_07.txt',ray_w,ray_a,format='(d,d)'
        ;readcol,'~/SORCE/data/overlap_uv_smooth_pos_50.txt',ray_w,ray_a,format='(d,d)'
        readcol,'~/SORCE/data/overlap_uv_smooth_pos_160.txt',ray_w,ray_a,format='(d,d)'
        ;readcol,'~/SORCE/data/overlap_uv_smooth_pos_20.txt',ray_w,ray_a,format='(d,d)'
        afact={wavelength:ray_w, a:ray_a}
    endif

    if size(kappa,/tname) ne 'STRUCT' then begin
        ;readcol,'~/SORCE/data/UVAB_453-1570_kappa_2011.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
        ;readcol,'~/SORCE/data/UVAB_453-1570_kappa_2011_15.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
        ;readcol,'~/SORCE/data/UVAB_453-1570_kappa_2011_mod_15.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
        ;readcol,'~/SORCE/data/UVAB_453-1570_kappa_2011_07.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
        ;readcol,'~/SORCE/data/UVAB_453-1570_kappa_2011_50.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
        ;readcol,'~/SORCE/data/UVAB_453-1570_kappa_2011_160.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
        ;readcol,'~/SORCE/data/UVAB_453-1200_kappa_2011_20.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
        ;readcol,'~/SORCE/data/UVAB_453-1200_kappa_2011_afact160_mod407.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
        ;readcol,'~/SORCE/data/UVAB_453-1200_kappa_2011_afact160.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
        ;readcol,'~/SORCE/data/UVAB_453-1200_kappa_2011_160_s54mod.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
        ;readcol,'~/SORCE/data/UVAB_453-1200_kappa_2011_160_s5455mod.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
        ;readcol,'~/SORCE/data/UVAB_453-1200_kappa_2011_160_s5455mod_lin.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
        ;readcol,'~/SORCE/data/UVAB_0-3900_kappa_2011_160.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
        readcol,'~/SORCE/data/UVAB_453-1570_kappa_2011_160_s5455_sol260nm.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
        kappa={x:pdkappa_W, y:pdkappa_k}
    endif

    if size(degcol,/tname) ne 'STRUCT' then begin
       ; get the surface representing the FFunction
        if instrumentModeId eq 43 then begin
           ;restore, '~/SORCE/data/sima_uv_degcol_v20.sav'
           ;restore, '~/SORCE/data/sima_uv_degcol_v20_15.sav'
           ;restore, '~/SORCE/data/sima_uv_degcol_v20_mod_15.sav'
           ;restore, '~/SORCE/data/sima_uv_degcol_v20_07.sav'
           ;restore, '~/SORCE/data/sima_uv_degcol_v20_07_2.sav'
           ;restore, '~/SORCE/data/sima_uv_degcol_v20_07_3.sav'
           ;restore, '~/SORCE/data/sima_uv_degcol_v20_50.sav'
           ;restore, '~/SORCE/data/sima_uv_degcol_v20_15_07.sav'
           ;restore, '~/SORCE/data/sima_uv_degcol_v20_07_15.sav'
           ;restore, '~/SORCE/data/sima_uv_degcol_v20_15_00.sav'
           ;restore, '~/SORCE/data/sima_uv_degcol_v20_07_453_1200.sav'
           ;restore, '~/SORCE/data/sima_uv_degcol_v20_20.sav'
           ;restore, '~/SORCE/data/sima_uv_degcol_v20_160_mod407.sav'
           restore, '~/SORCE/data/sima_uv_degcol_v20_afact160.sav'
           ;restore, '~/SORCE/data/sima_uv_degcol_v20_afact160_s54mod.sav'
           degcol=degcol_c54
           ;restore, '~/SORCE/data/sima_uv_degcol_v20_afact160_s5455mod.sav'
           ;restore, '~/SORCE/data/sima_uv_degcol_v20_afact160_s5455_kappaall.sav'
        endif else if instrumentModeId eq 47 then begin
           ;restore, '~/SORCE/data/simb_uv_degcol_v20.sav'
           ;restore, '~/SORCE/data/simb_uv_degcol_v20_15.sav'
           ;restore, '~/SORCE/data/simb_uv_degcol_v20_mod_15.sav'
           ;restore, '~/SORCE/data/simb_uv_degcol_v20_07.sav'
           ;restore, '~/SORCE/data/simb_uv_degcol_v20_07_2.sav'
           ;restore, '~/SORCE/data/simb_uv_degcol_v20_07_3.sav'
           ;restore, '~/SORCE/data/simb_uv_degcol_v20_20.sav'
           ;restore, '~/SORCE/data/simb_uv_degcol_v20_160_mod407.sav'
           restore, '~/SORCE/data/simb_uv_degcol_v20_afact160.sav'
           degcol=degcol_c55
           ;restore, '~/SORCE/data/simb_uv_degcol_v20_afact160_s54mod.sav'
           ;restore, '~/SORCE/data/simb_uv_degcol_v20_afact160_s5455mod.sav'
           ;restore, '~/SORCE/data/simb_uv_degcol_v20_afact160_s5455_kappaall.sav'
        endif
    endif

    ; process every SIM scan
    kappa_val=interpol(kappa.y,kappa.x,new_wave,/spl)
    afact_val=interpol(afact.a, afact.wavelength, new_wave,/spl)
    degcol_sd=gps2sd(degcol.t1/1d6)

    nspect = n_elements(simspect.timestamp20)
    delta_exp=dblarr(nspect) -1d

    parinfo=replicate({value:1.0d, fixed:0, limited:[1,1], limits:[-1d,2d], step:1d-3, tnside:2},1)

    ; use these to force FFunct=1 (raw solar exposure)
    ;restore,'~/SORCE/data/solarexp_54_55.sav'
    restore,'~/SORCE/data/solarexp_uv_54_55_adjusted_sol260nm.sav'
    if instrumentModeId eq 43 then solarexp=solar54 else solarexp=solar55
    solarexp_sd = gps2sd(solarexp.t1/1d6)


    ; adjust the corrected irradiance of SOLSTICE to match the "corrected" irradiance of SIM on the first day 
    ;degcol_val=interpol(solarexp.solar_exp,solarexp.t1,sd2gps(simspect.timestamp20[0])*1d6)
    ;res=my_adjust_uvexp_tnmin(0.0, simIrrad=median(simspect.v20_irrad[0:4]), degcol_val=degcol_val, $
    ;    solIrrad_val=median(solspect.v20_irrad[0:4]), afact_val=afact_val, kappa_val=kappa_val, corrIrrad=corrIrrad) 

    ;solspect.v20_irrad += (corrIrrad[0] - median(solspect.v20_irrad[0:4])
    solspect.v20_irrad += 0.013086d   ; estimated from the SIMB "corrected" irrad on day 80 compared to SOLSTICE corr on same day

    ; for each SIM spectra, get the corresponding value of the degradation column (slow)
    for i=0L,nspect-1L do begin
        if keyword_set(verbose) then print,i+1,nspect,format='(I," / ",I,"      ",$,%"\r")'
        simIrrad = SIMSPECT.V20_IRRAD[i]
        ;pos=where(degcol_sd le max(simspect.timestamp20[i]))
        ;degcol_val = poly(new_wave,degcol.coeffs[pos[-1],*])

; set FFunct=1 by using the solar exposure for our degCol
pos=where(solarexp_sd le max(simspect.timestamp20[i]),count)
degcol_val = solarexp[pos[-1]].solar_exp

        ; find the closest SOLSTICE data point (in time) - skip if more than 1 days away
        mn=min(abs(solspect.timestamp20 - simspect.timestamp20[i]),pos)
        if mn gt 1.0d then continue
        solirrad_val = solspect.V20_IRRAD[pos[0]]

        ;delta_exp[i] = AMOEBA(1d-8, P0=1d, scale=1d, FUNCTION_VALUE=fval, FUNCTION_NAME='my_adjust_uvexp_amoeba', nmax=maxiter)
        ;if keyword_set(verbose) then print,simspect.timestamp20[i], simIrrad,solirrad_val,delta_exp[i]

        functargs = {simIrrad:simIrrad, degcol_val:degcol_val, solIrrad_val:solIrrad_val, $
            afact_val:afact_val, kappa_val:kappa_val, my_verbose:my_verbose}
        coeffs = tnmin('my_adjust_uvexp_tnmin', functargs=functargs, bestmin=f0, status=status, $
            nfev=nfev, niter=niter, autoderivative=1, errmsg=errmsg, parinfo=parinfo, /quiet)
        delta_exp[i]=coeffs[0]
        
    endfor

    if keyword_set(verbose) then print,''
    p=where(delta_exp gt -1d,count)
    if count eq 0 then return,-1

    ; filter out the outliers
    if keyword_set(verbose) then print,' filtering out outliers ...'
    ;sset=bspline_iterfit(simspect.timestamp20[p],delta_exp[p],maxiter=0,requiren=60,bkspace=5,nord=3)
    ;delta_fit=bspline_valu(simspect.timestamp20[p],sset)
    ;resistant_mean,delta_exp[p] - delta_fit,5.0,mean,good=k0
    cc=robust_poly_fit(simspect.timestamp20[p],delta_exp[p],9,yfit) 
    resistant_mean,delta_exp[p]-yfit,8.0,mean,good=k0

    sset=bspline_iterfit(simspect.timestamp20[p[k0]],delta_exp[p[k0]],maxiter=0,requiren=60,bkspace=5,nord=3)
    delta_fit=bspline_valu(simspect.timestamp20,sset)

    ; define a new solar exposure record based on this fit
    new_solarexp=solarexp
    k=where(solarexp.solar_exp_orbit gt 0d)
    solarexp_val = solarexp[k].solar_exp * ( 1d + bspline_valu(solarexp_sd[k], sset) )
    ; get the new solar exposure per orbit by starting from the end and subtracting the previous one
    new_solarexp[k].solar_exp_orbit = [0d, reverse( (reverse(solarexp_val))[0:-2] - (reverse(solarexp_val))[1:-1]) * 86400d]
    new_solarexp.solar_exp = total(new_solarexp.solar_exp_orbit,/cum) / 86400d

    if NOT keyword_set(noplot) then begin
        plot_multi,simspect.timestamp20[p], delta_exp[p], simspect.timestamp20, delta_fit,/xst,/yst,psym=[-4,-3],$
            thick=[1.0,3.0],xtitle='Mission Day',ytitle='Delta Multiplier',charsize=1.4,$
            title='UV Degradation Column Adjustment @ '+strtrim(string(new_wave,format='(F0.2)'),2)+ 'nm'
    endif

    ; adjust the degcol and recalculate a set of coefficients for each day
    if keyword_set(verbose) then print,' adjusting new degradation column coefficients ...'
    new_degcol=degcol
    for i=0L, n_elements(degcol.t1)-1L do begin
        delta=bspline_valu(degcol_sd[i],sset)
        degcol_val = poly(degcol.wavelength,degcol.coeffs[i,*]) * delta
        coeffs = poly_fit(degcol.wavelength, degcol_val,3)
        new_degcol.coeffs[i,*] = coeffs
    endfor

    ; calculate a new SimCorrectedIrradiance with the updated degradation column
    if keyword_set(verbose) then print,' calculating a new SimCorrectedIrradiance ...'
    corrsim = simspect
    for i=0L,nspect-1L do begin
        ;pos=where(degcol_sd le max(corrsim.timestamp20[i]))
        ;degcol_val = poly(new_wave,new_degcol.coeffs[pos[-1],*])

pos=where(solarexp_sd le max(corrsim.timestamp20[i]),count)
degcol_val = solarexp[pos[-1]].solar_exp * (1d + delta_fit[i])

        prismDeg = (1d - afact_val)*exp(-abs(kappa_val*degcol_val)) + afact_val*exp(-abs(kappa_val*degcol_val)/2d)
        corrsim.V20_IRRAD[i] /= prismDeg
    endfor

    if NOT keyword_set(noplot) then begin
        plot_multi, corrsim.timestamp20, smooth(corrsim.v20_irrad,2), solspect.timestamp20, solspect.V20_IRRAD,/xst,/yst,$
            psym=[-3,-3],label=['SimCorrectedIrradiance','SOLSTICE'], xtitle='Mission Day',ytitle='Irradiance', $
            title='SimCorrectedIrradiance vs SOSLTICE @ '+strtrim(string(new_wave,format='(F0.2)'),2)+ 'nm'
    endif

    return,{wavelength:new_wave, timestamp:simspect.timestamp20, delta_exp:delta_exp, $
        delta_fit:delta_fit, degcol:new_degcol, solarexp:new_solarexp}

end
