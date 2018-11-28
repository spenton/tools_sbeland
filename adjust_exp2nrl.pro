;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Using the corrected irradiance from SORCE data, adjust the solar
;   exposure record (with the fixed kappa) to match the irradiance 
;   from the NRL2 SSI model data at the specifed wavelength.
;
; CALLING SEQUENCE:
;
; RETURNED VALUES:
;   Returns the same data structure that was provided for the exposur record
;   with the modified cummulative exposure data.
;
; INPUT PARAMETERS:
;
; OPTIONAL INPUT PARAMETERS:
;
; RETURNED PARAMETERS:
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
;
;------------------------------------------------------------------
function adjust_exp2nrl_amoeba, x
    common adjust_exp2nrl_common, subIrdA, subKappa, afact_val, nrl_irrad
    ; measured_irrad = SolarIrrad * exp(-abs(x[0]*solexp))
    ; SolarIrrad = measured_irrad / exp(-abs(x[0]*solexp))
    prismDeg_A = (1d - afact_val[0])*exp(-abs(x*subKappa[0])) + afact_val[0]*exp(-abs(x*subKappa[0])/2d)

    out_value = ABS(subIrdA/prismDeg_A - nrl_irrad)
    return, out_value*1d12
end
;*******************************************************************
function adjust_exp2nrl, instMode, indata, nrldata=nrldata, kappa=kappa, solarexp=solarexp, $
    afact=afact, wavelength=wavelength
    common adjust_exp2nrl_common

    if n_elements(nrldata) eq 0 then begin
        readcol,'~/SORCE/data/nrl2_ssi_280.csv',nrl_jd,ww,nrl_irrad,format='(d,d,d)'
        ; (days since 1610-01-01)
        nrl_sd = nrl_jd + (ymd2sd([1610d,1d,1d]))[0] + 0.5d
        nrl_spect = {timesd:temporary(nrl_sd), irradiance:temporary(nrl_irrad)}
    endif else nrl_spect=nrldata

    if n_elements(kappa) eq 0 then begin
        restore,'~/SORCE/data/simab_uv_kappa_2380.sav'
        ; use the kappa_smooth structure
        kappa=kappa_smooth
    endif

    if n_elements(solarexp) eq 0 then begin
        restore, '~/SORCE/data/solar5455_0_5050_1au.sav'
        if instMode lt 45 then solarexp=solar54 else solarexp=solar55
    endif

    if n_elements(afact) eq 0 then begin
        q1="SELECT * from CalibrationMetadata where calibrationTableName='SimRayPathDegradationParams'"
        q1=q1+" and instrumentModeId="+strtrim(string(instMode),2)+" order by version desc"
        query_database, q1,result,nrows
        calibId = result[0].calibrationSetId
        q2="select wavelength,singlePassAreaFraction as afact from SimRayPathDegradationParams where calibrationSetId="+strtrim(string(calibId),2)
        query_database,q2,afact,nrows
        ;if instMode eq 41 then afact_file = '~/SORCE/data/overlap_vis1_smooth.txt'
        ;if instMode eq 41 then afact_file = '~/SORCE/data/overlap_vis1_mod_407nm.txt'
        ;if instMode eq 43 then afact_file = '~/SORCE/data/overlap_uv_smooth_pos.txt'
        ;if instMode eq 44 then afact_file = '~/SORCE/data/overlap_ir_smooth.txt'
        ;if instMode eq 31 then afact_file = '~/SORCE/data/overlap_vis1_smooth.txt'
        ;readcol,afact_file,w0,a0,format='(d,d)'
        ;afact={wavelength:w0, afact:a0}
    endif

    if n_elements(wavelength) eq 0 then wavelength=280.5d

    ; check if indata is the full corrected irradiance dataset or 
    ; the timeseries for this wavelength
    p=where(strpos(tag_names(indata),'SPECT20') ge 0,count)
    if count ge 1 then begin
        ; get the timeseries for the requested wavelength
        res=compare_19_20(0,1,instMode,v20='',alldata=indata,plotwave=wavelength, /noplot)
        ; only process to the extend of the time overlap
        maxsd=min([max(res.timestamp20), max(nrl_spect.timesd)])
        p=where(res.timestamp20 le maxsd)
        sorce_spect={timesd:res.timestamp20[p], irradiance:res.v20_irrad[p]}
        val=interpol(nrl_spect.irradiance, nrl_spect.timesd, sorce_spect.timesd, /quad)
        nrl_spect = {timesd:sorce_spect.timesd, irradiance:val}
    endif else begin
        ; only process to the extend of the time overlap
        maxsd=min([max(indata.timestamp20), max(nrl_spect.timesd)])
        p=where(res.timestamp20 le maxsd)
        sorce_spect={timesd:indata.timestamp20[p], irradiance:indata.v20_irrad[p]}
        val=interpol(nrl_spect.irradiance, nrl_spect.timesd, sorce_spect.timesd, /quad)
        nrl_spect={timesd:sorce_spect.timesd, irradiance:val}
    endelse

    ; get the global offset between SORCE and NRL2 data around solar minimum
    pos = where(nrl_spect.timesd gt 220d and nrl_spect.timesd le 400d,count)
    if count gt 0 then nrl_med = median(nrl_spect.irradiance[pos]) else nrl_med=median(nrl_spect.irradiance)
    
    pos = where(sorce_spect.timesd gt 220d and sorce_spect.timesd le 400d,count)
    if count gt 0 then sorce_med = median(sorce_spect.irradiance[pos]) else sorce_med=median(sorce_spect.irradiance)

    nrl_spect.irradiance += (sorce_med - nrl_med)

    ; get the initial solar exposure on those days
    solarsd =  gps2sd(solarexp.t1/1d6)
    solarexp.solar_exp = total(solarexp.solar_exp_orbit,/cum)
    sorce_solexp = interpol(solarexp.solar_exp, solarsd, sorce_spect.timesd, /quad)
    sorce_solexp /= 86400d  ; convert to days from seconds
    
    ; get the Kappa used to calculate the corrected irradiance at this wavelength
    kv=dblarr(n_elements(kappa.timesd))
    for i=0L,n_elements(kappa.timesd)-1L do kv[i] = interpol(kappa.kappa[i,*], kappa.waves, wavelength, /quad)
    kappa_val = interpol(kv, kappa.timesd, sorce_spect.timesd, /quad)

    afact_wave = interpol(afact.afact, afact.wavelength, wavelength, /spl)
    afact_values= afact_wave[0]*(interpol(solarexp.oneau,solarsd,sorce_spect.timesd))^2d

    ; now adjust the solar exposure of the SORCE spectra to match the NRL data
    new_solexp=dblarr(n_elements(sorce_spect.timesd))
    for sd=0L, n_elements(sorce_spect.timesd)-1 do begin
        afact_val=afact_values[sd]
        subKappa = kappa_val[sd]
        prismDeg = (1d - afact_val[0]) * exp(-abs(subKappa * sorce_solexp[sd])) + afact_val[0] * exp(-abs(subKappa * sorce_solexp[sd] / 2d))
        subIrdA = sorce_spect.irradiance[sd] * prismDeg
        nrl_irrad=nrl_spect.irradiance[sd]
        cc = AMOEBA(1d-20, P0=1d, scale=1d, FUNCTION_NAME='adjust_exp2nrl_amoeba',ncalls=ncalls)
        if cc[0] ne -1 then begin
            new_solexp[sd] = abs(cc[0])
            ;print,sorce_spect.timesd, new_solexp[sd]
        endif else if ~keyword_set(quiet) then print,'Failed to converge for ',sorce_spect.timesd[sd]
    endfor

    ; get the new solar exposure per orbit
    ; filter out the points outside of 5 sigma from initial fit
    if instMode lt 45 then begin
        nnodes=40
        k=where(new_solexp gt 0d,count)
        if count eq 0 then return,-1
        sset=bspline_iterfit(sorce_spect.timesd[k], new_solexp[k], requiren=nnodes,bkspace=nnodes, nord=3)
        yfit=bspline_valu(sorce_spect.timesd, sset)
        resistant_mean, (new_solexp - yfit), 5.0,mean,goodvec=keep0
        sset=bspline_iterfit(sorce_spect.timesd[keep0], new_solexp[keep0], requiren=nnodes,bkspace=nnodes, nord=3)
        ; get the new cummulative solar exposure in seconds
        yfit_cum=bspline_valu(solarsd, sset) * 86400d
    endif else begin
        k=where(new_solexp gt 0d,count)
        coeff=robust_poly_fit(sorce_spect.timesd[k], new_solexp[k],4,/double)
        yfit_cum=poly(solarsd,coeff) * 86400d
    endelse
    ; only keep the region where we actually have NRL data
    p=where(solarsd gt nrl_spect.timesd[-1],count)
    if count gt 0 then begin
        yfit_cum[p] = yfit_cum[p[0]-1] + total(solarexp[p].solar_exp_orbit,/cum)
    endif

    new_solarexp = solarexp
    k=where(new_solarexp.solar_exp_orbit gt 0d,comp=cp)
    new_solarexp[k].solar_exp=yfit_cum[k]
    for i=0L,n_elements(cp)-1L do begin 
        pos=where(solarsd[k] lt solarsd[cp[i]],count) 
        if count eq 0 then new_solarexp[cp[i]].solar_exp=0d $
        else new_solarexp[cp[i]].solar_exp=new_solarexp[k[pos[-1]]].solar_exp 
    endfor
    new_solarexp.solar_exp_orbit = [0d,new_solarexp[1:-1].solar_exp - new_solarexp[0:-2].solar_exp]

    plt1=plot(solarsd,solarexp.solar_exp/86400d,xtitle='SORCE DAY',ytitle='Exposure (days)',color='red',$
              dimension=[1000,800],font_size=14,/xst,/yst,thick=3,title='SIMA Cumulative Solar Exposure')
    plt2=plot(solarsd,new_solarexp.solar_exp/86400d,/overplot,color='blue',thick=3)
    leg=legend(target=[plt1,plt2],position=[0.35,0.25],/norm,HORIZONTAL_ALIGNMENT='LEFT',shadow=0)
    leg[0].label='RAW Solar Exposure'
    leg[1].label='Adjusted Solar Exposure to NRL2SSI @ '+strtrim(string(wavelength,format='(F0.1)'),2)+' nm'

    return,new_solarexp

end
