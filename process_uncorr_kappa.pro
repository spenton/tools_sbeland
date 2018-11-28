;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Using the list of calibration files and the uncorrected Irradiance
;   (with wavelengths already established), process the uncorrected 
;   irradiance to a Calibrated state.
;
; CALLING SEQUENCE:
;
; RETURNED VALUES:
;   Returns the same data structure that was read in with the irradiance 
;   tag now corresponding to calibrated data.
;
; INPUT PARAMETERS:
;   startTime -
;      The lower time range for which data will be returned.
;      Expects times to be in mission days.
;   stopTime -
;      The upper time range for which data will be returned.
;      Expects times to be in mission days.
;   instrumentModeId -
;      Instrument mode to process.
;
; OPTIONAL INPUT PARAMETERS:
;   version -
;      Requests a specific version of the SimProfileIntegralCal, otherwise
;      pick the latest version valid for the specified time range.
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
; REVISION HISTORY:
;   Revision: $Id: process_uncorr_kappa.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;------------------------------------------------------------------
function process_uncorr_kappa, instrumentModeId, starttime=starttime, endtime=endtime, inspectrum=inspectrum, $
    solar54=solar54, solar55=solar55, kappa=kappa, no_detdeg=no_detdeg, afact=afact, degcol=degcol, no_aleph=no_aleph


    ; get the SimSolarExposureData for modes 31 and 32
    if n_elements(solar54) eq 0 or n_elements(solar55) eq 0 then begin
        ;restore,'~/SORCE/data/solarexp_54_55.sav'
        ;restore,'~/SORCE/data/solarexp_54_55_modified_v20.sav'
        ;restore,'~/SORCE/data/solarexp_54mod_55.sav'
        ;restore,'~/SORCE/data/solarexp_uv_54_55_mod.sav'
        ;restore,'~/SORCE/data/solarexp_uv_54_55_adjusted_sol260nm.sav'
        ;restore,'~/SORCE/data/solarexp_uv_54_55_mod2.sav'
        ;restore,'~/SORCE/data/solarexp_54_55_mod_407nm.sav'
        ;restore,'~/SORCE/data/solarexp_54_55_raw.sav'
        restore,'~/SORCE/data/solar5455_0_5050_1au.sav'
    endif

    if instrumentModeId eq 41 then  begin
        if size(inspectrum,/tname) ne 'STRUCT' then begin
            ;restore,'~/SORCE/data/sima_vis_uncorr_2011.sav'
            ;inspectrum = temporary(visa_uncorr_2011)
            ;restore,'~/SORCE/data/sima_vis_uncorr_2011_ccdshift.sav'
            ;inspectrum = temporary(visa_uncorr_2011_CCDSHIFT)
            restore,'~/SORCE/data/sima_vis_uncorr_2351.sav'
            inspectrum = temporary(visa_uncorr_2351)
        endif
        solarexp=solar54
    endif else if instrumentModeId eq 43 then begin
        if size(inspectrum,/tname) ne 'STRUCT' then begin
            ;restore,'~/SORCE/data/sima_uv_uncorr_2011.sav'
            ;inspectrum = temporary(UVA_UNCORR_2011)
            ;restore,'~/SORCE/data/sima_uv_uncorr_2011_ccdshift.sav'
            ;inspectrum = temporary(UVA_UNCORR_2011_CCDSHIFT)
            ;restore,'~/SORCE/data/sima_uv_uncorr_2351.sav'
            restore,'~/SORCE/data/sima_uv_uncorr_2375.sav'
            inspectrum = temporary(UVA_UNCORR_2375)
        endif
        solarexp=solar54
    endif else if instrumentModeId eq 44 then begin
        if size(inspectrum,/tname) ne 'STRUCT' then begin
            restore,'~/SORCE/data/sima_ir_uncorr_2011.sav'
            inspectrum = temporary(ira_uncorr_2011)
        endif
        solarexp=solar54
    endif else if instrumentModeId eq 45 then begin
        if size(inspectrum,/tname) ne 'STRUCT' then begin
            ;restore,'~/SORCE/data/simb_vis_uncorr_2011.sav'
            ;inspectrum = temporary(visb_uncorr_2011)
            ;restore,'~/SORCE/data/simb_vis_uncorr_2011_CCDSHIFT.sav'
            ;inspectrum = temporary(visb_uncorr_2011_CCDSHIFT)
            restore,'~/SORCE/data/simb_vis_uncorr_2371.sav'
            inspectrum = temporary(visb_uncorr_2371)
        endif
        solarexp=solar55
    endif else if instrumentModeId eq 47 then begin
        if size(inspectrum,/tname) ne 'STRUCT' then begin
            ;restore,'~/SORCE/data/simb_uv_uncorr_2011.sav'
            ;inspectrum = temporary(UVB_UNCORR_2011)
            ;restore,'~/SORCE/data/simb_uv_uncorr_2011_ccdshift.sav'
            ;inspectrum = temporary(UVB_UNCORR_2011_CCDSHIFT)
            ;restore,'~/SORCE/data/simb_uv_uncorr_2371.sav'
            restore,'~/SORCE/data/simb_uv_uncorr_2375.sav'
            inspectrum = temporary(UVB_UNCORR_2375)
        endif
        solarexp=solar55
    endif else if instrumentModeId eq 48 then begin
        if size(inspectrum,/tname) ne 'STRUCT' then begin
            restore,'~/SORCE/data/simb_ir_uncorr_2011.sav'
            inspectrum = temporary(irb_uncorr_2011)
        endif
        solarexp=solar55
    endif else if instrumentModeId eq 31 then begin
        if size(inspectrum,/tname) ne 'STRUCT' then begin
            restore,'~/SORCE/data/sima_esr_uncorr_2011.sav'
            inspectrum = temporary(ESRA_UNCORR_2011)
        endif
        solarexp=solar54
    endif else if instrumentModeId eq 32 then begin
        if size(inspectrum,/tname) ne 'STRUCT' then begin
            restore,'~/SORCE/data/simb_esr_uncorr_2011.sav'
            inspectrum = temporary(ESRB_UNCORR_2011)
        endif
        solarexp=solar55
    endif


    if n_elements(afact) eq 0 then begin
        ;if instrumentModeId eq 41 or instrumentModeId eq 45 then begin
        ;    ;readcol,'~/SORCE/data/overlap_vis1_smooth.txt',ray_w,ray_a,format='(d,d)'
        ;    ;readcol,'~/SORCE/data/overlap_vis1_v20.txt',ray_w,ray_a,format='(d,d)'
        ;    ;readcol,'~/SORCE/data/overlap_vis1a_v20_2.txt',ray_w,ray_a,format='(d,d)'
        ;    ;readcol,'~/SORCE/data/overlap_vis1a_v20_reduced55.txt',ray_w,ray_a,format='(d,d)'
        ;    ;readcol,'~/SORCE/data/overlap_vis1_smooth_mod_v20.txt',ray_w,ray_a,format='(d,d)'
        ;    readcol,'~/SORCE/data/overlap_vis1_mod_407nm.txt',ray_w,ray_a,format='(d,d)'
        ;    ;readcol,'~/SORCE/data/overlap_vis1_smooth_mod2.txt',ray_w,ray_a,format='(d,d)'
        ;    afact={wavelength:ray_w, a:ray_a}
        ;endif else if instrumentModeId eq 43 or instrumentModeId eq 47 then begin
        ;    ;readcol,'~/SORCE/data/overlap_uv_smooth_pos.txt',ray_w,ray_a,format='(d,d)'
        ;    ;readcol,'~/SORCE/data/overlap_uv_smooth_pos_15.txt',ray_w,ray_a,format='(d,d)'
        ;    ;readcol,'~/SORCE/data/overlap_uv_smooth_pos_07.txt',ray_w,ray_a,format='(d,d)'
        ;    ;readcol,'~/SORCE/data/overlap_uv_smooth_pos_50.txt',ray_w,ray_a,format='(d,d)'
        ;    readcol,'~/SORCE/data/overlap_uv_smooth_pos_160.txt',ray_w,ray_a,format='(d,d)'
        ;    afact={wavelength:ray_w, a:ray_a}
        ;endif else if instrumentModeId eq 44 or instrumentModeId eq 48 then begin
        ;    ;readcol,'~/SORCE/data/overlap_ir.txt',ray_w,ray_a,format='(d,d)'
        ;    readcol,'~/SORCE/data/overlap_esr_v20_mod2.txt',ray_w,ray_a,format='(d,d)'
        ;    afact={wavelength:ray_w, a:ray_a}
        ;endif else if instrumentModeId eq 31 or instrumentModeId eq 32 then begin
        ;    readcol,'~/SORCE/data/overlap_esr_fit.txt',ray_w,ray_a,format='(d,d)'
        ;    ;readcol,'~/SORCE/data/overlap_esr_v20.txt',ray_w,ray_a,format='(d,d)'
        ;    ;readcol,'~/SORCE/data/overlap_esr_v20_mod2.txt',ray_w,ray_a,format='(d,d)'
        ;    afact={wavelength:ray_w, a:ray_a}
        ;endif
        q1="SELECT * from CalibrationMetadata where calibrationTableName='SimRayPathDegradationParams'"
        q1=q1+" and instrumentModeId="+strtrim(string(instrumentModeId),2)+" order by version desc"
        query_database, q1,result,nrows
        calibId = result[0].calibrationSetId
        q2="select wavelength,singlePassAreaFraction as afact from SimRayPathDegradationParams where calibrationSetId="+strtrim(string(calibId),2)
        query_database,q2,afact,nrows
    endif

    if size(kappa,/tname) ne 'STRUCT' then begin
        if instrumentModeId eq 43 or instrumentModeId eq 47 then begin
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
            ;readcol,'~/SORCE/data/UVAB_453-1570_kappa_2011_160_s5455_sol260nm.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
            ;readcol,'~/SORCE/data/UVAB_453-1570_kappa_2011_160_s5455mod2.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
            restore,'~/SORCE/data/uv_kappa_2371_oneau.sav'
            ; TODO
            kappa=kappa_uv_smooth
        endif else if instrumentModeId eq 44 or instrumentModeId eq 48 then begin
            ; use ESR kappa for IR
            ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2011.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
            ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2011_uv_aligned.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
            readcol,'~/SORCE/data/ESR_kappa_453_1570_2011_mod2.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
            kappa={x:pdkappa_w, y:pdkappa_k}
        endif else if instrumentModeId eq 41 or instrumentModeId eq 45 then begin
            ; for the VIS diode, we use the same Kappa as the ESR
            ;readcol,'~/SORCE/data/ESR_kappa_453_1570_1004.txt',pdkappa_w, pdkappa_k,format='(d,d)'
            ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2000.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
            ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2001.txt',pdkappa_w, pdkappa_k,format='(d,d)'
            ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2002.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
            ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2011.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
            ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2011_reduced55.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
            ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2011_mod_v20.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
            ;readcol,'~/SORCE/data/VISAB_453-1570_kappa_2011_s5455mod2.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
            ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2011_mod_407nm.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
            restore,'~/SORCE/data/vis_kappa_2371.sav'
            ; TODO
            kappa=kappa_vis_smooth
        endif else if instrumentModeId eq 31 or instrumentModeId eq 32 then begin
            ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2011_mod_407nm.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
            ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2011_mod2.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
            readcol,'~/SORCE/data/ESR_kappa_453_1570_2011_uv_aligned.txt',pdkappa_w, pdkappa_k,format='(d,d)'
            kappa={x:pdkappa_W, y:pdkappa_k}
        endif
    endif

    if NOT keyword_set(no_detdeg) then begin
        if instrumentModeId eq 41 or instrumentModeId eq 45 then begin
            ;restore,'/Users/sbeland/SORCE/data/simab_vis_deg_coeffs_afact48.sav'
            ;restore,'/Users/sbeland/SORCE/data/simab_vis_deg_coeffs_2011_new_raypath.sav'
            ;restore,'/Users/sbeland/SORCE/data/simab_vis_deg_coeffs_mod_v20.sav'
            restore,'/Users/sbeland/SORCE/data/simab_vis_deg_coeffs_mod_407nm.sav'
            detdeg_coeffs = temporary(visab_coeffs)
        endif else if instrumentModeId eq 44 or instrumentModeId eq 48 then begin
            restore,'~/SORCE/data/simab_ir_deg_coeffs_2011.sav'
            ;restore,'~/SORCE/data/simab_ir_deg_coeffs_2011_1000_1350nm.sav'
            detdeg_coeffs = temporary(irab_coeffs)
        endif else begin
            ;detector degradation for VIS diode only
            restore,'/Users/sbeland/SORCE/data/simab_vis_deg_coeffs_mod_407nm.sav'
            detdeg_coeffs = temporary(visab_coeffs)*0d
            detdeg_coeffs[0,0]=1d
        endelse
    endif

    
    if NOT keyword_set(no_aleph) then begin
        if instrumentModeId eq 41 then begin
            readcol,'~/SORCE/data/sima_vis_aleph_2011.txt',ww,aa,format='(d,d)'
            aleph = {wavelength:ww, aleph:aa}
        endif else if instrumentModeId eq 43 then begin
            readcol,'~/SORCE/data/sima_uv_aleph_2011.txt',ww,aa,format='(d,d)'
            aleph = {wavelength:ww, aleph:aa}
        endif else if instrumentModeId eq 44 then begin
            readcol,'~/SORCE/data/sima_ir_aleph_2011.txt',ww,aa,format='(d,d)'
            aleph = {wavelength:ww, aleph:aa}
        endif else if instrumentModeId eq 45 then begin
            readcol,'~/SORCE/data/simb_vis_aleph_2011.txt',ww,aa,format='(d,d)'
            aleph = {wavelength:ww, aleph:aa}
        endif else if instrumentModeId eq 47 then begin
            readcol,'~/SORCE/data/simb_uv_aleph_2011.txt',ww,aa,format='(d,d)'
            aleph = {wavelength:ww, aleph:aa}
        endif else if instrumentModeId eq 48 then begin
            readcol,'~/SORCE/data/simb_ir_aleph_2011.txt',ww,aa,format='(d,d)'
            aleph = {wavelength:ww, aleph:aa}
        endif else if instrumentModeId eq 31 or instrumentModeId eq 32 then begin
            aleph = {wavelength:[200d,2500d], aleph:[1d,1d]}
        endif
    endif


    if n_elements(starttime) eq 0 then starttime=0d
    if n_elements(endtime) eq 0 then endtime=5400d
    t0=sd2gps(starttime)*1d6
    t1=sd2gps(endtime)*1d6
    ; the uv input structure is padded at the end with 0d timestamps - look before the end
    ;if instrumentModeId eq 41 or instrumentModeId eq 45 then pos=-4 else pos=-8
    ;sppos = where(inspectrum.spect20.timestamp[pos] ge t0 and inspectrum.spect20.timestamp[pos] lt t1, nspect)
    sppos=[]
    outspect=[]
    if ptr_valid(inspectrum.spect20[0]) then begin
        for i=0L,n_elements(inspectrum.spect20)-1 do begin 
            p=where((*inspectrum.spect20[i]).timestamp ge t0 and (*inspectrum.spect20[i]).timestamp lt t1,c) 
            if c gt 0 then outspect = [outspect, ptr_new((*inspectrum.spect20[i]))]
        endfor
        nspect=n_elements(outspect)
    endif else begin
        for i=0L,n_elements(inspectrum.spect20)-1 do begin 
            pos=where(inspectrum.spect20[i].timestamp gt 0d,count) 
            if count gt 0 then begin 
                p=where(inspectrum.spect20[i].timestamp[pos] ge t0 and inspectrum.spect20[i].timestamp[pos] lt t1,c) 
                if c gt 0 then sppos=[sppos,i] 
            endif 
        endfor
        nspect=n_elements(sppos)
        outspect = inspectrum.spect20[sppos]
    endelse
    p=where(strpos(tag_names(solarexp),'ONEAU') ge 0,oneau_flag)
    if solarexp[-1].solar_exp gt 1000d then scale=86400d else scale=1d

    ; process one spectra at a time
    for sp=0L, nspect-1L do begin
        if ptr_valid(outspect[sp]) then begin
            time_gps = (*outspect[sp]).timestamp
            sd0 = gps2sd(time_gps/1d6)
            wv0 = (*outspect[sp]).wavelength
            nwpos=n_elements(wv0)
        endif else begin
            wpos = where(finite(outspect[sp].wavelength) eq 1 and outspect[sp].wavelength gt 0d,nwpos)
            max_sd = max(outspect[sp].timestamp[wpos])
            time_gps = outspect[sp].timestamp[wpos]
            sd0 = gps2sd(time_gps/1d6)
            wv0 = outspect[sp].wavelength[wpos]
        endelse
        if nwpos eq 0 then continue

        detdeg=dblarr(nwpos) + 1.0d
        print,'Processing data from day: '+strtrim(string(sd0[0],format='(F0.2)'),2),format='("     ",a,"   ",$,%"\r")'
        ; only diode degradation with VIS 
        if not keyword_set(no_detdeg) then begin
            if instrumentModeId eq 41 or instrumentModeId eq 45 or instrumentModeId eq 44 or instrumentModeId eq 48 then begin
                if ptr_valid(outspect[sp]) then begin
                    for i=0L,nwpos-1 do detdeg[i]=(sfit_poly(sd0[i], (*outspect[sp]).wavelength[i], detdeg_coeffs)) < 1.0d
                endif else begin
                    for i=0L,nwpos-1 do detdeg[i]=(sfit_poly(sd0[i], outspect[sp].wavelength[wpos[i]], detdeg_coeffs)) < 1.0d
                endelse
            endif
        endif

        ; get the solar exposure (in days) at the beginning of this scan
        pos=where(solarexp.t1 le max(time_gps),count)
        if count gt 0 then begin
            solexp_val = replicate(solarexp[pos[-1]].solar_exp/scale, nwpos)
        endif else begin
            solexp_val = replicate(0d, nwpos)
        endelse
        
        ; get the kappas for these wavelengths for all days in the kappa structure
        ; them interpolate for these sd0
        kv=dblarr(n_elements(kappa.timesd),n_elements(wv0))
        for i=0L,n_elements(kappa.timesd)-1L do kv[i,*] = interpol(kappa.kappa[i,*], kappa.waves, wv0, /quad)
        kappa_val=dblarr(n_elements(wv0))
        for i=0L,n_elements(wv0)-1L do kappa_val[i] = interpol(kv[*,i], kappa.timesd, sd0[i], /quad)

        afact_wave = interpol(afact.afact, afact.wavelength, wv0)
        if (oneau_flag) eq 0 then afact_val=afact_wave[0]  else $
            afact_val= afact_wave*(interpol(solarexp.oneau,solarexp.t1,time_gps))^2 

        prismDeg = ((1d - afact_val)*exp(-abs(kappa_val*solexp_val)) + afact_val*exp(-abs(kappa_val*solexp_val/2d)))

        ;  TESTING A DIFFERENT KAPPA FUNCTION
        ;prismDeg = exp(-abs((2d - afact_val)* kappa_val *  solexp_val))


        if ptr_valid(outspect[sp]) then begin
            (*outspect[sp]).irradiance /= (prismDeg * detdeg)
            if size(aleph,/tname) eq 'STRUCT' then begin
                value = interpol(aleph.aleph, aleph.wavelength, (*outspect[sp].wavelength))
                (*outspect[sp]).irradiance *= value
            endif
        endif else begin
            outspect[sp].irradiance[wpos] /= (prismDeg * detdeg)
            if size(aleph,/tname) eq 'STRUCT' then begin
                value = interpol(aleph.aleph, aleph.wavelength, outspect[sp].wavelength[wpos])
                outspect[sp].irradiance[wpos] *= value
            endif
        endelse

    endfor
    print,""

    return, {spect20:outspect}

end
