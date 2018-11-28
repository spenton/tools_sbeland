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
;   Revision: $Id: process_uncorr.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;------------------------------------------------------------------
function process_uncorr, instrumentModeId, starttime=starttime, endtime=endtime, inspectrum=inspectrum, $
    solar54=solar54, solar55=solar55, kappa=kappa, no_detdeg=no_detdeg, afact=afact, degcol=degcol, no_aleph=no_aleph


    ; get the SimSolarExposureData for modes 31 and 32
    if n_elements(solar54) eq 0 or n_elements(solar55) eq 0 then begin
        ;restore,'~/SORCE/data/solarexp_54_55.sav'
        ;restore,'~/SORCE/data/solarexp_54_55_modified_v20.sav'
        ;restore,'~/SORCE/data/solarexp_54_55_mod_407nm.sav'
        ;restore,'~/SORCE/data/solarexp_54mod_55.sav'
        ;restore,'~/SORCE/data/solarexp_uv_54_55_mod.sav'
        ;restore,'~/SORCE/data/solarexp_uv_54_55_adjusted_sol260nm.sav'
        restore,'~/SORCE/data/solarexp_uv_54_55_mod2.sav'
    endif

    if size(inspectrum,/tname) ne 'STRUCT' then begin
        if instrumentModeId eq 41 then  begin
            ;restore,'~/SORCE/data/sima_vis_uncorr_2011.sav'
            ;inspectrum = temporary(visa_uncorr_2011)
            restore,'~/SORCE/data/sima_vis_uncorr_2011_ccdshift.sav'
            inspectrum = temporary(visa_uncorr_2011_CCDSHIFT)
        endif else if instrumentModeId eq 43 then begin
            ;restore,'~/SORCE/data/sima_uv_uncorr_2011.sav'
            ;inspectrum = temporary(UVA_UNCORR_2011)
            restore,'~/SORCE/data/sima_uv_uncorr_2011_ccdshift.sav'
            inspectrum = temporary(UVA_UNCORR_2011_CCDSHIFT)
        endif else if instrumentModeId eq 44 then begin
            restore,'~/SORCE/data/sima_ir_uncorr_2011.sav'
            inspectrum = temporary(ira_uncorr_2011)
        endif else if instrumentModeId eq 45 then begin
            ;restore,'~/SORCE/data/simb_vis_uncorr_2011.sav'
            ;inspectrum = temporary(visb_uncorr_2011)
            restore,'~/SORCE/data/simb_vis_uncorr_2011_CCDSHIFT.sav'
            inspectrum = temporary(visb_uncorr_2011_CCDSHIFT)
        endif else if instrumentModeId eq 47 then begin
            ;restore,'~/SORCE/data/simb_uv_uncorr_2011.sav'
            ;inspectrum = temporary(UVB_UNCORR_2011)
            restore,'~/SORCE/data/simb_uv_uncorr_2011_ccdshift.sav'
            inspectrum = temporary(UVB_UNCORR_2011_CCDSHIFT)
        endif else if instrumentModeId eq 48 then begin
            restore,'~/SORCE/data/simb_ir_uncorr_2011.sav'
            inspectrum = temporary(irb_uncorr_2011)
        endif else if instrumentModeId eq 31 then begin
            restore,'~/SORCE/data/sima_esr_uncorr_2011.sav'
            inspectrum = temporary(ESRA_UNCORR_2011)
        endif else if instrumentModeId eq 32 then begin
            restore,'~/SORCE/data/simb_esr_uncorr_2011.sav'
            inspectrum = temporary(ESRB_UNCORR_2011)
        endif
    endif

    if n_elements(afact) eq 0 then begin
        if instrumentModeId eq 41 or instrumentModeId eq 45 then begin
            ;readcol,'~/SORCE/data/overlap_vis1_smooth.txt',ray_w,ray_a,format='(d,d)'
            ;readcol,'~/SORCE/data/overlap_vis1_v20.txt',ray_w,ray_a,format='(d,d)'
            ;readcol,'~/SORCE/data/overlap_vis1a_v20_2.txt',ray_w,ray_a,format='(d,d)'
            ;readcol,'~/SORCE/data/overlap_vis1a_v20_reduced55.txt',ray_w,ray_a,format='(d,d)'
            ;readcol,'~/SORCE/data/overlap_vis1_smooth_mod_v20.txt',ray_w,ray_a,format='(d,d)'
            readcol,'~/SORCE/data/overlap_vis1_mod_407nm.txt',ray_w,ray_a,format='(d,d)'
            ;readcol,'~/SORCE/data/overlap_vis1_smooth_mod2.txt',ray_w,ray_a,format='(d,d)'
            afact={wavelength:ray_w, a:ray_a}
        endif else if instrumentModeId eq 43 or instrumentModeId eq 47 then begin
            ;readcol,'~/SORCE/data/overlap_uv_smooth_pos.txt',ray_w,ray_a,format='(d,d)'
            ;readcol,'~/SORCE/data/overlap_uv_smooth_pos_15.txt',ray_w,ray_a,format='(d,d)'
            ;readcol,'~/SORCE/data/overlap_uv_smooth_pos_07.txt',ray_w,ray_a,format='(d,d)'
            ;readcol,'~/SORCE/data/overlap_uv_smooth_pos_50.txt',ray_w,ray_a,format='(d,d)'
            readcol,'~/SORCE/data/overlap_uv_smooth_pos_160.txt',ray_w,ray_a,format='(d,d)'
            afact={wavelength:ray_w, a:ray_a}
        endif else if instrumentModeId eq 44 or instrumentModeId eq 48 then begin
            ;readcol,'~/SORCE/data/overlap_ir.txt',ray_w,ray_a,format='(d,d)'
            readcol,'~/SORCE/data/overlap_esr_v20_mod2.txt',ray_w,ray_a,format='(d,d)'
            afact={wavelength:ray_w, a:ray_a}
        endif else if instrumentModeId eq 31 or instrumentModeId eq 32 then begin
            readcol,'~/SORCE/data/overlap_esr_fit.txt',ray_w,ray_a,format='(d,d)'
            ;readcol,'~/SORCE/data/overlap_esr_v20.txt',ray_w,ray_a,format='(d,d)'
            ;readcol,'~/SORCE/data/overlap_esr_v20_mod2.txt',ray_w,ray_a,format='(d,d)'
            afact={wavelength:ray_w, a:ray_a}
        endif
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
            readcol,'~/SORCE/data/UVAB_453-1570_kappa_2011_160_s5455mod2.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
            ;readcol,'~/SORCE/data/UVAB_453-1570_kappa_2011_160_s5455_sol260nm.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
            kappa={x:pdkappa_W, y:pdkappa_k}
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
            readcol,'~/SORCE/data/ESR_kappa_453_1570_2011_mod_407nm.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
            ;readcol,'~/SORCE/data/VISAB_453-1570_kappa_2011_s5455mod2.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
            kappa={x:pdkappa_W, y:pdkappa_k}
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

    if n_elements(degcol) eq 0 then begin
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
           ;restore, '~/SORCE/data/sima_uv_degcol_v20_afact160.sav'
           ;restore, '~/SORCE/data/sima_uv_degcol_v20_afact160_s54mod.sav'
           ;restore, '~/SORCE/data/sima_uv_degcol_v20_afact160_s5455mod.sav'
           ;restore, '~/SORCE/data/sima_uv_degcol_v20_afact160_s5455_kappaall.sav'
           ;restore, '~/SORCE/data/sima_uv_degcol_v20_afact160_s54mod_adjusted.sav'
           ;restore, '~/SORCE/data/sima_uv_degcol_v20_constant_adjusted_2195.sav'
           ;restore, '~/SORCE/data/sima_uv_degcol_v20_afact15_adjusted_2195.sav'
           ;restore, '~/SORCE/data/sima_uv_degcol_v20_afact160_s5455mod2.sav'
           ;restore, '~/SORCE/data/sima_uv_degcol_v20_afact160_s5455mod2_nopoly.sav'
           ;restore, '~/SORCE/data/sima_uv_degcol_v20_s5455_mod2_tst.sav'
           ;restore, '~/SORCE/data/sima_uv_degcol_v20_afact160_s5455mod2_poly.sav'
           ;restore, '~/SORCE/data/sima_uv_degcol_v20_afact160_s5455mod2_2D.sav'
           ;restore, '~/SORCE/data/sima_uv_degcol_v20_afact160_s5455mod2_2D_smooth.sav'
           restore,'~/SORCE/data/sima_uv_degcol_v20_afact160_s5455mod2_2D_allw.sav'
           degcol=degcol_c54
        endif else if instrumentModeId eq 47 then begin
           ;restore, '~/SORCE/data/simb_uv_degcol_v20.sav'
           ;restore, '~/SORCE/data/simb_uv_degcol_v20_15.sav'
           ;restore, '~/SORCE/data/simb_uv_degcol_v20_mod_15.sav'
           ;restore, '~/SORCE/data/simb_uv_degcol_v20_07.sav'
           ;restore, '~/SORCE/data/simb_uv_degcol_v20_07_2.sav'
           ;restore, '~/SORCE/data/simb_uv_degcol_v20_07_3.sav'
           ;restore, '~/SORCE/data/simb_uv_degcol_v20_20.sav'
           ;restore, '~/SORCE/data/simb_uv_degcol_v20_160_mod407.sav'
           ;restore, '~/SORCE/data/simb_uv_degcol_v20_afact160.sav'
           ;restore, '~/SORCE/data/simb_uv_degcol_v20_afact160_s54mod.sav'
           ;restore, '~/SORCE/data/simb_uv_degcol_v20_afact160_s5455mod.sav'
           ;restore, '~/SORCE/data/simb_uv_degcol_v20_afact160_s5455_kappaall.sav'
           ;restore, '~/SORCE/data/simb_uv_degcol_v20_afact160_s5455mod2.sav'
           ;restore, '~/SORCE/data/simb_uv_degcol_v20_afact160_s5455mod2_nopoly.sav'
           ;restore, '~/SORCE/data/simb_uv_degcol_v20_s5455_mod2_tst.sav'
           ;restore, '~/SORCE/data/simb_uv_degcol_v20_afact160_s5455mod2_poly.sav'
           ;restore, '~/SORCE/data/simb_uv_degcol_v20_afact160_s5455mod2_2D.sav'
           restore, '~/SORCE/data/simb_uv_degcol_v20_afact160_s5455mod2_2D_smooth.sav'
           degcol=degcol_c55
        endif else if instrumentModeId eq 41 then begin
           ;restore, '~/SORCE/data/simab_vis_ffunc_mod_407.sav'
           restore, '~/SORCE/data/sima_vis_degcol_v20_mod_407.sav'
           ;restore, '~/SORCE/data/sima_vis_degcol_v20_s5455mod2.sav'
           degcol=degcol_c54
        endif else if instrumentModeId eq 45 then begin
           restore, '~/SORCE/data/simb_vis_degcol_v20_mod_407.sav'
           ;restore, '~/SORCE/data/simb_vis_degcol_v20_s5455mod2.sav'
           degcol=degcol_c55
        endif else if instrumentModeId eq 31 then begin
           ;restore, '~/SORCE/data/sima_vis_degcol_v20_mod_407.sav'
           ;degcol=degcol_c54

                ;restore,'~/SORCE/data/solarexp_54_55_mod_407nm.sav'
                p=where(solar54.solar_exp_orbit gt 0d, count)
                q=uniq(solar54[p].t1)
                count=n_elements(q)
                coeffs=dblarr(count,4)
                coeffs[*,0] = total(solar54[p[q]].solar_exp_orbit,/cum)/86400d
                degcol={wavelength:[200d, 820d, 1400d, 1800d, 3000d], t1:solar54[p[q]].t1, coeffs:coeffs}

        endif else if instrumentModeId eq 32 then begin
           restore, '~/SORCE/data/simb_vis_degcol_v20_mod_407.sav'
           degcol=degcol_c55
        endif else if instrumentModeId eq 44 or instrumentModeId eq 48 then begin
           ; there is no ffunct, assumed uniformly 1
           ; but we need degcol_c54 and degcol_c55
           ;restore, '~/SORCE/data/simab_ir_degcol.sav'
           p=where(solar54.solar_exp_orbit gt 0d, count)
           q=uniq(solar54[p].t1)
           count=n_elements(q)
           coeffs=dblarr(count,4)
           coeffs[*,0] = total(solar54[p[q]].solar_exp_orbit,/cum)/86400d
           degcol={wavelength:[800d, 1000d, 1400d, 1800d], t1:solar54[p[q]].t1, coeffs:coeffs}
        endif
    endif
    p=where(strpos(tag_names(degcol),'DEGCOL') ge 0,degcol_2D)
    
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
    for i=0L,n_elements(inspectrum.spect20)-1 do begin 
        if ptr_valid(inspectrum.spect20[0]) then begin
            p=where((*inspectrum.spect20[i]).timestamp ge t0 and (*inspectrum.spect20[i]).timestamp lt t1,c) 
            if c gt 0 then sppos=[sppos,i] 
        endif else begin
            pos=where(inspectrum.spect20[i].timestamp gt 0d,count) 
            if count gt 0 then begin 
                p=where(inspectrum.spect20[i].timestamp[pos] ge t0 and inspectrum.spect20[i].timestamp[pos] lt t1,c) 
                if c gt 0 then sppos=[sppos,i] 
            endif 
        endelse
    endfor
    nspect=n_elements(sppos)

    if ptr_valid(inspectrum.spect20[0]) then begin
        ;create a new array of pointers since we can't modify them 
        ;inplace without modifying the input as well
        outspect=ptrarr(nspect)
        for i=0L,nspect-1L do outspect[i]=PTR_NEW((*inspectrum.spect20[sppos[i]]))
    endif else begin
        outspect = inspectrum.spect20[sppos]
    endelse
    p=where(strpos(tag_names(solar54),'ONEAU') ge 0,oneau_flag)
    n54=n_elements(solar54.solar_exp)
    n55=n_elements(solar55.solar_exp)
    sd54=gps2sd(solar54.t1/1d6)
    sd55=gps2sd(solar55.t1/1d6)

    ; process one spectra at a time
    for sp=0L, nspect-1L do begin
        if ptr_valid(outspect[sp]) then begin
            ; the pointer only has valid data
            outwave = (*outspect[sp]).wavelength
            nwpos=n_elements(outwave)
            outtime = (*outspect[sp]).timestamp
            sd = gps2sd(outtime/1d6)
            detdeg=dblarr(nwpos) + 1.0d
            print,'Processing data from day: '+strtrim(string(sd[0],format='(F0.2)'),2),format='("     ",a,"   ",$,%"\r")'
        endif else begin
            wpos = where(finite(outspect[sp].wavelength) eq 1 and outspect[sp].wavelength gt 0d,nwpos)
            if nwpos eq 0 then continue

            outwave=outspect[sp].wavelength[wpos]
            outtime=outspect[sp].timestamp[wpos]
            sd = gps2sd(outtime/1d6)
            detdeg=dblarr(nwpos) + 1.0d
            print,'Processing data from day: '+strtrim(string(sd[0],format='(F0.2)'),2),format='("     ",a,"   ",$,%"\r")'
            ; only diode degradation with VIS 
        endelse

        if not keyword_set(no_detdeg) then begin
            if instrumentModeId eq 41 or instrumentModeId eq 45 or instrumentModeId eq 44 or instrumentModeId eq 48 then $
                for i=0L,nwpos-1 do detdeg[i]=(sfit_poly(sd[i], outwave[i], detdeg_coeffs)) < 1.0d
        endif

        kappa_val = interpol(kappa.y, kappa.x, outwave)
        afact_wave = interpol(afact.a, afact.wavelength, outwave)
        if instrumentModeID eq 41 then begin
            if (oneau_flag) eq 0 then afact_val=afact_wave[0]  else $
                afact_val= afact_wave*(interpol(solar54.oneau,solar54.t1,outtime))^2 
            if degcol_2D eq 1 then begin
                pos=where(degcol.t1 le max(outtime),count)
                degCol_val = interpol(degcol.degcol[pos[-1],*], degcol.wavelength, outwave)
            endif else begin
                pos=where(degcol.t1 le max(outtime),count)
                degCol_val = poly(outwave,degcol.coeffs[pos[-1],*])
            endelse

        endif else if instrumentModeID eq 45 then begin
            if (oneau_flag) eq 0 then afact_val=afact_wave[0]  else $
                afact_val= afact_wave*(interpol(solar55.oneau,solar55.t1,outtime))^2
            if degcol_2D eq 1 then begin
                pos=where(degcol.t1 le max(outtime),count)
                degCol_val = interpol(degcol.degcol[pos[-1],*], degcol.wavelength, outwave)
            endif else begin
                pos=where(degcol.t1 le max(outtime),count)
                degCol_val = poly(outwave,degcol.coeffs[pos[-1],*])
            endelse

        endif else if instrumentModeID eq 43 then begin
            if (oneau_flag) eq 0 then afact_val=afact_wave[0]  else $
                afact_val= afact_wave*(interpol(solar54.oneau,solar54.t1,outtime))^2 
            if degcol_2D eq 1 then begin
                pos=where(degcol.t1 le max(outtime),count)
                degCol_val = interpol(degcol.degcol[pos[-1],*], degcol.wavelength, outwave)
            endif else begin
                pos=where(degcol.t1 le max(outtime),count)
                degCol_val = poly(outwave,degcol.coeffs[pos[-1],*])
            endelse

        endif else if instrumentModeID eq 47 then begin
            if (oneau_flag) eq 0 then afact_val=afact_wave[0]  else $
                afact_val= afact_wave*(interpol(solar55.oneau,solar55.t1,outtime))^2
            if degcol_2D eq 1 then begin
                pos=where(degcol.t1 le max(outtime),count)
                degCol_val = interpol(degcol.degcol[pos[-1],*], degcol.wavelength, outwave)
            endif else begin
                pos=where(degcol.t1 le max(outtime),count)
                degCol_val = poly(outwave,degcol.coeffs[pos[-1],*])
            endelse

        endif else if instrumentModeID eq 44 then begin
            if (oneau_flag) eq 0 then afact_val=afact_wave[0]  else $
                afact_val= afact_wave*(interpol(solar54.oneau,solar54.t1,outtime))^2 
            ; no ffunc: degCol comes directly from the cumulative solar exposure for that orbit
            ;pos=where(degcol.t1 le max(outspect[sp].timestamp),count)
            pos=where(solar54.t1 le max(outtime),count)
            degCol_val=replicate(solar54[pos[-1]].solar_exp,nwpos)

        endif else if instrumentModeID eq 48 then begin
            if (oneau_flag) eq 0 then afact_val=afact_wave[0]  else $
                afact_val= afact_wave*(interpol(solar55.oneau,solar55.t1,outtime))^2
            ; no ffunc: degCol comes directly from the cumulative solar exposure for that orbit
            ;pos=where(degcol.t1 le max(outspect[sp].timestamp),count)
            pos=where(solar55.t1 le max(outtime),count)
            degCol_val=replicate(solar55[pos[-1]].solar_exp,nwpos)

        endif else if instrumentModeID eq 31 then begin
            if (oneau_flag) eq 0 then afact_val=afact_wave[0]  else $
                afact_val= afact_wave*(interpol(solar54.oneau,solar54.t1,outtime[wpos]))^2
            ;afact_val=afact_wave[0]
            if degcol_2D eq 1 then begin
                pos=where(degcol.t1 le max(outtime),count)
                degCol_val = interpol(degcol.degcol[pos[-1],*], degcol.wavelength, outwave)
            endif else begin
                pos=where(degcol.t1 le max(outtime),count)
                degCol_val = poly(outwave,degcol.coeffs[pos[-1],*])
            endelse

        endif else if instrumentModeID eq 32 then begin
            if (oneau_flag) eq 0 then afact_val=afact_wave[0]  else $
                afact_val= afact_wave*(interpol(solar55.oneau,solar55.t1,outtime))^2
            if degcol_2D eq 1 then begin
                pos=where(degcol.t1 le max(outtime),count)
                degCol_val = interpol(degcol.degcol[pos[-1],*], degcol.wavelength, outwave)
            endif else begin
                pos=where(degcol.t1 le max(outtime),count)
                degCol_val = poly(outwave,degcol.coeffs[pos[-1],*])
            endelse

        endif

        prismDeg = ((1d - afact_val)*exp(-abs(kappa_val*degCol_val)) + afact_val*exp(-abs(kappa_val*degCol_val)/2d))

        if ptr_valid(outspect[sp]) then (*outspect[sp]).irradiance /= (prismDeg * detdeg)  $
        else outspect[sp].irradiance[wpos] /= (prismDeg * detdeg)

        if size(aleph,/tname) eq 'STRUCT' then begin
            value = interpol(aleph.aleph, aleph.wavelength, outspect[sp].wavelength[wpos])
            if ptr_valid(outspect[sp]) then (*outspect[sp]).irradiance[wpos]*=value  $
            else outspect[sp].irradiance[wpos]*=value
        endif

    endfor
    print,""

    return, {spect20:outspect}

end
