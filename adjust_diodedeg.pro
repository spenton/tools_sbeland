;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Using the data structure from compare_19_20 with the uncorrected irradiance and
;   the resulting ESR and diode Kappa, we estimate a value of the diode degradation 
;   to make the ESR and diode data the same at each wavelength throughout the specified 
;   time range.
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
;   instrumentModeId -
;      Id of photo-diode to process.
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
;    The ESR and diode data is expected to be the uncorrectedIrradiance (with the 1AU
;    correction applied) as in SORCE_SIM_V20 version 1003 (and 1004 for the ESR).
;
; REVISION HISTORY:
;   Revision: $Id: adjust_diodedeg.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;------------------------------------------------------------------
function adjust_diodedeg, starttime, stoptime, instrumentModeId, $
    gps=gps, missionDays=missionDays, julianDays=julianDays, $
    pddata=pddata, esrdata=esrdata, waverange=waverange, noplot=noplot, $
    solar54=solar54, solar55=solar55, yalign=yalign, surfacefit=surfacefit, $
    sfit_coeffs=sfit_coeffs, wavezero=wavezero, contourplot=contourplot

    modes = [41,43,44,45,47,48]
    match,modes,instrumentModeId,suba,subb,count=count
    if count eq 0 then begin
        print,'Error: wrong instrumentModeId (expecting '+modes,+')'
        return,-1
    endif

    if keyword_set(gps) then begin
       ; user specified time in gps
       t0 = gps2sd(startTime/1.d6)
       t1 = gps2sd(stopTime/1.d6)
    endif else if keyword_set(julianDays) then begin
       ; user specified time in julian days
       t0 = jd2sd(startTime)
       t1 = jd2gsd(stopTime)
    endif else begin
       ; user specified timetags in mission days (default)
       t0 = startTime
       t1 = stopTime
    endelse

    if size(esrdata,/tname) ne 'STRUCT' then begin
        ; get the ESR data (version 2000 contains only the table scan data)
        if instrumentModeId lt 45 then begin
            ;restore,'~/SORCE/data/sima_esr_uncorr_2011.sav'
            restore,'~/SORCE/data/sima_esr_uncorr_2377.sav'   ; 2377 contains only ESR table scans and IRscans
            esrindata=ESRA_UNCORR_2377
            ;restore,'~/SORCE/data/sima_esr_uncorr_2376.sav'    ; 2376 contains only ESRFull scans
            ;esrindata=ESRA_UNCORR_2376
            if ptr_valid(esrindata.SPECT20[0]) then begin
                nsegs=n_elements(esrindata.SPECT20)
                esd=dblarr(nsegs)
                for i=0L,nsegs-1L do esd[i]=(*esrindata.SPECT20[i]).timestamp[0]
                p=where(esd ge sd2gps(t0)*1d6 and esd le sd2gps(t1)*1d6,count)
            endif else begin
                p=where(esrindata.spect20[*].timestamp[0] ge sd2gps(t0)*1d6 and esrindata.spect20[*].timestamp[0] le sd2gps(t1)*1d6,count)
            endelse
            if count eq 0 then begin
                print,'Error: no ESR data within requested time range'
                return,-1
            endif
            esrdata=esrindata.spect20[p]

        endif else begin
            ;restore,'~/SORCE/data/simb_esr_uncorr_2011.sav'
            restore,'~/SORCE/data/simb_esr_uncorr_2377.sav'
            if ptr_valid(ESRB_UNCORR_2377.SPECT20[0]) then begin
                nsegs=n_elements(ESRB_UNCORR_2377.SPECT20)
                esd=dblarr(nsegs)
                for i=0L,nsegs-1L do esd[i]=(*ESRB_UNCORR_2377.SPECT20[i]).timestamp[0]
                p=where(esd ge sd2gps(t0)*1d6 and esd le sd2gps(t1)*1d6,count)
            endif else begin
                p=where(ESRB_UNCORR_2377.spect20[*].timestamp[0] ge sd2gps(t0)*1d6 and ESRB_UNCORR_2377.spect20[*].timestamp[0] le sd2gps(t1)*1d6,count)
            endelse
            if count eq 0 then begin
                print,'Error: no ESR data within requested time range'
                return,-1
            endif
            esrdata=ESRB_UNCORR_2377.spect20[p]

        endelse
    endif

    ; get the latest ESR Kappa
    ;readcol,'~/SORCE/data/ESR_kappa_453_1570_1004.txt',esrkappa_w, esrkappa_k,format='(d,d)'
    ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2000.txt',esrkappa_w, esrkappa_k,format='(d,d)'
    ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2001.txt',esrkappa_w, esrkappa_k,format='(d,d)'
    ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2002.txt',esrkappa_w, esrkappa_k,format='(d,d)'
    ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2011_mod_v20.txt',esrkappa_w, esrkappa_k,format='(d,d)'

    ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2011.txt',esrkappa_w, esrkappa_k,format='(d,d)'
    ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2011_mod_407nm.txt',esrkappa_w, esrkappa_k,format='(d,d)'
    readcol,'~/SORCE/data/ESR_kappa_453_1570_2011_uv_aligned.txt',esrkappa_w, esrkappa_k,format='(d,d)'

    ; get the latest ESR raypath
    ;readcol,'~/SORCE/data/overlap_esr_fit.txt',esrpath_w, tmp, esrpath_a, format='(d,d,d)'
    ;readcol,'~/SORCE/data/overlap_esr_v20.txt',esrpath_w, esrpath_a, format='(d,d)'
    readcol,'~/SORCE/data/overlap_esr_v20_mod2.txt',esrpath_w, esrpath_a, format='(d,d)'

    if instrumentModeId eq 41 then begin
        ;readcol,'~/SORCE/data/overlap_vis1.txt',pd_aw,pd_apath,format='(d,d)'
        ;readcol,'~/SORCE/data/overlap_vis1_fit.txt',pd_aw,pd_apath,format='(d,d)'
        ;readcol,'~/SORCE/data/overlap_vis1_fit_2000.txt',pd_aw,pd_apath,format='(d,d)'
        ;readcol,'~/SORCE/data/overlap_vis1a_v20_2.txt',pd_aw,pd_apath,format='(d,d)'
        ;readcol,'~/SORCE/data/overlap_vis1_smooth_mod_v20.txt',pd_aw,pd_apath,format='(d,d)'
        readcol,'~/SORCE/data/overlap_vis1_mod_407nm.txt',pd_aw,pd_apath,format='(d,d)'

        if size(pddata,/tname) ne 'STRUCT' then begin
            ;restore,'~/SORCE/data/sima_vis_1003_uncorr.sav'
            ;restore,'~/SORCE/data/sima_vis_2000_uncorr.sav'
            ;restore,'~/SORCE/data/sima_vis_uncorr_2011.sav'
            restore,'~/SORCE/data/sima_vis_uncorr_2381.sav'
            if ptr_valid(VISA_UNCORR_2381.SPECT20[0]) then begin
                nsegs=n_elements(VISA_UNCORR_2381.SPECT20)
                esd=dblarr(nsegs)
                for i=0L,nsegs-1L do esd[i]=(*VISA_UNCORR_2381.SPECT20[i]).timestamp[0]
                p=where(esd ge sd2gps(t0)*1d6 and esd le sd2gps(t1)*1d6,count)
            endif else begin
                p=where(VISA_UNCORR_2381.spect20[*].timestamp[0] ge sd2gps(t0)*1d6 and VISA_UNCORR_2381.spect20[*].timestamp[0] le sd2gps(t1)*1d6,count)
            endelse
            if count eq 0 then begin
                print,'Error: no VIS data within requested time range'
                return,-1
            endif
            pddata=VISA_UNCORR_2381.spect20[p]
        endif
        ; for the VIS diode, we use the same Kappa as the ESR
        ;readcol,'~/SORCE/data/ESR_kappa_453_1570_1004.txt',pdkappa_w, pdkappa_k,format='(d,d)'
        ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2001.txt',pdkappa_w, pdkappa_k,format='(d,d)'
        ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2002.txt',pdkappa_w, pdkappa_k,format='(d,d)'
        ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2011.txt',pdkappa_w, pdkappa_k,format='(d,d)'
        ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2011_mod_v20.txt',pdkappa_w, pdkappa_k,format='(d,d)'
        readcol,'~/SORCE/data/ESR_kappa_453_1570_2011_mod_407nm.txt',pdkappa_w, pdkappa_k,format='(d,d)'

        ; select a wavelength range to try to avoid diode degradation
        wave0=310d
        wave1=980d
        yrange=[330d,980d]
    endif else if instrumentModeId eq 43 then begin
        ;readcol,'~/SORCE/data/overlap_uv.txt',pd_aw,pd_apath,format='(d,d)'
        ;readcol,'~/SORCE/data/overlap_uv_fit.txt',pd_aw,pd_apath,format='(d,d)'
        readcol,'~/SORCE/data/overlap_uv_smooth_pos_160.txt',pd_aw,pd_apath,format='(d,d)'
        if size(pddata,/tname) ne 'STRUCT' then begin
            ;restore,'~/SORCE/data/sima_uv_1003_uncorr.sav'
            restore,'~/SORCE/data/sima_uv_uncorr_2380.sav'
            if ptr_valid(UVA_UNCORR_2380.SPECT20[0]) then begin
                nsegs=n_elements(UVA_UNCORR_2380.SPECT20)
                esd=dblarr(nsegs)
                for i=0L,nsegs-1L do esd[i]=(*UVA_UNCORR_2380.SPECT20[i]).timestamp[0]
                p=where(esd ge sd2gps(t0)*1d6 and esd le sd2gps(t1)*1d6,count)
            endif else begin
                p=where(UVA_UNCORR_2380.spect20[*].timestamp[0] ge sd2gps(t0)*1d6 and UVA_UNCORR_2380.spect20[*].timestamp[0] le sd2gps(t1)*1d6,count)
            endelse
            if count eq 0 then begin
                print,'Error: no UV data within requested time range'
                return,-1
            endif
            pddata=UVA_UNCORR_2380.spect20[p]
        endif
        ; for the UV diode, we need a different Kappa since the ESR doesn't go blue enough
        ;readcol,'~/SORCE/data/UVAB_degradation_453-1570_kappa_ffunct_V70.txt',pdkappa_w, pdkappa_k,format='(d,d)'
        readcol,'~/SORCE/data/UVAB_453-1570_kappa_2011_160_s5455mod2.txt',pdkappa_w, pdkappa_k,format='(d,d)'
        esrkappa_w=pdkappa_w
        esrkappa_k=pdkappa_k
        wave0=210d
        wave1=308d
        yrange=[wave0,wave1]
    endif else if instrumentModeId eq 44 then begin
        ;readcol,'~/SORCE/data/overlap_ir.txt',pd_aw,pd_apath,format='(d,d)'
        readcol,'~/SORCE/data/overlap_esr_v20_mod2.txt',pd_aw,pd_apath,format='(d,d)'
        if size(pddata,/tname) ne 'STRUCT' then begin
            ;restore,'~/SORCE/data/sima_ir_1003_uncorr.sav'
            ;restore,'~/SORCE/data/sima_ir_uncorr_2011.sav'
            restore,'~/SORCE/data/sima_ir_uncorr_2377.sav'
            if ptr_valid(IRA_UNCORR_2377.SPECT20[0]) then begin
                nsegs=n_elements(IRA_UNCORR_2377.SPECT20)
                esd=dblarr(nsegs)
                for i=0L,nsegs-1L do esd[i]=(*IRA_UNCORR_2377.SPECT20[i]).timestamp[0]
                p=where(esd ge sd2gps(t0)*1d6 and esd le sd2gps(t1)*1d6,count)
            endif else begin
                p=where(IRA_UNCORR_2377.spect20[*].timestamp[0] ge sd2gps(t0)*1d6 and IRA_UNCORR_2377.spect20[*].timestamp[0] le sd2gps(t1)*1d6,count)
            endelse
            if count eq 0 then begin
                print,'Error: no IR data within requested time range'
                return,-1
            endif
            pddata=IRA_UNCORR_2377.SPECT20[p]
        endif
        ; for the IR diode, we use the same Kappa as the ESR
        ;readcol,'~/SORCE/data/ESR_kappa_453_1570_1004.txt',pdkappa_w, pdkappa_k,format='(d,d)'
        ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2011.txt',pdkappa_w, pdkappa_k,format='(d,d)'
        ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2011_mod_407nm.txt',pdkappa_w, pdkappa_k,format='(d,d)'
        readcol,'~/SORCE/data/ESR_kappa_453_1570_2011_uv_aligned.txt',pdkappa_w, pdkappa_k,format='(d,d)'
        
        ; select a wavelength range to check the diode degradation
        wave0=900d
        wave1=1400d
        yrange=[wave0,wave1]
    endif else if instrumentModeId eq 45 then begin
        ;readcol,'~/SORCE/data/overlap_vis1.txt',pd_aw,pd_apath,format='(d,d)'
        ;readcol,'~/SORCE/data/overlap_vis1_fit.txt',pd_aw,pd_apath,format='(d,d)'
        ;readcol,'~/SORCE/data/overlap_vis1_fit_2000.txt',pd_aw,pd_apath,format='(d,d)'
        ;readcol,'~/SORCE/data/overlap_vis1b_v20_2.txt',pd_aw,pd_apath,format='(d,d)'
        ;readcol,'~/SORCE/data/overlap_vis1_smooth_mod_v20.txt',pd_aw,pd_apath,format='(d,d)'
        readcol,'~/SORCE/data/overlap_vis1_mod_407nm.txt',pd_aw,pd_apath,format='(d,d)'

        if size(pddata,/tname) ne 'STRUCT' then begin
            ;restore,'~/SORCE/data/simb_vis_1003_uncorr.sav'
            ;restore,'~/SORCE/data/simb_vis_2000_uncorr.sav'
            restore,'~/SORCE/data/simb_vis_uncorr_2011.sav'
            p=where(VISB_UNCORR_2011.spect20[*].timestamp[0] ge sd2gps(t0)*1d6 and VISB_UNCORR_2011.spect20[*].timestamp[0] le sd2gps(t1)*1d6,count)
            if count eq 0 then begin
                print,'Error: no VIS data within requested time range'
                return,-1
            endif
            pddata=VISB_UNCORR_2011.spect20[p]
        endif
        ; for the VIS diode, we use the same Kappa as the ESR
        ;readcol,'~/SORCE/data/ESR_kappa_453_1570_1004.txt',pdkappa_w, pdkappa_k,format='(d,d)'
        ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2002.txt',pdkappa_w, pdkappa_k,format='(d,d)'
        ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2011.txt',pdkappa_w, pdkappa_k,format='(d,d)'
        ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2011_mod_v20.txt',pdkappa_w, pdkappa_k,format='(d,d)'
        readcol,'~/SORCE/data/ESR_kappa_453_1570_2011_mod_407nm.txt',pdkappa_w, pdkappa_k,format='(d,d)'

        ; select a wavelength range to try to avoid diode degradation
        wave0=310d
        wave1=980d
        yrange=[330d,980d]
    endif else if instrumentModeId eq 47 then begin
        ;readcol,'~/SORCE/data/overlap_uv.txt',pd_aw,pd_apath,format='(d,d)'
        readcol,'~/SORCE/data/overlap_uv_fit.txt',pd_aw,pd_apath,format='(d,d)'
        if size(pddata,/tname) ne 'STRUCT' then begin
            restore,'~/SORCE/data/simb_uv_1003_uncorr.sav'
            p=where(uvb_1003.spect20[*].timestamp[0] ge sd2gps(t0)*1d6 and uvb_1003.spect20[*].timestamp[0] le sd2gps(t1)*1d6,count)
            if count eq 0 then begin
                print,'Error: no UV data within requested time range'
                return,-1
            endif
            pddata=uvb_1003.spect20[p]
        endif
        ; for the UV diode, we need a different Kappa since the ESR doesn't go blue enough
        readcol,'~/SORCE/data/UVAB_degradation_453-1570_kappa_ffunct_V70.txt',pdkappa_w, pdkappa_k,format='(d,d)'
        esrkappa_w=pdkappa_w
        esrkappa_k=pdkappa_k
        wave0=260d
        wave1=308d
        yrange=[wave0,wave1]
    endif else if instrumentModeId eq 48 then begin
        readcol,'~/SORCE/data/overlap_ir.txt',pd_aw,pd_apath,format='(d,d)'
        if size(pddata,/tname) ne 'STRUCT' then begin
            ;restore,'~/SORCE/data/simb_ir_1003_uncorr.sav'
            restore,'~/SORCE/data/simb_ir_uncorr_2011.sav'
            p=where(irb_uncorr_2011.spect20[*].timestamp[0] ge sd2gps(t0)*1d6 and irb_uncorr_2011.spect20[*].timestamp[0] le sd2gps(t1)*1d6,count)
            if count eq 0 then begin
                print,'Error: no UV data within requested time range'
                return,-1
            endif
            pddata=irb_uncorr_2011.spect20[p]
        endif
        ; for the IR diode, we use the same Kappa as the ESR
        ;readcol,'~/SORCE/data/ESR_kappa_453_1570_1004.txt',pdkappa_w, pdkappa_k,format='(d,d)'
        readcol,'~/SORCE/data/ESR_kappa_453_1570_2011.txt',pdkappa_w, pdkappa_k,format='(d,d)'
        ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2011_mod_407nm.txt',pdkappa_w, pdkappa_k,format='(d,d)'
        ; select a wavelength range to try to avoid diode degradation
        wave0=900d
        wave1=1400d
        yrange=[wave0,wave1]
    endif


    ; use the modified solar exposure record  SBeland 2014-06-05
    ;if n_elements(solar54) eq 0 then restore,'~/SORCE/data/solarexp_54_55_modified_v20.sav'
    ;if n_elements(solar54) eq 0 then restore,'~/SORCE/data/solarexp_54_55_mod_407nm.sav'
    if n_elements(solar54) eq 0 then restore,'~/SORCE/data/solarexp_54_55.sav'
    tags=tag_names(solar54)
    kk=where(strpos(tags,'ONEAU') ge 0,count)
    if count ge 1 then solar54_oneau=1 else solar54_oneau=[]

    if n_elements(solar54) eq 0 and instrumentModeId lt 45 then begin
        if keyword_set(verbose) then print,'   getting solar54 ...'
        q1 = 'SELECT orbitStartTimeGps t0, orbitEndTimeGps t1, cumulativeSolarExpHrtOut solar_exp, '+$
             'solarExposureHrtOut as solar_exp_orbit from SimSolarExposureData where instrumentModeId=54'
        query_database, q1, solar54, info
    endif else if n_elements(solar55) eq 0 and instrumentModeId ge 45 then begin
        if keyword_set(verbose) then print,'   getting solar55 ...'
        q1 = 'SELECT orbitStartTimeGps t0, orbitEndTimeGps t1, cumulativeSolarExpHrtOut solar_exp, '+$
             'solarExposureHrtOut as solar_exp_orbit from SimSolarExposureData where instrumentModeId=55'
        query_database, q1, solar55, info
    endif

    ; use 15 esrdata scans to identify the list of wavelengths (by prism positions)
    p0= long((n_elements(esrdata) / 2.0 - 8.0) > 0)
    p1=(p0+15L)<(n_elements(esrdata)-1L)
    if ptr_valid(esrdata[0]) then begin
        esr_pos=[]
        for i=p0,p1 do esr_pos=[esr_pos, (*esrdata[i]).prismposition]
        p=where(esr_pos gt 0d)
        yhist=histogram(esr_pos[p],bin=1.0,location=xhist)
        k=where(yhist ge 1.0,count)
        if count eq 0 then begin
            print,'No valid ESR data to process'
            return,-1
        endif
        esr_prismpos = xhist[k]
    endif else begin
        p=where(esrdata[p0:p1].prismposition gt 0d)
        yhist=histogram((esrdata[p0:p1].prismposition)[p],bin=1.0,location=xhist)
        k=where(yhist ge 1.0,count)
        if count eq 0 then begin
            print,'No valid ESR data to process'
            return,-1
        endif
        esr_prismpos = xhist[k]
    endelse
    if instrumentModeId lt 45 then $
        esr_prismwaves = ccd2lambda(31, esr_prismpos, 17.2d ,tref=17.2d)  $
    else $
        esr_prismwaves = ccd2lambda(32, esr_prismpos, 17.2d ,tref=17.2d)  
    
    if n_elements(waverange) eq 2 then wave0=min(waverange,max=wave1)
    p=where(esr_prismwaves ge wave0 and esr_prismwaves le wave1,count)
    esr_prismwaves=esr_prismwaves[p]
    esr_prismpos=esr_prismpos[p]
    pd_a_factor=dblarr(n_elements(esr_prismwaves))
    detdeg=dblarr(n_elements(esr_prismwaves))
    outdata = ptrarr(n_elements(esr_prismwaves))

    if ptr_valid(pddata[0]) then begin
        pdtimestamp=dblarr(n_elements(pddata))
        for i=0L,n_elements(pddata)-1L do pdtimestamp[i]=(*pddata[i]).timestamp[0]
    endif

    ; loop over every ESR table scan wavelength within our range
    for w=0L,n_elements(esr_prismpos)-1 do begin
        ; get the list of ESR points
        print,'processing prism position '+strtrim(string(esr_prismpos[w],format='(I)'),2)+'  ('+strtrim(string(esr_prismwaves[w],format='(F0.2)'),2)+'nm)'
        if ptr_valid(esrdata[0]) then begin
            esrwaves=[]
            esrirrad=[]
            esrsd=[]
            for i=0L,n_elements(esrdata)-1L do begin
                p=where((*esrdata[i]).prismposition eq esr_prismpos[w] and (*esrdata[i]).timestamp ge sd2gps(t0)*1d6 and $
                        (*esrdata[i]).timestamp le sd2gps(t1)*1d6,count_esr)
                if count_esr eq 0 then continue
                esrwaves=[esrwaves,(*esrdata[i]).wavelength[p]]
                esrirrad=[esrirrad,(*esrdata[i]).irradiance[p]]
                esrsd=[esrsd,(*esrdata[i]).timestamp[p]]
            endfor
            count_esr=n_elements(esrirrad)
        endif else begin
            p=where(esrdata.prismposition eq esr_prismpos[w] and esrdata.timestamp ge sd2gps(t0)*1d6 and esrdata.timestamp le sd2gps(t1)*1d6,count_esr)
            if count_esr eq 0 then begin
                print,'No ESR data for ',esr_prismwaves[w]
                continue
            endif
            esrwaves=(esrdata.wavelength)[p]
            esrirrad=(esrdata.irradiance)[p]
            esrsd=(esrdata.timestamp)[p]
        endelse

        ; remove spikes in ESR data
        if n_elements(esrsd) lt 4 then continue
        ;coeff=robust_poly_fit(esrsd, esrirrad, 3, /double)
        ;resistant_mean,(esrirrad-poly(esrsd,coeff)), 4.0, mean, good=keep
        ;if n_elements(keep) ne n_elements(esrsd) then begin
        ;    print,'removed '+strtrim(string(n_elements(esrsd)-n_elements(keep)),2)+' ESR points ...'
        ;    count_esr=n_elements(keep)
        ;    esrwaves=esrwaves[keep]
        ;    esrirrad=esrirrad[keep]
        ;    esrsd=esrsd[keep]
        ;endif

        kappaesr = interpol(esrkappa_k, esrkappa_w, esrwaves, /lsq)
        kappapd = interpol(pdkappa_k, pdkappa_w, esrwaves, /lsq)
        a_esr = interpol(esrpath_a, esrpath_w, esrwaves, /lsq)
        a_pd = interpol(pd_apath,pd_aw, esrwaves, /lsq)

        ; scale the raypath by the 1AU factor
        if n_elements(solar54_oneau) gt 0 then begin
            if instrumentModeId lt 45 then begin
                a_esr = a_esr[*] * ((interpol(solar54.oneau, solar54.t1, esrsd))^2)[*]
            endif else begin
                a_esr = a_esr[*] * ((interpol(solar55.oneau, solar55.t1, esrsd))^2)[*]
            endelse
        endif 

        if instrumentModeId lt 45 then $
            solexpesr=interpol(solar54.solar_exp,solar54.t1,esrsd) / 86400d $
        else $
            solexpesr=interpol(solar55.solar_exp,solar55.t1,esrsd) / 86400d
; !! not needed now
;        solexpesr /= 86400d
        correctedESRirrad = esrirrad / ( (1d - a_esr) * exp(-kappaesr * solexpesr) + a_esr * exp(-kappaesr * solexpesr / 2d) )

        pdirrad=dblarr(count_esr)
        pdsd=dblarr(count_esr)
        new_pdirrad=dblarr(count_esr)
        delta_irrad=dblarr(count_esr)
        mnw=min(esrwaves,max=mxw)

        ; get diode data for the same wavelength and timestamp
        for t=0L,count_esr-1 do begin
            if ptr_valid(pddata[0]) then begin
                mn=min(abs(pdtimestamp - esrsd[t]),pp)
                p = where((*pddata[pp]).wavelength gt 0d)
                mn=min(((*pddata[pp]).wavelength)[p], max=mx)
                if (mn gt mxw or mx lt mnw) then continue
                s=sort(((*pddata[pp]).wavelength)[p])
                pdirrad[t] = interpol(((*pddata[pp]).irradiance)[p[s]], ((*pddata[pp]).wavelength)[p[s]], esrwaves[t],/spline)
                mn=min(abs(((*pddata[pp]).timestamp)[p] - esrsd[t]),post)
                pdsd[t]=((*pddata[pp]).timestamp)[post]
            endif else begin
                mn=min(abs(pddata.timestamp - esrsd[t]),pp)
                pos=array_indices(pddata.timestamp,pp)
                ; interpolate the pd spectra to get the uncorrectedIrradiance at the same wavelength as ESR
                p = where(pddata[pos[1]].wavelength gt 0d)
                mn=min((pddata[pos[1]].wavelength)[p], max=mx)
                if (mn gt mxw or mx lt mnw) then continue
                s=sort((pddata[pos[1]].wavelength)[p])
                pdirrad[t] = interpol((pddata[pos[1]].irradiance)[p[s]], (pddata[pos[1]].wavelength)[p[s]], esrwaves[t],/spline)
                mn=min(abs((pddata[pos[1]].timestamp)[p] - esrsd[t]),post)
                pdsd[t]=(pddata[pos[1]].timestamp)[post]
            endelse
        endfor

        k=where(pdsd gt 0d,count)
        if count lt 2 then begin
            print,'no diode data for ',esr_prismwaves[w]
            continue
        endif else if count ne n_elements(pdsd) then begin
            correctedESRirrad=correctedESRirrad[k]
            kappapd=kappapd[k]
            pdirrad=pdirrad[k]
            pdsd=pdsd[k]
            esrsd=esrsd[k]
            esrirrad=esrirrad[k]
            esrwaves=esrwaves[k]
            a_pd=a_pd[k]
        endif

        ; remove spikes in diode data
        ;coeff=robust_poly_fit(pdsd, pdirrad, 3, /double)
        ;resistant_mean,(pdirrad-poly(pdsd,coeff)), 4.0, mean, good=keep
        ;if n_elements(keep) ne n_elements(pdsd) then begin
        ;    print,'removed '+strtrim(string(n_elements(pdsd)-n_elements(keep)),2)+' PD points ...'
        ;    correctedESRirrad=correctedESRirrad[keep]
        ;    kappapd=kappapd[keep]
        ;    pdirrad=pdirrad[keep]
        ;    pdsd=pdsd[keep]
        ;    esrsd=esrsd[keep]
        ;    esrirrad=esrirrad[keep]
        ;    esrwaves=esrwaves[keep]
        ;    a_pd=a_pd[keep]
        ;endif

        ; scale the raypath by the 1AU factor
        if n_elements(solar54_oneau) gt 0 then begin
            if instrumentModeId lt 45 then begin
                a_pd = a_pd[*] * ((interpol(solar54.oneau, solar54.t1, pdsd))^2)[*]
            endif else begin
                a_pd = a_pd[*] * ((interpol(solar55.oneau, solar55.t1, pdsd))^2)[*]
            endelse
        endif 

        ; now find the diode Kappa at his wavelength that will match the ESR trend
        if instrumentModeId lt 45 then $
            solexppd=interpol(solar54.solar_exp,solar54.t1,pdsd) / 86400d $
        else $
            solexppd=interpol(solar55.solar_exp,solar55.t1,pdsd) / 86400d
; !! not needed any more
;        solexppd /= 86400d

        correctedPDirrad = pdirrad / ( (1d - a_pd) * exp(-kappapd * solexppd) + a_pd * exp(-kappapd * solexppd / 2d) )
        p1=min([n_elements(correctedPDirrad), n_elements(correctedESRirrad)]) -1L
        correctedPDirrad -= median(correctedPDirrad[0:p1]-correctedESRirrad[0:p1])
        detdeg = correctedPDirrad / correctedESRirrad
        p=where(finite(detdeg) eq 1,complement=cp,count)
        if count eq 0 then continue
        temp=detdeg[p]
        ;align_irrad,pdsd[p], temp, /gps,/trend
        align_irrad,pdsd[p], temp, /gps
        ; fit a plynomial
        coeff=robust_poly_fit(pdsd[p],temp,2,/double)
        ; force the degradation to be 0.0 on SD=0
        detdeg_fit=poly(pdsd,coeff)
        ; fit n exponential
        ;coeff = exponential_fit(gps2sd(pdsd[p]/1d6),temp)
        ;detdeg_fit = my_expo(gps2sd(pdsd/1d6), coeff)
        detdeg[p]=temp
        if n_elements(cp) gt 0 then detdeg[cp]=detdeg_fit[cp]
        yoffset=poly(sd2gps(0d)*1d6,coeff)-1d
        outdata[w] = ptr_new({wavelength:esrwaves, sd:pdsd, detdeg:detdeg-yoffset[0], coeff:[coeff[0]-yoffset,coeff[1:*]]})
   
        if not keyword_set(noplot) then begin
            title='Diode degradation ESR vs Diode '+strtrim(string(instrumentModeId),2)+' @ '+strtrim(string(esr_prismwaves[w],format='(F10.1)'),2)
            pdsd=gps2sd(pdsd/1d6)
            esrsd=gps2sd(esrsd/1d6)
            if keyword_set(yalign) then begin
                pdirrad -= median(pdirrad-esrirrad)
            endif
            correctedPDirrad /= detdeg_fit
           ; plot_multi,esrsd,esrirrad,pdsd,pdirrad,esrsd,correctedESRirrad,pdsd,correctedPDirrad,$
           ;     /xst,/yst, xtitle='Mission Day',ytitle='Irradiance',title=title,/bottom, $
           ;     color=[3,4,5,8], psym=[-4,-4,-4,-4], thick=[1.0,1.0,1.0,1.0],$
           ;     label=['ESR Uncorrected','Diode Uncorrected','ESR Corrected','Diode Corrected']
            plot_multi,pdsd,detdeg-yoffset[0],pdsd,detdeg_fit-yoffset[0], /xst,/yst, xtitle='Mission Day',ytitle='Diode Degradation',$
                title=title,/bottom, color=[9,10], psym=[-4,-3], thick=[1.0,2.0], yrange=[0.95,1.05], $
                label=['Diode Degradation','2nd order fit']
        endif

    endfor

    p=where(ptr_valid(outdata) eq 1,count)
    if count ne n_elements(outdata) then outdata=outdata[p]

    if keyword_set(surfacefit) then begin
        ;create the arrays from the coefficients obtained above
        wv=[]
        gps=[]
        degfit=[]
        degall=[]
        for i=0,n_elements(outdata)-1 do begin
            gps=[gps,(*outdata[i]).sd]
            wv=[wv,(*outdata[i]).wavelength]
            degfit=[degfit,poly((*outdata[i]).sd, (*outdata[i]).coeff)]
            degall=[degall,(*outdata[i]).detdeg]
        endfor
        degall = (degall<1.01d)>0.97d
        title='Diode Degradation Mode='+strtrim(string(instrumentModeId),2)
        ;srf=surface(degfit,gps2sd(gps/1d6),wv,irregular=grid,xtitle='Mission Day',ytitle='Wavelength',$
        ;   ztitle='Fractional Degradation',title=title)
        srf=surface(degall,gps2sd(gps/1d6),wv,irregular=grid,xtitle='Mission Day',ytitle='Wavelength',$
           ztitle='Fractional Degradation',title=title)

        ; get the polynomial fit
        data=dblarr(3,n_elements(gps))
        ; change the axis to mission day from microseconds since sfit doesn't handle very large numbers well
        ;
        ; WE'LL NEED TO MODIFY THE COEFFICIENTS TO WORK WITH MICROSECONDS IN JAVA CODE !!!
        ;
        data[0,*]=gps2sd(gps/1d6) & data[1,*]=wv & data[2,*]=degall
        res=sfit(data, 3, /irregular, kx=sfit_coeffs)

        ; do a 3-sigma cleanup of the data
        if n_elements(wavezero) gt 0 then begin
            ; force the degradation to nothing (1.0) for wavelengths smaller than provided wavezero
            ; and for mission days smaller than 200 days
            k=where(data[1,*] lt wavezero[0],count)
            if count gt 0 then data[2,k]=1d
        endif
        sigma3=stddev(degall-res) * 2d
        k=where(abs(degall-res) gt sigma3,count,complement=cp)
        if count gt 0 then res=sfit(data[*,cp], 3, /irregular, kx=sfit_coeffs)
        srf=surface(res,gps2sd(gps[cp]/1d6),wv[cp],irregular=grid,/overplot,color=[120,120,200],yrange=yrange,zrange=[0.98d,1.01d])
        
        if keyword_set(contourplot) then begin
            ct = COLORTABLE(72, /reverse)
            ; smoothed data
;            ctr = contour(res,gps2sd(gps[cp]/1d6),wv[cp], $
            ; raw data
            ctr = contour(degall,gps2sd(gps/1d6), wv, $
              /FILL, RGB_TABLE=ct, MARGIN=[0.15,0.3,0.1,0.1], $
              ; XRANGE=[min_temp1, max_temp1], YRANGE=[min_temp2, max_temp2], $
              ;dimensions=[8192,2978], $
              N_LEVELS=20, $
              TITLE='Degradation factor', $
              XTITLE='Mission day', $
              YTITLE='Wavelength (nm)')
            cb = COLORBAR(TITLE='Degradation factor')
        endif
    endif
    
    ; srf=surface(kappa_smooth.kappa,kappa_smooth.timesd,kappa_smooth.waves)

    return, outdata

end
