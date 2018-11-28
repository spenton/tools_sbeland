;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Using the fixed raypath from the ESR data, adjust the raypath for the  
;   requested photodiode so that the time series of the two detectors match
;   for the requested wavelengths.  The new rayPath value at each wavelength is 
;   ajusted so that the differences in the new "calibratedIrradiance" 
;   between ESR and diode is constant for the specified time range.
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
;    We use the SimUncorrectedIrradiance from version development which has no 
;    degradation but has 1AU correction applied.
;    Be careful to limit the wavelengths to avoid the diode degradation.
;    We assume the files with the time series data have been obtained
;    using the program compare_19_20.pro and saved in an IDL savefile.
;
;
; REVISION HISTORY:
;   Revision: $Id: adjust_pdesr_raypath.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;------------------------------------------------------------------
function func_getpdesr_slope, x, df, esrirrad=esrirrad, esrsd=esrsd, solexppd=solexppd, $
    pdirrad=pdirrad, new_pdirrad=new_pdirrad, delta_irrad=delta_irrad, kappa=kappa, $
    a_pd=a_pd, pdsd=pdsd, detectorDeg=detectorDeg

    ; calculate and minimize the slope of the differences between ESR and PD

    ;prismTransDegPD = (1d - a_pd) * exp(-kappa * solexppd) + a_pd * exp(-kappa * solexppd * 0.5d)
    prismTransDegPD = (1d - x[0]) * exp(-kappa * solexppd) + x[0] * exp(-kappa * solexppd / 2d)
    if n_elements(detectorDeg) eq 0 then detectorDeg=1d 
    new_pdirrad = pdirrad /  (prismTransDegPD * detectorDeg)
    delta_irrad = esrirrad - interpol(new_pdirrad, pdsd, esrsd)
    coeff = ladfit(esrsd, delta_irrad)
    ;print,'coeffs[1] = ',coeff[1]
    return, abs(coeff[1]*1d6)


    ;if n_elements(pdsd) gt 3 then begin
    ;    coeffs_pd = robust_poly_fit(pdsd, new_pdirrad, 2d, /double)
    ;endif else begin
    ;    coeffs_pd = ladfit(pdsd, new_pdirrad, /double)
    ;endelse
    ;if n_elements(esrsd) gt 3 then begin
    ;    coeffs_esr = robust_poly_fit(esrsd, esrirrad, 2d, /double)
    ;endif else begin
    ;    coeffs_esr = ladfit(esrsd, esrirrad, /double)
    ;endelse

    ;if n_elements(new_pdirrad) ge n_elements(esrirrad) then begin
    ;    fitA = poly(pdsd,coeffs_pd)
    ;    fitB = poly(pdsd,coeffs_esr)
    ;    delta_irrad = fitB - fitA
    ;    coeff = ladfit(pdsd, delta_irrad, /double)
    ;endif else begin
    ;    fitA = poly(esrsd,coeffs_pd)
    ;    fitB = poly(esrsd,coeffs_esr)
    ;    delta_irrad = fitB - fitA
    ;    coeff = ladfit(esrsd, delta_irrad, /double)
    ;endelse

    ;minimize the slope of the differences
    ;print,'coeffs = ',coeff
    ;return, abs(coeff[1]*1d14)

end
;*******************************************************************

function adjust_pdesr_raypath, starttime, stoptime, instrumentModeId, $
    gps=gps, missionDays=missionDays, julianDays=julianDays, $
    esrdata=esrdata, pddata=pddata, waverange=waverange, noplot=noplot, $
    solar54=solar54, solar55=solar55, yalign=yalign, pddeg=pddeg, alignobc=alignobc


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

    wavea = [    258.0d, 261.2,  265.2,  279.5,  280.4,  283.5,   288,    293.5,  304.5,  319,    $
                 333.5,  342,    355,    369,    375,    384.5,   395,    407,    431,    471,    $
                 480,    488,    519.5,  567,    594,    659,     710,    758,    806,    859,    $
                 867.5,  892.5,  901,    964.5, 1013,   1063.5,  1116,   1170,   1213.5, 1279.5, $
                1356.5, 1411,   1475.5, 1497,   1549.5, 1591,    1621.5, 1692,   1741,   1817,   $
                1882,   1909,   1944.5, 1988.5, 2014.2, 2106.5,  2171.2, 2287.5, 2397.5, 2502.5]

    prismposa = [4075L,   5050,   6250,  10000,  10232,  10916,  11950,  13075,  15100,  17425,  $
                 19375,  20350,  21775,  23125,  23650,  24400,  25150,  25975,  27400,  29275,  $
                 29650,  29950,  31000,  32275,  32875,  34075,  34825,  35425,  35950,  36475,  $
                 36550,  36775,  36850,  37375,  37750,  38125,  38500,  38875,  39175,  39625,  $
                 40150,  40525,  40975,  41125,  41500,  41800,  42025,  42550,  42925,  43525,  $
                 44050,  44275,  44575,  44950,  45175,  46000,  46600,  47725,  48850,  49975]

    waveb = [   278.40d, 279.37, 282.32, 284.17, 287.00, 292.43, 303.26, 317.69, 331.90, 339.90, $
                346.51,  352.87, 366.86, 372.82, 381.93, 391.83, 394.95, 403.72, 427.36, 466.35, $
                475.55,  483.34, 513.85, 559.75, 585.75, 649.07, 698.42, 744.69, 790.95, 843.00, $
                850.93,  875.37, 883.81, 945.85, 993.56, 1043.66, 1095.83, 1149.36, 1192.92, 1258.99, $
                1336.15, 1391.06, 1456.26, 1477.65, 1530.81, 1572.82, 1603.79, 1674.86, 1724.38, 1801.46, $
                1866.81, 1894.24, 1930.42, 1974.75, 2000.84, 2094.10, 2159.37, 2276.61, 2387.84, 2493.59]

    prismposb = [52388L,  52156,  51473,  51056,  50441,  49318,  47296,  44975,  43029,  42055,  $
                 41307L,  40633,  39286,  38762,  38014,  37265,  37041,  36442,  35020,  33149,  $
                 32775L,  32475,  31428,  30156,  29557,  28360,  27612,  27014,  26490,  25966,  $
                 25891L,  25667,  25592,  25069,  24695,  24321,  23946,  23572,  23273,  22824,  $
                 22301L,  21927,  21478,  21329,  20955,  20655,  20431,  19907,  19533,  18935,  $
                 18412L,  18188,  17888,  17514,  17290,  16467,  15869,  14748,  13626,  12504]

    if size(esrdata,/tname) ne 'STRUCT' then begin
        ; get the ESR data (version 1004 contains only the table scan data)
        if instrumentModeId lt 45 then begin
            ;restore,'~/SORCE/data/sima_esr_1004_uncorr.sav'
            ;restore,'~/SORCE/data/sima_esr_2001_uncorr.sav'
            restore,'~/SORCE/data/sima_esr_uncorr_2011.sav'
            p=where(ESRA_UNCORR_2011.spect20[*].timestamp[0] ge sd2gps(t0)*1d6 and ESRA_UNCORR_2011.spect20[*].timestamp[0] le sd2gps(t1)*1d6,count)
            if count eq 0 then begin
                print,'Error: no ESR data within requested time range'
                return,-1
            endif
            esrdata=ESRA_UNCORR_2011.spect20[p]
        endif else begin
            ;restore,'~/SORCE/data/simb_esr_2000_uncorr.sav'
            restore,'~/SORCE/data/simb_esr_uncorr_2011.sav'
            p=where(ESRB_UNCORR_2011.spect20[*].timestamp[0] ge sd2gps(t0)*1d6 and ESRB_UNCORR_2011.spect20[*].timestamp[0] le sd2gps(t1)*1d6,count)
            if count eq 0 then begin
                print,'Error: no ESR data within requested time range'
                return,-1
            endif
            esrdata=ESRB_UNCORR_2011.spect20[p]
        endelse
    endif

    ; get the latest ESR Kappa
    ;readcol,'~/SORCE/data/ESR_kappa_453_1570_1004.txt',esrkappa_w, esrkappa_k,format='(d,d)'
    ;readcol,'~/SORCE/data/ESRAB_453-1570_kappa_raypath_raytrace_1003.txt',esrkappa_w, esrkappa_k,format='(d,d)'
    ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2000.txt',esrkappa_w, esrkappa_k,format='(d,d)'
    ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2011.txt',esrkappa_w, esrkappa_k,format='(d,d)'
    ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2011_reduced55.txt',esrkappa_w, esrkappa_k,format='(d,d)'
    ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2011_mod_v20.txt',esrkappa_w, esrkappa_k,format='(d,d)'
    ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2011_mod_407nm.txt',esrkappa_w, esrkappa_k,format='(d,d)'

    ; get the latest ESR raypath
    ;readcol,'~/SORCE/data/overlap_esr_fit.txt',esrpath_w, tmp, esrpath_a, format='(d,d,d)'
    readcol,'~/SORCE/data/overlap_esr_v20.txt',esrpath_w, esrpath_a, format='(d,d)'

    if instrumentModeId eq 41 then begin
        ;readcol,'~/SORCE/data/overlap_vis1.txt',pd_aw,pd_apath,format='(d,d)'
        ;readcol,'~/SORCE/data/overlap_vis1_fit.txt',pd_aw,pd_apath,format='(d,d)'
        readcol,'~/SORCE/data/overlap_vis1_smooth.txt',pd_aw,pd_apath,format='(d,d)'
        if size(pddata,/tname) ne 'STRUCT' then begin
            ;restore,'~/SORCE/data/sima_vis_1003_uncorr.sav'
            ;restore,'~/SORCE/data/sima_vis_2000_uncorr.sav'
            restore,'~/SORCE/data/sima_vis_uncorr_2011.sav'
            p=where(VISA_UNCORR_2011.spect20[*].timestamp[0] ge sd2gps(t0)*1d6 and VISA_UNCORR_2011.spect20[0].timestamp[0] le sd2gps(t1)*1d6,count)
            if count eq 0 then begin
                print,'Error: no VIS data within requested time range'
                return,-1
            endif
            pddata=VISA_UNCORR_2011.spect20[p]
        endif
        ; for the VIS diode, we use the same Kappa as the ESR
        ;readcol,'~/SORCE/data/ESR_kappa_453_1570_1004.txt',pdkappa_w, pdkappa_k,format='(d,d)'
        ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2000.txt',pdkappa_w, pdkappa_k,format='(d,d)'
        ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2011_reduced55.txt',pdkappa_w, pdkappa_k,format='(d,d)'
        ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2011_mod_v20.txt',pdkappa_w, pdkappa_k,format='(d,d)'
        readcol,'~/SORCE/data/ESR_kappa_453_1570_2011_mod_407nm.txt',esrkappa_w, esrkappa_k,format='(d,d)'
        readcol,'~/SORCE/data/ESR_kappa_453_1570_2011_mod_407nm.txt',pdkappa_w, pdkappa_k,format='(d,d)'
        ; select a wavelength range to try to avoid diode degradation
        wave0=310d
        wave1=970d
        histbin=1.0
    endif else if instrumentModeId eq 43 then begin
        ;readcol,'~/SORCE/data/overlap_uv.txt',pd_aw,pd_apath,format='(d,d)'
        readcol,'~/SORCE/data/overlap_uv_smooth_pos.txt',pd_aw,pd_apath,format='(d,d)'
        if size(pddata,/tname) ne 'STRUCT' then begin
            ;restore,'~/SORCE/data/sima_uv_1003_uncorr.sav'
            restore,'~/SORCE/data/sima_uv_uncorr_2011.sav'
            p=where(UVA_UNCORR_2011.spect20[*].timestamp[0] ge sd2gps(t0)*1d6 and UVA_UNCORR_2011.spect20[*].timestamp[0] le sd2gps(t1)*1d6,count)
            if count eq 0 then begin
                print,'Error: no UV data within requested time range'
                return,-1
            endif
            pddata=UVA_UNCORR_2011.spect20[p]
        endif
        ; for the UV diode, we need a different Kappa since the ESR doesn't go blue enough
        ;readcol,'~/SORCE/data/UVAB_degradation_453-1570_kappa_ffunct_V70.txt',pdkappa_w, pdkappa_k,format='(d,d)'
        ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2011_mod_407nm.txt',pdkappa_w, pdkappa_k,format='(d,d)'
        readcol,'~/SORCE/data/UVAB_453-1570_kappa_2011.txt',pdkappa_w, pdkappa_k,format='(d,d)'
        readcol,'~/SORCE/data/UVAB_453-1570_kappa_2011.txt',esrkappa_w, esrkappa_k,format='(d,d)'
        ; select a wavelength range to try to avoid diode degradation
        wave0=210d
        wave1=308d
        histbin=0.5
    endif else if instrumentModeId eq 44 then begin
        readcol,'~/SORCE/data/overlap_ir.txt',pd_aw,pd_apath,format='(d,d)'
        if size(pddata,/tname) ne 'STRUCT' then begin
            restore,'~/SORCE/data/sima_ir_1003_uncorr.sav'
            p=where(ira_1003.spect20[*].timestamp[0] ge sd2gps(t0)*1d6 and ira_1003.spect20[*].timestamp[0] le sd2gps(t1)*1d6,count)
            if count eq 0 then begin
                print,'Error: no UV data within requested time range'
                return,-1
            endif
            pddata=ira_1003.spect20[p]
            delvar,ira_1003
        endif
        ; for the IR diode, we use the same Kappa as the ESR
        readcol,'~/SORCE/data/ESR_kappa_453_1570_1004.txt',pdkappa_w, pdkappa_k,format='(d,d)'
        ; select a wavelength range to try to avoid diode degradation
        wave0=1000d
        wave1=1400d
    endif else if instrumentModeId eq 45 then begin
        ;readcol,'~/SORCE/data/overlap_vis1.txt',pd_aw,pd_apath,format='(d,d)'
        ;readcol,'~/SORCE/data/overlap_vis1_fit.txt',pd_aw,pd_apath,format='(d,d)'
        readcol,'~/SORCE/data/overlap_vis1_smooth.txt',pd_aw,pd_apath,format='(d,d)'
        if size(pddata,/tname) ne 'STRUCT' then begin
            ;restore,'~/SORCE/data/simb_vis_1003_uncorr.sav'
            ;restore,'~/SORCE/data/simb_vis_2000_uncorr.sav'
            restore,'~/SORCE/data/simb_vis_uncorr_2011.sav'
            p=where(VISB_UNCORR_2011.spect20[*].timestamp[0] ge sd2gps(t0)*1d6 and VISB_UNCORR_2011.spect20[0].timestamp[0] le sd2gps(t1)*1d6,count)
            if count eq 0 then begin
                print,'Error: no VIS data within requested time range'
                return,-1
            endif
            pddata=VISB_UNCORR_2011.spect20[p]
        endif
        ; for the VIS diode, we use the same Kappa as the ESR
        ;readcol,'~/SORCE/data/ESR_kappa_453_1570_1004.txt',pdkappa_w, pdkappa_k,format='(d,d)'
        ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2000.txt',pdkappa_w, pdkappa_k,format='(d,d)'
        ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2011_mod_v20.txt',pdkappa_w, pdkappa_k,format='(d,d)'
        readcol,'~/SORCE/data/ESR_kappa_453_1570_2011_mod_407nm.txt',pdkappa_w, pdkappa_k,format='(d,d)'
        ; select a wavelength range to try to avoid diode degradation
        wave0=340d
        wave1=970d
    endif else if instrumentModeId eq 47 then begin
        ;readcol,'~/SORCE/data/overlap_uv.txt',pd_aw,pd_apath,format='(d,d)'
        readcol,'~/SORCE/data/overlap_uv_smooth_pos.txt',pd_aw,pd_apath,format='(d,d)'
        if size(pddata,/tname) ne 'STRUCT' then begin
            ;restore,'~/SORCE/data/simb_uv_1003_uncorr.sav'
            restore,'~/SORCE/data/simb_uv_uncorr_2011.sav'
            ;p=where(uvb_1003.spect20[*].timestamp[0] ge sd2gps(t0)*1d6 and uvb_1003.spect20[*].timestamp[0] le sd2gps(t1)*1d6,count)
            p=where(UVB_UNCORR_2011.spect20[*].timestamp[0] ge sd2gps(t0)*1d6 and UVB_UNCORR_2011.spect20[*].timestamp[0] le sd2gps(t1)*1d6,count)
            if count eq 0 then begin
                print,'Error: no UV data within requested time range'
                return,-1
            endif
            pddata=UVB_UNCORR_2011.spect20[p]
        endif
        ; for the UV diode, we need a different Kappa since the ESR doesn't go blue enough
        ;readcol,'~/SORCE/data/UV_kappa_453_1570_1003.txt',pdkappa_w, pdkappa_k,format='(d,d)'
        readcol,'~/SORCE/data/UVAB_453-1570_kappa_2011.txt',pdkappa_w, pdkappa_k,format='(d,d)'
        readcol,'~/SORCE/data/UVAB_453-1570_kappa_2011.txt',esrkappa_w, esrkappa_k,format='(d,d)'
        esrkappa_w=pdkappa_w
        esrkappa_k=pdkappa_k
        ; select a wavelength range to try to avoid diode degradation
        wave0=210d
        wave1=308d
    endif else if instrumentModeId eq 48 then begin
        readcol,'~/SORCE/data/overlap_ir.txt',pd_aw,pd_apath,format='(d,d)'
        if size(pddata,/tname) ne 'STRUCT' then begin
            restore,'~/SORCE/data/simb_ir_1003_uncorr.sav'
            p=where(irb_1003.spect20[*].timestamp[0] ge sd2gps(t0)*1d6 and irb_1003.spect20[*].timestamp[0] le sd2gps(t1)*1d6,count)
            if count eq 0 then begin
                print,'Error: no UV data within requested time range'
                return,-1
            endif
            pddata=irb_1003.spect20[p]
            delvar,irb_1003
        endif
        ; for the IR diode, we use the same Kappa as the ESR
        readcol,'~/SORCE/data/ESR_kappa_453_1570_1004.txt',pdkappa_w, pdkappa_k,format='(d,d)'
        ; select a wavelength range to try to avoid diode degradation
        wave0=1000d
        wave1=1400d
    endif


    ; use the modified solar exposure record  SBeland 2014-06-05
    ;if n_elements(solar54) eq 0 then restore,'~/SORCE/data/solarexp_54_55_modified_v20.sav'
    if n_elements(solar54) eq 0 then restore,'~/SORCE/data/solarexp_54_55_mod_407nm.sav'


    ;if n_elements(solar54) eq 0 and instrumentModeId lt 45 then begin
    ;    if keyword_set(verbose) then print,'   getting solar54 ...'
    ;    q1 = 'SELECT orbitStartTimeGps t0, orbitEndTimeGps t1, cumulativeSolarExpHrtOut solar_exp, '+$
    ;         'solarExposureHrtOut as solar_exp_orbit from SimSolarExposureData where instrumentModeId=54'
    ;    query_database, q1, solar54, info
    ;endif else if n_elements(solar55) eq 0 and instrumentModeId ge 45 then begin
    ;    if keyword_set(verbose) then print,'   getting solar55 ...'
    ;    q1 = 'SELECT orbitStartTimeGps t0, orbitEndTimeGps t1, cumulativeSolarExpHrtOut solar_exp, '+$
    ;         'solarExposureHrtOut as solar_exp_orbit from SimSolarExposureData where instrumentModeId=55'
    ;    query_database, q1, solar55, info
    ;endif

    ; starting parameters for minimizing function
    parinfo=replicate({value:0.0d, fixed:0, limited:[0,0], limits:[0d,0d], step:1d-5, tnside:2},1)

    if instrumentModeId ge 45 then begin
        wavea=waveb
        prismposa=prismposb
    endif

    if n_elements(waverange) eq 2 then wave0=min(waverange,max=wave1)
    p=where(wavea ge wave0 and wavea le wave1,count)
    wavea=wavea[p]
    prismposa=prismposa[p]
    pd_a_factor=dblarr(n_elements(wavea))
    detdeg=dblarr(n_elements(wavea))

    if keyword_set(pddeg) and (instrumentModeId eq 41 or instrumentModeId eq 45) then begin
        ; for now we only have data for mode 41
        ;query='SELECT * from SimDiodeDegCalPolynomialCoeffs where calibrationSetId=1159'
        ;query_database,query,data
        ;db_pddeg=reform(data[*].coeffvalue)
        ;restore,'/Users/sbeland/SORCE/data/simab_vis_deg_coeffs_afact48.sav'
        restore,'/Users/sbeland/SORCE/data/simab_vis_deg_coeffs_mod_v20.sav'
    endif

    ; will plot only every skip plot since lineplot is limited to 24 plots
    skip = round(float(n_elements(wavea)) / 12.0)

    ; loop over every ESR table scan wavelength within our range
    for w=0L,n_elements(wavea)-1 do begin
        ; get the list of ESR points
        print,'processing ',wavea[w]
        ;p=where(abs(esrdata.wavelength - wavea[w]) lt histbin and esrdata.timestamp ge sd2gps(t0)*1d6 and esrdata.timestamp le sd2gps(t1)*1d6,count)
        p=where(esrdata.prismposition eq prismposa[w] and esrdata.timestamp ge sd2gps(t0)*1d6 and esrdata.timestamp le sd2gps(t1)*1d6,count)
        if count eq 0 then begin
            print,'No ESR data for ',wavea[w],prismposa[w]
            continue
        endif

        esrwaves=(esrdata.wavelength)[p]
        esrirrad=(esrdata.irradiance)[p]
        esrsd=(esrdata.timestamp)[p]

        ; remove spikes in ESR data
        ;if n_elements(esrsd) lt 4 then continue
        ;coeff=robust_poly_fit(esrsd, esrirrad, 3, /double)
        ;resistant_mean,(esrirrad-poly(esrsd,coeff)), 4.0, mean, good=keep
        ;if n_elements(keep) ne n_elements(esrsd) then begin
        ;    print,'removed '+strtrim(string(n_elements(esrsd)-n_elements(keep)),2)+' ESR points ...'
        ;    p=where(abs(esrwaves[keep] - wavea[w]) lt 2.0d,count)
        ;    esrwaves=esrwaves[keep[p]]
        ;    esrirrad=esrirrad[keep[p]]
        ;    esrsd=esrsd[keep[p]]
        ;endif

        kappaesr = interpol(esrkappa_k, esrkappa_w, esrwaves, /lsq)
        kappapd = interpol(pdkappa_k, pdkappa_w, esrwaves, /lsq)
        a_esr = interpol(esrpath_a, esrpath_w, esrwaves, /lsq)
        a_pd = interpol(pd_apath,pd_aw, esrwaves, /lsq)
        if instrumentModeId lt 45 then $
            solexpesr=interpol(solar54.solar_exp,solar54.t1,esrsd) $
        else $
            solexpesr=interpol(solar55.solar_exp,solar55.t1,esrsd)
        solexpesr /= 86400d
        ; scale the raypath by the oneau
        ;if where(strpos(tag_names(solar55), 'ONEAU') eq 0) ge 0 then begin
        ;    a_esr = a_esr[0] * (interpol(solar54.oneau, solar54.t1, esrsd))^2d
        ;endif
        correctedESRirrad = esrirrad / ( (1d - a_esr) * exp(-kappaesr * solexpesr) + a_esr * exp(-kappaesr * solexpesr / 2d) )

        pdirrad=dblarr(n_elements(esrwaves))
        pdsd=pdirrad
        new_pdirrad=pdirrad
        delta_irrad=pdirrad
        mnw=min(esrwaves,max=mxw)

        ; get diode data for the same wavelength and timestamp
        for t=0L,count-1 do begin
            mn=min(abs(pddata.timestamp - esrsd[t]),pp)
            pos=array_indices(pddata.timestamp,pp)
            ; interpolate the pd spectra to get the uncorrectedIrradiance at the same wavelength as ESR
            p = where(pddata[pos[1]].wavelength gt 0d)
            mn=min(pddata[pos[1]].wavelength[p], max=mx)
            if (mn gt mxw or mx lt mnw) then continue
            s=sort((pddata[pos[1]].wavelength)[p])
            pdirrad[t] = interpol((pddata[pos[1]].irradiance)[p[s]], (pddata[pos[1]].wavelength)[p[s]], esrwaves[t],/spline)
            mn=min(abs((pddata[pos[1]].timestamp)[p] - esrsd[t]),post)
            pdsd[t]=(pddata[pos[1]].timestamp)[post]
        endfor

        k=where(pdsd gt 0d,count)
        if count lt 2 then begin
            print,'no diode data for ',wavea[w]
            continue
        endif else if count ne n_elements(pdsd) then begin
            correctedESRirrad=correctedESRirrad[k]
            kappapd=kappapd[k]
            pdirrad=pdirrad[k]
            pdsd=pdsd[k]
            esrsd=esrsd[k]
            esrirrad=esrirrad[k]
            a_pd=a_pd[k]
        endif

        ; now find the diode Kappa at his wavelength that will match the ESR trend
        if instrumentModeId lt 45 then $
            solexppd=interpol(solar54.solar_exp,solar54.t1,pdsd) $
        else $
            solexppd=interpol(solar55.solar_exp,solar55.t1,pdsd)
        solexppd /= 86400d
        parinfo[0].value=median(a_pd)
        ;parinfo[0].fixed=1
        ;parinfo[1].value=1.0d
        ;parinfo[1].fixed=1

        if keyword_set(pddeg) then begin
            detectorDeg = dblarr(n_elements(esrwaves))
            for i=0,n_elements(esrwaves)-1 do detectorDeg[i] = (sfit_poly(gps2sd(pdsd[i]/1d6), esrwaves[i], visab_coeffs)) < 1.0d
        endif else detectorDeg=1d

        ; scale the raypath by the oneau
        ;if where(strpos(tag_names(solar55), 'ONEAU') eq 0) ge 0 then begin
        ;    a_pd = a_pd[0] * (interpol(solar54.oneau, solar54.t1, pdsd))^2d
        ;endif

        if keyword_set(alignobc) then begin
            align_irrad, esrsd, correctedESRirrad, /mission
            align_irrad, pdsd, pdirrad, /mission
        endif


        functargs = {esrirrad:correctedESRirrad, esrsd:gps2sd(esrsd/1d6), solexppd:solexppd, kappa:kappapd, $
            pdirrad:pdirrad, new_pdirrad:new_pdirrad, delta_irrad:delta_irrad, a_pd:a_pd, $
            pdsd:gps2sd(pdsd/1d6), detectorDeg:detectorDeg}
        coeffs = tnmin('func_getpdesr_slope', functargs=functargs, bestmin=f0, status=status, $
            nfev=nfev, niter=niter, autoderivative=1, errmsg=errmsg, parinfo=parinfo)

        pd_a_factor[w]=coeffs[0]
        print,'detdeg[-1]=',detectorDeg[-1]

        ; get the delta_irrad with the best fit pd raypath
        res=func_getpdesr_slope(pd_a_factor[w], df, esrirrad=correctedESRirrad, esrsd=esrsd, solexppd=solexppd, $
            pdirrad=pdirrad, new_pdirrad=new_pdirrad, delta_irrad=delta_irrad, kappa=kappapd, $
            a_pd=a_pd, pdsd=pdsd, detectorDeg=detectorDeg)

        if not keyword_set(noplot) then begin
            title='ESR vs Diode'+' @ '+strtrim(string(wavea[w],format='(F10.1)'),2)+' raypath='+strtrim(string(pd_a_factor[w],format='(F0.4)'),2)
            title=title+' diodeDeg='+strtrim(string(detectorDeg[-1],format='(F0.4)'),2)
            pdsd=gps2sd(pdsd/1d6)
            esrsd=gps2sd(esrsd/1d6)
            if keyword_set(yalign) then new_pdirrad -= median(new_pdirrad-correctedESRirrad)
            if keyword_set(yalign) then pdirrad -= median(pdirrad-esrirrad)
            if keyword_set(alignobc) then begin
                align_irrad, esrsd, esrirrad, /mission
                align_irrad, esrsd, correctedESRirrad, /mission
                align_irrad, pdsd, pdirrad, /mission
                align_irrad, pdsd, new_pdirrad, /mission
            endif
            plot_multi,esrsd,esrirrad,pdsd,pdirrad,esrsd,correctedESRirrad,pdsd,new_pdirrad,/xst,/yst,xtitle='Mission Day',ytitle='Irradiance', $
                label=['ESR Uncorrected','Diode Uncorrected','ESR Corrected','Diode Corrected'],title=title,/bottom
            if (w mod skip) eq 0 or w eq n_elements(wavea)-1 then begin
                lineplot, esrsd,correctedESRirrad, title='ESR Corrected  '+strtrim(string(wavea[w],format='(F10.1)'),2), psym=-4
                lineplot, pdsd,new_pdirrad, title='Diode Corrected  '+strtrim(string(wavea[w],format='(F10.1)'),2),psym=-4
            endif
        endif

    endfor

    p=where(pd_a_factor ne 0.0d,count)
    if count eq 0 then return,-1

    wavea=wavea[p]
    ;pd_a_factor=abs(pd_a_factor[p])
    pd_a_factor=pd_a_factor[p]
    detdeg=detdeg[p]

    if not keyword_set(noplot) then begin
        if instrumentModeId eq 41 then inst='VISA'
        if instrumentModeId eq 43 then inst='UVA'
        if instrumentModeId eq 44 then inst='IRA'
        if instrumentModeId eq 45 then inst='VISB'
        if instrumentModeId eq 47 then inst='UVB'
        if instrumentModeId eq 48 then inst='IRB'

        ;k=[1,2,3,7,8,9,10,12,13,14,15,17,18,19,20]
        k=lindgen(n_elements(wavea))
        ;cc=robust_poly_fit(wavea[k],pd_a_factor[k],2,yfit)
        if n_elements(wavea) lt 6 then cc=[median(pd_a_factor),0d] else $
            cc=robust_poly_fit(wavea,pd_a_factor,4,yfit)
        if instrumentModeId eq 41 or instrumentModeId eq 45 then begin
            p=where(wavea[k] lt 500d)
            coef=poly_fit(wavea[k[p]], pd_a_factor[k[p]],1)
            delta = poly(420d,coef) - interpol(pd_apath,pd_aw,420d)
            plot_multi,pd_aw,pd_apath,wavea,pd_a_factor,wavea[k],pd_a_factor[k],wavea,poly(wavea,cc),pd_aw,pd_apath+delta,/xst,/yst,$
                xtitle='Wavelength (nm)',ytitle='Raypath',title=inst+' Raypath obtained from ESR geometric raypath ('+$
                strtrim(string(t0,format='(I)'),2)+' to '+strtrim(string(t1,format='(I)'),2)+')',$
                label=[inst+' Geometric Raypath',inst+' Optimized Raypath',inst+' Optimized Raypath CLEAN','4th ordr Fit',$
                'Geo Raypath +'+strtrim(string(delta,format='(F8.6)'),2)],/bottom,$
                psym=[-3,-4,-4,-3,-3],thick=[4.0,1.0,1.0,3.0,4.0],color=[2,3,4,5,8]
        endif else begin
            plot_multi,pd_aw,pd_apath,wavea,pd_a_factor,wavea[k],pd_a_factor[k],wavea,poly(wavea,cc),/xst,/yst,$
                xtitle='Wavelength (nm)',ytitle='Raypath',title=inst+' Raypath obtained from ESR geometric raypath ('+$
                strtrim(string(t0,format='(I)'),2)+' to '+strtrim(string(t1,format='(I)'),2)+')',$
                label=[inst+' Geometric Raypath',inst+' Optimized Raypath',inst+' Optimized Raypath CLEAN','4th ordr Fit'],/bottom,$
                psym=[-4,-4,-4,-3],thick=[4.0,1.0,1.0,3.0]
        endelse

        ;plot_multi,wavea,detdeg,/xst,/yst,$
        ;    xtitle='Wavelength (nm)',ytitle='Diode Degradation',title=inst+' diode degradation from ESR '
    endif

    return,{wavelength:wavea, a_factor:pd_a_factor, detdeg:detdeg}

end
