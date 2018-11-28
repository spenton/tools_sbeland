;+
; NAME:   GET_SINGLE_KAPPA
;
; PURPOSE: 
;    This routine compares one spectrum of uncorrected irradiance from SIMB 
;    and a spectrum from SIMA taken at about the same time and determines 
;    the value of Kappa using the solar exposure and the raypath values for 
;    each channel.
;
; CALLING SEQUENCE:
;    result=get_single_kappa(specta, spectb, afact=afact, $
;           solar54=solar54, solar55=solar55, waves=waves, t0=t0, t1=t1)
;
; RETURNED VALUE:
;    result - data structure containing the value of Kappa for each wavelength of the SIMB spectrum
; INPUT PARAMETERS:
;    specta - data structure containing the spectra from SIMA (saved from compare_19_20.pro)
;    spectb - data structure containing the spectra from SIMB (saved from compare_19_20.pro)
;
; OPTIONAL INPUT PARAMETERS:
;    afact - data structure containing the raypath for the specific detector (same for Sima and SIMB)
;    solar54 - data structure containing the solar exposure records for SIMA
;    solar55 - data structure containing the solar exposure records for SIMB
;    waves - array of wavelengths for which to calculate a value of Kappa
;    t0 - value specifying the start time to extract Kappa
;    t1 - value specifying the stop time to extract Kappa
;
; OPTIONAL INPUT KEYWORDS:
;    missionDays - indicates trange was specified in SORCE Mission Days (default)
;    gps - indicates trange was specified in GPS
;    julianDays - indicates trange was specified in Julian Days
;
; OUTPUT PARAMETERS:
;    none
;
; OPTIONAL OUTPUT PARAMETERS:
;    NONE
;
; EXAMPLE:  
;
; COMMON BLOCKS:
;    NONE
;
; NOTES:
;
; REVISION HISTORY:
;   Revision: $Id: get_single_kappa.pro,v 1.2 2018/11/26 16:16:02 spenton Exp $
;-
;*******************************************************************
function get_single_kappa_amoeba, x
    common single_kappa_common, subSolA, subIrdA, subSolB, subIrdB, afact_val, afactb_val
    ; measured_irrad = SolarIrrad * exp(-abs(x[0]*solexp))
    ; SolarIrrad = measured_irrad / exp(-abs(x[0]*solexp))
    prismDeg_A = (1d - afact_val[0])*exp(-abs(x*subSolA[0])) + afact_val[0]*exp(-abs(x*subSolA[0])/2d)
    prismDeg_B = (1d - afactb_val[0])*exp(-abs(x*subSolB[0])) + afactb_val[0]*exp(-abs(x*subSolB[0])/2d)

;  TESTING A DIFFERENT KAPPA FUNCTION
    ;prismDeg_A = exp(-abs((2d - afact_val[0])* x *subSolA[0]))
    ;prismDeg_B = exp(-abs((2d - afact_val[0])* x *subSolB[0]))

    out_value = ABS(subIrdA/prismDeg_A - subIrdB/prismDeg_B)

    ;TODO add the diode temperature correction 

    ;out_value = subIrdB/subIrdA - prismDeg_B/prismDeg_A
    ;print,x,abs(out_value)*1d6
    return, out_value*1d12
end
;*****************************************************************************
FUNCTION GET_SINGLE_KAPPA, specta, spectb, instMode, solar54=solar54, solar55=solar55, $
    afacta=afact, afactb=afactb, t0=t0, t1=t1, waves=waves, all_kappa=all_kappa, $
    gps=gps, julianDays=julianDays, missionDays=missionDays, $
    kappa_smooth=kappa_smooth, regular=regular, quiet=quiet, bspline=bspline

    common single_kappa_common

    ; get the SimSolarExposureData for modes 54 and 55
    if n_elements(solar54) eq 0 then begin
        q1 = 'SELECT orbitStartTimeGps t0, orbitEndTimeGps t1, cumulativeSolarExpHrtOut solar_exp, '+$
             'solarExposureHrtOut as solar_exp_orbit from SimSolarExposureData where instrumentModeId=54'
        query_database, q1, solar54, info
        solar54.solar_exp = total(solar54.solar_exp_orbit, /cum)
    endif

    if n_elements(solar55) eq 0 then begin
        q1 = 'SELECT orbitStartTimeGps t0, orbitEndTimeGps t1, cumulativeSolarExpHrtOut solar_exp, '+$
             'solarExposureHrtOut as solar_exp_orbit from SimSolarExposureData where instrumentModeId=55'
        query_database, q1, solar55, info
        solar55.solar_exp = total(solar55.solar_exp_orbit, /cum)
    endif

    get_oneau=0
    tnames = tag_names(solar54)
    p=where(strpos(tnames,"ONEAU") ge 0, count54)
    if count54 eq 0 then get_oneau=1
    tnames = tag_names(solar55)
    p=where(strpos(tnames,"ONEAU") ge 0, count55)
    if count55 eq 0 then get_oneau=1
    
    if get_oneau then begin
        ; get the 1-AU correction
        print,'  getting SUNOBSERVERDISTANCECORRECTION from database ...'
        q1="SELECT microsecondsSinceGpsEpoch, sunObserverDistanceCorrection FROM SolarDistAndDopplerFixedStep2 "+$
           "WHERE instrumentModeId=50 and version=7"
        query_database, q1, solardist, info
        if count54 eq 0 then begin
            oneau54 = interpol(solardist.SUNOBSERVERDISTANCECORRECTION, solardist.MICROSECONDSSINCEGPSEPOCH, solar54.(1))
            append_tag,solar54,'oneau',oneau54,/slim
        endif
        if count55 eq 0 then begin
            oneau55 = interpol(solardist.SUNOBSERVERDISTANCECORRECTION, solardist.MICROSECONDSSINCEGPSEPOCH, solar55.(1))
            append_tag,solar55,'oneau',oneau55,/slim
        endif
    endif

    query_database,/reset
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

    if n_elements(afactb) eq 0 then afactb=afact

    if n_elements(t0) gt 0 and n_elements(t1) gt 0 then begin
        ; t0 and t1 were specified -> get the list of SIMB spectra within this time
        if keyword_set(julianDays) then begin
            st0=jd2sd(t0)
            st1=jd2sd(t1)
        endif else if keyword_set(gps) then begin
            st0=gps2sd(t0/1d6)
            st1=gps2sd(t1/1d6)
        endif else begin
            st0=t0
            st1=t1
        endelse
        if PTR_VALID(spectb.spect20[0]) then begin
            ; we have an array of pointers from the new compare_19_20
            posb=[]
            for i=0L,n_elements(spectb.spect20)-1 do begin
                if (*spectb.spect20[i]).timestamp[0] ge sd2gps(st0)*1d6 and $
                   (*spectb.spect20[i]).timestamp[0] le sd2gps(st1)*1d6  then posb=[posb,i]
            endfor
            nspect=n_elements(posb)
        endif else begin
            posb = where(spectb.spect20[*].timestamp[5] ge sd2gps(st0)*1d6 and $
                         spectb.spect20[*].timestamp[5] le sd2gps(st1)*1d6,nspect)
        endelse
        if nspect eq 0 then begin
            print,'Error: no SIMB spectra within specified timerange'
            return,-1
        endif
    endif else if n_elements(t0) gt 0 then begin
        ; only t0 was specified 
        if keyword_set(julianDays) then begin
            st0=jd2sd(t0)
        endif else if keyword_set(gps) then begin
            st0=gps2sd(t0/1d6)
        endif else begin
            st0=t0
        endelse
        if PTR_VALID(spectb.spect20[0]) then begin
            ; we have an array of pointers from the new compare_19_20
            st1=0d
            posb=[]
            for i=0L,n_elements(spectb.spect20)-1 do begin
                if (*spectb.spect20[i]).timestamp[0] ge sd2gps(st0)*1d6 then posb=[posb,i]
                st1 = st1 > max((*spectb.spect20[i]).timestamp)
            endfor
            st1=gps2sd(st1/1d6)
            nspect=n_elements(posb)
        endif else begin
            st1 = gps2sd(max(spectb.spect20[*].timestamp)/1d6)
            posb = where(spectb.spect20[*].timestamp[5] ge sd2gps(st0)*1d6,nspect)
        endelse
        if nspect eq 0 then begin
            print,'Error: no SIMB spectra within specified timerange'
            return,-1
        endif
    endif else if n_elements(t1) gt 0 then begin
        ; only t1 was specified 
        if keyword_set(julianDays) then begin
            st1=jd2sd(t1)
        endif else if keyword_set(gps) then begin
            st0=gps2sd(t0/1d6)
        endif else begin
            st1=t1
        endelse
        if PTR_VALID(spectb.spect20[0]) then begin
            ; we have an array of pointers from the new compare_19_20
            st0 = 1d99
            posb=[]
            for i=0L,n_elements(spectb.spect20)-1 do begin
                if (*spectb.spect20[i]).timestamp[0] le sd2gps(st1)*1d6 then posb=[posb,i]
                st0 = st0 < min((*spectb.spect20[i]).timestamp)
            endfor
            st0=gps2sd(st0/1d6)
            nspect=n_elements(posb)
        endif else begin
            k=where(spectb.spect20.timestamp gt 1d6)
            st0 = gps2sd(min((spectb.spect20.timestamp)[k])/1d6)
            posb = where(spectb.spect20[*].timestamp[5] le sd2gps(st1)*1d6,nspect)
        endelse
        if nspect eq 0 then begin
            print,'Error: no SIMB spectra within specified timerange'
            return,-1
        endif
    endif else begin
        nspect = n_elements(spectb.spect20)
        posb = lindgen(nspect)
    endelse

    if n_elements(waves) eq 0 then begin
        q1='select stdWavelength as wavelength from SimStandardWavelengths where version=3'
        q1=q1+' and instrumentModeId='+strtrim(string(instMode),2)
        query_database,q1,refw,nrows
        waves=refw.wavelength
        if instMode eq 43 then waves=waves[where(waves le 308.5)]
    endif

    ; loop over every SIMB spectra
    wid=10
    window,wid,xsize=1200,ysize=600

    if solar54[-1].solar_exp gt 1d4 then day_scale=86400d else day_scale=1d
    ; if the kappas were passed in, do not re-calculate them
    if n_elements(all_kappa) eq 0 then begin
        ; calculate a Kappa vs wavelength for each day
        timeA0=dblarr(n_elements(specta.spect20))
        timeB0=dblarr(nspect)
        if ptr_valid(spectb.spect20[0]) then begin
            for i=0L,nspect-1L do timeB0[i]=gps2sd(max((*spectb.spect20[posb[i]]).timestamp)/1d6)
            for i=0L,n_elements(specta.spect20)-1L do timeA0[i]=gps2sd(max((*specta.spect20[i]).timestamp)/1d6)
        endif else begin
            for i=0L,nspect-1L do timeB0[i]=gps2sd(max(spectb.spect20[posb[i]].timestamp)/1d6)
            timeA0=gps2sd(spectA.spect20[*].timestamp[10]/1d6)
        endelse
        all_kappa = []
        if ~keyword_set(quiet) then print,'Processing ',nspect,' spectra ...'
        for sp=0L,nspect-1 do begin
            ; extract the valid values for this SIMB spectra
            ;if ~keyword_set(quiet) then print,sp+1,format='("   ",i0,$,%"\r")'
            print,sp+1,timeB0[sp],format='("   ",i0,"   ",F8.2)'
            if ptr_valid(spectb.spect20[0]) then begin
                k=where((*spectb.spect20[posb[sp]]).wavelength gt 0d,count)
                if count lt 5 then continue
                iradB=(*spectb.spect20[posb[sp]]).irradiance[k]
                waveB=(*spectb.spect20[posb[sp]]).wavelength[k]
                timeB=(*spectb.spect20[posb[sp]]).timestamp[k]
            endif else begin
                k=where(spectb.spect20[posb[sp]].wavelength gt 0d,count)
                if count lt 5 then continue
                iradB=spectb.spect20[posb[sp]].irradiance[k]
                waveB=spectb.spect20[posb[sp]].wavelength[k]
                timeB=spectb.spect20[posb[sp]].timestamp[k]
            endelse
            s=sort(waveB)
            iradB=iradB[s]
            waveB=waveB[s]
            timeB=timeB[s]
            ; find the SIMA spectrum closest in time
            if ptr_valid(spectb.spect20[0]) then begin
                mn=min(abs(timeA0 - timeB0[sp]),pos)
                k=where((*spectA.spect20[pos]).wavelength gt 0d,count)
                iradA=(*spectA.spect20[pos]).irradiance[k]
                waveA=(*spectA.spect20[pos]).wavelength[k]
                timeA=(*spectA.spect20[pos]).timestamp[k]
            endif else begin
                mn=min(abs(specta.spect20[*].timestamp[10] - timeB[sp]),pos)
                k=where(spectA.spect20[pos].wavelength gt 0d,count)
                iradA=spectA.spect20[pos].irradiance[k]
                waveA=spectA.spect20[pos].wavelength[k]
                timeA=spectA.spect20[pos].timestamp[k]
            endelse
            ; TODO ignore data separated by more than one day
            s=sort(waveA)
            iradA=iradA[s]
            waveA=waveA[s]
            timeA=timeA[s]
            ; get range of overlapping wavelengths
            minw=max([min(waveB),min(waveA)])
            maxw=min([max(waveB),max(waveA)])

            ;check if a list of wavelengths were specified
            if n_elements(waves) gt 0 then begin
                if min(waveB) gt min(waves) or max(waveB) lt max(waves) then begin
                    iradB = interpol(iradB, waveB, waves, /lsq) 
                    timeB = interpol(timeB, waveB, waves)
                endif else begin
                    iradB = interpol(iradB, waveB, waves, /spl)
                    timeB = interpol(timeB, waveB, waves)
                endelse
                waveB = waves
            endif
            ; use the same wavelength grid for SIMA and SIMB
            ; check if we need to extrapolate
            if min(waveB) lt min(waveA) or max(waveB) gt max(waveA) then $
                iradA = interpol(iradA, waveA, waveB, /lsq) $
            else $
                iradA = interpol(iradA, waveA, waveB, /spl)
            waveA = waveB
            kappa=dblarr(n_elements(waveB)) +1d6
            ; get the solar exposure for this orbit
            pa=where(solar54.t0 lt timeB[0])
            subSolA = solar54[pa[-1]].solar_exp / day_scale
            pb=where(solar55.t0 lt timeB[0])
            subSolB = solar55[pb[-1]].solar_exp / day_scale

            afact_wave = interpol(afact.afact, afact.wavelength, waveB, /spl)
            afact_values= afact_wave*(interpol(solar54.oneau,solar54.t1,timeB[0]))^2d
            afactb_wave = interpol(afactb.afact, afactb.wavelength, waveB, /spl)
            afactb_values= afactb_wave*(interpol(solar54.oneau,solar54.t1,timeB[0]))^2d

            ; fit one wavelength at a time
            for w=0L, n_elements(waveB)-1 do begin
                subIrdA = iradA[w]
                subIrdB = iradB[w]
                afact_val=afact_values[w]
                afactb_val=afactb_values[w]
                cc = AMOEBA(1d-20, P0=1d, scale=1d, FUNCTION_NAME='get_single_kappa_amoeba',ncalls=ncalls)
                if cc[0] ne -1 then begin
                    kappa[w] = abs(cc[0])
                    ;print,waveB[w], kappa[w]
                endif else if ~keyword_set(quiet) then print,'Failed to converge for ',sp,waveB[w]
            endfor
            answer=''
            if count gt 0 then begin
                ;if gps2sd(timeb[sp]/1d6) lt 4050 and ~keyword_set(bspline) then begin
                if 1 and instMode eq 43 then begin

                    maxw=306.0 
                    nnodes=20
                    nord=4
                    k=where(kappa lt 1d6, count)
                    p = where(waveB[k] ge minw and waveB[k] le maxw,count2)
                    if count2 eq 0 then p=lindgen(n_elements(k))
                    ; ignore data outside of 5 sigma from a 3rd order polynomial fit
                    coeff=robust_poly_fit(waveB[k[p]], kappa[k[p]], 6d)
                    resistant_mean, (kappa[k[p]] - poly(waveB[k[p]],coeff)), 5.0,mean,goodvec=keep0
                    sset=bspline_iterfit(waveB[k[p[keep0]]], kappa[k[p[keep0]]], requiren=nnodes,bkspace=nnodes, nord=nord)
                    yfit=bspline_valu(waveB, sset)
                    if instMode eq 43 then begin
                        aa=where(waveb gt 200d and waveb lt 300d)
                        cc=linfit(waveB[aa],yfit[aa])
                        if cc[1] gt 0d or max(yfit) lt 1d-4 then continue
                    endif

                    ; for data in DOOP, use the broad spline since we don't seem to have a 
                    ; perfect wavelength alignment which seems to introduce a lot of
                    ; variations in Kappa vs wavelength in one spectra
                    ;g=where(waveb[k[p[keep0]]] le 225,comp=cg)
                    ;sset=bspline_iterfit(waveB[k[p[keep0[g]]]], kappa[k[p[keep0[g]]]], requiren=5,bkspace=5, nord=3)
                    ;yfit1=bspline_valu(waveB, sset)
                    ;sset=bspline_iterfit(waveB[k[p[keep0[cg]]]], kappa[k[p[keep0[cg]]]], requiren=1,bkspace=1, nord=9)
                    ;yfit2=bspline_valu(waveB, sset)
                    ;p=where(waveb le 225, comp=cp)
                    ;yfit[p]=yfit1[p]
                    ;yfit[cp]=yfit2[cp]
                    ; test a smoother kappa for UV 
                    sset=bspline_iterfit(waveB[k[p[keep0]]], kappa[k[p[keep0]]], requiren=5,bkspace=5, nord=3)
                    yfit=bspline_valu(waveB, sset)
                    ; use a finer bspline for wavelengths greater than 240.0nm to capture bumps at 280 and 285
                    kk=where(waveB[k[p[keep0]]] ge 240.0d)
                    jj=where(waveB lt 240.0d)
                    sset=bspline_iterfit([waveB[jj],waveB[k[p[keep0[kk]]]]], [yfit[jj],kappa[k[p[keep0[kk]]]]], requiren=1,bkspace=1, nord=3)
                    yfit=bspline_valu(waveB, sset)

                    ;q=where(waveb ge 220.0,count)
                    ;if count gt 0 then yfit[q]=yfit2[q]
                    yrange=[0.001,0.008]
                endif else if instMode eq 41 then begin

                    minw=324.0
                    maxw=957.0
                    nnodes=30
                    nord=3
                    k=where(kappa lt 1d6, count)
                    p = where(waveB[k] ge minw and waveB[k] le maxw,count2)
                    if count2 eq 0 then p=lindgen(n_elements(k))
                    ; ignore data outside of 5 sigma from a 3rd order polynomial fit
                    coeff=robust_poly_fit(waveB[k[p]], kappa[k[p]], 6d)
                    resistant_mean, (kappa[k[p]] - poly(waveB[k[p]],coeff)), 5.0,mean,goodvec=keep0
                    sset=bspline_iterfit(waveB[k[p[keep0]]], kappa[k[p[keep0]]], requiren=nnodes,bkspace=nnodes, nord=nord)
                    yfit=bspline_valu(waveB, sset) > 0d

                    ; get a finer bspline for wavelengths above 410nm
                    ;sset2=bspline_iterfit(waveB[k[p[keep0]]], kappa[k[p[keep0]]], requiren=10,bkspace=10, nord=nord)
                    ;p=where(waveB gt 405.0)
                    ;yfit2=bspline_valu(waveB, sset2)
                    ;yfit[p]=yfit2[p]
                    yrange=[0.0,0.002]
                    ;yrange=[min(yfit,max=m)*0.75,m*1.25]
                endif else if instMode eq 31 then begin
                    minw=258.0
                    maxw=2900.0
                    nnodes=30
                    nord=3
                    k=where(kappa lt 1d6, count)
                    p = where(waveB[k] ge minw and waveB[k] le maxw,count2)
                    if count2 eq 0 then p=lindgen(n_elements(k))
                    ; ignore data outside of 5 sigma from a 3rd order polynomial fit
                    coeff=robust_poly_fit(waveB[k[p]], kappa[k[p]], 6d)
                    resistant_mean, (kappa[k[p]] - poly(waveB[k[p]],coeff)), 5.0,mean,goodvec=keep0
                    sset=bspline_iterfit(waveB[k[p[keep0]]], kappa[k[p[keep0]]], requiren=nnodes,bkspace=nnodes, nord=nord)
                    yfit=bspline_valu(waveB, sset)

                    yrange=[0.0,0.004]
                endif
                plot,waveB[k],kappa[k],yrange=yrange,/xst,/yst,psym=-4, xtitle='Wavelength (nm)',ytitle='Kappa',charsize=1.4
                oplot,waveB,yfit,color=2, thick=2.0
                xyouts,0.70,0.70,strtrim(string(gps2sd(timeB[0]/1d6),format='(F0.2)'),2),/normal,color=2,charsize=1.5
                ;read,'SKIP (s) or Keep :', answer
                ;if strupcase(answer) ne 'S' then $
                all_kappa = [all_kappa, ptr_new({time:median(timeB), wavelength:waveB, kappa:kappa, kappa_fit:yfit},/no_copy)]
            endif

        endfor
    endif   ; all_kappa

    ; smooth the data IN TIME by fitting a bspline for each wavelength
    ; first create a 2D array of kappas
    timesd=dblarr(n_elements(all_kappa))
    trange=[0d,1d6]
    for i=0,n_elements(all_kappa)-1 do timesd[i]=(*all_kappa[i]).time
    timesd=gps2sd(timesd/1d6)
    if instMode eq 41 then begin
        trange=[40d, 6000d]
    endif else if instMode eq 43 then begin
        trange=[40d, 6000d]
    endif else if instMode eq 31 then begin
        trange=[40d, 3000d]
    endif
    tpos=where(timesd ge trange[0] and timesd le trange[1],count)
    if count gt 0 then timesd=timesd[tpos] else tpos=lindgen(n_elements(timesd))

    kappa_smooth=dblarr(n_elements(timesd), n_elements(waves))
    ;for d=0L,n_elements(timesd)-1L do begin
    ;    coeff=robust_poly_fit((*all_kappa[tpos[d]]).wavelength, (*all_kappa[tpos[d]]).kappa_fit, 6d)
    ;    resistant_mean, ((*all_kappa[tpos[d]]).kappa_fit-poly((*all_kappa[tpos[d]]).wavelength,coeff)), 5.0,mean,goodvec=kk
    ;    ;val=interpol((*all_kappa[tpos[d]]).kappa_fit[p], (*all_kappa[tpos[d]]).wavelength[p], waves,/lsq)
    ;    sset=bspline_iterfit((*all_kappa[tpos[d]]).wavelength[kk],(*all_kappa[tpos[d]]).kappa_fit[kk],requiren=10,bkspace=10,nord=3)
    ;    val=(bspline_valu(waves,sset)) > 1d-7
    ;    kappa_smooth[d,*]=val
    ;endfor

    ;use the previous fit instead of a much coarser fit
    for d=0L,n_elements(timesd)-1L do kappa_smooth[d,*]=(*all_kappa[d]).kappa_fit

    ; fit a bspline for each wavelength as a function of time
    s=sort(timesd)
    timesd=timesd[s]
    kappa_smooth=kappa_smooth[s,*]

    if n_elements(timesd) lt 5 then begin
        kappa_smooth = {timesd:timesd, waves:waves, kappa:kappa_smooth}
        return, all_kappa
    endif

    wid=10
    window,wid,xsize=1200,ysize=600
    answer = ''

    for w=0L,n_elements(waves)-1L do begin
        ;fit the data in log space
        ;td=where(timesd ge 180d, comp=ctd)
        ;ctd=where(timesd le 460d)
        ;coeff=robust_poly_fit(timesd[td], alog10(kappa_smooth[td,w]), 6d)
        ;yfit=10d^poly(timesd, coeff)
        
        ;p=where(kappa_smooth[ctd,w] gt 0d,count)
        ;coeff2=ladfit(timesd[ctd[p]], kappa_smooth[ctd[p],w])
        ;yfit2 = poly(timesd[ctd], coeff2)
        ;if count ge 7 then begin
        ;    coeff2=robust_poly_fit(timesd[ctd[p]], kappa_smooth[ctd[p],w], 6d) 
        ;endif else if count gt 2 then begin
        ;    coeff2=robust_poly_fit(timesd[ctd[p]], kappa_smooth[ctd[p],w], count-1)
        ;endif else coeff2=ladfit(timesd[ctd],kappa_smooth[ctd,w])
        ;yfit2=10d^poly(timesd[ctd], coeff)
        ;resistant_mean, (kappa_smooth[ctd,w] - poly(timesd[ctd],coeff2)), 6.0,mean,goodvec=keep0
        ;p=where(kappa_smooth[keep0,w] gt 0.0, comp=cp, count)
        ;if cp[0] ne -1 then keep0=keep0[p]
        ;sset=bspline_iterfit(timesd[ctd[keep0]], kappa_smooth[ctd[keep0],w], requiren=2,bkspace=2, nord=3)
        ;td2=where(timesd le 400d)
        ;yfit2 = bspline_valu(timesd[td2], sset)
        ;yfit2 = bspline_valu(timesd[ctd], sset)

        ;sset=bspline_iterfit(timesd, kappa_smooth[*,w], requiren=60,bkspace=60, nord=3)
        ;resistant_mean, (kappa_smooth[*,w] - bspline_valu(timesd,sset)), 6.0,mean,goodvec=keep0
        ;sset=bspline_iterfit(timesd[keep0], kappa_smooth[keep0,w], requiren=20,bkspace=20, nord=3)
        ;yfit=bspline_valu(timesd,sset)
        ;p=where(kappa_smooth[*,w] gt 0d)

        if instmode eq 41 then begin
            q=where(timesd gt 118d)                                                                                         
            sset=bspline_iterfit(timesd[q], kappa_smooth[q,w], requiren=20,bkspace=20, nord=3)                                  
            resistant_mean, (kappa_smooth[*,w]) - bspline_valu(timesd,sset), 6.0,mean,goodvec=keep0
            p=where(timesd[keep0] gt 118d)
            sset=bspline_iterfit(timesd[keep0[p]], kappa_smooth[keep0[p],w], requiren=20,bkspace=20, nord=5)
            yfit = bspline_valu(timesd,sset)
            xrange=[0d,5100]
            yrange=[-0.0001,0.0020]
        endif else if instMode eq 43 then begin
            q=where(timesd gt 160,comp=cq)                                                                                         
            sset=bspline_iterfit(timesd[q], alog10(kappa_smooth[q,w]), requiren=40,bkspace=40, nord=3)                                  
            resistant_mean, (alog10(kappa_smooth[q,w]) - bspline_valu(timesd[q],sset)), 6.0,mean,goodvec=keep0
            nnodes=20
            sset=bspline_iterfit(timesd[q[keep0]], alog10(kappa_smooth[q[keep0],w]), requiren=nnodes,bkspace=nnodes, nord=3)
            yfit = 10d^bspline_valu(timesd,sset)
            sset=bspline_iterfit(timesd, [kappa_smooth[cq,w],yfit[q]], requiren=2,bkspace=2, nord=3)
            yfit = bspline_valu(timesd,sset)
            xrange=[0d,6000]
            yrange=[0.001,0.008]
        endif else if instMode eq 31 then begin
            q=where(timesd gt 100)                                                                                         
            sset=bspline_iterfit(timesd[q], alog10(kappa_smooth[q,w]), requiren=10,bkspace=10, nord=3)                                  
            resistant_mean, (alog10(kappa_smooth[*,w]) - bspline_valu(timesd,sset)), 6.0,mean,goodvec=keep0
            nnodes=10
            sset=bspline_iterfit(timesd[keep0], alog10(kappa_smooth[keep0,w]), requiren=nnodes,bkspace=nnodes, nord=3)
            yfit = 10d^bspline_valu(timesd,sset)
            xrange=[0d,3000]
            yrange=[-0.0002,0.003]
        endif
        
        plot,timesd,kappa_smooth[*,w],yrange=yrange,xrange=xrange,/xst,/yst,psym=-4, xtitle='!3Mission Days!X',$
            ytitle='!3Kappa!X',title='!3SORCE SIM Kappa!X',charsize=2.0, font=-1
        oplot,timesd,yfit,color=2,thick=2.0
        ;oplot,timesd[ctd],yfit2,color=2
        xyouts,0.70,0.70,'!3'+strtrim(string(waves[w],format='(F0.2," nm")'),2)+'!X',/normal,color=2,charsize=2.0, font=-1
        wait, 0.05
        if n_elements(kappa_smooth[*,w]) eq n_elements(yfit) then begin
            kappa_smooth[*,w] = yfit
            ;kappa_smooth[ctd,w] = yfit2
        endif
    endfor


    ; fit a bspline for each day
    if instMode eq 41 then begin
        ;for d=0L,n_elements(timesd)-1L do begin
        ;    sset=bspline_iterfit(waves, kappa_smooth[d,*], requiren=20,bkspace=20, nord=3)
        ;    kappa_smooth[d,*] = bspline_valu(waves, sset)
        ;endfor
    endif else if instMode eq 43 then begin
        ; for the UV, the fit for the days prior to day 300 are not very good
        ; we also need to bspline only to wavelength 306.3 and make sure we have no Kappa<0.0
        p=where(timesd lt 400d,count)
        k=where(waves lt 306.3)
        for d=0L,count-1L do begin
            sset=bspline_iterfit(waves[k],kappa_smooth[d,k],requi=2,bkspace=2,nord=3) 
            val=bspline_valu(waves,sset)
            kappa_smooth[d,*]=val
        endfor
    endif else if instMode eq 31 then begin
        ; set the value of Kappa to its minimum where it starts increasing again 
        ; (the Kappa should NOT increase at long wave)
        for d=0L,n_elements(timesd)-1L do begin
            mn=min(kappa_smooth[d,*],pos)
            if pos ne (n_elements(kappa_smooth[d,*])-1L) then kappa_smooth[d,pos:-1] = kappa_smooth[d,pos]
            sset=bspline_iterfit(waves, kappa_smooth[d,*], requiren=10,bkspace=10, nord=3)
            kappa_smooth[d,*] = bspline_valu(waves, sset)
        endfor
    endif
    s=sort(timesd)
    timesd=timesd[s]
    kappa_smooth=kappa_smooth[s,*]
    ; limit the valid range of kappa (true for every detector)
    kappa_smooth = (kappa_smooth > 0.0d) < 0.01d

    ; TODO
    ; to extrapolate further out in the mission, we can get a new low order bspline and extend the time
    ; sset=bspline_iterfit(kappa_smooth.timesd,kappa_smooth.kappa[*,100],requiren=10,bkspace=10, nord=2)
    ; new_timesd=findgen(5000d - timesd[-1])+timesd[-1]+0.5
    ; yfit = bspline_valu(new_timesd,sset)
    ; timesd=[timesd, new_timesd]
    ; kappa=[kappa,yfit]

    if keyword_set(regular) then begin
        ; if regular is set, resample the kappa_smooth array with linear xgrid and ygrid
        ; (for much easier interpolation using interpole)
        ; with the new regularly spaced array, we can interpolate for a specific day using:
        ;
        ; sd0=specta.spect20[i].timestamp
        ; wv0=specta.spect20[i].wavelength
        ; xslope = (max(kappa.timesd)-min(kappa.timesd)) / (n_elements(kappa.timesd)-1.)
        ; yslope = (max(kappa.waves)-min(kappa.waves)) / (n_elements(kappa.waves)-1.)
        ; x2 = (sd0 - min(kappa.timesd)) / xslope
        ; y2 = (wv0 - min(kappa.waves)) / yslope
        ; kappa_val = interpolate(kappa.kappa, x2, y2, /cubic)
        ;
        ;  OR by simply calling inter2d
        ;
        ; sd0=specta.spect20[i].timestamp
        ; wv0=specta.spect20[i].wavelength
        ; kappa_val=interp2d(kappa.kappa, kappa.timesd, kappa.waves, sd0, wv0, /double, /cubic)
        ;
        nx=n_elements(timesd)
        ny=n_elements(waves)
        nxny=[nx,ny]
        gs = [(max(timesd)-min(timesd))/(nxny[0]-1), (max(waves)-min(waves))/(nxny[1]-1)]
        x0 = rebin(temporary(timesd), nx, ny, /sample)
        y0 = transpose(rebin(temporary(waves), ny, nx, /sample))
        triangulate, x0, y0, tr,bound
        zz = trigrid(x0,y0, kappa_smooth, tr, gs, extrapolate=bound, /quintic,xgrid=xgrid,ygrid=ygrid)
        kappa_smooth={timesd:xgrid, waves:ygrid, kappa:zz}
    endif else begin
        kappa_smooth={timesd:timesd, waves:waves, kappa:kappa_smooth}
    endelse

    ; srf=surface(kappa_smooth.kappa,kappa_smooth.timesd,kappa_smooth.waves)

    return, all_kappa

END

