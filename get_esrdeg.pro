;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Extract the degradation information from the ESR table scans matching
;   the data for the requested wavelength taken at the "same" time for 
;   ESRA and ESRB.  This is done for two consecutive table scans of ESRB
;   to get a running slope of the degradation with mission time.
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
;    The updated version wull use the fully processed ESR data and remove the prism degradation
;    found in the table SimPrismTransDegradation without having to run a special version with
;    no degradation.
;
;
; REVISION HISTORY:
;   Revision: $Id: get_esrdeg.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;------------------------------------------------------------------
function get_esrdeg, starttime, stoptime, wavelength, version=version, bin55=bin55, step55=step55, $
    dburl=dburl, dbdriver=dbdriver, user=user, password=password, noplot=noplot, $
    gps=gps, missionDays=missionDays, julianDays=julianDays, inspectrum=inspectrum, $
    solar54=solar54, solar55=solar55, kappa_coeff=kappa_coeff, smooth=smooth, alignobc=alignobc

    ; list of wavelengths from ESRA table scan
    wavea = [    261.2d, 265.2, 279.5,  280.4,  283.5,  288,    293.5,  304.5,  319,    333.5,  342,    355,    369, $
                 375,    384.5,   395,    407,    431,    471,    480,    488,    519.5,  567, $
                 594,      659,   710,    758,    806,    859,    867.5,  892.5,  901,    964.5, $
                1013,   1063.5,  1116,   1170,   1213.5, 1279.5, 1356.5, 1411,   1475.5, 1497, $
                1549.5,   1591,1621.5,   1692,   1741,   1817,   1882,   1909,   1944.5, 1988.5, $
                2014.25,2106.5,2171.2, 2287.5, 2397.5, 2502.5]

    prismposa = [ 5050L,   6250,  10000,  10232,  10916,  11950,  13075,  15100,  17425,  19375, 20350, 21775, 23125, $
                 23650,   24400,  25150,  25975,  27400,  29275,  29650,  29950,  31000,  32275, $
                 32875,   34075,  34825,  35425,  35950,  36475,  36550,  36775,  36850,  37375, $
                 37750,   38125,  38500,  38875,  39175,  39625,  40150,  40525,  40975,  41125, $
                 41500,   41800,  42025,  42550,  42925,  43525,  44050,  44275,  44575,  44950, $
                 45175,   46000,  46600,  47725,  48850,  49975]

    ; list of wavelengths from ESRB table scan
    waveb = [   278.4d,   279.3,  282.3,    284,    287,  292.5,  303.3,  317.7,    332,    340,    346.5,  $
                   353,     367,    373,    382,    392,    395,    404,    428,    467,    476,    $
                   484,     515,    561,    587,    650,  699.5,    746,    792,  844.5,  852.5,  $
                   877,   885.5,    948,    996, 1046.5,   1099, 1153.5, 1197.5,   1264,   1341,  $
                  1396,    1459,   1480, 1533.5,   1576,   1607, 1678.5, 1728.5, 1801.5, 1866.5, $
                  1894,  1930.5, 1974.5, 2000.5,   2094,   2159, 2276.5, 2387.5, 2493.5]

    prismposb = [52388L,  52156,  51473,  51056,  50441,  49318,  47296,  44975,  43029,  42055, 41307, $
                 40633,   39286,  38762,  38014,  37265,  37041,  36442,  35020,  33149,  32775, $
                 32475,   31428,  30156,  29557,  28360,  27612,  27014,  26490,  25966,  25891, $
                 25667,   25592,  25069,  24695,  24321,  23946,  23572,  23273,  22824,  22301, $
                 21927,   21478,  21329,  20955,  20655,  20431,  19907,  19533,  18935,  18412, $
                 18188,   17888,  17514,  17290,  16467,  15869,  14748,  13626,  12504]

    ; we'll match the following prism positions between A and B (corresponding to closests wavelengths)
    matchingpos = [[10000L,52388L],[10000, 52156], [10232, 51473], [10916, 51056], [11950, 50441], [13075, 49318], $
                   [15100, 47296], [17425, 44975], [19375, 43029], [20350, 42055], [21775, 40633], $
                   [23125, 39286], [23650, 38762], [24400, 38014], [25150, 37041], [25975, 36442], $
                   [27400, 35020], [29275, 33149], [29650, 32775], [29950, 32475], [31000, 31428], $
                   [32275, 30156], [32875, 29557], [34075, 28360], [34825, 27612], [35425, 27014], $
                   [35950, 26490], [36475, 25891], [36550, 25667], [36775, 25592], [37375, 25069], $
                   [37750, 24695], [38125, 24321], [38500, 23946], [38875, 23572], [39175, 23273], $
                   [39625, 22824], [40150, 22301], [40525, 21927], [40975, 21478], [41125, 21329], $
                   [41500, 20955], [41800, 20655], [42025, 20431], [42550, 19907], [42925, 19533], $
                   [43525, 18935], [44050, 18412], [44275, 18188], [44575, 17888], [44950, 17514], $
                   [45175, 17290], [46000, 16467], [46600, 15869], [47725, 14748], [48850, 13626], $
                   [49975, 12504]]


    obctimes = [0.0d, 1570.0d, 2173d, 2455d, 2803d, 5000d]

    if keyword_set(missionDays) then begin
       ; user specified time in mission (sorce) days
       t0 = sd2gps(startTime)*1.d6
       t1 = sd2gps(stopTime)*1.d6
    endif else if keyword_set(julianDays) then begin
       ; user specified time in julian days
       t0 = jd2gps(startTime)*1.d6
       t1 = jd2gps(stopTime)*1.d6
    endif else begin
       ; user specified timetags in gps microseconds
       t0 = startTime
       t1 = stopTime
    endelse

    if n_elements(version) eq 0 then version=70
    if n_elements(user) eq 0 then     user='sbeland'
    if n_elements(password) eq 0 then password='sbeland1'
    if n_elements(dburl) eq 0 then    dbUrl='jdbc:sybase:Tds:lasp-db-dev:4100/SORCE_SIM_V19'
    if n_elements(dbdriver) eq 0 then dbDriver='com.sybase.jdbc3.jdbc.SybDriver'
    ; if n_elements(kappa_coeff) lt 3 then kappa_coeff=[0.0095408430d, -0.0070394634d, 1.7073941d-05]
    ; new Kappa coefficients from 20 iterations over Kappa and F-Function
    ;if n_elements(kappa_coeff) lt 3 then kappa_coeff=[0.023980000d,   -0.0070074030d,   3.7123102d-05]
    ;if n_elements(kappa_coeff) lt 3 then kappa_coeff=[0.0098435447d, -0.0070887388d, 1.0d-05]
    if n_elements(kappa_coeff) lt 3 then kappa_coeff=[0.0096235817d,   -0.0071923836d,   1.7092305d-05]
    if n_elements(bin55) eq 0 then bin55=1

    ; get the SimSolarExposureData for modes 31 and 32
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
        if count54 eq 0 then begin
            oneau55 = interpol(solardist.SUNOBSERVERDISTANCECORRECTION, solardist.MICROSECONDSSINCEGPSEPOCH, solar55.(1))
            append_tag,solar55,'oneau',oneau55,/slim
        endif
        ; apply the 1AU correction (the factor 4 was determine empirically)
        solar54.solar_exp=total(solar54.solar_exp_orbit,/cum) / (1d +(solar54.oneau-1d)/4d)
        solar55.solar_exp=total(solar55.solar_exp_orbit,/cum) / (1d +(solar55.oneau-1d)/4d)
    endif

    if size(inspectrum,/tname) ne 'STRUCT' then begin
        ; get all of the ESR data within the specified timerange and only keep around the requested wavelength
        ;esrA=get_science_product(['Wavelength','SimCalibratedIrradiance'],t0,t1,31,/gps,$
        ;esrA=get_science_product(['SimProfileIntegral','SimCalibratedIrradiance','SimConvertedDataNumbers'],t0,t1,31,/gps,$
        esrA=get_science_product(['Wavelength','SimCorrectedIrradiance','PhaseDetectedDN','SimPrismTransDegradation'],t0,t1,31,/gps,$
            version=version,user=user,password=password,dburl=dburl,dbdriver=dbdriver)
        if size(esrA,/tname) eq 'STRUCT' then begin
            ; remove the prism transmission degradation
            esrA.irradiance *= esrA.prismTransDegradation
        endif else begin
            ; skip the prismtransdeg which will be missing if we are looking at a No-Degradation version
            esrA=get_science_product(['Wavelength','SimCorrectedIrradiance','PhaseDetectedDN'],t0,t1,31,/gps,$
                version=version,user=user,password=password,dburl=dburl,dbdriver=dbdriver)
        endelse

        ;esrB=get_science_product(['Wavelength','SimCalibratedIrradiance'],t0,t1,32,/gps,$
        ;esrB=get_science_product(['SimProfileIntegral','SimCalibratedIrradiance','SimConvertedDataNumbers'],t0,t1,32,/gps,$
        esrB=get_science_product(['Wavelength','SimCorrectedIrradiance','PhaseDetectedDN','SimPrismTransDegradation'],t0,t1,32,/gps,$
            version=version,user=user,password=password,dburl=dburl,dbdriver=dbdriver)
        if size(esrB,/tname) eq 'STRUCT' then begin
            ; remove the prism transmission degradation
            esrB.irradiance *= esrB.prismTransDegradation
        endif else begin
            esrB=get_science_product(['Wavelength','SimCorrectedIrradiance','PhaseDetectedDN'],t0,t1,32,/gps,$
                version=version,user=user,password=password,dburl=dburl,dbdriver=dbdriver)
        endelse
        inspectrum={sp31:esrA, sp32:esrB}
    endif

    ; find closest wavelength from the table scan
    mn=min(abs(waveb-wavelength),posb)
    matchpos = where(matchingpos[1,*] eq prismposb[posb],count)
    wb = where(inspectrum.sp32.prismposition eq matchingpos[1,matchpos[0]],count)
    if count eq 0 then begin
        print,'Error: no ESRB data found at prism position ',matchingpos[1,matchpos]
        print,'       corresponding to wavelength ',waveb[posb]
        return,-1
    endif
    ; find closest wavelength from the table scan for ESRA
    posa=where(prismposa eq matchingpos[0,matchpos[0]])
    wa = where(inspectrum.sp31.prismposition eq matchingpos[0,matchpos[0]],count)
    if count eq 0 then begin
        print,'Error: no ESRA data found at prism position ',matchingpos[0,matchpos]
        print,'       corresponding to wavelength ',wavea[posa]
        return,-1
    endif
    print,'processing closest table scan wavelength for ESRA,B '+ $
        strtrim(string(wavea[posa],format='(F0.2)'),2)+', '+strtrim(string(waveb[posb],format='(F0.2)'),2)


    irrad32=[]
    time32=[]
    wave32=[]
    solarexp32=[]
    p=where(strpos(tag_names(inspectrum.sp32),'PRISMTRANSDEGRADATION') ge 0,prismTrans)
    while n_elements(wb) gt 0 do begin
        ; look for multiple points within 1 day (from ESRFullScan) and average them
        p=where(abs(inspectrum.sp32[wb].(0) - inspectrum.sp32[wb[0]].(0)) lt 86400d6,count,comp=comp)
        if prismTrans eq 0 then $
            irrad32=[irrad32,mean(inspectrum.sp32[wb[p]].irradiance)] $
        else $
            irrad32=[irrad32,mean(inspectrum.sp32[wb[p]].irradiance * inspectrum.sp32[wb[p]].prismTransDegradation)]
        time32=[time32,mean(gps2sd(inspectrum.sp32[wb[p]].MICROSECONDSSINCEGPSEPOCH/1d6))]
        wave32=[wave32,mean(inspectrum.sp32[wb[p]].wavelengthRef)]
        q=where(solar55.t1 le sd2gps(time32[-1])*1d6,count)
        solarexp32=[solarexp32, solar55[q[-1]].solar_exp/86400d]
        if comp[0] ne -1 then wb=wb[comp] else break
    endwhile


    irrad31=[]
    time31=[]
    wave31=[]
    solarexp31=[]
    p=where(strpos(tag_names(inspectrum.sp31),'PRISMTRANSDEGRADATION') ge 0,prismTrans)
    while n_elements(wa) gt 0 do begin
        ; look for multiple points within 1 day  (from ESRFullScan) and average them
        p=where(abs(inspectrum.sp31[wa].(0) - inspectrum.sp31[wa[0]].(0)) lt 2d*3600d6,count,comp=comp)
        if prismTrans eq 0 then $
            irrad31=[irrad31,mean(inspectrum.sp31[wa[p]].irradiance)] $
        else $
            irrad31=[irrad31,mean(inspectrum.sp31[wa[p]].irradiance * inspectrum.sp31[wa[p]].prismTransDegradation)]
        time31=[time31,mean(gps2sd(inspectrum.sp31[wa[p]].MICROSECONDSSINCEGPSEPOCH/1d6))]
        wave31=[wave31,mean(inspectrum.sp31[wa[p]].wavelengthRef)]
        q=where(solar54.t1 le sd2gps(time31[-1])*1d6,count)
        solarexp31=[solarexp31, solar54[q[-1]].solar_exp/86400d]
        if comp[0] ne -1 then wa=wa[comp] else break
    endwhile



    ; if alignobc was specified, we align each segment by moving the latter piece
    ; in irradiance to match the average extrapolated points from both segments.
    if keyword_set(alignobc) then begin
        align_irrad, time31, irrad31, /mission
        align_irrad, time32, irrad32, /mission
    endif

    ; if smooth was flagged, apply a bspline smoothing (from SDSS idlutils/bspline_iterfit.pro)
    if keyword_set(smooth) then begin
        ; fit a 2nd order curve to data and remove the outliers
        ; ESRA
        coeff31=robust_poly_fit(time31,irrad31,2,yfit,/double)
        resistant_mean,(irrad31-yfit),3.0,mean,goodvec=keep0
        sset=bspline_iterfit(time31[keep0],irrad31[keep0],maxiter=0,requiren=10,bkspace=5)
        yfit31=bspline_valu(time31,sset)
        ; ESRB
        coeff32=robust_poly_fit(time32,irrad32,2,yfit,/double)
        resistant_mean,(irrad32-yfit),3.0,mean,goodvec=keep0
        sset=bspline_iterfit(time32[keep0],irrad32[keep0],maxiter=0,requiren=10,bkspace=5)
        yfit32=bspline_valu(time32,sset)
    endif else begin
        yfit31=-1
        yfit32=-1
    endelse


    ; loop through the list of ESRB data points and match a pair
    ; of points from ESRB to a pair of points from ESRA closest in
    ; clock time (to have the most similar solar irradiance when 
    ; comparing A & B)
    kappa = kappa_coeff[0] * exp(kappa_coeff[1] * wave32) + kappa_coeff[2]
    f_factor=dblarr((n_elements(wave32)-bin55)>1) -1d6
    time_f=f_factor
    
    ; force the first F_factor_scaled to be 1.0 (F/Kappa=1)  (every other value will be relative to this time)
    f_factor[0]=1d * kappa[0]
    time_f[0]=time32[0]

    for j=1L,n_elements(f_factor)-1L do begin
        f_factor[j]=f_factor[j-1]
        if n_elements(step55) eq 0 then begin
            posb0=lindgen(bin55)+j
            posb1=posb0+1
        endif else begin
            posb0=lindgen(bin55)+j
            p=where(time32 ge time32[j]+step55,count)
            if count eq 0 then continue
            posb1=lindgen(bin55)+p[0]
        endelse
        mn=min(time32[posb0],max=mx)
        ; pad the start and stop time by up to one day 
        posa0 = where(time31 ge (mn-1.5d/24d) and time31 le (mx+1.5d/24d),count)
        if count eq 0 then posa0 = where(time31 ge (mn-12d/24d) and time31 le (mx+12d/24d),count)
        if count eq 0 then posa0 = where(time31 ge (mn-1d) and time31 le (mx+1d),count)
        if count eq 0 then begin
            print,'   no matching data for ESRA for time: ',mn,mx
            continue
        endif
        mn=min(time32[posb1],max=mx)
        ; average the yfit32 over +/- 15 days
        posa1 = where(time31 ge (mn-1.5d/24d) and time31 le (mx+1.5d/24d),count)
        if count eq 0 then posa1 = where(time31 ge (mn-12d/24d) and time31 le (mx+12d/24d),count)
        if count eq 0 then posa1 = where(time31 ge (mn-1d) and time31 le (mx+1d),count)
        if count eq 0 then begin
            print,'   no matching data for ESRA for time: ',mn,mx
            continue
        endif
        ;f_factor1[j-1] = alog(irrad32[j]/irrad32[j-1]) - alog(irrad31[p1]/irrad31[p0])
        ;f_factor1[j-1] /= ((solarexp32[j-1]-solarexp32[j]) - (solarexp31[p0] - solarexp31[p1]))
        if keyword_set(smooth) then begin
            ;f_factor[j] = alog(mean(yfit32[posb1])/mean(yfit32[posb0])) - alog(mean(yfit31[posa1])/mean(yfit31[posa0]))
            ;f_factor[j] /= ((mean(solarexp32[posb0])-mean(solarexp32[posb1])) - (mean(solarexp31[posa0]) - mean(solarexp31[posa1])))
            f_factor[j] = alog(mean(yfit32[posb1])/mean(yfit32[posb0])) - alog(mean(yfit31[posa1])/mean(yfit31[posa0]))
            f_factor[j] -= f_factor[j-1] * (mean(solarexp32[posb0]) - mean(solarexp31[posa0]))
            f_factor[j] /= (mean(solarexp31[posa1]) - mean(solarexp32[posb1])) 
        endif else begin
            ;f_factor[j] = alog(mean(irrad32[posb1])/mean(irrad32[posb0])) - alog(mean(irrad31[posa1])/mean(irrad31[posa0]))
            ;f_factor[j] /= ((mean(solarexp32[posb0])-mean(solarexp32[posb1])) - (mean(solarexp31[posa0]) - mean(solarexp31[posa1])))
            f_factor[j] = alog(mean(irrad32[posb1])/mean(irrad32[posb0])) - alog(mean(irrad31[posa1])/mean(irrad31[posa0]))
            f_factor[j] -= f_factor[j-1] * (mean(solarexp32[posb0]) - mean(solarexp31[posa0]))
            f_factor[j] /= (mean(solarexp31[posa1]) - mean(solarexp32[posb1])) 
        endelse
    endfor
    p=where(f_factor gt -1d6,count)
    if count eq 0 then begin
        print,'No matching ESRA and ESRB data within 2 days of each other'
        return,-1
    endif
    f_factor=f_factor[p]
    time_f=time32[p]
    kappa=kappa[p]

    coeff=robust_poly_fit(time_f,f_factor,2,yfit,/double)
    resistant_mean,f_factor-poly(time_f,coeff),3.0,mean,goodvec=keep0
    ;keep0=lindgen(n_elements(f_factor))
    if not keyword_set(noplot) then begin
        label='F_Factor @ '+strtrim(string(mean(wave32[keep0]),format='(F0.2)'),2)+'nm  '
        label=label+'Kappa='+strtrim(string(mean(kappa),format='(E0.4)'),2)
        ;if n_elements(step55) gt 0 then label=label+'   Step='+strtrim(string(step55,format='(F0.1)'),2)+' days'
        ;plot_multi,time_f[keep0],f_factor[keep0]/kappa[keep0],psym=-4,title=label, xtitle='Mission Day',ytitle='F/Kappa',$
        ;    yrange=[0,2],xrange=[0,4000], charsize=1.4
        if n_elements(step55) gt 0 then title='Step='+strtrim(string(step55,format='(F0.1)'),2)+' days' else title=''
        title += '  '+strtrim(string(mean(wave32[keep0]),format='(F0.2)'),2)+'nm  InstMode=31'
        ;lineplot,time_f[keep0],f_factor[keep0]/kappa[keep0],psym=-3,ptitle='F_Factor', xtitle='Mission Day',ytitle='F/Kappa',$
        lineplot,time_f[keep0],f_factor[keep0]*kappa[keep0],psym=-3,ptitle='F_Factor', xtitle='Mission Day',ytitle='F*Kappa',$
            yrange=[0,2],xrange=[0,4000], charsize=1.4,title=title
    endif

    if keyword_set(smooth) then begin
        return, {time31:time31, time32:time32, irrad31:irrad31, irrad32:irrad32, yfit31:yfit31, yfit32:yfit32, $
                 solarexp31:solarexp31, solarexp32:solarexp32, time_f:time_f, f_factor:f_factor, f_factor_scaled:f_factor/kappa} 
    endif else begin
        return, {time31:time31, time32:time32, irrad31:irrad31, irrad32:irrad32,$
                 solarexp31:solarexp31, solarexp32:solarexp32, time_f:time_f, f_factor:f_factor, f_factor_scaled:f_factor/kappa}
    endelse

end
