;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Compare the ESRA and ESRB for a list of wavelengths and during the 
;   specified time range and determine the value of Tau (degradation). 
;   We calculate a change in the Kappa as a function of mission time 
;   AND wavelength.  If a Kappa function is provided, we measure the 
;   change of this reference Kappa as a function of mission day (and
;   wavelength).
;   The associated raypath AFact is assumed constant for the time period 
;   specified (adjustments to the AFact will be left for later processing).
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
;   kappa -
;      2D Array containing the wavelength and the corresponding reference
;      Kappa value.
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
;    We are using a version of the processed ESR data where the SimUncorrectedIrradiance
;    has been 1AU corrected (version 1004 of SORCE_SIM_V20 for example).
;
;    We use the following expression to represent the degradation (at a specific wavelength and mission day)
;    where afact is the raypath:
;
;    SolarIrrad = MeasuredIrrad / ((1d - afact)*exp(-abs(Kappa * FFunc * solexp)) + afact*exp(-abs(Kappa * FFunc * solexp31)*0.5))
;
; REVISION HISTORY:
;   Revision: $Id: get_esr_tau.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;------------------------------------------------------------------
;
function mytnmin_esrtau, x, df, solexp31=solexp31, ird31=ird31, solexp32=solexp32, ird32=ird32, $
    time31=time31, time32=time32, new_ird31=new_ird31, new_ird32=new_ird32, afact31=afact31, afact32=afact32
    ; measured_irrad = SolarIrrad * exp(-abs(x[1]*solexp))
    ; SolarIrrad = measured_irrad / exp(-abs(x[1]*solexp))
    new_ird31 = ird31 / ((1d - afact31)*exp(-abs(x[0]*solexp31)) + afact31*exp(-abs(x[0]*solexp31)*0.5d))
    new_ird32 = ird32 / ((1d - afact32)*exp(-abs(x[0]*solexp32)) + afact32*exp(-abs(x[0]*solexp32)*0.5d))

    ; since the ESR Table scan data is at different times for A and B
    ; do a 2nd order fit to each iMode and then use the diff at same time
    if n_elements(time31) gt 3 then begin
        coeffs31 = robust_poly_fit(time31, new_ird31, 2d, /double)
    endif else begin
        coeffs31 = ladfit(time31, new_ird31, /double)
    endelse
    if n_elements(time31) gt 3 then begin
        coeffs32 = robust_poly_fit(time32, new_ird32, 2d, /double)
    endif else begin
        coeffs32 = ladfit(time32, new_ird32, /double)
    endelse
    if n_elements(new_ird31) ge n_elements(new_ird32) then begin
        fit31 = poly(time31,coeffs31)
        fit32 = poly(time31,coeffs32)
    endif else begin
        fit31 = poly(time32,coeffs31)
        fit32 = poly(time32,coeffs32)
    endelse
    F = fit32 - fit31
    ;out_value = abs(mean(F))
    out_value = stddev(F)
    ;out_value = MEANABSDEV(F, /median)
    ;out_value = robust_sigma(F, /zero)
    ;out_value = total(abs(f))
    ;print,x[0],x[1],out_value
    return, out_value
end
;*******************************************************************

function get_esr_tau, starttime, stoptime, wavelength, version=version, bin55=bin55, step55=step55, $
    dburl=dburl, dbdriver=dbdriver, user=user, password=password, noplot=noplot, $
    gps=gps, missionDays=missionDays, julianDays=julianDays, inspectrum=inspectrum, $
    solar54=solar54, solar55=solar55, kappa=kappa, smooth=smooth, alignobc=alignobc, afact=afact

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
    matchingpos = [[10000, 52156], [10232, 51473], [10916, 51056], [11950, 50441], [13075, 49318], $
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

    if n_elements(version) eq 0 then version=1004
    if n_elements(user) eq 0 then     user='sbeland'
    if n_elements(password) eq 0 then password='sbeland1'
    if n_elements(dburl) eq 0 then    dbUrl='jdbc:sybase:Tds:lasp-db-dev:4100/SORCE_SIM_V20'
    if n_elements(dbdriver) eq 0 then dbDriver='com.sybase.jdbc3.jdbc.SybDriver'
    if n_elements(bin55) eq 0 then bin55=1
    if n_elements(wavelength) eq 0 then begin
        ; get the list of wavea from matchingpos
        match, matchingpos[0,*], prismposa, suba, subb
        wavelength = wavea[subb]
    endif
    if n_elements(kappa) eq 0 then begin
        readcol,'~/SORCE/data/ESR_kappa_453_1570_1004.txt', wave_fit, kappa_fit,format='(d,d)'
    endif
    if n_elements(afact) eq 0 then begin
        readcol,'~/SORCE/data/overlap_esr_fit.txt',wave,a0,a1,format='(d,d,d)'
        afact={wavelength:wave, a:a1}
    endif

    ; get the SimSolarExposureData for modes 31 and 32
    if n_elements(solar54) eq 0 then begin
        q1 = 'SELECT orbitStartTimeGps t0, orbitEndTimeGps t1, cumulativeSolarExpHrtOut solar_exp, '+$
             'solarExposureHrtOut as solar_exp_orbit from SimSolarExposureData where instrumentModeId=54'
        query_database, q1, solar54, info
    endif
    ; organize the solar exposure for the beginning and end of orbit (solar_exp is for end of orbit)
    time54 = [solar54.t0, solar54.t1]
    solarexp54=[0.0,solar54[0:-2].solar_exp, solar54.solar_exp]
    s=sort(time54)
    time54=time54[s]
    solarexp54=solarexp54[s]

    if n_elements(solar55) eq 0 then begin
        q1 = 'SELECT orbitStartTimeGps t0, orbitEndTimeGps t1, cumulativeSolarExpHrtOut solar_exp, '+$
             'solarExposureHrtOut as solar_exp_orbit from SimSolarExposureData where instrumentModeId=55'
        query_database, q1, solar55, info
    endif
    ; organize the solar exposure for the beginning and end of orbit (solar_exp is for end of orbit)
    time55 = [solar55.t0, solar55.t1]
    solarexp55=[0.0,solar55[0:-2].solar_exp, solar55.solar_exp]
    s=sort(time55)
    time55=time55[s]
    solarexp55=solarexp55[s]

    if size(inspectrum,/tname) ne 'STRUCT' then begin
        esrA=get_science_product(['Wavelength','SimUncorrectedIrradiance','SimConvertedDataNumbers'],t0,t1,31,/gps,$
            version=version,user=user,password=password,dburl=dburl,dbdriver=dbdriver)
        esrB=get_science_product(['Wavelength','SimCorrectedIrradiance','SimConvertedDataNumbers'],t0,t1,32,/gps,$
            version=version,user=user,password=password,dburl=dburl,dbdriver=dbdriver)
        inspectrum={sp31:esrA, sp32:esrB}
    endif

    ; we can get the corresponding Tau (Kappa*FFunc*solarExp) for each wavelength
    ; for the timerange by either measuring kappa for all wavelengths for a timeslice
    ; and doing a bspline of the kappa
    ; OR
    ; measuring the kappa at the specific wavelength for all timeslices and smoothing
    ; the FFUNC.
    ; Do we bspline Kappa or do we bspline FFUNC? 
    ; Or do we generate a 3D data grid and smooth that? 

    ; loop over the SimB data and find the corresponding SimA data (loop over wavelengths)
    parinfo=replicate({value:0.0d, fixed:0, limited:[0,0], limits:[0d,0d], step:0d, tnside:2},1)
    parinfo[0].value=1.0d-4

    data_wave=dblarr(n_elements(inspectrum.sp32)*2)
    data_time=data_wave
    data_tau=data_wave

    p0=0L
    p1=1L
    for w=0L,n_elements(wavelength)-1 do begin
        ; get all the data points from ESRB at the matching wavelength
        mn=min(abs(waveb-wavelength[w]),posb)
        matchpos = where(matchingpos[1,*] eq prismposb[posb],count)
        wb = where(inspectrum.sp32.prismposition eq matchingpos[1,matchpos[0]],count)
        s=sort(inspectrum.sp32[wb].MICROSECONDSSINCEGPSEPOCH)
        wb=wb[s]
        print,'   processing ',strtrim(string(median(inspectrum.sp32[wb].wavelengthref),format='(F0.2)'),2)+'nm ...'

        ; find the corresponding data from ESRA
        wa = where(inspectrum.sp31.prismposition eq matchingpos[0,matchpos[0]],count)
        s=sort(inspectrum.sp31[wa].MICROSECONDSSINCEGPSEPOCH)
        wa=wa[s]

        ; work over a pair of points and calculate a Tau to aligned ESRA and ESRB
        for pb=0L,n_elements(wb)-2 do begin
            ; look for points within one day and average
            pos=where(abs(inspectrum.sp32[wb[pb]].(0) - inspectrum.sp32[wb].(0)) lt 86400d6,count)
            if count gt 0 then begin
                time32 = mean(inspectrum.sp32[wb[pos]].MICROSECONDSSINCEGPSEPOCH)
                irrad32 = mean(inspectrum.sp32[wb[pos]].IRRADIANCE)
                wave32 = mean(inspectrum.sp32[wb[pos]].wavelengthref)
            endif else begin
                time32 = inspectrum.sp32[wb[pb]].MICROSECONDSSINCEGPSEPOCH
                irrad32 = inspectrum.sp32[wb[pb]].IRRADIANCE
                wave32 = inspectrum.sp32[wb[pb]].wavelengthref
            endelse
            pos=where(abs(inspectrum.sp32[wb[pb+1]].(0) - inspectrum.sp32[wb].(0)) lt 86400d6,count)
            if count gt 0 then begin
                time32 = [time32,mean(inspectrum.sp32[wb[pos]].MICROSECONDSSINCEGPSEPOCH)]
                irrad32 = [irrad32,mean(inspectrum.sp32[wb[pos]].IRRADIANCE)]
                wave32 = [wave32,mean(inspectrum.sp32[wb[pos]].wavelengthref)]
            endif else begin
                time32 = [time32,inspectrum.sp32[wb[pb+1]].MICROSECONDSSINCEGPSEPOCH]
                irrad32 = [irrad32,inspectrum.sp32[wb[pb+1]].IRRADIANCE]
                wave32 = [wave32,inspectrum.sp32[wb[pb+1]].wavelengthref]
            endelse

            ; try interpolating the SimA data for the same time (instead of taking the closest point)
            time31 = time32
            irrad31 = interpol(inspectrum.sp31[wa].irradiance, inspectrum.sp31[wa].MICROSECONDSSINCEGPSEPOCH, time31)
            pos=where(inspectrum.sp31[wa].MICROSECONDSSINCEGPSEPOCH gt time31,count)
            wave31 = [inspectrum.sp31[wa[pos[0]-1]].wavelengthref, inspectrum.sp31[wa[pos[0]]].wavelengthref]
            ;pos=where(abs(inspectrum.sp31[wa].(0) - time32[0]) lt 86400d6,count)
            ;if count gt 0 then begin
            ;    time31 = mean(inspectrum.sp31[wa[pos]].MICROSECONDSSINCEGPSEPOCH)
            ;    irrad31 = mean(inspectrum.sp31[wa[pos]].IRRADIANCE)
            ;    wave31 = mean(inspectrum.sp31[wa[pos]].wavelengthref)
            ;endif else begin
            ;    ; get closest data point
            ;    mn=min(abs(inspectrum.sp31[wa].(0) - time32[0]),pos)
            ;    time31 = inspectrum.sp31[wa[pos]].MICROSECONDSSINCEGPSEPOCH
            ;    irrad31 = inspectrum.sp31[wa[pos]].IRRADIANCE
            ;    wave31 = inspectrum.sp31[wa[pos]].wavelengthref
            ;endelse
            ;pos=where(abs(inspectrum.sp31[wa].(0) - time32[1]) lt 86400d6,count)
            ;if count gt 0 then begin
            ;    time31 = [time31,mean(inspectrum.sp31[wa[pos]].MICROSECONDSSINCEGPSEPOCH)]
            ;    irrad31 = [irrad31, mean(inspectrum.sp31[wa[pos]].IRRADIANCE)]
            ;    wave31 = [wave31, mean(inspectrum.sp31[wa[pos]].wavelengthref)]
            ;endif else begin
            ;    mn=min(abs(inspectrum.sp31[wa].(0) - time32[1]),pos)
            ;    time31 = [time31, inspectrum.sp31[wa[pos]].MICROSECONDSSINCEGPSEPOCH]
            ;    irrad31 = [irrad31, inspectrum.sp31[wa[pos]].IRRADIANCE]
            ;    wave31 = [wave31, inspectrum.sp31[wa[pos]].wavelengthref]
            ;endelse

            afact31=interpol(afact.a, afact.wavelength, wave31, /spline)
            ;afact32=interpol(afact.a, afact.wavelength, wave32, /spline)
            afact32=interpol(afact.a, afact.wavelength, inspectrum.sp32[wb[pb:pb+1]].wavelengthref, /spline)

            solarexp31 = interpol(solarexp54, time54, time31)
            solarexp31/=86400d
            solarexp32 = interpol(solarexp55, time55, time32)
            solarexp32/=86400d

            functargs = {solexp31:solarexp31, ird31:irrad31, solexp32:solarexp32, ird32:irrad32, $
                         time31:time31, time32:time32, afact31:afact31, afact32:afact32}

            coeffs = tnmin('mytnmin_esrtau', functargs=functargs, bestmin=f0, status=status, $
                nfev=nfev, niter=niter, autoderivative=1, errmsg=errmsg, parinfo=parinfo, /quiet)

            data_wave[[p0,p1]] = inspectrum.sp32[wb[pb:pb+1]].wavelengthref
            data_time[[p0,p1]] = time32
            data_tau[[p0,p1]]  = [coeffs[0],coeffs[0]]
            p0=p1+1L
            p1=p0+1L
        endfor

    endfor

    p=where(data_time eq 0.0,comp=cp,count)
    if count gt 0 then begin
        data_wave=data_wave[cp]
        data_time=data_time[cp]
        data_tau=data_tau[cp]
    endif

    return,{MICROSECONDSSINCEGPSEPOCH:data_time, wavelength:data_wave, tau:data_tau}

end
