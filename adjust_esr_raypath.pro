;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Using an existing version of the ESR data, modify the rayPath factor 
;   and apply a new one to the time series for ESRA and ESRB at all
;   table scan wavelengths.  The new rayPath value at each wavelength is 
;   ajusted so that the differences in the new "calibratedIrradiance" 
;   between ESRA and ESRB is constant for the specified time range.
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
;    We use the SimUncorrectedIrradiance as a starting point and apply the 1AU,
;    the degradationColumn, and the kappa value.
;
;
; REVISION HISTORY:
;   Revision: $Id: adjust_esr_raypath.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;------------------------------------------------------------------
function func_getslope, x, df, sda=sda, sdb=sdb, kappa_A=kappa_A, kappa_b=kappa_b, degcol_a=degcol_a, degcol_b=degcol_b, $
    irrad_a=irrad_a, irrad_b=irrad_b, new_irrad_a=new_irrad_a, new_irrad_b=new_irrad_b, delta_irrad=delta_irrad

    ; we will align the delta_irradiance per OBC event segment 
    ;obctimes = sd2gps([0.0d, 441.5, 1570.0d, 2173d, 2455d, 2803d, 5000d])*1d6

    ; quiet down the warnings in poly_fit
    prismTransDegA = (1d - x[0]) * exp(-kappa_a * degcol_a) + x[0] * exp(-kappa_a * degcol_a / 2d)
    prismTransDegB = (1d - x[0]) * exp(-kappa_b * degcol_b) + x[0] * exp(-kappa_b * degcol_b / 2d)
    new_irrad_a = irrad_a / prismTransDegA
    new_irrad_b = irrad_b / prismTransDegB
    
    if n_elements(sda) gt 3 then begin
        coeffs_a = robust_poly_fit(sda, new_irrad_a, 2d, /double)
    endif else begin
        coeffs_a = ladfit(sda, new_irrad_a, /double)
    endelse
    if n_elements(sdb) gt 3 then begin
        coeffs_b = robust_poly_fit(sdb, new_irrad_b, 2d, /double)
    endif else begin
        coeffs_b = ladfit(sdb, new_irrad_b, /double)
    endelse
    ;fitB=new_irrad_b
    ;fitA = interpol(new_irrad_a, sda, sdb)
    ;fitA = spline(sda, new_irrad_a, sdb, /double)
    ;delta_irrad = fitB - fitA
    ;coeff = ladfit(sdb, delta_irrad, /double)

    if n_elements(new_irrad_a) ge n_elements(new_irrad_b) then begin
        fitA = poly(sda,coeffs_a)
        fitB = poly(sda,coeffs_b)
        delta_irrad = fitB - fitA
        coeff = ladfit(sda, delta_irrad, /double)
    endif else begin
        fitA = poly(sdb,coeffs_a)
        fitB = poly(sdb,coeffs_b)
        delta_irrad = fitB - fitA
        coeff = ladfit(sdb, delta_irrad, /double)
    endelse

    ;print,x[0],coeff[1]
    return, abs(coeff[1]*1d6)

end
;*******************************************************************

function adjust_esr_raypath, starttime, stoptime, version=version, $
    dburl=dburl, dbdriver=dbdriver, user=user, password=password, $
    gps=gps, missionDays=missionDays, julianDays=julianDays, $
    esrdata=esrdata, solar54=solar54, solar55=solar55, kappa=kappa, $
    order=order, solardist=solardist

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


    waveb=waveb[0:-2]
    prismposb=prismposb[0:-2]
    matchingpos=matchingpos[*,0:-2]
    minwave=min(wavea,max=maxwave)

    ;obctimes = gps2sd([0.0d, 441.5, 1570.0d, 2173d, 2455d, 2803d, 5000d])*1d6

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

    ;if n_elements(version) eq 0 then version=16
    ;if n_elements(user) eq 0 then     user='sbeland'
    ;if n_elements(password) eq 0 then password='sbeland1'
    ;if n_elements(dburl) eq 0 then    dbUrl='jdbc:sybase:Tds:lasp-db-dev:4100/SORCE_SIM_V19'
    ;if n_elements(dbdriver) eq 0 then dbDriver='com.sybase.jdbc3.jdbc.SybDriver'

    nwaves=n_elements(matchingpos[0,*])
    outdata = ptrarr(nwaves)

    if size(esrdata,/tname) ne 'STRUCT' then begin
        ; use the SimCorrectedIrradiance to avoid the aleph factor
        ;table_list=['Wavelength','SimUncorrectedIrradiance','SimConvertedDataNumbers']
        ;table_list=['Wavelength','SimCalibratedIrradiance','SimConvertedDataNumbers']
        ;esra=get_science_product(table_list,t0,t1,31,/gps, version=version,user=user,password=password,dburl=dburl,dbdriver=dbdriver)
        ;esrb=get_science_product(table_list,t0,t1,32,/gps, version=version,user=user,password=password,dburl=dburl,dbdriver=dbdriver)
        ;esrdata={esra:temporary(esra), esrb:temporary(esrb)}
        restore,'~/SORCE/data/simab_esr_uncorr_2011.sav'
    endif

    if n_elements(kappa) eq 0 then begin
        ;query_database, 'SELECT * FROM SimPrismDegKappaCal WHERE calibrationSetId=1051', kappa
        ;kappa.x=10d ^ kappa.x
        readcol,'~/SORCE/data/ESR_kappa_453_1570_2011.txt',kw,kk,format='(d,d)'
        kappa={x:kw, y:kk}
    endif

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

    if n_elements(solardist) eq 0 then begin
        print,'  getting SUNOBSERVERDISTANCECORRECTION from database ...'
        q1="SELECT microsecondsSinceGpsEpoch, sunObserverDistanceCorrection FROM SolarDistAndDopplerFixedStep2 "+$
           "WHERE instrumentModeId=50 and version=7"
        query_database, q1, solardist, info, database='SORCE_NEW'
    endif

    readcol,'~/SORCE/data/overlap_esr_v20.txt',ww,aa,format='(d,d)'
    afact={x:ww, y:aa}

    ; reform the solar exposure array to be able to interpolate between start and end of orbit only
    ; the start of orbit has no previous solar exposure since end of previous orbit
    ; the new time will have both the start and end of the orbit and we'll interpolate with the "wrong" 
    ; assumption that the exposure is evenly distributed throughout the orbit
    time54 = [solar54.t0, solar54.t1]
    solarExp54=[solar54.solar_exp_orbit*0d, solar54.solar_exp_orbit]
    s=sort(time54)
    time54=time54[s]
    solarExp54=solarExp54[s]
    solaDist54=interpol(solardist.SUNOBSERVERDISTANCECORRECTION, solardist.MICROSECONDSSINCEGPSEPOCH, time54)

    ; reform the solar exposure array to be able to interpolate between start and end of orbit only
    ; the start of orbit has no previous solar exposure since end of previous orbit
    time55 = [solar55.t0, solar55.t1]
    solarExp55=[solar55.solar_exp_orbit*0d, solar55.solar_exp_orbit]
    s=sort(time55)
    time55=time55[s]
    solarExp55=solarExp55[s]
    solaDist55=interpol(solardist.SUNOBSERVERDISTANCECORRECTION, solardist.MICROSECONDSSINCEGPSEPOCH, time55)

    solarExp54=total(solarExp54,/cum) /86400d; / (1d +(solaDist54-1d)/4d)
    solarExp55=total(solarExp55,/cum) /86400d; / (1d +(solaDist55-1d)/4d)

    ; process one wavelength at a time
    parinfo=replicate({value:0.0d, fixed:0, limited:[0,0], limits:[0d,1d], step:1d-8, tnside:2},1)

    tposA = where(esrdata.esra.(0) ge t0 and esrdata.esra.(0) le t1, count)
    if count eq 0 then begin
        print,'Error: no ESRA data within requested time range'
        return,-1
    endif

    tposB = where(esrdata.esrb.(0) ge t0 and esrdata.esrb.(0) le t1, count)
    if count eq 0 then begin
        print,'Error: no ESRB data within requested time range'
        return,-1
    endif


    a_factor=[]
    mean_wave=[]
    ; for each ESR point wavelength
    for w=0,nwaves-1 do begin
        pos=matchingpos[1,w]
        p=where(prismposb eq matchingpos[1,w],count)
        if count eq 0 then begin
            print,'Warning:  no ESR data for prismposition / wavelength: ',matchingpos[1,w]
            continue
        endif
        wavelength=waveb[p[0]]

        print,'processing '+strtrim(string(wavelength,format='(F10.1)'),2)+' nm'
        esrapos = where(esrdata.esra[tposA].prismposition eq matchingpos[0,w],counta)
        if counta eq 0 then begin
            print,'Warning:  no ESRA data for prismposition / wavelength: ',matchingpos[0,w]
            continue
        endif

        esrbpos = where(esrdata.esrb[tposB].prismposition eq matchingpos[1,w],countb)
        if countb eq 0 then begin
            print,'Warning:  no ESRB data for prismposition / wavelength: ',matchingpos[1,w]
            continue
        endif

        timetag_a  = esrdata.esra[tposA[esrapos]].(0)
        wave_a     = esrdata.esra[tposA[esrapos]].wavelengthref
        esra_irrad = esrdata.esra[tposA[esrapos]].irradiance

        timetag_b  = esrdata.esrb[tposB[esrbpos]].(0)
        wave_b     = esrdata.esrb[tposB[esrbpos]].wavelengthref
        esrb_irrad = esrdata.esrb[tposB[esrbpos]].irradiance

        solexp_a = interpol(solarExp54, time54, timetag_a)
        solexp_b = interpol(solarExp55, time55, timetag_b)
        kappa_pa = interpol(kappa.y, kappa.x, wave_a, /spline)
        kappa_pb = interpol(kappa.y, kappa.x, wave_b, /spline)

        parinfo[0].value=interpol(afact.y, afact.x, median([wave_a,wave_b]))

        ; find the best multiplation factor for the rayPath at this wavelength
        functargs = {sda:gps2sd(timetag_a/1d6), sdb:gps2sd(timetag_b/1d6), kappa_a:kappa_pa, kappa_b:kappa_pb, $
                     degcol_a:solexp_a, degcol_b:solexp_b, irrad_a:esra_irrad, irrad_b:esrb_irrad}
        coeffs = tnmin('func_getslope', functargs=functargs, bestmin=f0, status=status, $
            nfev=nfev, niter=niter, autoderivative=1, errmsg=errmsg, parinfo=parinfo)
        a_factor=[a_factor,(coeffs[0])]
        mean_wave=[mean_wave, median([wave_a,wave_b])]

        ; get the delta_irrad from linear fit after correction
        res=func_getslope(coeffs[0], df, sda=gps2sd(timetag_a/1d6), sdb=gps2sd(timetag_b/1d6), kappa_a=kappa_pa, kappa_b=kappa_pb, $
            degcol_a=solexp_a, degcol_b=solexp_b, irrad_a=esra_irrad, irrad_b=esrb_irrad, new_irrad_A=new_irradA, $
            new_irrad_B=new_irradB, delta_irrad=delta_irrad)

        if w mod 2 then begin
            resistant_mean,delta_irrad,5.0,mean,goodvec=k
            title='ESRA vs ESRB'+' @ '+strtrim(string(mean_wave[w],format='(F10.1)'),2)
            lineplot,gps2sd(timetag_a[k]/1d6),delta_irrad[k],title=title,xtitle='Mission Day',ytitle='Delta Irradiance',$
                ptitle='Differences in Irradiance between ESRA and ESRB',nsum=3
        endif

    endfor  ; loop for each ESR wavelength


    if n_elements(a_factor) eq 0 then return,-1 

    p=where(abs(a_factor) lt 50d)
    coeff=robust_poly_fit(mean_wave[p],a_factor[p], 4,/double)
    resistant_mean,a_factor[p]-poly(mean_wave[p],coeff),3.0,mean,good=k0
    sset=bspline_iterfit(mean_wave[p[k0]],a_factor[p[k0]],maxiter=0,requiren=10,bkspace=5,nord=4) 
    afit0=bspline_valu(mean_wave,sset)
    resistant_mean,a_factor[p]-afit0[p],3.0,mean,good=q0
    ww=dindgen((maxwave-minwave+10d)) + minwave -5d
    if n_elements(order) eq 0 then begin
        ; perform a bspline fit by default
        sset=bspline_iterfit(mean_wave[p[q0]],a_factor[p[q0]],maxiter=0,requiren=10,bkspace=5,nord=4)
        afit0=bspline_valu(ww,sset)
    endif else begin
        ; perform a polynomial fit
        order=order>1
        coeff=robust_poly_fit(mean_wave[p[q0]],a_factor[p[q0]], order,/double)
        afit0 = poly(ww, coeff)
        print,'order ',order,' coeff=[',coeff,']'
    endelse

    p=where(mean_wave ge 310d and mean_wave lt  1000d,count)
    if count ge 3 then cc=ladfit(mean_wave[p],a_factor[p]) else cc=ladfit(mean_wave,a_factor)
    plot_multi,mean_wave,a_factor, ww, afit0, ww,poly(ww,cc), /xst,/yst,xtitle='Wavelength (nm)', $
        ytitle='ESRA and ESRB RayPath A-Factor',$
        title='SimA vs SimB RayPath from '+strtrim(string(gps2sd(t0/1d6),format='(f0.1)'),2)+' to '+$
        strtrim(string(gps2sd(t1/1d6),format='(f0.1)'),2),charsize=1.4, psym=[-4,-3,-3],thick=[0,2,2],$
        label=['RayPath','Bspline Fit','Linear fit 310->1000nm']

    return,{wavelength:mean_wave, a_factor:a_factor, good:q0, wavefit:ww, a_fit:afit0}


end
