;+
; NAME:   GET_ESRKAPPA
;
; PURPOSE: 
;    This routine uses the ESRA and ESRB spectrum for the
;    first part of the mission to determine a Kappa vs wavelength
;
; CALLING SEQUENCE:
;    result=get_esrkappa(inspectrum=inspectrum, solar54=solar54, solar55=solar55)
;
; OPTIONAL INPUT PARAMETERS:
;    none
;
; OPTIONAL INPUT KEYWORDS:
;    none 
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
;   Revision: $Id: get_esrkappa.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;-
;*****************************************************************************
FUNCTION GET_ESRKAPPA, timerange=timerange, inspectrum=inspectesr, solar54=solar54, solar55=solar55, noplot=noplot, $
    ffunct=ffunct, afact=afact, verbose=verbose, version=version, rfact=rfact

    if n_elements(timerange) eq 0 then timerange=[453d,1570d]
    t0=min(timerange,max=t1)
    ;t0=453d
    ;t1=1570d
    ;t0=1570d
    ;t1=2173d
    ; wavea = [261.2d, 265.2, 279.5, 280.5, 283.4, 285.3, 288.2, 293.7, $
    ; 304.6, 319.2, 333.6, 341.8, 348.5, 355, 369.2, 375.2, 384.5, 394.6, $
    ; 397.8, 406.7, 430.9, 470.8, 480.2, 488.2, 519.6, 566.8, 593.7, 659.1, $
    ; 710.2, 758.1, 805.8, 859.4, 867.5, 892.7, 901.2, 964.8, 1013.3, 1064, $
    ; 1116.4, 1170.2, 1205.8, 1221.8, 1238.6, 1272, 1289, 1306, 1323, 1340, $
    ; 1356.4, 1373, 1390, 1406.5, 1423, 1440, 1456, 1472.8, 1489, 1505, 1521.5, $
    ; 1537.7, 1554, 1569.6, 1585.4, 1601, 1617, 1632.3, 1648, 1663.5, 1678.5, 1694, $
    ; 1709, 1724, 1738.8, 1753.6, 1768.4, 1783.1, 1797.7, 1812, 1826.5, 1841, 1855, $
    ; 1869, 1883, 1897, 1911, 1924.5, 1938, 1951.5, 1965, 1978.6, 1992, 2005.2, 2018.4, $
    ; 2031.2, 2044.2, 2070, 2083, 2095.5, 2108.2, 2120.6, 2133, 2145.5, 2158, 2170, 2183.4, $
    ; 2194.5, 2206.6, 2218.5, 2230.5, 2242.2, 2254, 2265.6, 2277.4, 2289, 2300.5, 2312, $
    ; 2323.4, 2334.7, 2346, 2357.2, 2368.4, 2379.4, 2390.6, 2401.5, 2412.4, 2423.2, 2434.2, $
    ; 2445, 2455.7, 2466.2, 2477, 2487.4, 2498]

    ;1549.6, 1591, $
    ;1621.7, 1692, 1741, 1817.2, 1881.9, 1909, 1944.7, 1988.5, 1992.1, $
    ;2277.5, 2300.56, 2323.5, 2357.3, 2390.6, 2412.5]

    ;wavea = [288, 293.5, 304.5, 319, 333, 342, 355, 369, 375, 385, 395, $
    ;        407, 431, 471, 481, 489, 519, 568, 593, 662, 713, 761, 808]

    ; these correpond to the ESRMode (complete table scan)
    ;wavea = [261.2, 265.2, 279.5, 280.5, 283.4, 285.3, 288.2, 293.7, 304.6, 319.2, 333.7, $
    ;         341.8, 348.5, 354.9, 369.2, 375.2, 384.5, 394.6, 397.8, 406.7, 430.9, 470.8, $
    ;         480.2, 488.2, 519.6, 566.8, 593.7, 659.1, 710.2, 758.1, 805.9, 859.4, 867.5, $
    ;         892.7, 901.3, 964.8, 1013.3, 1064.0, 1116.4, 1170.2, 1213.8, 1279.7, 1356.6, $
    ;         1475.7, 1497.0, 1549.6, 1591.1, 1621.7, 1692.0, 1740.9, 1817.2, 1881.9, 1909.1, $
    ;         1944.7, 1988.5, 2502.5]

    ;waveb = [278.4, 279.4, 282.3, 284.2, 287.0, 292.4, 303.3, 317.7, 331.9, 339.9, 346.5, $
    ;         352.9, 366.9, 372.8, 382.0, 391.9, 395.0, 403.8, 427.4, 466.4, 475.6, 483.4, $
    ;         513.9, 559.8, 585.8, 649.2, 698.6, 744.9, 791.1, 843.2, 851.1, 875.6, 884.0, $
    ;         946.1, 993.8, 1043.8, 1096.0, 1149.5, 1193.1, 1259.2, 1336.3, 1391.3, 1456.5, $
    ;         1477.8, 1531.0, 1573.0, 1604.0, 1675.0, 1724.5, 1801.6, 1867.0, 1894.4, 1930.6, $
    ;         1974.8, 2493.6]

    ; these are a subset of ESRA table scan that matches well with ESR table scan wavelengths
    ; [wavelength, prismposition]

    wavea = [279.5, 283.4, $
             ;285.3, 288.2, $
             293.7, 304.6, 319.2, 333.7, 341.8, 348.5, 354.9, $
             369.2, 375.2, 384.5, 394.6, $
             406.7, $
             430.9, 470.8, 480.2, 488.2, 519.6, 566.8, $
             593.7, 659.1, 710.2, 758.1, 805.9, 859.4, 867.5, 892.7, 901.3, 964.8, 1013.3, $
             1064.0, 1116.4, 1170.2, 1213.8, 1356.6, 1475.7, 1497.0, 1549.6, 1591.1, 1621.7, $
             1692.0, 1740.9, 1817.2, 1881.9, 1909.1, 1944.7, 1988.5, 2502.5]

    waveposa = [[279.5, 10000], $
                [283.4, 10916], $
;                [285.3, 11334], $
;                [288.2, 11950], $
                [293.6, 13075], $
                [304.6, 15100], $
                [319.2, 17425], $
                [333.7, 19375], $
                [341.8, 20350], $
                [348.5, 21100], $
                [354.9, 21775], $
                [369.2, 23125], $
                [375.2, 23650], $
                [384.5, 24400], $
                [394.6, 25150], $
                [406.7, 25975], $
                [430.9, 27400], $
                [470.8, 29275], $
                [480.2, 29650], $
                [488.2, 29950], $
                [519.6, 31000], $
                [566.9, 32275], $
                [593.7, 32875], $
                [659.1, 34075], $
                [710.2, 34825], $
                [758.1, 35425], $
                [805.9, 35950], $
                [859.4, 36475], $
                [867.5, 36550], $
                [892.7, 36775], $
                [901.3, 36850], $
                [964.8, 37375], $
                [1013.3, 37750], $
                [1064.0, 38125], $
                [1116.4, 38500], $
                [1170.2, 38875], $
                [1213.8, 39175], $
                [1356.6, 40150], $
                [1475.7, 40975], $
                [1497.0, 41125], $
                [1549.6, 41500], $
                [1591.1, 41800], $
                [1621.7, 42025], $
                [1692.0, 42550], $
                [1740.9, 42925], $
                [1817.2, 43525], $
                [1881.9, 44050], $
                [1909.1, 44275], $
                [1944.7, 44575], $
                [1988.5, 44950], $
                [2502.5, 49975]]

    waveb = [279.4, 282.3, $
             ;284.2, 287.0, $
             292.4, 303.3, 317.7, 331.9, 339.9, 346.5, 352.9, $
             366.9, 372.8, 382.0, 395.0, $
             403.8, $
             427.4, 466.4, 475.6, 483.4, 513.9, 559.8, $
             585.8, 649.2, 698.6, 744.9, 791.1, 843.2, 875.6, 884.0, 946.1, 993.8, 1043.8, $
             1096.0, 1149.5, 1193.1, 1259.2, 1336.3, 1456.5, 1477.8, 1531.0, 1573.0, 1604.0, $
             1675.0, 1724.5, 1801.6, 1867.0, 1894.4, 1930.6, 1974.8, 2493.6]

    waveposb = [[279.4, 52156], $
                [282.3, 51473], $
;                [284.2, 51056], $
;                [287.0, 50441], $
                [292.4, 49318], $
                [303.3, 47296], $
                [317.7, 44975], $
                [331.9, 43029], $
                [339.9, 42055], $
                [346.5, 41307], $
                [352.9, 40633], $
                [366.9, 39286], $
                [372.8, 38762], $
                [382.0, 38014], $
                [395.0, 37041], $
                [403.8, 36442], $
                [427.4, 35020], $
                [466.4, 33149], $
                [475.6, 32775], $
                [483.4, 32475], $
                [513.9, 31428], $
                [559.8, 30156], $
                [585.8, 29557], $
                [649.2, 28360], $
                [698.6, 27612], $
                [744.9, 27014], $
                [791.1, 26490], $
                [843.2, 25966], $
                [875.6, 25667], $
                [884.0, 25592], $
                [946.1, 25069], $
                [993.8, 24695], $
                [1043.8, 24321], $
                [1096.0, 23946], $
                [1149.5, 23572], $
                [1193.1, 23273], $
                [1259.2, 22824], $
                [1336.3, 22301], $
                [1456.5, 21478], $
                [1477.8, 21329], $
                [1531.0, 20955], $
                [1573.0, 20655], $
                [1604.0, 20431], $
                [1675.0, 19907], $
                [1724.5, 19533], $
                [1801.6, 18935], $
                [1867.0, 18412], $
                [1894.4, 18188], $
                [1930.6, 17888], $
                [1974.8, 17514], $
                [2493.6, 12504]]

    if n_elements(afact) eq 0 then begin
        ;readcol,'~/SORCE/data/overlap_esr_fit.txt',esrpath_w, tmp, esrpath_a, format='(d,d,d)'
        readcol,'~/SORCE/data/overlap_esr_v20.txt',esrpath_w, esrpath_a, format='(d,d)'
        afact = {wavelength:esrpath_w, a:esrpath_a}; +0.15d}
    endif

    if n_elements(solar54) eq 0 then begin
        restore,'~/SORCE/data/solarexp_54_55_mod_407nm.sav'
    endif

    esrkappa=dblarr(n_elements(wavea))
    ; loop over every wavelengths and find a kappa for that wavelength (this will take a while)
    ; we need to pick a version of the processed data that has no degradation applied (v15)
    for i=0,n_elements(wavea)-1 do begin 
        if keyword_set(verbose) then print,wavea[i]
        ;res=tnmin_esrdeg_all(t0,t1,prismposa=waveposa[1,i], prismposb=waveposb[1,i], /mission,coeffs=coeff,inspectrum=inspectesr, ffunct=ffunct, $
        res=min_esrdeg_all(t0,t1,prismposa=waveposa[1,i], prismposb=waveposb[1,i], /mission,coeffs=coeff,inspectrum=inspectesr, ffunct=ffunct, $
            solar54=solar54, solar55=solar55, afact=afact, version=version, verbose=verbose, rfact=rfact, /submean) ;, /oneau) 
            ; using version 1003 of SORCE_SIM_V20 SimUnccorectedIrradiance is already 1AU corrected
        esrkappa[i]=coeff 
    endfor

    ; cleanup the kappa (remove negative values)
    esrkappa=abs(esrkappa)>1.0d-12
    ; zero out the degradation for wavelengths longer than 1300nm
    kk=where(esrkappa lt 0.50d,count)
    esrkappa=esrkappa[kk]
    wavea=wavea[kk]
    ;esrkappa[kk] = 1d-12

    ;k = where(wavea ge 290d)
    k = where(wavea ge 200d)
    ; fit a bspline
    ww=[2100d, 2200d, 2300d, 2400d]
    p0=where(wavea lt 2100d)
    p1=where(wavea gt 2400d)
    tmp=interpol(esrkappa[[p0[-1],p1[0]]],wavea[[p0[-1],p1[0]]],ww)

    ;sset=bspline_iterfit(wavea[k],esrkappa[k],maxiter=10,requiren=10,bkspace=5,nord=5)
    sset=bspline_iterfit([wavea[p0],ww,wavea[p1]],[esrkappa[p0],tmp,esrkappa[p1]],requiren=10,bkspace=5,nord=6)
    ; define a fine grid wavelength
    ww=dindgen(2275) +225d
    kappa_fit=bspline_valu(ww,sset)

    guess=[0.01d, -0.007d, 1.0d-5]
    k0=where(wavea ge 290d and wavea le 900d)
    ;dd=exponential_fit(wavea[k0],esrkappa[k0],guess=guess, fita=[1,1,1])
    dd=exponential_fit(wavea[k],esrkappa[k],guess=guess, fita=[1,1,1])
    yfit = (dd[0]*exp(dd[1]*wavea)+dd[2])
    ;resistant_mean,esrkappa-yfit,4.0,mean,good=keep0
    ;if n_elements(keep0) lt n_elements(wavea) and n_elements(keep0) gt n_elements(wavea)*0.5 then begin
    ;    dd=exponential_fit(wavea[keep0],esrkappa[keep0],guess=guess, fita=[1,1,1])
    ;    sset=bspline_iterfit(wavea[keep0],esrkappa[keep0],maxiter=10,requiren=10,bkspace=5,nord=4)
    ;    kappa_fit=bspline_valu(ww,sset)
    ;endif

    ; plot the data
    if NOT keyword_set(noplot) then begin
        yfit = (dd[0]*exp(dd[1]*ww)+dd[2])
        label3='Exp Fit '+strtrim(string(dd[0]),2)+', '+strtrim(string(dd[1]),2)+', '+strtrim(string(dd[2]),2)
        title="ESRA, B Degradation model Kappa ["+strtrim(string(timerange,format='(F0.2,"->",F0.2,"]")'),2)
        plot_multi,wavea,esrkappa,ww,kappa_fit,ww,yfit,/xst,/yst,title=title, $
            xtitle="Wavelength (nm)",ytitle="Kappa",charsize=1.4,label=["ESR_Kappa","BSpline Fit",label3],$
            psym=[-4,-3,-3],thick=[1.0,2.0,2.0]
    endif

    return, {wavelength:wavea, kappa:esrkappa, x:ww, y:kappa_fit, exp_coeff:dd}

END

