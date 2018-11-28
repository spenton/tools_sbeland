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
;    The table used for the irradiance is expected to have the final degradation
;    correction applied to the ESR data but no raypath or diode degradation applied to
;    the photodiode data (as in version 15 of our development database for version 19).
;
;
; REVISION HISTORY:
;   Revision: $Id: compare_esrab.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;------------------------------------------------------------------
function compare_esrab, starttime, stoptime, version=version, $
    dburl=dburl, dbdriver=dbdriver, user=user, password=password, $
    gps=gps, missionDays=missionDays, julianDays=julianDays, $
    esrdata=esrdata, percent=percent, alignobc=alignobc

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


    waveb=waveb[4:-1]
    prismposb=prismposb[4:-1]
    matchingpos=matchingpos[*,4:-1]

    obctimes = [0.0d, 441.5, 1570.0d, 2173d, 2455d, 2803d, 5000d]

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

    if n_elements(version) eq 0 then version=21
    if n_elements(user) eq 0 then     user='sbeland'
    if n_elements(password) eq 0 then password='sbeland1'
    if n_elements(dburl) eq 0 then    dbUrl='jdbc:sybase:Tds:lasp-db-dev:4100/SORCE_SIM_V19'
    if n_elements(dbdriver) eq 0 then dbDriver='com.sybase.jdbc3.jdbc.SybDriver'

    nwaves=n_elements(matchingpos[0,*])
    outdata = ptrarr(nwaves)

    if size(esrdata,/tname) ne 'STRUCT' then begin
        ; use the SimCorrectedIrradiance to avoid the aleph factor
        table_list=['Wavelength','SimCorrectedIrradiance','SimConvertedDataNumbers','SimPrismTransDegradation']
        ;table_list=['Wavelength','SimCalibratedIrradiance','SimConvertedDataNumbers']
        esra=get_science_product(table_list,t0,t1,31,/gps, version=version,user=user,password=password,dburl=dburl,dbdriver=dbdriver)
        esrb=get_science_product(table_list,t0,t1,32,/gps, version=version,user=user,password=password,dburl=dburl,dbdriver=dbdriver)
        esrdata={esra:temporary(esra), esrb:temporary(esrb)}
    endif

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
        esrapos = where(esrdata.esra.prismposition eq matchingpos[0,w],counta)
        if counta eq 0 then begin
            print,'Warning:  no ESRA data for prismposition / wavelength: ',matchingpos[0,w]
            continue
        endif

        esrbpos = where(esrdata.esrb.prismposition eq matchingpos[1,w],countb)
        if countb eq 0 then begin
            print,'Warning:  no ESRB data for prismposition / wavelength: ',matchingpos[1,w]
            continue
        endif

        ; loop through every esr value for this prismpos
        timetag_a=[]
        timetag_b=[]
        esra_irrad=[]
        esrb_irrad=[]
        delta_irrad=[]
        wave_a=[]
        print,'   looking at '+strtrim(string(count),2)+' observations ...'
        for i=0L,countb-1 do begin
            ; find the ESRA data closest in time to this ESRB
            mn=min(abs(esrdata.esrb[esrbpos[i]].(0) - esrdata.esra.(0)),pos)
            esrb_irrad = [esrb_irrad, esrdata.esrb[esrbpos[i]].irradiance]
            esra_irrad = [esra_irrad, esrdata.esra[pos].irradiance]
            delta_irrad = [delta_irrad, esra_irrad[-1]-esrb_irrad[-1]]
            timetag_b=[timetag_b,esrdata.esrb[esrbpos[i]].(0)]
            timetag_a=[timetag_a,esrdata.esra[pos].(0)]
            wave_a=[wave_a, esrdata.esra[pos].wavelengthref]
        endfor

        if keyword_set(alignobc) then begin
            pos0 = where(timetag_b ge sd2gps(453.0d)*1d6 and timetag_b le sd2gps(1570.0d)*1d6,count)
            resistant_mean,delta_irrad[pos0],3.0,avg0
            for k=0,n_elements(obctimes)-2 do begin
                align_pos = where(timetag_b ge sd2gps(obctimes[k])*1d6 and timetag_b le sd2gps(obctimes[k+1])*1d6,count)
                if count eq 0 then continue
                resistant_mean,delta_irrad[align_pos],3.0,avg1
                delta_irrad[align_pos] -= (avg1-avg0)
                resistant_mean,delta_irrad[align_pos],3.0,avg0
            endfor
        endif

        outdata[w] = ptr_new({timetag_a:timetag_a, timetag_b:timetag_b, wavelength_a:wave_a, $
            wavelength_b:esrdata.esrb[esrbpos].wavelengthref, esra_irrad:esra_irrad, esrb_irrad:esrb_irrad, $
            delta_irrad:delta_irrad})

        if not keyword_set(noplot) then begin
            title='Irradiance ESRA-ESRB  @ '+strtrim(string(wavelength,format='(F10.1)'),2)
            resistant_mean,delta_irrad,3.0,mean,goodvec=k
            if keyword_set(percent) then begin
                lineplot,gps2sd(timetag_b[k]/1d6),delta_irrad[k]/esra_irrad[k]*100d,title=title,xtitle='Mission Day',$
                    ytitle='Delta Irradiance (Percent)', ptitle='Differences in Irradiance between ESRA and ESRB',nsum=nsum
            endif else begin
                lineplot,gps2sd(timetag_b[k]/1d6),delta_irrad[k],title=title,xtitle='Mission Day',ytitle='Delta Irradiance',$
                    ptitle='Differences in Irradiance between ESRA and ESRB',nsum=nsum
            endelse

        endif

    endfor  ; loop for each ESR wavelength

    if not keyword_set(noplot) then begin
        ; plot the corresponding surface in percent changes in irradiance
        xout=gps2sd((*outdata[0]).(1)/1d6)
        yout=[]
        x=[]
        y=[]
        z=[]
        for i=0,nwaves-1 do begin
            resistant_mean,(*outdata[i]).delta_irrad,3.0,mean,goodvec=k
            x=[x,gps2sd((*outdata[i]).(1)[k]/1d6)]
            y=[y,((*outdata[i]).wavelength_b[k])]
            if keyword_set(percent) then $
                z=[z,((*outdata[i]).delta_irrad[k])/((*outdata[i]).esra_irrad[k])*100d] $
            else $
                z=[z,(*outdata[i]).delta_irrad[k]]
            yout=[yout,median((*outdata[i]).wavelength_b[k])]
        endfor
        
        grid=griddata(x,y,z,smooth=3,xout=xout,yout=yout,/grid)
        if keyword_set(percent) then ztitle='Percent Delta Irradiance' else ztitle='Delta Irradiance'
        title='Irradiance ESRA-ESRB'
        isurface,grid,xout,yout,xtitle='Mission Day',ytitle='Wavelength (nm)',ztitle=ztitle, view_title=title 
    endif



    return,outdata

end
