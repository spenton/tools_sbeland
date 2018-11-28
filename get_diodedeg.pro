;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Extract the degradation information from the photo-diodes compared
;   to the ESR table scans matching the data for the requested wavelength 
;   taken at the "same" time for ESR and diode.
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
;    the photodiode data (as in version 21 of our development database for version 19).
;
;
; REVISION HISTORY:
;   Revision: $Id: get_diodedeg.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;------------------------------------------------------------------
function get_diodedeg, starttime, stoptime, instrumentModeId, version=version, $
    gps=gps, missionDays=missionDays, julianDays=julianDays, $
    dburl=dburl, dbdriver=dbdriver, user=user, password=password, $
    inspectrum=inspectrum, esrdata=esrdata, percent=percent, alignobc=alignobc, $
    irradCol=irradCol, degColumns=degColumns, noplot=noplot

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

    obctimes = sd2gps([0.0d, 441.5, 1570.0d, 2173d, 2455d, 2803d, 5000d])*1d6

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

    ;if n_elements(version) eq 0 then version=21
    if n_elements(version) eq 0 then version=1003
    if n_elements(user) eq 0 then     user='sbeland'
    if n_elements(password) eq 0 then password='sbeland1'
    if n_elements(dburl) eq 0 then    dbUrl='jdbc:sybase:Tds:lasp-db-dev:4100/SORCE_SIM_V20'
    if n_elements(dbdriver) eq 0 then dbDriver='com.sybase.jdbc3.jdbc.SybDriver'

    ; limit the wavelength range according to the specified instrumentModeId
    if instrumentModeId eq 41 then begin
        instrument='SIMA'
        pos = where(wavea ge 310.0 and wavea lt 960.0)
        waves=wavea[pos]
        prismpos=prismposa[pos]
        esrmode=31
        nsum=3
    endif else if instrumentModeId eq 43 then begin
        instrument='SIMA'
        pos = where(wavea ge 200.0 and wavea lt 308.0)
        waves=wavea[pos]
        prismpos=prismposa[pos]
        esrmode=31
        nsum=3
    endif else if instrumentModeId eq 44 then begin
        instrument='SIMA'
        pos = where(wavea ge 846.0 and wavea lt 1670.0)
        waves=wavea[pos]
        prismpos=prismposa[pos]
        esrmode=31
        nsum=3
    endif else if instrumentModeId eq 45 then begin
        instrument='SIMB'
        pos = where(waveb ge 318.0 and waveb lt 960.0)
        waves=waveb[pos]
        prismpos=prismposb[pos]
        esrmode=32
        nsum=0
    endif else if instrumentModeId eq 47 then begin
        instrument='SIMB'
        pos = where(waveb ge 200.0 and waveb lt 308.0)
        waves=waveb[pos]
        prismpos=prismposb[pos]
        esrmode=32
        nsum=0
    endif else if instrumentModeId eq 48 then begin
        instrument='SIMB'
        pos = where(waveb ge 846.0 and waveb lt 1670.0)
        waves=waveb[pos]
        prismpos=prismposb[pos]
        esrmode=32
        nsum=0
    endif else begin
        print,'Error: wrong instrumentModeId (should be between 41 and 48)'
        return,-1
    endelse
    nwaves=n_elements(waves)
    outdata = ptrarr(nwaves)

    ; the order here is important for the query_database call later
    if size(irradCol,/tname) ne 'STRING' then begin
        irradColumns=['Wavelength','SimCorrectedIrradiance','SimConvertedDataNumbers']
    endif else begin
        irradColumns=['Wavelength',irradCol[0],'SimConvertedDataNumbers']
    endelse
    if keyword_set(degColumns) then irradColumns=[irradColumns,'SimDegradationColumn','SimPrismTransDegradation']

    query_database,/reset
    if size(esrdata,/tname) ne 'STRUCT' then begin
        esrdata=get_science_product(irradColumns,t0,t1,esrmode,/gps,$
            version=version,user=user,password=password,dburl=dburl,dbdriver=dbdriver)
        if size(esrdata,/tname) ne 'STRUCT' then begin
            print,'Error: no ESR data found for those dates: ',strtrim(string(gps2sd([t0,t1]/1d6)),2)
            return,-1
        endif
    endif

    if size(inspectrum,/tname) ne 'STRUCT' then begin
        ; extract the spectrum from SolarQuickScan24 during this time range 
        ; which are closest to the ESR data
        plans = get_sorce_plan(t0, t1, /gps,instrument=instrument,activity='SolarQuickScan24')
        plan_pos=[]
        for i=0L,n_elements(esrdata)-1L do begin
            mn=min(abs(esrdata[i].(0)/1d6 - plans.starttime),pos)
            ; only keep if diode scan is within 2 days of ESR scan
            if mn lt 2d*86400d then plan_pos=[plan_pos,pos]
        endfor
        ; only keep the unique plans
        s=sort(plan_pos)
        q=uniq(plan_pos[s])
        plans=plans[plan_pos[s[q]]]
        nspect = n_elements(plans)
        print,'   '+strtrim(string(nspect),2)+' spectrum to process ...'
        diode_data=ptrarr(nspect)
        diode_starttime=dblarr(nspect)

        ; get all of the SolarQuickScan24 data
        ; using get_science_spectra is very slow here because of the amount of pd data
        ; do a direct query instead
        for pl=0, nspect-1L do begin
            gps0 = plans[pl].starttime*1d6
            gps1 = plans[pl].stoptime*1d6
            inst = max(instrumentModeId)
            if not keyword_set(degColumns) then begin
                q1='SELECT w.microsecondsSinceGpsEpoch,w.instrumentModeId,w.version,w.wavelengthRef,ird.irradiance '+$
                   'FROM '+irradColumns[0]+' w, '+irradColumns[1]+' ird where '+$
                   'w.version='+strtrim(string(version),2)+' and w.version=ird.version and '+$
                   'w.instrumentModeId='+strtrim(string(instrumentModeId),2)+' and '+$ 
                   'w.instrumentModeId=ird.instrumentModeId and '+$
                   'w.microsecondsSinceGpsEpoch>='+strtrim(string(ulong64(gps0)),2)+' and '+$
                   'w.microsecondsSinceGpsEpoch<='+strtrim(string(ulong64(gps1)),2)+' and '+$
                   'w.microsecondsSinceGpsEpoch=ird.microsecondsSinceGpsEpoch'
            endif else begin
                q1='SELECT w.microsecondsSinceGpsEpoch,w.instrumentModeId,w.version,w.wavelengthRef,ird.irradiance,'+$
                   'dc.degradationColumn, pd.prismTransDegradation FROM '+irradColumns[0]+' w, '+$
                   irradColumns[1]+' ird, '+irradColumns[3]+' dc, '+irradColumns[4]+' pd where '+$
                   'w.version='+strtrim(string(version),2)+' and w.version=ird.version and '+$
                   'dc.version=ird.version and pd.version=ird.version and '+$
                   'w.instrumentModeId='+strtrim(string(instrumentModeId),2)+' and '+$ 
                   'w.instrumentModeId=ird.instrumentModeId and '+$
                   'dc.instrumentModeId=ird.instrumentModeId and pd.instrumentModeId=ird.instrumentModeId and '+$
                   'w.microsecondsSinceGpsEpoch>='+strtrim(string(ulong64(gps0)),2)+' and '+$
                   'w.microsecondsSinceGpsEpoch<='+strtrim(string(ulong64(gps1)),2)+' and '+$
                   'w.microsecondsSinceGpsEpoch=ird.microsecondsSinceGpsEpoch and '+$
                   'dc.microsecondsSinceGpsEpoch=ird.microsecondsSinceGpsEpoch and '+$
                   'pd.microsecondsSinceGpsEpoch=ird.microsecondsSinceGpsEpoch'
            endelse
            query_database, q1, res, user=user,password=password,dburl=dburl,dbdriver=dbdriver
            if size(res,/tname) ne 'STRUCT' then begin
                print,'Error: no spectra found for this plan ',plans[pl]
                continue
            endif
            diode_data[pl]= ptr_new(res)
            diode_starttime[pl]= res[0].(0)
        endfor
        keep = where(ptr_valid(diode_data),count)
        if count eq 0 then begin
            print,'Error:  no valid diode data found'
            return,-1
        endif
        inspectrum = {starttime:diode_starttime[keep], diode:diode_data[keep]}
    endif

    ; for each ESR point wavelength
    for w=0,n_elements(waves)-1 do begin
        print,'processing '+strtrim(string(waves[w],format='(F10.1)'),2)+' nm'
        esrpos = where(esrdata.prismposition eq prismpos[w],count)
        if count eq 0 then begin
            print,'Warning:  no ESR data for prismposition / wavelength: ',prismpos[w],waves[w]
            continue
        endif

        ; loop through every esr value for this prismpos
        timetag=[]
        pd_irrad=[]
        esr_irrad=[]
        delta_irrad=[]
        timetag_pd=[]
        esr_degCol=[]
        pd_degCol=[]
        esr_prismDeg=[]
        pd_prismDeg=[]
        print,'   looking at '+strtrim(string(count),2)+' observations ...'
        for i=0L,count-1 do begin
            ; find the diode scan the closest in time to this esr and interpolate the irradiance
            mn=min(abs(esrdata[esrpos[i]].(0) - inspectrum.starttime),pos)
            if not ptr_valid(inspectrum.diode[pos]) then continue
            irrad = interpol((*inspectrum.diode[pos]).irradiance, (*inspectrum.diode[pos]).wavelengthref, $
                esrdata[esrpos[i]].wavelengthref, /spline)
            pd_irrad = [pd_irrad, irrad]
            d_irrad = esrdata[esrpos[i]].irradiance - irrad
            delta_irrad = [delta_irrad, d_irrad]
            timetag=[timetag,esrdata[esrpos[i]].(0)]
            mn=min(abs(esrdata[esrpos[i]].wavelengthref - (*inspectrum.diode[pos]).wavelengthref),tpos)
            t_pd = (*inspectrum.diode[pos])[tpos].(0)
            timetag_pd = [timetag_pd, t_pd]
            if keyword_set(degColumns) then begin
                pd_degCol=[pd_degCol, (*inspectrum.diode[pos])[tpos].degradationColumn]
                esr_degCol=[esr_degCol, esrdata[esrpos[i]].degradationColumn]
                pd_prismDeg=[pd_prismDeg, (*inspectrum.diode[pos])[tpos].prismTransDegradation]
                esr_prismDeg=[esr_prismDeg, esrdata[esrpos[i]].prismTransDegradation]
            endif
        endfor
        if keyword_set(degColumns) then begin
            outdata[w] = ptr_new({timetag:timetag, wavelength:esrdata[esrpos].wavelengthref, timetag_pd:timetag_pd, $
                esr_irrad:esrdata[esrpos].irradiance, pd_irrad:pd_irrad, delta_irrad:delta_irrad, pd_degcol:pd_degcol,$
                esr_degcol:esr_degcol, pd_prismDeg:pd_prismDeg, esr_prismDeg:esr_prismDeg})
        endif else begin
            outdata[w] = ptr_new({timetag:timetag, wavelength:esrdata[esrpos].wavelengthref, timetag_pd:timetag_pd, $
                esr_irrad:esrdata[esrpos].irradiance, pd_irrad:pd_irrad, delta_irrad:delta_irrad})
        endelse

        if keyword_set(alignobc) then begin
            align_irrad, timetag, delta_irrad, /gps
;            for i=0,n_elements(obctimes)-3 do begin
;                pleft = where(timetag ge obctimes[i] and timetag lt obctimes[i+1],count)
;                if count eq 0 then continue
;                pright = where(timetag ge obctimes[i+1] and timetag lt obctimes[i+2],count)
;                if count eq 0 then continue
;                p0=n_elements(pleft)-21 > 0
;                cleft=robust_poly_fit(timetag[pleft[p0:-2]], delta_irrad[pleft[p0:-2]], 1, /double)
;                p1=n_elements(pright)-1 <20 
;                cright=robust_poly_fit(timetag[pright[1:p1]], delta_irrad[pright[1:p1]], 1, /double)
;                v0=poly(timetag[pright[0]],cleft)
;                v1=poly(timetag[pright[0]],cright)
;                delta_irrad[pright] += (v0-v1)
;            endfor
            (*outdata[w]).delta_irrad = delta_irrad
        endif

        if not keyword_set(noplot) then begin
            title='ESR vs '+strtrim(string(instrumentModeId,format='(I0)'),2)+$
                ' @ '+strtrim(string(waves[w],format='(F10.1)'),2)
            resistant_mean,delta_irrad,5.0,mean,goodvec=k
            if keyword_set(percent) then begin
                lineplot,gps2sd(timetag[k]/1d6),delta_irrad[k]/esrdata[esrpos[k]].irradiance*100d,title=title,xtitle='Mission Day',$
                    ytitle='Delta Irradiance (Percent)', ptitle='Differences in Irradiance between ESR and Diode ('+$
                    strtrim(string(instrumentModeId,format='(I0)'),2)+')', nsum=nsum
            endif else begin
                lineplot,gps2sd(timetag[k]/1d6),delta_irrad[k],title=title,xtitle='Mission Day',ytitle='Delta Irradiance',$
                    ptitle='Differences in Irradiance between ESR and Diode ('+strtrim(string(instrumentModeId,format='(I0)'),2)+')',nsum=nsum
            endelse

        endif

    endfor  ; loop for each ESR wavelength

    if not keyword_set(noplot) then begin
        ; plot the corresponding surface in percent changes in irradiance
        xout=gps2sd((*outdata[0]).(0)/1d6)
        yout=[]
        x=[]
        y=[]
        z=[]
        for i=0,nwaves-1 do begin
            resistant_mean,(*outdata[i]).delta_irrad,3.0,mean,goodvec=k
            x=[x,gps2sd((*outdata[i]).(0)[k]/1d6)]
            y=[y,(*outdata[i]).wavelength[k]]
            if keyword_set(percent) then $
                z=[z,((*outdata[i]).delta_irrad[k])/((*outdata[i]).esr_irrad[k])*100d] $
            else $
                z=[z,(*outdata[i]).delta_irrad[k]] 
            yout=[yout,median((*outdata[i]).wavelength[k])]
        endfor
        grid=griddata(x,y,z,smooth=3,xout=xout,yout=yout,/grid)
        if keyword_set(percent) then ztitle='Percent Delta Irradiance' else ztitle='Delta Irradiance'
        title='ESR vs '+strtrim(string(instrumentModeId,format='(I0)'),2)
        isurface,grid,xout,yout,xtitle='Mission Day',ytitle='Wavelength (nm)',ztitle=ztitle, view_title=title 
    endif



    return,outdata

end
