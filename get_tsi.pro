;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Read each SolarQuickScan24 spectra within the requested timerange 
;   and INTEGRATE the Solar Spectra Irradiance over the desired wavelength 
;   range (to compare with the TSI time series data).
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
;
; REVISION HISTORY:
;   Revision: $Id: get_tsi.pro,v 1.2 2018/11/26 16:16:02 spenton Exp $
;------------------------------------------------------------------
function get_tsi, starttime, stoptime, sima=sima, simb=simb, wrange=wrange, $
    gps=gps, missionDays=missionDays, julianDays=julianDays, dailyavg=dailyavg, $
    dburl=dburl, dbdriver=dbdriver, user=user, password=password, align=align, $
    version=version, uv_version=uv_version, vis_version=vis_version, ir_version=ir_version, $
    level3=level3, alldata=alldata, refwaves=refwaves, mincount=mincount, noplot=noplot

    if keyword_set(gps) then begin
       ; user specified time in mission (sorce) days
       t0 = gps2sd(startTime/1.d6)
       t1 = gps2sd(stopTime/1.d6)
    endif else if keyword_set(julianDays) then begin
       ; user specified time in julian days
       t0 = jd2sd(startTime)
       t1 = jd2sd(stopTime)
    endif else begin
       ; user specified timetags in gps microseconds
       t0 = startTime
       t1 = stopTime
    endelse

    if n_elements(version) eq 0 then version=2204
    if n_elements(uv_version) eq 0 then uv_version=version
    if n_elements(vis_version) eq 0 then vis_version=version
    if n_elements(ir_version) eq 0 then ir_version=version
    if keyword_set(sima) then simb=0
    if keyword_set(simb) then sima=0
    if n_elements(sima) eq 0 or n_elements(simb) eq 0 then begin
        ; get sima by default
        sima=1
        simb=0
    endif

    if n_elements(wrange) ne 2 then wrange=[240.0d, 2400.0d]
    if n_elements(mincount) eq 0 then mincount=10

    if keyword_set(level3) then begin
        if size(alldata,/tname) ne 'STRUCT' then alldata = {uv:ptr_new(), vis:ptr_new(), ir:ptr_new(), esr:ptr_new()}
        if n_elements(wrange) ne 2 then wrange=[240.0d, 2400.0d]

        ; for level3 data the queries are more straightforward
        if n_elements(user) eq 0 then     user='sorce'
        if n_elements(password) eq 0 then password='sorcedb'
        if n_elements(dburl) eq 0 then    dbUrl='jdbc:sybase:Tds:lasp-db1:4100/SORCE'
        if n_elements(dbdriver) eq 0 then dbDriver='com.sybase.jdbc3.jdbc.SybDriver'

        if wrange[0] lt 309.d and not ptr_valid(alldata.uv) then begin
            print,'Getting UV data ...'
            iMode=43
            w1=max([wrange[0],240d])
            w2=min([wrange[1],308d])
            sql = 'select * from Level3SolarSpectra s, Level3SolarSpectraData d ' + $
                ' where s.spectraId=d.spectraId ' + $
                ' and s.version=' + strtrim(version,2) + $
                ' and s.instrumentModeId= ' + strtrim(iMode,2) + $
                ' and s.timeSpanInHours=24' + $
                ' and s.nominalGpsTimetag > ' + strtrim(string(ulong64(sd2gps(t0)*1d6)),2) + $
                ' and s.nominalGpsTimetag < ' + strtrim(string(ulong64(sd2gps(t1)*1d6)),2) + $
                ' and d.minWavelength >= ' + string(w1) + $
                ' and d.minWavelength <= ' + string(w2)
            query_database, sql, uv_data, nrows, database=database, dbUrl=dbUrl, dbDriver=dbDriver, user=user,password=password
            alldata.uv = ptr_new(uv_data,/no_copy)
        endif

        if wrange[0] lt 950.001  and wrange[1] ge 310.0 and not ptr_valid(alldata.vis) then begin
            print,'Getting VIS data ...'
            iMode=41
            w1=max([wrange[0],310d])
            w2=min([wrange[1],950d])
            sql = 'select * from Level3SolarSpectra s, Level3SolarSpectraData d ' + $
                ' where s.spectraId=d.spectraId ' + $
                ' and s.version=' + strtrim(version,2) + $
                ' and s.instrumentModeId= ' + strtrim(iMode,2) + $
                ' and s.timeSpanInHours=24' + $
                ' and s.nominalGpsTimetag > ' + strtrim(string(ulong64(sd2gps(t0)*1d6)),2) + $
                ' and s.nominalGpsTimetag < ' + strtrim(string(ulong64(sd2gps(t1)*1d6)),2) + $
                ' and d.minWavelength >= ' + string(w1) + $
                ' and d.minWavelength <= ' + string(w2)
            query_database, sql, vis_data, nrows, database=database, dbUrl=dbUrl, dbDriver=dbDriver, user=user,password=password
            alldata.vis = ptr_new(vis_data,/no_copy)
        endif

        if wrange[0] le 1600.001  and wrange[1] gt 950.001 and not ptr_valid(alldata.ir) then begin
            print,'Getting IR data ...'
            iMode=44
            w1=max([wrange[0],950d])
            w2=min([wrange[1],1600d])
            sql = 'select * from Level3SolarSpectra s, Level3SolarSpectraData d ' + $
                ' where s.spectraId=d.spectraId ' + $
                ' and s.version=' + strtrim(version,2) + $
                ' and s.instrumentModeId= ' + strtrim(iMode,2) + $
                ' and s.timeSpanInHours=24' + $
                ' and s.nominalGpsTimetag > ' + strtrim(string(ulong64(sd2gps(t0)*1d6)),2) + $
                ' and s.nominalGpsTimetag < ' + strtrim(string(ulong64(sd2gps(t1)*1d6)),2) + $
                ' and d.minWavelength > ' + string(w1) + $
                ' and d.minWavelength <= ' + string(w2)
            query_database, sql, ir_data, nrows, database=database, dbUrl=dbUrl, dbDriver=dbDriver, user=user,password=password
            alldata.ir = ptr_new(ir_data,/no_copy)
        endif

        if wrange[1] gt 1600.0 and not ptr_valid(alldata.esr) then begin
            print,'Getting ESR data ...'
            iMode=31
            w1=max([wrange[0],1600d])
            w2=min([wrange[1],2400d])
            sql = 'select * from Level3SolarSpectra s, Level3SolarSpectraData d ' + $
                ' where s.spectraId=d.spectraId ' + $
                ' and s.version=' + strtrim(version,2) + $
                ' and s.instrumentModeId= ' + strtrim(iMode,2) + $
                ' and s.timeSpanInHours=24' + $
                ' and s.nominalGpsTimetag > ' + strtrim(string(ulong64(sd2gps(t0)*1d6)),2) + $
                ' and s.nominalGpsTimetag < ' + strtrim(string(ulong64(sd2gps(t1)*1d6)),2) + $
                ' and d.minWavelength >= ' + string(w1) + $
                ' and d.minWavelength <= ' + string(w2)
            query_database, sql, esr_data, nrows, database=database, dbUrl=dbUrl, dbDriver=dbDriver, user=user,password=password
            alldata.esr = ptr_new(esr_data,/no_copy)
        endif

        ; if align is specified, re-aligned the irradiance at OBC one wavelength at a time
        if keyword_set(align) then begin
            if ptr_valid(alldata.uv) gt 0 then begin
                s=sort((*alldata.uv).minwavelength)
                q=uniq((*alldata.uv)[s].minwavelength)
                p=where((*alldata.uv)[s[q]].minwavelength gt 0.0)
                waves = (*alldata.uv)[s[q[p]]].minwavelength
                for i=0L,n_elements(waves)-1L do begin
                    p=where((*alldata.uv).minwavelength eq waves[i],count)
                    if count eq 0 then continue
                    irrad = (*alldata.uv)[p].irradiance
                    timetag = (*alldata.uv)[p].NOMINALGPSTIMETAG
                    align_irrad, timetag, irrad ,/gps
                    (*alldata.uv)[p].irradiance = irrad
                endfor
            endif
            if ptr_valid(alldata.vis) gt 0 then begin
                s=sort((*alldata.vis).minwavelength)
                q=uniq((*alldata.vis)[s].minwavelength)
                p=where((*alldata.vis)[s[q]].minwavelength gt 0.0)
                waves = (*alldata.vis)[s[q[p]]].minwavelength
                for i=0L,n_elements(waves)-1L do begin
                    p=where((*alldata.vis).minwavelength eq waves[i],count)
                    if count eq 0 then continue
                    irrad = (*alldata.vis)[p].irradiance
                    timetag = (*alldata.vis)[p].NOMINALGPSTIMETAG
                    align_irrad, timetag, irrad ,/gps
                    (*alldata.vis)[p].irradiance = irrad
                endfor
            endif
            if ptr_valid(alldata.ir) gt 0 then begin
                s=sort((*alldata.ir).minwavelength)
                q=uniq((*alldata.ir)[s].minwavelength)
                p=where((*alldata.ir)[s[q]].minwavelength gt 0.0)
                waves = (*alldata.ir)[s[q[p]]].minwavelength
                for i=0L,n_elements(waves)-1L do begin
                    p=where((*alldata.ir).minwavelength eq waves[i],count)
                    if count eq 0 then continue
                    irrad = (*alldata.ir)[p].irradiance
                    timetag = (*alldata.ir)[p].NOMINALGPSTIMETAG
                    align_irrad, timetag, irrad ,/gps
                    (*alldata.ir)[p].irradiance = irrad
                endfor
            endif
            if ptr_valid(alldata.esr) gt 0 then begin
                s=sort((*alldata.esr).minwavelength)
                q=uniq((*alldata.esr)[s].minwavelength)
                p=where((*alldata.esr)[s[q]].minwavelength gt 0.0)
                waves = (*alldata.esr)[s[q[p]]].minwavelength
                for i=0L,n_elements(waves)-1L do begin
                    p=where((*alldata.esr).minwavelength eq waves[i],count)
                    if count eq 0 then continue
                    irrad = (*alldata.esr)[p].irradiance
                    timetag = (*alldata.esr)[p].NOMINALGPSTIMETAG
                    align_irrad, timetag, irrad ,/gps
                    (*alldata.esr)[p].irradiance = irrad
                endfor
            endif
        endif


        ; integrate the data for each day
        ndays = ceil(t1) - floor(t0) + 1d
        timestamp=dblarr(ndays)
        tsi=dblarr(ndays)

        gt0=sd2gps(floor(t0))*1d6

        print,'Doing daily integration ... '
        for i=0L, ndays-1L do begin
            print,i+1,double(i+1)/double(ndays+1)*100d,format='("   ",i0,"  (",i0," % )",$,%"\r")'
            ttime=[]
            waves=[]
            irrad=[]
            if ptr_valid(alldata.esr) gt 0 then begin
                pesr=where((*alldata.esr).NOMINALGPSTIMETAG ge gt0+i*86400d6 and $
                    (*alldata.esr).NOMINALGPSTIMETAG lt gt0+(i+1d)*86400d6 and $
                    (*alldata.esr).minwavelength ge wrange[0] and $
                    (*alldata.esr).minwavelength le wrange[1],count)
                if count eq 0 then continue
                ttime=[ttime,(*alldata.esr)[pesr].NOMINALGPSTIMETAG]
                waves=[waves,(*alldata.esr)[pesr].minwavelength]
                irrad=[irrad,(*alldata.esr)[pesr].irradiance]
            endif
            if ptr_valid(alldata.uv) gt 0 then begin
                puv=where((*alldata.uv).NOMINALGPSTIMETAG ge gt0+i*86400d6 and $
                    (*alldata.uv).NOMINALGPSTIMETAG lt gt0+(i+1d)*86400d6 and $
                    (*alldata.uv).minwavelength ge wrange[0] and $
                    (*alldata.uv).minwavelength le wrange[1],count)
                if count eq 0 then continue
                ttime=[ttime,(*alldata.uv)[puv].NOMINALGPSTIMETAG]
                waves=[waves,(*alldata.uv)[puv].minwavelength]
                irrad=[irrad,(*alldata.uv)[puv].irradiance]
            endif
            if ptr_valid(alldata.vis) gt 0 then begin
                pvis=where((*alldata.vis).NOMINALGPSTIMETAG ge gt0+i*86400d6 and $
                    (*alldata.vis).NOMINALGPSTIMETAG lt gt0+(i+1d)*86400d6 and $
                    (*alldata.vis).minwavelength ge wrange[0] and $
                    (*alldata.vis).minwavelength le wrange[1],count)
                if count eq 0 then continue
                ttime=[ttime,(*alldata.vis)[pvis].NOMINALGPSTIMETAG]
                waves=[waves,(*alldata.vis)[pvis].minwavelength]
                irrad=[irrad,(*alldata.vis)[pvis].irradiance]
            endif
            if ptr_valid(alldata.ir) gt 0 then begin
                pir=where((*alldata.ir).NOMINALGPSTIMETAG ge gt0+i*86400d6 and $
                    (*alldata.ir).NOMINALGPSTIMETAG lt gt0+(i+1d)*86400d6 and $
                    (*alldata.ir).minwavelength ge wrange[0] and $
                    (*alldata.ir).minwavelength le wrange[1],count)
                if count eq 0 then continue
                ttime=[ttime,(*alldata.ir)[pir].NOMINALGPSTIMETAG]
                waves=[waves,(*alldata.ir)[pir].minwavelength]
                irrad=[irrad,(*alldata.ir)[pir].irradiance]
            endif
            timestamp[i]=median(ttime)
            s=sort(waves)
            q=uniq(waves[s])
            tsi[i]=int_tabulated(waves[s[q]], irrad[s[q]], /double, /sort)
        endfor

        print,''
        p =where(timestamp gt 0d and tsi gt 0d,count)
        if count eq 0 then return,-1
        timestamp=timestamp[p]
        tsi=tsi[p]
        if keyword_set(align) then align_irrad,timestamp, tsi, /gps
        if NOT keyword_set(noplot) then $
            lineplot,gps2sd(timestamp/1d6), tsi, xtitle='Mission Day', ytitle='Integrated Irradiance', $
                ptitle='SimA Level3 Integrated Irradiance from '+strcompress(string(wrange[0],format='(F0.1)')+$
                ', '+string(wrange[1],format='(F0.1)')), charsize=1.4,psym=-4,title='TSI version='+strtrim(string(version),2)

        return,{timestamp:timestamp, tsi:tsi} 

    endif


    if size(alldata,/tname) eq 'STRUCT' then begin
        valid_tag = where(strpos(tag_names(alldata),'SPECT20') ge 0)
        if valid_tag[0] ne -1 then nplans=n_elements(alldata.spect20) else nplans=n_elements(alldata)
    endif else begin
        if n_elements(user) eq 0 then     user='sbeland'
        if n_elements(password) eq 0 then password='sbeland1'
        if n_elements(dburl) eq 0 then    dbUrl='jdbc:sybase:Tds:lasp-db-dev:4100/SORCE_SIM_V22'
        if n_elements(dbdriver) eq 0 then dbDriver='com.sybase.jdbc3.jdbc.SybDriver'
        if n_elements(wrange) ne 2 then wrange=[240.0d, 1600.0d]


        ;jstmt = fjava_get_jdbc_statement(user=user, password=password, dbdriver=dbdriver, dburl=dburl) 
        ;oJavaDbExchange = OBJ_NEW("oJava_DbExchange", jstmt)


; for after day 3570, the planning database is not being populated.
; we could rely on the prism position, shutter, half_cycle and gain_mode_hc
; to identify specific activities
;
; q1="select sampleVtcw 'timetag',prismPosition 'dn', prismPosition 'eu' from SimAScienceSamples where "+$
;    "prismPosition=4038 and sampleVtcw>=727401613000000 and sampleVtcw<1035849616000000 order by sampleVtcw"
; query_database,q1,prismPos,nrows,database="SORCE_L1S"
;
; q1="select sampleVtcw 'timetag',prismPosition 'dn', shutterPosition 'eu' from SimAScienceSamples where "+$
;    "prismPosition=4038 and sampleVtcw>=727401613000000 and sampleVtcw<1035849616000000 order by sampleVtcw"
; query_database,q1,shutterPos,nrows,database="SORCE_L1S"
;
; q1="select TMD.Value 'dn', TMD.Value 'eu',  TMD.SCT_VTCW 'timetag' from TManalog TMD, SORCE_L1..TMnames TMN where "+$
;    "TMD.TMID = TMN.TMID and TMN.tlmId=101952 and TMD.SCT_VTCW >= 766540813000000 and TMD.SCT_VTCW < 766713613000000 "+$
;    "order by TMD.SCT_VTCW"
; query_database,q1,gain_mode_hc_a,nrows,database="SORCE_L1A"
;
; q1="select TMD.Value 'dn', TMD.Value 'eu',  TMD.SCT_VTCW 'timetag' from TManalog TMD, SORCE_L1..TMnames TMN where "+$
;    "TMD.TMID = TMN.TMID and TMN.tlmId=100580 and TMD.SCT_VTCW >= 766540813000000 and TMD.SCT_VTCW < 766713613000000 "+$
;    "order by TMD.SCT_VTCW"
; query_database,q1,half_cycle_a,nrows,database="SORCE_L1A"
;
; q1="select TMD.Value 'dn', TMD.Value 'eu',  TMD.SCT_VTCW 'timetag' from TManalog TMD, SORCE_L1..TMnames TMN where "+$
;    "TMD.TMID = TMN.TMID and TMN.tlmId=100038 and TMD.SCT_VTCW >= 766540813000000 and TMD.SCT_VTCW < 766713613000000 "+$
;    "order by TMD.SCT_VTCW"
; query_database,q1,hrt_out_a,nrows,database="SORCE_L1A"
;
; q1="select TMD.Value 'dn', TMD.Value 'eu',  TMD.SCT_VTCW 'timetag' from TManalog TMD, SORCE_L1..TMnames TMN where "+$
;    "TMD.TMID = TMN.TMID and TMN.tlmId=100331 and TMD.SCT_VTCW >= 766540813000000 and TMD.SCT_VTCW < 766713613000000 "+$
;    "order by TMD.SCT_VTCW"
; query_database,q1,hrt_in_a,nrows,database="SORCE_L1A"


        print,' getting the list of plans ...'
        plans0=[]
        plans1=[]
        if wrange[0] lt 1600d then $
            plans0=get_sorce_plan(t0, t1, sima=sima, simb=simb, /mission, activity='SolarQuickScan24')
        if wrange[1] gt 1600d then $
            plans1=get_sorce_plan(t0, t1, sima=sima, simb=simb, /mission, activity='SolarIRScan')

        plans=[plans0,plans1]
        s=sort(plans.(0))
        plans=plans[s]
        nplans=n_elements(plans)
        irradColumns=['Wavelength','SimCalibratedIrradiance']
        ;irradColumns=['Wavelength','SimCorrectedIrradiance']
        ;irradColumns=['Wavelength','SimUncorrectedIrradiance']

        query1 = 'SELECT c.microsecondsSinceGpsEpoch,w.wavelengthRef,c.irradiance '+$
            'FROM '+irradColumns[0]+' as w, '+irradColumns[1]+' as c '
        query1 += ' where w.version=c.version'
        query1 += ' and w.instrumentModeId=c.instrumentModeId'

        query_database,/reset
    endelse
    print,'    '+strtrim(string(nplans),2)+' plans to process'

    tsi=dblarr(nplans)
    timestamp=dblarr(nplans)-1
    if n_elements(refwaves) gt 0 then begin
       p1=where(refwaves ge wrange[0] and refwaves le wrange[1],count)
       if count gt 0 then my_refwaves=refwaves[p1] else my_refwaves=reffwaves
    endif

    for pl=0L, nplans-1L do begin
        print,pl+1,double(pl+1)/double(nplans+1)*100d,format='("   ",i0,"  (",i0," % )",$,%"\r")'
        waves=[]
        irrad=[]
        ttime=[]
        if size(alldata,/tname) eq 'STRUCT' then begin
            ; use the provided structure
            if valid_tag[0] ge 0 then begin
                if ptr_valid(alldata.(valid_tag)[pl]) then begin
                    p = where((*alldata.(valid_tag)[pl]).wavelength ge wrange[0] and (*alldata.(valid_tag)[pl]).wavelength le wrange[1], count)
                    if count ge mincount then begin
                        waves=[waves,(*alldata.(valid_tag)[pl]).wavelength[p]]
                        irrad=[irrad,(*alldata.(valid_tag)[pl]).irradiance[p]]
                        ttime=[ttime,mean((*alldata.(valid_tag)[pl]).timestamp[p])]
                    endif
                endif else begin
                    p = where(alldata.(valid_tag)[pl].wavelength ge wrange[0] and alldata.(valid_tag)[pl].wavelength le wrange[1], count)
                    if count ge mincount then begin
                        waves=[waves,alldata.(valid_tag)[pl].wavelength[p]]
                        irrad=[irrad,alldata.(valid_tag)[pl].irradiance[p]]
                        ttime=[ttime,mean(alldata.(valid_tag)[pl].timestamp[p])]
                    endif
                endelse
            endif else begin
                if ptr_valid(alldata[pl]) then begin
                    p = where((*alldata[pl]).wavelength ge wrange[0] and (*alldata[pl]).wavelength le wrange[1], count)
                    if count ge mincount then begin
                        waves=[waves,(*alldata[pl]).wavelength[p]]
                        irrad=[irrad,(*alldata[pl]).irradiance[p]]
                        ttime=[ttime,mean((*alldata[pl]).timestamp[p])]
                    endif
                endif else begin
                    p = where(alldata[pl].wavelength ge wrange[0] and alldata[pl].wavelength le wrange[1], count)
                    if count ge mincount then begin
                        waves=[waves,alldata[pl].wavelength[p]]
                        irrad=[irrad,alldata[pl].irradiance[p]]
                        ttime=[ttime,mean(alldata[pl].timestamp[p])]
                    endif
                endelse
            endelse

        endif else begin
            ; query the database for the data
            print,pl+1,double(pl+1)/double(nplans+1)*100d,format='("   ",i0,"  (",i0," % )",$,%"\r")'
            query3 = ' and c.microsecondsSinceGpsEpoch>='+strtrim(string(ulong64(sd2gps(plans[pl].starttime)*1d6)),2)
            query3 += ' and c.microsecondsSinceGpsEpoch<='+strtrim(string(ulong64(sd2gps(plans[pl].stoptime)*1d6)),2)
            query3 += ' and c.microsecondsSinceGpsEpoch=w.microsecondsSinceGpsEpoch order by c.microsecondsSinceGpsEpoch'
            if wrange[0] lt 309.d then begin
                if sima then $
                    query2 = ' and c.instrumentModeId=43' $
                else $
                    query2 = ' and c.instrumentModeId=47'
                query2 += ' and c.version='+strtrim(string(uv_version),2)
                ;spectuv=oJavaDbExchange->getAllValues(query1+query2+query3)
                query_database, query1+query2+query3, spectuv, nrows, dburl=dburl, user=user, password=password, dbdriver=dbdriver
                ;spectuv=get_science_product(irradColumns,plans[pl].starttime, plans[pl].stoptime, 43,/mission,$
                ;    version=version,user=user,password=password,dburl=dburl,dbdriver=dbdriver)
                if size(spectuv,/tname) ne 'STRUCT' then begin
                    print,'Error: no UV data found for those dates: ',$
                        strtrim(string([plans[pl].starttime, plans[pl].stoptime]),2)
                    continue
                endif
                w1=max([wrange[0],210d])
                w2=min([wrange[1],310d])
                p = where(spectuv.wavelengthref ge w1 and spectuv.wavelengthref le w2, count)
                if count gt 0 then begin
                    waves=[waves,spectuv[p].wavelengthref]
                    irrad=[irrad,spectuv[p].irradiance]
                    ttime=[ttime,mean(spectuv[p].(0))]
                endif
            endif

            if wrange[0] le 950.0  and wrange[1] ge 310.0 then begin
                if sima then $
                    query2 = ' and c.instrumentModeId=41' $
                else $
                    query2 = ' and c.instrumentModeId=45'
                query2 += ' and c.version='+strtrim(string(vis_version),2)
                ;spectvis=oJavaDbExchange->getAllValues(query1+query2+query3)
                query_database, query1+query2+query3, spectvis, nrows, dburl=dburl, user=user, password=password, dbdriver=dbdriver
                ;spectvis=get_science_product(irradColumns,plans[pl].starttime, plans[pl].stoptime, 41,/mission,$
                ;    version=version,user=user,password=password,dburl=dburl,dbdriver=dbdriver)
                if size(spectvis,/tname) ne 'STRUCT' then begin
                    print,'Error: no VIS data found for those dates: ',$
                        strtrim(string([plans[pl].starttime, plans[pl].stoptime]),2)
                    continue
                endif
                w1=max([wrange[0],210d])
                w2=min([wrange[1],950d])
                p = where(spectvis.wavelengthref ge w1 and spectvis.wavelengthref le w2, count)
                if count gt 0 then begin
                    waves=[waves,spectvis[p].wavelengthref]
                    irrad=[irrad,spectvis[p].irradiance]
                    ttime=[ttime,mean(spectvis[p].(0))]
                endif
            endif

            if wrange[1] gt 950.001 and wrange[0] le 1600.001 then begin
                if sima then $
                    query2 = ' and c.instrumentModeId=44' $
                else $
                    query2 = ' and c.instrumentModeId=48' 
                query2 += ' and c.version='+strtrim(string(ir_version),2)
                ;spectir=oJavaDbExchange->getAllValues(query1+query2+query3)
                query_database, query1+query2+query3, spectir, nrows, dburl=dburl, user=user, password=password, dbdriver=dbdriver
                ;spectir=get_science_product(irradColumns,plans[pl].starttime, plans[pl].stoptime, 44,/mission,$
                ;    version=version,user=user,password=password,dburl=dburl,dbdriver=dbdriver)
                if size(spectir,/tname) ne 'STRUCT' then begin
                    print,'Error: no IR data found for those dates: ',$
                        strtrim(string([plans[pl].starttime, plans[pl].stoptime]),2)
                    continue
                endif
                w1=max([wrange[0],950d])
                w2=min([wrange[1],1600d])
                p = where(spectir.wavelengthref ge w1 and spectir.wavelengthref le w2, count)
                if count gt 0 then begin
                    waves=[waves,spectir[p].wavelengthref]
                    irrad=[irrad,spectir[p].irradiance]
                    ttime=[ttime,mean(spectir[p].(0))]
                endif
            endif

            if wrange[0] gt 1600.0 then begin
                if sima then $
                    query2 = ' and c.instrumentModeId=31' $
                else $
                    query2 = ' and c.instrumentModeId=32' 
                query2 += ' and c.version='+strtrim(string(ir_version),2)
                ;spectesr=oJavaDbExchange->getAllValues(query1+query2+query3)
                query_database, query1+query2+query3, spectesr, nrows, dburl=dburl, user=user, password=password, dbdriver=dbdriver
                ;spectir=get_science_product(irradColumns,plans[pl].starttime, plans[pl].stoptime, 44,/mission,$
                ;    version=version,user=user,password=password,dburl=dburl,dbdriver=dbdriver)
                if size(spectesr,/tname) ne 'STRUCT' then begin
                    print,'Error: no ESR data found for those dates: ',$
                        strtrim(string([plans[pl].starttime, plans[pl].stoptime]),2)
                    continue
                endif
                w1=max([wrange[0],1600d])
                w2=min([wrange[1],2400d])
                p = where(spectesr.wavelengthref ge w1 and spectesr.wavelengthref le w2, count)
                if count gt 0 then begin
                    waves=[waves,spectesr[p].wavelengthref]
                    irrad=[irrad,spectesr[p].irradiance]
                    ttime=[ttime,mean(spectesr[p].(0))]
                endif
            endif

        endelse

        if n_elements(waves) gt 2 then begin
            timestamp[pl]=mean(ttime)
            s=sort(waves)
            q=uniq(waves[s])
            if n_elements(my_refwaves) eq 0 then begin
                tsi[pl]=int_tabulated(waves[s[q]], irrad[s[q]])
            endif else begin
                ; a fixed set of wavelengths was provided, interpolate to this grid before integrating
                val=interpol( irrad[s[q]], waves[s[q]], my_refwaves, /lsq)
                tsi[pl]=int_tabulated(my_refwaves, val)
            endelse
        endif

    endfor

    print,''
    p =where(timestamp gt 0d and tsi gt 0d,count)
    if count eq 0 then return,-1
    timestamp=gps2sd(timestamp[p]/1d6)
    tsi=tsi[p]
    if keyword_set(dailyavg) then begin
        ; handle the data from DOOP mode where we should combine the spectra of 3 consecutive orbits
        ; to make sure we have a complete spectra
        ;k=where(timestamp gt 4040d,count)
        ;if count gt 0 then begin
        ;    deltat=3d*90d*60d/86400d  ; (3-90 minutes orbites in days)
        ;    p=where((timestamp[k[1]:k[-1]] - timestamp[k[0]:k[-2]]) gt deltat,count)
        ;    stop
        ;endif




        ; average the tsi values within a single day
        result=histogram(timestamp,binsize=1.0,min=floor(min(timestamp)),max=ceil(max(timestamp)),location=xloc,revers=r)
        new_tsi=dblarr(n_elements(result))
        new_timestamp=dblarr(n_elements(result))
        for i=0L,n_elements(result)-1L do begin
            IF R[i] NE R[i+1] THEN BEGIN
                new_tsi[i] = MEAN(tsi[R[R[i] : R[i+1]-1]])
                new_timestamp[i] = MEAN(timestamp[R[R[i] : R[i+1]-1]])
            ENDIF
        endfor
        p=where(new_timestamp gt 0d)
        timestamp=new_timestamp[p]
        tsi=new_tsi[p]
    endif
    if keyword_set(align) then align_irrad,timestamp, tsi, /mission
    if NOT keyword_set(noplot) then $
        lineplot,timestamp, tsi, xtitle='Mission Day', ytitle='Integrated Irradiance', $
            ptitle='SimA Integrated Irradiance from '+strcompress(string(wrange[0],format='(F0.1)')+$
            ', '+string(wrange[1],format='(F0.1)')), charsize=1.4,psym=-4,title='TSI version='+strtrim(string(version),2)

    return,{timestamp:timestamp, tsi:tsi} 

end
