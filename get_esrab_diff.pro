;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Calculates the differences in SimCalibratedIrradiance between SimA and SimB
;   for the specified timerange and wavelength range. The program extracts the
;   data from the ESR table scans.
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
;   Revision: $Id: get_esrab_diff.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;------------------------------------------------------------------
function get_esrab_diff, starttime, stoptime, modea, modeb, $
    version=version, gps=gps, missionDays=missionDays, julianDays=julianDays, $
    dburl=dburl, dbdriver=dbdriver, user=user, password=password

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

    if n_elements(version) eq 0 then version=95
    if n_elements(user) eq 0 then     user='sbeland'
    if n_elements(password) eq 0 then password='sbeland1'
    if n_elements(dburl) eq 0 then    dbUrl='jdbc:sybase:Tds:lasp-db-dev:4100/SORCE_SIM_V20'
    if n_elements(dbdriver) eq 0 then dbDriver='com.sybase.jdbc3.jdbc.SybDriver'

    jstmt = fjava_get_jdbc_statement(user=user, password=password, dbdriver=dbdriver, dburl=dburl) 
    oJavaDbExchange = OBJ_NEW("oJava_DbExchange", jstmt)

    print,' getting the list of plans ...'
    ; only get the plans from SimB since SimA and SimB are taken at the same time

    case modea of
        31: begin
            activityA='SolarESRMode' 
            sima=1
            simb=0
            end
        32: begin
            activityA='SolarESRMode' 
            sima=0
            simb=1
            end
        44: begin
            activityA='SolarQuickScan24' 
            sima=1
            simb=0
            end
    endcase

    case modeb of
        31: begin
            activityB='SolarESRMode' 
            sima=1
            simb=0
            end
        32: begin
            activityB='SolarESRMode' 
            sima=0
            simb=1
            end
        44: begin
            activityB='SolarQuickScan24' 
            sima=1
            simb=0
            end
    endcase

    plansA=get_sorce_plan(t0, t1, sima=sima, simb=simb, /mission, activity=activityA)
    nplansA=n_elements(plansA)
    plansB=get_sorce_plan(t0, t1, sima=sima, simb=simb, /mission, activity=activityB)
    nplansB=n_elements(plansB)

    print,'    '+strtrim(string(nplansA),2)+' plans to process'

    irradColumns=['Wavelength','SimCalibratedIrradiance','SimConvertedDataNumbers']

    query1 = 'SELECT c.microsecondsSinceGpsEpoch,w.wavelengthRef,c.irradiance,p.prismPosition '+$
        'FROM '+irradColumns[0]+' as w, '+irradColumns[1]+' as c, '+irradColumns[2]+' as p '
    query1 += ' where w.version=c.version and p.version=c.version '
    query1 += ' and w.instrumentModeId=c.instrumentModeId and p.instrumentModeId=c.instrumentModeId'

    query_database,/reset
    spectA=[]
    spectB=[]

    ; get SimA data
    for pl=0L, nplansA-1L do begin
        print,pl+1,double(pl+1)/double(nplansA+1)*100d,format='("   ",i0,"  (",i0," % )",$,%"\r")'
        waves=[]
        irrad=[]
        ttime=[]
        query2 = ' and c.instrumentModeId='+strtrim(string(modea),2)+' and c.version='+strtrim(string(version),2)

        query3 = ' and c.microsecondsSinceGpsEpoch>='+strtrim(string(ulong64(sd2gps(plansA[pl].starttime)*1d6)),2)
        query3 += ' and c.microsecondsSinceGpsEpoch<='+strtrim(string(ulong64(sd2gps(plansA[pl].stoptime)*1d6)),2)
        query3 += ' and w.microsecondsSinceGpsEpoch=c.microsecondsSinceGpsEpoch '
        query3 += ' and p.microsecondsSinceGpsEpoch=c.microsecondsSinceGpsEpoch order by c.microsecondsSinceGpsEpoch'
        
        tmp=oJavaDbExchange->getAllValues(query1+query2+query3)
        if size(tmp,/tname) ne 'STRUCT' then begin
            print,'Error: no SimA data found for those dates: ',$
                strtrim(string([plansA[pl].starttime, plansA[pl].stoptime]),2)
            continue
        endif
        spectA=[spectA, tmp]
    endfor

    ; get SimB data
    for pl=0L, nplansB-1L do begin
        print,pl+1,double(pl+1)/double(nplansB+1)*100d,format='("   ",i0,"  (",i0," % )",$,%"\r")'
        waves=[]
        irrad=[]
        ttime=[]
        query2 = ' and c.instrumentModeId='+strtrim(string(modeB),2)+' and c.version='+strtrim(string(version),2)

        query3 = ' and c.microsecondsSinceGpsEpoch>='+strtrim(string(ulong64(sd2gps(plansB[pl].starttime)*1d6)),2)
        query3 += ' and c.microsecondsSinceGpsEpoch<='+strtrim(string(ulong64(sd2gps(plansB[pl].stoptime)*1d6)),2)
        query3 += ' and w.microsecondsSinceGpsEpoch=c.microsecondsSinceGpsEpoch '
        query3 += ' and p.microsecondsSinceGpsEpoch=c.microsecondsSinceGpsEpoch order by c.microsecondsSinceGpsEpoch'
        
        tmp=oJavaDbExchange->getAllValues(query1+query2+query3)
        if size(tmp,/tname) ne 'STRUCT' then begin
            print,'Error: no SimB data found for those dates: ',$
                strtrim(string([plansB[pl].starttime, plansB[pl].stoptime]),2)
            continue
        endif
        spectB=[spectB, tmp]
    endfor

    ; find the unique prism positions for SimB and match them with SimA
    minwa=min(spectA.wavelengthref, max=maxwa)
    minwb=min(spectB.wavelengthref, max=maxwb)
    minwave=max([minwa,minwb])-1d
    maxwave=min([maxwa,maxwb])+1d
    p=where(spectA.wavelengthref ge minwave and spectA.wavelengthref le maxwave,count)
    if count eq 0 then begin
        print,'Error:  NO overlap in wavelengths'
        return,-1
    endif
    spectA=spectA[p]

    p=where(spectB.wavelengthref ge minwave and spectB.wavelengthref le maxwave,count)
    if count eq 0 then begin
        print,'Error:  NO overlap in wavelengths'
        return,-1
    endif
    spectB=spectB[p]

    sb=sort(spectB.prismposition)
    qb=uniq(spectB[sb].prismposition)
    prismPosB = spectB[sb[qb]].prismPosition
    prismPosA = PrismPosB * 0d
    for i=0,n_elements(prismPosB)-1L do begin
        mn=min(abs(spectB[sb[qb[i]]].wavelengthref - spectA.wavelengthref),pos)
        prismPosA[i] = spectA[pos[0]].prismPosition
    endfor

    ; go through the list of SimB wavelengths and match SimA closest in time
    refwave=[]
    wdiff=[]
    nptsa=n_elements(spectA)
    for i=0,n_elements(prismPosB)-1L do begin
        pa=where(spectA.prismPosition eq prismPosA[i], counta)
        pb=where(spectB.prismPosition eq prismPosB[i], countb)
        diff=[]
        for j=0,countb-1 do begin
            ; find the closest spectA in time at this position
            mn=min(abs(spectB[pb[j]].MICROSECONDSSINCEGPSEPOCH - spectA[pa].MICROSECONDSSINCEGPSEPOCH),pos)
            if mn gt 86400d6 then continue
            ; interpolate if modeb is IR
            if modea eq 44 then begin
                k=where(spectA.MICROSECONDSSINCEGPSEPOCH eq spectA[pa[pos]].MICROSECONDSSINCEGPSEPOCH)
                irradA = interpol(spectA[(k-3)>0:(k+3)<nptsa-1].irradiance, spectA[(k-3)>0:(k+3)<nptsa-1].wavelengthref, spectB[pb[j]].wavelengthref)
            endif else begin
                irradA=spectA[pa[pos]].irradiance
            endelse
            diff = [diff, (spectB[pb[j]].irradiance - irradA) / irradA]
        endfor

        refWave=[refWave, mean(spectB[pb].wavelengthref)]
        wdiff = [wdiff, stddev(diff)]
    endfor

    plot_multi,refwave,wdiff,xtitle='Wavelength (nm)',ytitle='STDDEV Fractional Difference in Irradiance', $
        title='SimA-SimB Irradiance for modes '+strtrim(string(modea),2)+', '+strtrim(string(modeb),2)+$
        ' from '+strtrim(string(t0,format='(F7.2)'),2)+' to '+strtrim(string(t1,format='(F7.2)'),2)+$
        ' version='+strtrim(string(version),2), charsize=1.4,psym=-4,/xst,/yst

    return, {refwave:refwave, wdiff:wdiff}

end
