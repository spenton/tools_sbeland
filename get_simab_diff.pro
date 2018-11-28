;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Calculates the differences in SimCalibratedIrradiance between SimA and SimB
;   for the specified timerange and wavelength range. The SimB spectra is 
;   interpolated to the SimA wavelengths to make sure we're comparing the same wavelength.
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
;   instrumentModeId -
;      Specify the detector to process (only specify SimA and program will also select SimB)
;
; OPTIONAL INPUT PARAMETERS:
;   version -
;      Requests a specific version of the SimProfileIntegralCal, otherwise
;      pick the latest version valid for the specified time range.
;
;   wave0 -
;      Starting wavelength to compare
;
;   wave1 -
;      Ending wavelength to compare
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
;   Revision: $Id: get_simab_diff.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;------------------------------------------------------------------
function get_simab_diff, starttime, stoptime, instrumentModeId, wave0=wave0, wave1=wave1, $
    version=version, gps=gps, missionDays=missionDays, julianDays=julianDays, $
    dburl=dburl, dbdriver=dbdriver, user=user, password=password, dowave=dowave, wavefrac=wavefrac

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

    if n_elements(version) eq 0 then version=70
    if n_elements(user) eq 0 then     user='sbeland'
    if n_elements(password) eq 0 then password='sbeland1'
    if n_elements(dburl) eq 0 then    dbUrl='jdbc:sybase:Tds:lasp-db-dev:4100/SORCE_SIM_V20'
    if n_elements(dbdriver) eq 0 then dbDriver='com.sybase.jdbc3.jdbc.SybDriver'
    if n_elements(wave0) eq 0 then wave0=0d
    if n_elements(wave1) eq 0 then wave1=2600d

    jstmt = fjava_get_jdbc_statement(user=user, password=password, dbdriver=dbdriver, dburl=dburl) 
    oJavaDbExchange = OBJ_NEW("oJava_DbExchange", jstmt)

    if instrumentModeId lt 31 or instrumentModeId gt 48 then begin
        print,'program only works for photo-diode (modes 41 -> 48)'
        return,-1
    endif
    modea = instrumentModeId
    if instrumentModeId eq 31 then modeb=modea+1 else modeb=modea+4

    print,' getting the list of plans ...'
    ; only get the plans from SimB since SimA and SimB are taken at the same time
    if instrumentModeId eq 31 then activity='SolarIRScan' else activity='SolarQuickScan24'
    plans=get_sorce_plan(t0, t1, /simb, /mission, activity=activity)
    nplans=n_elements(plans)
    print,'    '+strtrim(string(nplans),2)+' plans to process'

    tsi=dblarr(nplans)
    timestamp=dblarr(nplans)-1

    irradColumns=['Wavelength','SimCalibratedIrradiance']
    ;irradColumns=['Wavelength','SimCorrectedIrradiance']

    query1 = 'SELECT c.microsecondsSinceGpsEpoch,w.wavelengthRef,c.irradiance '+$
        'FROM '+irradColumns[0]+' as w, '+irradColumns[1]+' as c '
    query1 += ' where w.version=c.version'
    query1 += ' and w.instrumentModeId=c.instrumentModeId'

    query_database,/reset
    outdata = replicate({starttime:0d, stoptime:0d, mean:0d, median:0d, mdev:0d, stdv:0d, min:0d, max:0d}, nplans)
    spectra=ptrarr(nplans)
    minwave=1d10
    maxwave=-1d10

    for pl=0L, nplans-1L do begin
        print,pl+1,double(pl+1)/double(nplans+1)*100d,format='("   ",i0,"  (",i0," % )",$,%"\r")'
        waves=[]
        irrad=[]
        ttime=[]
        query3 = ' and c.microsecondsSinceGpsEpoch>='+strtrim(string(ulong64(sd2gps(plans[pl].starttime)*1d6)),2)
        query3 += ' and c.microsecondsSinceGpsEpoch<='+strtrim(string(ulong64(sd2gps(plans[pl].stoptime)*1d6)),2)
        query3 += ' and c.microsecondsSinceGpsEpoch=w.microsecondsSinceGpsEpoch order by c.microsecondsSinceGpsEpoch'
        
        ; get the SimA data
        query2 = ' and c.instrumentModeId='+strtrim(string(modea),2)+' and c.version='+strtrim(string(version),2)
        spectA=oJavaDbExchange->getAllValues(query1+query2+query3)
        if size(spectA,/tname) ne 'STRUCT' then begin
            print,'Error: no SimA data found for those dates: ',$
                strtrim(string([plans[pl].starttime, plans[pl].stoptime]),2)
            continue
        endif

        ; get the SimB data
        query2 = ' and c.instrumentModeId='+strtrim(string(modeb),2)+' and c.version='+strtrim(string(version),2)
        spectB=oJavaDbExchange->getAllValues(query1+query2+query3)
        if size(spectB,/tname) ne 'STRUCT' then begin
            print,'Error: no SimB data found for those dates: ',$
                strtrim(string([plans[pl].starttime, plans[pl].stoptime]),2)
            continue
        endif

        if keyword_set(dowave) then begin
            if n_elements(refwave) eq 0 then begin
                refwave=spectA.wavelengthref
                p=where(refwave ge wave0 and refwave le wave1,count)
                refwave=refwave[p]
                minwave=refwave[0]
                maxwave=refwave[-1]
            endif
            irradA = interpol(spectA.irradiance, spectA.wavelengthref, refwave, /spline)
            irradB = interpol(spectB.irradiance, spectB.wavelengthref, refwave, /spline)
        endif else begin
            irradA = spectA.irradiance
            irradB = interpol(spectB.irradiance, spectB.wavelengthref, spectA.wavelengthref, /spline)
            p=where(spectA.wavelengthref ge wave0 and spectA.wavelengthref le wave1,count)
            minwave=min([spectA.wavelengthref,minwave])
            maxwave=max([spectA.wavelengthref,maxwave])
        endelse

        if count gt 0 then begin
            diff = (irradA[p] - irradB[p]) / irradA[p]
            if keyword_set(dowave) then spectra[pl] = ptr_new(diff)
            stats = moment(diff, mean=mean_diff, mdev=mdev_diff, sdev=sdev_diff)
            min_diff=min(diff,max=max_diff)
            outdata[pl].starttime = plans[pl].starttime
            outdata[pl].stoptime = plans[pl].stoptime
            outdata[pl].mean= mean(diff)
            outdata[pl].median= median(diff)
            outdata[pl].mdev= meanabsdev(diff)
            outdata[pl].stdv= stddev(diff)
            outdata[pl].min= min_diff
            outdata[pl].max= max_diff
        endif

    endfor

    print,''
    p =where(outdata.starttime gt 0d,count)
    if count eq 0 then return,-1
    outdata=outdata[p]
    plot_multi,outdata.starttime, outdata.stdv, xtitle='Mission Day', ytitle='STDEV of Difference in Irradiance', $
        title='SimA-SimB Irradiance for modes '+strtrim(string(modea),2)+', '+strtrim(string(modeb),2)+$
        ' from '+strcompress(string(minwave,format='(F0.1)')+', '+string(maxwave,format='(F0.1)'))+$
        ' version='+strtrim(string(version),2), charsize=1.4,psym=-4

    if keyword_set(dowave) then begin
        spectra=spectra[p]
        wdiff=dblarr(n_elements(refwave))
        value=dblarr(n_elements(spectra))
        for j=0,n_elements(refwave)-1 do begin
            value*=0d
            for i=0,n_elements(spectra)-1 do begin
                value[i]=(*spectra[i])[j]
            endfor
            wdiff[j]=stddev(value)
        endfor
        wavefrac={wavelength:refwave, fracdiff:wdiff}
        plot_multi,refwave,wdiff,xtitle='Wavelength (nm)',ytitle='STDDEV Fractional Difference in Irradiance', $
            title='SimA-SimB Irradiance for modes '+strtrim(string(modea),2)+', '+strtrim(string(modeb),2)+$
            ' from '+strtrim(string(t0,format='(F7.2)'),2)+' to '+strtrim(string(t1,format='(F7.2)'),2)+$
            ' version='+strtrim(string(version),2), charsize=1.4,psym=-4,/xst,/yst
    endif

    return,outdata

end
