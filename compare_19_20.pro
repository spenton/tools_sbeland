;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Compare the SimCalibratedIrradiance results from version 19 and 20
;   of the production code.
;
; CALLING SEQUENCE:
;   result = COMPARE_19_20(instrumentModeId, path20=path20)
;
; INPUT PARAMETERS:
;   instrumentModeId -
;      The instrument mode of interest:
;      31	SIM_A	ESR 
;      32	SIM_B	ESR
;      41	SIM_A	VIS1
;      42	SIM_A	VIS2
;      43	SIM_A	UV
;      44	SIM_A	IR
;      45	SIM_B	VIS1
;      46	SIM_B	VIS2
;      47	SIM_B	UV
;      48	SIM_B	IR
;
; OPTIONAL INPUT PARAMETERS:
;   PATH20 - 
;      If specified, will look in that directory for the V20 files.
;   gps - 
;      If set, indicates that input times are to be interpreted as
;      microseconds since Jan 6, 2080 midnight UT.  (default)
;   missionDays - 
;      If set, indicates that input times are to be interpreted as
;      days elapsed since launch, i.e. mission days.  If only the startTime
;      argument contains a non-zero value (stopTime is not supplied), then
;      all data for the specified mission day number will be retrieved.
;   julianDays - 
;      If set, indicates that the input times are to be interpreted as
;      julian day numbers.
;
; RETURNED PARAMETERS:
;   A structure with the spectrum and their differences.
;
; OPTIONAL OUTPUT PARAMETERS:
;      NONE
;
; EXAMPLE:  
;
; COMMON BLOCKS:
;      NONE
;
; REVISION HISTORY:
;   Revision: $Id: compare_19_20.pro,v 1.2 2018/11/26 16:16:02 spenton Exp $
;   Copied from compare_19_20.pro
;-
;
function compare_19_20, startTime, stopTime, instrumentModeId, v19=v19, v20=v20, $
         gps=gps, missionDays=missionDays, julianDays=julianDays, alldata=alldata, $
         plotwave=plotwave, skip19=skip19, irrad_table=irrad_table,irrad_column=irrad_column, $
         prismposition=prismposition, closest=closest, align=align, dettemp=dettemp, $
         user=user, password=password, dburl=dburl, dbdriver=dbdriver, degradation=degradation, $
         noplot=noplot, prismpos_v20=prismpos_v20, avg_doop=avg_doop, esrtable=esrtable, $
         psym=psym, irrad_version=irrad_version

    ; if alldata was passed, simply use it with plotwave and skip the database extraction
    if size(alldata,/tname) ne 'STRUCT' then  begin

        nrows=0
        ; validate the instrumentModeId
        if n_elements(instrumentModeId) eq 0 then begin
            doc_library,'get_sim_spectra'
            print,''
            print,'Missing instrumentModeId '
            return,-1
        endif else begin
            ; verify we got the right instrument
            valid_inst = [31,32,41,42,43,44,45,46,47,48]
            instrumentModeId = fix(instrumentModeId)
            if (where(valid_inst eq instrumentModeId))[0] lt 0 then  begin
                doc_library,'compare_19_20'
                print,''
                print,'Invalid instrumentModeId was provided'
                return,-1
            endif
        endelse

        ; define the telemetry items according to the specified mode
        case instrumentModeId of
            31: begin
                ; SIMA ESR
                instrument='SIMA'
            end
            41: begin
                ; SIMA VIS1
                instrument='SIMA'
            end
            42: begin
                ; SIMA VIS2
                instrument='SIMA'
            end
            43: begin
                ; SIMA UV
                instrument='SIMA'
            end
            44: begin
                ; SIMA IR
                instrument='SIMA'
            end
            32: begin
                ; SIMB ESR
                instrument='SIMB'
            end
            45: begin
                ; SIMB VIS1
                instrument='SIMB'
            end
            46: begin
                ; SIMB VIS2
                instrument='SIMB'
            end
            47: begin
                ; SIMB UV
                instrument='SIMB'
            end
            48: begin
                ; SIMB IR
                instrument='SIMB'
            end
        endcase

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

        if n_elements(user) eq 0 then user='sorce'
        if n_elements(password) eq 0 then password='sorcedb'
        if n_elements(dburl) eq 0 then dbUrl='jdbc:sybase:Tds:lasp-db-dev:4100/SORCE_SIM_DEV'
        if n_elements(dbdriver) eq 0 then dbDriver='com.sybase.jdbc3.jdbc.SybDriver'

        if n_elements(v19) eq 0 then v19='0'
        if n_elements(v20) eq 0 then begin
            query = 'SELECT MAX(version) v20 from SimCalibratedIrradiance where instrumentModeId='+strtrim(string(instrumentModeId),2)
            query = query + ' and microsecondsSinceGpsEpoch>='+strtrim(string(ulong64(t0)),2)
            query = query + ' and microsecondsSinceGpsEpoch<='+strtrim(string(ulong64(t1)),2)
            query_database,query, res, nrows, user=user, password=password, dburl=dburl, dbdriver=dbdriver
            if nrows eq 0 then begin
                print,'Error: no max(version) for SimCalibratedIrradiance found for the requested date/version'
                return,-1
            endif
            v20 = res[0].v20
        endif
        
        if n_elements(prismpos_v20) eq 0 then prismpos_v20=v20
        if n_elements(irrad_version) eq 0 then irrad_version=v20

        ; get the spectra for version 20
        print,'Extracting version 20 ...'
        if v20 eq 17 then begin
            ; if the version was specified as 17 -> read from the production database 
            ; with the right tables to extract from
            if n_elements(irrad_column) eq 0 then irrad_col='c.irradiance' else irrad_col='c.'+irrad_column
            if n_elements(irrad_table) eq 0 then irrad_table='SimCalibratedIrradiance'
            select_query = 'SELECT c.microsecondsSinceGpsEpoch,w.wavelengthRef,w.prismTransDegradation,'+irrad_col+',p.prismPosition'
            from_query = ' FROM '+irrad_table[0]+' as c, SimConvertedDataNumbers as p, SimProfileIntegral as w'
            where_query = ' where  w.instrumentModeId='+strtrim(string(instrumentModeId),2)
            where_query += ' and c.instrumentModeId=w.instrumentModeId'
            where_query += ' and p.instrumentModeId=w.instrumentModeId'
            where_query += ' and w.version='+strtrim(string(v20),2)
            where_query += ' and c.version='+strtrim(string(irrad_version),2)
            where_query += ' and p.version='+strtrim(string(prismpos_v20),2)
            where_query += ' and w.microsecondsSinceGpsEpoch>='+strtrim(string(ulong64(t0)),2)
            where_query += ' and w.microsecondsSinceGpsEpoch<='+strtrim(string(ulong64(t1)),2)
            where_query += ' and c.microsecondsSinceGpsEpoch=w.microsecondsSinceGpsEpoch'
            where_query += ' and p.microsecondsSinceGpsEpoch=w.microsecondsSinceGpsEpoch'
            order_query = ' order by c.microsecondsSinceGpsEpoch'
            query=select_query + from_query + where_query + order_query
            query_database, query, res, nrows, user='sorce', password='sorcedb', dbdriver=dbdriver, $
                dburl='jdbc:sybase:Tds:lasp-db1:4100/SORCE_SIM_V17'
        endif else begin
            if n_elements(irrad_column) eq 0 then irrad_col='c.irradiance' else irrad_col='c.'+irrad_column
            if n_elements(irrad_table) eq 0 then irrad_table='SimCalibratedIrradiance'
            select_query = 'SELECT c.microsecondsSinceGpsEpoch,w.wavelengthRef,w.correctedPrismPosition,'+irrad_col+',p.prismPosition'
            from_query = ' FROM '+irrad_table[0]+' as c, SimConvertedDataNumbers as p, Wavelength as w'
            where_query = ' where  w.instrumentModeId='+strtrim(string(instrumentModeId),2)
            where_query += ' and c.instrumentModeId=w.instrumentModeId'
            where_query += ' and p.instrumentModeId=w.instrumentModeId'
            where_query += ' and w.version='+strtrim(string(v20),2)
            where_query += ' and c.version='+strtrim(string(irrad_version),2)
            where_query += ' and p.version='+strtrim(string(prismpos_v20),2)
            where_query += ' and w.microsecondsSinceGpsEpoch>='+strtrim(string(ulong64(t0)),2)
            where_query += ' and w.microsecondsSinceGpsEpoch<='+strtrim(string(ulong64(t1)),2)
            where_query += ' and c.microsecondsSinceGpsEpoch=w.microsecondsSinceGpsEpoch'
            where_query += ' and p.microsecondsSinceGpsEpoch=w.microsecondsSinceGpsEpoch'
            order_query = ' order by c.microsecondsSinceGpsEpoch'
            ;print,query
            query=select_query + from_query + where_query + order_query
            if v20 eq 20 then begin
                ; if the version was specified as 20 -> read from the production database 
                ; where the officially released data resides
                dburl='jdbc:sybase:Tds:lasp-db1:4100/SORCE_SIM_V20'
            endif else if v20 eq 21 then begin
                ; if the version was specified as 21 -> read from the production database 
                ; where the officially released data resides
                dburl='jdbc:sybase:Tds:lasp-db1:4100/SORCE_SIM_V21'
            endif else if v20 eq 22 then begin
                ; if the version was specified as 22 -> read from the production database 
                ; where the officially released data resides
                dburl='jdbc:sybase:Tds:lasp-db1:4100/SORCE_SIM_V22'
            endif else if v20 eq 23 then begin
                ; if the version was specified as 22 -> read from the production database 
                ; where the officially released data resides
                dburl='jdbc:sybase:Tds:lasp-db1:4100/SORCE_SIM_V23'
            endif else if v20 eq 24 then begin
                ; if the version was specified as 22 -> read from the production database 
                ; where the officially released data resides
                dburl='jdbc:sybase:Tds:lasp-db1:4100/SORCE_SIM_V24'
            endif
            print,query
            query_database, query, res, nrows, user=user, password=password, dbdriver=dbdriver, dburl=dburl
        endelse

        if size(res,/tname) ne 'STRUCT' then begin
            print,'Error: no '+irrad_table+' found for the requested date/version'
            return,-1
        endif

        print,'   got '+strtrim(string(n_elements(res.(0))),2)+' entries '
        timerange20 = res.(0)
        wave20 = res.(1)
        irrad20 = res.(3)
        position20 = res.prismPosition
        if v20 ne 17 then correctedPosition20=res.correctedPrismPosition else correctedPosition20=dblarr(n_elements(wave20))

        if keyword_set(dettemp) then begin
            ; we need to interpolate the detector temperature
            query = 'SELECT microsecondsSinceGpsEpoch, temperature FROM DetectorTemperature '
            query += 'WHERE  instrumentModeId='+strtrim(string(instrumentModeId),2)
            query += ' and version='+strtrim(string(v20),2)
            query += ' and microsecondsSinceGpsEpoch>='+strtrim(string(ulong64(t0)),2)
            query += ' and microsecondsSinceGpsEpoch<='+strtrim(string(ulong64(t1)),2)
            query += ' order by microsecondsSinceGpsEpoch'
            query_database, query, dettemp20, nrows, user=user, password=password, dbdriver=dbdriver, dburl=dburl
        endif 

        if keyword_set(degradation) then begin
            ; remove the prism degradation from the irradiances
            query = 'SELECT * FROM SimPrismTransDegradation '
            query += 'where  instrumentModeId='+strtrim(string(instrumentModeId),2)
            query += ' and version='+strtrim(string(v20),2)
            query += ' and microsecondsSinceGpsEpoch>='+strtrim(string(ulong64(t0)),2)
            query += ' and microsecondsSinceGpsEpoch<='+strtrim(string(ulong64(t1)),2)
            query_database, query, pDeg, nrows, user=user, password=password, dbdriver=dbdriver, dburl=dburl
            if nrows eq 0 then begin
                print,'Error: no prismTransDegradation found for the requested date/version -> Returning data WITH degradation'
            endif else begin
                temp = interpol(pDeg.prismTransDegradation, pDeg.microsecondsSinceGpsEpoch, timerange20)
                irrad20 *= temp
            endelse
        endif 

        ; MODIFIED FOR DO-OP mode where we have partial scans over 2 to 3 consecutive orbits
        ; Process pre-DO-OP as before EXCEPT we will generate an array of pointers for each
        ; scan.
        ; SBeland 2016-11-27
        ;
        ck=where(timerange20 gt sd2gps(4040d)*1d6, count_doop,complement=compk)
        if compk[0] ne -1 then begin

            ; look for gaps of 200 seconds or more (48 hours for ESR to capture a complete FullScan 15 orbits)
            if instrumentModeId gt 32 or keyword_Set(esrtable) then  deltat=200d6 else deltat=48.0*3600d6
            p=where((timerange20[compk[1]:compk[-1]] - timerange20[compk[0]:compk[-2]]) gt deltat,count)
            spect20=ptrarr(count+1)
            p0=compk[0L]
            for j=0L,count do begin
                if j eq count then p1=compk[-1] else p1=compk[p[j]]
                if size(dettemp20,/tname) eq 'STRUCT' then begin
                    dettemp = interpol(dettemp20.temperature, dettemp20.microsecondsSinceGpsEpoch, timerange20[p0:p1])
                    spect20[j] = ptr_new($
                        {wavelength : wave20[p0:p1], $
                        irradiance : irrad20[p0:p1], $
                        timestamp : timerange20[p0:p1], $
                        prismPosition : position20[p0:p1], $
                        correctedPosition : correctedPosition20[p0:p1], $
                        dettemp : dettemp})
                endif else begin
                    spect20[j] = ptr_new($
                        {wavelength : wave20[p0:p1], $
                        irradiance : irrad20[p0:p1], $
                        timestamp : timerange20[p0:p1], $
                        prismPosition : position20[p0:p1], $
                        correctedPosition : correctedPosition20[p0:p1]})
                endelse
                p0=p1+1L
            endfor
        endif 

        if count_doop gt 0 then begin

            ; if p0 is not defined, no data before DOOP otherwise use the last value of p0
            if n_elements(p0) eq 0 then p0=0L

            ; look for gaps of 3 orbits and combine all the partial scans into one
            if instrumentModeId gt 32 then  deltat=3d*90d*60d6 else deltat=24.0*3600d6
            p=where((timerange20[ck[1]:ck[-1]] - timerange20[ck[0]:ck[-2]]) gt deltat,count)
            spect20_doop=ptrarr(count+1)
            if instrumentModeId eq 43 or instrumentModeId eq 47 then binsize=0.01d else binsize=0.1875d
            p0=ck[0L]
            for j=0L,count do begin
                if j eq count then p1=ck[-1] else p1=ck[p[j]]
                wavelength = wave20[p0:p1]
                irradiance = irrad20[p0:p1]
                timestamp = timerange20[p0:p1]
                prismPos= position20[p0:p1]
                correctedPosition = correctedPosition20[p0:p1]
                s=sort(wavelength)
                wavelength = wavelength[s]
                irradiance = irradiance[s]
                timestamp = timestamp[s]
                prismPos= prismPos[s]
                correctedPosition = correctedPosition[s]
                if keyword_set(avg_doop) then begin
                   ; average the data from the multiple scans within 3 orbits into one
                   yhis=histogram(wavelength,binsize=binsize,location=xhis,reverse=r)
                   ww=[] 
                   ird=[]
                   tt=[]
                   pp=[]
                   cp=[]
                   for k=0L,n_elements(xhis)-1L do begin 
                       if r[k] ne r[k+1] then begin 
                           ww = [ww,avg(wavelength[r[r[k]:r[k+1]-1]])] 
                           ird=[ird, avg(irradiance[r[r[k]:r[k+1]-1]])] 
                           tt=[tt, avg(timestamp[r[r[k]:r[k+1]-1]])] 
                           pp=[pp, avg(prismPos[r[r[k]:r[k+1]-1]])] 
                           cp=[cp, avg(correctedPosition[r[r[k]:r[k+1]-1]])] 
                       endif 
                   endfor
                   wavelength=ww
                   irradiance=ird
                   timestamp=tt
                   prismPos=pp
                   correctedPosition=cp
                endif

                if size(dettemp20,/tname) eq 'STRUCT' then begin
                    dettemp = interpol(dettemp20.temperature, dettemp20.microsecondsSinceGpsEpoch, timestamp)
                    spect20_doop[j] = ptr_new($
                        {wavelength : wavelength, $
                        irradiance : irradiance, $
                        timestamp : timestamp, $
                        prismPosition : prismPos, $
                        correctedPosition : correctedPosition, $
                        dettemp : dettemp})
                endif else begin
                    spect20_doop[j] = ptr_new($
                        {wavelength : wavelength, $
                        irradiance : irradiance, $
                        timestamp : timestamp, $
                        prismPosition : prismPos, $
                        correctedPosition : correctedPosition})
                endelse
                p0=p1+1L
            endfor
            if n_elements(spect20) gt 0 then spect20=[spect20,spect20_doop] else spect20=spect20_doop
        endif 


        ; get the spectra for version 19
        if NOT keyword_set(skip19) then begin
            print,'Extracting version 19 ...'
            if n_elements(irrad_column) eq 0 then irrad_col='c.irradiance' else irrad_col='c.'+irrad_column
            if n_elements(irrad_table) eq 0 then irrad_table='SimCalibratedIrradiance'
            select_query = 'SELECT c.microsecondsSinceGpsEpoch,w.wavelengthRef,w.correctedPrismPosition,'+irrad_col+',p.prismPosition'
            from_query = ' FROM '+irrad_table[0]+' as c, SimConvertedDataNumbers as p, Wavelength as w'
            where_query = ' where  w.instrumentModeId='+strtrim(string(instrumentModeId),2)
            where_query += ' and c.instrumentModeId=w.instrumentModeId'
            where_query += ' and p.instrumentModeId=w.instrumentModeId'
            where_query += ' and w.version='+strtrim(string(v19),2)
            where_query += ' and c.version=w.version'
            where_query += ' and p.version=w.version'
            where_query += ' and w.microsecondsSinceGpsEpoch>='+strtrim(string(ulong64(t0)),2)
            where_query += ' and w.microsecondsSinceGpsEpoch<='+strtrim(string(ulong64(t1)),2)
            where_query += ' and c.microsecondsSinceGpsEpoch=w.microsecondsSinceGpsEpoch'
            where_query += ' and p.microsecondsSinceGpsEpoch=w.microsecondsSinceGpsEpoch'
            order_query = ' order by c.microsecondsSinceGpsEpoch'
            query=select_query + from_query + where_query + order_query
            if v19 eq 19 then begin
                ; if the version was specified as 19 -> read from the production database 
                ; where the officially released data resides
                dburl='jdbc:sybase:Tds:lasp-db1:4100/SORCE_SIM_V19'
            endif else if v19 eq 20 then begin
                ; if the version was specified as 19 -> read from the production database 
                ; where the officially released data resides
                dburl='jdbc:sybase:Tds:lasp-db1:4100/SORCE_SIM_V20'
            endif
            query_database, query, res, nrows, user=user, password=password, dbdriver=dbdriver, dburl=dburl
            if size(res,/tname) ne 'STRUCT' then begin
                print,'Error: no SimCalibratedIrradiance found for the requested date/version'
                return,-1
            endif
            timerange19 = res.(0)
            wave19 = res.(1)
            irrad19 = res.(3)
            position19 = res.prismPosition
            if v19 ne 17 then correctedPosition19=res.correctedPrismPosition else correctedPosition19=dblarr(n_elements(wave19))

            if keyword_set(dettemp) then begin
                ; we need to interpolate the detector temperature
                query = 'SELECT microsecondsSinceGpsEpoch, temperature FROM DetectorTemperature '
                query += 'WHERE  instrumentModeId='+strtrim(string(instrumentModeId),2)
                query += ' and version='+strtrim(string(v19),2)
                query += ' and microsecondsSinceGpsEpoch>='+strtrim(string(ulong64(t0)),2)
                query += ' and microsecondsSinceGpsEpoch<='+strtrim(string(ulong64(t1)),2)
                query += ' order by microsecondsSinceGpsEpoch'
                query_database, query, dettemp19, nrows, user=user, password=password, dbdriver=dbdriver, dburl=dburl
            endif 

            if keyword_set(degradation) then begin
                ; remove the prism degradation from the irradiances
                query = 'SELECT * FROM SimPrismTransDegradation '
                query += 'where  instrumentModeId='+strtrim(string(instrumentModeId),2)
                query += ' and version='+strtrim(string(v19),2)
                query += ' and microsecondsSinceGpsEpoch>='+strtrim(string(ulong64(t0)),2)
                query += ' and microsecondsSinceGpsEpoch<='+strtrim(string(ulong64(t1)),2)
                query_database, query, pDeg, nrows, user=user, password=password, dbdriver=dbdriver, dburl=dburl
                if nrows eq 0 then begin
                    print,'Error: no prismTransDegradation found for the requested date/version -> Returning data WITH degradation'
                endif else begin
                    temp = interpol(pDeg.prismTransDegradation, pDeg.microsecondsSinceGpsEpoch, timerange19)
                    irrad19 *= temp
                endelse
            endif 

            ; organize the data as an array of structure organize per scan
            ; assuming each spectra has same number of elements
            ;p=where(wave19[0:-2]-wave19[1:-1] gt 0.0d,count)
            if instrumentModeId gt 32 then deltat=240d6 else deltat=24.0*3600d6
            p = where((timerange19[1:-1] - timerange19[0:-2]) gt deltat,count)
            ; count the number of points per spectra
            if count eq 0 then nwave19 = n_elements(wave19) else begin
                p=[0,p,n_elements(wave19)-1]
                nwave19=max(p[1:-1]-p[0:-2])
            endelse
            count_scans = count + 1L
            spect19 = replicate({timestamp:dblarr(nwave19),wavelength:dblarr(nwave19), dettemp:dblarr(nwave19), $
                irradiance:dblarr(nwave19), prismPosition:dblarr(nwave19), correctedPosition:dblarr(nwave19)}, count_scans)
            p0=0L
            for j=0,count_scans-1L do begin
                ; for the last scan point to end of array
                if j eq count then p1=n_elements(wave19)-1L else p1=p[j+1]
                ; check if we have more points then the size of the array
                if (p1-p0) ge n_elements(spect19[j].wavelength) then p2=p1-1 else p2=p1
                spect19[j].wavelength = wave19[p0:p2]
                spect19[j].irradiance = irrad19[p0:p2]
                spect19[j].timestamp = timerange19[p0:p2]
                spect19[j].prismPosition = position19[p0:p2]
                spect19[j].correctedPosition = correctedPosition19[p0:p2]
                if size(dettemp19,/tname) eq 'STRUCT' then $
                    spect19[j].dettemp = interpol(dettemp19.temperature, dettemp19.microsecondsSinceGpsEpoch, spect19[j].timestamp)
                p0=p1+1L
            endfor

            alldata = {spect19:temporary(spect19), spect20:temporary(spect20)}
        endif else alldata={spect20:temporary(spect20)}

    endif   ; alldata structure was passed as an input


    nscans20= n_elements(alldata.spect20)
    if max(strpos(tag_names(alldata),'SPECT19')) lt 0 then begin
        skip19=1
        nscans19=0
    endif else begin
        if not keyword_set(skip19) then begin
            skip19=0
            nscans19= n_elements(alldata.spect19)
        endif else begin
            skip19=1
            nscans19=0
        endelse
    endelse

    print,'   got '+strtrim(string(nscans20),2)+' scans '
    if n_elements(psym) eq 0 then psym=-4

    if n_elements(prismposition) gt 0 then begin
        ; plot the data from a specfic prismposition instead of interpolating 
        if skip19 eq 0 then begin
            ts_19_time=dblarr(nscans19)
            ts_19_irrad=dblarr(nscans19)
            ts_19_wave=dblarr(nscans19)
            tnames = tag_names(alldata.spect19)
            if max(strpos(tag_names(alldata.spect19),'PRISMPOSITION')) ge 0 then begin
                for i=0L,nscans19-1L do begin
                    p19 = where(alldata.spect19[i].prismposition eq prismposition[0],count)
                    if count gt 0 then begin
                        ts_19_irrad[i] = mean(alldata.spect19[i].irradiance[p19])
                        ts_19_time[i] = mean(alldata.spect19[i].timestamp[p19])
                        ts_19_wave[i] = mean(alldata.spect19[i].wavelength[p19])
                    endif 
                endfor

                p=where(finite(ts_19_irrad) eq 1 and ts_19_irrad gt 0d)
                ts_19_irrad=ts_19_irrad[p]
                ts_19_time=gps2sd(ts_19_time[p]/1d6)
                ts_19_wave=ts_19_wave[p]

                title='SIM_V19 version='+strtrim(string(v19),2)
                title+=' pos='+strtrim(string(min(prismposition[0]),format='(f0.1)'),2)
                if keyword_set(align) then align_irrad,ts_19_time,ts_19_irrad,/mission
                if NOT keyword_set(noplot) then $
                    lineplot,ts_19_time,ts_19_irrad,psym=psym,xtitle='Mission Day',ytitle='Irradiance',title=title,charsize=1.5
            endif
        endif

        ts_20_time=dblarr(nscans20)
        ts_20_irrad=dblarr(nscans20)
        ts_20_wave=dblarr(nscans20)
        tnames = tag_names(alldata.spect20)
        if max(strpos(tag_names(alldata.spect20),'PRISMPOSITION')) ge 0 then begin
			ww=[]
            for i=0L,nscans20-1L do begin
                p20 = where(alldata.spect20[i].prismposition eq prismposition[-1],count)
                if count gt 0 then begin
                    ts_20_irrad[i] = mean(alldata.spect20[i].irradiance[p20])
                    ts_20_time[i] = mean(alldata.spect20[i].timestamp[p20])
                    ts_20_wave[i] = mean(alldata.spect20[i].wavelength[p20])
					ww=[ww,alldata.spect20[i].wavelength[p20]]
                endif 
            endfor

            p=where(finite(ts_20_irrad) eq 1 and ts_20_irrad gt 0d)
            ts_20_irrad=ts_20_irrad[p]
            ts_20_time=gps2sd(ts_20_time[p]/1d6)
            ts_20_wave=ts_20_wave[p]

            title='SIM_V20 version='+strtrim(string(v20),2)
            title+=' pos='+strtrim(string(min(prismposition[-1]),format='(f0.1)'),2)
            if keyword_set(align) then align_irrad,ts_20_time,ts_20_irrad,/mission
            if NOT keyword_set(noplot) then $
                lineplot,ts_20_time,ts_20_irrad,psym=psym,xtitle='Mission Day',ytitle='Irradiance',title=title,charsize=1.5
        endif
        if n_elements(ww) eq 0 then return,-1
        if skip19 eq 0 then begin
            return, {plotwave:median(ww), v20_wave:ts_20_wave, v19_irrad:ts_19_irrad, v20_irrad:ts_20_irrad, $
                v19_wave:ts_19_wave, timestamp19:ts_19_time, timestamp20:ts_20_time }
        endif else begin
            return, {plotwave:median(ww), v20_wave:ts_20_wave, v20_irrad:ts_20_irrad, timestamp20:ts_20_time }
        endelse

    endif else if n_elements(plotwave) gt 0 then begin

        ; plot the irradiance for the specified wavelength for all spectra
        ; interpolating to match the requested wavelength
        if nscans19 gt 0 then ts_19_time=dblarr(nscans19)
        if nscans19 gt 0 then ts_19_irrad=dblarr(nscans19)
        if nscans19 gt 0 then ts_19_dtemp=dblarr(nscans19)
        ts_20_time=dblarr(nscans20)
        ts_20_irrad=dblarr(nscans20)
        ts_20_dtemp=dblarr(nscans20)
        ; if plotwave is an array of size two we assume we want to integrate
        ; the irradiance within that range 
        ; in this case interpolate at every nm (and divide by number of nm)
        if n_elements(plotwave) gt 1 then plotwave = dindgen(ceil(max(plotwave,min=mn)-mn)+1)+mn
        mnw=min(plotwave,max=mxw)
        nwave = double(n_elements(plotwave))
        for i=0L,nscans19-1L do begin
            print,i+1,double(i+1)/double(nscans19+1)*100d,format='("   ",i0,"  (",i0," % )",$,%"\r")'
            if skip19 eq 0 then begin
                temp=alldata.spect19[i].wavelength
                p=where(temp gt 200.0,count)
                mnsw=min(temp[p],max=mxsw)
                if mnsw le mnw and mxsw ge mxw then begin
                    diff = min(abs(temp[p] - plotwave),pos)
                    ts_19_time[i] = mean(alldata.spect19[i].timestamp[p[pos[0]]])
                    if keyword_set(closest) then begin
                        ts_19_irrad[i] = mean(alldata.spect19[i].irradiance[p[pos[0]]])
                    endif else begin
                        if count ge 3 then begin
                            s=sort(temp[p])
                            if instrumentModeId gt 32 then begin
                                y2=spl_init(temp[p[s]], alldata.spect19[i].irradiance[p[s]])
                                irrad = spl_interp(temp[p[s]], alldata.spect19[i].irradiance[p[s]], y2, plotwave)
                                ;irrad = interpol(alldata.spect19[i].irradiance[p[s]],temp[p[s]],plotwave,/lsq)
                            endif else begin 
                                irrad = interpol(alldata.spect19[i].irradiance[p[s]],temp[p[s]],plotwave,/spline)
                            endelse
                            ts_19_irrad[i] = total(irrad)
                        endif
                    endelse
                endif
            endif
        endfor

        for i=0L,nscans20-1L do begin
            print,i+1,double(i+1)/double(nscans20+1)*100d,format='("   ",i0,"  (",i0," % )",$,%"\r")'
            if ptr_valid(alldata.spect20[i]) then begin
                tempw = (*alldata.spect20[i]).wavelength
                tempt = (*alldata.spect20[i]).timestamp
                tempi = (*alldata.spect20[i]).irradiance
                if keyword_set(dettemp) then tempdt=(*alldata.spect20[i]).dettemp
            endif else begin
                tempw = alldata.spect20[i].wavelength
                tempt = alldata.spect20[i].timestamp
                tempi = alldata.spect20[i].irradiance
                if keyword_set(dettemp) then tempdt=alldata.spect20[i].dettemp
            endelse
            p=where(tempw gt 200.0, count)
            mnsw=min(tempw[p],max=mxsw)
            if mnsw le mnw and mxsw ge mxw then begin
                diff = min(abs(tempw[p] - plotwave),pos)
                ts_20_time[i] = tempt[p[pos[0]]]
                if keyword_set(dettemp) then ts_20_dtemp[i] = tempdt[p[pos[0]]]
                if keyword_set(closest) then begin
                    ts_20_irrad[i] = tempi[p[pos[0]]]
                endif else begin
                    if count lt 3 then continue
                    s=sort(tempw[p])
                    ; make sure the spectra fully covers the wavelength range requested
                    ; if in DO-OP mode, we have multiple scans within one group which
                    ; messes up the interpolation
                    ; when the end of timerange is in DO-OP mode, average the data from multiple scans
                    if max(tempt[p]) gt sd2gps(4040d)*1d6  and ~keyword_set(avg_doop) then begin
                        yhis=histogram(tempw[p[s]],binsize=0.01d,location=xhis,reverse=r)
                        ww=[] 
                        ird=[]
                        for k=0L,n_elements(xhis)-1L do begin 
                            if r[k] ne r[k+1] then begin 
                                ww = [ww,avg(tempw[p[s[r[r[k]:r[k+1]-1]]]])] 
                                ird=[ird, avg(tempi[p[s[r[r[k]:r[k+1]-1]]]])] 
                            endif 
                        endfor
                        irrad = interpol(ird,ww,plotwave,/lsq)
                    endif else begin
                        irrad = interpol(tempi[p[s]],tempw[p[s]],plotwave,/lsq)
                    endelse
                    ts_20_irrad[i] = total(irrad)
                endelse
            endif
        endfor

        p=where(finite(ts_20_irrad) eq 1 and ts_20_irrad gt 0d)
        ts_20_irrad=ts_20_irrad[p]
        ts_20_time=ts_20_time[p]
        if keyword_set(dettemp) then ts_20_dtemp = ts_20_dtemp[p]
        print,''

        if skip19 eq 0 then begin
            p=where(finite(ts_19_irrad) eq 1 and ts_19_irrad gt 0d)
            ts_19_irrad=ts_19_irrad[p]
            ts_19_time=ts_19_time[p]
            p=where(ts_19_irrad gt 0.0d and ts_19_irrad lt 10.0d*nwave,count)
            ts_19_irrad=ts_19_irrad[p]
            ts_19_time=gps2sd(ts_19_time[p]/1d6)
            if keyword_set(align) then align_irrad,ts_19_time,ts_19_irrad,/mission
            title='SIM_V19 version='+strtrim(string(v19),2)
            if nwave eq 1  then $
                title+=' wl=['+strtrim(string(min(plotwave),format='(f0.1)'),2)+'-'+strtrim(string(max(plotwave),format='(f0.1)'),2)+']'
            if NOT keyword_set(noplot) then $
                lineplot,ts_19_time,ts_19_irrad,psym=psym,xtitle='Mission Day',ytitle='Irradiance',title=title,charsize=1.5
        endif

        if n_elements(ts_20_irrad) gt 0 then begin
            ;p=where(ts_20_irrad gt 0.0d and ts_20_irrad lt 10.0d*nwave,count)
            p=where(finite(ts_20_irrad) eq 1 and ts_20_irrad gt 0d)
            ts_20_irrad=ts_20_irrad[p]
            ts_20_time=gps2sd(ts_20_time[p]/1d6)
            if keyword_set(align) then align_irrad,ts_20_time,ts_20_irrad,/mission
            title='SIM_V20 version='+strtrim(string(v20),2)
            if nwave eq 1  then $
                title+=' wl='+strtrim(string(plotwave,format='(f0.1)'),2) $
            else $
                title+=' wl=['+strtrim(string(min(plotwave),format='(f0.1)'),2)+'-'+strtrim(string(max(plotwave),format='(f0.1)'),2)+']'
            if NOT keyword_set(noplot) then $
                lineplot,ts_20_time,ts_20_irrad,psym=psym,xtitle='Mission Day',ytitle='Irradiance',title=title,charsize=1.5

            if skip19 eq 0 then begin
                if keyword_set(dettemp) then begin
                    dt=reform(alldata.spect19.dettemp,n_elements(alldata.spect19.dettemp))
                    ts=reform(alldata.spect19.timestamp,n_elements(alldata.spect19.timestamp))
                    s=sort(ts)
                    dettemp19=interpol(dt[s], ts[s], sd2gps(ts_19_time)*1d6)
                    dt=reform(alldata.spect20.dettemp,n_elements(alldata.spect20.dettemp))
                    ts=reform(alldata.spect20.timestamp,n_elements(alldata.spect20.timestamp))
                    s=sort(ts)
                    dettemp20=interpol(dt[s], ts[s], sd2gps(ts_20_time)*1d6)
                    return, {plotwave:plotwave, v20_wave:replicate(plotwave,n_elements(ts_20_irrad)), $
                        v19_irrad:ts_19_irrad, v20_irrad:ts_20_irrad, $
                        v19_wave:replicate(plotwave,n_elements(ts_19_irrad)), $
                        timestamp19:ts_19_time, timestamp20:ts_20_time, dettemp19:dettemp19, dettemp20:dettemp20 }
                endif else begin
                    return, {plotwave:plotwave, v20_wave:replicate(plotwave,n_elements(ts_20_irrad)), $
                        v19_irrad:ts_19_irrad, v20_irrad:ts_20_irrad, $
                        v19_wave:replicate(plotwave,n_elements(ts_19_irrad)), $
                        timestamp19:ts_19_time, timestamp20:ts_20_time}
                endelse
            endif else begin
                if keyword_set(dettemp) then begin
                    return, {plotwave:plotwave, v20_wave:replicate(plotwave,n_elements(ts_20_irrad)), $
                    v20_irrad:ts_20_irrad, timestamp20:ts_20_time, dettemp20:ts_20_dtemp }
                endif else begin
                    return, {plotwave:plotwave, v20_wave:replicate(plotwave,n_elements(ts_20_irrad)), $
                        v20_irrad:ts_20_irrad, timestamp20:ts_20_time }
                endelse
            endelse
        endif else begin
            print,'No Data to plot'
        endelse
    endif

    return, nscans20
 
end
