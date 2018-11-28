;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Compare the SimCalibratedIrradiance results from version 17 and 19
;   of the production code.
;
; CALLING SEQUENCE:
;   result = COMPARE_17_19(instrumentModeId, path19=path19)
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
;   PATH19 - 
;      If specified, will look in that directory for the V19 files.
;   gps - 
;      If set, indicates that input times are to be interpreted as
;      microseconds since Jan 6, 1980 midnight UT.  (default)
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
;   Revision: $Id: compare_17_19.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;   Changed to access the data from the database instead of ASCII files
;-
;
function compare_17_19, startTime, stopTime, instrumentModeId, v17=v17, v19=v19, $
         gps=gps, missionDays=missionDays, julianDays=julianDays, alldata=alldata, $
         plotwave=plotwave, skip17=skip17, irrad_table=irrad_table,irrad_column=irrad_column, $
         prismposition=prismposition, closest=closest, align=align, $
         user=user, password=password, dburl=dburl, dbdriver=dbdriver, degradation=degradation


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
                doc_library,'compare_17_19'
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

        if n_elements(user) eq 0 then user='sbeland'
        if n_elements(password) eq 0 then password='sbeland1'
        if n_elements(dburl) eq 0 then dbUrl='jdbc:sybase:Tds:lasp-db-dev:4100/SORCE_SIM_V19'
        if n_elements(dbdriver) eq 0 then dbDriver='com.sybase.jdbc3.jdbc.SybDriver'

        jstmt = fjava_get_jdbc_statement(user=user, password=password, dbdriver=dbdriver, dburl=dburl) 
        oJavaDbExchange = OBJ_NEW("oJava_DbExchange", jstmt)

        if n_elements(v17) eq 0 then v17='0'
        if n_elements(v19) eq 0 then begin
            query = 'SELECT MAX(version) v19 from SimCalibratedIrradiance where instrumentModeId='+strtrim(string(instrumentModeId),2)
            query = query + ' and microsecondsSinceGpsEpoch>='+strtrim(string(ulong64(t0)),2)
            query = query + ' and microsecondsSinceGpsEpoch<='+strtrim(string(ulong64(t1)),2)
            res=oJavaDbExchange->getAllValues(query)
            if size(res,/tname) ne 'STRUCT' then begin
                print,'Error: no max(version) for SimCalibratedIrradiance found for the requested date/version'
                return,-1
            endif
            v19 = res[0].v19
        endif

        ; get the spectra for version 19
        print,'Extracting version 19 ...'
        if n_elements(irrad_column) eq 0 then irrad_column='c.irradiance' else irrad_column='c.'+irrad_column
        if n_elements(irrad_table) eq 0 then begin
            query = 'SELECT c.microsecondsSinceGpsEpoch,w.wavelengthRef,'+irrad_column+',p.prismPosition '+$
                    'FROM SimCalibratedIrradiance as c, SimConvertedDataNumbers as p, Wavelength as w '
        endif else begin
            query = 'SELECT c.microsecondsSinceGpsEpoch,w.wavelengthRef,'+irrad_column+',p.prismPosition '+$
                    'FROM '+irrad_table[0]+' as c, SimConvertedDataNumbers as p, Wavelength as w '
        endelse
        query += 'where  w.instrumentModeId='+strtrim(string(instrumentModeId),2)
        query += ' and c.instrumentModeId=w.instrumentModeId'
        query += ' and p.instrumentModeId=w.instrumentModeId'
        query += ' and w.version='+strtrim(string(v19),2)
        query += ' and c.version=w.version'
        query += ' and p.version=w.version'
        query += ' and w.microsecondsSinceGpsEpoch>='+strtrim(string(ulong64(t0)),2)
        query += ' and w.microsecondsSinceGpsEpoch<='+strtrim(string(ulong64(t1)),2)
        query += ' and c.microsecondsSinceGpsEpoch=w.microsecondsSinceGpsEpoch'
        query += ' and p.microsecondsSinceGpsEpoch=w.microsecondsSinceGpsEpoch order by c.microsecondsSinceGpsEpoch'

        ; print,query
        res=oJavaDbExchange->getAllValues(query)
        if size(res,/tname) ne 'STRUCT' then begin
            print,'Error: no SimCalibratedIrradiance found for the requested date/version'
            return,-1
        endif

        print,'   got '+strtrim(string(n_elements(res.(0))),2)+' entries '

        timerange19 = res.(0)
        wave19 = res.(1)
        irrad19 = res.(2)
        position19 = res.prismPosition

        if keyword_set(degradation) then begin
            query = 'SELECT * FROM SimPrismTransDegradation '
            query += 'where  instrumentModeId='+strtrim(string(instrumentModeId),2)
            query += ' and version='+strtrim(string(v19),2)
            query += ' and microsecondsSinceGpsEpoch>='+strtrim(string(ulong64(t0)),2)
            query += ' and microsecondsSinceGpsEpoch<='+strtrim(string(ulong64(t1)),2)
            pDeg=oJavaDbExchange->getAllValues(query)
            if size(pDeg,/tname) ne 'STRUCT' then begin
                print,'Error: no prismTransDegradation found for the requested date/version -> Returning data WITH degradation'
            endif else begin
                temp = interpol(pDeg.prismTransDegradation, pDeg.microsecondsSinceGpsEpoch, timerange19)
                irrad19 *= temp
            endelse
        endif 

        ; organize the data as an array of structure organize per scan
        ; assuming each spectra has same number of elements
        ; p=where(wave19[0:-2]-wave19[1:-1] gt 0.0d,count)
        ; look for gaps of 240 seconds or more 
        p = where((timerange19[1:-1] - timerange19[0:-2]) gt 240.0d6,count)
        ; count the number of points per spectra
        if count eq 0 then nwave19 = n_elements(wave19) else begin
            p=[0,p,n_elements(wave19)-1]
            nwave19=max(p[1:-1]-p[0:-2])+1
        endelse
        count_scans = count + 1L
        spect19 = replicate({timestamp:dblarr(nwave19),wavelength:dblarr(nwave19), $
            irradiance:dblarr(nwave19), prismPosition:dblarr(nwave19)}, count_scans)
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
            p0=p1+1L
        endfor


        ; get the spectra for version 17
        if NOT keyword_set(skip17) then begin
            print,'Extracting version 17 ...'
            if n_elements(irrad_table) eq 0 then begin
                query = 'SELECT c.microsecondsSinceGpsEpoch,w.wavelengthRef,'+irrad_column+',p.prismPosition '+$
                        ' FROM SimCalibratedIrradiance as c, SimConvertedDataNumbers as p, '
            endif else begin
                query = 'SELECT c.microsecondsSinceGpsEpoch,w.wavelengthRef,'+irrad_column+',p.prismPosition '+$
                        ' FROM '+irrad_table[0]+' as c, SimConvertedDataNumbers as p, '
            endelse
            query += 'Wavelength as w where  w.instrumentModeId='+strtrim(string(instrumentModeId),2)
            query += ' and c.instrumentModeId=w.instrumentModeId'
            query += ' and p.instrumentModeId=w.instrumentModeId'
            query += ' and w.version='+strtrim(string(v17),2)
            query += ' and c.version=w.version'
            query += ' and p.version=w.version'
            query += ' and w.microsecondsSinceGpsEpoch>='+strtrim(string(ulong64(t0)),2)
            query += ' and w.microsecondsSinceGpsEpoch<='+strtrim(string(ulong64(t1)),2)
            query += ' and c.microsecondsSinceGpsEpoch=w.microsecondsSinceGpsEpoch'
            query += ' and p.microsecondsSinceGpsEpoch=w.microsecondsSinceGpsEpoch order by c.microsecondsSinceGpsEpoch'
            res=oJavaDbExchange->getAllValues(query)
            if size(res,/tname) ne 'STRUCT' then begin
                print,'Error: no SimCalibratedIrradiance found for the requested date/version'
                return,-1
            endif
            timerange17 = res.(0)
            wave17 = res.(1)
            irrad17 = res.(2)
            position17 = res.prismPosition

            if keyword_set(degradation) then begin
                query = 'SELECT * FROM SimPrismTransDegradation '
                query += 'where  instrumentModeId='+strtrim(string(instrumentModeId),2)
                query += ' and version='+strtrim(string(v17),2)
                query += ' and microsecondsSinceGpsEpoch>='+strtrim(string(ulong64(t0)),2)
                query += ' and microsecondsSinceGpsEpoch<='+strtrim(string(ulong64(t1)),2)
                pDeg=oJavaDbExchange->getAllValues(query)
                if size(pDeg,/tname) ne 'STRUCT' then begin
                    print,'Error: no prismTransDegradation found for the requested date/version -> Returning data WITH degradation'
                endif else begin
                    temp = interpol(pDeg.prismTransDegradation, pDeg.microsecondsSinceGpsEpoch, timerange17)
                    irrad17 *= pDeg
                endelse
            endif 

            ; organize the data as an array of structure organize per scan
            ; assuming each spectra has same number of elements
            ;p=where(wave17[0:-2]-wave17[1:-1] gt 0.0d,count)
            p = where((timerange17[1:-1] - timerange17[0:-2]) gt 240.0d6,count)
            ; count the number of points per spectra
            if count eq 0 then nwave17 = n_elements(wave17) else begin
                p=[0,p,n_elements(wave17)-1]
                nwave17=max(p[1:-1]-p[0:-2])
            endelse
            count_scans = count + 1L
            spect17 = replicate({timestamp:dblarr(nwave17),wavelength:dblarr(nwave17), $
            irradiance:dblarr(nwave17), prismPosition:dblarr(nwave17)}, count_scans)
            p0=0L
            for j=0,count_scans-1L do begin
                ; for the last scan point to end of array
                if j eq count then p1=n_elements(wave17)-1L else p1=p[j+1]
                ; check if we have more points then the size of the array
                if (p1-p0) ge n_elements(spect17[j].wavelength) then p2=p1-1 else p2=p1
                spect17[j].wavelength = wave17[p0:p2]
                spect17[j].irradiance = irrad17[p0:p2]
                spect17[j].timestamp = timerange17[p0:p2]
                spect17[j].prismPosition = position17[p0:p2]
                p0=p1+1L
            endfor

            alldata = {spect17:temporary(spect17), spect19:temporary(spect19)}
        endif else alldata={spect19:temporary(spect19)}

    endif   ; alldata structure was passed as an input


    nscans19= n_elements(alldata.spect19)
    if max(strpos(tag_names(alldata),'SPECT17')) lt 0 then begin
        skip17=1
        nscans17=0
    endif else begin
        if not keyword_set(skip17) then begin
            skip17=0
            nscans17= n_elements(alldata.spect17)
        endif else begin
            skip17=1
            nscans17=0
        endelse
    endelse

    print,'   got '+strtrim(string(nscans19),2)+' scans '

    if n_elements(prismposition) gt 0 then begin
        ; plot the data from a specfic prismposition instead of interpolating 
        ts_19_time=dblarr(nscans19)
        ts_19_irrad=dblarr(nscans19)
        if skip17 eq 0 then begin
            ts_17_time=dblarr(nscans17)
            ts_17_irrad=dblarr(nscans17)
            tnames = tag_names(alldata.spect17)
            if max(strpos(tag_names(alldata.spect17),'PRISMPOSITION')) ge 0 then begin
                for i=0L,nscans17-1L do begin
                    p17 = where(alldata.spect17[i].prismposition eq prismposition[0],count)
                    if count gt 0 then begin
                        ts_17_irrad[i] = mean(alldata.spect17[i].irradiance[p17])
                        ts_17_time[i] = mean(alldata.spect17[i].timestamp[p17])
                    endif 
                endfor

                p=where(abs(ts_17_irrad) gt 0d)
                ts_17_irrad=ts_17_irrad[p]
                ts_17_time=gps2sd(ts_17_time[p]/1d6)

                title='SIM_V17 version='+strtrim(string(v17),2)
                title+=' pos='+strtrim(string(min(prismposition[0]),format='(f0.1)'),2)
                if keyword_set(align) then align_irrad,ts_17_time,ts_17_irrad,/mission
                lineplot,ts_17_time,ts_17_irrad,psym=-4,xtitle='Mission Day',ytitle='Irradiance',title=title,charsize=1.5
            endif
        endif

        tnames = tag_names(alldata.spect19)
        if max(strpos(tag_names(alldata.spect19),'PRISMPOSITION')) ge 0 then begin
			ww=[]
            for i=0L,nscans19-1L do begin
                p19 = where(alldata.spect19[i].prismposition eq prismposition[-1],count)
                if count gt 0 then begin
                    ts_19_irrad[i] = mean(alldata.spect19[i].irradiance[p19])
                    ts_19_time[i] = mean(alldata.spect19[i].timestamp[p19])
					ww=[ww,alldata.spect19[i].wavelength[p19]]
                endif 
            endfor

            p=where(abs(ts_19_irrad) gt 0d)
            ts_19_irrad=ts_19_irrad[p]
            ts_19_time=gps2sd(ts_19_time[p]/1d6)

            title='SIM_V19 version='+strtrim(string(v19),2)
            title+=' pos='+strtrim(string(min(prismposition[-1]),format='(f0.1)'),2)
            if keyword_set(align) then align_irrad,ts_19_time,ts_19_irrad,/mission
            lineplot,ts_19_time,ts_19_irrad,psym=-4,xtitle='Mission Day',ytitle='Irradiance',title=title,charsize=1.5
        endif
        if n_elements(ww) eq 0 then return,-1
        if skip17 eq 0 then begin
            return, {plotwave:median(ww), v17_irrad:ts_17_irrad, v19_irrad:ts_19_irrad, timestamp17:ts_17_time, timestamp19:ts_19_time }
        endif else begin
            return, {plotwave:median(ww), v19_irrad:ts_19_irrad, timestamp19:ts_19_time }
        endelse

    endif else if n_elements(plotwave) gt 0 then begin

        ; plot the irradiance for the specified wavelength for all spectra
        ; interpolating to match the requested wavelength
        if nscans17 gt 0 then ts_17_time=dblarr(nscans17)
        if nscans17 gt 0 then ts_17_irrad=dblarr(nscans17)
        ts_19_time=dblarr(nscans19)
        ts_19_irrad=dblarr(nscans19)
        ; if plotwave is an array of size two we assume we want to integrate
        ; the irradiance within that range 
        ; in this case interpolate at every nm (and divide by number of nm)
        if n_elements(plotwave) gt 1 then plotwave = dindgen(ceil(max(plotwave,min=mn)-mn)+1)+mn
        mnw=min(plotwave,max=mxw)
        nwave = double(n_elements(plotwave))
        for i=0L,nscans17-1L do begin
            print,i+1,double(i+1)/double(nscans17+1)*100d,format='("   ",i0,"  (",i0," % )",$,%"\r")'
            if skip17 eq 0 then begin
                temp=alldata.spect17[i].wavelength
                p=where(temp gt 0.0,count)
                mnsw=min(temp[p],max=mxsw)
                if mnsw le mnw and mxsw ge mxw then begin
                    diff = min(abs(temp[p] - plotwave),pos)
                    ts_17_time[i] = mean(alldata.spect17[i].timestamp[p[pos[0]]])
                    if keyword_set(closest) then begin
                        ts_17_irrad[i] = mean(alldata.spect17[i].irradiance[p[pos[0]]])
                    endif else begin
                        if count ge 3 then begin
                            if instrumentModeId gt 32 then $
                                irrad = interpol(alldata.spect17[i].irradiance[p],temp[p],plotwave,/lsq) $
                            else $
                                irrad = interpol(alldata.spect17[i].irradiance[p],temp[p],plotwave)
                            ts_17_irrad[i] = total(irrad)/nwave
                        endif
                    endelse
                endif
            endif
        endfor

        for i=0L,nscans19-1L do begin
            print,i+1,double(i+1)/double(nscans19+1)*100d,format='("   ",i0,"  (",i0," % )",$,%"\r")'
            tempw = alldata.spect19[i].wavelength
            tempt = alldata.spect19[i].timestamp
            tempi = alldata.spect19[i].irradiance
            p=where(tempw gt 0.0, count)
            mnsw=min(tempw[p],max=mxsw)
            if mnsw le mnw and mxsw ge mxw then begin
                diff = min(abs(tempw[p] - plotwave),pos)
                ts_19_time[i] = tempt[p[pos[0]]]
                if keyword_set(closest) then begin
                    ts_19_irrad[i] = tempi[p[pos[0]]]
                endif else begin
                    if count lt 3 then continue
                    k0=0
                    k1=n_elements(p)-1L
                    if instrumentModeId gt 32 then $
                        irrad = interpol(tempi[p[(pos[0]-4)>k0:(pos[0]+4)<k1]],tempw[p[(pos[0]-4)>k0:(pos[0]+4)<k1]],plotwave,/lsq) $
                    else $
                        irrad = interpol(tempi[p[(pos[0]-4)>k0:(pos[0]+4)<k1]],tempw[p[(pos[0]-4)>k0:(pos[0]+4)<k1]],plotwave)
                    ts_19_irrad[i] = total(irrad)/nwave
                endelse
            endif
        endfor

        p=where(abs(ts_19_irrad) gt 0d)
        ts_19_irrad=ts_19_irrad[p]
        ts_19_time=ts_19_time[p]
        print,''

        if skip17 eq 0 then begin
            p=where(abs(ts_17_irrad) gt 0d)
            ts_17_irrad=ts_17_irrad[p]
            ts_17_time=ts_17_time[p]
            p=where(ts_17_irrad gt 0.0d and ts_17_irrad lt 10.0d*nwave,count)
            ts_17_irrad=ts_17_irrad[p]
            ts_17_time=gps2sd(ts_17_time[p]/1d6)
            if keyword_set(align) then align_irrad,ts_17_time,ts_17_irrad,/mission
            title='SIM_V17 version='+strtrim(string(v17),2)
            if nwave eq 1  then $
                title+=' wl=['+strtrim(string(min(plotwave),format='(f0.1)'),2)+'-'+strtrim(string(max(plotwave),format='(f0.1)'),2)+']'
            lineplot,ts_17_time,ts_17_irrad,psym=-4,xtitle='Mission Day',ytitle='Irradiance',title=title,charsize=1.5
        endif

        if n_elements(ts_19_irrad) gt 0 then begin
            ;p=where(ts_19_irrad gt 0.0d and ts_19_irrad lt 10.0d*nwave,count)
            p=where(abs(ts_19_irrad) gt 0d)
            ts_19_irrad=ts_19_irrad[p]
            ts_19_time=gps2sd(ts_19_time[p]/1d6)
            if keyword_set(align) then align_irrad,ts_19_time,ts_19_irrad,/mission
            title='SIM_V19 version='+strtrim(string(v19),2)
            if nwave eq 1  then $
                title+=' wl='+strtrim(string(plotwave,format='(f0.1)'),2) $
            else $
                title+=' wl=['+strtrim(string(min(plotwave),format='(f0.1)'),2)+'-'+strtrim(string(max(plotwave),format='(f0.1)'),2)+']'
            lineplot,ts_19_time,ts_19_irrad,psym=-4,xtitle='Mission Day',ytitle='Irradiance',title=title,charsize=1.5
            if skip17 eq 0 then begin
                return, {plotwave:plotwave, v17_irrad:ts_17_irrad, v19_irrad:ts_19_irrad, timestamp17:ts_17_time, timestamp19:ts_19_time }
            endif else begin
                return, {plotwave:plotwave, v19_irrad:ts_19_irrad, timestamp19:ts_19_time }
            endelse
        endif else begin
            print,'No Data to plot'
        endelse
    endif

    return, nscans19
 
end
