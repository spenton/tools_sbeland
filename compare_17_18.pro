;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Compare the SimCalibratedIrradiance results from version 17 and 18
;   of the production code.
;
; CALLING SEQUENCE:
;   result = COMPARE_17_18(t0, t1, instrumentModeId, /hrtout, /missionDays)
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
;   HRTOUT - 
;      If specified, will only return result for when the HRT was OUT.
;      By default, the program ignores the state of the HRT.
;   HRTIN - 
;      If specified, will only return result for when the HRT was IN.
;      By default, the program ignores the state of the HRT.
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
;   Revision: $Id: compare_17_18.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;-
;
function compare_17_18, startTime, stopTime, instrumentModeId, $
         gps=gps, missionDays=missionDays, julianDays=julianDays, $
         hrtout=hrtout, hrtin=hrtin, first=first, plotwave=plotwave, $
         skip17=skip17, alldata=alldata

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
                doc_library,'compare_17_18'
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
        if instrument eq 'SIMA' then begin
            rad_trap_in  = 20245L
            rad_trap_out = 20035L
        endif else begin
            rad_trap_in  = 20374L
            rad_trap_out = 20815L
        endelse

        ; convert the start and stop times to microsecondsSinceGpsEpoch
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

        ; prepare the database connection 
        stmt17 = fjava_get_jdbc_statement(user='sorce', password='sorcedb', dburl=dbURL, $
            dbdriver=dbdriver, server='sorce-db', database="SORCE")
        stmt18 = fjava_get_jdbc_statement(user='sorce', password='sorcedb', dburl=dbURL, $
            dbdriver=dbdriver, server='sorce-db', database="SORCE_SIM_V18")
        conn17 = OBJ_NEW("oJava_DbExchange", stmt17)
        conn18 = OBJ_NEW("oJava_DbExchange", stmt18)

        ; get the irradiance for version 17
        query1 = "SELECT microsecondsSinceGpsEpoch, irradiance FROM SimCalibratedIrradiance WHERE "
        query1 = query1+"instrumentModeId="+strtrim(string(instrumentModeId),2)
        query1 = query1+" AND microsecondsSinceGpsEpoch >= "+strtrim(string(ulong64(t0)),2)
        query1 = query1+" AND microsecondsSinceGpsEpoch <= "+strtrim(string(ulong64(t1)),2)
        query1 = query1+" AND version=17"
        query1 = query1+" ORDER BY microsecondsSinceGpsEpoch"
        res=conn17->getAllValues(query1)
        if size(res,/tname) ne 'STRUCT' then begin
            print,'Error: no SimCalibratedIrradiance found for the requested date/version'
            return,-1
        endif
        irrad17 = res.irradiance
        timerange1 = res.microsecondsSinceGpsEpoch
     
        ; get the wavelength for version 17
        query1 = "SELECT microsecondsSinceGpsEpoch, wavelengthRef FROM SimProfileIntegral WHERE "
        query1 = query1+"instrumentModeId="+strtrim(string(instrumentModeId),2)
        query1 = query1+" AND microsecondsSinceGpsEpoch >= "+strtrim(string(ulong64(t0)),2)
        query1 = query1+" AND microsecondsSinceGpsEpoch <= "+strtrim(string(ulong64(t1)),2)
        query1 = query1+" AND version=17"
        query1 = query1+" ORDER BY microsecondsSinceGpsEpoch"
        res=conn17->getAllValues(query1)
        if size(res,/tname) ne 'STRUCT' then begin
            print,'Error: no SimProfileIntegral found for the requested date/version'
            return,-1
        endif
        timerange2 = res.microsecondsSinceGpsEpoch
        wave17 = res.wavelengthref
        if n_elements(wave17) ne n_elements(irrad17) then begin
            ; find the matching times 
            match,timerange1,timerange2,sub1,sub2
            irrad17=irrad17[sub1]
            wave17=wave17[sub2]
            timerange17=timerange2[sub2]
        endif else timerange17=timerange2

        if keyword_set(first) then begin
            ; requested to only keep the first spectra found within the time span
            delta = abs(wave17[1]-wave17[0]) / 10.0
            p = where(abs(wave17[1:-1]-wave17[0]) lt delta,count)
            if count gt 0 then begin
                pos=lindgen(p[0])
                wave17=wave17[pos]
                irrad17=irrad17[pos]
                timerange17=timerange17[pos]
                t0=timerange[0]
                t1=timerange[-1]
            endif
        endif

        ; get the irradiance for version 18
        query1 = "SELECT microsecondsSinceGpsEpoch, irradiance FROM SimCalibratedIrradiance WHERE "
        query1 = query1+"instrumentModeId="+strtrim(string(instrumentModeId),2)
        query1 = query1+" AND microsecondsSinceGpsEpoch >= "+strtrim(string(ulong64(t0)),2)
        query1 = query1+" AND microsecondsSinceGpsEpoch <= "+strtrim(string(ulong64(t1)),2)
        query1 = query1+" AND version=18"
        query1 = query1+" ORDER BY microsecondsSinceGpsEpoch"
        res=conn18->getAllValues(query1)
        if size(res,/tname) ne 'STRUCT' then begin
            print,'Error: no SimCalibratedIrradiance found for the requested date/version'
            return,-1
        endif
        irrad18 = res.irradiance
        timerange1=res.microsecondsSinceGpsEpoch
     
        ; get the wavelength for version 18
        query1 = "SELECT microsecondsSinceGpsEpoch, wavelengthRef FROM SimProfileIntegral WHERE "
        query1 = query1+"instrumentModeId="+strtrim(string(instrumentModeId),2)
        query1 = query1+" AND microsecondsSinceGpsEpoch >= "+strtrim(string(ulong64(t0)),2)
        query1 = query1+" AND microsecondsSinceGpsEpoch <= "+strtrim(string(ulong64(t1)),2)
        query1 = query1+" AND version=18"
        query1 = query1+" ORDER BY microsecondsSinceGpsEpoch"
        res=conn18->getAllValues(query1)
        if size(res,/tname) ne 'STRUCT' then begin
            print,'Error: no SimProfileIntegral found for the requested date/version'
            return,-1
        endif
        wave18 = res.wavelengthref
        timerange2=res.microsecondsSinceGpsEpoch
        if n_elements(wave18) ne n_elements(irrad18) then begin
            ; find the matching times 
            match,timerange1,timerange2,sub1,sub2
            irrad18=irrad18[sub1]
            wave18=wave18[sub2]
            timerange18=timerange2[sub2]
        endif else timerange18=timerange2

     

        ; if the state of the hrt was specified then trim the data for these only
        if keyword_set(hrtin) or keyword_set(hrtout) then begin
           q1 = "select SCT_VTCW 'timetag', Value 'dn' from SORCE_L1.dbo.TMdiscrete where TMID="
           q2 = " AND SCT_VTCW BETWEEN "+strtrim(ulong64(t0),2)+" AND " + $
                strtrim(ulong64(t1),2)+" ORDER by SCT_VTCW"
           if keyword_set(hrtout) then begin
               q_str = q1+strtrim(rad_trap_out,2)+q2
               rad_trap_out=conn17->getAllValues(q_str)
               hrtout_dn  = dethin_data(rad_trap_out.dn, rad_trap_out.timetag, timerange17)
               keep = where(hrtout_dn gt 0,count)
               if count eq 0 then begin
                   print,'Error: no data found with HRTOUT during specified timerange'
                   return,-1
               endif
               wave17=wave17[keep]
               irrad17=irrad17[keep]
               hrtout_dn  = dethin_data(rad_trap_out.dn, rad_trap_out.timetag, timerange18)
               keep = where(hrtout_dn gt 0,count)
               if count eq 0 then begin
                   print,'Error: no data found with HRTOUT during specified timerange'
                   return,-1
               endif
               wave18=wave18[keep]
               irrad18=irrad18[keep]
           endif else begin
               q_str = q1+strtrim(rad_trap_in,2)+q2
               rad_trap_in=conn17->getAllValues(q_str)
               hrtin_dn  = dethin_data(rad_trap_in.dn, rad_trap_in.timetag, timerange17)
               keep = where(hrtout_dn gt 0,count)
               if count eq 0 then begin
                   print,'Error: no data found with HRTOUT during specified timerange'
                   return,-1
               endif
               wave17=wave17[keep]
               irrad17=irrad17[keep]
               hrtin_dn  = dethin_data(rad_trap_in.dn, rad_trap_in.timetag, timerange18)
               keep = where(hrtout_dn gt 0,count)
               if count eq 0 then begin
                   print,'Error: no data found with HRTOUT during specified timerange'
                   return,-1
               endif
               wave18=wave18[keep]
               irrad18=irrad18[keep]
           endelse
        endif

        ; organize the data as an array of structure organize per scan
        ; assuming each spectra has same number of elements
        p=where(wave17[0:-2]-wave17[1:-1] gt 0.0d,count)
        ; count the number of points per spectra
        if count eq 0 then nwave17 = n_elements(wave17) else nwave17=p[0]+1L
        count_scans = count + 1L
        spect17 = replicate({timestamp:dblarr(nwave17),wavelength:dblarr(nwave17), irradiance:dblarr(nwave17)}, count_scans)
        p0=0L
        for j=0,count_scans-1L do begin
            ; for the last scan point to end of array
            if j eq count then p1=n_elements(wave17)-1L else p1=p[j]
            ; check if we have more points then the size of the array
            if (p1-p0) ge n_elements(spect17[j].wavelength) then $
                p2=p0+n_elements(spect17[j].wavelength)-1 else p2=p1
            spect17[j].wavelength = wave17[p0:p2]
            spect17[j].irradiance = irrad17[p0:p2]
            spect17[j].timestamp = timerange17[p0:p2]
            p0=p1+1L
        endfor

        ; organize the data as an array of structure organize per scan
        ; assuming each spectra has same number of elements
        p=where(wave18[0:-2]-wave18[1:-1] gt 0.0d,count)
        ; count the number of points per spectra
        if count eq 0 then nwave18 = n_elements(wave18) else nwave18=p[0]+1L
        count_scans = count + 1L
        spect18 = replicate({timestamp:dblarr(nwave18),wavelength:dblarr(nwave18), irradiance:dblarr(nwave18)}, count_scans)
        p0=0L
        for j=0,count_scans-1L do begin
            ; for the last scan point to end of array
            if j eq count then p1=n_elements(wave18)-1L else p1=p[j]
            ; check if we have more points then the size of the array
            if (p1-p0) ge n_elements(spect18[j].wavelength) then $
                p2=p0+n_elements(spect18[j].wavelength)-1 else p2=p1
            spect18[j].wavelength = wave18[p0:p2]
            spect18[j].irradiance = irrad18[p0:p2]
            spect18[j].timestamp = timerange18[p0:p2]
            p0=p1+1L
        endfor
        alldata = {spect17:temporary(spect17), spect18:temporary(spect18)}
        nscans= min([n_elements(alldata.spect17),n_elements(alldata.spect18)])

    endif   ; alldata structure was passed as an input

    if max(strpos(tag_names(alldata),'SPECT17')) lt 0 then begin
        skip17=1
        nscans=n_elements(alldata.spect18) 
    endif else begin
        if not keyword_set(skip17) then begin
            skip17=0
            nscans= min([n_elements(alldata.spect17),n_elements(alldata.spect18)])
        endif else begin
            skip17=1
            nscans=n_elements(alldata.spect18) 
        endelse
    endelse

    if n_elements(plotwave) gt 0 then begin
        ; plot the irradiance for the specified wavelength for all spectra
        ; interpolating to match the requested wavelength
        ts_17_time=[]
        ts_17_irrad=[]
        ts_18_time=[]
        ts_18_irrad=[]
        ; if plotwave is an array of size two we assume we want to integrate
        ; the irradiance within that range 
        ; in this case interpolate at every nm (and divide by number of nm)
        if n_elements(plotwave) gt 1 then plotwave = dindgen(ceil(max(plotwave,min=mn)-mn)+1)+mn
        nwave = double(n_elements(plotwave))
        for i=0L,nscans-1L do begin
            if skip17 eq 0 then begin
                p=where(alldata.spect17[i].wavelength gt 0.0,count)
                if count ge 3 then begin
                    irrad = interpol(alldata.spect17[i].irradiance[p],alldata.spect17[i].wavelength[p],plotwave,/lsq)
                    irrad = total(irrad)/nwave
                    ts_17_irrad = [ts_17_irrad,irrad]
                    diff = min(abs(alldata.spect17[i].wavelength - plotwave),pos)
                    ts_17_time = [ts_17_time,alldata.spect17[i].timestamp[pos[0]]]
                endif
            endif
            p=where(alldata.spect18[i].wavelength gt 0.0, count)
            if count lt 3 then continue
            irrad = interpol(alldata.spect18[i].irradiance[p],alldata.spect18[i].wavelength[p],plotwave,/lsq)
            irrad = total(irrad)/nwave
            ts_18_irrad = [ts_18_irrad,irrad]
            diff = min(abs(alldata.spect18[i].wavelength - plotwave),pos)
            ts_18_time = [ts_18_time,alldata.spect18[i].timestamp[pos[0]]]
        endfor
        if skip17 eq 0 then begin
            p=where(ts_17_irrad gt 0.0d and ts_17_irrad lt 10.0d*nwave,count)
            ts_17_irrad=ts_17_irrad[p]
            ts_17_time=ts_17_time[p]
            title='SIM_V17'
            if nwave eq 1  then $
                title+=' wl=['+strtrim(string(min(plotwave),format='(f0.1)'),2)+'-'+strtrim(string(max(plotwave),format='(f0.1)'),2)+']'
            lineplot,gps2sd(ts_17_time/1d6),ts_17_irrad,psym=-4,xtitle='Mission Day',ytitle='Irradiance',title=title,charsize=1.5
        endif
        p=where(ts_18_irrad gt 0.0d and ts_18_irrad lt 10.0d*nwave,count)
        ts_18_irrad=ts_18_irrad[p]
        ts_18_time=ts_18_time[p]
        title='SIM_V18'
        if nwave eq 1  then $
            title+=' wl='+strtrim(string(plotwave,format='(f0.1)'),2) $
        else $
            title+=' wl=['+strtrim(string(min(plotwave),format='(f0.1)'),2)+'-'+strtrim(string(max(plotwave),format='(f0.1)'),2)+']'
        lineplot,gps2sd(ts_18_time/1d6),ts_18_irrad,psym=-4,xtitle='Mission Day',ytitle='Irradiance',title=title,charsize=1.5
        if skip17 eq 0 then begin
            return, {plotwave:plotwave, v17_irrad:ts_17_irrad, v18_irrad:ts_18_irrad, timestamp17:ts_17_time, timestamp18:ts_18_time }
        endif else begin
            return, {plotwave:plotwave, v18_irrad:ts_18_irrad, timestamp18:ts_18_time }
        endelse


    endif

 
    return, nscans

end
