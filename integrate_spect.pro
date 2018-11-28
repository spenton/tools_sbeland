;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Integrate the SimCalibratedIrradiance results for the specified
;   timerange, instrument mode and wavelength range.
;
; CALLING SEQUENCE:
;   integrated_spectra = INTEGRATE_SPECT(starttime, endtime, instrumentModeId, wrange)
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
;   wrange - 
;      A 2 element array with the lowest and highest wavelength to integrate.
;
; OPTIONAL INPUT PARAMETERS:
;   development -
;      If specified will use the DEVELOPMENT database
;   production -
;      If specified will use the PRODUCTION database (default)
;   version -
;      Specifies the processed data version from the database. If not specified,
;      will use the highest numbered version.
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
;   A structure with the integrated spectra. The wavelength scale will be 
;   taken from the first spectra or from the provided wavelength input array.
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
;   Revision: $Id: integrate_spect.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;-
;
function integrate_spect, startTime, stopTime, instrumentModeId, wrange, version=version, $
         gps=gps, missionDays=missionDays, julianDays=julianDays, plans=plans, indata=indata, $
         development=development, production=production, noplot=noplot, verbose=verbose, dayavg=dayavg

    nrows=0
    ; validate the instrumentModeId
    if n_elements(instrumentModeId) eq 0 or n_elements(wrange) lt 2 then begin
        doc_library,'integrate_spect'
        print,''
        print,'Missing instrumentModeId '
        return,-1
    endif else begin
        ; verify we got the right instrument
        valid_inst = [31,32,41,42,43,44,45,46,47,48]
        instrumentModeId = fix(instrumentModeId)
        if (where(valid_inst eq instrumentModeId))[0] lt 0 then  begin
            doc_library,'integrate_spect'
            print,''
            print,'Invalid instrumentModeId was provided'
            return,-1
        endif
    endelse

    ; define the telemetry items according to the specified mode
    case instrumentModeId of
        31: begin
            ; SimA ESR
            instrument='SimA'
        end
        41: begin
            ; SimA VIS1
            instrument='SimA'
        end
        42: begin
            ; SimA VIS2
            instrument='SimA'
        end
        43: begin
            ; SimA UV
            instrument='SimA'
        end
        44: begin
            ; SimA IR
            instrument='SimA'
        end
        32: begin
            ; SimB ESR
            instrument='SimB'
        end
        45: begin
            ; SimB VIS1
            instrument='SimB'
        end
        46: begin
            ; SimB VIS2
            instrument='SimB'
        end
        47: begin
            ; SimB UV
            instrument='SimB'
        end
        48: begin
            ; SimB IR
            instrument='SimB'
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

    if keyword_set(verbose) then begin
        print,'Processing mode='+strtrim(string(instrumentmodeid),2)+' from '+$
        strtrim(string(gps2sd(t0/1d6),format='(F0.2)'),2)+' to '+strtrim(string(gps2sd(t1/1d6),format='(F0.2)'),2)+$
        ' covering wavelengths ['+strtrim(string(wrange[0]),2)+', '+strtrim(string(wrange[1]),2)+']'
    endif

    if size(indata,/tname) ne 'STRUCT' then begin
        if keyword_set(development) then begin
            if keyword_set(verbose) then print,'   using the development database ...'
            wave_table='Wavelength'
            dbUrl='jdbc:sybase:Tds:lasp-db-dev:4100/SORCE_SIM_V19'
            dbDriver='com.sybase.jdbc3.jdbc.SybDriver'
            user='sbeland'
            password='sbeland1'
            jstmt = fjava_get_jdbc_statement(user=user, password=password, dbdriver=dbdriver, dburl=dburl) 
            oJavaDbExchange = OBJ_NEW("oJava_DbExchange", jstmt)
        endif else begin
            if keyword_set(verbose) then print,'   using the production database ...'
            wave_table='SimProfileIntegral'
            dbUrl='jdbc:sybase:Tds:lasp-db1:4100/SORCE'
            dbDriver='com.sybase.jdbc3.jdbc.SybDriver'
            user='sorce'
            password='sorcedb'
            jstmt = fjava_get_jdbc_statement(user=user, password=password, dbdriver=dbdriver, dburl=dburl) 
            oJavaDbExchange = OBJ_NEW("oJava_DbExchange", jstmt)
        endelse

        if n_elements(version) eq 0 then begin
            if keyword_set(verbose) then print,'   looking for version to extract ...'
            query = 'SELECT MAX(version) version from SimCalibratedIrradiance where instrumentModeId='+strtrim(string(instrumentModeId),2)
            query = query + ' and microsecondsSinceGpsEpoch>='+strtrim(string(ulong64(t0)),2)
            query = query + ' and microsecondsSinceGpsEpoch<='+strtrim(string(ulong64(t1)),2)
            res=oJavaDbExchange->getAllValues(query)
            if size(res,/tname) ne 'STRUCT' then begin
                print,'Error: no max(version) for SimCalibratedIrradiance found for the requested date/version'
                return,-1
            endif
            version = res[0].version
        endif
        if keyword_set(verbose) then print,'   extracting data from version='+strtrim(string(version),2)+' ...'

        ; get the spectra for this version
        query = 'SELECT c.microsecondsSinceGpsEpoch as timestamp,w.wavelengthRef as wavelength,c.irradiance FROM SimCalibratedIrradiance as c, '
        query += wave_table+' as w where  c.instrumentModeId='+strtrim(string(instrumentModeId),2)
        query += ' and w.instrumentModeId='+strtrim(string(instrumentModeId),2)
        query += ' and w.version='+strtrim(string(version),2)
        query += ' and c.version='+strtrim(string(version),2)
        query += ' and c.microsecondsSinceGpsEpoch>='+strtrim(string(ulong64(t0)),2)
        query += ' and c.microsecondsSinceGpsEpoch<='+strtrim(string(ulong64(t1)),2)
        query += ' and c.microsecondsSinceGpsEpoch=w.microsecondsSinceGpsEpoch order by c.microsecondsSinceGpsEpoch'
        ; print,query
        indata=oJavaDbExchange->getAllValues(query)
        if size(indata,/tname) ne 'STRUCT' then begin
            print,'Error: no SimCalibratedIrradiance found for the requested date/version'
            return,-1
        endif
    endif

    ; keep only the data within our wavelength range (plus some)
    wrange=[min(wrange,max=m),m]
    p=where(indata.WAVELENGTH ge wrange[0]-10d and indata.WAVELENGTH le wrange[1]+10d,count)
    if count eq 0 then begin
        print,'No data withing requested wavelength range: ',wrange
        return,-1
    endif

    s=sort(indata[p].TIMESTAMP)
    timerange19 = indata[p[s]].TIMESTAMP
    irrad19 = indata[p[s]].IRRADIANCE
    wave19 = indata[p[s]].WAVELENGTH
    nwave=n_elements(wave19)

    if instrumentModeId ne 31 and instrumentModeId ne 32 then begin
        ; for the photo-diodes, use the SolarQuickScan24 as our indicator
        if n_elements(plans) eq 0 then $
            plans=get_sorce_plan(t0,t1,/gps,instrument=instrument,activity='SolarQuickScan24')
        pl0 = plans.starttime * 1d6
        pl1 = plans.stoptime * 1d6
        p0=[]
        p1=[]
        q=0L
        ; get the positions in the large array for each scan
        for i=0,n_elements(plans)-1L do begin
            p=where(timerange19 ge pl0[i] and timerange19 le pl1[i],count)
            if count eq 0 then continue
            p0=[p0,p[0]]
            p1=[p1,p[-1]]
            if p1[-1] ge nwave then break
        endfor
    endif else begin
        ; organize the data as an array of structure organize per scan
        ; assuming each spectra has same number of elements
        ; p=where(wave19[0:-2]-wave19[1:-1] gt 0.0d,count)
        ; look for gaps of 240 seconds or more 
        p = where((timerange19[1:-1] - timerange19[0:-2]) gt 240.0d6,count)
        ; count the number of points per spectra
        if count eq 0 then begin
            p0=0L
            p1=n_elements(wave19)-1L
        endif else begin
            p0=[0L,p]
            p1=[p+1L,n_elements(wave19)-1]
        endelse
    endelse
    count_scans = n_elements(p0)
    integrated_irrad = dblarr(count_scans)
    timestamp = dblarr(count_scans)


    for j=0,count_scans-1L do begin
        ; for the last scan point to end of array
        if keyword_set(verbose) then print,'   '+strtrim(string(j+1),2)+' / '+strtrim(string(count_scans),2)
        ; only get uniq wavelength points
        timestamp[j] = timerange19[p0[j]]
        pos=lindgen(p1[j]-p0[j]+1)+p0[j]
        mn = min(wave19[pos],max=mx)
        if mn gt wrange[0] then begin
            print,gps2sd(timestamp[j]/1d6),' min wavelength is larger than requested wrange (no extrapolation)  '+$
                '['+strtrim(string(mn,format='(F0.2)'),2)+', '+strtrim(string(mx,format='(F0.2)'),2)+']'
            continue
        endif else if mx lt wrange[1] then begin
            print,gps2sd(timestamp[j]/1d6),' max wavelength is smaller than requested wrange (no extrapolation)  '+$
                '['+strtrim(string(mn,format='(F0.2)'),2)+', '+strtrim(string(mx,format='(F0.2)'),2)+']'
            continue
        endif
        s=sort(wave19[pos])
        added_irrad=interpol(irrad19[pos[s]], wave19[pos[s]], wrange)
        wave=[wrange[0],wave19[pos[s]],wrange[1]]
        irrad=[added_irrad[0],irrad19[pos[s]],added_irrad[1]]
        q=uniq(wave)
        p=where(wave[q] ge wrange[0] and wave[q] le wrange[1],count)
        if count gt 0 then integrated_irrad[j] = int_tabulated(wave[q[p]], irrad[q[p]])
    endfor

    resistant_mean,integrated_irrad,4.0,mean,goodvec=keep0
    integrated_irrad=integrated_irrad[keep0]
    timestamp=timestamp[keep0]
    sd = gps2sd(timestamp/1d6)

    if keyword_set(dayavg) then begin
        minh = floor(min(sd,max=mx))
        maxh = floor(mx+1.0)
        hist=histogram(sd,binsize=1.0,min=minh,max=maxh,location=day,reverse=rev)
        new_irrad=dblarr(n_elements(day))
        for i=0,n_elements(day)-1 do begin
            if rev[i] ne rev[i+1] then begin
                new_irrad[i]=mean(INTEGRATED_IRRAD[rev[rev[i]:rev[i+1]-1]])
            endif
        endfor
        p=where(new_irrad gt 0.0,count)
        if count gt 0 then begin
            new_irrad=new_irrad[p]
            day=day[p]
        endif
        resistant_mean,new_irrad,4.0,mean,goodvec=keep0
        integrated_irrad=new_irrad[keep0]
        sd=day[keep0]+0.5d
        timestamp=sd2gps(sd)*1d6
    endif

    if NOT keyword_set(noplot) then begin
        title='iMode='+strtrim(string(instrumentModeId),2)+' version='+strtrim(string(version),2)+' wrange=['
        title = title+strtrim(string(wrange[0],format='(F0.2)'),2)+', '+strtrim(string(wrange[1],format='(F0.2)'),2)+']'
        if keyword_set(dayavg) then title=title+' Daily Average'
        lineplot,sd, integrated_irrad, psym=-4,title=title,xtitle='Mission Day', ytitle='Integrated Irradiance'
    endif


    return, {timesd:sd, timestamp:timestamp, integrated_irrad:integrated_irrad}
 
end
