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
;   Revision: $Id: int_spectra.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;-
;
function int_spectra, startTime, stopTime, instrumentModeId, wrange, version=version, $
         gps=gps, missionDays=missionDays, julianDays=julianDays, indata=indata, $
         development=development, production=production, noplot=noplot, verbose=verbose

    nrows=0
    ; validate the instrumentModeId
    if n_elements(instrumentModeId) eq 0 or n_elements(wrange) lt 2 then begin
        doc_library,'int_spectra'
        print,''
        print,'Missing instrumentModeId '
        return,-1
    endif else begin
        ; verify we got the right instrument
        valid_inst = [31,32,41,42,43,44,45,46,47,48]
        instrumentModeId = fix(instrumentModeId)
        if (where(valid_inst eq instrumentModeId))[0] lt 0 then  begin
            doc_library,'int_spectra'
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
            ;wave_table='SimProfileIntegral'
            wave_table='Wavelength'
            ;dbUrl='jdbc:sybase:Tds:lasp-db1:4100/SORCE'
            ;dbUrl='jdbc:sybase:Tds:lasp-db1:4100/SORCE'
            dbUrl='jdbc:sybase:Tds:lasp-db1:4100/SORCE_SIM_V20'
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
        query = 'SELECT c.microsecondsSinceGpsEpoch,w.wavelengthRef,c.irradiance FROM SimCalibratedIrradiance as c, '
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

    timerange19 = indata.MICROSECONDSSINCEGPSEPOCH
    irrad19 = indata.IRRADIANCE
    wave19 = indata.WAVELENGTHREF

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
    integrated_irrad = dblarr(count_scans)
    timestamp = dblarr(count_scans)
    p0=0L
    p1=-1L
    wrange=[min(wrange,max=m),m]
    for j=0,count_scans-1L do begin
        ; for the last scan point to end of array
        if keyword_set(verbose) then print,'   '+strtrim(string(j+1),2)+' / '+strtrim(string(count_scans),2)
        p0=p1+1L
        if j eq count then p1=n_elements(wave19)-1L else p1=p[j+1]
        ; only get uniq wavelength points
        pos = lindgen(p1-p0+1L) + p0
        timestamp[j] = timerange19[pos[0]]
        mn = min(wave19[pos],max=mx)
        if mn gt wrange[0] then begin
            print,gps2sd(timestamp[j]/1d6),' min wavelength is larger than requested wrange (no extrapolation)'
            continue
        endif else if mx lt wrange[1] then begin
            print,gps2sd(timestamp[j]/1d6),' max wavelength is smaller than requested wrange (no extrapolation)'
            continue
        endif
        added_irrad=interpol(irrad19[pos], wave19[pos], wrange, /spline)
        ; only keep the data from within the requested wavelength range
        k = where(wave19[pos] ge wrange[0] and wave19[pos] le wrange[1])
        wave=[wave19[pos[k]],wrange]
        irrad=[irrad19[pos[k]],added_irrad]
        s=sort(wave)
        q=uniq(wave[s])
        integrated_irrad[j] = int_tabulated(wave[s[q]], irrad[s[q]])
    endfor

    pos = where(integrated_irrad gt 0.0,count)
    if NOT keyword_set(noplot) then begin
        if count eq 0 then begin
            print,'No data to plot'
        endif else begin
            title='iMode='+strtrim(string(instrumentModeId),2)+' version='+strtrim(string(version),2)+' wrange=['
            title = title+strtrim(string(wrange[0],format='(F0.2)'),2)+', '+strtrim(string(wrange[1],format='(F0.2)'),2)+']'
            lineplot,gps2sd(timestamp[pos]/1d6), integrated_irrad[pos], psym=-4,title=title,xtitle='Mission Day', $
                ytitle='Integrated Irradiance'
        endelse
    endif


    return, {timestamp:timestamp, integrated_irrad:integrated_irrad}
 
end
