;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Measures and apply a baseline correction to the raw ESR data from the
;   dark measurements for the specified instrumentModeId and time span.
;
; CALLING SEQUENCE:
;   spectra = ESR_BASELINE(t0, t1, instrumentModeId, /missionDays)
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
;      31	SIM_A	ESR   -> not yet implemented
;      32	SIM_B	ESR   -> not yet implemented
;
; OPTIONAL INPUT PARAMETERS:
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
;   A structure with the corrected ESR data.
;
; OPTIONAL OUTPUT PARAMETERS:
;      NONE
;
; EXAMPLE:  
;
; COMMON BLOCKS:
;      NONE
;
; NOTES:
;  You may want to run get_sorce_plan before hand to make sure an activity
;  with observations in the desired mode was performed during the time
;  span of interest.
;
; REVISION HISTORY:
;   Revision: $Id: esr_baseline.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;-
;
function esr_baseline, startTime, stopTime, instrumentModeId, $
         gps=gps, missionDays=missionDays, julianDays=julianDays, $
         verbose=verbose, noplot=noplot

    nrows=0
    ; validate the instrumentModeId
    if n_elements(instrumentModeId) eq 0 then begin
        doc_library,'esr_baseline'
        print,''
        print,'Missing instrumentModeId '
        return,-1
    endif else begin
        ; verify we got the right instrument
        valid_inst = [31,32]
        instrumentModeId = fix(instrumentModeId)
        if (where(valid_inst eq instrumentModeId))[0] lt 0 then  begin
            doc_library,'esr_baseline'
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
            detectorId = 100058L
            shutterId = 102121L
            detectorTemp=100170L
            NOMINAL_DARK = 42000d
        end
        32: begin
            ; SIMB ESR
            instrument='SIMB'
            detectorId = 100048L
            shutterId = 100197L
            detectorTemp=100939L
            NOMINAL_DARK = 35800d
        end
    endcase

    SHT_OPEN = 0
    SHT_CLOSE = 1

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

    if keyword_set(verbose) then $
        print,'Processing '+strtrim(string(instrumentModeId),2)+' from ',gps2sd(t0/1d6),' to ',gps2sd(t1/1d6)

    ; get the list of prism positions for the specified time range (+/- 2 minutes to get darks)
    get_sorce_telemetry, tlm_data, info, t0-120.0d6, t1+120.0d6, /gps, tlmId=[shutterId,detectorId]
    if size(tlm_data,/tname) ne 'STRUCT' then begin
        print,'No data found for the specified instrumentId and time range'
        return,-1
    endif else if (*tlm_data.(0)).n_science eq 0 then begin
        print,'No data found for the specified instrumentId and time range'
        return,-1
    endif

    ; the shutter timing is off by 1.0 second (too early) - add 1 second to shutter
    (*tlm_data.(0)).science.timetag += 1d6

    ; dethin the shutter position and set to ESR data timetag
    shutter_pos = dethin_data((*tlm_data.(0)).science.dn, (*tlm_data.(0)).science.timetag, (*tlm_data.(1)).science.timetag)

    ; we only care about the ESR data from SolarESRMode, SolarIRScan and ESRFullScan1 to ESRFullScan15
    activities = ['SolarESRMode','SolarIRScan','ESRFullScan1','ESRFullScan2','ESRFullScan3','ESRFullScan4',$
                  'ESRFullScan5','ESRFullScan6','ESRFullScan7','ESRFullScan8','ESRFullScan9','ESRFullScan10',$
                  'ESRFullScan11','ESRFullScan12','ESRFullScan13','ESRFullScan14','ESRFullScan15']
    plans = get_sorce_plan(t0,t1,/gps, instrument=instrument, activity=activities)

    if size(plans,/tname) ne 'STRUCT' then begin
        print, 'No valid activities during specified time range'
        return,-1
    endif

    ; join plans that are closer appart than 30 seconds
    if n_elements(plans) gt 1 then begin
        pos = where(abs(plans[1:-1].starttime - plans[0:-2].stoptime) lt 30,count)
        if count gt 0 then begin
            plans[pos].stoptime = plans[pos+1].stoptime
            plans[pos+1].starttime = -1d
            p=where(plans.starttime gt 0d)
            plans=plans[p]
        endif
    endif

    outdata=replicate({timetag:0d, dn:0L, dn_corr:0D}, (*tlm_data.(1)).n_science)
    outdata.timetag = (*tlm_data.(1)).science.timetag
    outdata.dn = (*tlm_data.(1)).science.dn
    outdata.dn_corr = (*tlm_data.(1)).science.dn

    ; for each planned activity, we fit a bspline through the dark data and remove
    ; this baseline the ESR data for this padded time span
    for pl=0L,n_elements(plans)-1 do begin
        pt0 = plans[pl].starttime*1d6 - 200d6
        pt1 = plans[pl].stoptime *1d6 + 200d6
        pl_pos=where(outdata.timetag ge pt0 and outdata.timetag le pt1,allcount)
        dark_pos=where(outdata.timetag ge pt0 and outdata.timetag le pt1 and shutter_pos eq SHT_CLOSE,dark_count)
        if dark_count eq 0 then continue

        ; perform the bspline on the dark data
        ;sset=bspline_iterfit(outdata[dark_pos].timetag,outdata[dark_pos].dn,maxiter=0,requiren=10,bkspace=5)
        ;baseline=bspline_valu(outdata[pl_pos].timetag,sset)
        ;outdata[pl_pos].dn_corr -= (baseline - NOMINAL_DARK)

        ; a global bspline over the whole data doesn't work very well
        ; instead, get the meadian value for each dark period and do a simple spline interpolation
        segpos=where(outdata[dark_pos[1:-1]].timetag - outdata[dark_pos[0:-2]].timetag gt 10d6,count)
        p0=0L
        tempx=[]
        tempy=[]
        for i=0L,count do begin
            if i lt count then p1=segpos[i] else p1=n_elements(dark_pos)-1L
            tempx = [tempx, mean(outdata[dark_pos[p0:p1]].timetag)]
            tempy = [tempy, median(outdata[dark_pos[p0:p1]].dn)]
            p0=p1+1
        endfor
        baseline = interpol(tempy, tempx, outdata[pl_pos].timetag, /quad)
        outdata[pl_pos].dn_corr -= (baseline - NOMINAL_DARK)
    endfor

    if NOT keyword_set(noplot) then begin
        lineplot,gps2sd(outdata.timetag/1d6), outdata.dn, title='ESR data',xtitle='Mission Day', ytitle='DN', charsize=1.25, ptitle='ESR Baseline Corrected Data'
        lineplot,gps2sd(outdata.timetag/1d6), outdata.dn_corr, title='ESR corrected data'
    endif

    return, outdata

end
