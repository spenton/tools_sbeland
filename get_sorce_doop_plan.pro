; docformat = 'rst'
;+
;   Return a structure containing information about all SimA
;   and/or SimB activities planned in the time range specified.
;
; :Author:
;   S. Beland
;   
; :Examples:
; ::  
;   plans = get_sorce_plan(453,500,/mission,/sima) 
;   ; get all activities for SimA from SORCE mission day 453 to 500 
;         
;   plans = get_sorce_plan(453,500, /sima, activity='SolarQuickScan24')
;   ; only return the list of SolarQuickScan24 activities for SimA   
;     
;   plans = get_sorce_plan(453,500, instrument='Tim')
;   ; get all the activities for the TIM instrument between day 453 and 500
;      
;   plans = get_sorce_plan(453,500, instrument='SolsticeA', activity='MUVSolarMiniScan')
;   ; get the list of MUVSolarMiniScan for the TIM instrument between day 453 and 500
;      
;   plans = get_sorce_plan(453,470,/sima,act='ESRFullScan%')
;   ; get the list of plans for SimA with activity name matching ESRFullScan*
;     
;   plans = get_sorce_plan(453,470,/sima,act='ESRFullScan1[3-5]%')
;   ; get the list of plans for SimA with activity names ESRFullScan13, ESRFullScan14, ESRFullScan15
;      
;   ; get the list of FUVStellar scans where the targetName was the 'moon'
;   plans = get_sorce_plan(julday(1,1,2006), julday(1,1,2012), /jul, instrument='SolsticeB', $
;                          activity='FUVStellarScan', paramName='targetName', paramValue='moon')
;   ; and then get the ABComparison activities within 5 days to one of these observations
;   solar_plan = get_sorce_plan(plans[k].starttime-5.0, plan2[k].starttime+5.0, /julianDays, $
;                          instrument='SolsticeB', activity='ABComparison')
;
;   ; get the list of "Safe Mode Config" commanded
;   pl=get_sorce_plan(4227d,4228d,/mission,instrument='SORCE',activity='Safe Mode Config’)
;      
;   result = GET_SORCE_PLAN(t0, t1, /simA, /mission, activity="SolarAlignment5", version="1.10") 
;-

;+
; :Params:
;   startTime : in, required
;      The lower time range for which data will be returned.
;      May be specified as Julian days, gps microseconds, 
;      VMS time string or SORCE mission day numbers (the default).
;
;   stopTime : in, required
;      The upper time range for which data will be returned.
;      May be specified as Julian days, gps microseconds, 
;      VMS time string or SORCE mission day numbers (the default).
;
; :Keywords:
;   simA - in, optional
;      Specifies which SIM channel to retrieve. Zero, one or both 
;      instruments can be requested. If neither are selected, it defaults
;      to both simA and simB.
;      
;   simB - in, optional
;      Specifies which SIM channel to retrieve. Zero, one or both
;      instruments can be requested. If neither are selected, it defaults
;      to both simA and simB.
;         
;   instrument - in, optional
;      Specifies the exact instrument to extract the plan: SIMA or SIMB
;      
;   gps - in, optional
;      If set, indicates that input times are to be interpreted as
;      microseconds since Jan 6, 1980 midnight UT.  (default)
;      
;   missionDays - in, required
;      If set, indicates that input times are to be interpreted as
;      days elapsed since launch, i.e. mission days.  If only the startTime
;      argument contains a non-zero value (stopTime is not supplied), then
;      all data for the specified mission day number will be retrieved.
;      
;   julianDays - in, optional
;      If set, indicates that the input times are to be interpreted as
;      julian day numbers.
;      
;   noconvert - in, optional
;      If set, will NOT convert the returned start and stop times to the same
;      format as the input starttime and stoptime and leave them in the default
;      format returned from the planning database (Sybase date/time format).
;      
; :Returns:
;   An array of structures whose fields correspond to the returned
;   database column names. The format for the returned timestamps
;   will be the same as the input timestamp.
;
; :Notes:
;   In DOOP mode we can only identify SolarQuickScan24, SolarESRMode and SolarIRScan 
;   at this time.  These are identified by SQS24, esrScan and irScan.
;-
function get_sorce_doop_plan, starttime, stoptime, activity=activity, sima=sima, $
             simb=simb, gps=gps, missionDays=missionDays, vms=vms, $
             julianDays=julianDays, noconvert=noconvert,verbose=verbose

    ; will extract data for SIM_A by default
    instrument='SimA'
    if keyword_set(simA) then instrument='SimA'
    if keyword_set(simB) then instrument='SimB'

    if instrument eq 'SimA' then begin
        sima=1
        simb=0
        externalElement='sim_a'
        scienceTable='SimAScienceSamples'
    endif else begin
        externalElement='sim_b'
        sima=0
        simb=1
        scienceTable='SimBScienceSamples'
    endelse

    if keyword_set(vms) then begin
        ; user specified timetags as VMS strings 
        t0 = (string(startTime))
        t1 = (string(stopTime))
    endif else if keyword_set(julianDays) then begin
        ; user specified time in julian days
        t0 = (jd2vms(startTime))
        t1 = (jd2vms(stopTime))
    endif else if keyword_set(gps) then begin
        ; user specified timetags in gps microseconds
        t0 = (gps2vms(startTime/1d6))
        t1 = (gps2vms(stopTime/1d6))
    endif else begin
        ; user specified time in mission (sorce) days
        missionDays=1
        t0 = (sd2vms(startTime))
        t1 = (sd2vms(stopTime))
    endelse
    startSD=vms2sd(t0)
    stopSD=vms2sd(t1)
 
    if n_elements(version) gt 0 then version=strtrim(version,2)

    ; since we can't rely on the planning database in hybrid mode, figure out which experiment
    ; was run during each orbit
    ; pad the start and stop time to make sure we get at least a complete orbit

    if keyword_set(verbose) then print,'   extracting Solar Periods ...'

    solarPeriod=get_solar_period(startSD-0.0625, stopSD+0.0625, /mission)
    keep=where(solarperiod.stoptime_gps ge sd2gps(startSD)*1d6 and solarperiod.starttime_gps le sd2gps(stopSD)*1d6,count)
    plans=[]
    if keyword_set(verbose) then print,'   determine the experiment for the '+strtrim(string(count),2)+' orbits found ...'
    for i=0,count-1 do begin
        ; since we only care about the prismposition, avoid the overhead of get_sorce_telemetry by calling query_database directly
        ;get_sorce_telemetry, prismpos, info, solarperiod[keep[i]].starttime_gps, solarperiod[keep[i]].stoptime_gps, $
        ;    itemName='position', externalElement=externalElement
        sql = "select sampleVtcw 'timetag',prismPosition 'eu' from "+scienceTable+" where "
        sql += " sampleVtcw >= "+strtrim(string(solarperiod[keep[i]].starttime_gps),2)
        sql += " and sampleVtcw < "+strtrim(string(solarperiod[keep[i]].stoptime_gps),2)
        sql += " order by sampleVtcw"
        query_database, sql, prismpos, nrows, database = "SORCE_L1S"
        if nrows gt 0 then begin
            ; we found some telemetry for this orbit - figure out which activity was run
            act = determine_experiment(prismpos.eu, sima=sima, simb=simb)
            if size(act,/tname) eq 'STRUCT' then begin
                keepit=0
                if n_elements(activity) then begin
                    if (act[0].type).contains('SQS24') and activity.contains('solarquickscan',/fold) then keepit=1
                    if (act[0].type).contains('esrScan') and activity.contains('SolarESRMode',/fold) then keepit=1
                    if (act[0].type).contains('irScan') and activity.contains('SolarIRScan',/fold) then keepit=1
                endif else keepit=1
                if keepit then begin
                    plan = {starttime:gps2sd(solarperiod[keep[i]].starttime_gps/1d6), stoptime:gps2sd(solarperiod[keep[i]].stoptime_gps/1d6),$
                        ACTIVITYNAME:act[0].type, PLANNINGNAME:instrument, COMMANDVERSION:''}
                    plans = [plans, plan]
                endif
            endif
        endif
    endfor

    nplans=n_elements(plans)
    if (nplans eq 0) then return, 0L

    if keyword_set(verbose) then print,'   found '+strtrim(string(nplans),2)+' matching plans'

    ; convert the names to the standard one
    p=where(plans.activityname eq 'SQS24A' or plans.activityname eq 'SQS24B',count)
    if count gt 0 then plans[p].activityname = 'SolarQuickScan24'

    p=where(plans.activityname eq 'esrScanA' or plans.activityname eq 'esrScanB', count)
    if count gt 0 then plans[p].activityname = 'SolarESRMode'

    p=where(plans.activityname eq 'irScanA' or plans.activityname eq 'irScanB',count)
    if count gt 0 then plans[p].activityname = 'SolarIRScan'


    ; convert the time formats
    if keyword_set(noconvert) then begin
        plans.starttime=sd2gps(plans.starttime)*1d6
        plans.stoptime=sd2gps(plans.stoptime)*1d6
    endif else begin
        if keyword_set(julianDays) then begin
            ; user specified time in julian days
            plans.starttime=sd2jd(jdbc2vms(plans.starttime))
            plans.stoptime=sd2jd(jdbc2vms(plans.stoptime))
        endif else if keyword_set(gps) then begin
            ; user specified timetags in gps microseconds
            plans.starttime=sd2gps(jdbc2vms(plans.starttime))
            plans.stoptime=sd2gps(jdbc2vms(plans.stoptime))
        endif else if keyword_set(vms) then begin
            ; user specified timetags in gps microseconds
            plans.starttime=sd2vms(jdbc2vms(plans.starttime))
            plans.stoptime=sd2vms(jdbc2vms(plans.stoptime))
        endif
    endelse

    return, plans

end
