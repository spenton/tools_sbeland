;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Calculate the solar exposure within the specified time range for
;   the requested instrument (SIMA or SIMB) speicifically for the
;   HYBRID mode.  This mode will very often have gaps in the telemetry
;   resulting from intermittant TDRS contacts. 
;   In HYBRID mode, a single experiment is scheduled per orbit and repeated
;   until a Safe Mode is executed which will force the shutter closed in
;   preparation to entering the eclipse portion of the orbit.
;   The planning database will contain the scheduled experiment (with many
;   potential start time), and the SolarPeriod corresponding to the time
;   ranges when the sun is visible, as well as the times when Safe Mode
;   is triggered.  
;
; CALLING SEQUENCE:
;   out_data = CALC_SIM_EXP_hybrid(startTime, stopTime, simA=simA, simB=simB)
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
; OPTIONAL INPUT KEYWORDS:
;   simA, simB - 
;      Specifies which SIM channel to retreive.  These are mutually
;      exclusive and defaults to simA.
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
; OUTPUT PARAMETERS:
;   out_data -
;      The resulting solar exposure time for the specified instrument calculated
;      per orbit.  The structure will have the following format:
;         instrument:  string - (sim_a or sim_b)
;         boxwidth:    float - (0.0 indicates boxwidth was not defined)
;         hrtpos:      string - ('IN', 'OUT', or 'NA')
;         starttime:   double array - gps timestamp (microsecs) of start of each orbit (when entering day)
;         endtime:     double array - gps timestamp (microsecs) of end of each orbit (when entering night)
;         insun_time:  double array - total time in_sun for each orbit (seconds)
;         gap_sci:     double array - total time with missing science telemetry (seconds) (gaps larger than 2 sec)
;         gap_sc:      double array - total time with missing spacecraft telemetry (seconds) (fssunangle0)
;         gap_sc_open: double array - total time with missing spacecraft telemetry with shutter open
;         gap_hrt:     double array - total time with missing HRT tlm (when power_on and shutter_open)
;         same_hrt:    double array - total time with HRT_IN=HRT_OUT (seconds)
;         power_off:   double array - total time during orbit when instrument's power is off
;         hrt_from_plan:double array - total time when HRT state was obtained from the plan database 
;                      (power_on and shutter_open)
;         safehold:    double array - total time during orbit when SC is in safehold/contingency mode
;         abs_gap:     double array - total absolute gap time (combined gaps per time slice)
;                      This data only reflects the times when we were not able to determine the
;                      state of the HRT, SHUTTER or found gaps in the spacecraft telemetry.
;                      These gaps come from:
;                        - missing science data for more than 2.0 seconds (and not in safehold),
;                        - missing SC HK (no fssunagle tlm) when not in safehold with power_on and shutter_open
;                        - HRT indicates both IN and OUT
;         solar_exp:   double array - total amount of solar exposure time (seconds)
;
; OPTIONAL OUTPUT PARAMETERS:
;      NONE
;
; EXAMPLE:  
;
; COMMON BLOCKS:
;      NONE
;
; NOTE:
;    This program returns the solar exposure per orbit only. The orbits are defined in
;    SolarPeriod.  The orbits for which the solar exposure will be calculated are the 
;    ones starting within the specified time range. If the specified starttime is in 
;    the middle of a solar period, the following orbit will defined the actual start.
;    The end of the time period considered will correspond to the end of the last 
;    orbit which started within the time range requested.  This way we only include
;    complete orbits.
;    In the Hybrid mode, only a single type of experiment can be run in a single
;    orbit.
;
; REVISION HISTORY:
;   2014.11.30  SBeland
;   Revision: $Id: calc_sim_exp_hybrid.pro,v 1.2 2018/11/26 16:16:02 spenton Exp $
;-
;*****************************************************************************
FUNCTION calc_sim_exp_hybrid,start_time,stop_time, sima=sima, simb=simb, verbose=verbose, $
                       missionDays=missionDays, gps=gps, julianDays=julianDays, $
                       user=user, password=password, dbdriver=dbdriver, dbUrl=dbUrl

; do sima by default
if keyword_set(sima) then begin
    sima = 1
    simb = 0
    externalElement='sim_a'
    instrument = "SimA"
endif else if keyword_set(simb) then begin
    sima = 0
    simb = 1
    externalElement='sim_b'
    instrument = "SimB"
endif else begin
    sima = 1
    simb = 0
    externalElement='sim_a'
    instrument = "SimA"
endelse

if keyword_Set(gps) then begin
    sd0 = gps2sd(start_time/1d6)
    sd1 = gps2sd(stop_time/1d6)
endif else if keyword_set(julianDays) then begin
    sd0 = jd2sd(start_time)
    sd1 = jd2sd(stop_time)
endif else begin
    sd0=start_time
    sd1=stop_time
endelse

; seconds when shutter is closed between 2 consecutive SolarQuickScan24
sqs24_between_scans = 211.0d
; delay from start of solar orbit to start of first scan
sqs24_delay_scans = 1200.0d

;Get data on when the solar-period is for SORCE - this defines the part of each orbit when spacecraft sees the sun.
if keyword_set(verbose) then print,' extracting the SolarPeriod ...'
; add 90 minutes to the end to make sure we get a complete orbit at the end
query_database,/reset
solarPeriod = get_sorce_plan(sd0, sd1+1.5/24.0, activity='SolarPeriod', /events, user=user, password=password, dbUrl=dbUrl, dbdriver=dbdriver)

solarPeriod_gps0=sd2gps(solarperiod.starttime)*1d6   ; gps micro-seconds
solarPeriod_gps1=sd2gps(solarperiod.stoptime)*1d6    ; gps micro-seconds

; only keep the orbits that START within the requested time range
keep = where(solarPeriod_gps0 ge sd2gps(sd0)*1d6 and solarPeriod_gps0 lt sd2gps(sd1)*1d6, norbits)
if norbits eq 0 then begin
    print,'  Warning:  no orbit found within requested time range'
    return,-1
endif
solarPeriod = solarPeriod[keep]
solarPeriod_gps0=SolarPeriod_gps0[keep]
solarPeriod_gps1=SolarPeriod_gps1[keep]

; initialize the output structure - keep the same data structure BEFORE HYBRID for compatibility
tmpstr = {instrument:instrument, boxwidth:0.0, hrtpos:'OUT', $
         starttime:0.0d, endtime:0.0d, insun_time:0.0d, $
         gap_sci:0.0d, gap_sc:0.0d, gap_sc_open:0.0d, $
         gap_hrt:0.0d, same_hrt:0.0d, $
         gap_power:0.0d, power_off:0.0d, $
         abs_gap:0.0d, hrt_from_plan:0.0d, $
         time_from_plan:0.0d, time_std_from_plan:0.0d, $
         safehold:0.0d, solar_exp:0.0d } 

out_data = replicate(tmpstr,norbits)
out_data.starttime = solarPeriod_gps0
out_data.endtime = solarPeriod_gps1
out_data.insun_time = (solarPeriod_gps1-solarPeriod_gps0)/1.0d6

items = ['position', 'shutter_pos']     

; process only a subset of orbits (~2 days) at a time to limit the number and size of queries
orbs2proc = 30.0
ngrp = ceil(float(norbits)/orbs2proc)

for orbgrp=0L,ngrp-1L do begin
   
    orb0 = orbgrp*orbs2proc
    orb1 = ((orbgrp+1)*orbs2proc-1L) < (norbits-1L)
    ; stretch the start and end of sunpresence by 10 minutes to cover delay in telemetry
    orb_start = solarPeriod_gps0[orb0] - 600.0d6
    orb_end   = solarPeriod_gps1[orb1] + 600.0d6

    if keyword_set(verbose) then print,' extracting the telemetry for orbits '+strtrim(string(orb0,format='(I)'),2)+$
        ' to '+strtrim(string(orb1,format='(I)'),2)+' ...'
    ; get the telemetry for this chunk of time
    query_database,/reset
    get_sorce_telemetry, sim_data, info, externalElement=externalElement, itemName=items, orb_start, orb_end,/gps
    ; make sure we have some telemetry data to look at
    if (*sim_data.(0)).n_science eq 0 then begin
        ; no telemetry at all -> we can't tell anything about the solarExp
        print,'  WARNING: no telemetry was found for the requested timerange'
        continue
    endif

    ; get the times when the spacecraft was commanded to safe-mode, forcing the shutter to close
    if keyword_set(verbose) then print,' extracting the Safe Mode Config plans ...'
    safe_mode_plan  = get_sorce_plan(orb_start,orb_end,/gps,instrument='SORCE',activity="Safe Mode Config", $
        user=user, password=password, dbUrl=dbUrl, dbdriver=dbdriver)
    ; check if the orbit info was available (if not we can't do anytghing - not in a valid part of hybrid mode)
    if size(SAFE_MODE_PLAN,/tname) eq 'STRUCT' then begin
        safe_mode_plan.starttime *= 1d6   ; gps micro-seconds
        safe_mode_plan.stoptime *= 1d6    ; gps micro-seconds
    endif else begin
        ; with no safe_mode entry in the database - simply assume the orbit ends with solarPeriod
        safe_mode_plan = replicate({starttime:0d, stoptime:0d}, n_elements(solarPeriod_gps1))
        ; when the safe_mode_plan is available, it is schedule 5 minutes before the end of solarPeriod
        ; but early in the hybrid mode (when safe_mode_plan is absent) - it's unclear
        safe_mode_plan.starttime = solarPeriod_gps1 - 300d6
        safe_mode_plan.stoptime  = solarPeriod_gps1 
    endelse

    ; now that we have extracted the data from the database, process one orbit at a time
    for orb=orb0, orb1 do begin
        if keyword_set(verbose) then print,'   processing orbit ('+strtrim(string(orb+1,format='(I)'),2)+') starting at '+$
            jd2iso(gps2jd(solarPeriod_gps0[orb]/1d6))
        pos = where((*sim_data.(0)).science.timetag ge solarPeriod_gps0[orb] and  $
                    (*sim_data.(0)).science.timetag le solarPeriod_gps1[orb],count)
        ; if no telemetry in this orbit - assume no experiments were scheduled and shutter was closed
        if count eq 0 then continue

        ; assume shutter closes 4.0 seconds after first safeMode if present
        k=where(safe_mode_plan.stoptime gt solarPeriod_gps0[orb] and safe_mode_plan.stoptime le solarPeriod_gps1[orb],count)
        if count gt 0 then begin
            end_orbit=safe_mode_plan[k[0]].starttime + 4d6
        endif else begin
            ; assume the safe mode command was sent 5 minutes before entering eclipse (like all other times)
            end_orbit=solarPeriod_gps1[orb] - 300d6
        endelse 

        ; in some cases the fisrt safe_mode is not executed (they are repeated until one catches)
        ; look for the last telemetry with shutter open and use that if larger than end_orbit
        ; FOR SOME WEIRD REASON, POSITION is an array of STRUCT while SHUTTER_POS is a structure of arrays ???
        p=where((*sim_data.(1)).science.timetag[pos] le solarPeriod_gps1[orb] and (*sim_data.(1)).science.eu[pos] eq 0)
        end_orbit = max( [((*sim_data.(1)).science.timetag[pos[p[-1]]])+2d6, end_orbit] )

        time_pos = (*sim_data.(0)).science[pos].timetag
        positions = (*sim_data.(0)).science[pos].eu
        ; shutter telemetry is synchronized with prism telemetry 
        ; but returned structure is different (is this because of IdlDBInterface 2.4?)
        shutter = (*sim_data.(1)).science.eu[pos]
        
        ;determine which experiment was run during this orbit
        experiment = determine_experiment(positions, sima=sima, simb=simb, validpos=validpos)
        ;verify we have data to process
        if size(experiment,/tname) ne 'STRUCT' then continue

        if keyword_set(verbose) then print,'   experiment: ',experiment.type
        if strpos(strupcase(experiment.type), 'SQS24') ge 0 then begin
            ; processing SolarQuickScan24 data
            ; get the largest prismposition and figure out the actual exposure
            maxexp=max(experiment.data.time)
            q = where(experiment.data.position eq max(validpos))
            cumul_time=experiment.data[q].time 
            ; see if the beginning of this experiment is after the start of the orbit
            k=(where(positions eq max(validpos)))[0]
            start_scan = time_pos[k] - cumul_time*1d6
            end_scan = start_scan + maxexp *1d6
            ; determine if this is the first or second scanfor this orbit 
            time0 = start_scan - (experiment.time_between_scans + maxexp + experiment.time_delay_scans)*1d6
            if time0 gt solarPeriod_gps0[orb] then begin
                ; we are dealing with the second scan - assume a complete first scan was fully performed
                out_data[orb].solar_exp += maxexp
                ; count the exposure time for this scan up to the last prism position
                out_data[orb].solar_exp += cumul_time
                out_data[orb].solar_exp += (end_orbit - time_pos[k])/1d6
            endif else begin
                ; we are dealing with the first scan
                ; make sure it ends before safeMode
                if start_scan + maxexp*1d6 lt end_orbit then begin
                    out_data[orb].solar_exp += maxexp
                    ; account for the second scan (if any)
                    start_scan += (maxexp + experiment.time_between_scans)*1d6
                    end_scan = (start_scan + maxexp*1d6) < (end_orbit)
                    ; check that telemetry (if any) agrees with a second scan
                    p=where(time_pos ge start_scan and time_pos le end_scan,npos)
                    if npos gt 0 then begin
                        q=where(shutter[p] eq 0,nopen, comp=nclose)
                        ; look for closed shutter before end_scan
                        p0=where(time_pos[p[nclose]] gt time_pos[p[q[-1]]],c)
                        if c gt 0 then end_scan = time_pos[p[nclose[p0[0]]]]
                        if nopen gt 0 then out_data[orb].solar_exp += (end_scan - start_scan)/1d6
                    endif else begin
                        ; we have no telemetry -> assume a second scan
                        out_data[orb].solar_exp += (end_scan - start_scan)/1d6
                    endelse
                endif else begin
                    end_scan = (start_scan + maxexp*1d6) < (end_orbit)
                    out_data[orb].solar_exp += (end_scan - start_scan)/1d6
                endelse
            endelse

        endif else if strpos(strupcase(experiment.type), 'IRSCAN') ge 0 then begin
            ; processing SolarIRScan data
            ; get the last prism transition and figure out the actual exposure
            maxexp=max(experiment.data.time)
            diff = positions[1:-1] - positions[0:-2]
            cycle_time = experiment.data[-1].time - experiment.data[-2].time
            p=where(diff eq experiment.data[1].position-experiment.data[0].position, count)
            if count gt 0 then begin
                q = where(experiment.data.position eq positions[p[-1]])
                cumul_time=experiment.data[q].time 
                start_scan = time_pos[p[-1]] - cumul_time * 2d6  ; to account for the times with shutter closed
                end_scan = start_scan + maxexp *1d6
                ; count the time for complete cycles
                out_data[orb].solar_exp += cumul_time
                ; now count the time from the last fractional cycle
                out_data[orb].solar_exp += ((end_orbit - time_pos[p[-1]]) / 2d6)
            endif else begin
            endelse

        endif else if strpos(strupcase(experiment.type), 'ESRSCAN') ge 0 then begin
            ; processing SolarESRMode1-8 data
            cycle_time = experiment.data[-1].time - experiment.data[-2].time
            ; look for a transition in prism position within about 2 second
            pos=where(abs(positions[1:-1] - positions[0:-2]) gt 0 and positions[1:-1] gt 0 and $
                      (time_pos[1:-1] - time_pos[0:-2]) lt 2d6, count)
            if count gt 0 then begin
                ; figure out how much exposure up to that transition
                pos1=where(experiment.data.position eq positions[pos[0]+1])
                out_data[orb].solar_exp += experiment.data[pos1[0]].time
                ; now calculate the number of full and partial cycles up to end_orbit
                t1 = time_pos[pos[0]+1]
                num_full_cycle = floor((end_orbit - t1) / (cycle_time * 2d6))
                out_data[orb].solar_exp += (num_full_cycle * cycle_time)
                t1 += num_full_cycle * cycle_time * 2d6
                ; figure out the remainder
                fraction = floor((end_orbit - t1)/1d6 / cycle_time)
                if fraction eq 0 then begin
                    out_data[orb].solar_exp += ((end_orbit - t1)/1d6 < (50.0001d))
                endif else begin
                    out_data[orb].solar_exp += fraction * cycle_time / 2d
                    remain = (end_orbit - t1)/1d6 - cycle_time
                    if remain gt 0.0 then out_data[orb].solar_exp += (remain<(cycle_time/2d))
                endelse
            endif

        endif

    endfor

endfor

return, out_data

END
