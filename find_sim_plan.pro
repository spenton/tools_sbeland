;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Return a structure containing information about all SimA
;   and/or SimB activities planned in the time range specified.
;   Instaed of querying the planning database (like get_sorce_pla.pro),
;   this routines looks at the telemetry to identify a subset of
;   plans.
;
; CALLING SEQUENCE:
;   result = FIND_SIM_PLAN(t0, t1, /simA, activity="SolarQuickScan24")
;
; INPUT PARAMETERS:
;   startTime -
;      The lower time range for which data will be returned.
;      May be specified as Julian days, gps microseconds, 
;      VMS time string or SORCE mission day numbers (the default).
;
;   stopTime -
;      The upper time range for which data will be returned
;      May be specified as Julian days, gps microseconds, 
;      VMS time string or SORCE mission day numbers (the default).
;
; OPTIONAL INPUT KEYWORDS:
;   simA, simB - 
;      Specifies which SIM channel to retreive. Zero, one or both 
;      instruments can be requested. If neither are selected, it defaults
;      to both simA and simB.
;   instrument -
;      Specifies the exact instrument to extract the plan: sim_a or sim_b
;   gps - 
;      If set, indicates that input times are to be interpreted as
;      microseconds since Jan 6, 1980 midnight UT.
;   missionDays - 
;      If set, indicates that input times are to be interpreted as
;      days elapsed since launch, i.e. mission days.  If only the startTime
;      argument contains a non-zero value (stopTime is not supplied), then
;      all data for the specified mission day number will be retrieved.
;      (default)
;   julianDays - 
;      If set, indicates that the input times are to be interpreted as
;      julian day numbers.
;   activity_name -
;      If set, returns only the activities whose name matches the provided 
;      name.  If not provided or an empty string, will return all activities
;      for the requested instrument. Wild cards are accepted (in the SQL sense).
;      The % character is used as the wildcard.  An activity of 'SolarQuickScan%'
;      will return all activities starting with SolarQuickScan.
;   noconvert -
;      If set, will NOT convert the returned start and stop times to the same
;      format as the input starttime and stoptime and leave them in the default
;      format returned from the planning database (Sybase date/time format).
;
; RETURNED PARAMETERS:
;   An array of structures whose fields correspond to the returned
;   database column names. The format for the returned timestamps
;   will be the same as the input timestamp.
;
; OPTIONAL OUTPUT PARAMETERS:
;      NONE
;
; EXAMPLES:  
;   ; get all activities for SimA from day 453 to 500
;   plans = get_sorce_plan(453,500,/mission,/sima)  
;
;   ; only return the list of SolarQuickScan24 activities for SimA
;   plans = get_sorce_plan(453,500,/mission,/sima, activity='SolarQuickScan24')  
;
;   ; get all the activities for the TIM instrument between day 453 and 500
;   plans = get_sorce_plan(453,500,/mission, instrument='Tim')
;
;   ; get the list of MUVSolarMiniScan for the TIM instrument between day 453 and 500
;   plans = get_sorce_plan(453,500,/mission, instrument='SolsticeA', activity='MUVSolarMiniScan')
;
; COMMON BLOCKS:
;      NONE
;
; NOTE:
;    Based and expanded from Chris Pankrtaz get_sim_planned_times.pro
;
; REVISION HISTORY:
;   Revision: $Id: find_sim_plan.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;-
;*****************************************************************************************
function find_sim_plan, starttime, stoptime, activity=activity, sima=sima, $
             simb=simb, gps=gps, missionDays=missionDays, $
             julianDays=julianDays, noconvert=noconvert

    ; will extract data for SIM_A by default
    if n_elements(instrument) eq 0 then instrument='SimA'
    if keyword_set(simB) then instrument='SimB'
    if keyword_set(simA) then instrument='SimA'
    if instrument eq 'SimB' then ext_elem='sim_b' else ext_elem='sim_a'

    if keyword_set(julianDays) then begin
        ; user specified time in julian days
        t0 = (jd2sd(startTime))
        t1 = (jd2sd(stopTime))
    endif else if keyword_set(gps) then begin
        ; user specified timetags in gps microseconds
        t0 = (gps2sd(startTime/1d6))
        t1 = (gps2sd(stopTime/1d6))
    endif else begin
        ; user specified time in mission (sorce) days
        missionDays=1
        t0 = startTime
        t1 = stopTime
    endelse
 
    get_sorce_telemetry, data,info, t0, t1, missionDays=missionDays, gps=gps, julianDays=julianDays, $
        externalElement=ext_elem, itemName=['shutter_pos','position']

    ; only consider when shutter is open
    is_open = where((*data.(0)).science.eu eq 0)
    prism_pos = ((*data.(1)).science.eu)[is_open]
    delta_prism = abs(prism_pos[0:-2] - prism_pos[1:-1])

    ; look for SolarQuickScan24 (steps of 38 subpixels)
    k=where(delta_prism eq 38d,count)
    if count eq 0 and strupcase(activity) eq 'SOLARQUICKSCAN24' then return,-1
    
    ; find the end of the SolarQuickScan24
    delta_time = ((*data.(1)).science.timetag)[is_open[k[1:-1]]] - ((*data.(1)).science.timetag)[is_open[k[0:-2]]]
    ; looking for a difference larger than 10 seconds
    q=where(delta_time gt 10d6,count)
    if count eq 0 and strupcase(activity) eq 'SOLARQUICKSCAN24' then return,-1
    s0=(*data.(0))[is_open[k[q]]].science.timetag

    ; find the end of the SolarQuickScan24
    q=where(delta_time gt 10d6,count)
    if count eq 0 and strupcase(activity) eq 'SOLARQUICKSCAN24' then return,-1
    ; add 60 seconds to the end 
    s1=((*data.(1)).science.timetag)[is_open[k[q]]] + 60d6
    s0=((*data.(1)).science.timetag)[is_open[k[q]]] - 1515d6

    out_data = replicate({starttime:0d, stoptime:0d, activityname:'SolarQuickScan24', planningname:instrument, commandversion:''}, n_elements(s0))
    out_data.starttime=gps2sd(s0/1d6)
    out_data.stoptime=gps2sd(s1/1d6)

    return, out_data
    

end
