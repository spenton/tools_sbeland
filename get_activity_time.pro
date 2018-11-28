;+
; Author: Stephane Beland
;
; PURPOSE: 
;   The program identifies each activityName/version (read from
;   an input file) that was executed during the specified time 
;   range and calculates the solar exposure for each of these 
;   for four different cases:
;      HRTIN_FOV
;      HRTIN_10arc
;      HRTOUT_FOV
;      HRTOUT_10arc
;
; CALLING SEQUENCE:
;   result = GET_ACTIVITY_TIME(t0, t1, /mission, actfile="sima_clean_rlog.txt")
;
; INPUT PARAMETERS:
;   startTime -
;      The lower time range for which data will be returned.
;      May be specified as Julian days, gps microseconds, 
;      VMS time string (the default) or SORCE mission day numbers.
;
;   stopTime -
;      The upper time range for which data will be returned
;      May be specified as Julian days, gps microseconds, 
;      VMS time string (the default) or SORCE mission day numbers.
;
; OPTIONAL INPUT KEYWORDS:
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
;   actfile -
;      Will read the list of activityName and version numbers from this ASCII
;      file.  If not specified, the program will promt the user for the file.
;   activity -
;      If provided with version, only extract the specific activity.  The name
;      should be from dbp.PlanningSummary.commandLink:  'SimA_SolarQuickScan24'
;   version -
;      If provided with activity, only extract the specific activity/version.
;   maxplan -
;      Specifies the maximum number of plans to extract to calculate the stdev.
;      Defaults to 10.
;
; RETURNED PARAMETERS:
;   An array of structures.
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
;
; REVISION HISTORY:
;   Revision: $Id: get_activity_time.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;-
;
function get_activity_time, starttime, stoptime, actfile=actfile, $
            activity=activity, version=version, hrtstate=hrtstate, maxplan=maxplan, $
            gps=gps, missionDays=missionDays, julianDays=julianDays, verbose=verbose

    if n_elements(actfile) gt 0 then begin
        readcol,actfile,activity,version,ci_date_vms,ci_date_gps,hrtstate,format='(A,A,A,D,A)',delimiter=','
    endif else if n_elements(activity) gt 0 and n_elements(version) gt 0 then begin
        activity = strtrim(activity,2)
        version = strtrim(string(version),2)
    endif else begin
        actfile = dialog_pickfile(title='Select the file with the list of Activities')
        if strlen(actfile) eq 0 then return,-1
        readcol,actfile,activity,version,ci_date_vms,ci_date_gps,hrtstate,format='(A,A,A,D,A)',delimiter=','
    endelse

    ; loop through the whole list of activity/version
    nact=n_elements(activity)

    out_data = replicate({activity:'', version:'',$
                          ci_date_vms:'', ci_date_gps:0.0d, hrtstate:'', $
                          time_in_fov:0.0d,  time_in_fov_std:0.0d, $
                          time_in_10:0.0d,   time_in_10_std:0.0d, $
                          time_out_fov:0.0d, time_out_fov_std:0.0d, $
                          time_out_10:0.0d,  time_out_10_std:0.0d}, nact)
 
    if n_elements(maxplan) eq 0 then maxplan=10 else maxplan = maxplan > 1

    out_data[*].activity = activity[*]
    out_data[*].version = version[*]
    out_data[*].hrtstate = hrtstate[*]
    if n_elements(ci_date_vms) eq nact then out_data[*].ci_date_vms = ci_date_vms[*] 
    if n_elements(ci_date_gps) eq nact then out_data[*].ci_date_gps = ci_date_gps[*]

    for i=0L,nact-1L do begin
        ; look through the Planning summary and find when this specific activity/version was executed
        planningName=strmid(activity[i],0,4)  ; SimA or SimB
        activityName = strmid(activity[i],5)  ; SolarQuickScan24, ...
        commandVersion=strtrim(version[i],2)  ; version string "1.10"

        ; get the start/stop times of scheduled activity/version from the database
        print,'Getting plans for '+planningName+'/'+activityName+', version '+commandVersion
        plan = get_sorce_plan(startTime, stopTime, missionDays=missionDays, gps=gps, $
               julianDays=julianDays, instrument=planningName, activity=activityName, $
               version=commandVersion)

        if size(plan,/tname) ne 'STRUCT' then continue

        print,'   found '+strtrim(string(n_elements(plan)),2)+' matching plans'
        simA='SimA'
        simB='SimB'
        if strupcase(planningName) eq 'SIMA' then simB=''
        if strupcase(planningName) eq 'SIMB' then simA=''
        ; limit the number of plans to extract to maxplan and disperse the ones we use
        if n_elements(plan) gt maxplan then begin
           stp = floor(n_elements(plan) / maxplan)
           diff = floor(n_elements(plan) - (maxplan-1.0)*stp -1.0)
           p=lindgen(maxplan)*stp + diff
           plan=plan[p]
        endif
        t0=vms2gps(jdbc2vms(plan.starttime)) * 1d6
        t1=vms2gps(jdbc2vms(plan.stoptime)) * 1d6
        time_in_fov=[]
        time_in_10=[]
        time_out_fov=[]
        time_out_10=[]
        for k=0L,n_elements(plan)-1L do begin
            print,'   Extracting the SOLAR_EXP ',k+1
            res_in_10 = CALC_SIM_EXP(t0[k], t1[k], /gps, simA=simA, simB=simB, boxwidth=20.0, $
                                     hrtpos='IN', /orbitForce, /silent)
            res_in_FOV = CALC_SIM_EXP(t0[k], t1[k], /gps, simA=simA, simB=simB, boxwidth=0.0, $
                                     hrtpos='IN', /orbitForce, /silent)
            res_out_10 = CALC_SIM_EXP(t0[k], t1[k], /gps, simA=simA, simB=simB, boxwidth=20.0, $
                                     hrtpos='OUT', /orbitForce, /silent)
            res_out_FOV = CALC_SIM_EXP(t0[k], t1[k], /gps, simA=simA, simB=simB, boxwidth=0.0, $
                                     hrtpos='OUT', /orbitForce, /silent)
            time_in_10=[time_in_10,res_in_10[0].solar_exp]
            time_in_fov=[time_in_fov,res_in_fov[0].solar_exp]
            time_out_10=[time_out_10,res_out_10[0].solar_exp]
            time_out_fov=[time_out_fov,res_out_fov[0].solar_exp]
            if keyword_set(verbose) then $
                print,gps2sd(t0[k]/1d6),'  '+hrtstate,res_in_fov[0].solar_exp,res_in_10[0].solar_exp,res_out_fov[0].solar_exp,res_out_10[0].solar_exp
        endfor
        out_data[i].time_in_10  = median(time_in_10)
        if n_elements(time_in_10) gt 1 then out_data[i].time_in_10_std = stdev(time_in_10)
        out_data[i].time_in_fov = median(time_in_fov)
        if n_elements(time_in_fov) gt 1 then out_data[i].time_in_fov_std = stdev(time_in_fov)
        out_data[i].time_out_10  = median(time_out_10)
        if n_elements(time_out_10) gt 1 then out_data[i].time_out_10_std = stdev(time_out_10)
        out_data[i].time_out_fov = median(time_out_fov)
        if n_elements(time_out_fov) gt 1 then out_data[i].time_out_fov_std = stdev(time_out_fov)
    endfor

    print,'             '
    return, out_data

end
