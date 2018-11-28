;+
; NAME:   GET_EXP_FROM_PLAN
;
; PURPOSE: 
;   Look at the SIM_HRT solar exposure data and fill the missing
;   science packets with the solar exposure estimated from the 
;   planned activity database. For now, the data for the whole 
;   orbit will be replaced by the estimated times from the plan
;   database.
;
; CALLING SEQUENCE:
;   out_data = GET_EXP_FROM_PLAN(inFile, outFile=outFile, planFile=planFile, minGap=minGap)
;
; INPUT PARAMETERS:
;   inFile -
;      The file name of the IDL save file with the SIM_EXP structure.
;
; OPTIONAL INPUT KEYWORDS:
;   outFile -
;      The name of the new IDL save file with the updated data.
;
;   planFile - 
;      The name of the ASCII file with the list of activities and
;      their solar exposure times. If the file is not specified,
;      the program will use the data from the system variable !SIM_PLAN_HRT
;      if defined, otherwise, will prompt the user for the right file.
;
;   minGap -
;      The minimum amount of gap time (in seconds) to consider.
;      Defaults to 300.0 seconds (any gap_sci smaller than this will be ignored).
;
; OUTPUT PARAMETERS:
;   out_data -
;      Will be a copy of the inpu SIM_EXP with the solar exposures updated for
;      orbits with missing science packets.
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
;   2012.03.05  SBeland
;   Revision: $Id: get_exp_from_plan.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;-
;*****************************************************************************
FUNCTION GET_EXP_FROM_PLAN, starttime, stoptime, instrument=instrument, hrtpos=hrtpos, boxwidth=boxwidth, $
         gps=gps, missionDays=missionDays, julianDays=julianDays, verbose=verbose

    defsysv,'!SIM_PLANS_HRT',exists=exists
    if exists eq 0 then begin
        ; res=def_sim_plans()
        res=insert_sim_plans()
        if size(res,/tname) ne 'STRUCT' then begin
            print,'No !SIM_PLANS_HRT was found'
            return,-1
        endif
       defsysv,'!SIM_PLANS_HRT',temporary(res) 
    endif

   if keyword_set(missionDays) then begin
      ; user specified time in mission (sorce) days
      t0 = ulong64(sd2gps(startTime)*1.d6)
      t1 = ulong64(sd2gps(stopTime)*1.d6)
   endif else if keyword_set(julianDays) then begin
      ; user specified time in julian days
      t0 = ulong64(jd2gps(startTime)*1.d6)
      t1 = ulong64(jd2gps(stopTime)*1.d6)
   endif else begin
      ; user specified timetags in gps microseconds
      t0 = ulong64(startTime)
      t1 = ulong64(stopTime)
   endelse

   ; defaults to SimA
   if n_elements(instrument) eq 0 then instrument="SIMA" else instrument=strupcase(strtrim(instrument,2))
   ; hrt defaults to NA
   if n_elements(hrtpos) eq 0 then hrtpos='NA'
   ; defautls the boxwidth to FOV
   if n_elements(boxwidth) eq 0 then boxwidth=0.0

   out_data = {time_from_plan:(t1-t0)/1d6, solar_exp:0.0d, solar_exp_stdev:0.0d}
 
   ; get the plans for this time range
   plan = get_sorce_plan(t0, t1, instrument=instrument, /gps)
   if size(plan,/tname) ne 'STRUCT' then begin
       ; no plans found for this orbit 
       ; it is safe to assume nothing was going on during this time
       return, out_data
   endif 

   sim_plans_hrt_gps = double(!SIM_PLANS_HRT.ci_date_gps)

   ; loop through every plan in this orbit and update the structure
   ;plan.activityName=strupcase(plan.planningName)+'_'+strupcase(plan.activityName)
   plan.activityName=plan.planningName+'_'+plan.activityName
   plan_duration = (plan.stoptime - plan.starttime)

   for i=0L, n_elements(plan)-1L do begin
       ; loop over every planned activity to determine the state of HRT
       p=where(!SIM_PLANS_HRT.commandLink eq plan[i].activityName and $
               !SIM_PLANS_HRT.version eq plan[i].commandVersion ,count)
       if count eq 0 then begin
           ; check if only the version is missing in which case use the plan from
           ; version before orbit start time
           if plan[i].commandVersion ne '' then begin
               ; no matching activityName was found
               if keyword_set(verbose) then print, '   No matching activityName/Version ('+ plan[i].activityName+'/'+$
                   plan[i].commandVersion+')'
               continue 
           endif
           p=where(!SIM_PLANS_HRT.commandLink eq plan[i].activityName and $
                   sim_plans_hrt_gps lt plan[i].starttime ,count)
           if count gt 0 then p=p[-1] else begin
               ; no matching activityName was found
               if keyword_set(verbose) then print, '   No matching activityName/Version ('+ plan[i].activityName+'/'+$
                   plan[i].commandVersion+')'
               continue 
           endelse

       endif

       ; if only a partial activity is run during the requested time
       ; we'll ratio the solar_exp for this activity by the fraction of
       ; the activity requested (not perfect but gets us pretty close)
       if t0 gt plan[i].starttime or t1 lt plan[i].stoptime then begin
           factor = (min([plan[i].stoptime,t1/1d6]) - max([plan[i].starttime,t0/1d6])) / plan_duration[i]
       endif else factor=1.0d

       if hrtpos eq 'IN' then begin
           if boxwidth eq 0.0 then begin
               out_data.solar_exp+=(!SIM_PLANS_HRT.inFov[p[0]]*factor)
               out_data.solar_exp_stdev+=!SIM_PLANS_HRT.inFovStd[p[0]]
           endif else begin
               out_data.solar_exp+=(!SIM_PLANS_HRT.in10[p[0]]*factor)
               out_data.solar_exp_stdev+=!SIM_PLANS_HRT.in10Std[p[0]]
           endelse
       endif else if hrtpos eq 'OUT' then begin
           if boxwidth eq 0.0 then begin
               out_data.solar_exp+=(!SIM_PLANS_HRT.outFov[p[0]]*factor)
               out_data.solar_exp_stdev+=!SIM_PLANS_HRT.outFovStd[p[0]]
           endif else begin
               out_data.solar_exp+=(!SIM_PLANS_HRT.out10[p[0]]*factor)
               out_data.solar_exp_stdev+=!SIM_PLANS_HRT.out10Std[p[0]]
           endelse
       endif else begin
           if boxwidth eq 0.0 then begin
               out_data.solar_exp+=(!SIM_PLANS_HRT.outFov[p[0]]*factor)
               out_data.solar_exp_stdev+=!SIM_PLANS_HRT.outFovStd[p[0]]
               out_data.solar_exp+=(!SIM_PLANS_HRT.inFov[p[0]]*factor)
               out_data.solar_exp_stdev+=!SIM_PLANS_HRT.inFovStd[p[0]]
           endif else begin
               out_data.solar_exp+=(!SIM_PLANS_HRT.out10[p[0]]*factor)
               out_data.solar_exp_stdev+=!SIM_PLANS_HRT.out10Std[p[0]]
               out_data.solar_exp+=(!SIM_PLANS_HRT.in10[p[0]]*factor)
               out_data.solar_exp_stdev+=!SIM_PLANS_HRT.in10Std[p[0]]
           endelse
       endelse
   endfor
   return,out_data

END 

