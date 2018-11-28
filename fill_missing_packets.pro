;+
; NAME:   FILL_MISSING_PACKETS
;
; PURPOSE: 
;   Look at the SIM_HRT solar exposure data and fill the missing
;   science packets with the solar exposure estimated from the 
;   planned activity database. For now, the data for the whole 
;   orbit will be replaced by the estimated times from the plan
;   database.
;
; CALLING SEQUENCE:
;   out_data = FILL_MISSING_PACKETS(inFile, outFile=outFile, planFile=planFile, minGap=minGap)
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
;   Revision: $Id: fill_missing_packets.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;-
;*****************************************************************************
FUNCTION FILL_MISSING_PACKETS, infile, outFile=outFile, minGap=minGap

   ; check for the inFile
   if n_elements(infile) eq 0 then begin
       infile = dialog_pickfile(title='Select the SIM Solar Exposure IDL save file')
       if strlen(infile) eq 0 then return,-1
       infile=infile[0]
   endif
   restore, inFile ; expecting the structure name to be SIM_EXP
   if size(sim_exp,/tname) ne 'STRUCT' then begin
       print,'file '+inFile+' does not conatain expected structure SIM_EXP'
       return,-1
   endif

   defsysv,'!SIM_PLANS_HRT',exists=exists
   if exists eq 0 then begin
       res=def_sim_plans()
       if res ne 1 then begin
           print,'No valid !SIM_PLANS_HRT found'
           return,-1
       endif
   endif

   if n_elements(minGap) eq 0 then minGap = 300.0d

   ; determine which instrument, HRT and pointing details the inFile corresponds to
   instrument=sim_exp[0].instrument
   hrtpos = sim_exp[0].hrtpos
   boxwidth = sim_exp[0].boxwidth

   ; look for every orbit with a sim_exp.gap_sci larger than minGap
   out_data=sim_exp
   gappos = where(out_data.gap_sci ge minGap,gapcount)
   if gapcount eq 0 then begin
       print,'No gaps larger that minGap ('+strtrim(string(minGap),2)+') found'
       return,-1
   endif

   gapstr=strtrim(string(gapcount),2)
   ; loop for every orbit with a gap
   for orb=0,gapcount-1 do begin
       ; get the plans for this orbit
       plan = get_sorce_plan(out_data[gappos[orb]].starttime, out_data[gappos[orb]].endtime, $
               instrument=instrument, /gps)
       if size(plan,/tname) ne 'STRUCT' then begin
           ; no plans found for this orbit 
           ; it is safe to assume nothing was going on during this time
           print, 'orbit '+strtrim(string(orb+1),2)+'/'+gapstr+'  processing 0 plans'
           out_data[gappos[orb]].gap_sci = 0.0d
           out_data[gappos[orb]].gap_sc = 0.0d
           out_data[gappos[orb]].gap_hrt = 0.0d
           out_data[gappos[orb]].gap_power = 0.0d
           out_data[gappos[orb]].abs_gap = 0.0d
           out_data[gappos[orb]].hrt_from_plan = out_data[gappos[orb]].insun_time
           out_data[gappos[orb]].solar_exp = 0.0d
           continue  
       endif 
       print, 'orbit '+strtrim(string(orb+1),2)+'/'+gapstr+$
           '  processing '+strtrim(n_elements(plan),2)+' plans'

       ; loop through every plan in this orbit and update the structure
       solar_exp=0.0D
       solar_exp_stdev=0.0D
       plan.activityName=strupcase(plan.planningName)+'_'+strupcase(plan.activityName)
       for i=0L, n_elements(plan)-1L do begin
           ; loop over every planned activity to determine the state of HRT
           p=where(!SIM_PLANS_HRT.activityName eq plan[i].activityName and $
                   !SIM_PLANS_HRT.activityVersion eq plan[i].commandVersion ,count)
           if count eq 0 then begin
               ; no matching activityName was found
               print, '   No matching activityName/Version ('+ plan[i].activityName+'/'+plan[i].commandVersion+')'
               continue 
           endif

           if hrtpos eq 'IN' then begin
               if boxwidth eq 0.0 then begin
                   solar_exp+=!SIM_PLANS_HRT.inFov[p[0]]
                   solar_exp_stdev+=!SIM_PLANS_HRT.inFovStd[p[0]]
               endif else begin
                   solar_exp+=!SIM_PLANS_HRT.in10[p[0]]
                   solar_exp_stdev+=!SIM_PLANS_HRT.in10Std[p[0]]
               endelse
           endif else if hrtpos eq 'OUT' then begin
               if boxwidth eq 0.0 then begin
                   solar_exp+=!SIM_PLANS_HRT.outFov[p[0]]
                   solar_exp_stdev+=!SIM_PLANS_HRT.outFovStd[p[0]]
               endif else begin
                   solar_exp+=!SIM_PLANS_HRT.out10[p[0]]
                   solar_exp_stdev+=!SIM_PLANS_HRT.out10Std[p[0]]
               endelse
           endif else begin
               if boxwidth eq 0.0 then begin
                   solar_exp+=!SIM_PLANS_HRT.outFov[p[0]]
                   solar_exp_stdev+=!SIM_PLANS_HRT.outFovStd[p[0]]
                   solar_exp+=!SIM_PLANS_HRT.inFov[p[0]]
                   solar_exp_stdev+=!SIM_PLANS_HRT.inFovStd[p[0]]
               endif else begin
                   solar_exp+=!SIM_PLANS_HRT.out10[p[0]]
                   solar_exp_stdev+=!SIM_PLANS_HRT.out10Std[p[0]]
                   solar_exp+=!SIM_PLANS_HRT.in10[p[0]]
                   solar_exp_stdev+=!SIM_PLANS_HRT.in10Std[p[0]]
               endelse
           endelse
       endfor
       out_data[gappos[orb]].gap_sci = 0.0d
       out_data[gappos[orb]].gap_sc = 0.0d
       out_data[gappos[orb]].gap_hrt = 0.0d
       out_data[gappos[orb]].gap_power = 0.0d
       out_data[gappos[orb]].abs_gap = solar_exp_stdev
       out_data[gappos[orb]].hrt_from_plan = out_data[gappos[orb]].insun_time
       out_data[gappos[orb]].solar_exp = solar_exp
   endfor

   return,out_data

END 

