;+
; NAME:   CALC_SIM_EXP
;
; PURPOSE: 
;   Calculates the total time where either SIM instrument were exposed
;   to solar radiation, per orbit within the provided time range.
;
; CALLING SEQUENCE:
;   out_data = CALC_SIM_EXP(startTime, stopTime, simA=simA, simB=simB, 
;                           center=center, hrtpos=hrtpos)
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
;   boxwidth -
;      If specified, will only use the data when the instrument is pointing 
;      within a box of "width" arcminutes centered on the instrument FOV.
;      If not specified (or set to 0.0), the program will use the data as 
;      soon as the pointing falls inside the SIM FOV (2.25 x 2.5°).
;   hrtpos -
;      Specifies the state of the Hard Radiation Trap to consider (IN or OUT).
;      If not specified, will use data when the HRT is in either position.
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
;   orbitForce -
;      If this flag is set and the program can't find a complete orbit
;      within the specified timerange, we will assume one orbit starting 
;      at startTime and ending at stopTime. This is usefull when the 
;      specified time range is a fraction of an orbit.
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
;    If the time range starts when the satellite is in the sun, the partial orbit is ignored
;    and the following full orbit will be the new startTime.  The same is true for the 
;    stopTime falling when sun is present. The last partial orbit will be ignored.
;    It is possible to not have any data returned if the startTime and stopTime don't
;    span a full orbit.
;    An orbit is defined by the telemetry mu$sunpresence (tlmId=8025). For SIM, the orbit
;    spans a transition from IN_SUN to OUT_OF_SUN to next similar transition.
;
; REVISION HISTORY:
;   2011.11.03  SBeland
;   Revision: $Id: calc_sim_exp.pro,v 1.2 2018/11/26 16:16:02 spenton Exp $
;-
;*****************************************************************************
FUNCTION FIX_CHANGES, hrt, shutter, t0, t1
    
    pos = where(hrt.timetag ge t0 and hrt.timetag le t1,count)
    if count le 1 then return,{timetag:hrt.timetag, dn:hrt.dn}

    ; if no shutter information simply return
    if size(shutter,/tname) ne 'STRUCT' then return,{timetag:hrt.timetag, dn:hrt.dn}

    hrt_dn=hrt[pos].dn
    hrt_trap_time=hrt[pos].timetag
    delta_hrt = hrt_dn[0:-2] - hrt_dn[1:-1]
    hrt_change = where(delta_hrt ne 0, n_hrt_change)

    if n_hrt_change eq 0 then return, {timetag:hrt_trap_time, dn:hrt_dn}

    hrt_change_times = dblarr(2, n_hrt_change)
    hrt_change_times[0,*] = hrt_trap_time[hrt_change]   ;before change
    hrt_change_times[1,*] = hrt_trap_time[hrt_change+1] ;after change
    
    delta_shutter = shutter[0:-2].dn - shutter[1:-1].dn
    shutter_change = where(delta_shutter ne 0, n_shutter_change)
    shutter_change_time = shutter[shutter_change].timetag

    found_it=0
    for i=0L,n_hrt_change-1L do begin
        pos_sc = where(shutter_change_time gt hrt_change_times[0,i]  $
                  and shutter_change_time lt hrt_change_times[1,i], nsctime)
        if nsctime eq 0 then continue ; no shutter change during this HRT change
        ; get the time of the last shutter change
        sctime = shutter_change_time[pos_sc[-1]]
        found_it=1
        new_hrt_sample = hrt_dn[hrt_change[i]+1]
        hrt_dn=[temporary(hrt_dn),new_hrt_sample]
        hrt_trap_time=[temporary(hrt_trap_time),sctime]
    endfor

    if found_it then begin
        pos = sort(hrt_trap_time)
        hrt_dn=hrt_dn[pos]
        hrt_trap_time=hrt_trap_time[pos]
    endif

    return, {timetag:hrt_trap_time, dn:hrt_dn}

END
;*****************************************************************************
FUNCTION CALC_SIM_EXP, startTime, stopTime, simA=simA, simB=simB,  orbits=orbits, $
         minSciGap=minSciGap, boxwidth=boxwidth, hrtpos=hrtpos, silent=silent, $
         gps=gps, missionDays=missionDays, julianDays=julianDays, $
         savefile=savefile, missing_plans=missing_plans, orbitForce=orbitForce, indata=indata

   ; will extract data for SIM_A by default
   instrument = "SimA"
   if (keyword_set(simB) and NOT keyword_set(simA)) then begin
      instrument = "SimB"
   endif

   if n_elements(boxwidth) eq 0 then boxwidth=0.0 else boxwidth=abs(boxwidth)
   if n_elements(hrtpos) eq 0 then begin
       hrtpos='NA'
   endif else begin
       hrtpos = strupcase(strtrim(hrtpos,2))
       if strpos(hrtpos,'IN') lt 0 and strpos(hrtpos,'OUT') lt 0 then begin
           print,'Error: hrtpos should only have values of "IN" or "OUT"'
        return,-1
       endif
   endelse

   if size(indata,/tname) eq 'STRUCT' then begin
       out_data=indata
       HRTPOS = out_data[0].HRTPOS
       boxwidth= out_data[0].boxwidth
       instrument = out_data[0].instrument
       GOTO,CLEANUP
   endif

   SIMA_FSS0_CENTER = -3.8420044 ; fssunangle0 value when centered in FOV (arcmin)
   SIMA_FSS1_CENTER = -1.0706056 ; fssunangle1 value when centered in FOV (arcmin)
   ;
   ;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ;
   ;  NEED TO VERIFY THE VALUES OF FSSUNAGLE0 AND FSSUNANGLE1 for SIM_B when centered
   ;
   ;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ;
   SIMB_FSS0_CENTER = -3.8420044 ; fssunangle0 value when centered in FOV (arcmin)
   SIMB_FSS1_CENTER = -1.0706056 ; fssunangle1 value when centered in FOV (arcmin)
   FSS0_FOV = 2.25 * 60.0    ; FOV in arcmin
   FSS1_FOV = 2.50* 60.0     ; FOV in arcmin

   if instrument eq "SimA" then begin
       TMID = {fssunangle0:    7871L, $
               fssunangle1:    7872L, $ 
               fssunanglevld:  7749L, $
               sunpresence:    8025L, $ 
               acsmd:          7926L, $
               shutter:       21019L, $
               rad_trap_in:   20245L, $
               rad_trap_out:  20035L, $ 
               uv_temp:       20592L, $
               uv_tlm:       100046L, $
               pri_power:     20242L  $
       }
       tlmid = {fssunangle0:   7871L, $
               fssunangle1:    7872L, $ 
               fssunanglevld:  7749L, $
               sunpresence:    8025L, $ 
               acsmd:          7926L, $
               shutter:      102121L, $
               rad_trap_in:  100331L, $
               rad_trap_out: 100038L, $ 
               uv_temp:      101069L, $
               uv_tlm:       100046L, $
               pri_power:    100328L  $
       }
       fss0_center = SIMA_FSS0_CENTER
       fss1_center = SIMA_FSS1_CENTER
       uv_min_signal = -10.0
   endif else begin
       TMID = {fssunangle0:    7871L, $
               fssunangle1:    7872L, $ 
               fssunanglevld:  7749L, $
               sunpresence:    8025L, $ 
               acsmd:          7926L, $
               shutter:       20177L, $
               rad_trap_in:   20374L, $
               rad_trap_out:  20815L, $ 
               uv_temp:       20816L, $
               uv_tlm:       100057L, $
               pri_power:     20243L  $
       }
       tlmid = {fssunangle0:   7871L, $
               fssunangle1:    7872L, $ 
               fssunanglevld:  7749L, $
               sunpresence:    8025L, $ 
               acsmd:          7926L, $
               shutter:      100197L, $
               rad_trap_in:  100582L, $
               rad_trap_out: 101579L, $ 
               uv_temp:      101580L, $
               uv_tlm:       100046L, $
               pri_power:    100329L  $
       }
       fss0_center = SIMB_FSS0_CENTER
       fss1_center = SIMB_FSS1_CENTER
       uv_min_signal = -10.0
   endelse

   ; prepare to filter the data by distance from center
   if boxwidth eq 0.0 then begin
       min_fss0 = fss0_center - FSS0_FOV/2.0
       max_fss0 = min_fss0 + FSS0_FOV
       min_fss1 = fss1_center - FSS1_FOV/2.0
       max_fss1 = min_fss1 + FSS1_FOV
   endif else begin
       min_fss0 = fss0_center - boxwidth/2.0
       max_fss0 = min_fss0 + boxwidth
       min_fss1 = fss1_center - boxwidth/2.0
       max_fss1 = min_fss1 + boxwidth
   endelse

   ; initialize the output structure
   tmpstr = {instrument:instrument, boxwidth:boxwidth, hrtpos:hrtpos, $
             starttime:0.0d, endtime:0.0d, insun_time:0.0d, $
             gap_sci:0.0d, gap_sc:0.0d, gap_sc_open:0.0d, $
             gap_hrt:0.0d, same_hrt:0.0d, $
             gap_power:0.0d, power_off:0.0d, $
             abs_gap:0.0d, hrt_from_plan:0.0d, $
             time_from_plan:0.0d, time_std_from_plan:0.0d, $
             safehold:0.0d, solar_exp:0.0d } 

   if n_elements(orbits) eq 0 then $
       orbits = get_sim_orbits(startTime, stopTime, simA=simA, simB=simB, silent=silent, $
                     gps=gps, missionDays=missionDays, julianDays=julianDays)

   if size(orbits,/tname) ne 'STRUCT' then begin
      if not keyword_set(orbitForce) then begin
          ; there was an error in determining the orbits
          print,'Error: no valid complete orbit for the specified time range'
          return,-1
      endif
      ; force the provided times as one orbit
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
      orbits={instrument:instrument, starttime:t0, endtime:t1, insun_time:(t1-t0)/1d6}
   endif 

   norbits=n_elements(orbits)
   missing_plans=PTRARR(norbits)
   out_data = replicate(tmpstr,norbits)
   insun_t0 = orbits.starttime
   insun_t1 = orbits.endtime
   out_data.starttime = insun_t0
   out_data.endtime = insun_t1
   out_data.insun_time = (insun_t1-insun_t0)/1.0d6

   ; get the conversion coeffs for the fssunangles to transform from dn to eu
   q1 = "select c0,c1,c2,c3,c4,c5,c6,c7 from SORCE_CT_SC.dbo.TelemetryAnalogConversions where tlmId="
   query_database, q1+strtrim(TMID.fssunangle0,2), coeffs, nrows, user='sorce', password='sorcedb', server='sorce-db', database='SORCE'
   fss0_coeffs = [coeffs.c0, coeffs.c1, coeffs.c2, coeffs.c3, coeffs.c4, coeffs.c5, coeffs.c6, coeffs.c7]

   query_database, q1+strtrim(TMID.fssunangle1,2), coeffs, nrows
   fss1_coeffs = [coeffs.c0, coeffs.c1, coeffs.c2, coeffs.c3, coeffs.c4, coeffs.c5, coeffs.c6, coeffs.c7]

   ; select the correct database for the shutter depending on the instrument
   if instrument eq 'SimB' then begin
       shutter_db='SORCE_L1S.dbo.SimBScienceSamples'
   endif else begin
       shutter_db='SORCE_L1S.dbo.SimAScienceSamples'
   endelse

   if not keyword_set(silent) then print,' looping over the orbits ...'

   ; process only a subset of orbits at a time to limit the number and size of queries
   orbs2proc = 30.0
   ngrp = ceil(float(norbits)/orbs2proc)

   for orbgrp=0L,ngrp-1L do begin
       
       orb0 = orbgrp*orbs2proc
       orb1 = ((orbgrp+1)*orbs2proc-1L) < (norbits-1L)
       ; stretch the start and end of sunpresence by 10 minutes to cover delay in telemetry
       orb_start = insun_t0[orb0] - 600.0d6
       orb_end = insun_t1[orb1]  + 600.0d6
       if not keyword_set(silent) then print,' querying from '+gps2la(orb_start/1d6)+$
                                             ' to '+gps2la(orb_end/1d6)
       ; get_sorce_telemetry, data, info, insun_t0[orb0], insun_t1[orb1], tlmId=TMID

       ; extract the fssunangle information
       if not keyword_set(silent) then print,format='($,"  fssunangle0 ...")'
       q1 = "select SCT_VTCW 'timetag', Value 'dn', Value 'eu' from SORCE_L1A.dbo.TManalog where TMID="
       q2 = " AND SCT_VTCW BETWEEN "+strtrim(ulong64(orb_start),2)+" AND " + $
            strtrim(ulong64(orb_end),2)+" ORDER by SCT_VTCW"
       q_str = q1+strtrim(TMID.fssunangle0,2)+q2
       query_database, q_str, fss0, nrows
       if size(fss0,/tname) ne 'STRUCT' then begin
           if not keyword_set(silent) then print,format='($,"  no fssunangle0 found -> skipping")'
           ; continue
           ; instead of skipping assume it's pointing dead on for now
           fss0={timetag:orb_start, dn:0.0d, eu:0.0d}
       endif
       fss0.eu = poly(fss0.dn, fss0_coeffs) * 60.0d  ; eu is now in arcmin

       if not keyword_set(silent) then print,format='($,"  fssunangle1 ...")'
       q_str = q1+strtrim(TMID.fssunangle1,2)+q2
       query_database, q_str, fss1, nrows
       if size(fss1,/tname) ne 'STRUCT' then begin
           if not keyword_set(silent) then print,format='($,"  no fssunangle1 found -> skipping")'
           ; continue
           fss1={timetag:orb_start, dn:0.0d, eu:0.0d}
       endif
       fss1.eu = poly(fss1.dn, fss1_coeffs) * 60.0d  ; eu is now in arcmin

       if not keyword_set(silent) then print,format='($,"  fssunanglevld ...")'
       q1 = "select SCT_VTCW 'timetag', Value 'dn' from SORCE_L1.dbo.TMdiscrete where TMID="
       q_str = q1+strtrim(TMID.fssunanglevld,2)+q2
       query_database, q_str, fssv, nrows
       if size(fssv,/tname) ne 'STRUCT' then begin
           if not keyword_set(silent) then print,format='($,"  no fssunanglevld found -> skipping")'
           ; continue
           fssv={timetag:orb_start, dn:0.0d}
       endif

       ; extract the radiation trap status
       if not keyword_set(silent) then print,format='($,"  rad_trap_in ...")'
       q1 = "select SCT_VTCW 'timetag', Value 'dn' from SORCE_L1.dbo.TMdiscrete where TMID="
       q_str = q1+strtrim(TMID.rad_trap_in,2)+q2
       query_database, q_str, rad_trap_in, nrows

       if not keyword_set(silent) then print,format='($,"  rad_trap_out ...")'
       q_str = q1+strtrim(TMID.rad_trap_out,2)+q2
       query_database, q_str, rad_trap_out, nrows

       ; extract the shutter position information
       if not keyword_set(silent) then print,format='($,"  shutter_pos ...")'
       q_str = "select sampleVtcw 'timetag',shutterPosition 'dn' from "+shutter_db+" where " + $
            "sampleVtcw BETWEEN "+strtrim(ulong64(orb_start),2)+" AND " + $
            strtrim(ulong64(orb_end),2)+" ORDER by sampleVtcw"
       query_database, q_str, shutter, nrows

       ; extract the state of the primary power
       if not keyword_set(silent) then print,format='("  pri_power ...")'
       q1 = "select SCT_VTCW 'timetag', Value 'dn' from SORCE_L1.dbo.TMdiscrete where TMID="
       q2 = " AND SCT_VTCW BETWEEN "+strtrim(ulong64(orb_start),2)+" AND " + $
            strtrim(ulong64(orb_end),2)+" ORDER by SCT_VTCW"
       q_str = q1+strtrim(string(TMID.pri_power),2)+q2
       query_database, q_str, pri_power, nrows

       ; process each orbit at a time
       for orb=orb0,orb1 do begin
           if (not keyword_set(silent) and (round(orb+1) MOD 5) eq 0) or orb eq orb1 then $
               print,format='($, TL1,"   orbit '+strtrim(string(long(orb)+1),2)+string(13b)+'")'
           t0=insun_t0[orb]
           t1=insun_t1[orb]
           ; use a 0.1 second sampling bin
           ref_time = dindgen((t1-t0)/1d5+1.0)*1d5+t0
           bintime = (ref_time[1] - ref_time[0])/1d6
           ; keep everything at first and remove time ranges as we go
           gap_sci = bytarr(n_elements(ref_time))
           gap_power = gap_sci * 0b
           gap_hrt = gap_sci * 0b
           same_hrt = gap_sci * 0b
           hrt_from_plan = gap_sci * 0b
           gap_sc = gap_sci * 0b
           gap_sc_open = gap_sci * 0b
           gap_abs = gap_sci * 0b
           fss_keep = gap_sci * 0b

           ; make sure the power is on during at least part of the orbit
           ; the safehold times will be detected here and power for those times will be set to OFF
           power_count=0
           if size(pri_power,/tname) eq 'STRUCT' then $
               power_pos = where(pri_power.timetag ge t0 and pri_power.timetag le t1,power_count)
           if size(pri_power,/tname) ne 'STRUCT' or power_count le 2 then begin
               ; if "orbit" is really small, we may not have any power signal
               ; look for shutter info in this case
               if size(shutter,/tname) eq 'STRUCT' and (t1-t0) lt 300d6 then begin
                   tick_pos = where(shutter.timetag ge t0 and shutter.timetag le t1,tick_count)
                   if tick_count gt 0 then power=replicate(1b,n_elements(ref_time))
               endif else begin
                   ; no power status for this orbit -> safehold ???
                   ; extract the state of the AC Task: ACS Mode
                   safe_time=get_safehold_times(t0,t1, tlmid=TMID.acsmd,ref_time=ref_time)
                   if safe_time[0] ne -1.0  then begin
                       gap_power=temporary(safe_time)
                       out_data[orb].safehold = TOTAL(gap_power) * bintime
                       ; assume power is OFF when in safehold and ON otherwise
                       power = ABS(gap_power-1)
                   endif else begin
                       ; no acsmd info was found to tell if we're in safehold or not
                       ; assume we're not and keep going with power on
                       gap_power[*]=1b
                       power=gap_power
                   endelse
               endelse
           endif else begin
               ; the power telemetry seems normal - look for gaps 
               ; dethin may hide the gaps showing either on or off - depending on joining states
               power = byte(dethin_data(pri_power[power_pos].dn,pri_power[power_pos].timetag, ref_time))
               delta_t = 362.0d6 
               power_time=[t0, pri_power[power_pos].timetag, t1]
               gappos = WHERE(abs(power_time[0:-2] - power_time[1:-1]) GE delta_t,gapcount)
               if gapcount gt 0 then begin
                   ; look for a safehold
                   safe_time=get_safehold_times(t0,t1, tlmid=TMID.acsmd,ref_time=ref_time)
                   if safe_time[0] ne -1.0 then begin
                       gap_power=temporary(safe_time)
                       out_data[orb].safehold = TOTAL(gap_power) * bintime
                       ; mark power as OFF during safehold
                       power = ABS(gap_power-1)
                   endif 
                   for i=0,gapcount-1 do begin
                       ; loop over the gaps and only count times when power tlm is missing
                       p=where(ref_time gt power_time[gappos[i]] and $
                               ref_time lt power_time[gappos[i]+1],count)
                       if count gt 0 then begin
                           ; the gap_power is either from safehold or missing packet
                           gap_power[p] = 1b
                       endif
                   endfor
               endif
           endelse

           ; mark times when power is turned off
           p=where(power eq 0,count)
           if count eq n_elements(power) then begin
               ; the power is off for the whole orbit -> skip
               out_data[orb].power_off = out_data[orb].insun_time
               continue
           endif else if count gt 0 then begin
               out_data[orb].power_off += (count * bintime)
           endif

           ; look at the state of the shutter
           tick_count=0
           if size(shutter,/tname) eq  'STRUCT' then $
               tick_pos = where(shutter.timetag ge t0 and shutter.timetag le t1,tick_count)
           if size(shutter,/tname) ne  'STRUCT' or tick_count le 1 then begin
               ; no science telemetry in this time range (expecting a STRUCTURE)
               ; the missing sci tlm could be from a safehold or missing packets
               safed_time=get_safehold_times(t0,t1, tlmid=TMID.acsmd)
               if safed_time gt 0.0 then begin
                   ; since we have a safehold situation and NO science telemetry
                   ; we can't tell which part of the plan was executed
                   ; simply skip this whole orbit  and mark it as a gap_sci
                   out_data[orb].safehold = safed_time > 0.0
                   out_data[orb].gap_sci = out_data[orb].insun_time
                   ; subtract from the gap, the time when the instrument is off
                   out_data[orb].abs_gap = out_data[orb].insun_time - total(abs(power-1))*bintime
                   continue
               endif

               ; not in safehold/contingency so we must be missing all the science packets

               ; get the plans and times for this orbit
               result = get_hrt_plan(t0, t1, instrument=instrument, hrtpos=hrtpos, boxwidth=boxwidth, /gettimes)
               if size(result,/tname) eq 'STRUCT' then begin
                   out_data[orb].solar_exp = result.solar_exp
                   ; out_data[orb].abs_gap = result.solar_exp_stdev 
                   out_data[orb].time_from_plan = result.solar_exp
                   out_data[orb].time_std_from_plan = result.solar_exp_stdev
                   out_data[orb].hrt_from_plan = result.solar_exp
               endif else begin
                   ; no planned activity for this time range - assume shutter closed
                   out_data[orb].solar_exp = 0.0
               endelse
               continue
           endif else begin
               ; looks like a normal coverage
               shutter_pos = dethin_data(shutter[tick_pos].dn, shutter[tick_pos].timetag, ref_time)
               ; look for gaps in the sci telemetry (dethin will hide them)
               delta_t = 2.0d6   ; look for gaps longer than 2 seconds 
               delta_t1=180.0d6
               gappos = WHERE(abs(shutter[tick_pos[0:-2]].timetag- shutter[tick_pos[1:-1]].timetag) GT delta_t,gapcount)
               if gapcount gt 0 then begin
                   ; we have gaps 
                   ; for gaps smaller than 3 minutes look at the plan runnign during this time
                   gappos = WHERE(abs(shutter[tick_pos[0:-2]].timetag- shutter[tick_pos[1:-1]].timetag) GT delta_t AND $ 
                       abs(shutter[tick_pos[0:-2]].timetag - shutter[tick_pos[1:-1]].timetag) LE delta_t1,gapcount)
                   for i=0,gapcount-1 do begin
                       p=where(ref_time gt shutter[tick_pos[gappos[i]]].timetag and $
                               ref_time lt shutter[tick_pos[gappos[i]+1]].timetag and $
                               power eq 1,count)
                       result = get_exp_from_plan(shutter[tick_pos[gappos[i]]].timetag, shutter[tick_pos[gappos[i]+1]].timetag, /gps, $
                           instrument=instrument, hrtpos=hrtpos,boxwidth=boxwidth)
                       if result.solar_exp eq 0.0 then begin
                           ; shutter was closed during this time - fill the gap
                           shutter_pos[p]=1b
                           gap_sci[p]=0b
                       endif else begin
                           ; shutter is open during some of this time
                           ; since we can't tell exactly what fraction is missing, mark as gap and don't count in solar_exp
                           gap_sci[p]=1b
                           shutter_pos[p]=1b
                       endelse
                   endfor

                   ;for gaps larger than 3  minutes - mark as real gaps
                   gappos = WHERE(abs(shutter[tick_pos[0:-2]].timetag-shutter[tick_pos[1:-1]].timetag) GT delta_t1, gapcount)
                   for i=0,gapcount-1 do begin
                       ; loop over the gaps and only count times when shutter is missing AND power is ON
                       p=where(ref_time gt shutter[tick_pos[gappos[i]]].timetag and $
                               ref_time lt shutter[tick_pos[gappos[i]+1]].timetag and $
                               power eq 1,count)
                       if count gt 0 then begin
                           ; since we don't know the state of the shutter in the gaps
                           ; we won't count them in the solar_exp (assume shutter closed)
                           ; and the gap_sci will be accounted for in abs_gap below
                           ; only when the power is ON (take care of partial safehold)
                           gap_sci[p]=1b
                           shutter_pos[p]=1b
                       endif
                   endfor
               endif
           endelse

           ; a shutter value of 0 indicates an open shutter (1 is close)
           p=where(shutter_pos eq 0,count,complement=comp)
           if count eq 0 then begin
               ; shutter is always closed for this orbit -> skip
               ; out_data[orb].gap_sci = (TOTAL(gap_sci) * bintime)
               ; out_data[orb].abs_gap = (TOTAL(gap_sci) * bintime)
               ; continue
           endif 

           hrtin_count = 0L
           if size(rad_trap_in,/tname) eq 'STRUCT' then $
               hrtin_pos = where(rad_trap_in.timetag ge t0 and rad_trap_in.timetag le t1, hrtin_count)
           hrtout_count = 0L
           if size(rad_trap_out,/tname) eq 'STRUCT' then $
               hrtout_pos = where(rad_trap_out.timetag ge t0 and rad_trap_out.timetag le t1, hrtout_count)
           if hrtpos ne 'NA' then begin
              if size(rad_trap_in,/tname) ne 'STRUCT' or size(rad_trap_out,/tname) ne 'STRUCT' or $
                  hrtin_count le 2 or hrtout_count le 2 then begin
                  ; we're missing the data for the state of the HRT - get it from the plan
                  print,'Warning: missing some radiation trap information - looking at plan'
                  hrt_state = get_hrt_plan(t0,t1,instrument=instrument)
                  if size(hrt_state,/tname) ne 'STRUCT' then begin
                      ; no planned activities were found for the timerange
                      ; count only the times when the shutter is open
                      out_data[orb].gap_hrt= TOTAL(ABS(shutter_pos-1)) * bintime
                      out_data[orb].abs_gap = TOTAL(ABS(shutter_pos-1)) * bintime
                      continue
                  endif
                  hrtin = hrt_state.hrtin
                  hrtout = hrt_state.hrtout
                  hrtin_dn  = dethin_data(hrtin.dn, hrtin.timetag, ref_time)
                  hrtout_dn = dethin_data(hrtout.dn, hrtout.timetag, ref_time)
                  ; only account for HRT gaps when the power is on and shutter is open
                  hrt_from_plan[*] = power AND ABS(shutter_pos-1)
                  gap_hrt=hrt_from_plan 
              endif else begin
                  ; look for holes in the HK telemetry
                  new_hrtin  = fix_changes(rad_trap_in[hrtin_pos], shutter, t0, t1)
                  hrtin_dn  = dethin_data(new_hrtin.dn, new_hrtin.timetag, ref_time)
                  new_hrtout  = fix_changes(rad_trap_out[hrtout_pos], shutter, t0, t1)
                  hrtout_dn = dethin_data(new_hrtout.dn, new_hrtout.timetag, ref_time)
                  delta_t = 362.0d6 
                  rad_pos=where(rad_trap_in.timetag ge t0 and rad_trap_in.timetag le t1,count)
                  if count gt 0 then begin
                      temp_t = [t0, rad_trap_in[rad_pos].timetag, t1]
                      t_gappos = WHERE(abs(temp_t[0:-2]- temp_t[1:-1]) GT delta_t,t_gapcount)
                      if t_gapcount gt 0 then begin
                          for i=0,t_gapcount-1 do begin
                              p=where(ref_time ge temp_t[t_gappos[i]] and $
                                      ref_time le temp_t[t_gappos[i]+1],count)
                              if count gt 0 then gap_hrt[p]=power[p] AND ABS(shutter_pos[p]-1)
                          endfor
                          ; get the HRT state from the plan database for this orbit
                          hrt_state = get_hrt_plan(t0,t1,instrument=instrument)
                          if size(hrt_state,/tname) eq 'STRUCT' then begin
                              p=where(gap_hrt eq 1,count)
                              ; only count when power is ON and shutter open
                              hrt_from_plan[p]=power[p] AND ABS(shutter_pos[p]-1)
                              if count gt 0 then begin
                                  in_dn  = dethin_data(hrt_state.hrtin.dn, hrt_state.hrtin.timetag, ref_time)
                                  out_dn  = dethin_data(hrt_state.hrtout.dn, hrt_state.hrtout.timetag, ref_time)
                                  hrtin_dn[p] = in_dn[p]
                                  hrtout_dn[p] = out_dn[p]
                              endif
                          endif else begin
                              ; no hrt state from the planned activities - skip this orbit
                              out_data[orb].gap_hrt = TOTAL(ABS(shutter_pos-1) AND power) * bintime
                              out_data[orb].abs_gap = TOTAL(ABS(shutter_pos-1) AND power) * bintime
                          endelse
                      endif
                  endif 
              endelse

              ; look for times when rad_trap_in and rad_trap_out are both in the same state
              ; only when the power is actually on and shutter is open
              same_pos = where(hrtin_dn eq hrtout_dn and power eq 1.0 and shutter_pos eq 0,hrt_bad_count)
              if hrt_bad_count gt 0 then begin
                  same_hrt[same_pos] = 1b
                  print,'Warning: rad_trap_in and rad_trap_out showed same state for '+$
                         strtrim(string(total(same_hrt)*bintime,format='(f10.1)'),2)+' seconds (get plan)'
                  ; get the state of the HRT from the plan to try to fix this
                  ; since instances when the rad_trap_out is sticky (when state of rad_trap_in changes)
                  hrt_state = get_hrt_plan(t0,t1,instrument=instrument)
                  if size(hrt_state,/tname) ne 'STRUCT' then begin
                      ; no planned activities were found for the timerange
                      ; skip this whole orbit
                      ; count only the times when the shutter is open
                      out_data[orb].gap_hrt= TOTAL(ABS(shutter_pos-1)) * bintime
                      out_data[orb].same_hrt = TOTAL(same_hrt) * bintime
                      out_data[orb].abs_gap = TOTAL(ABS(shutter_pos-1)) * bintime
                      continue
                  endif
                  new_hrtin  = fix_changes(hrt_state.hrtin, shutter, t0, t1)
                  hrtin_pos = where(new_hrtin.timetag ge t0 and new_hrtin.timetag le t1, hrtin_count)
                  in_dn  = dethin_data(new_hrtin.dn[hrtin_pos], new_hrtin.timetag[hrtin_pos], ref_time)
                  hrtin_dn[same_pos] = in_dn[same_pos]
                  new_hrtout = fix_changes(hrt_state.hrtout, shutter, t0, t1)
                  hrtout_pos = where(new_hrtout.timetag ge t0 and new_hrtout.timetag le t1, hrtout_count)
                  out_dn = dethin_data(new_hrtout.dn[hrtout_pos], new_hrtout.timetag[hrtout_pos], ref_time)
                  hrtout_dn[same_pos] = out_dn[same_pos]
              endif
           endif

           ; if we're here, the SC is not in safehold/contingency mode
           pos0 = where(fss0.timetag ge ref_time[0] and fss0.timetag le ref_time[-1],count0)
           pos1 = where(fss1.timetag ge ref_time[0] and fss1.timetag le ref_time[-1],count1)
           posv = where(fssv.timetag ge ref_time[0] and fssv.timetag le ref_time[-1],countv)
           if count0 lt 2 or count1 lt 2 or countv lt 2 then begin
               ; no SC data for this whole orbit (no fssangle info either)
               ; we will assume that the spacecraft is pointing dead-on
               ; and will account for the gap_sc in the uncertainties
               gap_sc[*] = 1b
               gap_sc_open = power AND ABS(shutter_pos-1)
               fss_keep[*] = 1b
           endif else begin
               ; only keep the data within the boxwidth
               fss0_eu = dethin_data(fss0[pos0].eu, fss0[pos0].timetag, ref_time)
               fss_valid = dethin_data(fssv[posv].dn, fssv[posv].timetag, ref_time)
               fss1_eu = dethin_data(fss1[pos1].eu, fss1[pos1].timetag, ref_time)
               fss_pos = where(fss_valid eq 1 and fss0_eu ge min_fss0 and fss0_eu le max_fss0 and $
                            fss1_eu ge min_fss1 and fss1_eu le max_fss1, fss_count)
               if fss_count gt 0 then fss_keep[fss_pos] = 1b

               ; look for gaps in the spacecraft HK telemetry (fssunangle)
               delta_t = 60.0d6
               fss0_time = [ref_time[0], fss0[pos0].timetag, ref_time[-1]]
               sc_gappos = WHERE(abs(fss0_time[0:-2] - fss0_time[1:-1]) GE delta_t,sc_gapcount)
               if sc_gapcount gt 0 then begin
                   ; for now just find out which plan was scheduled
                   plan = get_sorce_plan(insun_t0[orb], insun_t1[orb], /gps, instrument=instrument)
                   missing_plans[orb] = ptr_new(plan)
                   for i=0L,sc_gapcount-1L do begin
                       ; remove the gaps from the counted solar exposures
                       p=where(ref_time gt fss0_time[sc_gappos[i]] and $
                               ref_time lt fss0_time[sc_gappos[i]+1], count)
                       if count gt 0 then gap_sc[p]=1b
                       gap_sc_open[p] = byte(gap_sc[p]) AND byte(ABS(shutter_pos[p]-1))
                   endfor
               endif
           endelse

           if hrtpos eq 'OUT' then begin
               keep = where(gap_sci eq 0 and power eq 1 and shutter_pos eq 0 and fss_keep eq 1 and $
                            hrtout_dn eq 1, count)
           endif else if hrtpos eq 'IN' then begin
               keep = where(gap_sci eq 0 and power eq 1 and shutter_pos eq 0 and fss_keep eq 1 and $
                            hrtin_dn eq 1, count)
           endif else begin
               keep = where(gap_sci eq 0 and power eq 1 and shutter_pos eq 0 and fss_keep eq 1, count)
           endelse
           out_data[orb].solar_exp = count * bintime
           out_data[orb].gap_sci = total(gap_sci) * bintime
           out_data[orb].gap_sc = total(gap_sc) * bintime
           out_data[orb].gap_sc_open = total(gap_sc_open) * bintime
           out_data[orb].gap_hrt = total(gap_hrt) * bintime
           out_data[orb].same_hrt = total(same_hrt) * bintime
           out_data[orb].gap_power = total(gap_power) * bintime
           out_data[orb].hrt_from_plan = total(hrt_from_plan) * bintime
           ; ignore gaps when power is OFF (no sci telemetry comes in then) and when no sunvalid
           ; out_data[orb].abs_gap = total((gap_sci OR gap_sc_open) AND power AND fss_keep) * bintime

           ; I'm removing the gap_sc_open from the abs_gap since this overestimates the errors when
           ; only the spacecraft packets are missing and everything else is nominal
           out_data[orb].abs_gap = total(gap_sci AND power AND fss_keep) * bintime

           ; if a minSciGap was provided, fill the whole orbit with exposure data from
           ; the planned activity database
           if n_elements(minSciGap) gt 0 then begin
               if out_data[orb].gap_sci ge minSciGap then begin
                   ; get the solar exposure information for this WHOLE orbit from the plan
                   result = get_exp_from_plan(insun_t0[orb], insun_t1[orb], /gps, instrument=instrument, hrtpos=hrtpos,boxwidth=boxwidth)
                   if size(result,/tname) eq 'STRUCT' then begin
                       out_data[orb].gap_sci = 0.0d
                       out_data[orb].gap_sc = 0.0d
                       out_data[orb].gap_hrt = 0.0d
                       out_data[orb].gap_power = 0.0d
                       out_data[orb].abs_gap = 0.0d
                       out_data[orb].solar_exp = result.solar_exp
                       out_data[orb].time_from_plan = result.time_from_plan 
                       out_data[orb].time_std_from_plan = result.solar_exp_stdev 
                   endif

               endif
           endif

       endfor  ; orbit loop
       if not keyword_set(silent) then print,''

   endfor  ; orbit group loop

   ; clean up the data using the information we got from visual inspection
   ; of the data and inspection of the mission logs indicating the presence of
   ; dropped packets in the communication with the spececraft or the ground system.
   ; In these cases, we assume the plan was exectued as scheduled.

   ; first ignore any data prior to mission day 10.0
   p=where(out_data.endtime le 728248333020000d,count)
   if count gt 0 then begin
       out_data[p].solar_exp=0.0
       out_data[p].gap_sci=0.0
       out_data[p].abs_gap=0.0
   endif


   ;return,out_data

   CLEANUP:

   query_database,/dbclose

   ; clear the abs_gap and solar_exp when we know the instrument was turned off
   ; properly after a safehold (with missing telemetry)
   ; and when the shutter was closed when safehold was triggered - no gaps in exposure
   pos = where(out_data.starttime ge 863647476079997d and out_data.endtime le 863651640560004d,count)
   if count gt 0 then begin
       ; HRT is IN here
       out_data[pos].solar_exp=0.0
       out_data[pos].abs_gap=0.0
   endif
   ;pos=23181L
   pos = where(out_data.starttime ge 910020518030000d and out_data.endtime le 910024638050000d,count)
   if count gt 0 then begin
       ; HRT is IN here
       out_data[pos].solar_exp=0.0
       out_data[pos].abs_gap=0.0
   endif
   ;pos=32074L
   pos = where(out_data.starttime ge 915290177030000d and out_data.endtime le 915293977050000d,count)
   if count gt 0 then begin
       out_data[pos].solar_exp=0.0
       out_data[pos].abs_gap=0.0
   endif
   ;pos=[32156L,32157L]
   pos = where(out_data.starttime ge 915808827030000d and out_data.endtime le 915818347040000d,count)
   if count gt 0 then begin
       out_data[pos].solar_exp=0.0
       out_data[pos].abs_gap=0.0
   endif
   ;pos=[36247L,36248L]
   pos = where(out_data.starttime ge 939671597030000d and out_data.endtime le 939681107050000d,count)
   if count gt 0 then begin
       out_data[pos].solar_exp=0.0
       out_data[pos].abs_gap=0.0
   endif
   ;pos= 41366 -> 41434
   pos = where(out_data.starttime ge 969535647040000d and out_data.endtime le 969941652040000d,count)
   if count gt 0 then begin
       out_data[pos].solar_exp=0.0
       out_data[pos].abs_gap=0.0
   endif
   ;pos=41890L
   pos = where(out_data.starttime ge 972595492030000d and out_data.endtime le 972623881150000d,count)
   if count gt 0 then begin
       out_data[pos].solar_exp=0.0
       out_data[pos].abs_gap=0.0
   endif
   ;pos=[41913L, 41914L, 41915L]
   pos = where(out_data.starttime ge 972752838030000d and out_data.endtime le 972768248050000d,count)
   if count gt 0 then begin
       out_data[pos].solar_exp=0.0
       out_data[pos].abs_gap=0.0
   endif
   ;pos=42716L
   pos = where(out_data.starttime ge 977432638030000d and out_data.endtime le 977483108150000d,count)
   if count gt 0 then begin
       out_data[pos].solar_exp=0.0
       out_data[pos].abs_gap=0.0
   endif
   ;pos=[42751L, 42752L, 42753L]
   pos = where(out_data.starttime ge 977683194030000d and out_data.endtime le 977698574050000d,count)
   if count gt 0 then begin
       out_data[pos].solar_exp=0.0
       out_data[pos].abs_gap=0.0
   endif
   ;pos=[46505L, -> 46517
   pos = where(out_data.starttime ge 999639214030000d and out_data.endtime le 999712834050000d,count)
   if count gt 0 then begin
       out_data[pos].solar_exp=0.0
       out_data[pos].abs_gap=0.0
   endif

   ; this region between 84.9 and 85.91 had an under-voltage event which triggered a shut down
   ; (with closing the shutter as one of the first commands)
   ; the first part of the first orbit was ok (leave the solar_exp the same but 0 the abs_gap)
   pos = where(out_data.starttime ge 734738551040000d and out_data.endtime le 734742771040000d,count)
   if count gt 0 then begin
       out_data[pos].abs_gap=0.0
   endif
   pos = where(out_data.starttime ge 734744391040000d and out_data.endtime le 734824551040000d,count)
   if count gt 0 then begin
       out_data[pos].solar_exp=0.0
       out_data[pos].abs_gap=0.0
   endif


   if HRTPOS eq 'OUT' and instrument eq 'SimA' then begin
       ; the door was opened on mission day 18.87: anything before count as if
       ; HRT was in (since door is same material as HRT)
       pos = where(out_data.starttime lt 729043213000000d,count)
       if count gt 0 then begin
           ; door is still closed (act as if HRT was IN)
           out_data[pos].gap_sci = 0.0d
           out_data[pos].gap_sc = 0.0d
           out_data[pos].gap_hrt = 0.0d
           out_data[pos].gap_power = 0.0d
           out_data[pos].solar_exp=0.0
           out_data[pos].abs_gap=0.0
       endif

       ; add the gaps from the safeholds when we know the shutter was open
       ; when the safehold occured (there's no way to be certain the instrument
       ; was shutdown cleanly when a safehold occurs so we need to account for those)
       ;pos=23181L
       pos = where(out_data.starttime ge 910020518030000d and out_data.endtime le 910024638050000d,count)
       if count gt 0 then begin
           ; HRT is IN here
           out_data[pos].solar_exp=0.0
           out_data[pos].abs_gap=0.0
       endif

       ; shutter was closed when safehold was triggered - no gaps in exposure
       ;pos=32041L
       pos = where(out_data.starttime ge 915033688030000d and out_data.endtime le 915102769150000d,count)
       if count gt 0 then begin
           exposure=1460.0d
           out_data[pos].solar_exp=exposure
           out_data[pos].abs_gap=0.0
           out_data[pos].safehold=out_data[pos].insun_time - exposure
       endif

       ; shutter was opened when safehold was triggered - but rad_trap was in
       ;sd=1571.86
       pos = where(out_data.starttime ge 863210469040000d and out_data.endtime le 863221921204971d,count)
       if count gt 0 then begin
           out_data[pos].abs_gap=0.0
           out_data[pos].safehold=out_data[pos].gap_sci
       endif

   endif else if HRTPOS eq 'IN' and instrument eq 'SimA' then begin

       ;pos=lindgen(15)+1229L -> HRT was OUT for this - leave everything as 0.0

       ;pos=23181L
       pos = where(out_data.starttime ge 910020518030000d and out_data.endtime le 910024638050000d,count)
       if count gt 0 then begin
           ; HRT is IN here
           out_data[pos].solar_exp=1451.0 + 2116.8
           out_data[pos].abs_gap=148.0 + 3820.0
       endif

       ;pos=32041L
       pos = where(out_data.starttime ge 915033688030000d and out_data.endtime le 915102769150000d,count)
       if count gt 0 then begin
           ; HRT is OUT for this time
           exposure=0.0d
           out_data[pos].solar_exp=exposure
           out_data[pos].abs_gap=0.0
           out_data[pos].safehold=out_data[pos].insun_time - exposure
       endif

   endif else if instrument eq 'SimB' then begin

       ; the door was opened on mission day 18.87: anything before count as if
       ; HRT was in (since door is same material as HRT)
       pos = where(out_data.starttime lt 729043213000000d,count)
       if count gt 0 and HRTPOS eq 'OUT' then begin
           ; door is still closed
           out_data[pos].gap_sci = 0.0d
           out_data[pos].gap_sc = 0.0d
           out_data[pos].gap_hrt = 0.0d
           out_data[pos].gap_power = 0.0d
           out_data[pos].solar_exp=0.0
           out_data[pos].abs_gap=0.0
       endif

       ; on SD=24.4, SimB shutter was left open for many orbits with HRTIN (missing telemetry)
       ; and SimA and SimB were "locked" for several hours
       ; change the exposure times and remove abs_gap  and do not use the results from the plan
       pos=[]
       pos0 = where(out_data.starttime ge 729467601030000d and out_data.endtime le 729535581050000d,count)
       if count gt 0 then pos=pos0
       pos0 = where(out_data.starttime ge 729467601030000d and out_data.endtime le 729535581050000d,count)
       if count gt 0 then begin
           if HRTPOS eq 'IN' then $
               out_data[pos].solar_exp=out_data[pos].insun_time $
           else $
               out_data[pos].solar_exp=0.0
           out_data[pos].abs_gap=0.0
           out_data[pos].time_from_plan=0.0
           out_data[pos].time_std_from_plan=0.0
           out_data[pos].hrt_from_plan=0.0
       endif

       ;pos=23181L -> SimB shutter was closed before  safehold
       pos = where(out_data.starttime ge 910020518030000d and out_data.endtime le 910024638050000d,count)
       if count gt 0 then begin
           out_data[pos].solar_exp=0.0
           out_data[pos].abs_gap=0.0
       endif

       ;pos=32041L  -> SimB shutter was closed before  safehold
       pos = where(out_data.starttime ge 915033688030000d and out_data.endtime le 915102769150000d,count)
       if count gt 0 then begin
           out_data[pos].solar_exp=0.0
           out_data[pos].abs_gap=0.0
       endif

       ; shutter was closed during the CalSweepXXXX activities
       ;sd=582.5
       pos = where(out_data.starttime ge 777734568903998d and out_data.endtime le 777767869192013d,count)
       if count gt 0 then begin
           out_data[pos].abs_gap=0.0
           out_data[pos].safehold=out_data[pos].gap_sci
       endif

       ; shutter was closed when safehold was triggered
       ;sd=1571.86
       pos = where(out_data.starttime ge 863210469040000d and out_data.endtime le 863221921204971d,count)
       if count gt 0 then begin
           out_data[pos].abs_gap=0.0
           out_data[pos].safehold=out_data[pos].gap_sci
       endif

   endif
   
   ; correct for orbits where we know from the mission ops notes that there was a 
   ; problem with the download of the data but NOT with the execution of the plan 
   ; use the times from the plan and remove the uncertainties.
   pos=[]
   ; yd=2003056.9
   pos1=where(out_data.starttime ge 730244081030000d and out_data.endtime le 730247891050000d,count)
   if count gt 0 then pos=[pos,pos1]

   ; yd=2003151.8
   pos1=where(out_data.starttime ge 738444939030000d and out_data.endtime le 738460354050000d,count)
   if count gt 0 then pos=[pos,pos1]

   ; yd=2003224.2
   pos1=where(out_data.starttime ge 744695959030000d and out_data.endtime le 744700139050000d,count)
   if count gt 0 then pos=[pos,pos1]

   ; yd=2003253.1
   pos1=where(out_data.starttime ge 747199579030000d and out_data.endtime le 747203389050000d,count)
   if count gt 0 then pos=[pos,pos1]

   ; yd=2004204.8
   pos1=where(out_data.starttime ge 774556069030000d and out_data.endtime le 774571429050000d,count)
   if count gt 0 then pos=[pos,pos1]

   ; yd=2004221.7
   pos1=where(out_data.starttime ge 776020009030000d and out_data.endtime le 776099549050000d,count)
   if count gt 0 then pos=[pos,pos1]

   ; yd=2004267.8
   pos1=where(out_data.starttime ge 780003618630000d and out_data.endtime le 780095039370000d,count)
   if count gt 0 then pos=[pos,pos1]

   ; yd=2005121.6
   pos1=where(out_data.starttime ge 798993199030000d and out_data.endtime le 798996928950000d,count)
   if count gt 0 then pos=[pos,pos1]

   ; yd=2006192.1
   pos1=where(out_data.starttime ge 836620311900000d and out_data.endtime le 836629850490000d,count)
   if count gt 0 then pos=[pos,pos1]

   ; yd=2006194.9
   pos1=where(out_data.starttime ge 836865109030000d and out_data.endtime le 836938739150000d,count)
   if count gt 0 then pos=[pos,pos1]

   ; yd=2009082.5
   pos1=where(out_data.starttime ge 921847417030000d and out_data.endtime le 921868797050000d,count)
   if count gt 0 then pos=[pos,pos1]

   ; yd=2011186.7
   pos1=where(out_data.starttime ge 993923614030000d and out_data.endtime le 994003294050000d,count)
   if count gt 0 then pos=[pos,pos1]

    for orb=0,n_elements(pos)-1 do begin
        result = get_exp_from_plan(out_data[pos[orb]].starttime, out_data[pos[orb]].endtime, $
            /gps, instrument=instrument, hrtpos=hrtpos,boxwidth=boxwidth)
        if size(result,/tname) eq 'STRUCT' then begin
            out_data[pos[orb]].gap_sci = 0.0d
            out_data[pos[orb]].gap_sc = 0.0d
            out_data[pos[orb]].gap_hrt = 0.0d
            out_data[pos[orb]].gap_power = 0.0d
            out_data[pos[orb]].abs_gap = 0.0d
            out_data[pos[orb]].solar_exp = result.solar_exp
            out_data[pos[orb]].time_from_plan = result.time_from_plan 
            out_data[pos[orb]].time_std_from_plan = result.solar_exp_stdev 
        endif
    endfor




   if n_elements(savefile) ge 1 then begin
      sim_exp = temporary(out_data)
      save,sim_exp,file=savefile
      return,sim_exp
   endif

   return,out_data

END 

