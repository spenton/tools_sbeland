;+
; NAME:   GET_ORBITS  
;
; PURPOSE: 
;   Determines the start and stop times of each orbit within a
;   specified time range.
;
; CALLING SEQUENCE:
;   out_data = GET_ORBITS(startTime, stopTime, simA=simA, simB=simB) 
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
;
; OUTPUT PARAMETERS:
;   out_data -
;      The resulting solar exposure time for the specified instrument calculated
;      per orbit.  The structure will have the following format:
;         instrument:  string - (sim_a or sim_b)
;         starttime:   double array - gps timestamp (microsecs) of start of each orbit (when entering day)
;         endtime:     double array - gps timestamp (microsecs) of end of each orbit (when entering night)
;         insun_time:  double array - total time in_sun for each orbit (seconds)
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
;   2012.01.26  SBeland
;   Revision: $Id: get_sim_orbits.pro,v 1.2 2018/11/26 16:16:02 spenton Exp $
;-
;*****************************************************************************
FUNCTION GET_SIM_ORBITS, startTime, stopTime, simA=simA, simB=simB,  $
         gps=gps, missionDays=missionDays, julianDays=julianDays, silent=silent

   ; will extract data for SIM_A by default
   instrument = "SimA"
   if (keyword_set(simB) and NOT keyword_set(simA)) then begin
      instrument = "SimB"
   endif

   if instrument eq "SimA" then begin
       TMID = {fssunangle0:    7871L, $
               fssunangle1:    7872L, $ 
               fssunanglevld:  7749L, $
               sunpresence:    8025L, $ 
               sunpresent:     7464L $ 
       }
       tlmid = {fssunangle0:   7871L, $
               fssunangle1:    7872L, $ 
               fssunanglevld:  7749L, $
               sunpresence:    8025L, $ 
               sunpresent:     7464L $ 
       }
   endif else begin
       TMID = {fssunangle0:    7871L, $
               fssunangle1:    7872L, $ 
               fssunanglevld:  7749L, $
               sunpresence:    8025L, $ 
               sunpresent:     7464L $ 
       }
       tlmid = {fssunangle0:   7871L, $
               fssunangle1:    7872L, $ 
               fssunanglevld:  7749L, $
               sunpresence:    8025L, $ 
               sunpresent:     7464L $ 
       }
   endelse

   ; get the start and stop time of every orbits within the requested time range
   if not keyword_set(silent) then print,' querying the database ...'
   ; get_sorce_telemetry, data, info, startTime, stopTime, TMID=sunpresence, gps=gps, $
   ;     missionDays=missionDays, julianDays=julianDays
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

   if not keyword_set(silent) then print,' determining the number of orbits ...'
   ;jstmt = fjava_get_jdbc_statement(user='sorce', password='sorcedb', server='sorce-db', database='SORCE')
   ;sorceDbExchange = OBJ_NEW("oJava_DbExchange", jstmt)
   ; q1="select Value 'dn', SCT_VTCW 'timetag' from SORCE_L1.dbo.TMdiscrete where TMID="+strtrim(TMID.sunpresence,2)+$
   ; q1="select Value 'dn', SCT_VTCW 'timetag' from SORCE_L1.dbo.TMdiscrete where TMID="+strtrim(TMID.fssunanglevld,2)+$
   ; the mu$fssunanglevld sometimes stops showing up (when mu is off) which messes up the orbits
   ; use ac$sunpresent instead
   q1="select Value 'dn', SCT_VTCW 'timetag' from SORCE_L1.dbo.TMdiscrete where TMID="+strtrim(TMID.sunpresent,2)+$
      " and SCT_VTCW BETWEEN "+strtrim(t0,2)+" AND "+strtrim(t1,2)+ " ORDER by SCT_VTCW"
   ;sunpresence=sorceDbExchange->getAllValues(q1)
   query_database,q1,sunpresence,nrows

   ; for mu$sunpresence
   ; IN_SUN=0
   ; OUT_SUN=1
   ; for ac$fssunanglevld or ac$sunpresent
   IN_SUN=1
   OUT_SUN=0
   
   ; get times for when SC enters sunpresence
   n_sunp = n_elements(SUNPRESENCE)
   if n_sunp lt 2 or size(SUNPRESENCE,/tname) ne 'STRUCT' then begin
       if not keyword_set(silent) then print,'No orbits found during the specified time range'
       return,-1
   endif

   ; set the sunpresence=2 to sunpresence=1 for early in the mission
   p=where(sunpresence.dn eq 2,count)
   if count gt 0 then sunpresence[p].dn=1

   ; remove the entries where the state is undefined
   p = where(SUNPRESENCE.dn ne OUT_SUN and SUNPRESENCE.dn ne IN_SUN,badcount,complement=comp)
   if badcount gt 0 then begin
       if not keyword_set(silent) then print,'Warning: '+strtrim(badcount,2)+' bad mu$sunpresence'
       SUNPRESENCE=temporary(SUNPRESENCE[comp])
       n_sunp -= badcount
   endif
   ; instead of removing the UNKNOWN sunpresence, force them to OUT 
   ; (otherwise the orbits get somewhat weird at the start of the mission)
   ; sunpresence.dn = sunpresence.dn < 1

   ; look for transitions from OUT to IN (start of orbit)
   p0=where(SUNPRESENCE[0L:n_sunp-2L].dn eq OUT_SUN AND $
            SUNPRESENCE[1L:n_sunp-1L].dn eq IN_SUN,ntrans0)
   if ntrans0 eq 0 then return,-1
   p0+=1L

   ; get times for when SC leaves sunpresence
   p1=where(SUNPRESENCE[0L:n_sunp-2L].dn eq IN_SUN AND $
            SUNPRESENCE[1L:n_sunp-1L].dn eq OUT_SUN,ntrans1)
   if ntrans1 eq 0 then return,-1
   p1+=1L

   start_pos=[]
   end_pos=[]
   i=0L
   while 1 do begin
       ; find the end of this orbit (longer than 20 minutes)
       k=where(sunpresence[p1].timetag gt sunpresence[p0[i]].timetag,count)
       if count eq 0 then break
       ; get the end of orbit just before a full orbit from this one
       k=where(sunpresence.timetag lt sunpresence[p0[i]].timetag+4800d6 and sunpresence.dn eq IN_SUN,count)
       if count eq 0 then break
       p=where(sunpresence[p1].timetag gt sunpresence[k[-1]].timetag,count)
       if count eq 0 then break
       start_pos=[start_pos,p0[i]]
       end_pos=[end_pos,p1[p[0]]]
       i=where(sunpresence[p0].timetag gt sunpresence[end_pos[-1]].timetag,count)
       if count eq 0 then break
       i=i[0]
   endwhile

   insun_t0 = SUNPRESENCE[start_pos].timetag
   insun_t1 = SUNPRESENCE[end_pos].timetag

   norbits = n_elements(insun_t0)
   if not keyword_set(silent) then print,' found '+strtrim(string(norbits),2)+' orbits ...'

   ; verify that the end time of every orbit is after its start time
   diff = where((insun_t1-insun_t0) le 0.0,count)
   if count gt 0 then begin
       print,'Error in determining orbits: '+strtrim(count,2)+' bad orbits '
       print,'     first bad orbit: ',insun_t0[diff[0]],insun_t1[diff[0]]
       return,-1
   endif

   ; initialize the output structure
   out_data = replicate( {instrument:instrument, starttime:0.0d, endtime:0.0d, insun_time:0.0d} , norbits)
   out_data.starttime = insun_t0
   out_data.endtime = insun_t1
   out_data.insun_time = (insun_t1-insun_t0)/1.0d6


   return,out_data

END 

