;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Return a structure containing information about missing 
;   science packets from either SimA or SimB from SORCE_L0:dbo.Packets
;   database table.
;
; CALLING SEQUENCE:
;   result = MISSING_PACKETS(startTime, stopTime, simA=sima, SimB=simb, /mission)
;
; INPUT PARAMETERS:
;   startTime -
;      The lower time range for which data will be returned.
;      May be specified as Julian days, gps microseconds (the default) 
;      or SORCE mission day numbers.
;
;   stopTime -
;      The upper time range for which data will be returned
;      May be specified as Julian days, gps microseconds (the default) 
;      or SORCE mission day numbers.
;
; OPTIONAL INPUT KEYWORDS:
;   simA, simB - 
;      Specifies which SIM channel to retreive. 
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
;   convert -
;      If set, will convert the returned start and stop times to the same
;      format as the input starttime and stoptime.
;   noSafe - 
;      If set, will check if the missing packets are due to a 
;      safehold/contingency conditions and will ignore these.
;
; RETURNED PARAMETERS:
;   An array of structures with the startTime of the previous and 
;   next science packet (this defines an approximate time range 
;   where th emissing packet would.
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
;    There are a lot of packets in the database, so the query can be very
;    slow if the requested time range is large: YOU'VE BEEN WARNED !
;-
;
;*****************************************************************************
FUNCTION GET_SAFEHOLD_TIMES, t0, t1, sorceDbExchange, tlmid, ref_time=ref_time
    ; return the amount of time the spacecraft was in safehold between t0 and t1

    q1 = "select SCT_VTCW 'timetag', Value 'dn' from SORCE_L1.dbo.TMdiscrete where TMID="
    q2 = " AND SCT_VTCW BETWEEN "+strtrim(ulong64(t0),2)+" AND " + $
        strtrim(ulong64(t1),2)+" ORDER by SCT_VTCW"
    q_str = q1+strtrim(string(tlmid),2)+q2
    acsmd=sorceDbExchange->getAllValues(q_str)
    ; check if empty
    if size(acsmd,/tname) eq 'OBJREF' then return,-1.0d

    ; acsmd samples at 300 seconds
    ; count the time when the SC is in normal operation
    if n_elements(ref_time) eq 0 then begin
        ; just return the total time in safe/contingency
        p=where(acsmd.dn ge 4,count,complement=comp)
        if count gt 0 then $
            return, ((t1 - t0)/1d6 - count* 300.0) > 0.0 $
        else $
            return,(t1 - t0)/1d6 > 0.0d
    endif else begin
        ; return an array with a value of 1 for when in safe/contingency
        p=where(acsmd.dn lt 4, count,comp=comp)
        if count eq 0 then return,bytarr(n_elements(ref_time))  ; always in normal mode
        if count eq n_elements(dn) then return,bytarr(n_elements(ref_time))+1b  ; always in safe/contingency mode
        acsmd.dn[p] = 1b     ; indicates safehold mode
        acsmd.dn[comp] = 0b  ; indicates normal mode
        safehold = dethin_data(acsmd.dn,acsmd.timetag, ref_time)
        return, safehold
    endelse

END
;*****************************************************************************
FUNCTION GET_POWERON_TIMES, t0, t1, sorceDbExchange, tlmid
    ; return the amount of time the spacecraft was in safehold between t0 and t1

    q1 = "select SCT_VTCW 'timetag', Value 'dn' from SORCE_L1.dbo.TMdiscrete where TMID="
    q2 = " AND SCT_VTCW BETWEEN "+strtrim(ulong64(t0),2)+" AND " + $
         strtrim(ulong64(t1),2)+" ORDER by SCT_VTCW"
    q_str = q1+strtrim(string(tlmid),2)+q2
    pri_power=sorceDbExchange->getAllValues(q_str)
    ; check if empty
    if size(pri_power,/tname) eq 'OBJREF' then return,-1.0d

    ; power normally samples at 300 seconds
    ; count the time when the instrument is ON
    if n_elements(ref_time) eq 0 then begin
        ; just return the total time when instrument is ON
        p=where(pri_power.dn gt 0,count)
        if count gt 0 then $
            return, count* 300.0 $
        else $
            return, 0.0d
    endif else begin
        ; return an array with a value of 1 for when instrument is ON
        p=where(pri_power.dn lt 1, count,comp=comp)
        if count eq 0 then return,bytarr(n_elements(ref_time))  ; always OFF
        if count eq n_elements(dn) then return,bytarr(n_elements(ref_time))+1b  ; always ON
        power = dethin_data(pri_power.dn,pri_power.timetag, ref_time)
        return, power
    endelse

END
;*****************************************************************************


function missing_packets, starttime, stoptime, sima=sima, simb=simb, $
             gps=gps, missionDays=missionDays, julianDays=julianDays, $
             convert=convert, showAll=showAll

    ; will extract data for SIM_A by default
    if keyword_set(simA) then instrument='SimA'
    if keyword_set(simB) then instrument='SimB'
    if n_elements(instrument) eq 0 then instrument='SimA'
    ; APID for the corresponding sim science packet
    if instrument eq 'SimA' then begin
        apid=498 
        acsmd=7926
        pri_power=20242
    endif else begin
        apid=514
        acsmd=7926
        pri_power=20243
    endelse

    if keyword_set(missionDays) then begin
        ; user specified time in mission (sorce) days
        t0 = sd2gps(startTime)
        t1 = sd2gps(stopTime)
    endif else if keyword_set(julianDays) then begin
        ; user specified time in julian days
        t0 = jd2gps(startTime)
        t1 = jd2gps(stopTime)
    endif else begin
        t0 = double(starttime)
        t1 = double(stoptime)
    endelse

    ; convert from seconds to micro seconds
    t0*=1d6
    t1*=1d6
 
    ; q1 = 'SELECT t1.gapID as startOfGroup, MIN(t2.gapID) ', $
    ;     'as endOfGroup FROM ', $
    ;     '(SELECT gapID FROM gaps tbl1 ', $
    ;     'WHERE NOT EXISTS(SELECT * FROM gaps tbl2 ', $
    ;     'WHERE tbl1.gapID = tbl2.gapID + 1)) t1 ', $
    ;     'INNER JOIN (SELECT gapID FROM gaps tbl1 ', $
    ;     'WHERE NOT EXISTS(SELECT * FROM gaps tbl2 ', $
    ;     'WHERE tbl2.gapID = tbl1.gapID + 1)) t2 ', $
    ;     'ON t1.gapID <= t2.gapID', $
    ;     'GROUP BY t1.gapID'
      
  
    jstmt = fjava_get_jdbc_statement(user='sorce', password='sorcedb', server='sorce-db', database='SORCE')
    sorceDbExchange = OBJ_NEW("oJava_DbExchange", jstmt)
    sql='SELECT SCT_VTCW, APID, SSCount, Length FROM SORCE_L0.dbo.Packets where SCT_VTCW>='+$
        strtrim(ulong64(t0),2)+' and SCT_VTCW<'+strtrim(ulong64(t1),2)+' and APID='+strtrim(apid,2)
    ;print, sql
    data=sorceDbExchange->getAllValues(sql)
 
    ; If no rows were returned, return a simple long 0
    nrows=n_elements(data)
    if nrows eq 0 or size(data,/tname) ne 'STRUCT' then return, 0L

    ; look for missing packets 
    ; SSCount is out of sequence indicates a hole
    ; SSCount is a 15bit counter so rolls back to 0 after 16383
    diff = abs(data[0:-2].sscount - data[1:-1].sscount)
    pos = where(diff gt 1 and diff ne 16383,count)
    if count eq 0 then return,0L

    ; create a structure with the the time of the previous and following packet
    starttime = dblarr(count)
    endtime = dblarr(count)
    SSCount0 = intarr(count)
    SSCount1 = intarr(count)
    k=0L
    BITMASK = UINT(16383)

    for i=0,count-1 do begin
        starttime[k] = data[pos[i]].SCT_VTCW+1.0
        endtime[k]   = data[pos[i]+1].SCT_VTCW-1.0
        p0 = (data[pos[i]].SSCount + 1L) AND BITMASK
        p1 = (data[pos[i]+1].SSCount - 1L) AND BITMASK
        SSCount0[k] = p0
        if p1 lt p0 and p1 lt BITMASK then begin
            ; we'll have to make a new entry starting at 0
            SSCount1[k] = 16383
            SSCount0 = [SSCount0[0:k],0,SSCount0[k+1:*]]
            SSCount1 = [SSCount1[0:k],p1,SSCount1[k+1:*]]
            starttime = [starttime[0:k],starttime[k],starttime[k+1:*]]
            endtime = [endtime[0:k],endtime[k],endtime[k+1:*]]
            k+=2
        endif else begin
            SSCount1[k] = p1
            k+=1
        endelse
    endfor

    if not keyword_set(showall) then begin
       ; remove from the list the entries when the SC was in safehold
       ; or when the instrument was turned off
       for i=0,count-1 do begin
           safe_time = get_safehold_times(starttime[i],endtime[i],sorceDbExchange,acsmd)
           pwr_on = get_poweron_times(starttime[i],endtime[i],sorceDbExchange,pri_power)
           if safe_time gt 0.0 or pwr_on eq 0.0 then begin
               ; mark this entry to be removed
               SSCount0[i] = -1
           endif
       endfor
       p=where(SSCount0 lt 0,count,complement=comp)
       if count gt 0 then begin
           starttime=starttime[comp]
           endtime=endtime[comp]
           SSCount0=SSCount0[comp]
           SSCount1=SSCount1[comp]
       endif
    endif

    out_data = {starttime:starttime, endtime:endtime, apid:replicate(apid,n_elements(sscount0)),$
        sscount0:sscount0, sscount1:sscount1}
    return, out_data

    if not keyword_set(convert) then return, data

    ; convert the time formats
    if keyword_set(missionDays) then begin
        ; user specified time in mission (sorce) days
        out_data=replicate({starttime:0.0d, stoptime:0.0d, activityname:''}, n_elements(data))
        out_data.starttime=vms2sd(jdbc2vms(data.starttime))
        out_data.stoptime=vms2sd(jdbc2vms(data.stoptime))
        out_data.activityname = data.activityname
    endif else if keyword_set(julianDays) then begin
        ; user specified time in julian days
        out_data=replicate({starttime:0.0d, stoptime:0.0d, activityname:''}, n_elements(data))
        out_data.starttime=vms2jd(jdbc2vms(data.starttime))
        out_data.stoptime=vms2jd(jdbc2vms(data.stoptime))
        out_data.activityname = data.activityname
    endif else if keyword_set(gps) then begin
        ; user specified timetags in gps microseconds
        out_data=replicate({starttime:0.0d, stoptime:0.0d, activityname:''}, n_elements(data))
        out_data.starttime=vms2gps(jdbc2vms(data.starttime))
        out_data.stoptime=vms2gps(jdbc2vms(data.stoptime))
        out_data.activityname = data.activityname
    endif else begin
        ; user specified timetags in gps microseconds
        out_data=TEMPORARY(data)
        out_data.starttime=jdbc2vms(out_data.starttime)
        out_data.stoptime=jdbc2vms(out_data.stoptime)
    endelse
 

    return, out_data

end
