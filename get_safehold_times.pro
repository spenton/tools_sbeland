;+
; NAME:   GET_SAFEHOLD_TIMES
;
; PURPOSE: 
;    Get the time when the SIM instrument was in safehold within
;    the requested time range.
;
; CALLING SEQUENCE:
;
; INPUT PARAMETERS:
;
; OPTIONAL INPUT KEYWORDS:
;
; OUTPUT PARAMETERS:
;
; OPTIONAL OUTPUT PARAMETERS:
;
; EXAMPLE:  
;
; COMMON BLOCKS:
;
; NOTE:
;
; REVISION HISTORY:
;   2011.11.03  SBeland
;   Revision: $Id: get_safehold_times.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;-
;*****************************************************************************
FUNCTION GET_SAFEHOLD_TIMES, t0, t1, tlmid=tlmid, dbex=dbex, ref_time=ref_time
    ; return the amount of time the spacecraft was in safehold between t0 and t1

    ;if n_elements(dbex) eq 0 then begin
    ;    ; the database connection needs to be extablished
    ;    jstmt = fjava_get_jdbc_statement(user='sorce', password='sorcedb', server='sorce-db', database='SORCE')
    ;    dbex = OBJ_NEW("oJava_DbExchange", jstmt)
    ;endif

    ; use the ACSMD telemetry id by default to determine a safehold condition
    if n_elements(tlmid) eq 0 then tlmid=7926L

    q1 = "select SCT_VTCW 'timetag', Value 'dn' from SORCE_L1.dbo.TMdiscrete where TMID="
    q2 = " AND SCT_VTCW BETWEEN "+strtrim(ulong64(t0),2)+" AND " + $
        strtrim(ulong64(t1),2)+" ORDER by SCT_VTCW"
    q_str = q1+strtrim(string(tlmid),2)+q2
    ;acsmd=dbex->getAllValues(q_str)
    ; if already connected the login info is ignored (version 2.4)
    query_database, q_str, acsmd, nrows, user='sorce', password='sorcedb', server='sorce-db', database='SORCE'
    ; check if empty
    if size(acsmd,/tname) eq 'OBJREF' then return,-1.0d

    ; acsmd samples at 300 seconds
    ; count the time when the SC is in normal operation
    if n_elements(ref_time) eq 0 then begin
        ; just return the total time in safe/contingency
        dn=acsmd.dn
        p=where(dn lt 4, count,comp=comp)
        if count eq 0 then return,0.0  ; always in normal mode
        if count eq n_elements(dn) then return,(t1 - t0)/1d6 > 0.0d ; always in safe mode
        ; count the time in safehold
        dn[p] = 1b     ; indicates safehold mode
        dn[comp] = 0b  ; indicates normal mode
        dn[(p-1)>0] = 1b   ; the previous (dn ge 4) will indicate start of safehold
        safe_time = 0.0
        flag=0
        for i=0,n_elements(dn)-1 do begin
           if dn[i] eq 1 and flag eq 0 then begin
               st0=acsmd[i].timetag
               flag=1
           endif else if dn[i] eq 0 and flag eq 1 then begin
               safe_time += ((acsmd[i].timetag - st0)/1d6)
               flag=0
           endif
        endfor
        if flag eq 1 then safe_time += ((acsmd[-1].timetag - st0)/1d6)
        return, safe_time
    endif else begin
        ; return an array with a value of 1 for when in safe/contingency
        dn=acsmd.dn
        p=where(dn lt 4, count,comp=comp)
        if count eq 0 then return,bytarr(n_elements(ref_time))  ; always in normal mode
        if count eq n_elements(dn) then return,bytarr(n_elements(ref_time))+1b  ; always in safe/contingency mode
        dn[p] = 1b     ; indicates safehold mode
        dn[comp] = 0b  ; indicates normal mode
        dn[(p-1)>0] = 1b   ; the previous (dn ge 4) will indicate start of safehold
        safehold = dethin_data(dn,acsmd.timetag, ref_time)
        return, safehold
    endelse

END


