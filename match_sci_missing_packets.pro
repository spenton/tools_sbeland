;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Loops through the list of missing packets and look if the science 
;   telemetry is actually missing from that time period.
;
; CALLING SEQUENCE:
;   sci_tlm = MATCH_SCI_MISSING_PACKETS(infile=infile)
;
; INPUT PARAMETERS:
;   infile -
;      The name of the IDL save file with the missing packets structure.
;      (Generated from MISSING_PACKETS.PRO)
;
; RETURNED PARAMETERS:
;   An array with the time (seconds) for which science telemetry was found 
;   for each missing packet time range in the input file.
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
;-
;
function match_sci_missing_packets, infile=infile

    if n_elements(infile) eq 0 then begin
        infile=dialog_pickfile(filter='*.sav')
        ; check if we canceled
        if strlen(infile) eq 0 then return,0
    endif

    sObj=OBJ_NEW('IDL_Savefile',infile)
    sNames = sObj->Names()
    sObj->Restore, sNames[0]

    ; make a local copy of the restored data (since we don't know the name of the variable)
    data = SCOPE_VARFETCH(sNames[0], level=0)

    ; the expected structure is:
    ;     STARTTIME DOUBLE Array
    ;     ENDTIME   DOUBLE Array
    ;     APID      INT    Array
    ;     SSCOUNT0  INT    Array
    ;     SSCOUNT1  LONG   Array
    if size(data,/tname) ne 'STRUCT' then begin
        print,'Error: unrecognized data structure'
        help,data
        return,-1
    endif

    ; select the correct database for the shutter depending on the instrument
    if data.apid[0] eq 498 then begin
        instrument='SimA'
        shutter_db='SORCE_L1S.dbo.SimAScienceSamples'
    endif else if data.apid[0] eq 514 then begin
        instrument='SimB'
        shutter_db='SORCE_L1S.dbo.SimBScienceSamples'
    endif else begin
        print,'Error: unrecognized APID='+strtrim(data.apid[0],2)
        return,-1
    endelse
        
    n_missing = n_elements(data.starttime)
    out_data=lonarr(n_missing)

    jstmt = fjava_get_jdbc_statement(user='sorce', password='sorcedb', $
                                     server='sorce-db', database='SORCE')
    sorceDb = OBJ_NEW("oJava_DbExchange", jstmt)
 
    for i=0L,n_missing-1L do begin
        if ((i+1L) MOD 5) eq 0 or i eq n_missing-1 then $
           print,format='($, TL1,"   '+strtrim(string(i+1),2)+' of '+strtrim(string(n_missing),2)+string(13b)+'     ")'
        q_str = "select sampleVtcw 'timetag',shutterPosition 'dn' from "+shutter_db+" where " + $
             "sampleVtcw BETWEEN "+strtrim(ulong64(data.starttime[i]),2)+" AND " + $
             strtrim(ulong64(data.endtime[i]),2)+" ORDER by sampleVtcw"
        shutter=sorceDb->getAllValues(q_str)
        if size(shutter,/tname) ne 'STRUCT' then continue
        if n_elements(shutter) le 1 then continue

        ; shutter is sampled at 1.0 Hz or missing when instrument is off 
        ; or when the packet is missing
        diff = ABS(shutter[0:-2].timetag - shutter[1:-1].timetag)/1d6
        p=where(diff ge 0.0 and diff le 1.1,count)
        if count gt 0 then out_data[i]=count
    endfor

    print,''
    return, out_data

end
