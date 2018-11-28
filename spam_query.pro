; $Id: spam_query.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;+
; NAME:
;   spam_query.pro
;
; PURPOSE:
;   Return an array of Static Packet Monitor (SPaM) response data 
;   corresponding to GCI lockups, for the specified instrument 
;   and time period.
;
; INPUT PARAMETERS:
;   start_time - 
;      The lower bound of the time range for which data are requested. 
;      The default format is laTime (2008/012-01:23:45).
;
;   stop_time  - 
;      The upper bound of the time range for which data are requested, 
;      The default format is laTime (2008/012-01:23:45).
;
; OPTIONAL INPUT KEYWORDS:
;   instrument - 
;      The instrument for which data are requested.  Permitted values are:
;      SIM_A, SIM_B, SOLSTICE_A, SOLSTICE_B, TIM, XPS, ALL. If not provided,
;      will return information for all SIM instruments.
;   gps - 
;      If set, indicates that input times are to be interpreted as
;      microseconds since (1980/01/06 midnight UT)
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
;      The resulting structure will contain information on when a GCI lockup    
;      happened for the requested instrument.  These are determined from a flag
;      set by the MU as a non-response to a ping by the instrument.
;      RESPONSE_TIME_GPS  DOUBLE     8.1484838e+08     GPS when inst was restarted
;      RESPONSE_TIME_LA   STRING  '2005/305-02:46:10'  LA when inst was restarted
;      LOCKUP_TIME_GPS    DOUBLE     8.1484695e+08     GPS when lockup started
;      LOCKUP_TIME_LA     STRING  '2005/305-02:22:14'  LA  when lockup started
;      LOCKUP_DURATION    LONG            1435         Lockup duration (seconds)
;      PACKET             STRING  'SIM_A_HOUSEKEEPING' Telemetry that caused the lockup
;      RESPONSE_COUNT     LONG             670         SPM response count for this APID
;
; OPTIONAL OUTPUT PARAMETERS:
;      NONE
;
; COMMON BLOCKS:
;      NONE
;
; AUTHOR:
;   Stephane Beland 2011/12/07 based on Chris Pankratz get_spam_responses.pro
;   which was based on David Gathright's spam_query.cgi perl script
;   (in svn MODS/web/projects/sorce/web)
;
; VERSION:
;   $Id: spam_query.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;
; NOTES:
;
; USAGE EXAMPLE:
;   la1 = '2008/030-04:42:57'
;   la2 = '2008/044-04:42:57'
;   results = spam_query(la1, la2)
;   results = spam_query(1082.0, 1083.0, /mission, instr='sim_a')
;
;
;-
;*****************************************************************************
;
function spam_query, start_time, stop_time, instrument=instrument, DEBUG=debug, $
        gps=gps, missionDays=missionDays, julianDays=julianDays

    if n_elements(instrument) eq 0 then instrument_name='all' else $
        instrument_name = strlowcase(string(instrument))

    case instrument_name of
      'all': st_apid='ax'
      'sim_a': st_apid='ax1f'
      'sim_b': st_apid='ax20'
      'solstice_a': st_apid='ax1d'
      'solstice_b': st_apid='ax1e'
      'tim': st_apid='ax21'
      'xps': st_apid='ax1c'
      else:  begin
              doc_library, 'spam_query'
              return, -1
             end
    endcase

    if not keyword_set(debug) then debug=0

    ; convert the input time to a gps microsecond string
    if keyword_set(missionDays) then begin
        sgps1 = string(sd2gps(start_time) * 1.d6, format='(i18)')
        sgps2 = string(sd2gps(stop_time) * 1.d6, format='(i18)')
    endif else if keyword_set(julianDays) then begin
        sgps1 = string(jd2gps(start_time) * 1.d6, format='(i18)')
        sgps2 = string(jd2gps(stop_time) * 1.d6, format='(i18)')
    endif else if keyword_set(gps) then begin
        sgps1 = string(start_time * 1.d6, format='(i18)')
        sgps2 = string(stop_time * 1.d6, format='(i18)')
    endif else begin
        ; default expected format is laTime
        sgps1 = string(la2gps(start_time) * 1.d6, format='(i18)')
        sgps2 = string(la2gps(stop_time) * 1.d6, format='(i18)')
    endelse

    if (debug) then print, sgps1, sgps2

    ; prepare the database connection
    jstmt = fjava_get_jdbc_statement(user=name, password=password, dburl=dbURL, $
        dbdriver=dbdriver, server=server, database="SORCE")
    oJavaDbExchange = OBJ_NEW("oJava_DbExchange", jstmt)

    ; query for spam busy times
    query1 = "SELECT td.SCT_VTCW, tn.Mnemonic FROM SORCE_L1..TMnames tn, SORCE_L1..TMdiscrete td " + $
      "WHERE tn.Subsystem = 'spm' AND tn.Mnemonic LIKE '"+st_apid+"%_rsp_busy' " + $
      "AND td.Value = 1 AND td.TMID = tn.TMID AND td.SCT_VTCW BETWEEN " + $
      sgps1 + " and " + sgps2 + " ORDER BY td.SCT_VTCW, tn.Mnemonic"

    if (debug) then print, query1

    responses=oJavaDbExchange->getAllValues(query1)
    nresponses=n_elements(responses)

    rec = {response_time_gps:0.0d0, $
                     response_time_la:'', $
                     lockup_time_gps:0.0d0, $
                     lockup_time_la:'', $
                     lockup_duration:0L, $
                     packet:'', $
                     response_count:0L $
                     }

    if size(responses,/tname) eq 'OBJREF' then return,rec

    return_result = replicate(rec, nresponses)

    ; cycle through each spam response
    for i = 0, nresponses - 1 do begin
       apid = strmid(responses[i].mnemonic, 0, 5)
       vtcw = responses[i].sct_vtcw
       svtcw = string(ulong64(vtcw))
       
       query2 = "SELECT tn.Mnemonic, ta.Value FROM SORCE_L1..TMnames tn, SORCE_L1A..TManalog ta " + $
           "WHERE tn.Subsystem = 'spm' AND tn.Mnemonic IN ('" + apid + "_rsp_ct', '" + $
           apid + "_rsp_start', '" + apid + "_pkt_age', 'det_sec' ) " + $
           "AND ta.SCT_VTCW = " + svtcw + " AND tn.TMID = ta.TMID order by Mnemonic"
           if (debug) then print, query2
       
       response_data=oJavaDbExchange->getAllValues(query2)
       nrows=n_elements(response_data)
       
       if (debug) then help, /st, response_data
       
       ; determine the name of the packet that became stale   
       case apid of
          "ax211": packet_name = "TIM_HOUSEKEEPING"
          "ax212": packet_name = "TIM_SCIENCE"
          "ax213": packet_name = "TIM_GCI_ERROR"
          "ax215": packet_name = "TIM_DSP_CMD_ACCEPT"
          "ax216": packet_name = "TIM_DSP_CMD_REJECT"
          "ax217": packet_name = "TIM_DSP_MEM_DUMP"
          "ax1c1": packet_name = "XPS_HOUSEKEEPING"
          "ax1c2": packet_name = "XPS_SCIENCE"
          "ax1c3": packet_name = "XPS_GCI_ERROR"
          "ax1f1": packet_name = "SIM_A_HOUSEKEEPING"
          "ax1f2": packet_name = "SIM_A_SCIENCE"
          "ax1f3": packet_name = "SIM_A_GCI_ERROR"
          "ax1f5": packet_name = "SIM_A_DSP_CMD_ACCEPT"
          "ax1f6": packet_name = "SIM_A_DSP_CMD_REJECT"
          "ax1f7": packet_name = "SIM_A_DSP_MEM_DUMP"
          "ax201": packet_name = "SIM_B_HOUSEKEEPING"
          "ax202": packet_name = "SIM_B_SCIENCE"
          "ax203": packet_name = "SIM_B_GCI_ERROR"
          "ax205": packet_name = "SIM_B_DSP_CMD_ACCEPT"
          "ax206": packet_name = "SIM_B_DSP_CMD_REJECT"
          "ax207": packet_name = "SIM_B_DSP_MEM_DUMP"
          "ax1d1": packet_name = "SOLSTICE_A_HOUSEKEEPING"
          "ax1d2": packet_name = "SOLSTICE_A_SCIENCE"
          "ax1d3": packet_name = "SOLSTICE_A_GCI_ERROR"
          "ax1d4": packet_name = "SOLSTICE_A_DSP_PSR"
          "ax1d5": packet_name = "SOLSTICE_A_DSP_CMD_ACCEPT"
          "ax1d6": packet_name = "SOLSTICE_A_DSP_CMD_REJECT"
          "ax1d7": packet_name = "SOLSTICE_A_DSP_MEM_DUMP"
          "ax1e1": packet_name = "SOLSTICE_B_HOUSEKEEPING"
          "ax1e2": packet_name = "SOLSTICE_B_SCIENCE"
          "ax1e3": packet_name = "SOLSTICE_B_GCI_ERROR"
          "ax1e4": packet_name = "SOLSTICE_B_DSP_PSR"
          "ax1e5": packet_name = "SOLSTICE_B_DSP_CMD_ACCEPT"
          "ax1e6": packet_name = "SOLSTICE_B_DSP_CMD_REJECT"
          "ax1e7": packet_name = "SOLSTICE_B_DSP_MEM_DUMP"
       endcase

       tags = tag_names(response_data)

       return_result[i].response_time_gps = vtcw/1d6
       return_result[i].response_time_la = gps2la(vtcw/1.d6)
       return_result[i].lockup_duration = response_data[where(STRMATCH(response_data.mnemonic, '*pkt_age'))].value
       return_result[i].lockup_time_gps = response_data[where(STRMATCH(response_data.mnemonic, '*det_sec'))].value - $
                                          return_result[i].lockup_duration
       return_result[i].lockup_time_la = gps2la(return_result[i].lockup_time_gps)
       return_result[i].packet = packet_name
       return_result[i].response_count = response_data[where(STRMATCH(response_data.mnemonic, '*rsp_ct'))].value
       
    endfor

    return, return_result
end
