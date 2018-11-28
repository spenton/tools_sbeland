;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Get the Distance and Doppler correction from the database for a requested
;   date or date range.
;
; CALLING SEQUENCE:
;   result = get_1au(starttime, stoptime, /gps)
;
; INPUT PARAMETERS:
;   startTime -
;      The start time for the query
;
;   stopTime -
;      The stop time for the query. If not specified, will return the
;      correction following starttime.
;
; OPTIONAL INPUT PARAMETERS:
;   instrumentModeId -
;      The instrument we want to get the 1AU for.  Defaults to SORCE.
;   gps - 
;      If set, indicates that input times are to be interpreted as
;      microseconds since Jan 6, 1980 midnight UT.
;   missionDays - 
;      If set, indicates that input times are to be interpreted as
;      days elapsed since launch, i.e. mission days. (DEFAULT) 
;   julianDays - 
;      If set, indicates that the input times are to be interpreted as
;      julian day numbers.
;
; RETURNED PARAMETERS:
;   A structure with the cross-correlation residuals for a list of CCD
;   offset positions and the best estimate of the offset (by fitting a Gaussian).
;
; OPTIONAL OUTPUT PARAMETERS:
;      NONE
;
; EXAMPLE:  
;
; COMMON BLOCKS:
;      NONE
;
; NOTES:
;
; REVISION HISTORY:
;   Revision: $Id: get_1au.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;-
;
function get_1au, startTime, stopTime, instrumentModeId=instrumentModeId, verbose=verbose, $
         gps=gps, missionDays=missionDays, julianDays=julianDays

    if n_elements(stopTime) eq 0 then stopTime=startTime
    if keyword_set(gps) then begin
       ; user specified time ingps microseconds
       t0 = startTime
       t1 = stopTime
    endif else if keyword_set(julianDays) then begin
       ; user specified time in julian days
       t0 = jd2gps(startTime)*1.d6
       t1 = jd2gps(stopTime)*1.d6
    endif else begin
       ; user specified timetags in mission day
       t0 = sd2gps(startTime)*1.d6
       t1 = sd2gps(stopTime)*1.d6
    endelse

    ; is start and stop times are the same add 1 minute to stop to make sure we get something back
    if stopTime eq startTime then stopTime+=90d6
    if n_elements(instrumentModeId) eq 0 then instrumentModeId=50

    ; correct for 1AU and Doppler
    ; the instrumentId for SORCE is 50 and we're at version=6  (querying for the max(version) is VERY SLOW)
    query2="SELECT microsecondsSinceGpsEpoch,sunObserverDistanceCorrection,sunObserverDopplerFactor "
    query2 = query2+"FROM SolarDistAndDopplerFixedStep2 WHERE instrumentModeId="+strtrim(string(instrumentModeId),2)+" and version=7 and "
    query2 = query2+"microsecondsSinceGpsEpoch>="+strtrim(string(ulong64(t0)),2)+" and "
    query2 = query2+"microsecondsSinceGpsEpoch<="+strtrim(string(ulong64(t1)),2)
    if keyword_set(verbose) then print,query2
    query_database, query2, solarDist, info
    if size(solarDist,/tname) ne 'STRUCT' then return,-1
    
    return, solarDist

end
