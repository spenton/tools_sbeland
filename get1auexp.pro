;+
; Author: Stephane Beland
;
; PURPOSE: 
;
; CALLING SEQUENCE:
;
; INPUT PARAMETERS:
;
; OPTIONAL INPUT PARAMETERS:

; RETURNED PARAMETERS:
;
; OPTIONAL OUTPUT PARAMETERS:

; EXAMPLE:  
;
; COMMON BLOCKS:
;      NONE
;
; NOTES:

;
; REVISION HISTORY:
;   Revision: $Id: get1auexp.pro,v 1.2 2018/11/26 16:16:02 spenton Exp $
;------------------------------------------------------------------

function get1auexp, outfile=outfile

    ; get the SimSolarExposureData for SimA and SimB
    q1 = 'SELECT orbitStartTimeGps t0, orbitEndTimeGps t1, solarExposureHrtOut solarexp, '+$
         'cumulativeSolarExpHrtOut cumsolarexp from SimSolarExposureData where instrumentModeId=55'
    query_database, q1, solar55, info

    data=replicate({t0:0d, t1:0d, solarExp:0d, cumSolarExp:0d, oneAU:1d, correctedCumSolarExp:0d}, n_elements(solar55))
    data.t0=solar55.t0
    data.t1=solar55.t1
    data.solarExp=solar55.solarexp
    data.cumSolarExp=solar55.cumsolarexp

    ; unfortunately we can't use a single sql query to extract all of the 
    ; data since the timestamps for orbits in SimSOlarExposureData don't match the
    ; timestamps in the SolarDistance table

    ; THAT'S A LOT OF QUERIES !!!

    ; the instrumentId for SORCE is 50 and we're at version=6  (querying for the max(version) is VERY SLOW)
    query1="SELECT * FROM SolarDistAndDopplerFixedStep2 WHERE instrumentModeId=50 and version=7 and "
    cr=string(13b)
    for i=0L, n_elements(solar55)-1L do begin
        if (i MOD 10) eq 0 then print,format='($, I, A,TL100)',i,cr
        query2 = query1+" microsecondsSinceGpsEpoch>="+strtrim(string(ulong64(solar55[i].t0-120d6)),2)
        query2 = query2+" and microsecondsSinceGpsEpoch<="+strtrim(string(ulong64(solar55[i].t0+120d6)),2)
        query_database, query2, res, nrows
        if nrows gt 0 then data[i].oneau = res[0].sunObserverDistanceCorrection
    endfor

    correctedCumSolarExp = TOTAL(data.solarExp * data.oneAU, /cum)
    data.correctedCumSolarExp = correctedCumSolarExp

    print,''
    return,data

end
