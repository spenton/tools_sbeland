;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Get the RAW SOLSTICE telemetry and format it in a structure similar
;   to what get_sim_spectra.pro would return.
;
; CALLING SEQUENCE:
;   spectra = GET_SOLSTICE_spectra(t0, t1, externalElement, /missionDays)
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
;   externalElement -
;      solstice_a or solstice_b
;
; OPTIONAL INPUT PARAMETERS:
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
; RETURNED PARAMETERS:
;   A structure with the wavelength and irradiance (in raw dn).
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
;  You may want to run get_sorce_plan before hand to make sure an activity
;  with observations in the desired mode was performed during the time
;  span of interest.
;
; REVISION HISTORY:
;   Revision: $Id: get_solstice_spectra.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;-
;
function get_solstice_spectra, startTime, stopTime, instrumentModeId, $
         gps=gps, missionDays=missionDays, julianDays=julianDays

    nrows=0
    ; if externalElement, sola and solb are omitted, default to sola

    filter_scale=10d

    if instrumentModeId eq 7 then begin
        ; FUV A
        externalElement='solstice_a'
        detector='detector_b'
        pmtTemp='g_pmt_pc_temp'
    endif else if instrumentModeId eq 9 then begin
        ; MUV A
        externalElement='solstice_a'
        detector='detector_a'
        pmtTemp='f_pmt_pc_temp'
    endif else if instrumentModeId eq 11 then begin
        ; FUV B
        externalElement='solstice_b'
        detector='detector_b'
        pmtTemp='g_pmt_pc_temp'
    endif else if instrumentModeId eq 13 then begin
        ; MUV B
        externalElement='solstice_b'
        detector='detector_a'
        pmtTemp='f_pmt_pc_temp'
    endif
    items=[detector,'grat_pos','int_time',pmtTemp,'f_filter1_in','f_filter2_in']

    get_sorce_telemetry,soldata,info,startTime, stopTime, gps=gps, missionDays=missionDays, $
        julianDays=julianDays, externalElement=externalElement, item=items

    query_database,/reset

    ; form the data structure corresponding to GET_SIM_SPECTRA
    outdata=replicate({timetag:0d, ccdpos:0d, ccdpos_uncorr:0d, wavelength:0d, dn:0d, dn_tempcorr:0d, $
        irradiance:0d, prismtemp:0d, detectortemp:0d, filter1_in:0B, filter2_in:0B}, $
        n_elements((*soldata.SOLSTICE_A$DETECTOR_A).science))

    outdata.timetag = (*soldata.SOLSTICE_A$DETECTOR_A).science.timetag
    ; transform dn from counts to count rate (counts-per-seconds)  (int_time is in milli-sec)
    outdata.dn = (*soldata.SOLSTICE_A$DETECTOR_A).science.dn / ( (*soldata.SOLSTICE_A$INT_TIME).science.dn*1d-3)
    outdata.dn_tempcorr = outdata.dn
    outdata.ccdpos = (*soldata.SOLSTICE_A$grat_pos).science.dn
    outdata.ccdpos_uncorr = outdata.ccdpos
    outdata.wavelength = grt2lambda(outdata.ccdpos)
    outdata.detectortemp = interpol((*soldata.SOLSTICE_A$F_PMT_PC_TEMP).housekeeping.eu, $
        (*soldata.SOLSTICE_A$F_PMT_PC_TEMP).housekeeping.timetag, outdata.timetag)
    outdata.filter1_in = (*soldata.SOLSTICE_A$F_FILTER1_IN).science.eu
    outdata.filter2_in = (*soldata.SOLSTICE_A$F_FILTER2_IN).science.eu
    ; scale the data back when the ND1 filters are in place
    p=where(outdata.filter1_in eq 1,count)
    if count gt 0 then outdata[p].dn *=filter_scale
    p=where(outdata.filter2_in eq 1,count)
    if count gt 0 then outdata[p].dn *=filter_scale


    return, outdata

end
