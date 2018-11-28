;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Within the known time range, track the diode responsivity with its
;   temperature.  The time range corresponds to the recovery time after
;   the OBC anomaly in 2007 when the instrument was slowly warming up.
;
; CALLING SEQUENCE:
;   result = GET_DIODE_tempcorr(instrumentModeId, starttime=starttime, stoptime=stoptime)
;
; INPUT PARAMETERS:
;   instrumentModeId -
;      41 for VIS1, 43 for UV or 44 for IR (SimA)
;
; OPTIONAL INPUT PARAMETERS:
;   startTime -
;      The lower time range for which data will be returned.
;      May be specified as julian days, gps microseconds, or
;      SORCE mission day numbers.
;
;   stopTime -
;      The upper time range for which data will be returned
;      May be specified as julian days, gps microseconds, or
;      SORCE mission day numbers.
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
;   
;
; OPTIONAL OUTPUT PARAMETERS:
;      NONE
;
; EXAMPLE:  
;
; COMMON BLOCKS:
;      NONE
;
; REVISION HISTORY:
;   Revision: 
;-
;
function get_diode_tempcorr_hybrid, instrumentModeId, version=version, alldata=alldata, noplot=noplot

    if instrumentModeId lt 41 or instrumentModeId gt 45 then begin
        print,'Error: invalid instrumentModeId (should be between 41 and 48)'
        return,-1
    endif

    if n_elements(version) eq 0 then version=24

    ; check if we already have the data on hand
    if size(alldata,/tname) ne 'POINTER' then begin

        ; get the list of plans to process and the corresponding ccdshift to apply to each spectra
        file='~/SORCE/data/sima_ir_plans_ccdshift_4110_4200.sav'
        restore,file

        nplans=n_elements(plans)
        ;if size(alldata,/tname) ne 'POINTER' then begin

            alldata=ptrarr(nplans)
            for i=0L,nplans-1 do begin
                ; read each spectra one at a time
                fitCoeff = reform(ccd_shifts.ccdfit[i,*])
                spect = get_sim_spectra(plans[i].starttime, plans[i].stoptime, instrumentModeId, /mission, version=version, $
                    fitCoeff=fitCoeff, profile_data=profile_data) ;, /no1au, /nodark)
                alldata[i]=ptr_new(spect,/no_copy)
            endfor

        ;endif
    endif


    return, alldata
   
end
