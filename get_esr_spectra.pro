;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Return a structure containing the solar spectra measured from the
;   specified ESR and time span using the data already in the database.
;
; CALLING SEQUENCE:
;   spectra = GET_ESR_spectra(t0, t1, instrumentModeId, /missionDays, version=19)
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
;   instrumentModeId -
;      The instrument mode of interest:
;      31	SIM_A	ESR   
;      32	SIM_B	ESR  
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
;   A structure with the wavelength and  raw dn.
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
;   Revision: $Id: get_esr_spectra.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;-
;
function get_esr_spectra, startTime, stopTime, instrumentModeId, $
         gps=gps, missionDays=missionDays, julianDays=julianDays, $
         user=user, password=password, server=server, database=database, $
         no_duplicate=no_duplicate, version=version

    nrows=0
    ; validate the instrumentModeId
    if n_elements(instrumentModeId) eq 0 then begin
        doc_library,'get_sim_spectra'
        print,''
        print,'Missing instrumentModeId '
        return,-1
    endif else begin
        ; verify we got the right instrument
        valid_inst = [31,32]
        instrumentModeId = fix(instrumentModeId)
        if (where(valid_inst eq instrumentModeId))[0] lt 0 then  begin
            doc_library,'get_sim_spectra'
            print,''
            print,'Invalid instrumentModeId was provided'
            return,-1
        endif
    endelse

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

    if n_elements(version) eq 0 then version=19
    if n_elements(database) eq 0 then database='SORCE_SIM_V20'
    if n_elements(user) eq 0 then user='sbeland'
    if n_elements(password) eq 0 then password='sbeland1'
    if n_elements(dburl) eq 0 then dbUrl='jdbc:sybase:Tds:lasp-db-dev:4100/SORCE_SIM_V20'
    if n_elements(dbdriver) eq 0 then dbDriver='com.sybase.jdbc3.jdbc.SybDriver'



    q1="SELECT microsecondsSinceGpsEpoch as 'timetag', temperature FROM PrismDriveTemperature where "
    q1=q1+" instrumentModeId="+strtrim(string(instrumentModeId),2)+" AND "
    q1=q1+" version="+strtrim(string(version),2)+" AND "
    q3=" microsecondsSinceGpsEpoch>="+strtrim(ulong64(t0),2)+" AND microsecondsSinceGpsEpoch<="+strtrim(ulong64(t1),2)
    q3=q3+" ORDER by microsecondsSinceGpsEpoch ASC"
    query_database,q1+q3, temp_data, user=user, password=password, server=server, database=database

    q2="SELECT microsecondsSinceGpsEpoch as 'timetag', prismPosition as 'ccdpos', prismPosition as 'ccdpos_uncorr', "
    q2=q2+"realComponent as 'DN' FROM PhaseDetectedDN WHERE "
    q2=q2+" instrumentModeId="+strtrim(string(instrumentModeId),2)+" AND "
    q2=q2+" version="+strtrim(string(version),2)+" AND "
    query_database,q2+q3, spect_data, user=user, password=password, server=server, database=database
    
    if size(spect_data,/tname) ne 'STRUCT' then return,-1

    prismtemp=interpol(temp_data.temperature,temp_data.timetag, spect_data.timetag)
    ; for the expected structure
    spect=replicate({timetag:0d, ccdpos:0d, ccdpos_uncorr:0d, wavelength:0d, dn:0d, prismtemp:0d}, n_elements(spect_data))
    spect.timetag = spect_data.timetag
    spect.ccdpos = spect_data.ccdpos
    spect.ccdpos_uncorr = spect_data.ccdpos_uncorr
    spect.dn = spect_data.dn
    spect.prismtemp = prismtemp
    spect.wavelength = ccd2lambda(instrumentModeId, spect.ccdpos, spect.prismtemp)

    if NOT keyword_set(no_duplicate) then return, spect

    ; clean up the spectra of duplicates
    diff=spect[1:-1].ccdpos - spect[0:-2].ccdpos
    p =where(diff le 0d,count)
    if count gt 0 then begin
        ; remove duplicates
        spect[p].ccdpos=-1d
        spect[p+1].ccdpos=-1d
        k=where(spect.ccdpos gt 0d)
        spect=spect[k]
    endif

    s=sort(spect.ccdpos)
    q=uniq(spect[s].ccdpos)
    if n_elements(q) ne n_elements(s) then begin
        spect=spect[s[q]]
        s=sort(spect.timetag)
        spect=spect[s]
    endif

    return, spect

end
