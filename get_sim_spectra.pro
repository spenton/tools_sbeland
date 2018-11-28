;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Return a structure containing the solar spectra measured from the
;   specified instrumentModeId and time span.
;
; CALLING SEQUENCE:
;   spectra = GET_SIM_spectra(t0, t1, instrumentModeId, /missionDays)
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
;      31	SIM_A	ESR   -> not yet implemented
;      32	SIM_B	ESR   -> not yet implemented
;      41	SIM_A	VIS1
;      42	SIM_A	VIS2
;      43	SIM_A	UV
;      44	SIM_A	IR
;      45	SIM_B	VIS1
;      46	SIM_B	VIS2
;      47	SIM_B	UV
;      48	SIM_B	IR   -> not yet implemented
;
; OPTIONAL INPUT PARAMETERS:
;   notempCorr - 
;      If specified, will skip the prism temperature correction when 
;      assigning a wavelength, returning the interpolated value only.
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
;   fitCoeff -
;      If provided, the provided polynomial coefficients are used to calculate
;      the ccd spot positions before wavelengths are calculated.
;   profileDate -
;      If provided, the SimProfileIntegralCal used will have an effectiveDate 
;      less than this date (IN GPS microseconds). Allows to use the ProfileIntegral
;      from the early part of the mission in the later part - ignoring embedded
;      CCD shift offsets).
;   tsis - 
;      If set, will use the TSIS Sellmeier coefficients to determine the wavelength
;      change with temperature.
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
;   Revision: $Id: get_sim_spectra.pro,v 1.4 2018/11/26 16:16:02 spenton Exp $
;-
;
function get_sim_spectra, startTime, stopTime, instrumentModeId, $
         gps=gps, missionDays=missionDays, julianDays=julianDays, $
         noTempCorr=noTempCorr, noDark=noDark, version=version, $
         fitCoeff=fitCoeff, noCcdCorr=noCcdCorr, no1AU=no1AU, dn_only=dn_only, $
         profileDate=profileDate, verbose=verbose, profile_data=profile_data, $
         solarDist=solarDist, uncorr=uncorr, dettempcorr=dettempcorr, tsis=tsis

    run_t0 = systime(/sec)

    nrows=0
    ; validate the instrumentModeId
    if n_elements(instrumentModeId) eq 0 then begin
        doc_library,'get_sim_spectra'
        print,''
        print,'Missing instrumentModeId '
        return,-1
    endif else begin
        ; verify we got the right instrument
        valid_inst = [41,42,43,44,45,46,47,48]
        instrumentModeId = fix(instrumentModeId)
        if (where(valid_inst eq instrumentModeId))[0] lt 0 then  begin
            doc_library,'get_sim_spectra'
            print,''
            print,'Invalid instrumentModeId was provided'
            return,-1
        endif
    endelse

    ; define the telemetry items according to the specified mode
    case instrumentModeId of
        31: begin
            ; SIMA ESR
            instrument='SIMA'
            instrumentName='sim_a'
            prismId = 101862L
            prismTempId = 101909L
            detectorId = 100058L
            shutterId = 102121L
            pos0=4000L
            pos1=55590L
            intg_period=101395L
            detectorTemp=100170L
        end
        41: begin
            ; SIMA VIS1
            instrument='SIMA'
            instrumentName='sim_a'
            prismId = 101862L
            prismTempId = 101909L
            detectorId = 100043L
            shutterId = 102121L
            pos0=4120L
            pos1=25690L
            intg_period=101395L
            detectorTemp=100778L
        end
        42: begin
            ; SIMA VIS2
            instrument='SIMA'
            instrumentName='sim_a'
            prismId = 101862L
            prismTempId = 101909L
            detectorId = 100055L
            shutterId = 102121L
            pos0=4000L
            pos1=22650L
            intg_period=101395L
            detectorTemp=100963L
        end
        43: begin
            ; SIMA UV
            instrument='SIMA'
            instrumentName='sim_a'
            prismId = 101862L
            prismTempId = 101909L
            detectorId = 100046L
            shutterId = 102121L
            pos0=5900L
            pos1=50322L
            intg_period=101395L
            detectorTemp=101069L
        end
        44: begin
            ; SIMA IR
            instrument='SIMA'
            instrumentName='sim_a'
            prismId = 101862L
            prismTempId = 101909L
            detectorId = 100052L
            shutterId = 102121L
            pos0=16740L
            pos1=22850L
            intg_period=101395L
            detectorTemp=102112L
        end
        32: begin
            ; SIMB ESR
            instrument='SIMB'
            instrumentName='sim_b'
            prismId = 100478L
            prismTempId = 102212L
            detectorId = 100048L
            shutterId = 100197L
            pos0=6840L
            pos1=56000L
            intg_period=101647L
            detectorTemp=100939L
        end
        45: begin
            ; SIMB VIS1
            instrument='SIMB'
            instrumentName='sim_b'
            prismId = 100478L
            prismTempId = 102212L
            detectorId = 100044L
            shutterId = 100197L
            pos0=36730L
            pos1=56000L
            intg_period=101647L
            detectorTemp=101901L
        end
        46: begin
            ; SIMB VIS2
            instrument='SIMB'
            instrumentName='sim_b'
            prismId = 100478L
            prismTempId = 102212L
            detectorId = 100056L
            shutterId = 100197L
            pos0=40100L
            pos1=56000L
            intg_period=101647L
            detectorTemp=102062L
        end
        47: begin
            ; SIMB UV
            instrument='SIMB'
            instrumentName='sim_b'
            prismId = 100478L
            prismTempId = 102212L
            detectorId = 100057L
            shutterId = 100197L
            pos0=12330L
            pos1=55870L
            intg_period=101647L
            detectorTemp=101580L
        end
        48: begin
            ; SIMB IR
            instrument='SIMB'
            instrumentName='sim_b'
            prismId = 100478L
            prismTempId = 102212L
            detectorId = 100047L
            shutterId = 100197L
            pos0=39610L
            pos1=45600L
            intg_period=101647L
            detectorTemp=100698L
        end
    endcase

    SHT_OPEN = 0
    SHT_CLOSE = 1

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

    if keyword_set(verbose) then $
        print,'Processing '+strtrim(string(instrumentModeId),2)+' from ',gps2sd(t0/1d6),' to ',gps2sd(t1/1d6)

    ; get the list of prism positions for the specified time range (+/- 2 minutes to get darks)
    query_database,/reset
    get_sorce_telemetry, tlm_data, info, t0-120.0d6, t1+120.0d6, /gps, $
        tlmId=[prismId,prismTempId,detectorId,shutterId,intg_period,detectorTemp]
    if size(tlm_data,/tname) ne 'STRUCT' then begin
        print,'No data found for the specified instrumentId and time range'
        return,-1
    endif else if (*tlm_data.(0)).n_science eq 0 then begin
        print,'No data found for the specified instrumentId and time range'
        return,-1
    endif

    run_t1=systime(/sec)
    if keyword_set(verbose) then print,string(run_t1-run_t0,format='("Telem read ",F0.3," seconds")')

    ; convert the prism temperature
    if size((*tlm_data.(1)).housekeeping, /tname) eq 'STRUCT' then begin
        ; get_sorce_telemetry uses a 5th order polynomial coefficient to convert the DN to degree C 
        ; and this doesn't match the Steinhart-Hart model that well.
        ; Use the same method the Java code uses to convert prism_temp DN to degrees
        query="SELECT * FROM dbo.TimSimThermistorCalibration where thermistorItemName='prism_drive_temp'"
        query+=" and instrumentName='"+instrumentName+"' and version=1"
        query_database,/reset
        query_database, query, javaCal
        if size(javaCal,/tname) ne 'STRUCT' then begin
           ; no database access - force values
           javaCal={resistance:19600d, referenceVoltage:7.17d, $
                    a:0.0011292410d, b:0.00023410770, c:8.7754680d-08, zeroCelsius:273.15d}
        endif
        vout = -20d/2d^16 * (*tlm_data.(1)).housekeeping.dn
        rt=vout * javaCal.resistance / (javaCal.referenceVoltage - vout)
        tmp = 1d / (javaCal.a + javaCal.b * alog(rt) + javaCal.c * alog(rt)^3d) - javaCal.zerocelsius
        (*tlm_data.(1)).housekeeping.eu = tmp
    endif

    ; convert the detector temperature (miniature thermistors)
    if size((*tlm_data.(5)).housekeeping, /tname) eq 'STRUCT' then begin
        query="SELECT * FROM dbo.TimSimThermistorCalibration where thermistorItemName='vis_1_temp'"
        query+=" and instrumentName='"+instrumentName+"' and version=1"
        query_database,/reset
        query_database, query, javaCal
        if size(javaCal,/tname) ne 'STRUCT' then begin
           ; no database access - force values
           javaCal={resistance:23340d, referenceVoltage:2.0d, gain:7.98d, $
                    a:0.0011292410d, b:0.00023410770, c:8.7754680d-08, zeroCelsius:273.15d}
        endif
        vout = -20d/2d^16 * (*tlm_data.(5)).housekeeping.dn
        rt = javaCal.resistance * (javaCal.referenceVoltage * javaCal.gain + vout) 
        rt /= (javaCal.referenceVoltage * javaCal.gain - vout)
        tmp = 1d / (javaCal.a + javaCal.b * alog(rt) + javaCal.c * alog(rt)^3d) - javaCal.zerocelsius
        (*tlm_data.(5)).housekeeping.eu = tmp
    endif

    ; average the photo-diode array data 
    valid_pos = where((*tlm_data.(0)).science.dn ge pos0 and (*tlm_data.(0)).science.dn le pos1 and $
                      (*tlm_data.(0)).science.timetag ge t0-10d6 and (*tlm_data.(0)).science.timetag le t1+10d6 , ndata)
    if ndata eq 0 then begin
        print,'No valid prism position for the specified instrumentId and time range'
        return,-1
    endif

    outdata = replicate({timetag:0.0d, ccdpos:0.0d, ccdpos_uncorr:0.0d, wavelength:0.0d, dn:0.0d, $
                         dn_tempcorr:0.0d, irradiance:0.0d, prismTemp:-999d, detectorTemp:0.0d},ndata)

    ; average out the dns during integration period
    delta_t = median((*tlm_data.(0)).science[valid_pos[1:-1]].timetag - (*tlm_data.(0)).science[valid_pos[0:-2]].timetag)
    for pos=0L,ndata-1L do begin
        ; prism position is reported every 1.0 second
        tt0 = (*tlm_data.(0)).science[valid_pos[pos]].timetag
        if pos lt ndata-1L then tt1 = (*tlm_data.(0)).science[valid_pos[pos+1]].timetag else tt1=tt0+delta_t
        if tt1-tt0 gt delta_t*2 then tt1=tt0+delta_t ;new
        ; diode signal is reported every 0.1 second and we skip the first 4 samples (settling time)
        p=where((*tlm_data.(2)).science.timetag ge tt0 and (*tlm_data.(2)).science.timetag lt tt1 , count)
        if count ge 5 then begin
            ; skip the first 4 points after new ccdpos and average the remaining (6 points)
            outdata[pos].dn = mean((*tlm_data.(2)).science[p[4:-1]].dn)
            ; define the timetag as the middle of the sample
            outdata[pos].timetag = mean((*tlm_data.(2)).science[p[4:-1]].timetag)
        endif else begin
            outdata[pos].dn = mean((*tlm_data.(2)).science[p].dn)
            outdata[pos].timetag = mean((*tlm_data.(2)).science[p].timetag)
        endelse
    endfor

    ; look for dark (shutter closed) before and after our time range
    dark_dn=[]
    dark_timetag=[]
    if NOT keyword_set(noDark) then begin
        ; get the shutter dethinned to the detector cadence
        if keyword_set(verbose) then print,' calculating darks ...'
        shutter_pos = dethin_data((*tlm_data.(3)).science.dn, (*tlm_data.(3)).science.timetag, (*tlm_data.(2)).science.timetag)
        darkpos0 = where(shutter_pos eq SHT_CLOSE and (*tlm_data.(3)).science.timetag le t0, dcount0)
        darkpos1 = where(shutter_pos eq SHT_CLOSE and (*tlm_data.(3)).science.timetag ge t1, dcount1)
        ; use 30 seconds of dark exposure on both ends of the exposure
        if dcount0 gt 0  and dcount1 gt 0 then begin
            ; print,'Darks found at both ends of activity'
            tt1 = (*tlm_data.(2)).science[darkpos0[-1]].timetag
            tt0 = tt1 - 30.0d6 
            p0 = where((*tlm_data.(2)).science.timetag ge tt0 and (*tlm_data.(2)).science.timetag lt tt1 and $
                       shutter_pos eq SHT_CLOSE, np0)
            if np0 gt 0 then begin
                ; make sure the dark is ok because we found times when it looks like
                ; the shutter was open with high counts for darks 
                ; (even if telemetry indicates shutter closed)
                resistant_mean,(*tlm_data.(2)).science[p0].dn,2.0,mean_val,sigma
                if sigma lt 1.0 then begin
                    dark_dn=[dark_dn,double(median((*tlm_data.(2)).science[p0].dn))]
                endif else dark_dn=[dark_dn,0.0d]
                dark_timetag=[dark_timetag,mean((*tlm_data.(2)).science[p0].timetag)]
            endif
            tt0 = (*tlm_data.(2)).science[darkpos1[0]].timetag
            tt1 = tt0 + 30.0d6
            p1 = where((*tlm_data.(2)).science.timetag ge tt0 and (*tlm_data.(2)).science.timetag lt tt1 and $
                       shutter_pos eq SHT_CLOSE, np1)
            if np1 gt 0 then begin
                resistant_mean,(*tlm_data.(2)).science[p1].dn,2.0,mean_val,sigma
                if sigma lt 1.0 then begin
                    dark_dn=[dark_dn,double(median((*tlm_data.(2)).science[p1].dn))]
                endif else dark_dn=[dark_dn,0.0d]
                dark_timetag=[dark_timetag,mean((*tlm_data.(2)).science[p1].timetag)]
            endif
            p=where(dark_dn eq 0.0, comp=c,count)
            if count gt 0 then dark_dn[p[0]]=dark_dn[c[0]]
        endif else if dcount0 gt 0 then begin
            ; print,'Darks found after activity only'
            tt1 = (*tlm_data.(2)).science[darkpos0[-1]].timetag
            tt0 = tt1 - 30.0d6
            p0 = where((*tlm_data.(2)).science.timetag ge tt0 and (*tlm_data.(2)).science.timetag lt tt1 and $
                       shutter_pos eq SHT_CLOSE, np0)
            if np0 gt 0 then begin
                resistant_mean,(*tlm_data.(2)).science[p0].dn,2.0,mean_val,sigma
                if sigma lt 1.0 then begin
                    value = double(median((*tlm_data.(2)).science[p0].dn))
                endif else value=0.0d
                dark_dn=[value,value]
                dark_timetag=[mean((*tlm_data.(2)).science[p0].timetag),t1]
            endif
        endif else if dcount1 gt 0 then begin
            ; print,'Darks found before activity only'
            tt0 = (*tlm_data.(2)).science[darkpos1[0]].timetag
            tt1 = tt0 + 30.0d6
            p1 = where((*tlm_data.(2)).science.timetag ge tt0 and (*tlm_data.(2)).science.timetag lt tt1 and $
                       shutter_pos eq SHT_CLOSE, np1)
            if np1 gt 0 then begin
                resistant_mean,(*tlm_data.(2)).science[p1].dn,2.0,mean_val,sigma
                if sigma lt 1.0 then begin
                    value = double(median((*tlm_data.(2)).science[p1].dn))
                endif else value=0.0d
                dark_dn=[value, value]
                dark_timetag=[t0,dark_timetag,mean((*tlm_data.(2)).science[p1].timetag)]
            endif
        endif else begin
            ; no dark values found: set to 0.0
            print,'No DARK found -> using 0.0'
            dark_dn=[0.0d,0.0d]
            dark_timetag=[t0,t1]
        endelse
        if keyword_set(verbose) then print,' subtracting darks ...'
        darks = interpol(dark_dn,dark_timetag,outdata.timetag)
        outdata.dn -= darks
    endif

    if keyword_set(verbose) then print,"reading ccdCorr ..."
    run_t2=systime(/sec)

    outdata.dn_tempcorr = outdata.dn
    outdata.ccdpos_uncorr = (*tlm_data.(0)).science[valid_pos].dn
    outdata.ccdpos = outdata.ccdpos_uncorr
    if n_elements(fitCoeff) gt 0 and NOT keyword_set(noCcdCorr) then begin
        ; coefficients were provided -> calculate new CCD positions
        if keyword_set(verbose) then print,' applying provided CCD shift coefficients ...'
        offset = poly(outdata.ccdpos, fitCoeff)
        ;outdata.ccdpos += offset
        outdata.ccdpos -= offset   ; with the ccdshift from SimCcdShiftCal version=2
    endif else if NOT keyword_set(noCcdCorr) then begin
        ; determine the intg_period
        ; we EXPECT it to be constant for the duration of the specified time range
        ; only care about the intg_period before day 450
        mn=min((*tlm_data.(4)).housekeeping.dn,max=mx)
        if gps2sd(t0/1d6) lt 450.0 AND mn ne mx then begin
            print,'Error: intg_period is NOT constant throughout the specified time range'
            return,-1
        endif
        intg_period = (*tlm_data.(4)).housekeeping[0].dn

        if keyword_set(verbose) then print,' calculating the CCD shifts  ...'

        ; grab the highest version effective on the date of this spectra
        query1 = "SELECT calibrationSetId,version FROM CalibrationMetadata WHERE"
        query1 = query1+" calibrationTableName='SimCcdShiftCal' AND instrumentModeId="+strtrim(string(instrumentModeId),2)
        query1 = query1+" AND effectiveDate <= '"+strcompress(jd2syb(gps2jd(t0/1d6)))+"' ORDER BY version desc"
        ; for some unknown IDL reason, the string doesn't work until we "print" it
        query1=string(query1,format='(a)')
        query_database, query1, res, nrows
        if nrows eq 0 then begin
            print,'Error: no CalibrationMetadata found for SimCcdShiftCal for the requested instruemntModeId'
            return,-1
        endif
        calibrationSetId=res[0].calibrationSetId
        if keyword_set(verbose) then $
            print,'  using SimCcdShiftCal iMode='+strtrim(string(instrumentModeId),2)+$
            ' version='+strtrim(string(res[0].version),2)+', calibrationSetId='+strtrim(string(calibrationSetId),2)

        query2 ="SELECT * FROM SimCcdShiftCal WHERE calibrationSetId="+strtrim(string(calibrationSetId),2)
        query2+=" and effectiveDateGps <= "+string(LONG64(t1))
        query2+=" ORDER BY effectiveDateGps asc"
        query_database, query2, res, nrows
        if nrows eq 0 then begin
            print,'Error: no SimCcdShiftCal found for modeId / intg_period = ',instrumentModeId, intg_period
            return,-1
        endif

        if nrows eq 0 then begin
            print,'Error: No matching entry in the SimCcdShift table'
            return,-1
        endif
        ; use the latest valid entry
        p=-1
        ; calculate the CCD spot position correction
        ccdCorr = poly(outdata.ccdpos, [res[p].c0, res[p].c1,res[p].c2,res[p].c3])
        ; with the version 20 of the ccdshift, the correction should be subtracted
        ;outdata.ccdpos += ccdCorr
        outdata.ccdpos -= ccdCorr

        if keyword_set(verbose) then $
            print,'  using SimCcdShiftCal: intg_period='+strtrim(string(intg_period),2)+$
            ', coeffs=',res[p].c0, res[p].c1,res[p].c2,res[p].c3
    endif

    if keyword_set(verbose) then print,string(run_t2-run_t1,format='("ccdCorr ",F0.3," seconds")')

    if keyword_set(verbose) then print,"reading profileIntegralCal ..."
    run_t3=systime(/sec)


    ; add the prism temperature to structure
    if size((*tlm_data.(1)).housekeeping, /tname) eq 'STRUCT' then begin
        prismTemp = interpol((*tlm_data.(1)).housekeeping.eu, (*tlm_data.(1)).housekeeping.timetag, outdata.timetag, /spline)
        outdata.prismTemp = prismTemp
    endif

    ; get the profileIntegralCal data for this instrumentMode
    if size(profile_data,/tname) ne 'STRUCT' then begin
        query1 = "SELECT calibrationSetId,version FROM CalibrationMetadata WHERE "
        query1 = query1+"calibrationTableName='SimProfileIntegralCal' AND instrumentModeId="+strtrim(string(instrumentModeId),2)
        if n_elements(profileDate) ne 0 then effectiveDate=profileDate else effectiveDate=t0
        query1 = query1+" AND effectiveDate <= '"+jd2syb(gps2jd(effectiveDate/1d6))+"'"
        if keyword_set(version) then query1=query1+" AND version<="+strtrim(string(version),2)
        query1 = query1+" ORDER BY version DESC"
        if keyword_set(verbose) then print,'query="'+query1+'"'
        query_database, query1[0], res, nrows, database='SORCE'
        if nrows eq 0 then begin
            print,'Error: no CalibrationMetadata found for the requested date/version'
            return,-1
        endif
        calibrationSetId=res[0].calibrationSetId
        if keyword_set(verbose) then $
            print,'  using SimProfileIntegralCal iMode='+strtrim(string(instrumentModeId),2)+$
            ' version='+strtrim(string(res[0].version),2)+', calibrationSetId='+strtrim(string(calibrationSetId),2)

        ; now get the SimProfileIntegralCal entry for this calibrationSetId
        query2="SELECT * FROM SimProfileIntegralCal WHERE calibrationSetId="+strtrim(string(calibrationSetId),2)
        query_database, query2, profile_data, nrows, database='SORCE'

        if nrows eq 0 then begin
            if keyword_set(verbose) then $
                print,'NO SimProfileIntegralCal entries found for the calibrationSetId='+strtrim(string(calibrationSetId),2)
            return,outdata
        endif
    endif

    if keyword_set(verbose) then print,string(run_t3-run_t2,format='("profileIntegralCal ",F0.3," seconds")')
    run_t4=systime(/sec)

    spotRef = profile_data.x
    ; the reference temperature is really the same for every row inSimProfileIntegralCal
    prismRefTemp = mean(profile_data.y0)
    wavelengthRef = profile_data.y4
    wavelengthRefTCoef = profile_data.y18

    ; interpolate the profileIntegralCal to get the wavelength for each ccdpos
    outdata.wavelength = interpol(wavelengthRef, spotRef, outdata.ccdpos,/spline)
    if (*tlm_data.(5)).n_housekeeping gt 0 then $
        outdata.detectorTemp = interpol((*tlm_data.(5)).housekeeping.eu, (*tlm_data.(5)).housekeeping.timetag, outdata.timetag, /spline)

    if keyword_set(verbose) then print,"reading 1AU ..."

    if NOT keyword_set(no1au) then begin
        if keyword_set(verbose) then $
            print,'  using SolarDistAndDopplerFixedStep2 iMode=50, version=7'
        ; correct for 1AU and Doppler
        ; the instrumentId for SORCE is 50 and we're at version=6  (querying for the max(version) is VERY SLOW)
        if size(solarDist,/tname) ne 'STRUCT' then begin
            query2="SELECT * FROM SolarDistAndDopplerFixedStep2 WHERE instrumentModeId=50 and version=7 and "
            query2 = query2+"microsecondsSinceGpsEpoch>="+strtrim(string(ulong64(t0-0.5d6*86400d)),2)+" and "
            query2 = query2+"microsecondsSinceGpsEpoch<="+strtrim(string(ulong64(t1+0.5d6*86400d)),2)
            query_database, query2, solarDist, nrows, database='SORCE'
        endif
        sunObserverDistance = interpol(solarDist.sunObserverDistanceCorrection, $
            solarDist.microsecondsSinceGpsEpoch, outdata.timetag,/spline)
        outdata.dn_tempcorr /= sunObserverDistance
    endif

    if keyword_set(verbose) then print,string(run_t4-run_t3,format='("1AU ",F0.3," seconds")')

    if NOT keyword_set(noTempCorr) then begin
        ; perform the wavelength correction due to prism temperature
        if (*tlm_data.(1)).n_housekeeping eq 0 then begin
            print,'Error: no prism_drive_temp found in telemetry'
            return,-1
        endif

        if keyword_set(verbose) then print,'  applying prismTempCorr ...'

        ; test the instrument model to get the shift of wavelength with temp instead of Y18
        ;outdata.wavelength = ccd2lambda(instrumentModeId, outdata.ccdpos, outdata.prismTemp, tsis=tsis)

        dndt = interpol(wavelengthRefTCoef, spotRef, outdata.ccdpos, /spline)
        p=where(outdata.prismTemp gt -999d,npts)
        if npts gt 0 then begin
            outdata.wavelength = outdata.wavelength + (outdata.prismTemp - prismRefTemp) * dndt
        endif else begin
            print,'No PrismTemp found for this scan - SKIPPING TEMPCORR'
        endelse

    endif

    if keyword_set(verbose) then print,"reading diode tempCorr ..."

    ; apply the new diode temperature correction to DN
    ; get the data (not in the DB yet)

    if size(dettempcorr,/tname) ne 'STRUCT' then begin
        case instrumentModeId of
            41: begin
                ;infile = file_search("$SORCE_DATA","sima_vis_tempcorr_3.txt",/expand_env,/full)
                ;infile = file_search("$SORCE_DATA","sima_vis_tempcorr_20.txt",/expand_env,/full)
                infile = file_search("$SORCE_DATA","sima_vis_tempcorr_6.txt",/expand_env,/full)
                readcol,infile[0],wv,dndt,detreftemp,format='(D,D,D)',/silent
                tempcorr = interpol(dndt,wv, outdata.wavelength, /lsq)
                outdata.dn_tempcorr = outdata.dn_tempcorr * (1d + (detreftemp[0] - outdata.detectorTemp) * tempcorr)
            end
            43: begin
                ;infile = file_search("$SORCE_DATA","sima_uv_tempcorr_3.txt",/expand_env,/full)
                ;infile = file_search("$SORCE_DATA","sima_vis_tempcorr_20.txt",/expand_env,/full)
                infile = file_search("$SORCE_DATA","sima_vis_tempcorr_6.txt",/expand_env,/full)
                readcol,infile[0],wv,dndt,detreftemp,format='(D,D,D)',/silent
                tempcorr = interpol(dndt,wv, outdata.wavelength, /lsq)
                outdata.dn_tempcorr = outdata.dn_tempcorr * (1d + (detreftemp[0] - outdata.detectorTemp) * tempcorr)
            end
            44: begin
                ;infile = file_search("$SORCE_DATA","sima_ir_tempcorr_3.txt",/expand_env,/full)
                ;infile = file_search("$SORCE_DATA","sima_ir_tempcorr_21.txt",/expand_env,/full)
                infile = file_search("$SORCE_DATA","sima_ir_tempcorr_24.txt",/expand_env,/full)
                ;infile = file_search("$SORCE_DATA","sima_ir_tempcorr_v24_2.txt",/expand_env,/full)
                ;infile = file_search("$SORCE_DATA","sima_ir_tempcorr_v24_7.txt",/expand_env,/full)
                ;infile = file_search("$SORCE_DATA","sima_ir_tempcorr_25.txt",/expand_env,/full)
                readcol,infile[0],wv,dndt,detreftemp,format='(D,D,D)',/silent
                tempcorr = interpol(dndt,wv, outdata.wavelength, /lsq)

                ;infile = file_search("$SORCE_DATA","sima_ir_tempcorr_v24_3.txt",/expand_env,/full)
                ;readcol,infile[0],wv,detreftemp,c0,c1,c2,c3,format='(D,D,D,D,D,D)',/silent
                ;; we need to loop for every wavelength since the diode temperature changes during the scan
                ;tempcorr=dblarr(n_elements(outdata))
                ;for i=0L,n_elements(outdata.wavelength)-1 do begin
                ;    dt=(detreftemp[0] - outdata[i].detectorTemp)
                ;    p0=(where(wv le outdata[i].wavelength))[-1]
                ;    p1=(where(wv gt outdata[i].wavelength))[0]
                ;    temp0=poly(dt,[c0[p0],c1[p0],c2[p0],c3[p0]])
                ;    temp1=poly(dt,[c0[p1],c1[p1],c2[p1],c3[p1]])
                ;    tempcorr[i]=interpol([temp0,temp1],[wv[p0],wv[p1]], outdata[i].wavelength)
                ;endfor


                ; the file below contains a 5th order fit to the alog(-deltadndt)
                ;infile = file_search("$SORCE_DATA","sima_ir_tempcorr_v24_6.txt",/expand_env,/full)
                ;readcol,infile[0],wv,detreftemp,c0,c1,c2,c3,c4,c5,format='(D,D,D,D,D,D,D,D)',/silent

                ; we need to loop for every wavelength since the diode temperature changes during the scan
                ;tempcorr=dblarr(n_elements(outdata))
                ;for i=0L,n_elements(outdata.wavelength)-1 do begin
                ;    dt=(detreftemp[0] - outdata[i].detectorTemp)
                ;    p0=(where(wv le outdata[i].wavelength))[-1]
                ;    p1=(where(wv gt outdata[i].wavelength))[0]
                ;    temp0=poly(dt,[c0[p0],c1[p0],c2[p0],c3[p0],c4[p0],c5[p0]])
                ;    temp1=poly(dt,[c0[p1],c1[p1],c2[p1],c3[p1],c4[p1],c5[p1]])
                ;    temp=interpol([temp0,temp1],[wv[p0],wv[p1]], outdata[i].wavelength)
                ;    tempcorr[i] = -1d * exp(temp)
                ;endfor

                outdata.dn_tempcorr = outdata.dn_tempcorr * (1d + (detreftemp[0] - outdata.detectorTemp) * tempcorr)
            end
            45: begin
                ;infile = file_search("$SORCE_DATA","sima_vis_tempcorr_3.txt",/expand_env,/full)
                ;infile = file_search("$SORCE_DATA","sima_vis_tempcorr_20.txt",/expand_env,/full)
                infile = file_search("$SORCE_DATA","sima_vis_tempcorr_6.txt",/expand_env,/full)
                readcol,infile[0],wv,dndt,detreftemp,format='(D,D,D)',/silent
                tempcorr = interpol(dndt,wv, outdata.wavelength, /lsq)
                outdata.dn_tempcorr = outdata.dn_tempcorr * (1d + (detreftemp[0] - outdata.detectorTemp) * tempcorr)
            end
            47: begin
                ;infile = file_search("$SORCE_DATA","sima_uv_tempcorr_3.txt",/expand_env,/full)
                ;infile = file_search("$SORCE_DATA","sima_vis_tempcorr_20.txt",/expand_env,/full)
                infile = file_search("$SORCE_DATA","sima_vis_tempcorr_6.txt",/expand_env,/full)
                readcol,infile[0],wv,dndt,detreftemp,format='(D,D,D)',/silent
                tempcorr = interpol(dndt,wv, outdata.wavelength, /lsq)
                outdata.dn_tempcorr = outdata.dn_tempcorr * (1d + (detreftemp[0] - outdata.detectorTemp) * tempcorr)
            end
            48: begin
                ;infile = file_search("$SORCE_DATA","sima_ir_tempcorr_3.txt",/expand_env,/full)
                infile = file_search("$SORCE_DATA","sima_ir_tempcorr_21.txt",/expand_env,/full)
                readcol,infile[0],wv,dndt,detreftemp,format='(D,D,D)',/silent
                tempcorr = interpol(dndt,wv, outdata.wavelength, /lsq)
                outdata.dn_tempcorr = outdata.dn_tempcorr * (1d + (detreftemp[0] - outdata.detectorTemp) * tempcorr)
            end
        endcase
        dettempcorr = {infile:infile, wv:wv, dndt:dndt, detreftemp:detreftemp}
    endif else begin
        tempcorr = interpol(dettempcorr.dndt, dettempcorr.wv, outdata.wavelength, /lsq)
        outdata.dn_tempcorr = outdata.dn_tempcorr * (1d + (dettempcorr.detreftemp[0] - outdata.detectorTemp) * tempcorr)
    endelse
   
    if keyword_set(verbose) then print,'  applying diodeTempCorr from '+dettempcorr.infile+' ...'
    run_t5=systime(/sec)
    if keyword_set(verbose) then print,string(run_t5-run_t4,format='("diodeTempCorr ",F0.3," seconds")')

    if keyword_set(dn_only) then begin
        if keyword_set(verbose) then print,string(run_t5-run_t0,format='("Total time ",F0.3," seconds")')
        return, outdata
    endif

    ; now calculate the irradiance
    ; dnData = (AveragedDN - DarkDN) / ADCGain /feedBackResistance
    ; UncorrectedIrradiance = dnData / sensitivityIntegral / sunObserverDistance / Doppler 
    ; correctedIrradiance = unCorrectedIrradiance / DetectorDeg / PrismDeg 
    ; calibratedIrradiance = correctedIrradiance * alephFactorValue
    ; 
    ; for now we skip the Doppler, DetectorDeg and PrismDeg
    ;
    ; get the calibrationSetId for the ADCGain 
    query1 = "SELECT calibrationSetId,version FROM CalibrationMetadata WHERE "
    query1 = query1+"calibrationTableName='SimDiodeAdcGainCalibrationData' AND instrumentModeId="+strtrim(string(instrumentModeId),2)
    query1 = query1+" AND effectiveDate <= '"+jd2syb(gps2jd(t0/1d6))+"'"
    if keyword_set(version) then query1=query1+" AND version<="+strtrim(string(version),2)
    query1 = query1+" ORDER BY version DESC"
    query_database, query1[0], res, nrows, database='SORCE'
    if nrows eq 0 then begin
        print,'Error: no CalibrationMetadata found for the requested date/version'
        return,-1
    endif
    calibrationSetId=res[0].calibrationSetId
    if keyword_set(verbose) then $
        print,'  using SimDiodeAdcGainCalibrationData iMode='+strtrim(string(instrumentModeId),2)+$
        ' version='+strtrim(string(res[0].version),2)+', calibrationSetId='+strtrim(string(calibrationSetId),2)

    query2="SELECT * FROM SimDiodeAdcGainCalibrationData WHERE calibrationSetId="+strtrim(string(calibrationSetId),2)
    query_database, query2, adcgain_data, nrows, database='SORCE'
    adcgain = adcgain_data.diodeadcgain

    query1 = "SELECT calibrationSetId,version FROM CalibrationMetadata WHERE "
    query1 = query1+"calibrationTableName='SimDiodeFeedbackResistanceCal' AND instrumentModeId="+strtrim(string(instrumentModeId),2)
    query1 = query1+" AND effectiveDate <= '"+jd2syb(gps2jd(t0/1d6))+"'"
    if keyword_set(version) then query1=query1+" AND version<="+strtrim(string(version),2)
    query1 = query1+" ORDER BY version DESC"
    query_database, query1[0], res, nrows, database='SORCE'
    if nrows eq 0 then begin
        print,'Error: no CalibrationMetadata found for the requested date/version'
        return,-1
    endif
    calibrationSetId=res[0].calibrationSetId
    if keyword_set(verbose) then $
        print,'  using SimDiodeFeedbackResistanceCal iMode='+strtrim(string(instrumentModeId),2)+$
        ' version='+strtrim(string(res[0].version),2)+', calibrationSetId='+strtrim(string(calibrationSetId),2)

    query2="SELECT * FROM SimDiodeFeedbackResistanceCal WHERE calibrationSetId="+strtrim(string(calibrationSetId),2)
    query_database, query2, feedbackR_data, nrows, database='SORCE'
    feedbackResistance = feedbackR_data.resistanceInOhms
    outdata.irradiance = outdata.dn_tempcorr / adcGain / feedbackResistance

    ; get the sensitivityIntegral = y8 + y10*(prismTemp - y0) + y11*(degColumn - y1) + y12*(detTemp - y2)
    y0 = prismRefTemp                                                    ; reference prism temperature
    y1 = interpol(profile_data.y1, wavelengthRef, outdata.wavelength,/spline)      ; profile integral
    y2 = mean(profile_data.y2)                                           ; reference detector temperature
    y8 = interpol(profile_data.y8, wavelengthRef, outdata.wavelength,/spline)      ; profile integral
    y10 = interpol(profile_data.y10, wavelengthRef, outdata.wavelength,/spline)    ; profIntegralPrismTCoef
    y11 = interpol(profile_data.y11, wavelengthRef, outdata.wavelength,/spline)    ; profIntegralDegrCoef 
    y12 = interpol(profile_data.y12, wavelengthRef, outdata.wavelength,/spline)    ; profIntegralDetTCoef

    ; the diode temperature correction is now applied to the dark subtracted dn
    ; and the y8 has now been interpolated in wavelengthh space - no need to add a correction for prismTemp twice
    sensitivityProfile = y8 ;+y11*(degColumn-y1); + y10*(outdata.prismTemp - y0) + y12*(outdata.detectorTemp-y2)
    outdata.irradiance /= sensitivityProfile

    ; if uncorr was specified do not apply the prism degradation correction
    if keyword_set(uncorr) then return, outdata


    ; starting with version 21, we changed the degColumn to a 2D polynomial fit
    ; prismDeg = (1-a)*exp(-Kappa*degCol) + a*exp(-Kappa*degCol/2)
    ; correctedIrrad = uncorrectedIrrad / prismDeg

    ; get the rayPath for this instrunemnetModeId
    query1 = "SELECT calibrationSetId,version FROM CalibrationMetadata WHERE "
    query1 = query1+"calibrationTableName='SimRayPathDegradationParams' AND instrumentModeId="+strtrim(string(instrumentModeId),2)
    query1 = query1+" ORDER BY version desc"
    query_database, query1, res, nrows, database='SORCE'
    if nrows eq 0 then begin
        print,'Error: no CalibrationMetadata found for the requested date/version'
        return,outdata
    endif
    calibrationSetId=res[0].calibrationSetId
    query1 = "SELECT * from SimRayPathDegradationParams where calibrationSetId="+strtrim(string(calibrationSetId),2)
    query_database, query1, rayPathCal, nrows, database='SORCE'

    ; get the Kappa for this instrunemnetModeId
    query1 = "SELECT calibrationSetId,version FROM CalibrationMetadata WHERE "
    query1 = query1+"calibrationTableName='SimPrismDegKappaCal' AND instrumentModeId="+strtrim(string(instrumentModeId),2)
    query1 = query1+" ORDER BY version desc"
    query_database, query1, res, nrows, database='SORCE'
    if nrows eq 0 then begin
        print,'Error: no CalibrationMetadata found for the requested date/version'
        return,outdata
    endif
    calibrationSetId=res[0].calibrationSetId
    query1 = "SELECT * from SimPrismDegKappaCal where calibrationSetId="+strtrim(string(calibrationSetId),2)
    query_database, query1, kappaCal, nrows, database='SORCE'

    ; get the degradationColumn for this instrunemnetModeId
    query1 = "SELECT calibrationSetId,version FROM CalibrationMetadata WHERE "
    query1 = query1+"calibrationTableName='SimPrismDegColCalWave2D' AND instrumentModeId="+strtrim(string(instrumentModeId),2)
    query1 = query1+" ORDER BY version desc"
    query_database, query1, res, nrows, database='SORCE'
    if nrows eq 0 then begin
        print,'Error: no CalibrationMetadata found for the requested date/version'
        return,outdata
    endif
    calibrationSetId=res[0].calibrationSetId
    query1 = "SELECT * from SimPrismDegColCalWave2D where calibrationSetId="+strtrim(string(calibrationSetId),2)+" and "
    query2 = "orbitEndGpsMicroseconds>="+strtrim(string(ulong64(t0)),2)+" and "
    query2 = query2+"orbitEndGpsMicroseconds<="+strtrim(string(ulong64(t1)),2)
    query_database, query1+query2, degColumnCal, nrows, database='SORCE'
    if nrows eq 0 then begin
        ; get the scan before t0 and use that
        qq="SELECT top 1 * FROM dbo.SimPrismDegColCalWave2D where calibrationSetId="+strtrim(string(calibrationSetId),2)+" and "
        qq=qq+"orbitEndGpsMicroseconds-"+string(LONG64(t0))+"<0 order by orbitEndGpsMicroseconds desc"
        query_database, qq, temp, nrows, database='SORCE'
        query2 = "orbitEndGpsMicroseconds="+strtrim(string(ulong64(temp.ORBITENDGPSMICROSECONDS)),2)
        query_database, query1+query2, degColumnCal, nrows, database='SORCE'
    endif

    kappa = interpol(kappaCal.Y, 10^kappaCal.x, outdata.wavelength, /lsq)
    degCol = interpol(degColumnCal.DEGRADATIONCOLUMN, degColumnCal.wavelength, outdata.wavelength, /lsq)
    raypath=interpol(rayPathCal.SINGLEPASSAREAFRACTION, rayPathCal.wavelength, outdata.wavelength, /lsq)
    if ~ keyword_set(no1au) then raypath *= sunObserverDistance ^ 2d

    prismDeg = (1d - raypath) * exp(-1d * kappa * degCol) + raypath * exp(-1d * kappa * degCol / 2d)
    outdata.irradiance /= prismDeg


; TODO: add the aleph for the calibrated irradiance

    return, outdata

end
