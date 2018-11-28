;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Align spectra with a reference spectra by changing the ccd subpixel
;   position and minimizing the residuals (optionally of the first or second
;   derivatives). We're comparing the irradiance for each spectra.
;
; CALLING SEQUENCE:
;   result = find_ccdoffset(t0, t1, instrumentModeId, refSpect=refSpect, /missionDays)
;
; INPUT PARAMETERS:
;   startTime -
;      The start time of a spectral scan.
;      May be specified as julian days, gps microseconds, or
;      SORCE mission day numbers.
;
;   stopTime -
;      The stop time of a spectral scan.
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
;  You may want to run get_sorce_plan before hand to make sure an activity
;  with observations in the desired mode was performed during the time
;  span of interest.
;  It is assumed that the CCD offset is constant for the whole scan.
;  This is NOT the case for data taken with a saturated CCD spot image
;  (pre mission 450).
;
; REVISION HISTORY:
;   Revision: $Id: find_ccdoffset.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;-
;
function find_ccdoffset, startTime, stopTime, instrumentModeId, $
         gps=gps, missionDays=missionDays, julianDays=julianDays, $
         refT0=refT0, refT1=refT1, noTempCorr=noTempCorr, noDark=noDark, $
         version=version, noCcdCorr=noCcdCorr, no1AU=no1AU, verbose=verbose, $
         deriv_order=deriv_order, noplot=noplot, refSpect=refSpect

    nrows=0
    if n_elements(deriv_order) eq 0 then deriv_order=2
    ; validate the instrumentModeId
    if n_elements(instrumentModeId) eq 0 then begin
        doc_library,'find_ccdoffset'
        print,''
        print,'Missing instrumentModeId '
        return,-1
    endif else begin
        ; verify we got the right instrument
        valid_inst = [41,42,43,44,45,46,47,48]
        instrumentModeId = fix(instrumentModeId)
        if (where(valid_inst eq instrumentModeId))[0] lt 0 then  begin
            doc_library,'find_ccdoffset'
            print,''
            print,'Invalid instrumentModeId was provided'
            return,-1
        endif
    endelse

    SHT_CLOSE=1
    SHT_OPEN=0

    ; define the telemetry items according to the specified mode
    case instrumentModeId of
        31: begin
            ; SIMA ESR
            instrument='SIMA'
            detector='esr'
            channel='a'
            spot='nominal'
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
            detector='vis'
            channel='a'
            spot='nominal'
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
            detector='vis'
            channel='a'
            spot='nominal'
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
            detector='uv'
            channel='a'
            spot='uv'
            prismId = 101862L
            prismTempId = 101909L
            detectorId = 100046L
            shutterId = 102121L
            pos0=6590L
            pos1=49890L
            intg_period=101395L
            detectorTemp=101069L
        end
        44: begin
            ; SIMA IR
            instrument='SIMA'
            detector='ir'
            channel='a'
            spot='nominal'
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
            detector='esr'
            channel='b'
            spot='nominal'
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
            detector='vis'
            channel='b'
            spot='nominal'
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
            detector='vis'
            channel='b'
            spot='nominal'
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
            detector='uv'
            channel='b'
            spot='uv'
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
            detector='ir'
            channel='b'
            spot='nominal'
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

    if keyword_set(missionDays) then begin
       ; user specified time in mission (sorce) days
       t0 = sd2gps(startTime)*1.d6
       t1 = sd2gps(stopTime)*1.d6
       if n_elements(refT0) ne 0 then rT0=sd2gps(refT0)*1d6
       if n_elements(refT1) ne 0 then rT1=sd2gps(refT1)*1d6
    endif else if keyword_set(julianDays) then begin
       ; user specified time in julian days
       t0 = jd2gps(startTime)*1.d6
       t1 = jd2gps(stopTime)*1.d6
       if n_elements(refT0) ne 0 then rT0=jd2gps(refT0)*1d6
       if n_elements(refT1) ne 0 then rT1=jd2gps(refT1)*1d6
    endif else begin
       ; user specified timetags in gps microseconds
       t0 = startTime
       t1 = stopTime
    endelse

    if n_elements(refT0) eq 0 or n_elements(refT1) eq 0 then begin
        rT0=sd2gps(453.67840d)*1d6
        rT1=sd2gps(453.69524d)*1d6
    endif

    ; get the reference spectra
    if keyword_set(verbose) then print,'Extracting the reference spectra ...'
    if size(refSpect,/tname) ne 'STRUCT' then $
        refSpect = get_sim_spectra(rT0, rT1, instrumentModeId, /gps, version=version, /noccdCorr, $
            noTempCorr=noTempCorr, nodark=nodark)

    if instrumentModeId eq 43 or instrumentModeId eq 47 then begin 
        wpos = where(refspect.wavelength ge 260.0)
    endif else if instrumentModeId eq 41 or instrumentModeId eq 45 then begin
        wpos = where(refspect.wavelength le 550.0)
    endif else if instrumentModeId eq 44 or instrumentModeId eq 48 then begin
        wpos = where(refspect.wavelength le 1630.0)
    endif

    ; now get the requested spectra from the raw DN 


    ; get the new diode temperature correction to DN (not in the DB yet)
    case instrumentModeId of
        41: begin
            infile = file_search("$SORCE_DATA","sima_vis_tempcorr_3.txt",/expand_env,/full)
            readcol,infile[0],wv,c0,c1,c2,c3,format='(D,D,D,D,D)',/silent
        end
        43: begin
            infile = file_search("$SORCE_DATA","sima_uv_tempcorr_3.txt",/expand_env,/full)
            readcol,infile[0],wv,c0,c1,c2,c3,format='(D,D,D,D,D)',/silent
        end
        44: begin
            infile = file_search("$SORCE_DATA","sima_ir_tempcorr_3.txt",/expand_env,/full)
            readcol,infile[0],wv,c0,c1,c2,c3,format='(D,D,D,D,D)',/silent
        end
        45: begin
            infile = file_search("$SORCE_DATA","sima_vis_tempcorr_3.txt",/expand_env,/full)
            readcol,infile[0],wv,c0,c1,c2,c3,format='(D,D,D,D,D)',/silent
        end
        47: begin
            infile = file_search("$SORCE_DATA","sima_uv_tempcorr_3.txt",/expand_env,/full)
            readcol,infile[0],wv,c0,c1,c2,c3,format='(D,D,D,D,D)',/silent
        end
        48: begin
            infile = file_search("$SORCE_DATA","sima_ir_tempcorr_3.txt",/expand_env,/full)
            readcol,infile[0],wv,c0,c1,c2,c3,format='(D,D,D,D,D)',/silent
        end
    endcase



    ; prepare the database connection 
    jstmt = fjava_get_jdbc_statement(user=name, password=password, dburl=dbURL, $
        dbdriver=dbdriver, server=server, database="SORCE")
    oJavaDbExchange = OBJ_NEW("oJava_DbExchange", jstmt)

    ; get the list of prism positions for the specified time range (+/- 2 minutes to get darks)
    if keyword_set(verbose) then print,'Extracting the spectra to compare ...'
    get_sorce_telemetry, tlm_data, info, t0-120.0d6, t1+120.0d6, /gps, $
        tlmId=[prismId,prismTempId,detectorId,shutterId,intg_period,detectorTemp]
    if size(tlm_data,/tname) ne 'STRUCT' then begin
        print,'No data found for the specified instrumentId and time range'
        return,-1
    endif else if (*tlm_data.(0)).n_science eq 0 then begin
        print,'No data found for the specified instrumentId and time range'
        return,-1
    endif

    ; average the photo-diode array data 
    valid_pos = where((*tlm_data.(0)).science.dn ge pos0 and (*tlm_data.(0)).science.dn le pos1 and $
                      (*tlm_data.(0)).science.timetag ge t0 and (*tlm_data.(0)).science.timetag le t1, ndata)
    if ndata eq 0 then begin
        print,'No valid prism position for the specified instrumentId and time range'
        return,-1
    endif

    outdata = replicate({timetag:0.0d, ccdpos:0.0d, ccdpos_uncorr:0.0d, wavelength:0.0d, dn:0.0d, $
                         dn_tempcorr:0.0d, irradiance:0.0d, prismTemp:0.0d, detectorTemp:0.0d},ndata)

    ; average out the dns during integration period
    if keyword_set(verbose) then print,'Average out the dns during integration period ...'
    delta_t = max((*tlm_data.(0)).science[valid_pos[1:-1]].timetag - (*tlm_data.(0)).science[valid_pos[0:-2]].timetag)
    for pos=0L,ndata-1L do begin
        tt0 = (*tlm_data.(0)).science[valid_pos[pos]].timetag
        if pos lt ndata-1L then tt1 = (*tlm_data.(0)).science[valid_pos[pos+1]].timetag else tt1=tt0+delta_t
        p=where((*tlm_data.(2)).science.timetag ge tt0 and (*tlm_data.(2)).science.timetag lt tt1 , count)
        if count lt 5 then continue
        ; skip the first 4 points after new ccdpos and average the remaining (6 points)
        outdata[pos].dn = mean((*tlm_data.(2)).science[p[4:-1]].dn)
        ; define the timetag as the middle of the sample
        outdata[pos].timetag = mean((*tlm_data.(2)).science[p[4:-1]].timetag)
    endfor

    ; look for dark (shutter closed) before and after our time range
    dark_dn=[]
    dark_timetag=[]
    if NOT keyword_set(noDark) then begin
        ; get the shutter dethinned to the detector cadence
        if keyword_set(verbose) then print,'Subtracting the darks ...'
        shutter_pos = dethin_data((*tlm_data.(3)).science.dn, (*tlm_data.(3)).science.timetag, (*tlm_data.(2)).science.timetag)
        darkpos0 = where(shutter_pos eq SHT_CLOSE and (*tlm_data.(2)).science.timetag le t0, dcount0)
        darkpos1 = where(shutter_pos eq SHT_CLOSE and (*tlm_data.(2)).science.timetag ge t1, dcount1)
        ; use 30 seconds of dark exposure on both ends of the exposure
        if dcount0 gt 0  and dcount1 gt 0 then begin
            ; print,'Darks found at both ends of activity'
            tt1 = (*tlm_data.(2)).science[darkpos0[-1]].timetag
            tt0 = tt1 - 30.0d6 
            p0 = where((*tlm_data.(2)).science.timetag ge tt0 and (*tlm_data.(2)).science.timetag lt tt1 and $
                       shutter_pos eq SHT_CLOSE, np0)
            if np0 gt 1 then begin
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
            if np0 gt 1 then begin
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
            if np1 gt 1 then begin
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
        darks = interpol(dark_dn,dark_timetag,outdata.timetag)
        outdata.dn -= darks
    endif

    outdata.ccdpos_uncorr = (*tlm_data.(0)).science[valid_pos].dn
    outdata.ccdpos = outdata.ccdpos_uncorr

    if (*tlm_data.(1)).n_housekeeping eq 0 then begin
        print,'Error: no prismTemp found for the requested date/version'
        return,-1
    endif
    if (*tlm_data.(5)).n_housekeeping eq 0 then begin
        print,'Error: no detectorTemp found for the requested date/version'
        return,-1
    endif

    outdata.prismTemp = interpol((*tlm_data.(1)).housekeeping.eu, (*tlm_data.(1)).housekeeping.timetag, outdata.timetag, /quad)
    outdata.detectorTemp = interpol((*tlm_data.(5)).housekeeping.eu, (*tlm_data.(5)).housekeeping.timetag, outdata.timetag, /quad)

    ; get the calibrationSetId for the ADCGain 
    if keyword_set(verbose) then print,'Retrieving the profile integral ...'
    query1 = "SELECT calibrationSetId FROM CalibrationMetadata WHERE "
    query1 = query1+"calibrationTableName='SimDiodeAdcGainCalibrationData' AND instrumentModeId="+strtrim(string(instrumentModeId),2)
    query1 = query1+" AND effectiveDate <= '"+jd2syb(gps2jd(t0/1d6))+"'"
    if keyword_set(version) then query1=query1+" AND version<="+strtrim(string(version),2)
    query1 = query1+" ORDER BY version DESC"
    res=oJavaDbExchange->getAllValues(query1[0])
    if size(res,/tname) ne 'STRUCT' then begin
        print,'Error: no CalibrationMetadata found for the requested date/version'
        return,-1
    endif
    calibrationSetId=res[0].calibrationSetId

    query2="SELECT * FROM SimDiodeAdcGainCalibrationData WHERE calibrationSetId="+strtrim(string(calibrationSetId),2)
    adcgain_data=oJavaDbExchange->getAllValues(query2)
    adcgain = adcgain_data.diodeadcgain

    query1 = "SELECT calibrationSetId FROM CalibrationMetadata WHERE "
    query1 = query1+"calibrationTableName='SimDiodeFeedbackResistanceCal' AND instrumentModeId="+strtrim(string(instrumentModeId),2)
    query1 = query1+" AND effectiveDate <= '"+jd2syb(gps2jd(t0/1d6))+"'"
    if keyword_set(version) then query1=query1+" AND version<="+strtrim(string(version),2)
    query1 = query1+" ORDER BY version DESC"
    res=oJavaDbExchange->getAllValues(query1[0])
    if size(res,/tname) ne 'STRUCT' then begin
        print,'Error: no CalibrationMetadata found for the requested date/version'
        return,-1
    endif
    calibrationSetId=res[0].calibrationSetId
    query2="SELECT * FROM SimDiodeFeedbackResistanceCal WHERE calibrationSetId="+strtrim(string(calibrationSetId),2)
    feedbackR_data=oJavaDbExchange->getAllValues(query2)
    feedbackResistance = feedbackR_data.resistanceInOhms
    outdata.irradiance = outdata.dn / adcGain / feedbackResistance

    ; the irradiance is defined as:
    ; dnData = (AveragedDN - DarkDN) / ADCGain /feedBackResistance
    ; UncorrectedIrradiance = dnData / sensitivityIntegral
    ; correctedIrradiance = unCorrectedIrradiance / sunObserverDistance / Doppler / DetectorDeg / PrismDeg 
    ; calibratedIrradiance = correctedIrradiance * alephFactorValue
    ; 
    ; for now we skip the Doppler, DetectorDeg and PrismDeg, alephFactorValue
    ;


    ; get the profileIntegralCal data for this instrumentMode
    query1 = "SELECT calibrationSetId FROM CalibrationMetadata WHERE "
    query1 = query1+"calibrationTableName='SimProfileIntegralCal' AND instrumentModeId="+strtrim(string(instrumentModeId),2)
    ; we will ALWAYS use the same SimProfileIntegralCal ID as the reference spectra
    query1 = query1+" AND effectiveDate <= '"+jd2syb(gps2jd(refSpect[0].timetag/1d6))+"'"
    if keyword_set(version) then query1=query1+" AND version<="+strtrim(string(version),2)
    query1 = query1+" ORDER BY version DESC"
    res=oJavaDbExchange->getAllValues(query1[0])
    if size(res,/tname) ne 'STRUCT' then begin
        print,'Error: no CalibrationMetadata found for the requested date/version'
        return,-1
    endif
    calibrationSetId=res[0].calibrationSetId
    ; now get the SimProfileIntegralCal entry for this calibrationSetId
    query2="SELECT * FROM SimProfileIntegralCal WHERE calibrationSetId="+strtrim(string(calibrationSetId),2)
    profile_data=oJavaDbExchange->getAllValues(query2)
    if size(profile_data,/tname) ne 'STRUCT' then begin
        print,'Error: no SimProfileIntegralCal found for the requested date/version'
        return,-1
    endif
    spotRef = profile_data.x
    ; the reference temperature is really the same for every row inSimProfileIntegralCal
    prismRefTemp = mean(profile_data.y0)
    wavelengthRef = profile_data.y4
    bandCenterPrismTCoef = profile_data.y13 

    ; the degColumn comes from SimPrismDegColumnCalTable
    ; get the calibrationSetId for the degColumn 
    if instrumentModeId le 44 then simMode=54 else simMode=55   ; entry is for 54 and 55 for simA and simB
    query1 = "SELECT calibrationSetId FROM CalibrationMetadata WHERE "
    query1 = query1+"calibrationTableName='SimPrismDegColumnCalTable' AND instrumentModeId="+strtrim(string(simMode),2)
    query1 = query1+" AND effectiveDate <= '"+jd2syb(gps2jd(t0/1d6))+"'"
    if keyword_set(version) then query1=query1+" AND version<="+strtrim(string(version),2)
    query1 = query1+" ORDER BY version DESC"
    res=oJavaDbExchange->getAllValues(query1[0])
    if size(res,/tname) ne 'STRUCT' then begin
        print,'Error: no CalibrationMetadata found for SimPrismDegColumnCalTable and the requested date/version'
        print,query1
        return,-1
    endif
    calibrationSetId=res[0].calibrationSetId
    ; now get the degColumn entry for this calibrationSetId
    ; the X column is in seconds and not microSeconds
    query2="SELECT * FROM SimPrismDegColumnCalTable WHERE calibrationSetId="+strtrim(string(calibrationSetId),2)
    degColumn_data=oJavaDbExchange->getAllValues(query2)
    degColumn = interpol(degColumn_data.y, degColumn_data.x, outdata.timetag/1d6, /quad)


    ; correct for 1AU and Doppler
    ; the instrumentId for SORCE is 50 and we're at version=6  (querying for the max(version) is VERY SLOW)
    query2="SELECT * FROM SolarDistAndDopplerFixedStep2 WHERE instrumentModeId=50 and version=7 and "
    query2 = query2+"microsecondsSinceGpsEpoch>="+strtrim(string(ulong64(t0)),2)+" and "
    query2 = query2+"microsecondsSinceGpsEpoch<="+strtrim(string(ulong64(t1)),2)
    solarDistDop=oJavaDbExchange->getAllValues(query2)
    if size(solarDistDop,/tname) ne 'STRUCT' then return,-1
    if n_elements(solarDistDop.sunObserverDistanceCorrection) le 3 then $
        sunObserverDistance = interpol(solarDistDop.sunObserverDistanceCorrection, $
            solarDistDop.microsecondsSinceGpsEpoch, outdata.timetag) $
    else $ 
        sunObserverDistance = interpol(solarDistDop.sunObserverDistanceCorrection, $
            solarDistDop.microsecondsSinceGpsEpoch, outdata.timetag,/quad)


    ; we assume the separation between dwell points has the same number of sub-pixel
    pix_steps = refSpect[1].ccdpos - refSpect[0].ccdpos

    ; apply all the wavelength related corrections and measure the cross-correlation
    ; with our refSpect 
    ; get the initial offset from the raw dn
    ; interpolate to the refspect grid in case some points are missing
    int_dn = interpol(outdata.dn,outdata.ccdpos,refspect.ccdpos,/quad)
    if deriv_order eq 0 then begin
        cross_correlate, refSpect.dn, int_dn, offset0, corr, width=30
    endif else if deriv_order eq 1 then begin
        cross_correlate, DERIV(refSpect.ccdpos,refSpect.dn), DERIV(refSpect.ccdpos,int_dn), offset0, corr, width=30
    endif else if deriv_order eq 2 then begin
        cross_correlate, DERIV(refSpect.ccdpos,DERIV(refSpect.ccdpos,refSpect.dn)), $
            DERIV(refSpect.ccdpos,DERIV(refSpect.ccdpos,int_dn)), offset0, corr, width=30
    endif else begin
        cross_correlate, DERIV(refSpect.ccdpos,DERIV(refSpect.ccdpos,DERIV(refSpect.ccdpos,refSpect.dn))), $
            DERIV(refSpect.ccdpos,DERIV(refSpect.ccdpos,DERIV(refSpect.ccdpos,int_dn))), offset0, corr, width=30
    endelse
    ; we'll be testing for 20 sub-pixels around the nominal offset
    if keyword_set(verbose) then print,'Initial offset = '+strtrim(string(offset0*pix_steps),2)
    offset_list = (dindgen(40)/1.0 -20)*1.0 + offset0*pix_steps
    offsets = offset_list * 0.0d
    fit_goodness = offsets * 0.0d

    tempcorr = dblarr(n_elements(outdata.wavelength))
    ;lineplot,refSpect.wavelength,refspect.irradiance,title='Reference Spectra',psym=-3
    ;lineplot,refSpect.ccdpos,refspect.dn,title='Reference Spectra',psym=-3
    ;lineplot,outdata.ccdpos,outdata.dn,title='Base Spectra',psym=-3
    mn=min(refSpect.irradiance,max=mx,pos)
    norm_ref_irrad = 1.0d - (1.0d - 0.001d) * (mx-refSpect.irradiance) / (mx-mn)
    if deriv_order eq 0 then begin
        deriv_ref_irrad = norm_ref_irrad
    endif else if deriv_order eq 1 then begin
        deriv_ref_irrad = DERIV(refSpect.wavelength,norm_ref_irrad)
    endif else if deriv_order eq 2 then begin
        deriv_ref_irrad = DERIV(refSpect.wavelength,DERIV(refSpect.wavelength,norm_ref_irrad))
    endif else begin
        deriv_ref_irrad = DERIV(refSpect.wavelength,DERIV(refSpect.wavelength,DERIV(refSpect.wavelength,norm_ref_irrad)))
    endelse

    ; use this to test the TSIS prism temperature correction
    ; dldt= get_dwl_dt(detector,channel)

    for i=0, n_elements(offsets)-1 do begin
        outdata.ccdpos = outdata.ccdpos_uncorr + offset_list[i]
        outdata.wavelength = interpol(wavelengthRef, spotRef, outdata.ccdpos,/quad)
        if NOT keyword_set(noTempCorr) then begin
            dndt = interpol(bandCenterPrismTCoef, spotRef, outdata.ccdpos, /quad)
            ; test the TSIS prism temp correction for giggles
            ; dndt = interpol(dldt.dldt, dldt.wavelength, outdata.wavelength,/quad)
            outdata.wavelength = outdata.wavelength + (outdata.prismTemp - prismRefTemp) * dndt
            ; sim2_wavelength_to_subpixel_2d_raytrace, outdata.wavelength, spot, channel, detector, subpixel1, temp_prism=prismRefTemp
            ; sim2_subpixel_to_wavelength_2d_raytrace, subpixel1, spot, channel, detector, wave1, temp_prism=outdata.prismTemp
            ; outdata.wavelength=wave1
        endif

        ; get the correction for each wavelength by calculating the correction for the wavelengths
        ; smaller and greater than required and do a linear interpolation
        for j=0L,n_elements(outdata.wavelength)-1L do begin
            p0=where(wv le outdata[j].wavelength,count0)
            lopos=p0[-1]
            p1=where(wv ge outdata[j].wavelength,count1)
            hipos=p1[0]
            corrLo=poly(outdata[j].detectorTemp,[c0[lopos],c1[lopos],c2[lopos],c3[lopos]])
            corrHi=poly(outdata[j].detectorTemp,[c0[hipos],c1[hipos],c2[hipos],c3[hipos]])
            if count0 gt 0 and count1 gt 0 then begin
                ; interpolate for the desired wavelength
                tempcorr[j] = corrLo - (wv[hipos] - outdata[j].wavelength) * (corrHi - corrLo) / (wv[hipos] - wv[lopos])
            endif else if count0 gt 0 then begin
                ; use the lower wavelength value
                tempcorr[j]=corrLo
            endif else if count1 gt 0 then begin
                ; use the higher wavelength value
                tempcorr[j]=corrHi
            endif
        endfor

        new_irradiance = outdata.irradiance + tempcorr

        ; get the sensitivityIntegral = y8 + y10*(prismTemp - y0) + y11*(degColumn - y1) + y12*(detTemp - y2)
        y0 = prismRefTemp                                            ; reference prism temperature
        y1 = interpol(profile_data.y1, spotRef, outdata.ccdpos)      ; Prism Degradation column
        y2 = mean(profile_data.y2)                                   ; reference detector temperature
        y8 = interpol(profile_data.y8, spotRef, outdata.ccdpos)      ; profile integral
        y10 = interpol(profile_data.y10, spotRef, outdata.ccdpos)    ; profIntegralPrismTemperatureCoef
        y11 = interpol(profile_data.y11, spotRef, outdata.ccdpos)    ; profIntegralDegradationCoef 
        y12 = interpol(profile_data.y12, spotRef, outdata.ccdpos)    ; profIntegralDetTemperatureCoef

        ; the diode temperature correction is now applied to the dark subtracted dn
        sensitivityProfile = y8 + y10*(outdata.prismTemp - y0) +y11*(degColumn-y1); + y12*(outdata.detectorTemp-y2)
        new_irradiance /= sensitivityProfile
        new_irradiance /= sunObserverDistance

        ; get the cross-correlation and the quality of the fit after resampling to the same wavelength grid
        ;lineplot,outdata.ccdpos,outdata.dn,psym=-3,title='Calib Spectra '+strtrim(string(offset_list[i]),2)
        ;new_irradiance = interpol(new_irradiance, outdata.wavelength, refSpect.wavelength,/spline)
        ;new_irradiance = interpol(new_irradiance, outdata.wavelength, refSpect.wavelength,/lsq)
        new_irradiance = interpol(new_irradiance, outdata.wavelength, refSpect.wavelength,/quad)
        mn=min(new_irradiance,max=mx,pos)
        norm_new_irrad = 1.0d - (1.0d - 0.001d) * (mx-new_irradiance) / (mx-mn)
        ;lineplot,refSpect.wavelength,new_irradiance,psym=-3,title='Calib Spectra '+strtrim(string(offset_list[i]),2)
        ;lineplot,refSpect.wavelength,(refspect.irradiance-outdata.irradiance),title=strtrim(string(offset_list[i]),2)
        ; get the variance of the difference of the normalized polts
        if deriv_order eq 0 then begin
            deriv_irrad = norm_new_irrad
        endif else if deriv_order eq 1 then begin
            deriv_irrad = DERIV(refSpect.wavelength,norm_new_irrad)
        endif else if deriv_order eq 2 then begin
            deriv_irrad = DERIV(refSpect.wavelength,DERIV(refSpect.wavelength,norm_new_irrad))
        endif else begin
            deriv_irrad = DERIV(refSpect.wavelength,DERIV(refSpect.wavelength,DERIV(refSpect.wavelength,norm_new_irrad)))
        endelse

        if n_elements(wpos) eq 0 then begin
            cross_correlate, deriv_ref_irrad, deriv_irrad, offset, corr, width=31
        endif else begin
            cross_correlate, deriv_ref_irrad[wpos], deriv_irrad[wpos], offset, corr, width=31
        endelse
        offsets[i] = offset * pix_steps
        fit_goodness[i] = max(corr)
        ;lineplot,corr,psym=-3,title=strtrim(string(offset_list[i]),2)
        if keyword_set(verbose) then print,'Processing offset '+strtrim(string(offset_list[i]),2)+', '+$
            strtrim(string(offsets[i]),2)+', '+strtrim(string(fit_goodness[i]),2)

    endfor
    
    if not keyword_set(noplot) then $
        lineplot, offset_list, fit_goodness, psym=-4, charsize=1.2, xtitle='CCD Offset (subpixels)',$
        ytitle='Maximum Correlation',title='InstrumentModeId='+strtrim(string(instrumentModeId),2)+$
        ' from '+strtrim(string(gps2sd(t0/1d6)),2)+' to '+strtrim(string(gps2sd(t1/1d6)),2)+$
        ' DERIV='+strtrim(string(deriv_order),2)
    ; fine tune the peak by re-sampling to 0.1 pix
    mx=max(fit_goodness,pos)
    delta=(offset_list[1]-offset_list[0])
    xpos=(dindgen(600)/100.0)*delta + round(offset_list[pos]) - 3.0*delta
    ypos=interpol(fit_goodness,offset_list,xpos,/quad)
    mx=max(ypos,pos)
    
    return, [xpos[pos], max(fit_goodness)]

end
