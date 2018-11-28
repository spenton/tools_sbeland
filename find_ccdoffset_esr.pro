;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Align spectra with a reference spectra by changing the ccd subpixel
;   position and minimizing the residuals (optionally of the first or second
;   derivatives). We're comparing the irradiance for each spectra.
;
;   Wrapper for the find_ccdoffset routine specifically for the ESR data
;   using the SORCE development database version 0.
;
; CALLING SEQUENCE:
;   result = find_ccdoffset_esr(t0, t1, instrumentModeId, refSpect=refSpect, /missionDays)
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
;   Revision: $Id: find_ccdoffset_esr.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;-
;
function find_ccdoffset_esr, startTime, stopTime, instrumentModeId, $
         gps=gps, missionDays=missionDays, julianDays=julianDays, $
         refT0=refT0, refT1=refT1, noTempCorr=noTempCorr, $
         version=version, no1AU=no1AU, verbose=verbose, $
         deriv_order=deriv_order, noplot=noplot, refSpect=refSpect

    nrows=0
    if n_elements(deriv_order) eq 0 then deriv_order=2
    ; validate the instrumentModeId
    if n_elements(instrumentModeId) eq 0 then begin
        doc_library,'find_ccdoffset_esr'
        print,''
        print,'Missing instrumentModeId '
        return,-1
    endif else begin
        ; verify we got the right instrument
        valid_inst = [31,32]
        instrumentModeId = fix(instrumentModeId)
        if (where(valid_inst eq instrumentModeId))[0] lt 0 then  begin
            doc_library,'find_ccdoffset'
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
            detector='esr'
            channel='a'
            prismId = 101862L
            prismTempId = 101909L
            detectorId = 100058L
        end
        32: begin
            ; SIMB ESR
            instrument='SIMB'
            detector='esr'
            channel='b'
            prismId = 100478L
            prismTempId = 102212L
            detectorId = 100048L
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
        rT0=sd2gps(453.02066d)*1d6
        rT1=sd2gps(453.98990d)*1d6
    endif

    user='sbeland' 
    password='sbeland1'
    dbUrl='jdbc:sybase:Tds:lasp-db-dev:4100/SORCE_SIM_V19'
    dbDriver='com.sybase.jdbc3.jdbc.SybDriver' 

    if n_elements(version) eq 0 then version=0

    ; get the reference spectra
    if size(refSpect,/tname) ne 'STRUCT' then begin
        if keyword_set(verbose) then print,'Extracting the reference spectra ...'
        spect = get_science_product(['Wavelength','SimCalibratedIrradiance','PhaseDetectedDN'],rT0,rT1,$
            instrumentModeId, /gps, version=version,dburl=dburl,dbdriver=dbdriver,user=user,password=password)
        refSpect=replicate({timetag:0d, ccdpos:0d, wavelength:0d, dn:0d, irradiance:0d},n_elements(spect))
        s=sort(spect.wavelengthRef)
        spect=spect[s]
        refSpect.timetag = spect.MICROSECONDSSINCEGPSEPOCH
        refSpect.ccdpos = spect.PRISMPOSITION
        refSpect.wavelength = spect.WAVELENGTHREF
        refSpect.irradiance = spect.IRRADIANCE
        refSpect.dn = spect.REALCOMPONENT
    endif


    ; we assume the separation between dwell points has the same number of sub-pixel
    pix_steps = refSpect[1].ccdpos - refSpect[0].ccdpos

    mn=min(refSpect.irradiance,max=mx,pos)
    norm_ref_irrad = 1.0d - (1.0d - 0.001d) * (mx-refSpect.irradiance) / (mx-mn)
    if deriv_order eq 0 then begin
        deriv_ref_irrad = norm_ref_irrad
    endif else if deriv_order eq 1 then begin
        deriv_ref_irrad = DERIV(refSpect.wavelength,norm_ref_irrad)
    endif else begin
        deriv_ref_irrad = DERIV(refSpect.wavelength,DERIV(refSpect.wavelength,norm_ref_irrad))
    endelse

    ; get rid of spikes in the ref spectra
    resistant_mean,deriv_ref_irrad,9.0,mean,goodvec=keep0

    ; limit the wavelength coverage
    ;w0=320d & w1=550d
    w0=280d & w1=500d
    wpos = where(refspect[keep0].wavelength ge w0 and refSpect[keep0].wavelength le w1)
    wpos=keep0[wpos]

    rSpect=refSpect[wpos]
    deriv_ref_irrad=deriv_ref_irrad[wpos]

    ; get all of the ESRFullScan within the requested time range
    if instrumentModeId eq 31 then begin
        plans = get_sorce_plan(t0,t1,/gps,/sima,activity='ESRFullScan%')
    endif else if instrumentModeId eq 32 then begin
        plans = get_sorce_plan(t0,t1,/gps,/simb,activity='ESRFullScan%')
    endif
    p0=0
    p1=1
    t0=[] & t1=[]
    while p0 lt n_elements(plans) do begin
        ; group plans within 3 full days (microsec)
        p=where(plans[p0:-1].stoptime - plans[p0].starttime lt 86400d*3d,count)
        t0 = [t0,plans[p0].starttime*1d6]
        t1 = [t1,plans[p[-1]+p0].stoptime*1d6]
        p0=p0 + p[-1] + 1
    endwhile

    if n_elements(t0) eq 0 then begin
        print,'Error: no ESRFullScan found within the time range'
        return,-1
    endif

    ; prepare the database connection 
    jstmt = fjava_get_jdbc_statement(server=server, database="SORCE")
    oJavaDbExchange = OBJ_NEW("oJava_DbExchange", jstmt)

    ; get each spectra at a time and perform the wavelength correction due to
    ; prism temperature (and corresponding irradiance change)
    nspect =  n_elements(t0)
    result=replicate({starttime:0d, offset:0d, fit_goodness:0d}, nspect)

    for sp=0,nspect-1 do begin

        if keyword_set(verbose) then print,'Retrieving the spectra from  ',gps2sd(t0[sp]/1d6)
        spect = get_science_product(['Wavelength','SimCalibratedIrradiance','SimConvertedDataNumbers'],t0[sp],t1[sp],$
            instrumentModeId, /gps, version=version,dburl=dburl,dbdriver=dbdriver,user=user,password=password)

        outdata = replicate({timetag:0.0d, ccdpos:0.0d, ccdpos_uncorr:0.0d, wavelength:0.0d, dn:0.0d, $
                  dn_tempcorr:0.0d, irradiance:0.0d, prismTemp:0.0d, detectorTemp:0.0d},n_elements(spect))
        outdata.timetag = spect.MICROSECONDSSINCEGPSEPOCH
        outdata.ccdpos_uncorr = spect.PRISMPOSITION
        outdata.ccdpos = outdata.ccdpos_uncorr
        outdata.wavelength = spect.WAVELENGTHREF
        outdata.irradiance = spect.IRRADIANCE
        outdata.dn = spect.REALCOMPONENT
        outdata.dn_tempcorr = spect.REALCOMPONENT

        ; get the prism temperature for these timetags
        if keyword_set(verbose) then print,'Retrieving the prismTemp ...'
        get_sorce_telemetry,tlm_data,info,t0[sp]-120.0d6,t1[sp]+120.0d6,/gps,tlmId=[prismId,prismTempId]
        outdata.prismTemp = interpol((*tlm_data.(1)).housekeeping.eu, (*tlm_data.(1)).housekeeping.timetag, outdata.timetag) ;, /quad)

        ; correct for 1AU and Doppler
        ; the instrumentId for SORCE is 50 and we're at version=6  (querying for the max(version) is VERY SLOW)
        if keyword_set(verbose) then print,'Retrieving the sun-observer distance ...'
        query_database, /reset
        query2="SELECT microsecondsSinceGpsEpoch,sunObserverDistanceCorrection FROM SolarDistAndDopplerFixedStep2 "
        query2 = query2+"WHERE instrumentModeId=50 and version=7 and "
        query2 = query2+"microsecondsSinceGpsEpoch>="+strtrim(string(ulong64(t0[sp])),2)+" and "
        query2 = query2+"microsecondsSinceGpsEpoch<="+strtrim(string(ulong64(t1[sp])),2)
        solarDistDop=oJavaDbExchange->getAllValues(query2)
        if size(solarDistDop,/tname) ne 'STRUCT' then return,-1
        sunObserverDistance = interpol(solarDistDop.sunObserverDistanceCorrection, $
            solarDistDop.microsecondsSinceGpsEpoch, outdata.timetag) 

        ; the degColumn comes from SimPrismDegColumnCalTable
        ; get the calibrationSetId for the degColumn 
        if keyword_set(verbose) then print,'Retrieving the SimPrismDegColumnCalTable ...'
        if instrumentModeId eq 31 then simMode=54 else simMode=55   ; entry is for 54 and 55 for simA and simB
        query1 = "SELECT calibrationSetId FROM CalibrationMetadata WHERE "
        query1 = query1+"calibrationTableName='SimPrismDegColumnCalTable' AND instrumentModeId="+strtrim(string(simMode),2)
        query1 = query1+" AND effectiveDate <= '"+jd2syb(gps2jd(t0[sp]/1d6))+"'"
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
        degColumn = interpol(degColumn_data.y, degColumn_data.x, outdata.timetag/1d6) ;, /quad)

        ; get the profileIntegralCal data for this instrumentMode
        if keyword_set(verbose) then print,'Retrieving the profile integral ...'
        query1 = "SELECT calibrationSetId FROM CalibrationMetadata WHERE "
        query1 = query1+"calibrationTableName='SimProfileIntegralCal' AND instrumentModeId="+strtrim(string(instrumentModeId),2)
        ; we will ALWAYS use the same SimProfileIntegralCal ID as the reference spectra
        query1 = query1+" AND effectiveDate <= '"+jd2syb(gps2jd(rSpect[0].timetag/1d6))+"'"
        ;if keyword_set(version) then query1=query1+" AND version<="+strtrim(string(version),2)
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
        ; the reference temperature is really the same for every row in SimProfileIntegralCal
        prismRefTemp = mean(profile_data.y0)
        wavelengthRef = profile_data.y4
        bandCenterPrismTCoef = profile_data.y13 

        ; sort the outdata with wavelength now that we got the time dependant stuff
        s=sort(outdata.wavelength)
        outdata=outdata[s]
        sunObserverDistance=sunObserverDistance[s]

        ; interpolate to the rSpect grid in case some points are missing
        int_dn = interpol(outdata.irradiance,outdata.ccdpos,rSpect.ccdpos) ;,/quad)
        mn=min(int_dn,max=mx,pos)
        temp = 1.0d - (1.0d - 0.001d) * (mx-int_dn) / (mx-mn)
        if deriv_order eq 0 then begin
            deriv_irrad = temp
        endif else if deriv_order eq 1 then begin
            deriv_irrad = DERIV(rSpect.wavelength,temp)
        endif else begin
            deriv_irrad = DERIV(rSpect.wavelength,DERIV(rSpect.wavelength,temp))
        endelse
        cross_correlate, deriv_ref_irrad, deriv_irrad, offset0, corr, width=31
        if keyword_set(verbose) then print,'Initial offset = '+strtrim(string(offset0*pix_steps),2)
        ; we'll be testing for 20 sub-pixels around the nominal offset
        offset_list = dindgen(80) -40d + offset0*pix_steps
        offsets = offset_list * 0.0d
        fit_goodness = offsets * 0.0d

        tempcorr = dblarr(n_elements(outdata.wavelength))
        for i=0, n_elements(offsets)-1 do begin
            outdata.ccdpos = outdata.ccdpos_uncorr + offset_list[i]
            outdata.wavelength = interpol(wavelengthRef, spotRef, outdata.ccdpos) ;,/quad)
            if NOT keyword_set(noTempCorr) then begin
                dndt = interpol(bandCenterPrismTCoef, spotRef, outdata.ccdpos); , /quad)
                outdata.wavelength = outdata.wavelength + (outdata.prismTemp - prismRefTemp) * dndt
            endif

            new_irradiance = outdata.irradiance
            ;s=sort(outdata.wavelength)
            ;new_irradiance = outdata[s].dn

            ; get the sensitivityIntegral = y8 + y10*(prismTemp - y0) + y11*(degColumn - y1) + y12*(detTemp - y2)
            y0 = prismRefTemp                                            ; reference prism temperature
            y1 = interpol(profile_data.y1, spotRef, outdata.ccdpos)      ; Prism Degradation column
            y2 = mean(profile_data.y2)                                   ; reference detector temperature
            y8 = interpol(profile_data.y8, spotRef, outdata.ccdpos)      ; profile integral
            y10 = interpol(profile_data.y10, spotRef, outdata.ccdpos)    ; profIntegralPrismTemperatureCoef
            y11 = interpol(profile_data.y11, spotRef, outdata.ccdpos)    ; profIntegralDegradationCoef 
            y12 = interpol(profile_data.y12, spotRef, outdata.ccdpos)    ; profIntegralDetTemperatureCoef

            ; the diode temperature correction is now applied to the dark subtracted dn
            sensitivityProfile = y8 + y11*(degColumn-y1)
            ;new_irradiance /= sensitivityProfile
            ;new_irradiance /= sunObserverDistance

            ; get the cross-correlation and the quality of the fit after resampling to the same wavelength grid
            ;lineplot,outdata.ccdpos,outdata.dn,psym=-3,title='Calib Spectra '+strtrim(string(offset_list[i]),2)
            ;new_irradiance = interpol(new_irradiance, outdata.wavelength, rSpect.wavelength,/spline)
            ;new_irradiance = interpol(new_irradiance, outdata.wavelength, rSpect.wavelength,/lsq)
            temp = interpol(new_irradiance, outdata.wavelength, rSpect.wavelength) ;,/quad)
            mn=min(temp,max=mx,pos)
            p=where(finite([mn,mx]) eq 0,count)
            ; if interpol exploded do a linear interpol
            if count gt 0 then begin
                temp = interpol(new_irradiance, outdata.wavelength, rSpect.wavelength)
                mn=min(temp,max=mx,pos)
            endif
            new_irradiance=temp
            norm_new_irrad = 1.0d - (1.0d - 0.001d) * (mx-new_irradiance) / (mx-mn)
            ;lineplot,rSpect.wavelength,new_irradiance,psym=-3,title='Calib Spectra '+strtrim(string(offset_list[i]),2)
            ;lineplot,rSpect.wavelength,(rSpect.irradiance-outdata.irradiance),title=strtrim(string(offset_list[i]),2)
            ; get the variance of the difference of the normalized polts
            if deriv_order eq 0 then begin
                deriv_irrad = norm_new_irrad
            endif else if deriv_order eq 1 then begin
                deriv_irrad = DERIV(rSpect.wavelength,norm_new_irrad)
            endif else begin
                deriv_irrad = DERIV(rSpect.wavelength,DERIV(rSpect.wavelength,norm_new_irrad))
            endelse

            cross_correlate, deriv_ref_irrad, deriv_irrad, offset, corr, width=31
            offsets[i] = offset * pix_steps
            fit_goodness[i] = max(corr)
            ;lineplot,corr,psym=-3,title=strtrim(string(offset_list[i]),2)
            if keyword_set(verbose) then print,'Processing offset '+strtrim(string(offset_list[i]),2)+', '+$
                strtrim(string(offsets[i]),2)+', '+strtrim(string(fit_goodness[i]),2)

        endfor
        
        ; fine tune the peak by re-sampling to 0.1 pix
        mx=max(fit_goodness,pos)
        delta=(offset_list[1]-offset_list[0])
        xpos=(dindgen(600)/100.0)*delta + round(offset_list[pos]) - 3.0*delta
        ypos=interpol(fit_goodness,offset_list,xpos) ;,/quad)
        mx=max(ypos,pos)

        result[sp].starttime=t0[sp]
        result[sp].offset=xpos[pos]
        result[sp].fit_goodness=max(fit_goodness)

        if not keyword_set(noplot) then $
            lineplot, offset_list, fit_goodness, psym=-4, charsize=1.2, xtitle='CCD Offset (subpixels)',$
            ytitle='Maximum Correlation',title='ModeId='+strtrim(string(instrumentModeId),2)+$
            ' from '+strtrim(string(gps2sd(t0[sp]/1d6)),2)+' to '+strtrim(string(gps2sd(t1[sp]/1d6)),2)+$
            ' DERIV='+strtrim(string(deriv_order),2)+', '+strtrim(string(xpos[pos],format='(F0.1)'),2)
    endfor


    return, result

end
