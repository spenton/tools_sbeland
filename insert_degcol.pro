;+
; NAME:   INSERT_DEGCOL
;
; PURPOSE: 
;    This routine populates the SORCE/SimPrismDegColumnCalTable database table with the 
;    new degradation column values obtained from comparing the degradation of ESRA and B
;    and running get_global_degcol.pro.
;
; CALLING SEQUENCE:
;    insert_degcol, 54, /dbinsert
;
; OPTIONAL INPUT PARAMETERS:
;    none
;
; OPTIONAL INPUT KEYWORDS:
;    dbinsert - if set will insert the re-formatted data in the database
;
; OUTPUT PARAMETERS:
;    none
;
; OPTIONAL OUTPUT PARAMETERS:
;    NONE
;
; EXAMPLE:  
;
; COMMON BLOCKS:
;    NONE
;
; NOTES:
;
; REVISION HISTORY:
;   Revision: $Id: insert_degcol.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;-
;*****************************************************************************
PRO update_calibMetadata_degcol, instMode, calId, effDate, version=version, conn17=conn17, $
    verbose=verbose, dbinsert=insert

    ; make sure we have an entry in the CalibrationMetadata table for this calibrationSetId
    q1="SELECT * from CalibrationMetadata where calibrationSetId="+strtrim(string(CalId),2)
    if insert then begin
        query_database,/reset
        query_database,q1,data,nrows
        date_str=string(systime(/jul),format='(C(CYI,"-",CMOI2.2,"-",CDI2.2," ",CHI2.2,":",CMI2.2,":",CSI2.2))')
        if nrows eq 0 then begin
            ; insert a new row in the CalibrationMetadata table
            q1="INSERT INTO dbo.CalibrationMetadata(calibrationSetId, calibrationTableName, instrumentModeId, "
            q1+="version, releaseDate, effectiveDate, rationaleForNewVersion, responsibleIndividual, "
            q1+="dataDescription, sourceOfData, applicableReferences, additionalInformation) VALUES("
            q1=q1+strtrim(string(CalId),2)+", 'SimPrismDegColumnCalTable', "
            q1=q1+strtrim(string(instMode),2)+", "+strtrim(string(version),2)
            q1=q1+", '"+date_str+"', '"+effDate+"'"
            q1=q1+", 'New Degradation Column for version 20 from updated solar exposures, 1AU and F'"
            q1=q1+", 'Stephane Beland'"
            q1=q1+", 'Stores an arbitrary number of table values (x, y, yUncertainty) for a TabularFunction1d "
            q1=q1+"object representing the SIM-A column degradation, as a function of time (seconds since GPS epoch in this version)'"
            q1=q1+", 'Stephane processing and ESR to diode comparison'"
            q1=q1+", 'Generated from best solar_exp, 1AU correction and F function'"
            q1=q1+", 'none')"
            if keyword_set(verbose) then print,q1
            if insert then result=conn17->execute(q1)
        endif else begin
            q1="UPDATE CalibrationMetadata SET releaseDate='"+date_str+"'"
            q1=q1+", effectiveDate='"+effDate+"'"
            q1=q1+", version="+strtrim(string(version),2)
            q1=q1+" WHERE calibrationSetId="+strtrim(string(CalId),2)
            if keyword_set(verbose) then print,q1
            if insert then result=conn17->execute(q1)
        endelse
    endif

END

;-----------------------------------------------------------------------------
FUNCTION INSERT_DEGCOL, instrumentModeId, dbinsert=dbinsert, verbose=verbose, $
    solar54=solar54, solar55=solar55, lya=lya, version=version

    ;modes = [31, 32, 54, 55]
    modes = [54, 55]
    pos = where(modes eq instrumentModeId, count)
    if count eq 0 then begin
        print,'Error: invalid InstrumentModeId ',instrumentModeId,' [',mode,']'
        return,-1
    endif

    if n_elements(version) eq 0 then begin
        version=19
        newCalBase=1033 
    endif else if version eq 20 then begin
        newCalBase=1123-2
    endif

    case instrumentModeId of
        31: begin
            newCalId=newCalBase
            end
        32: begin
            newCalId=newCalBase+1
            end
        54: begin
            newCalId=newCalBase+2
            end
        55: begin
            newCalId=newCalBase+3
            end
    endcase

    ; get the solar exposure data from the database
    get_solar54=0
    get_oneau=0
    get_lya=0
    if size(solar54,/tname) ne 'STRUCT' then begin
        get_solar54=1
        get_oneau=1
        get_lya=1
    endif else begin
        tnames = tag_names(solar54)
        p=where(strpos(tnames,"ONEAU") ge 0, count)
        if count eq 0 then get_oneau=1
        p=where(strpos(tnames,"LYALPHA_CORR") ge 0, count)
        if count eq 0 then get_lya=1
    endelse

    if get_solar54 eq 1 then begin
        print,'  getting solar54 from database ...'
        q1 = 'SELECT orbitStartTimeGps t0, orbitEndTimeGps t1, solarExposureHrtOut as solar_exp_orbit, '+$
             'cumulativeSolarExpHrtOut as solar_exp, combinedGapsDuration as gaps from SimSolarExposureData ' + $
             'where instrumentModeId=54'
        query_database, q1, solar54, info
    endif

    get_solar55=0
    if size(solar55,/tname) ne 'STRUCT' then begin
        get_solar55=1
        get_oneau=1
        get_lya=1
    endif else begin
        tnames = tag_names(solar55)
        p=where(strpos(tnames,"ONEAU") ge 0, count)
        if count eq 0 then get_oneau=1
        p=where(strpos(tnames,"LYALPHA_CORR") ge 0, count)
        if count eq 0 then get_lya=1
    endelse

    if get_solar55 eq 1 then begin
        print,'  getting solar55 from database ...'
        q1 = 'SELECT orbitStartTimeGps t0, orbitEndTimeGps t1, solarExposureHrtOut as solar_exp_orbit, '+$
             'cumulativeSolarExpHrtOut as solar_exp, combinedGapsDuration as gaps from SimSolarExposureData ' + $
             'where instrumentModeId=55'
        query_database, q1, solar55, info
    endif

    ;time54 = (solar54.(0)+solar54.(1))/2d
    ;time55 = (solar55.(0)+solar55.(1))/2d
    ; the solar exposure time is valid at end of orbit
    time54 = solar54.(1)
    time55 = solar55.(1)

    ; get the 1-AU correction
    if get_oneau eq 1 then begin
        print,'  getting SUNOBSERVERDISTANCECORRECTION from database ...'
        q1="SELECT microsecondsSinceGpsEpoch, sunObserverDistanceCorrection FROM SolarDistAndDopplerFixedStep2 "+$
           "WHERE instrumentModeId=50 and version=7"
        query_database, q1, solardist, info
        oneau54 = interpol(solardist.SUNOBSERVERDISTANCECORRECTION, solardist.MICROSECONDSSINCEGPSEPOCH, time54)
        append_tag,solar54,'oneau',oneau54,/slim
        oneau55 = interpol(solardist.SUNOBSERVERDISTANCECORRECTION, solardist.MICROSECONDSSINCEGPSEPOCH, time55)
        append_tag,solar55,'oneau',oneau55,/slim
    endif

    if keyword_set(lya) and get_lya eq 1 then begin
        print,'  getting Lyman Alpha dosage from Level3 data ...'
        ; extract the Lyman alpha irradiance from SOLSTICE data normalizing around the solar minimum
        res=get_level3_spectrum(0d,5000d, 11, 24d, /mission, min_wave=121d, max_wave=122d, /released)
        ; normalized to irradiance from "quiet sun"
        quiet_time=sd2gps([2130d,2230d])*1d6
        p=where(res.(1) ge quiet_time[0] and res.(1) le quiet_time[1])
        avg_irrad = mean(res[p].irradiance)
        res.irradiance /= avg_irrad
        ; add the oneau to Lyman Alpha to have a more accurate "dosage" factor (in reality makes little difference)
        lyalpha_corr = interpol(res.irradiance, res.(1), time54)
        min_time=min(res.(1),pos)
        p=where(time54 lt min_time,count)
        if count gt 0 then lyalpha_corr[p]=res[pos].irradiance
        lyalpha_corr *= solar54.oneau
        append_tag,solar54,'LYALPHA_CORR',lyalpha_corr,/slim
        lyalpha_corr = interpol(res.irradiance, res.(1), time55)
        p=where(time54 lt min_time,count)
        if count gt 0 then lyalpha_corr[p]=res[pos].irradiance
        lyalpha_corr *= solar55.oneau
        append_tag,solar55,'LYALPHA_CORR',lyalpha_corr,/slim
    endif

    print,'  calculating the Degradation Column ...'
    if version eq 20 then begin
        ; the cumulative solar exposure is in days to go with our Kappa function
        corr54=1.0d
        corr55=1.0d

        ; add the "average" F_Function for all wavelengths (2014-04-30) as a function of mission day
        coeff_corr=[1.4005348d,  -0.00061723119d,   1.5586607d-07,   3.9341566d-11,  -1.3793075d-14]
        corr54 = poly(gps2sd(time54/1d6),coeff_corr)
        corr55 = poly(gps2sd(time55/1d6),coeff_corr)

        ;solar54.solar_exp=total(solar54.solar_exp_orbit*corr54,/cum)
        ;solar54.solar_exp/=86400d
        ;solar55.solar_exp=total(solar55.solar_exp_orbit*corr55,/cum)
        ;solar55.solar_exp/=86400d

        ; adding a 1AU contribution from what we determined analytically (4 seems to work best)
        solar54.solar_exp=total(solar54.solar_exp_orbit*corr54,/cum) / (1d +(solar54.oneau-1d)/4d)
        solar54.solar_exp/=86400d
        solar55.solar_exp=total(solar55.solar_exp_orbit*corr55,/cum) / (1d +(solar55.oneau-1d)/4d)
        solar55.solar_exp/=86400d

        ; these coefficients were obtained from the fractional differences between the integrated SIM
        ; data for the whole mission and the TIM data and fitting a 5th order polynomial
        ; apply the correction to SIMB as a test
        ;v20_coeff = [-0.00039542735d, 5.4600012d-06, -5.4198504d-09, 2.4542743d-12, -6.1814658d-16, 6.9183804d-20]
        ;v20_corr55 = poly(gps2sd(time55/1d6),v20_coeff)
        ;solar55.solar_exp *= (1d + v20_corr55/max(v20_corr55))   ; version 300
        ;solar55.solar_exp *= (1d - v20_corr55/max(v20_corr55))   ; version 302
        ;solar55.solar_exp *= (1d - v20_corr55 * 100d)   ; version 303
        ;solar55.solar_exp = (solar55.solar_exp - 4d) > 0d ; version 304

        ; version 305
        ; try adding the errors for 55
        ;solar55.solar_exp=total((solar55.solar_exp_orbit - solar55.gaps)*corr55,/cum) / (1d +(solar55.oneau-1d)/4d)
        ;solar55.solar_exp/=86400d
        ; version 306
        ;solar54.solar_exp=total((solar54.solar_exp_orbit + solar54.gaps)*corr54,/cum) / (1d +(solar54.oneau-1d)/4d)
        ;solar54.solar_exp/=86400d
        ;solar54.solar_exp = (solar54.solar_exp + 20d) > 0d ; version 307

        ; vesion 308 scaling solar54.solar_exp by alog(Irrad_sim/Irrad_tim)
        ;v20_coeff=[9.1201296d, -0.0078704640d, 3.0253863d-06, -2.8708480d-08, 4.0338162d-11, -2.1772839d-14, 5.2187623d-18, -4.6365170d-22]
        ;v20_corr54 = poly(gps2sd(time54/1d6),v20_coeff)
        ;solar54.solar_exp+=v20_corr54

        ; version 309 - simply add 10 days to solar54 to compare with version 308
        ;solar54.solar_exp+=12d

    endif else if keyword_set(lya) then begin
        ; when we take the Dosage into account, the time dependant F function is linear
        ; these coefficients come from get_global_degcol
        coeff_corr=[0.92418784d, -7.5035887d-05]
        ; coefficients were calculated in days 
        corr54 = poly(gps2sd(time54/1d6),coeff_corr)
        corr55 = poly(gps2sd(time55/1d6),coeff_corr)
        ; the cumulative solar exposure is in seconds
        solar54.solar_exp=total(solar54.solar_exp_orbit*solar54.LYALPHA_CORR*corr54,/cum) / (1d +(solar54.oneau-1d)/4d)
        solar54.solar_exp/=86400d
        solar55.solar_exp=total(solar55.solar_exp_orbit*solar55.LYALPHA_CORR*corr55,/cum) / (1d +(solar55.oneau-1d)/4d)
        solar55.solar_exp/=86400d
    endif else begin
        ; without the Dosage, the time dependant F function is a 2nd order polynomial (up to day 2805 at least)
        ; these coefficients come from get_global_degcol
        ;coeff_corr=[1.4907146d, -0.00069623171d, 1.6402941d-07]
        ; NEW F-FUNCTION coefficients from iterating over 20 times between Kappa and F-Function on ESR
        ;coeff_corr=[0.58673234d, -0.00027095508d,   6.3361721d-08]
        ;
        ; a new F-Function was obtained using a step55=365 to remove the 1AU effect and shown to be
        ; the same for ESR and VIS for wavelengths from 330 to 650nm
        coeff_corr=[1.5054297d, -0.00072414176d, 1.8150587d-07]

        ; coefficients were calculated in days 
        corr54 = poly(gps2sd(time54/1d6),coeff_corr)
        corr55 = poly(gps2sd(time55/1d6),coeff_corr)

        ; AS a test (SBeland 2013/07/18), we add the ratio of exposure of SimB/SimA to our F_Function (Version 69)
        ; we then get a new Kappa and raypath
        sum54 = total(solar54.solar_exp_orbit,/cum)
        sum55 = total(solar55.solar_exp_orbit,/cum)
        p=where(sum54 gt 0.0)
        ;corr54[p] += (sum55[p] / sum54[p])  ; Version 69
        ;corr55[p] += (sum55[p] / sum54[p])  ; Version 69
        ;corr54[p] *= (1d - sum55[p] / sum54[p])   ; Version 71
        ;corr55[p] *= (1d - sum55[p] / sum54[p])   ; Version 71

        ; the cumulative solar exposure is in days to go with our Kappa function
        solar54.solar_exp=total(solar54.solar_exp_orbit*corr54,/cum) / (1d +(solar54.oneau-1d)/4d)
        solar54.solar_exp/=86400d
        solar55.solar_exp=total(solar55.solar_exp_orbit*corr55,/cum) / (1d +(solar55.oneau-1d)/4d)
        solar55.solar_exp/=86400d

        ; As a different test, scale the solar_exp for SimB with the profile of the ratio B/A (Version 70)
        ;sum55 = interpol(solar55.solar_exp, solar55.t0, solar54.t0)
        ;p=where(solar54.solar_exp gt 0d)
        ;sum55[p] *= ((1d - sum55[p] / sum54[p]) > 0d)
        ;solar55.solar_exp = interpol(sum55, solar54.t0, solar55.t0)


    endelse


    stmt17 = fjava_get_jdbc_statement(user='sorce', password='sorcedb', server='sorce-db', database="SORCE")
    conn17 = OBJ_NEW("oJava_DbExchange", stmt17)
    if keyword_set(dbinsert) then insert=1 else insert=0


    ; first delete all entries in the table for this particular version (version 1 only for now)
    q1 = "DELETE FROM SimPrismDegColumnCalTable WHERE calibrationSetId="+strtrim(string(newcalId),2)
    if keyword_set(verbose) then print,q1
    if insert then result=conn17->execute(q1)

    ; make sure we have an entry in the CalibrationMetadata table for this calibrationSetId
    effectiveDate='2003-01-01 00:00:01.0'
    update_calibMetadata_degcol, instrumentModeId, newCalId, effectiveDate, $
        version=version, conn17=conn17, dbinsert=insert, verbose=verbose

    q1="SELECT * from CalibrationMetadata where calibrationSetId="+strtrim(string(newCalId),2)
    if insert then begin
        date_str=string(systime(/jul),format='(C(CYI,"-",CMOI2.2,"-",CDI2.2," ",CHI2.2,":",CMI2.2,":",CSI2.2))')
        query_database,q1,data,nrows
        if nrows eq 0 and version eq 20 then begin 
        endif else if nrows eq 0 then begin
            print,'enrty in CalibrationMetadata for calibrationSetId='+strtrim(strin(newCalId),2)+' does note exist -> quit'
            return,-1
        endif
        q1="UPDATE CalibrationMetadata SET releaseDate='"+date_str+"'"
        q1=q1+" WHERE calibrationSetId="+strtrim(string(newCalId),2)
        if keyword_set(verbose) then print,q1
        if insert then result=conn17->execute(q1)
    endif


    ; insert one row at a time
    q1="INSERT INTO dbo.SimPrismDegColumnCalTable(calibrationSetId, x, y, yUncertainty) VALUES("
    n54=n_elements(solar54)
    n55=n_elements(solar55)
    if instrumentModeId eq 31 or instrumentModeId eq 54 then begin
        for i=0L,n_elements(solar54)-1 do begin
            print,i+1,double(i+1)/double(n54+1)*100d,format='("   ",i0,"  (",i0," % )",$,%"\r")'
            str = strtrim(string(newCalId,format='(I)'),2)
            str = str+', '+strtrim(string(time54[i]/1d6),2)
            str = str+', '+strtrim(string(solar54[i].solar_exp),2)
            str = str+', 0.0'
            q2 = q1+str+')'
            if keyword_set(verbose) then print,q2
            if insert then result=conn17->execute(q2)
        endfor
        outdata = replicate({calibrationSetId:0L, X:0d, Y:0d}, n_elements(solar54))
        outdata.calibrationSetId = newCalId
        outdata.x = time54/1d6
        outdata.y = solar54.solar_exp
    endif else if instrumentModeId eq 32 or instrumentModeId eq 55 then begin
        for i=0L,n_elements(solar55)-1 do begin
            print,i+1,double(i+1)/double(n55+1)*100d,format='("   ",i0,"  (",i0," % )",$,%"\r")'
            str = strtrim(string(newCalId,format='(I)'),2)
            str = str+', '+strtrim(string(time55[i]/1d6),2)
            str = str+', '+strtrim(string(solar55[i].solar_exp),2)
            str = str+', 0.0'
            q2 = q1+str+')'
            if keyword_set(verbose) then print,q2
            if insert then result=conn17->execute(q2)
        endfor
        outdata = replicate({calibrationSetId:0L, X:0d, Y:0d}, n_elements(solar55))
        outdata.calibrationSetId = newCalId
        outdata.x = time55/1d6
        outdata.y = solar55.solar_exp
    endif
    print,''

    return, outdata

END





