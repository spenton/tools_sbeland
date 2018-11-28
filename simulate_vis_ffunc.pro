;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Simulate some data (starting with day 453 as reference) by adding an
;   exposure rate for A and B and applying our Kappa degradation and a F-function.
;   The returned structure can then be used as input to get_global_degcol.pro .
;
; CALLING SEQUENCE:
;   data = simulate_ffunc()
;
; INPUT PARAMETERS:
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
;   version -
;      Requests a specific version of the SimProfileIntegralCal, otherwise
;      pick the latest version valid for the specified time range.
;   order - 
;      Order of the polynomial correction.
;
; RETURNED PARAMETERS:
;
; OPTIONAL OUTPUT PARAMETERS:
;
; EXAMPLE:  
;
; COMMON BLOCKS:
;      NONE
;
; NOTES:
;    SimCalibratedIrradiance = SimCorrectedIrradiance * aleph
;    SimCorrectedIrradiance  = SimUncorrectedIrradiance / (prismTransDeg * detectorDeg)
;    SimPrismTransDegradation = (1 - rayA) * exp(-kappa*SimDegradationColumn) + rayA * exp(-kappa * SimDegradationColumn * rayB)
;    SimDegradationColumn = SolarExposure * F_Function    (a function of time and wavelength)
;
;       for ESR:
;          detectorDeg=1
;          SimUncorrectedIrradiance = PhaseDetectedDN / SunObserverDistance / ProfileIntegral
;
;       for diodes:
;          detectorDeg = c[0] * exp(-c[1] * (time - c[2]))
;          SimUncorrectedIrradiance = SimConvertedDataNumbers / diodeTempCorr / SunObserverDistance / ProfileIntegral
;          SimConvertedDataNumber = (DN - dark) / diodeADCGain / diodeFeedbackResistance
;
;
; REVISION HISTORY:
;   Revision: $Id: simulate_vis_ffunc.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;------------------------------------------------------------------
;
function simulate_vis_ffunc, version=version, inspectrum=inspectrum, $
         solar54=solar54, solar55=solar55, step55=step55,  order=order, $
         coeff_corr=coeff_corr, kappa=kappa, afact=afact, outdata=outdata, $
         fakedata=fakedata


    modeA=41
    modeB=45
    refT0 = 453.67840d
    refT1 = 453.69524d

    if n_elements(version) eq 0 then  version=2000
    if n_elements(user) eq 0 then     user='sbeland'
    if n_elements(password) eq 0 then password='sbeland1'
    if n_elements(dburl) eq 0 then    dbUrl='jdbc:sybase:Tds:lasp-db-dev:4100/SORCE_SIM_V20'
    if n_elements(dbdriver) eq 0 then dbDriver='com.sybase.jdbc3.jdbc.SybDriver'

    ;if the kappa is not provided, use the default exponent coefficients
    if n_elements(kappa) eq 0 then begin
        ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2002.txt',wave0,kappa0,format='(d,d)'
        readcol,'~/SORCE/data/ESR_kappa_453_1570_2000.txt',wave0,kappa0,format='(d,d)'
        kappa = {x:wave0, y:kappa0}
    endif

    if n_elements(afact) eq 0 then begin
        readcol,'~/SORCE/data/overlap_vis1.txt',wave0,singlePass,format='(d,d)' 
        ;singlePass *= 4.8d
        firstSurfaceDegradation = singlePass*0d + 0.5d
        afact={wavelength:wave0, singlePassAreaFraction:singlePass, firstSurfaceDegradation:firstSurfaceDegradation}
    endif

    ; if no F-Function was provided, use the default polynomial coefficients initially found for SORCE SIM
    ;if n_elements(coeff_corr) eq 0 then coeff_corr=[1.5054297d, -0.00072414176d, 1.8150587d-07]

    ; get the REAL SimSolarExposureData for modes 31 and 32
    if n_elements(solar54) eq 0 then begin
        q1 = 'SELECT orbitStartTimeGps t0, orbitEndTimeGps t1, cumulativeSolarExpHrtOut solar_exp, '+$
             'solarExposureHrtOut as solar_exp_orbit from SimSolarExposureData where instrumentModeId=54'
        query_database, q1, solar54, info
    endif

    if n_elements(solar55) eq 0 then begin
        q1 = 'SELECT orbitStartTimeGps t0, orbitEndTimeGps t1, cumulativeSolarExpHrtOut solar_exp, '+$
             'solarExposureHrtOut as solar_exp_orbit from SimSolarExposureData where instrumentModeId=55'
        query_database, q1, solar55, info
    endif

    get_oneau=0
    ;tnames = tag_names(solar54)
    ;p=where(strpos(tnames,"ONEAU") ge 0, count54)
    ;if count54 eq 0 then get_oneau=1
    ;tnames = tag_names(solar55)
    ;p=where(strpos(tnames,"ONEAU") ge 0, count55)
    ;if count55 eq 0 then get_oneau=1
    
    if get_oneau then begin
        ; get the 1-AU correction
        print,'  getting SUNOBSERVERDISTANCECORRECTION from database ...'
        q1="SELECT microsecondsSinceGpsEpoch, sunObserverDistanceCorrection FROM SolarDistAndDopplerFixedStep2 "+$
           "WHERE instrumentModeId=50 and version=7"
        query_database, q1, solardist, info
        if count54 eq 0 then begin
            oneau54 = interpol(solardist.SUNOBSERVERDISTANCECORRECTION, solardist.MICROSECONDSSINCEGPSEPOCH, solar54.(1))
            append_tag,solar54,'oneau',oneau54,/slim
        endif
        if count54 eq 0 then begin
            oneau55 = interpol(solardist.SUNOBSERVERDISTANCECORRECTION, solardist.MICROSECONDSSINCEGPSEPOCH, solar55.(1))
            append_tag,solar55,'oneau',oneau55,/slim
        endif
        ; apply the 1AU correction (the factor 4 was determine empirically)
        solar54.solar_exp=total(solar54.solar_exp_orbit,/cum);  / (1d +(solar54.oneau-1d)/4d)
        solar55.solar_exp=total(solar55.solar_exp_orbit,/cum);  / (1d +(solar55.oneau-1d)/4d)
    endif

    ; simulate an F-Function to apply to the solar exposure
    corr54=1d
    corr55=1d
    ;corr54 = poly(gps2sd(solar54.t1/1d6),coeff_corr)
    ;corr55 = poly(gps2sd(solar55.t1/1d6),coeff_corr)
    ;cum54 = total(poly(gps2sd(solar54.t1/1d6),coeff_corr),/cum)
    ;cum55 = total(poly(gps2sd(solar55.t1/1d6),coeff_corr),/cum)
    ;mx=max(cum54,min=mn)
    ;corr54 = (cum54 - mn)/(mx-mn)/10d + 1d
    ;mx=max(cum55,min=mn)
    ;corr55 = (cum55 - mn)/(mx-mn)/10d + 1d

    ; try a sinusoidal FFunction for testing
    corr54 = sin(gps2sd(solar54.t1/1d6) / 500d) + 1.0d
    corr55 = interpol(corr54, gps2sd(solar54.t1/1d6), gps2sd(solar55.t1/1d6))

    ;solar54.solar_exp=total(solar54.solar_exp_orbit*corr54,/cum) / (1d +(solar54.oneau-1d)/4d)
    solar54.solar_exp=total(solar54.solar_exp_orbit*corr54,/cum)
    solar54.solar_exp/=86400d
    ;solar55.solar_exp=total(solar55.solar_exp_orbit*corr55,/cum) / (1d +(solar55.oneau-1d)/4d)
    solar55.solar_exp=total(solar55.solar_exp_orbit*corr55,/cum)
    solar55.solar_exp/=86400d


    ; skip the making of the data from the reference spectra
    if keyword_set(fakedata) then begin
        ; get the reference spectra for SimA and SimB
        refSpectA=get_science_product(['Wavelength','SimCorrectedIrradiance'],reft0,reft1,modeA,$
            /mission, version=version, user=user,password=password,dburl=dburl,dbdriver=dbdriver)

        refSpectB=get_science_product(['Wavelength','SimCorrectedIrradiance'],reft0,reft1,modeB,$
            /mission, version=version, user=user,password=password,dburl=dburl,dbdriver=dbdriver)

        if size(refSpectA,/tname) ne 'STRUCT' or  size(refSpectB,/tname) ne 'STRUCT' then begin
            if keyword_set(verbose) then print,'No valid reference spectra found'
            return,-1
        endif

        ; redefine the structure to the format expected by get_global_degcol
        s=sort(refspectA.wavelengthRef)
        q=uniq(refSpectA[s].CORRECTEDPRISMPOSITION)
        refspectA=refspectA[s[q]]
        spectA={timestamp:refSpectA.MICROSECONDSSINCEGPSEPOCH, $
                wavelength:refspectA.wavelengthref, $
                irradiance:refSpectA.irradiance, $
                prismposition:refSpectA.CORRECTEDPRISMPOSITION, $
                CORRECTEDPOSITION:refSpectA.CORRECTEDPRISMPOSITION}

        s=sort(refspectB.wavelengthRef)
        q=uniq(refSpectB[s].CORRECTEDPRISMPOSITION)
        refspectB=refspectB[s[q]]
        spectB={timestamp:refspectB.MICROSECONDSSINCEGPSEPOCH, $
                wavelength:refspectB.wavelengthref, $
                irradiance:refspectB.irradiance, $
                prismposition:refspectB.CORRECTEDPRISMPOSITION, $
                CORRECTEDPOSITION:refspectB.CORRECTEDPRISMPOSITION}

        ; add a solar cycle to simA and simB data with minimum on day 2200
        x=dindgen(1100)/100 ;in years (full cycle)
        y=sin((x+2.25d)*2d * !dpi / 11.0) * 7.32d-4 / 2d + (7.32d-4/2d + 1.0d)  ; minimum aligned around year 6.0
        x*=365.25d       ; x now in ~days
        x=sd2gps(x)*1d6  ; x now in micronSecondsSinceGpsEpoch
        solcycle = interpol(y, x, spectA.timestamp, /spline)
        spectA.irradiance *= solcycle
        solcycle = interpol(y, x, spectB.timestamp, /spline)
        spectB.irradiance *= solcycle

        nptsA = n_elements(spectA)
        nptsB = n_elements(spectB)

        ; process the data over 3200 days 
        ; the first spectra starts on day 0 with both SimA and SimB
        ; SimA is taken every 7 days and SimB every 30 days
        ndays=4300d
        skipdaysA=7d
        skipDaysB=28d
        nspectA = ceil(ndays / skipDaysA)
        nspectB = ceil(ndays / skipDaysB)

        ; make data start at actual first time prism was exposed
        p=where(solar54.solar_exp_orbit gt 0d)
        t0A = min(solar54[p[0]].t1)
        spectA.timestamp -= (min(spectA.timestamp) - t0A)
        spectB.timestamp -= (min(spectB.timestamp) - t0A)

        sp41 = replicate(spectA[0], nptsA * nspectA)
        sp45 = replicate(spectB[0], nptsB * nspectB)

        ; get the kappa corresponding to each wavelengths for SimA and SimB
        kappa_a = interpol(kappa.y, kappa.x, spectA.wavelength, /spline)
        kappa_b = interpol(kappa.y, kappa.x, spectB.wavelength, /spline)

        afact_a = interpol(afact.singlePassAreaFraction, afact.wavelength, spectA.wavelength, /spline)
        afact_b = interpol(afact.singlePassAreaFraction, afact.wavelength, spectB.wavelength, /spline)

        for i=0,nspectA-1L do begin
            sp41[i*nptsA:(i+1)*nptsA-1L] = spectA
            sp41[i*nptsA:(i+1)*nptsA-1L].timestamp += double(i) * skipDaysA *86400d6
            degCol = interpol(solar54.solar_exp, solar54.t1, sp41[i*nptsA:(i+1)*nptsA-1L].timestamp)
            prismtransdegradation = (1d - afact_a) * exp ( -kappa_a *  degCol) + afact_a * exp(-kappa_a * degCol / 2d)
            sp41[i*nptsA:(i+1)*nptsA-1L].irradiance *= prismtransdegradation
        endfor

        for i=0L,nspectB-1L do begin
            sp45[i*nptsB:(i+1)*nptsB-1L] = spectB
            sp45[i*nptsB:(i+1)*nptsB-1L].timestamp += double(i) * skipDaysB *86400d6
            degCol = interpol(solar55.solar_exp, solar55.t1, sp45[i*nptsB:(i+1)*nptsB-1L].timestamp)
            prismtransdegradation = (1d - afact_b) * exp ( -kappa_b *  degCol) + afact_b * exp(-kappa_b * degCol / 2d)
            sp45[i*nptsB:(i+1)*nptsB-1L].irradiance *= prismtransdegradation
        endfor

        inspectrum = {simA:temporary(sp41), simB:temporary(sp45)}

    endif else begin
        restore,'~/SORCE/data/sima_vis_2000_calibrated.sav'
        restore,'~/SORCE/data/simb_vis_2000_calibrated.sav'

        ; evaluate a prismDegaradtion based on Kappa, Afact and FFunc and convert to UnCorrectedIrradiance
        for i=0L,n_elements(visa_2000.spect20)-1 do begin
            kappa_a = interpol(kappa.y, kappa.x, visa_2000.spect20[i].wavelength, /spline)
            afact_a = interpol(afact.singlePassAreaFraction, afact.wavelength, visa_2000.spect20[i].wavelength, /spline)
            degCol = interpol(solar54.solar_exp, solar54.t1, visa_2000.spect20[i].timestamp)
            prismtransdegradation = (1d - afact_a) * exp ( -kappa_a *  degCol) + afact_a * exp(-kappa_a * degCol / 2d)
            visa_2000.spect20[i].irradiance *= prismtransdegradation
        endfor

        for i=0L,n_elements(visb_2000.spect20)-1 do begin
            kappa_b = interpol(kappa.y, kappa.x, visb_2000.spect20[i].wavelength, /spline)
            afact_b = interpol(afact.singlePassAreaFraction, afact.wavelength, visb_2000.spect20[i].wavelength, /spline)
            degCol = interpol(solar55.solar_exp, solar55.t1, visb_2000.spect20[i].timestamp)
            prismtransdegradation = (1d - afact_b) * exp ( -kappa_b *  degCol) + afact_b * exp(-kappa_b * degCol / 2d)
            visb_2000.spect20[i].irradiance *= prismtransdegradation
        endfor

        inspectrum = {simA:temporary(visa_2000.spect20), simB:temporary(visb_2000.spect20)}

    endelse

    ; now use the get_global_degcol to see if we can retrieve the input F-Function above
    coeffs=get_global_degcol(0,4300,/mission,inspectrum=inspectrum, /visdiode, solar54=solar54, solar55=solar55, $
        kappa=kappa, outdata=outdata, step55=step55, order=order)

    return,coeffs

end


