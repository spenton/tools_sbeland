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
;   Revision: $Id: simulate_ffunc.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;------------------------------------------------------------------
;
function simulate_ffunc, modeA, version=version, inspectrum=inspectrum, $
         solar54=solar54, solar55=solar55, step55=step55,  order=order, $
         coeff_corr=coeff_corr, kappa_coeff=kappa_coeff, outdata=outdata


    if modeA eq 31 then begin
        ; full ESR scan
        refT0 = 469.21836d
        refT1 = 469.78784d
        modeB=32
    endif else begin
        refT0 = 453.67840d
        refT1 = 453.69524d
        modeB=modeA+4
    endelse

    if n_elements(version) eq 0 then  version=95
    if n_elements(user) eq 0 then     user='sbeland'
    if n_elements(password) eq 0 then password='sbeland1'
    if n_elements(dburl) eq 0 then    dbUrl='jdbc:sybase:Tds:lasp-db-dev:4100/SORCE_SIM_V19'
    if n_elements(dbdriver) eq 0 then dbDriver='com.sybase.jdbc3.jdbc.SybDriver'

    ;if the kappa is not provided, use the default exponent coefficients
    if n_elements(kappa_coeff) lt 3 then kappa_coeff=[0.0096235817d,   -0.0071923836d,   1.7092305d-05]

    ; if no F-Function was provided, use the default polynomial coefficients initially found for SORCE SIM
    if n_elements(coeff_corr) eq 0 then coeff_corr=[1.5054297d, -0.00072414176d, 1.8150587d-07]

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
    tnames = tag_names(solar54)
    p=where(strpos(tnames,"ONEAU") ge 0, count54)
    if count54 eq 0 then get_oneau=1
    tnames = tag_names(solar55)
    p=where(strpos(tnames,"ONEAU") ge 0, count55)
    if count55 eq 0 then get_oneau=1
    
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
        solar54.solar_exp=total(solar54.solar_exp_orbit,/cum) / (1d +(solar54.oneau-1d)/4d)
        solar55.solar_exp=total(solar55.solar_exp_orbit,/cum) / (1d +(solar55.oneau-1d)/4d)
    endif

    ; get the reference spectra for SimA and SimB
    refSpectA=get_science_product(['Wavelength','SimCorrectedIrradiance','PhaseDetectedDN','SimPrismTransDegradation'],reft0,reft1,modeA,$
        /mission, version=version, user=user,password=password,dburl=dburl,dbdriver=dbdriver)

    refSpectB=get_science_product(['Wavelength','SimCorrectedIrradiance','PhaseDetectedDN','SimPrismTransDegradation'],reft0,reft1,modeB,$
        /mission, version=version, user=user,password=password,dburl=dburl,dbdriver=dbdriver)

    if size(refSpectA,/tname) ne 'STRUCT' or  size(refSpectB,/tname) ne 'STRUCT' then begin
        if keyword_set(verbose) then print,'No valid reference spectra found'
        return,-1
    endif

    ; remove the prismTransDeg from the correctedIrradiance
    refSpectA.irradiance *= refSpectA.prismtransdegradation
    s=sort(refspectA.wavelengthRef)
    q=uniq(refSpectA[s].prismposition)
    refspectA=refspectA[s[q]]

    s=sort(refspectB.wavelengthRef)
    q=uniq(refSpectB[s].prismposition)
    refspectB=refspectB[s[q]]

    ; modify the irradiance of refSpectB to be the same as A (with the corresponding wavelengths)
    ; this will be our day 0 data
    refSpectB.irradiance = interpol(refSpectA.irradiance, refSpectA.wavelengthref, refSpectB.wavelengthref, /spline)

    nptsA = n_elements(refSpectA)
    nptsB = n_elements(refSpectB)

    ; process the data over 3200 days 
    ; the first spectra starts on day 0 with both SimA and SimB
    ; SimA is taken every 7 days and SimB every 30 days
    ndays=3200d
    skipdaysA=7d
    skipDaysB=28d
    nspectA = ceil(ndays / skipDaysA)
    nspectB = ceil(ndays / skipDaysB)

    ; make data start at actual first time prism was exposed
    p=where(solar54.solar_exp_orbit gt 0d)
    t0A = min(solar54[p[0]].t1)
    refSpectA.microsecondsSinceGPSepoch -= (min(refSpectA.microsecondsSinceGPSepoch) - t0A)
    refSpectB.microsecondsSinceGPSepoch -= (min(refSpectB.microsecondsSinceGPSepoch) - t0A)

    sp31 = replicate(refSpectA[0], nptsA * nspectA)
    sp32 = replicate(refSpectB[0], nptsB * nspectB)

    ; get the kappa corresponding to each wavelengths for SimA and SimB
    kappa_a = kappa_coeff[0] * exp(kappa_coeff[1] * refSpectA.wavelengthref) + kappa_coeff[2]
    kappa_b = kappa_coeff[0] * exp(kappa_coeff[1] * refSpectB.wavelengthref) + kappa_coeff[2]

    corr54 = poly(gps2sd(solar54.t1/1d6),coeff_corr)
    corr55 = poly(gps2sd(solar55.t1/1d6),coeff_corr)
    cum54 = total(poly(gps2sd(solar54.t1/1d6),coeff_corr),/cum)
    cum55 = total(poly(gps2sd(solar55.t1/1d6),coeff_corr),/cum)
    mx=max(cum54,min=mn)
    corr54 = (cum54 - mn)/(mx-mn)/10d + 1d
    mx=max(cum55,min=mn)
    corr55 = (cum55 - mn)/(mx-mn)/10d + 1d

    ; try a sinusoidal FFunction for testing
    ;corr54 = sin(gps2sd(solar54.t1/1d6) / 500d) + 1.5d
    ;corr55 = sin(gps2sd(solar55.t1/1d6) / 500d) + 1.5d


    ;solar54.solar_exp=total(solar54.solar_exp_orbit*corr54,/cum) / (1d +(solar54.oneau-1d)/4d)
    solar54.solar_exp=total(solar54.solar_exp_orbit,/cum) * corr54
    solar54.solar_exp/=86400d
    ;solar55.solar_exp=total(solar55.solar_exp_orbit*corr55,/cum) / (1d +(solar55.oneau-1d)/4d)
    solar55.solar_exp=total(solar55.solar_exp_orbit,/cum) * corr55
    solar55.solar_exp/=86400d

    for i=0,nspectA-1L do begin
        sp31[i*nptsA:(i+1)*nptsA-1L] = refSpectA
        sp31[i*nptsA:(i+1)*nptsA-1L].microsecondsSinceGPSepoch += double(i) * skipDaysA *86400d6
        degCol = interpol(solar54.solar_exp, solar54.t1, sp31[i*nptsA:(i+1)*nptsA-1L].microsecondsSinceGPSepoch)
        sp31[i*nptsA:(i+1)*nptsA-1L].prismtransdegradation = exp ( -kappa_a *  degCol)
        sp31[i*nptsA:(i+1)*nptsA-1L].irradiance *= sp31[i*nptsA:(i+1)*nptsA-1L].prismtransdegradation
    endfor

    for i=0L,nspectB-1L do begin
        sp32[i*nptsB:(i+1)*nptsB-1L] = refSpectB
        sp32[i*nptsB:(i+1)*nptsB-1L].microsecondsSinceGPSepoch += double(i) * skipDaysB *86400d6
        degCol = interpol(solar55.solar_exp, solar55.t1, sp32[i*nptsB:(i+1)*nptsB-1L].microsecondsSinceGPSepoch)
        sp32[i*nptsB:(i+1)*nptsB-1L].prismtransdegradation = exp ( -kappa_b *  degCol)
        sp32[i*nptsB:(i+1)*nptsB-1L].irradiance *= sp32[i*nptsB:(i+1)*nptsB-1L].prismtransdegradation
    endfor

    inspectrum = {sp31:temporary(sp31), sp32:temporary(sp32)}

    ; now use the get_global_degcol to see if we can retrieve the input F-Function above
    coeffs=get_global_degcol(0,3200,/mission,inspectrum=inspectrum, solar54=solar54, solar55=solar55, $
        kappa=kappa_coeff, outdata=outdata, step55=step55, order=order)

    return,coeffs

end


