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
;   Revision: $Id: simulate_vis_ffunc2.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;------------------------------------------------------------------
;
function simulate_vis_ffunc2, wavelength, inspectrum=inspectrum, step=step, kappa=kappa, $
    afact=afact, solar54=solar54, solar55=solar55, order=order


    if n_elements(step) eq 0 then step=3

    ;if the kappa is not provided, use the default exponent coefficients
    if n_elements(kappa) eq 0 then begin
        ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2002.txt',wave0,kappa0,format='(d,d)'
        readcol,'~/SORCE/data/ESR_kappa_453_1570_2011.txt',wave0,kappa0,format='(d,d)'
        s=sort(wave0)
        q=uniq(wave0[s])
        kappa = {x:wave0[s[q]], y:kappa0[s[q]]}
    endif

    if n_elements(afact) eq 0 then begin
        readcol,'~/SORCE/data/overlap_vis1_smooth.txt',wave0,singlePass,format='(d,d)' 
        s=sort(wave0)
        q=uniq(wave0[s])
        firstSurfaceDegradation = singlePass*0d + 0.5d
        afact={wavelength:wave0[s[q]], a:singlePass[s[q]], $
            singlePassAreaFraction:singlePass[s[q]], firstSurfaceDegradation:firstSurfaceDegradation[s[q]]}
    endif

    ; if no F-Function was provided, use the default polynomial coefficients initially found for SORCE SIM
    ;if n_elements(coeff_corr) eq 0 then coeff_corr=[1.5054297d, -0.00072414176d, 1.8150587d-07]

    ; reload the solar exposure record every time since we're modifying it here
    if size(solar54,/tname) ne 'STRUCT' or size(solar55,/tname) ne 'STRUCT' then begin
        restore,'~/SORCE/data/everyOrbitExpos_AB.sav'
    endif

    ; simulate an F-Function to apply to the solar exposure
    corr54=replicate(1d, n_elements(solar54.t1))
    corr55=replicate(1d, n_elements(solar55.t1))
    ;corr54 = poly(gps2sd(solar54.t1/1d6),coeff_corr)
    ;corr55 = poly(gps2sd(solar55.t1/1d6),coeff_corr)
    ;cum54 = total(poly(gps2sd(solar54.t1/1d6),coeff_corr),/cum)
    ;cum55 = total(poly(gps2sd(solar55.t1/1d6),coeff_corr),/cum)
    ;mx=max(cum54,min=mn)
    ;corr54 = (cum54 - mn)/(mx-mn)/10d + 1d
    ;mx=max(cum55,min=mn)
    ;corr55 = (cum55 - mn)/(mx-mn)/10d + 1d

    ; try a sinusoidal FFunction for testing
    corr54 = sin(gps2sd(solar54.t1/1d6) / 250d)*0.25d + 1.0d
    corr55 = interpol(corr54, gps2sd(solar54.t1/1d6), gps2sd(solar55.t1/1d6))

    solar54.solar_exp=total(solar54.solar_exp_orbit*corr54,/cum) / 86400d
    solar55.solar_exp=total(solar55.solar_exp_orbit*corr55,/cum) / 86400d


    if size(inspectrum,/tname) ne 'STRUCT' then begin
        ; simulate the data using the a starting uncorrected spectra
        restore,'~/SORCE/data/sima_vis_uncorr_2011.sav'           

        n54 = n_elements(solar54)-1L
        n55 = n_elements(solar55)-1L

        for i=0L,n_elements(visa_uncorr_2011.spect20)-1L do begin
            good = where(visa_uncorr_2011.spect20[i].wavelength gt 0d)
            s=sort(visa_uncorr_2011.spect20[i].wavelength[good])
            q=uniq(visa_uncorr_2011.spect20[i].wavelength[good[s]])
            good=good[s[q]]
            kk = interpol(kappa.y, kappa.x, visa_uncorr_2011.spect20[i].wavelength[good])
            aa = interpol(afact.singlePassAreaFraction, afact.wavelength, visa_uncorr_2011.spect20[0].wavelength[good])
            if i eq 0 then begin
                refa={wavelength:visa_uncorr_2011.spect20[i].wavelength[good], $
                    irradiance:visa_uncorr_2011.spect20[i].irradiance[good], $
                    timestamp:visa_uncorr_2011.spect20[i].timestamp[good]}
                pos = where(solar54.t0 lt refa.timestamp[0], count)
                initialSolarExp = solar54[pos[-1]].solar_exp
                ; apply the oneau to the raypath
                if where(strpos(tag_names(solar54), 'ONEAU') eq 0) ge 0 then aa *= (solar54[pos[-1]+1].oneau)^2
                prismDeg = (1d - aa) * exp(-kk * initialSolarExp) + aa * exp(-kk * initialSolarExp / 2d)
                ; remove the degradation from our reference spectra as if it had nor degrad
                refa.irradiance /= prismDeg
            endif 
            pos = where(solar54.t0 lt visa_uncorr_2011.spect20[i].timestamp[good[0]], count)
            solarExp = solar54[pos[-1]].solar_exp
            ; apply the oneau to the raypath
            pp = (pos[-1]+1) < n54
            if where(strpos(tag_names(solar54), 'ONEAU') eq 0) ge 0 then aa *= (solar54[pp].oneau)^2
            prismDeg = (1d - aa) * exp(-kk * solarExp) + aa * exp(-kk * solarExp / 2d)
            visa_uncorr_2011.spect20[i].irradiance[good] = spline(refa.wavelength, refa.irradiance, $
                visa_uncorr_2011.spect20[i].wavelength[good], 1d)

            ;visa_uncorr_2011.spect20[i].irradiance[good] = interpol(refa.irradiance, refa.wavelength, $
            ;    visa_uncorr_2011.spect20[i].wavelength[good])
            ; apply the corresponding degradation
            visa_uncorr_2011.spect20[i].irradiance[good] *= prismDeg
        endfor
        
        ; simulate the data using the a starting uncorrected spectra
        restore,'~/SORCE/data/simb_vis_uncorr_2011.sav'           

        for i=0L,n_elements(VISB_uncorr_2011.spect20)-1L do begin
            good = where(VISB_uncorr_2011.spect20[i].wavelength gt 0d)
            s=sort(visb_uncorr_2011.spect20[i].wavelength[good])
            q=uniq(visa_uncorr_2011.spect20[i].wavelength[good[s]])
            good=good[s[q]]
            kk = interpol(kappa.y, kappa.x, VISB_uncorr_2011.spect20[i].wavelength[good])
            aa = interpol(afact.singlePassAreaFraction, afact.wavelength, VISB_uncorr_2011.spect20[0].wavelength[good])
            if i eq 0 then begin
                refb={wavelength:VISB_uncorr_2011.spect20[i].wavelength[good], $
                    irradiance:VISB_uncorr_2011.spect20[i].irradiance[good], $
                    timestamp:VISB_uncorr_2011.spect20[i].timestamp[good]}
                ; use refa spectra for refb sampled at different wavelengths
                refb.irradiance = interpol(refa.irradiance, refa.wavelength, refb.wavelength)
            endif 
            pos = where(solar55.t0 lt VISB_uncorr_2011.spect20[i].timestamp[good[0]], count)
            solarExp = solar55[pos[-1]].solar_exp
            ; apply the oneau to the raypath
            pp = (pos[-1]+1) < n55
            if where(strpos(tag_names(solar55), 'ONEAU') eq 0) ge 0 then aa *= (solar55[pp].oneau)^2
            prismDeg = (1d - aa) * exp(-kk * solarExp) + aa * exp(-kk * solarExp / 2d)
            VISB_uncorr_2011.spect20[i].irradiance[good] = spline(refb.wavelength, refb.irradiance, $
                VISB_uncorr_2011.spect20[i].wavelength[good], 1d)
            ;VISB_uncorr_2011.spect20[i].irradiance[good] = interpol(refb.irradiance, refb.wavelength, $
            ;    VISB_uncorr_2011.spect20[i].wavelength[good])
            ; apply the corresponding degradation
            VISB_uncorr_2011.spect20[i].irradiance[good] *= prismDeg
        endfor

        inspectrum={sima:visa_uncorr_2011.spect20, simb:visb_uncorr_2011.spect20}

    endif

    for w=0L, n_elements(wavelength)-1L do begin
        res=tnmin_pdtau2(wavelength[w], inspectrum=inspectrum, coeffs=coeffs, solar54=solar54, solar55=solar55, $
            step=step, kappa=kappa, afact=afact, /no_detdeg, order=order, /noplot)

        tmp=label_date(date_format='%M %Y')
        trange=gps2sd(([min(res.sima.MICROSECONDSSINCEGPSEPOCH,max=mx),mx])/1d6)
        allsd=dindgen(ceil(trange[1])-floor(trange[0]))+floor(trange[0])
        wv=strtrim(string(wavelength,format='(F0.2)'),2)+'nm'

        plot_multi, res.fdeg_sd,res.fdeg, allsd,poly(allsd,res.coeffs),gps2sd(solar54.t1/1d6),corr54, $
            /xst,/yst, psym=[-4,-3,-3],thick=[1.0,4.0,4.0],label=['FDeg','Polynomial Fit','Injected FFunc'], $
            xtitle='Mission Days',ytitle='F*Kappa',title='Diode F*Kappa @ '+wv, charsize=1.4, yrange=[0d,2d]

        delta2=median(res.sima.new_irrad)-median(res.simb.new_irrad)
        sda=gps2sd(res.sima.(0)/1d6)
        sdb=gps2sd(res.simb.(0)/1d6)
        lineplot,sda,res.sima.irradiance, title='VIS A UnCorrected @ '+wv, psym=-3, nsum=2
        lineplot,sda,res.sima.new_irrad, title='VIS A Corrected FFunc poly @ '+wv, psym=-3, nsum=2
        lineplot,sdb,res.simb.irradiance, title='VIS B UnCorrected @ '+wv, psym=-3, nsum=2
        lineplot,sdb,res.simb.new_irrad+delta2, title='VIS B Corrected FFunc poly @ '+wv, psym=-4, /unzoom

    endfor


    return,coeffs

end


