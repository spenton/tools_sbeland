;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Fit a spectra to the SORCE SIM instrument model sectra by using the
;   AMOEBA routine to minimize the difference with the reference when 
;   adjusting GammaZ and PixelScale.  The routine
;   finds the coefficients of the best fit and returns the original spectra 
;   with the new wavelength scale.
;
; CALLING SEQUENCE:
;   spectra = ALIGN_SPECTRA3(t0, t1, instrumentModeId, /missionDays, refT0=t0, refT1=t1)
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
;      The instrument mode of interest are currently limited to:
;      41	SIM_A	VIS1
;      43	SIM_A	UV
;      44	SIM_A	IR
;      45	SIM_B	VIS1
;      47	SIM_B	UV
;      31	SIM_A	ESR
;      32	SIM_B	ESR
;
; OPTIONAL INPUT PARAMETERS:
;   spect - 
;      Structure containing the spectra to align. If provided, the
;      startTime, stopTime and instrumentModeId are ignored.
;   refSpect - 
;      Structure containg the reference spectra to align to.
;   refT0 -
;      StartTime of the reference spectra in the same units as T0.
;      If not specified, uses our "standard" reference time range.
;   refT1 -
;      StopTime of the reference spectra in the same units as T1.
;      If not specified, uses our "standard" reference time range.
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
;   quiet -
;      Do not print extra messages
;   noderiv -
;      By default, the program fits the derivative of the spectrum. This flag
;      will let the program use the original spectrum.
;
; RETURNED PARAMETERS:
;   A structure with the spectra aligned to the reference spectra in wavelength.
;
; OPTIONAL OUTPUT PARAMETERS:
;   coeffs -
;      The coefficients of the best polynomial fit.
;   status -
;      Returned STATUS of the mpfitfun routine.
;
; EXAMPLE:  
;   Adjust the wavelength scale of the SimB UV diode to the SimA UV diode for 
;   the SolarQuickScan24 on day 453.67:
;
;   IDL> T0=453.67840d
;   IDL> T1=453.69524d
;   IDL> dbtables=['SimProfileIntegral','SimCalibratedIrradiance']
;   IDL> res43=get_science_product(dbtables, T0, T1, 43, /mission, version=17)
;   IDL> res=align_spectra3(T0,T1,47,/mission,refspect=res43,coeffs=coeffs,status=status)
;   IDL> help,status
;   IDL> STATUS          INT       =        2
;   IDL> help,res
;   ** Structure <3877418>, 2 tags, length=18000, data length=18000, refs=1:
;      WAVELENGTH      DOUBLE    Array[1125]
;      IRRADIANCE      DOUBLE    Array[1125]
;
; COMMON BLOCKS:
;      NONE
;
; NOTES:
;  You may want to run get_sorce_plan before hand to make sure an activity
;  with observations in the desired mode was performed during the time
;  span of interest.
;
;
;   STILL A WORK IN PROGRESS:  WORKS FOR SOME DATA BUT NOT WITH OTHER
;   -----------------------------------------------------------------
;
; REVISION HISTORY:
;   Revision: $Id: align_spectra3.pro,v 1.2 2018/11/26 16:16:02 spenton Exp $
;-
;*******************************************************************
function my_amoeba3, coeffs
    common myamoeba3_common, sp, rsp, no_deriv, wave_range, use_rms, imode, dettempcorr

    if n_elements(coeffs) eq 1 then begin
        newx = ccd2lambda(iMode, sp.x, sp.prismtemp, gamz=coeffs[0])
    endif else begin
        newx = ccd2lambda(iMode, sp.x, sp.prismtemp, gamz=coeffs[0], pixsize=coeffs[1])
    endelse
    newy=interpol(sp.y,newx,rsp.x, /lsq)
    if size(dettempcorr,/tname) eq 'STRUCT' then begin
        ; detCorr is already resampled to refSpect wavelength space
        newy *= (1d + dettempcorr.tempcorr*(dettempcorr.refTemp - sp.detectorTemp))
    endif

    if no_deriv eq 0 then newy = deriv(rsp.x,newy)
    pos = where(rsp.x ge wave_range[0] and rsp.x le wave_range[1], count)
    if count eq 0 then return,0d

    if use_rms eq 1 then begin
        corr_value = sqrt(total((newy[pos] - rsp.y[pos])^2d)/count)
    endif else begin
        corr_value = abs(1d - abs(correlate(newy[pos], rsp.y[pos], /double))) * 1d10
    endelse

    ;print,coeffs, corr_value/1d10,format='(2(G,","),G)'
    return, corr_value

end
;*******************************************************************

function align_spectra3, startTime, stopTime, instrumentModeId, $
         gps=gps, missionDays=missionDays, julianDays=julianDays, $
         refT0=refT0, refT1=refT1, refSpect=refSpect, spect=spect, $
         version=version, quiet=quiet, coeffs=coeffs, $
         noderiv=noderiv, status=status, profile_data=profile_data, $
		 verbose=verbose, goodness=goodness, wrange=wrange, rms=rms, $
         ccdfit=ccdfit, no_plot=no_plot, no_dettemp=no_dettemp, $
         detCorr=detCorr, order=order, solardist=solardist

    common myamoeba3_common
    common SIM_DEF, sim_model

    if n_elements(version) eq 0 then version=20
    if keyword_set(noderiv) then no_deriv=1 else no_deriv=0
    if keyword_set(rms) then use_rms=1 else use_rms=0
	if n_elements(order) eq 0 then order=3
    if n_elements(wrange) eq 0 then wave_range=[200d,2400d] else wave_range=wrange
    goodness=0d
    coeffs=0d

    maxiter = 5000d

    ; check if the common block has been defined
    if n_elements(sim_model) eq 0 then get_sim_model

    iMode_pos=where(sim_model.imode eq instrumentModeId, count)
    if count eq 0 then return,-1 else iMode=instrumentModeId

    ; if a reference spectra was provided use that instead
    if size(refSpect,/tname) ne 'STRUCT' then begin
        if n_elements(refT0) eq 0 or n_elements(refT1) eq 0 then begin
            ; this is our default REFERENCE SolarQuickScan24 spectra 
            if instrumentModeId eq 31 or instrumentModeId eq 32 then begin
                ; full ESR scan
                refT0 = 453.02065d
                refT1 = 453.98991d
            endif else begin
                refT0 = 453.67840d
                refT1 = 453.69524d
            endelse
        endif
        if keyword_set(missionDays) then begin
           ; do nothing here since already in mission days
        endif else if keyword_set(julianDays) then begin
           ; user specified time in julian days
           refT0 = sd2jd(refT0)
           refT1 = sd2jd(refT1)
        endif else begin
           ; user specified time in gps days (default)
           refT0 = sd2gps(refT0)*1.d6
           refT1 = sd2gps(refT1)*1.d6
        endelse

        ; get the reference spectra
        if instrumentModeId gt 32 then begin
            ; use get_sim_spectra to get the diode data
            refSpect = get_sim_spectra(refT0, refT1, instrumentModeId, gps=gps, missionDays=missionDays, $
                julianDays=julianDays, version=version, /noCcdCorr, /dn_only, profile_data=profile_data, solardist=solardist)
        endif else begin
            ; extract the ESR data from the database for version 19
            refSpect = get_esr_spectra(refT0, refT1, instrumentModeId, gps=gps, missionDays=missionDays, $
                julianDays=julianDays, version=version, /noCcdCorr, database=database, /no_duplicate, solardist=solardist)
        endelse

        if size(refSpect,/tname) ne 'STRUCT' then begin
            if keyword_set(verbose) then print,'No valid reference spectra found'
            return,-1
        endif
        s=sort(refSpect.ccdpos)
        refSpect=refSpect[s]
        ; force the wavelength of the reference spectra to our model at the reference temperature
        refSpect.wavelength=ccd2lambda(instrumentModeId,refSpect.ccdpos, TREF_SILICA)
    endif

    ; get the spectra to fit to the reference
    if size(spect,/tname) ne 'STRUCT' then begin
        if instrumentModeId gt 32 then begin
            spect = get_sim_spectra(startTime, stopTime, instrumentModeId, gps=gps, missionDays=missionDays, $
                julianDays=julianDays, version=version, /noCcdCorr, /noTempCorr, /dn_only, $
                solarDist=solarDist,profile_data=profile_data)
        endif else begin
            spect = get_esr_spectra(startTime, stopTime, instrumentModeId, gps=gps, missionDays=missionDays, $
                julianDays=julianDays, version=version, /noCcdCorr, database=database, solarDist=solarDist,/no_duplicate)
        endelse
        if size(spect,/tname) ne 'STRUCT' then begin
            if keyword_set(verbose) then print,'No valid spectra found'
            return,-1
        endif
        s=sort(spect.ccdpos)
        spect=spect[s]
    endif

    ;refSpect.dn_tempcorr = refSpect.dn
    ;spect.dn_tempcorr = spect.dn

    ; match the same ccdpos between the 2 scans
    match, refSpect.ccdpos, spect.ccdpos, suba, subb
    if n_elements(suba) lt 10 then begin
        print,'Not enough matching points'
        return,-1
    endif

    if keyword_set(no_dettemp) then begin
        detCorr=-1
    endif else if (instrumentModeId eq 41 or instrumentModeId eq 45) and size(detCorr,/tname) ne 'STRUCT' then begin
        infile = file_search("$SORCE_DATA","sima_vis_tempcorr_20.txt",/expand_env,/full)
        ; the following file is the results of lab measurements in Feb/2014
        ;infile = file_search("$SORCE_DATA","SORCE_Si_drdt_fit.txt",/expand_env,/full)
        readcol,infile[0],wv,c0,tref,format='(D,D,D)',/silent
        corr=interpol(c0, wv, refSpect.wavelength, /lsq)
        refSpect.dn_tempcorr = refSpect.dn * (1d + corr*(tref[0] - refSpect.detectorTemp))
        detCorr={wavelength:refSpect[suba].wavelength, tempcorr:corr[suba], refTemp:tref[0]}
    endif else if (instrumentModeId eq 44 or instrumentModeId eq 48) and size(detCorr,/tname) ne 'STRUCT' then begin
        infile = file_search("$SORCE_DATA","sima_ir_tempcorr_21.txt",/expand_env,/full)
        ; the following file is the results of lab measurements in Feb/2014
        ;infile = file_search("$SORCE_DATA","SORCE_InGaAs_drdt_fit.txt",/expand_env,/full)
        readcol,infile[0],wv,c0,tref,format='(D,D,D)',/silent
        corr=interpol(c0, wv, refSpect.wavelength, /lsq)
        refSpect.dn_tempcorr = refSpect.dn * (1d + corr*(tref[0] - refSpect.detectorTemp))
        detCorr={wavelength:refSpect[suba].wavelength, tempcorr:corr[suba], refTemp:tref[0]}
    endif else detCorr=-1
    dettempcorr=detCorr


    if size(detCorr,/tname) eq 'STRUCT' then begin
        rsp={x:refSpect[suba].wavelength, y:refSpect[suba].dn_tempcorr}
        sp={x:spect[subb].ccdpos, y:spect[subb].dn, prismTemp:spect[subb].prismTemp, detectorTemp:spect[subb].detectorTemp}
    endif else begin
        rsp={x:refSpect[suba].wavelength, y:refSpect[suba].dn}
        sp={x:spect[subb].ccdpos, y:spect[subb].dn, prismTemp:spect[subb].prismTemp}
    endelse

    if no_deriv eq 0 then rsp.y=deriv(rsp.x, rsp.y)

    init_coef=[sim_model.gamz[iMode_pos], sim_model.pix[iMode_pos]]
    coeffs = AMOEBA(1.0d-8, P0=init_coef, scale=[1d-5,1d-5], FUNCTION_VALUE=fval, FUNCTION_NAME='my_amoeba3', nmax=maxiter)

    if use_rms eq 0 then fval /= 1d10
    goodness= 1d - fval
    if keyword_set(verbose) then print,'Final coefficients = ',coeffs
    if keyword_set(verbose) then print,'Fit goodness=',goodness

    if n_elements(coeffs) eq 1 and coeffs[0] eq -1 then begin
        print,'AMOEBA failed to converge'
        return,-1
    endif

    outspect=spect
    outspect.wavelength = ccd2lambda(iMode, outspect.ccdpos_uncorr, outspect.prismtemp, $
        gamz=coeffs[0], pixsize=coeffs[1]) 

    if size(detCorr,/tname) eq 'STRUCT' then begin
        corr=interpol(c0, wv, outspect.wavelength, /lsq)
        outspect.dn_tempcorr = outspect.dn * (1d + corr*(Tref[0] - outspect.detectorTemp))
    endif

    ; get the delta_ccdpos vs ccdpos (using a 3rd order polynomial fit)
    new_ccdpos = lambda2ccd(iMode, outspect.wavelength, outspect.prismtemp)
    delta_ccdpos = new_ccdpos - outspect.ccdpos_uncorr
    if order lt 2 then begin
        ccdfit = poly_fit(outspect.ccdpos_uncorr, delta_ccdpos, order, /double)
    endif else begin
        ccdfit = robust_poly_fit(outspect.ccdpos_uncorr, delta_ccdpos, order, /double)
    endelse
    if NOT keyword_set(no_plot) then begin
        plot_multi,outspect.ccdpos_uncorr, delta_ccdpos, outspect.ccdpos_uncorr, poly(outspect.ccdpos_uncorr, ccdfit), $
            /xst,/yst,xtitle='CCD Position (subpixels)',ytitle='Delta CCDPOS',title='Mode='+strtrim(string(instrumentModeId),2)+$
            ' Align_Spectra2 on '+strtrim(string(gps2sd(outspect[0].timetag/1d6),format='(F7.2)'),2), psym=[4,-3], $
            label=['Delta CCD Pos','Polynomial Fit'], charsize=1.4, thickness=[1,3]
    endif

	; define the goodness of the fit by the correlate
	goodness=0d
    ; find the overlap in wavelengths between refSpect and outspect
    refMin = min(refSpect.wavelength,max=refMax)
    spMin  = min(outspect.wavelength,max=spMax)
    mnwave=max([refMin,spMin,wave_range[0]])
    mxwave=min([refMax,spMax,wave_range[1]])
    refp=where(refSpect.wavelength ge mnWave and refSpect.wavelength le mxWave,refCount)
    if refCount gt 0 then begin
        if size(detCorr,/tname) eq 'STRUCT' then begin
            newy = interpol(outspect.dn_tempcorr, outspect.wavelength, refspect[refp].wavelength, /lsq)
            goodness=correlate(deriv(refspect[refp].wavelength,newy), deriv(refspect[refp].wavelength,refspect[refp].dn_tempcorr))
        endif else begin
            newy = interpol(outspect.dn, outspect.wavelength, refspect[refp].wavelength, /lsq)
            goodness=correlate(deriv(refspect[refp].wavelength,newy), deriv(refspect[refp].wavelength,refspect[refp].dn))
        endelse
	endif
 
    return, outspect

end



