;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Aligned the specified spectra to the Kurucz reference spectra resampled
;   and convolved with the SIM instrument profile.
;
; CALLING SEQUENCE:
;   spectra = MPFIT_SPECTRA(t0, t1, instrumentModeId, /missionDays, reffile='convolvechkur_uv.sav')
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
;   reffile - 
;      String containing the path of the file with the appropriate 
;      Kurucz spectra to align to.
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
;   quiet -
;      Do not print extra messages
;
; RETURNED PARAMETERS:
;   A structure with the spectra aligned to the reference spectra in wavelength.
;
; OPTIONAL OUTPUT PARAMETERS:
;   coeffs -
;      The coefficients of the best fit.
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
;
;
; REVISION HISTORY:
;   Revision: $Id: align2kurucz.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;-
;*******************************************************************

function my_amoebaK, coeffs
    common myamoebaK_common, sx, sy, rx, ry, prismtemp, no_deriv, wave_range, use_rms, imode

    ;newx = ccd2lambda(iMode, sx, prismtemp, gamz=coeffs[0], pixsize=coeffs[1], wedge=coeffs[2]) 
    ;newx = ccd2lambda(iMode, sx, prismtemp, gamz=coeffs[0], pixsize=coeffs[1], focal=coeffs[2])
    ;newx = ccd2lambda(iMode, sx, prismtemp, gamz=coeffs[0], pixsize=coeffs[1])
    newx = ccd2lambda(iMode, sx, prismtemp, gamz=coeffs[0])
    newy=interpol(sy,newx,rx, /lsq)

    if no_deriv eq 0 then newy = deriv(rx,newy)
    pos = where(rx ge wave_range[0] and rx le wave_range[1], count)
    if count eq 0 then return,0d

    if use_rms eq 1 then begin
        corr_value = sqrt(total((newy[pos] - ry[pos])^2d)/count)
    endif else begin
        corr_value = abs(1d - abs(correlate(newy[pos], ry[pos], /double))) * 1d10
    endelse

    ;print,corr_value,format='(A, $)'
    return, corr_value

end

;*******************************************************************

function align2kurucz, startTime, stopTime, instrumentModeId, $
         gps=gps, missionDays=missionDays, julianDays=julianDays, $
         version=version, wrange=wrange, reffile=reffile, rms=rms, $
         spect=spect, noderiv=noderiv, coeffs=coeffs, goodness=goodness

    common myamoebaK_common
    common SIM_DEF, sim_model

    if n_elements(version) eq 0 then version=19
    if n_elements(wrange) eq 0 then wave_range=[200d,2400d] else wave_range=wrange
    if n_elements(rms) eq 0 then use_rms=0 else use_rms=1
    if n_elements(noderiv) eq 0 then no_deriv=0 else no_deriv=1
    goodness=0d
    coeffs=0d

    maxiter = 5000d
    ; check if the common block has been defined
    if n_elements(sim_model) eq 0 then get_sim_model

    iMode_pos=where(sim_model.imode eq instrumentModeId, count)
    if count eq 0 then return,-1 else iMode=instrumentModeId

    if n_elements(reffile) eq 0 then begin
        reffile='/Users/sbeland/SORCE/data/convolvechkur.sav'
    endif
    restore, reffile
    
    ; get Kurucz irradiance 
    if instrumentModeId eq 31 or instrumentModeId eq 32 then begin
       waveK = chkur_esr.wavelength
       irradK = chkur_esr.irradiance
    endif else if instrumentModeId eq 41 or instrumentModeId eq 45 then begin
       waveK = chkur_vis.wavelength
       irradK = chkur_vis.irradiance
    endif else if instrumentModeId eq 43 or instrumentModeId eq 47 then begin
       waveK = chkur_uv.wavelength
       irradK = chkur_uv.irradiance
    endif else if instrumentModeId eq 44 or instrumentModeId eq 45 then begin
       waveK = chkur_ir.wavelength
       irradK = chkur_ir.irradiance
    endif

    ; check if the common block has been defined
    if n_elements(sim_model) eq 0 then get_sim_model

    iMode_pos=where(sim_model.imode eq instrumentModeId, count)
    if count eq 0 then return,-1 else iMode=instrumentModeId

    ; get the spectra to fit to the reference
    if size(spect,/tname) ne 'STRUCT' then begin
        spect = get_sim_spectra(startTime, stopTime, instrumentModeId, gps=gps, missionDays=missionDays, $
            julianDays=julianDays, version=version, /noCcdCorr)
        if size(spect,/tname) ne 'STRUCT' then begin
            if keyword_set(verbose) then print,'No valid spectra found'
            return,-1
        endif
    endif

    ; to use align_spectra2, we'll replace the ccdpos by wavelengths and dn by the irradiance
    outspect=spect
    outspect.wavelength = ccd2lambda(iMode, outspect.ccdpos_uncorr, outspect.prismtemp) 

    new_irradK = interpol(irradK, waveK, spect.wavelength, /spline)

    ;sx=outspect.wavelength
    sx=outspect.ccdpos_uncorr
    sy=outspect.irradiance
    prismtemp=outspect.prismtemp
    rx = outspect.wavelength
    if no_deriv eq 0 then ry=deriv(rx, new_irradK) else ry=new_irradK

    ;coeffs=[sim_model.gamz[iMode_pos], sim_model.pix[iMode_pos], sim_model.thp[iMode_pos]]
    ;coeffs=[sim_model.gamz[iMode_pos], sim_model.pix[iMode_pos], sim_model.foc[iMode_pos]]
    ;coeffs=[sim_model.gamz[iMode_pos], sim_model.thp[iMode_pos]]
    ;coeffs=[sim_model.gamz[iMode_pos], sim_model.pix[iMode_pos]]
    coeffs=[sim_model.gamz[iMode_pos]]

    coeffs = AMOEBA(1.0d-8, P0=coeffs, scale=[1d-3,1d-3], FUNCTION_VALUE=fval, FUNCTION_NAME='my_amoebaK', nmax=maxiter)

    goodness= 1d - fval;/1d10
    if keyword_set(verbose) then print,'Final coefficients = ',coeffs
    if keyword_set(verbose) then print,'Fit goodness=',goodness

    if n_elements(coeffs) eq 1 and coeffs[0] eq -1 then begin
        print,'AMOEBA failed to converge'
        return,-1
    endif

    outspect.wavelength = ccd2lambda(iMode, outspect.ccdpos_uncorr, outspect.prismtemp, $
        ;gamz=coeffs[0], pixsize=coeffs[1], wedge=coeffs[2] $
        ;gamz=coeffs[0], pixsize=coeffs[1] $
        gamz=coeffs[0]$
        ) 


	; define the goodness of the fit by the correlate
    new_irradK = interpol(irradK, waveK, outspect.wavelength, /spline)
    goodness=correlate( deriv(outspect.wavelength,outspect.irradiance), deriv(outspect.wavelength,new_irradK) )
 
    return, outspect

end



