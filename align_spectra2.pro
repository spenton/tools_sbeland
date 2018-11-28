;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Fit a spectra to a reference spectra by using the
;   AMOEBA routine to minimize the difference with the reference when 
;   adjusting the parameters to the spectrograph model.  The routine
;   finds the coefficients of the best fit and returns the original spectra 
;   with the new wavelength scale.
;
; CALLING SEQUENCE:
;   spectra = MPFIT_SPECTRA(t0, t1, instrumentModeId, /missionDays, refT0=t0, refT1=t1)
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
;   tsis -
;      If set, will use the TSIS Sellmeier coefficients to determine the wavelength
;      change with temperature.
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
;   IDL> res=align_spectra(T0,T1,47,/mission,refspect=res43,coeffs=coeffs,status=status)
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
;   Revision: $Id: align_spectra2.pro,v 1.2 2018/11/26 16:16:02 spenton Exp $
;-
;*******************************************************************
function my_amoeba2, coeffs
    common myamoeba2_common, sp, rsp, no_deriv, wave_range, use_rms, imode, dettempcorr, do_tsis
    common SIM_DEF, sim_model

    pos=where(sim_model.imode eq iMode, count)
    if n_elements(coeffs) eq 1 then begin
        gamz=(1d + coeffs[0])*sim_model.gamz[pos]
        newx = ccd2lambda(iMode, sp.x, sp.prismtemp, gamz=gamz[0], tsis=do_tsis)
    endif else if n_elements(coeffs) eq 2 then begin
        gamz=(1d + coeffs[0])*sim_model.gamz[pos]
        pix =(1d + coeffs[1])*sim_model.pix[pos]
        newx = ccd2lambda(iMode, sp.x, sp.prismtemp, gamz=gamz[0], pixsize=pix[0], tsis=do_tsis)
    endif else if n_elements(coeffs) eq 3 then begin
        gamz=(1d + coeffs[0])*sim_model.gamz[pos]
        yab =(1d + coeffs[1])*sim_model.yab[pos]
        cz =(1d + coeffs[2])*sim_model.cz[pos]
        newx = ccd2lambda(iMode, sp.x, sp.prismtemp, gamz=gamz[0], yab=yab[0], cz=cz[0], tsis=do_tsis)
    endif
    if size(dettempcorr,/tname) eq 'STRUCT' then begin
        ; detCorr is already resampled to refSpect wavelength space
        corr=interpol(dettempcorr.tempcorr, dettempcorr.wavelength, newx)
        ;newy = sp.y / (1d - corr*(dettempcorr.refTemp - sp.detectorTemp))
        ; WE NEED TO MULTIPLY (NOT DIVIDE)
        newy = sp.y * (1d - corr*(dettempcorr.refTemp - sp.detectorTemp))
    endif else newy=sp.y

    s=sort(newx)
    newx=newx[s]
    newy=newy[s]
    if no_deriv eq 0 then newy = deriv(newx,newy)
    y2=spl_init(newx, newy)
    good=where(finite(y2) eq 1,count)
    if count eq n_elements(newy) then newy=spl_interp(newx,newy, y2, rsp.x) else newy=interpol(newy, newx, rsp.x)
    pos = where(rsp.x ge wave_range[0] and rsp.x le wave_range[1] and finite(rsp.y) and finite(newy), count)
    if count eq 0 then return,1d

    if use_rms eq 1 then begin
        corr_value = sqrt(total((newy[pos] - rsp.y[pos])^2d))
    endif else begin
        corr_value = abs(1d - abs(correlate(newy[pos], rsp.y[pos], /double)))
    endelse

    ;print,coeffs, corr_value,n_elements(pos),n_elements(sp.x),n_elements(rsp.x),format='(2(G,","),G,3(",  ",I))'
    return, corr_value

end
;*******************************************************************

function align_spectra2, startTime, stopTime, instrumentModeId, $
         gps=gps, missionDays=missionDays, julianDays=julianDays, $
         refT0=refT0, refT1=refT1, refSpect=refSpect, spect=spect, $
         version=version, quiet=quiet, coeffs=coeffs, $
         noderiv=noderiv, status=status, profile_data=profile_data, $
		 verbose=verbose, goodness=goodness, wrange=wrange, rms=rms, $
         ccdfit=ccdfit, no_plot=no_plot, database=database, no_dettemp=no_dettemp, $
         detCorr=detCorr, order=order, solardist=solardist, no1au=no1au, tsis=tsis

    common myamoeba2_common
    common SIM_DEF, sim_model

    if keyword_set(tsis) then do_tsis=1 else do_tsis=0
    if n_elements(version) eq 0 then version=19
    if n_elements(database) eq 0 then database='SORCE_SIM_V19'
    if keyword_set(noderiv) then no_deriv=1 else no_deriv=0
    if keyword_set(rms) then use_rms=1 else use_rms=0
	if n_elements(order) eq 0 then order=3
    ; adjust the wavelength range for each detector
    if n_elements(wrange) gt 0 then begin
        wave_range=wrange
    endif else if instrumentmodeId eq 41 or instrumentModeId eq 45 then begin
        wave_range=[310d,950d]
    endif else if instrumentmodeId eq 43 or instrumentModeId eq 47 then begin
        wave_range=[220d,305d]
    endif else if instrumentmodeId eq 44 or instrumentModeId eq 48 then begin
        wave_range=[950d,1600d]
    endif else begin
        wave_range=[220d,2400d]
    endelse
    goodness=0d
    coeffs=0d

    maxiter = 5000d

    ; check if the common block has been defined
    if n_elements(sim_model) eq 0 then get_sim_model

    iMode_pos=where(sim_model.imode eq instrumentModeId, count)
    if count eq 0 then return,-1 else iMode=instrumentModeId

    ; if a reference spectra was provided use that instead
    if size(refSpect,/tname) ne 'STRUCT' then begin
        if keyword_set(verbose) then print,'getting the reference spectra to align to ...'
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
            refSpect = get_sim_spectra(refT0, refT1, instrumentModeId, gps=gps, missionDays=missionDays, no1au=no1au, /nodark, $
                julianDays=julianDays, version=version, /noCcdCorr, profile_data=profile_data, solardist=solardist, /uncorr, tsis=do_tsis)
        endif else begin
            ; extract the ESR data from the database for version 19
            refSpect = get_esr_spectra(refT0, refT1, instrumentModeId, gps=gps, missionDays=missionDays, no1au=no1au, $
                julianDays=julianDays, version=version, database=database, /no_duplicate, solardist=solardist, /uncorr)
        endelse

        if size(refSpect,/tname) ne 'STRUCT' then begin
            if keyword_set(verbose) then print,'No valid reference spectra found'
            return,-1
        endif
        ; use the instrument model to define our reference wavelength scale (at the actual prism temp)
        ;refSpect.wavelength = ccd2lambda(iMode, refSpect.ccdpos_uncorr, refSpect.prismtemp)
    endif

    ; get the spectra to fit to the reference
    if size(spect,/tname) ne 'STRUCT' then begin
        if keyword_set(verbose) then print,'getting the spectra to align ...'
        if instrumentModeId gt 32 then begin
            spect = get_sim_spectra(startTime, stopTime, instrumentModeId, gps=gps, missionDays=missionDays, no1au=no1au, /nodark, $
                julianDays=julianDays, version=version, /noCcdCorr, profile_data=profile_data, solarDist=solarDist, /uncorr, $
                /noTempCorr)
        endif else begin
            spect = get_esr_spectra(startTime, stopTime, instrumentModeId, gps=gps, missionDays=missionDays, no1au=no1au, /nodark, $
                julianDays=julianDays, version=version, database=database, solarDist=solarDist,/no_duplicate, /uncorr)
        endelse
        if size(spect,/tname) ne 'STRUCT' then begin
            if keyword_set(verbose) then print,'No valid spectra found'
            return,-1
        endif
    endif

    ;refSpect.dn_tempcorr = refSpect.dn
    ;spect.dn_tempcorr = spect.dn

    ; before missionDay 190, we had glint issues -> remove affected ccdpos
    if gps2sd(spect[0].timetag/1d6) le 190.0 then begin
        if instrumentModeId eq 41 then p=where(spect.ccdpos ge 12260L and spect.ccdpos lt 12500L, comp=cp)
        if instrumentModeId eq 43 then p=where(spect.ccdpos ge 34720L and spect.ccdpos lt 34920L, comp=cp)
        if instrumentModeId eq 45 then p=where(spect.ccdpos ge 50150L and spect.ccdpos lt 50380L, comp=cp)
        if instrumentModeId eq 47 then p=where(spect.ccdpos ge 27465L and spect.ccdpos lt 27700L, comp=cp)
        if n_elements(p) gt 0 then spect=spect[cp]
    endif
    p=where(spect.ccdpos gt 0,count)
    if count lt n_elements(spect) then spect=spect[p]

    if keyword_set(no_dettemp) then begin
        detCorr=-1
    endif else if (instrumentModeId eq 41 or instrumentModeId eq 45) and size(detCorr,/tname) ne 'STRUCT' then begin
        if keyword_set(verbose) then print,'getting the detector temparature correction ...'
        ;infile = file_search("$SORCE_DATA","sima_vis_tempcorr_20.txt",/expand_env,/full)
        infile = file_search("$SORCE_DATA","sima_vis_tempcorr_6.txt",/expand_env,/full)
        ; the following file is the results of lab measurements in Feb/2014
        ;infile = file_search("$SORCE_DATA","SORCE_Si_drdt_fit.txt",/expand_env,/full)
        readcol,infile[0],wv,c0,tref,format='(D,D,D)',/silent
        corr=interpol(c0, wv, refSpect.wavelength, /lsq)
        ;refSpect.dn_tempcorr = refSpect.dn / (1d - corr*(tref[0] - refSpect.detectorTemp))
        refSpect.dn_tempcorr = refSpect.dn * (1d - corr*(tref[0] - refSpect.detectorTemp))
        detCorr={wavelength:wv, tempcorr:c0, refTemp:tref}
    endif else if (instrumentModeId eq 44 or instrumentModeId eq 48) and size(detCorr,/tname) ne 'STRUCT' then begin
        if keyword_set(verbose) then print,'getting the detector temparature correction ...'
        ;infile = file_search("$SORCE_DATA","sima_ir_tempcorr_21.txt",/expand_env,/full)
        infile = file_search("$SORCE_DATA","sima_ir_tempcorr_24.txt",/expand_env,/full)
        ;infile = file_search("$SORCE_DATA","sima_ir_tempcorr_25.txt",/expand_env,/full)
        ; the following file is the results of lab measurements in Feb/2014
        ;infile = file_search("$SORCE_DATA","SORCE_InGaAs_drdt_fit.txt",/expand_env,/full)
        readcol,infile[0],wv,c0,tref,format='(D,D,D)',/silent
        corr=interpol(c0, wv, refSpect.wavelength, /lsq)
        ;refSpect.dn_tempcorr = refSpect.dn / (1d - corr*(tref[0] - refSpect.detectorTemp))
        refSpect.dn_tempcorr = refSpect.dn * (1d - corr*(tref[0] - refSpect.detectorTemp))
        detCorr={wavelength:wv, tempcorr:c0, refTemp:tref}
    endif else if  size(detCorr,/tname) ne 'STRUCT' then detCorr=-1
    dettempcorr=detCorr

    ; match the same ccdpos between the 2 scans
    match, refSpect.ccdpos, spect.ccdpos, suba, subb
    if n_elements(suba) lt 20 then begin
        print,'Not enough matching points -> interpolating'
        ;return,-1
        ; resample refspect to spect
        ; make sure we keep the data from only the latest first ccdpos and the first of the last ccdpos
        p=where(spect.ccdpos eq spect[0].ccdpos,count)
        if count gt 1 then spect[p[0:-2]].ccdpos = 0d
        p=where(spect.ccdpos eq spect[-1].ccdpos,count)
        if count gt 1 then spect[p[1:-1]].ccdpos = 0d
        if spect[1].dn le 0d then p=where(spect.ccdpos gt 0d and spect.dn lt 0d,count) else  p=where(spect.ccdpos gt 0d and spect.dn gt 0d,count)
        s=sort(spect[p].ccdpos)
        if (spect[2].ccdpos - spect[1].ccdpos) lt 0 then s=reverse(s)
        k=uniq(spect[p[s]].ccdpos)
        newRefSpect = Spect[p[s[k]]]
        newRefSpect.DN = interpol(Refspect.dn, Refspect.ccdpos, newRefSpect.ccdpos, /lsq)
        newRefSpect.DN_TEMPCORR = interpol(Refspect.dn_tmepcorr, Refspect.ccdpos, newRefSpect.ccdpos, /lsq)
        newRefSpect.PRISMTEMP = interpol(Refspect.PRISMTEMP, Refspect.ccdpos, newRefSpect.ccdpos, /lsq)
        newRefSpect.DETECTORTEMP = interpol(Refspect.DETECTORTEMP, Refspect.ccdpos, newRefSpect.ccdpos, /lsq)
        newRefSpect.TIMETAG = interpol(Refspect.TIMETAG, Refspect.ccdpos, newRefSpect.ccdpos, /lsq)
        ; get the new wavelengths for the reference spectra
        newRefSpect.wavelength = ccd2lambda(iMode, newRefSpect.ccdpos_uncorr, newRefSpect.prismtemp)
        match, newRefSpect.ccdpos, Spect.ccdpos, suba, subb
    endif else newRefSpect=RefSpect

    ; make sure the spectra is in increasing wavelength
    sa=sort(newrefSpect[suba].wavelength)
    sb=sort(Spect[subb].ccdpos)
    if instrumentModeId eq 32 or (instrumentModeId ge 45 and instrumentModeId le 48) then begin
        ; for SimB, the wavelength increases with decreasing CCDPOS
        sb=reverse(sb)
    endif 

    if size(detCorr,/tname) eq 'STRUCT' then begin
        rsp={x:newRefSpect[suba[sa]].wavelength, y:newRefSpect[suba[sa]].dn_tempcorr}
        sp={x:Spect[subb[sb]].ccdpos, y:Spect[subb[sb]].dn, prismTemp:Spect[subb[sb]].prismTemp, detectorTemp:Spect[subb[sb]].detectorTemp}
    endif else begin
        rsp={x:newRefSpect[suba[sa]].wavelength, y:newRefSpect[suba[sa]].dn}
        sp={x:Spect[subb[sb]].ccdpos, y:Spect[subb[sb]].dn, prismTemp:Spect[subb[sb]].prismTemp}
    endelse

    if no_deriv eq 0 then rsp.y=deriv(rsp.x, rsp.y)

    if keyword_set(verbose) then print,'Starting AMOEBA ...'
    ;init_coef=[sim_model.gamz[iMode_pos],sim_model.pix[iMode_pos]]
    ; changed to a fractional perturbation of gammaZ and pixsize
    ;init_coef=[0d, 0d, 0d]
    init_coef=[0d, 0d]
    coeffs = AMOEBA(1d-8, P0=init_coef, scale=0.02d, FUNCTION_VALUE=fval, FUNCTION_NAME='my_amoeba2', ncalls=ncalls)

    if use_rms eq 0 then fval /= 1d10
    goodness= 1d - fval
    if keyword_set(verbose) then print,'Final coefficients = ',coeffs
    if keyword_set(verbose) then print,'Fit goodness=',goodness

    if n_elements(coeffs) eq 1 and coeffs[0] eq -1 then begin
        print,'AMOEBA failed to converge'
        return,-1
    endif

    coeffs[0]=(1d + coeffs[0])*sim_model.gamz[iMode_pos]
    coeffs[1]=(1d + coeffs[1])*sim_model.pix[iMode_pos]
    ;coeffs[1] =(1d + coeffs[1])*sim_model.yab[iMode_pos]
    ;coeffs[2] =(1d + coeffs[2])*sim_model.cz[iMode_pos]

    outspect=spect
    ;outspect.wavelength = ccd2lambda(iMode, outspect.ccdpos_uncorr, outspect.prismtemp, $
    ;    ;gamz=coeffs[0], yab=coeffs[1], cz=coeffs[2])
    ;    gamz=coeffs[0], pixsize=coeffs[1])
    ;    ;gamz=coeffs[0])

    ;if size(detCorr,/tname) eq 'STRUCT' then begin
    ;    corr=interpol(dettempcorr.tempcorr, dettempcorr.wavelength, outspect.wavelength)
    ;    outspect.dn_tempcorr = outspect.dn / (1d - corr*(dettempcorr.refTemp[0] - outspect.detectorTemp))
    ;endif

    ; get the delta_ccdpos vs ccdpos (using a 3rd order polynomial fit)
    ref_ccdpos = lambda2ccd(iMode, refspect.wavelength, spect.prismTemp, tsis=do_tsis)
    new_ccdpos = lambda2ccd(iMode, refspect.wavelength, spect.prismTemp, gamz=coeffs[0], pixsize=coeffs[1], tsis=do_tsis)
    ;new_ccdpos = lambda2ccd(iMode, refspect.wavelength, spect.prismTemp, gamz=coeffs[0], pixsize=coeffs[1], cz=coeffs[2])
    ;new_ccdpos = lambda2ccd(iMode, refspect.wavelength, spect.prismTemp, gamz=coeffs[0], yab=coeffs[1], cz=coeffs[2])
    delta_ccdpos = new_ccdpos - ref_ccdpos
    if order lt 2 then begin
        ccdfit = poly_fit(ref_ccdpos, delta_ccdpos, order, /double)
    endif else begin
        ccdfit = robust_poly_fit(ref_ccdpos, delta_ccdpos, order, /double)
    endelse

    ; get the new ccdpos from the fit
    outspect.ccdpos = outspect.ccdpos_uncorr - poly(outspect.ccdpos_uncorr, ccdfit)
    outspect.wavelength = ccd2lambda(iMode, outspect.ccdpos, outspect.prismtemp, tsis=do_tsis)

    if NOT keyword_set(no_plot) then begin
        plot_multi,ref_ccdpos, delta_ccdpos, outspect.ccdpos_uncorr, poly(outspect.ccdpos_uncorr, ccdfit), $
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
    if refCount gt 2 then begin
        ss=sort(outspect.wavelength)
        rs=sort(refspect[refp].wavelength)
        if size(detCorr,/tname) eq 'STRUCT' then begin
            newy = interpol(outspect[ss].dn_tempcorr, outspect[ss].wavelength, refspect[refp[rs]].wavelength, /quad)
            goodness=correlate(deriv(refspect[refp[rs]].wavelength,newy), deriv(refspect[refp[rs]].wavelength,refspect[refp[rs]].dn_tempcorr))
        endif else begin
            newy = interpol(outspect[ss].dn, outspect[ss].wavelength, refspect[refp[rs]].wavelength, /quad)
            goodness=correlate(deriv(refspect[refp[rs]].wavelength,newy), deriv(refspect[refp[rs]].wavelength,refspect[refp[rs]].dn))
        endelse
	endif
 
    return, outspect

end



