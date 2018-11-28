;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Within the known time range, track the diode responsivity with its
;   temperature.  The time range corresponds to the recovery time after
;   the OBC anomaly in 2007 when the instrument was slowly warming up.
;
; CALLING SEQUENCE:
;   result = GET_DIODE_tempcorr(instrumentModeId, starttime=starttime, stoptime=stoptime)
;
; INPUT PARAMETERS:
;   instrumentModeId -
;      41 for VIS1, 43 for UV or 44 for IR (SimA)
;
; OPTIONAL INPUT PARAMETERS:
;   startTime -
;      The lower time range for which data will be returned.
;      May be specified as julian days, gps microseconds, or
;      SORCE mission day numbers.
;
;   stopTime -
;      The upper time range for which data will be returned
;      May be specified as julian days, gps microseconds, or
;      SORCE mission day numbers.
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
;   
;
; OPTIONAL OUTPUT PARAMETERS:
;      NONE
;
; EXAMPLE:  
;
; COMMON BLOCKS:
;      NONE
;
; REVISION HISTORY:
;   Revision: $Id: get_diode_tempcorr.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;-
;
function get_diode_tempcorr, instrumentModeId, starttime=startTime, stoptime=stopTime, $
         gps=gps, missionDays=missionDays, julianDays=julianDays, version=version, $
         alldata=alldata, noplot=noplot, outfile=outfile, gammaz=gammaz, pixsize=pixsize, $
         nodarks=nodarks, plans=plans

    ; convert the time range in mission days
    if n_elements(starttime) eq 0 or n_elements(stoptime) eq 0 then begin
        ;starttime = 1577.7970d
        ;starttime = 1577.9320d
        ;stoptime =  1578.4390d
        ;stoptime =  1578.0842d
        starttime=3620.3136d
        stoptime=3683.7917d
        missionDays=1
    endif
    ;stablestart=1578.3540d
    ;stablestop =1578.3709d
    ;stablestart=3656.2833d
    ;stablestop= 3656.3008d

    ;starttime=2470.18d
    ;stoptime=2569.65d
    ;stablestart=2477.2010d
    ;stablestop=2477.2178d

    ;starttime = 1577.950d
    ;stoptime =  1578.14d
    stablestart=1578.9444d
    stablestop= 1578.9613d

    starttime = 3591.0d
    ;stoptime =  3629.0d
    stoptime =  3680.0d
    ;stablestart=3602.7875d
    ;stablestop= 3602.8049d

    ; try in the DOOP time range
    ;stablestart=4179.85d
    ;stablestop= 4179.89d

    ; DOOP plans covering the Ir SOlarQuickScan24 (for IR diode)
    ;plans=[ {starttime:4165.8604d, stoptime:4165.9032d}, $
    ;     {starttime:4165.9277d, stoptime:4165.9705d}, $
    ;     {starttime:4167.8782d, stoptime:4167.9214d}, $
    ;     {starttime:4167.9454d, stoptime:4167.9887d}, $
    ;     {starttime:4168.8197d, stoptime:4168.8632d}, $
    ;     {starttime:4168.8870d, stoptime:4168.9305d}, $
    ;     ;{starttime:4169.8285d, stoptime:4169.8724d}, $
    ;     {starttime:4169.8958d, stoptime:4169.9397d}, $
    ;     ;{starttime:4170.8373d, stoptime:4170.8817d}, $
    ;     {starttime:4170.9046d, stoptime:4170.9490d}, $
    ;     {starttime:4171.8461d, stoptime:4171.8911d}, $
    ;     {starttime:4171.9133d, stoptime:4171.9584d}, $
    ;     {starttime:4172.8548d, stoptime:4172.9006d}, $
    ;     ;{starttime:4172.9221d, stoptime:4172.9679d}, $
    ;     ;{starttime:4173.8636d, stoptime:4173.9102d}, $
    ;     {starttime:4175.8811d, stoptime:4175.9100d}, $
    ;     {starttime:4175.9484d, stoptime:4175.9770d}, $
    ;     {starttime:4176.8227d, stoptime:4176.8520d}, $
    ;     {starttime:4176.8900d, stoptime:4176.9180d}, $
    ;     {starttime:4177.7644d, stoptime:4177.7930d}, $
    ;     {starttime:4177.8316d, stoptime:4177.8600d}, $
    ;     {starttime:4178.6389d, stoptime:4178.6670d}]
    ;starttime = min(plans.starttime, max=stoptime)

    if keyword_set(missionDays) then begin
       ; do nothing here since already in mission days
    endif else if keyword_set(julianDays) then begin
       ; user specified time in julian days
       starttime = jd2sd(starttime)
       stoptime = jd2sd(stoptime)
    endif else begin
       ; user specified time in gps micro seconds
       starttime = gps2sd(starttime/1.d6)
       stoptime = gps2sd(starttime/1.d6)
    endelse

    if instrumentModeId lt 41 or instrumentModeId gt 45 then begin
        print,'Error: invalid instrumentModeId (should be between 41 and 48)'
        return,-1
    endif

    if n_elements(version) eq 0 then version=23

    case instrumentModeId of
        41: begin
            ; SIMA VIS1
            sima=1
            simb=0
            instrument='SIMA'
            detector='VIS1'
            prismId = 101862L
            prismTempId = 101909L
            detectorId = 100043L
            shutterId = 102121L
            detectorTemp=100778L
            adcgain = -3276.8d
            feedbackR = 563500.0d
            smooth_factor=9
            ;wrange=[330d, 860d]
            ;wrange=[350d, 700d]
            wrange=[320d, 960d]
            ;gammaZ = 0.99546587765d
            ;pixSize=-0.00128866120453d
        end
        42: begin
            ; SIMA VIS2
            sima=1
            simb=0
            instrument='SIMA'
            detector='VIS2'
            prismId = 101862L
            prismTempId = 101909L
            detectorId = 100055L
            shutterId = 102121L
            detectorTemp=100963L
            adcgain = -3276.8d
            feedbackR = 563500.0d
            smooth_factor=9
        end
        43: begin
            ; SIMA UV
            sima=1
            simb=0
            instrument='SIMA'
            detector='UV'
            prismId = 101862L
            prismTempId = 101909L
            detectorId = 100046L
            shutterId = 102121L
            detectorTemp=101069L
            adcgain = -3276.8d
            feedbackR = 50090000.0d
            smooth_factor=1
            wrange=[212d, 307d]
        end
        44: begin
            ; SIMA IR
            sima=1
            simb=0
            instrument='SIMA'
            detector='IR'
            prismId = 101862L
            prismTempId = 101909L
            detectorId = 100052L
            shutterId = 102121L
            detectorTemp=102112L
            adcgain = 3276.8d
            feedbackR = 249500.0d
            smooth_factor=1
            ;wrange=[880d, 1570d]
            ;wrange=[1050d, 1500d]
            ;wrange=[920d,1650d]
            wrange=[850d,1670d]
            ;gammaZ = 0.99546587765d
            ;pixSize=-0.00128866120453d
        end
    endcase

    ; check if we already have the data on hand
    if size(alldata,/tname) ne 'STRUCT' then begin
        ; get our golden standard spectra from day 453
        sd453=get_sim_spectra(453.67840d, 453.69524d, instrumentModeId, version=version, /mission,/noCCDCorr,/dn_only, profile_data=profile_data)

        darks=[{timetag:0.0d, value:0.0d}, {timetag:sd2gps(8000d)*1d6, value:0.0d}]
        if ~ keyword_set(nodarks) then begin
            ; get the plans for the specified timerange for Dark
            dark_plans=get_sorce_plan(starttime, stoptime,/mission,sima=sima,simb=simb,activity='Dark')
            ndarks=n_elements(dark_plans)
            darks = replicate({timetag:0.0d, value:0.0d}, ndarks)

            ; get the average dark value for each Dark activity
            print,'Extracting '+strtrim(string(ndarks),2)+' darks ...'
            for i=0,ndarks-1 do begin
                get_sorce_telemetry, dark_data, info, dark_plans[i].starttime, dark_plans[i].stoptime, /mission, $
                    tlmId=[detectorId,shutterId,detectorTemp]
                if size(dark_data,/tname) ne 'STRUCT' then begin
                    print,'No dark data found for the specified instrumentId and time range'
                    return,-1
                endif else if (*dark_data.(0)).n_science eq 0 then begin
                    print,'No dark data found for the specified instrumentId and time range'
                    return,-1
                endif
                ; only keep the data when the shutter is closed
                darkpos = where((*dark_data.(1)).science.dn eq 1,count)
                ; and trim 5 seconds at both ends (Dark activity lasts about 30 seconds)
                t0 = (*dark_data.(1)).science.timetag[darkpos[0]] + 5.0d6
                t1 = (*dark_data.(1)).science.timetag[darkpos[-1]] - 5.0d6
                pos = where((*dark_data.(0)).science.timetag ge t0 and $
                            (*dark_data.(0)).science.timetag le t1, dcount0)
                if dcount0 eq 0 then pos=darkpos
                darks[i].timetag = median((*dark_data.(0)).science[pos].timetag)
                darks[i].value = median((*dark_data.(0)).science[pos].dn)
            endfor
        endif
        print,'Darks (dn) = ',darks.value

        ; get our stable spectrum to use as our baseline
        stable_spect = get_sim_spectra(stablestart, stablestop, instrumentModeId, /mission,/noCCDCorr,/dn_only, /nodark, profile_data=profile_data)
        ; get the dark value at each timestamp by linear interpolation with closest dark values
        dark_values = interpol(darks.value, darks.timetag, stable_spect.timetag)
        stable_spect.dn -= dark_values

        ; align the spectrum with SD=453 for wavelength reference
        stable_align=align_spectra2(0d,1d, instrumentModeId, /mission,version=version,coeffs=coeff,spect=stable_spect,refspect=sd453,$
            profile_data=profile_data, wrange=wrange, /no_dettemp, /no_plot) ;no_plot=noplot)
        ;stable_align=align_spectra(0d,1d, instrumentModeId, /mission,version=version,coeffs=coeff,spect=stable_spect,refspect=sd453,$
        ;    profile_data=profile_data, order=1); , wrange=wrange)
        ;stable_align=align_test(0d,1d, instrumentModeId, /mission,version=version,coeffs=coeff,spect=stable_spect,refspect=sd453,$
        ;    profile_data=profile_data, order=1);, wrange=wrange)
        if n_elements(gammaz) gt 0 and n_elements(pixsize) gt 0 then begin
            stable_align.wavelength = ccd2lambda(instrumentModeId, stable_align.ccdpos, stable_align.prismtemp, gamz=gammaz, pixsize=pixsize) 
        endif else if n_elements(gammaz) then begin
            stable_align.wavelength = ccd2lambda(instrumentModeId, stable_align.ccdpos, stable_align.prismtemp, gamz=gammaz)
        endif else if n_elements(pixsize) then begin
            stable_align.wavelength = ccd2lambda(instrumentModeId, stable_align.ccdpos, stable_align.prismtemp, pixsize=pixsize)
        endif
        stable_spect_wave = stable_align.wavelength
        stable_spect_dn = stable_align.dn
        nwave=n_elements(stable_spect_wave)

        ; get the plans for the specified timerange for SolarquickScan24
        if n_elements(plans) eq 0 then begin
            plans=get_sorce_plan(starttime, stoptime,/mission,sima=sima,simb=simb,activity='SolarQuickScan24')
            k=where(abs(plans.starttime - 1578.4890) gt 0.0001)
            plans=plans[k]
        endif
        nplans=n_elements(plans)

        alldata=replicate({deltadn_dt:dblarr(nwave), refTemp:dblarr(nwave), diodeTemp:dblarr(nwave), $
                           refDn:dblarr(nwave), diodeDn:dblarr(nwave), wavelength:dblarr(nwave), timetag:dblarr(nwave)},nplans)

        for i=0,nplans-1 do begin
            print, '   plan '+strtrim(string(i+1),2)+'/'+strtrim(string(nplans),2)
            ; for each planned SolarQuickScan24 extract the spectra and apply the dark correction
            spect = get_sim_spectra(plans[i].starttime, plans[i].stoptime, instrumentModeId, /mission,/noCCDCorr,/dn_only, /nodark, profile_data=profile_data)
            if size(spect,/tname) ne 'STRUCT' then continue

            ; get the dark value at each timestamp by linear interpolation with closest dark values
            dark_values = interpol(darks.value, darks.timetag, spect.timetag)
            spect.dn -= dark_values

            ; align the spectrum with SD=453 for wavelength reference
            spect_align=align_spectra2(0d,1d, instrumentModeId, /mission,version=version,coeffs=coeff,spect=spect,refspect=sd453,$
                profile_data=profile_data, wrange=wrange, /no_dettemp, /no_plot) ;no_plot=noplot)
            ;spect_align=align_spectra(0d,1d, instrumentModeId, /mission,version=version,coeffs=coeff,spect=spect,refspect=sd453,$
            ;    profile_data=profile_data, order=1);, wrange=wrange)
            ;spect_align=align_test(0d,1d, instrumentModeId, /mission,version=version,coeffs=coeff,spect=spect,refspect=sd453,$
            ;    profile_data=profile_data, order=1);, wrange=wrange)
            if n_elements(gammaz) gt 0 and n_elements(pixsize) gt 0 then begin
                spect_align.wavelength = ccd2lambda(instrumentModeId, spect_align.ccdpos, spect_align.prismtemp, gamz=gammaz, pixsize=pixsize)
            endif else if n_elements(gammaz) then begin
                spect_align.wavelength = ccd2lambda(instrumentModeId, spect_align.ccdpos, spect_align.prismtemp, gamz=gammaz)
            endif else if n_elements(pixsize) then begin
                spect_align.wavelength = ccd2lambda(instrumentModeId, spect_align.ccdpos, spect_align.prismtemp, pixsize=pixsize)
            endif
            spect_wave = spect_align.wavelength
            ; interpolate the dns in our stable spectrum wavelength grid
            spect_dn = interpol(spect.dn, spect_wave, stable_spect_wave,/lsq)
            timetag = interpol(spect.timetag, spect_wave, stable_spect_wave, /lsq)

            if NOT keyword_Set(noplot) then begin
                title='Spect '+strtrim(string(i),2)+' SD='+strtrim(string(plans[i].starttime),2)
                title=title+' T='+strtrim(string(median(spect.detectortemp)),2)
                if i eq 0 then begin
                    ; plot the reference spectra
                    lineplot,stable_spect_wave,deriv(stable_spect_wave,spect_dn),xtitle="Wavelength (nm)",$
                        ytitle="Derivative Corrected DN", ptitle=instrument+' '+detector,psym=-3,charsize=1.2, $
                        title='SD453',thick=1.5
                endif
                lineplot,stable_spect_wave,deriv(stable_spect_wave,spect_dn),xtitle="Wavelength (nm)",$
                    ytitle="Derivative Corrected DN", $
                    title=title, ptitle=instrument+' '+detector,psym=-3,charsize=1.2
            endif

            ; compare our dn values with the close by, temperature stable/wavelength aligned spectrum
            ; we define deltadn_dt as a fractional difference (change in responsivity)
            delta_t = stable_spect.detectortemp - spect.detectortemp
            deltadn_dt=delta_t*0d
            p=where(delta_t eq 0d, comp=cp,np)
            if np gt 0 then deltadn_dt[p]=0d
            deltadn_dt[cp] = (stable_spect_dn[cp] - spect_dn[cp]) / spect_dn[cp] / delta_t[cp]
            alldata[i].deltadn_dt = temporary(deltadn_dt)
            alldata[i].refTemp = stable_spect.detectortemp
            alldata[i].diodeTemp = spect.detectortemp
            alldata[i].refDN = stable_spect_dn
            alldata[i].diodeDN = spect_dn
            alldata[i].wavelength = stable_spect_wave
            alldata[i].timetag = timetag
            ; only keep data that have a diode temperature difference greater than 1C from stable_spect
            if abs(median(alldata[i].diodeTemp) - median(stable_spect.detectortemp)) lt 0.05d then alldata[i].refTemp*=0d
        endfor
        p=where(alldata.refTemp[0] ne 0d)
        alldata=alldata[p]
    endif

    ; process each wavelength to get the dndt and reference temperature (zero crossing)
    nwave = n_elements(alldata[0].wavelength)
    return_data = {diodeDetTempCoef:dblarr(nwave,4), diodeAmpDetTempCoef:dblarr(nwave,4), refTemp:dblarr(nwave), $
                   stableTemp:dblarr(nwave), wavelength:dblarr(nwave)}
    if instrumentModeId eq 41 then begin
        ; lets combine all the data together and filter with temperature and wavelength ranges
        ; lets overlay each of these at 1100nm (assuming there's no degradation there)
        ; this should take care of solar changes in irradiance
        newwave=[]
        newdndt=[]
        newx=dindgen((962d -309d)*10d)/10d + 309d
        mn=min(abs(newx - 550.0d),pos)
        for i=0,n_elements(alldata)-1 do begin
            p=where(alldata[i].deltadn_dt gt 0d,count)
            if count eq 0 then continue
            ;temp=interpol(alldata[i].deltadn_dt[p], alldata[i].wavelength[p], newx, /lsq)
            ;temp -= temp[pos]
            ;newwave=[newwave, newx]
            ;newdndt=[newdndt, temp]
            newwave=[newwave, alldata[i].wavelength[p]]
            newdndt=[newdndt, alldata[i].deltadn_dt[p]]
        endfor
        s=sort(newwave)
        newwave=newwave[s]
        newdndt=newdndt[s]
        p=where(newwave ge 435d)
        fit=exponential_fit(newwave[p],newdndt[p],yfit=yfit)
        newz=(fit[0] * exp(fit[1]*newx) + fit[2]) > 0d
        return_data={wavelength:newx, tempcorr:newz, reftemp:replicate(24.12d,n_elements(newx))}
        if NOT keyword_set(noplot) then begin
            plot_multi,newwave,newdndt,newx,newz,xtitle='Wavelength (nm)',ytitle='Fractional dDN/dT',$
            title='Mode='+strtrim(string(instrumentModeId),2)+' Temperature Correction',/xst,/yst,charsize=1.4,psym=[4,-3],thick=[1,3], $
            label=['dDN/dT','BSpline Fits']
        endif
        if n_elements(outfile) gt 0 then $
            forprint,return_data.wavelength, return_data.tempcorr, return_data.refTemp, format='(g,g,d8.2)',text=outfile
        return, return_data

    endif else if instrumentModeId eq 43 then begin
        mult=1d

    endif else if instrumentModeId eq 44 then begin
        ; lets overlay each of these at 1100nm (assuming there's no degradation there)
        ; this should take care of solar changes in irradiance
        newwave=[]
        newdndt=[]
        newx=dindgen(1700-850) + 850d
        mn=min(abs(newx - 1150d),pos)
        for i=0,n_elements(alldata)-1 do begin
            temp=interpol(alldata[i].deltadn_dt, alldata[i].wavelength, newx, /lsq)
            ;temp -= temp[pos]
            newwave=[newwave, newx]
            newdndt=[newdndt, temp]
        endfor
        s=sort(newwave)
        newwave=newwave[s]
        newdndt=newdndt[s]
        sset=bspline_iterfit(newwave,newdndt,maxiter=0,requiren=20,bkspace=20)
        newz=bspline_valu(newx,sset)
        return_data={wavelength:newx, tempcorr:newz, reftemp:replicate(median(stable_spect.detectortemp),n_elements(newx))}
        ;return_data={wavelength:newx, tempcorr:newz, reftemp:replicate(24.12d,n_elements(newx))}
        if NOT keyword_set(noplot) then begin
            plot_multi,newwave,newdndt,newx,newz,xtitle='Wavelength (nm)',ytitle='Fractional dDN/dT',$
            title='Mode='+strtrim(string(instrumentModeId),2)+' Temperature Correction',/xst,/yst,charsize=1.4,psym=[4,-3],thick=[1,3], $
            label=['dDN/dT','BSpline Fits']
        endif
        if n_elements(outfile) gt 0 then $
            forprint,return_data.wavelength, return_data.tempcorr, return_data.refTemp, format='(g,g,d8.2)',text=outfile
        return, return_data

    endif


    order=3
    ;return_data={wavelength:dblarr(nwave), tempcorr:dblarr(nwave,order+1), reftemp:replicate(24.12d,nwave)}
    return_data={wavelength:dblarr(nwave), tempcorr:dblarr(nwave,order+1), reftemp:replicate(median(stable_spect.detectortemp),nwave)}

    for i=0L,nwave-1L do begin
        ; performing a linear fit of dn vs T at each wavelength
        ; the slope will be our DeltaDN/DeltaT and our reference temperature will be 
        ; the temperature at which the DeltaDn is 0.0
        ;
        ; Since the production code corrects the responsivity as follows:
        ; 
        ; correctedProfileIntegral= uncorrectedProfileIntegral + tempcorr * (detTemp - refTemp)
        ;
        ; we will calculate the slope of (refDn - diodeDn) vs diodeTemp to get the sign right
        ;
        ; hopefully the robust_poly_fit will ignore bad data points
        if instrumentModeId eq 41 then begin
            ; for the VIS diode, responsivities behave linearly with all temperatures
            p = where(alldata.diodetemp[i] ge -100.0 and alldata.diodetemp[i] ne 0d and alldata.diodetemp[i] le 22.0)
            if alldata[p[0]].wavelength[i] le 400d then mult=0d else mult=1d
        endif else if instrumentModeId eq 43 then begin
            ; for the UV diode, only keep the data with temperature greater than +9.0 C
            ; anything colder shows a very different slope
            p = where(alldata.diodetemp[i] ge 9.0)
        endif else if instrumentModeId eq 44 then begin
            ; we fit a 5th order polynomial to the ALOG of the data (which provides the best fit)
            ; the correction will now have to be applied as follows:
            ; calulate the deltadn_dt for the wavelengths on either side of the actual wavelength for a specific deltaT
            ; perform a linear interpolation to get the actual deltadn_dt for the desired deltaT and wavelength
            ; p0 = (where(SimDiodeTempCorrectionCal.wavelength le spect.wavelength[i]))[-1]
            ; p1 = (where(SimDiodeTempCorrectionCal.wavelength gt spect.wavelength[i]))[0]
            ; deltadn_dt0 = -1d *  exp(poly((reftemp - spect.diodetemp[i]),cc[p0,*]))
            ; deltadn_dt1 = -1d *  exp(poly((reftemp - spect.diodetemp[i]),cc[p1,*]))
            ; deltadn_dt=interpol([deltadn_dt0, deltadn_dt1], [SimDiodeTempCorrectionCal.wavelength[p0], $
            ;            SimDiodeTempCorrectionCal.wavelength[p1], spect.wavelength[i])
            newx=(alldata[*].reftemp[i] - alldata[*].diodetemp[i])
            ;cc=robust_poly_fit(newx, alog((-1d *alldata.deltadn_dt[i])>0d), order,/double) 
            cc=robust_poly_fit(newx, alldata.deltadn_dt[i], order,/double) 
            return_data.wavelength[i]=alldata[0].wavelength[i]
            return_data.tempcorr[i,*] = cc[0:order]

            if NOT keyword_Set(noplot) then begin
               mx=max((-1d * alldata[*].deltadn_dt[i]))
               ;plot_multi, newx, (-1d * alldata[*].deltadn_dt[i]), newx, exp(poly(newx,cc)),/xst,/yst,$
               ;    yrange=[0d,mx],title=string(alldata[0].wavelength[i], format='(F7.2)')
               plot_multi, newx, alldata[*].deltadn_dt[i], newx, poly(newx,cc),/xst,/yst,$
                   title=string(alldata[0].wavelength[i], format='(F7.2)')
               b=''
               ;READ, B, PROMPT='hit Return to continue ... '
            endif

        endif


        ;; keep all data points for our 3rd order polynomial fit
        ;coeff=robust_poly_fit(alldata[p].diodetemp[i], alldata[p].deltadn_dt[i]*mult,3,sig,/double)
        ;; return_data.diodeDetTempCoef[i] = coeff[1]
        ;return_data.diodeDetTempCoef[i,*] = coeff
        ;coeff=robust_poly_fit(alldata[p].diodetemp[i], alldata[p].deltadn_dt[i]*mult/ADCGAIN/feedbackR,3,sig,/double)
        ;return_data.diodeAmpDetTempCoef[i,*] = coeff
        ;return_data.stableTemp[i] = alldata[0].refTemp[i]
        ;return_data.wavelength[i] = alldata[0].wavelength[i]
        ;if NOT keyword_Set(noplot) and (i mod 10) eq 0 then begin
        ;    plot_multi,alldata[p].diodetemp[i], alldata[p].deltadn_dt[i]*mult, $
        ;        alldata[p].diodetemp[i],poly(alldata[p].diodetemp[i],return_data.diodeDetTempCoef[i,*]), $
        ;        /xst,/yst, yrange=[-0.002,0.004],xtitle='Diode Temperature (C)',ytitle='Delta DN/DT',$
        ;        title='Mode='+strtrim(string(instrumentModeId),2)+'  Wavelength='+strtrim(string(alldata[0].wavelength[i]))
        ;endif
    endfor

    if instrumentModeId eq 44 then return, return_data

    ; The change in responsivity with temperature is well represented by a 3rd order polynomial
    ; for each wavelength (a different fit at each wavelength).
    ; On the other hand, the change in responsivity with wavelength is a smooth function but can not
    ; be well represented by a low order polynomial.

    ; calculate the data as dndt vs wavelength and perform a boxcar smoothing on that
    ; do this for the same set of diode temperatures
    ; then, for each wavelength, obtain a new 3rd order polynomial fit of the 
    ; diode temperature vs dndt (basically what we started with)
    diode_temp = alldata[*].diodetemp[0]
    ndtemp=n_elements(diode_temp)
    dndt = dblarr(nwave,ndtemp)
    dndtAmp = dblarr(nwave,ndtemp)
    for dt=0L,ndtemp-1L do begin
        for wn=0L,nwave-1L do begin
            dndt[wn,dt]=poly(diode_temp[dt],  return_data.diodeDetTempCoef[wn,*])
            dndtAmp[wn,dt]=poly(diode_temp[dt],  return_data.diodeAmpDetTempCoef[wn,*])
        endfor
        dndt[*,dt] = smooth(dndt[*,dt], smooth_factor)
        dndtAmp[*,dt] = smooth(dndtAmp[*,dt], smooth_factor)
    endfor
    
    for wn=0L,nwave-1L do begin
        coeff=robust_poly_fit(diode_temp, dndt[wn,*], 3,sig,/double)
        return_data.diodeDetTempCoef[wn,*] = coeff
        coeff=robust_poly_fit(diode_temp, dndtAmp[wn,*], 3,sig,/double)
        return_data.diodeAmpDetTempCoef[wn,*] = coeff
    endfor

    ; smooth out the dndt data
    return_data.refTemp = (round(mean(return_data.stableTemp)*1000.0d))/1000.0d

    return, return_data
   
end
