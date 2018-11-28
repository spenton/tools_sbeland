;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Converts the SimCorrectedIrradiance save files to a format similar to
;   Level3 data structure (daily average and fixed wavelength grid) to
;   be able to sue "LISIRD" to look at the data.
;
; CALLING SEQUENCE:
;
; RETURNED VALUES:
;   Returns an array of structures similar to what get_level3_spectrum returns
;   with the daily averages and re-sampled at a fixed wavelength grid.
;
; INPUT PARAMETERS:
;   startTime -
;      The lower time range for which data will be returned.
;      Expects times to be in mission days.
;   stopTime -
;      The upper time range for which data will be returned.
;      Expects times to be in mission days.
;   instrumentModeId -
;      Instrument mode to process.
;
; OPTIONAL INPUT PARAMETERS:
;   version -
;      Requests a specific version of the SimProfileIntegralCal, otherwise
;      pick the latest version valid for the specified time range.
;
; RETURNED PARAMETERS:
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

;
; REVISION HISTORY:
;   Revision: $Id: process_corr2level3.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;------------------------------------------------------------------
function process_corr2level3, instrumentModeId, starttime=starttime, endtime=endtime, inspectrum=inspectrum, $
    outfile=outfile

    if size(inspectrum,/tname) ne 'STRUCT' then begin
        if instrumentModeId eq 41 then  begin
            ;restore,'~/SORCE/data/sima_vis_corr_2011.sav'
            ;inspectrum = temporary(visa_corr_2011_CCDSHIFT)
            ;restore,'~/SORCE/data/sima_vis_corr_2011_mod_407_ccdshift.sav'
            restore,'~/SORCE/data/sima_vis_calib_2204.sav'
            inspectrum = temporary(visa_corr_2011)
        endif else if instrumentModeId eq 43 then begin
            ;restore,'~/SORCE/data/sima_uv_corr_2011_s5455_160_kappaALL.sav'
            restore,'~/SORCE/data/sima_uv_calib_2204.sav'
            inspectrum = temporary(uva_corr_2011_160_all)
        endif else if instrumentModeId eq 44 then begin
            ;restore,'~/SORCE/data/sima_ir_corr_2011.sav'
            restore,'~/SORCE/data/sima_ir_calib_2204.sav'
            inspectrum = temporary(IRA_CORR_2011_CCDSHIFT)
        endif else if instrumentModeId eq 31 then begin
            ;restore,'~/SORCE/data/sima_esr_corr_2011.sav'
            restore,'~/SORCE/data/sima_esr_calib_2204.sav'
            inspectrum = temporary(ESRA_CORR_2011_CCDSHIFT)
        endif
    endif

   
    if n_elements(starttime) eq 0 then starttime=0d
    if n_elements(endtime) eq 0 then endtime=5400d
    starttime=min([starttime,endtime],max=endtime)
    t0=sd2gps(starttime)*1d6
    t1=sd2gps(endtime)*1d6

    ; get the standard level3 wavelength from day 500
    res=get_level3_spectrum(500d,501d,instrumentModeId,24d,version=19)
    waves=(res.minwavelength + res.maxwavelength) /2d
    nwaves=n_elements(waves)

    sppos=[]
    mintime=1d99
    for i=0L,n_elements(inspectrum.spect20)-1 do begin 
        pos=where(inspectrum.spect20[i].timestamp gt 0d,count) 
        if count gt 0 then begin 
            mintime=min([mintime,inspectrum.spect20[i].timestamp[pos]])
            p=where(inspectrum.spect20[i].timestamp[pos] ge t0 and inspectrum.spect20[i].timestamp[pos] lt t1,c) 
            if c gt 0 then sppos=[sppos,i] 
        endif 
    endfor
    nspect=n_elements(sppos)
    ; first day to start counting
    mintime=max([mintime,t0])
    mintime=FLOOR(gps2sd(mintime/1d6))

    ; skip the first couple which are often at 0.0
    mysd=gps2sd(inspectrum.spect20[sppos].timestamp[2]/1d6)
    ; get the list of spectra for each day
    hist=histogram(mysd,binsize=1.0,min=mintime,location=xhist,reverse_ind=rind)

    tmp_str={NOMINALGPSTIMETAG:0d, TIMESPANINHOURS:24, INSTRUMENTMODEID:INSTRUMENTMODEID, $
        VERSION:2204, MINWAVELENGTH:0d, MAXWAVELENGTH:0d, IRRADIANCE:0d, IRRADIANCEUNCERTAINTY:0d, QUALITY:0d}

    p=where(hist gt 0d,ndays)
    outdata = replicate(tmp_str, ndays * nwaves)

    ; now process each day using the reverse indices
    daycount=0L
    for day=0L,n_elements(hist)-1L do begin
       if rind[day] eq rind[day+1] then continue
       ; get the indices of the spectrum to process for that day
       sind = rind[rind[day]:rind[day+1]-1]
       irrad=replicate(0d,nwaves)
       navg=replicate(0d,nwaves)
       for i=0,n_elements(sind)-1 do begin
           p=where(inspectrum.spect20[sind[i]].wavelength gt 0d,count)
           if count gt 0 then begin
               s=sort(inspectrum.spect20[sind[i]].wavelength[p])
               y2=spl_init(inspectrum.spect20[sind[i]].wavelength[p[s]], inspectrum.spect20[sind[i]].irradiance[p[s]])
               tmp = spl_interp(inspectrum.spect20[sind[i]].wavelength[p[s]], inspectrum.spect20[sind[i]].irradiance[p[s]], y2, waves)
               ;tmp=interpol(inspectrum.spect20[sind[i]].irradiance[p[s]], inspectrum.spect20[sind[i]].wavelength[p[s]], waves,/spl)
               k0=where(tmp gt 0d and tmp lt 5.0,count)
               if count eq 0 then continue
               irrad[k0] += tmp[k0]
               navg[k0] += 1d
           endif
       endfor

       ; create the daily average
       k=where(navg gt 0,count)
       if count eq 0 then continue
       irrad[k] /= navg[k]
       outdata[daycount*nwaves:(daycount+1)*nwaves-1L].NOMINALGPSTIMETAG = replicate(sd2gps(xhist[day]+0.5d)*1d6, nwaves)
       outdata[daycount*nwaves:(daycount+1)*nwaves-1L].MINWAVELENGTH = waves
       outdata[daycount*nwaves:(daycount+1)*nwaves-1L].MAXWAVELENGTH = waves
       outdata[daycount*nwaves:(daycount+1)*nwaves-1L].IRRADIANCE = irrad
       daycount+=1L
    endfor

    return,outdata

end
