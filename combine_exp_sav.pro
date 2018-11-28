pro combine_exp_sav, data=data, outfile=outfile, noprint=noprint, daily=daily
    ; this routine simply reads the 8 IDL save files
    ; from calc_sim_exp.pro and totals the results
    ; in a meaningful way

    ; the first day with solar exposure is 17.5 - ignore all safehold etc prior
    sdays=[17.0d, 453.00d, 491.02d, 673.02d, 852.02d, 1034.05d, 1223.00d, 1419.03d]
    gpsdays=sd2gps(sdays) * 1d6
    files = ['~/SORCE/data/test9/sima_hrtout_10arc.sav', $
             '~/SORCE/data/test9/sima_hrtin_10arc.sav', $
             '~/SORCE/data/test9/sima_hrtout_FOV.sav', $
             '~/SORCE/data/test9/sima_hrtin_FOV.sav', $
             '~/SORCE/data/test9/simb_hrtout_10arc.sav', $
             '~/SORCE/data/test9/simb_hrtin_10arc.sav', $
             '~/SORCE/data/test9/simb_hrtout_FOV.sav', $
             '~/SORCE/data/test9/simb_hrtin_FOV.sav']

    ;files =  '~/SORCE/data/test3/sima_hrtout_FOV_test3.sav'
    ndays=n_elements(sdays)-1L
    nfiles=n_elements(files)
    indata=PTRARR(nfiles)
    data = replicate({file:'', day:dblarr(ndays), TotalExposure:dblarr(ndays), $
                          totalgap_sci:dblarr(ndays), totalgap_hk:dblarr(ndays),$
                          totalgap_sc_open:dblarr(ndays), totalgap_hrt:dblarr(ndays), $
                          total_safehold:dblarr(ndays), totalgap_abs:dblarr(ndays), $
                          totalfrom_plan:dblarr(ndays)}, nfiles)
    
    sec_per_day = 24.0d * 3600.0d
    for fn=0L,nfiles-1L do begin
        data[fn].file = files[fn]
        restore,files[fn]
        for dn=0L,ndays-1L do begin
            p=where(sim_exp.endtime ge gpsdays[0] and sim_exp.endtime lt gpsdays[dn+1])
            data[fn].day[dn] = sdays[dn+1]
            data[fn].totalexposure[dn] = TOTAL(sim_exp[p].solar_exp) / sec_per_day
            data[fn].totalgap_sci[dn] = TOTAL(sim_exp[p].gap_sci) / sec_per_day
            data[fn].totalgap_sc_open[dn] = TOTAL(sim_exp[p].gap_sc_open) / sec_per_day
            data[fn].totalgap_hrt[dn] = TOTAL(sim_exp[p].gap_hrt) / sec_per_day
            data[fn].totalgap_abs[dn] = TOTAL(sim_exp[p].abs_gap) / sec_per_day
            data[fn].totalfrom_plan[dn] = TOTAL(sim_exp[p].hrt_from_plan) / sec_per_day
            data[fn].total_safehold[dn] = TOTAL(sim_exp[p].safehold) / sec_per_day
        endfor
        indata[fn] = PTR_NEW(sim_exp,/no_copy)
    endfor

    if n_elements(outfile) eq 0 and keyword_set(noprint) and not keyword_set(daily) then return

    ; print the data to file and/or screen
    header=$
'     DAY   TotalExposure   TotalGap_ABS  TotalGap_Sci  TotalGap_SC_Open  TotalGap_HRT Total_Safehold TotalFrom_Plan'

    if n_elements(outfile) gt 0 then begin
        openw,unit,outfile,/get_lun
        printf,unit,outfile
        printf,unit,systime()
        printf,unit,''
        for fn=0L,nfiles-1L do begin
            printf,unit,data[fn].file
            printf,unit, header
            for dn=0L,ndays-1L do begin
            printf,unit,format='(F10.2,7F13.4)',data[fn].day[dn], data[fn].totalexposure[dn], $
                data[fn].totalgap_abs[dn], data[fn].totalgap_sci[dn], $
                data[fn].totalgap_sc_open[dn], data[fn].totalgap_hrt[dn], $
                data[fn].total_safehold[dn], data[fn].totalfrom_plan[dn]
            endfor
            printf,unit,''
            printf,unit,''
        endfor
        free_lun, unit
    endif

    if keyword_set(daily) then begin
        ; calculate the daily cumulative solar_exposure and save it to a file
        ; for example:  for day sd=200.0 we accumulate all solar_exp that ended
        ; before 200.0.
        ; This means that at the start of day 200 (at 200.0) we have the correct
        ; solar_exp up to the start of that day
        print,'calculating cumulative daily exposure ...'
        sdhdr = '#     Day               JD         GPS      Cumulative Solar_exp  Abs_Gap'
        for fn=0,nfiles-1 do begin
            sdfile = file_dirname(files[fn])+PATH_SEP()+file_basename(files[fn],'.sav')+'_sd.txt'
            print, '   '+sdfile
            openw,unit,sdfile,/get_lun
            printf,unit,'# '+sdfile
            printf,unit,'# '+systime()
            printf,unit,'#'
            printf,unit,sdhdr
            sd0=gps2sd(min((*indata[fn]).starttime)/1d6)
            sd1=gps2sd(max((*indata[fn]).endtime)/1d6)
            nsd=ceil(sd1) - floor(sd0)
            sd = dindgen(nsd)+floor(sd0)
            endtime = floor(gps2sd((*indata[fn]).endtime/1d6))
            for i=1L,nsd-1L do begin
                p=max(where(endtime lt sd[i]))
                printf,unit, sd[i], sd2jd(sd[i]), ulong(sd2gps(sd[i])), $
                    total((*indata[fn])[0:p[0]].solar_exp)/86400.0d, $
                    total((*indata[fn])[0:p[0]].abs_gap)/86400.0d
            endfor
            free_lun, unit
        endfor
    endif


    if keyword_set(noprint) then return

    for fn=0L,nfiles-1L do begin
        print,data[fn].file
        print, header
        for dn=0L,ndays-1L do begin
            print,format='(F10.2,7F13.4)',data[fn].day[dn], data[fn].totalexposure[dn], $
                data[fn].totalgap_abs[dn], data[fn].totalgap_sci[dn], $
                data[fn].totalgap_sc_open[dn], data[fn].totalgap_hrt[dn], $
                data[fn].total_safehold[dn], data[fn].totalfrom_plan[dn]
        endfor
        print,''
        print,''
    endfor

    return
end
