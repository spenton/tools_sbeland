pro combine_exp_sav, data=data, outfile=outfile, noprint=noprint
    ; this routine simply reads the 8 IDL save files
    ; from calc_sim_exp.pro and totals the results
    ; in a meaningful way

    ; the first day with solar exposure is 17.5 - ignore all safehold etc prior
    sdays=[17.0d, 453.00d, 491.02d, 673.02d, 852.02d, 1034.05d, 1223.00d, 1419.03d]
    gpsdays=sd2gps(sdays) * 1d6
    files = ['~/SORCE/data/sima_hrtout_10arc.sav', $
             '~/SORCE/data/sima_hrtin_10arc.sav', $
             '~/SORCE/data/sima_hrtout_FOV.sav', $
             '~/SORCE/data/sima_hrtin_FOV.sav', $
             '~/SORCE/data/simb_hrtout_10arc.sav', $
             '~/SORCE/data/simb_hrtin_10arc.sav', $
             '~/SORCE/data/simb_hrtout_FOV.sav', $
             '~/SORCE/data/simb_hrtin_FOV.sav']

    ndays=n_elements(sdays)-1L
    nfiles=n_elements(files)
    data = replicate({file:'', day:dblarr(ndays), TotalExposure:dblarr(ndays), $
                          totalgap_sci:dblarr(ndays), totalgap_hk:dblarr(ndays),$
                          totalgap_sc:dblarr(ndays), totalgap_hrt:dblarr(ndays), $
                          totalgap_abs:dblarr(ndays), totalfrom_plan:dblarr(ndays)}, nfiles)
    
    sec_per_day = 24.0d * 3600.0d
    for fn=0L,nfiles-1L do begin
        data[fn].file = files[fn]
        restore,files[fn]
        for dn=0L,ndays-1L do begin
            p=where(sim_exp.endtime ge gpsdays[0] and sim_exp.endtime lt gpsdays[dn+1])
            data[fn].day[dn] = sdays[dn+1]
            data[fn].totalexposure[dn] = TOTAL(sim_exp[p].solar_exp) / sec_per_day
            data[fn].totalgap_sci[dn] = TOTAL(sim_exp[p].gap_sci_tlm) / sec_per_day
            data[fn].totalgap_sc[dn] = TOTAL(sim_exp[p].gap_sc_tlm) / sec_per_day
            data[fn].totalgap_hk[dn] = TOTAL(sim_exp[p].gap_hk_tlm) / sec_per_day
            data[fn].totalgap_hrt[dn] = TOTAL(sim_exp[p].gap_hrt_tlm) / sec_per_day
            data[fn].totalgap_abs[dn] = TOTAL(sim_exp[p].abs_gap) / sec_per_day
            data[fn].totalfrom_plan[dn] = TOTAL(sim_exp[p].hrt_from_plan) / sec_per_day
        endfor
    endfor

    if n_elements(outfile) eq 0 and keyword_set(noprint) then return

    ; print the data to file and/or screen
    header=$
'     DAY   TotalExposure   TotalGap_Sci TotalGap_Hk   TotalGap_SC  TotalGap_HRT TotalGap_ABS TotalFrom_Plan'

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
                data[fn].totalgap_sci[dn], data[fn].totalgap_hk[dn], $
                data[fn].totalgap_sc[dn], data[fn].totalgap_hrt[dn], $
                data[fn].totalgap_abs[dn], data[fn].totalfrom_plan[dn]
            endfor
            printf,unit,''
            printf,unit,''
        endfor
        free_lun, unit
    endif

    if keyword_set(noprint) then return

    for fn=0L,nfiles-1L do begin
        print,data[fn].file
        print, header
        for dn=0L,ndays-1L do begin
            print,format='(F10.2,7F13.4)',data[fn].day[dn], data[fn].totalexposure[dn], $
                data[fn].totalgap_sci[dn], data[fn].totalgap_hk[dn], $
                data[fn].totalgap_sc[dn], data[fn].totalgap_hrt[dn], $
                data[fn].totalgap_abs[dn], data[fn].totalfrom_plan[dn]
        endfor
        print,''
        print,''
    endfor

    return
end
