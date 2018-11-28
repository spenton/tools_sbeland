pro compare_jh_sb, jhfile=jhfile, sbfile=sbfile
    ; this routine simply reads the ascii file from Jerry Harder
    ; with the solar exposure results and plot them against the
    ; data from SBeland.
    ; The channel is verified from the file name.

    if n_elements(jhfile) eq 0 then begin
       jhfile = '~/SORCE/code/tools_harder/simAexposureTime_100429.txt'
       ; jhfile = '~/SORCE/code/tools_harder/simBexposureTime_100429.txt'
       ; jhfile=dialog_pickfile(title='Pick the file from Jerry',/must_exist,/read)
       if jhfile eq '' then begin
           print,'File selection was canceled'
           return
       endif
    endif

    if n_elements(sbfile) eq 0 then begin
       ; sbfile = '~/SORCE/data/sima_solar_exp_hrtout_fov_sd.txt'
       ; sbfile = '~/SORCE/data/simb_solar_exp_hrtout_fov_sd.txt'
       sbfile=dialog_pickfile(title='Pick the file from Stephane',/must_exist,/read)
       if sbfile eq '' then begin
           print,'File selection was canceled'
           return
       endif
    endif

    if strpos(strupcase(jhfile),'SIMA') ge 0 then begin
        jh_inst='SIMA'
    endif else if strpos(strupcase(jhfile),'SIMB') ge 0 then begin
        jh_inst='SIMB'
    endif else begin
        jh_inst = 'NA'
    endelse

    if strpos(strupcase(sbfile),'SIMA') ge 0 then begin
        sb_inst='SIMA'
    endif else if strpos(strupcase(sbfile),'SIMB') ge 0 then begin
        sb_inst='SIMB'
    endif else begin
        sb_inst = 'NA'
    endelse

    if jh_inst eq 'NA' or sb_inst eq 'NA' or (jh_inst ne sb_inst) then begin
        print,'JH_File and SB_File are for different channels'
        return
    endif

    readcol, jhfile, jh_sd, jh_jd, jh_gps,jh_exp
    readcol, sbfile, sb_sd, sb_jd, sb_gps, sb_exp

    mx=max([min(jh_sd),min(sb_sd)])
    p=where(jh_sd ge mx)
    jh_sd=jh_sd[p]
    jh_jd=jh_jd[p]
    jh_gps=jh_gps[p]
    jh_exp=jh_exp[p]

    p=where(sb_sd ge mx)
    sb_sd=sb_sd[p]
    sb_jd=sb_jd[p]
    sb_gps=sb_gps[p]
    sb_exp=sb_exp[p]

    diff_exp = sb_exp - jh_exp
    lineplot,jh_sd,jh_exp, title=jh_inst+' - JH',ptitle=jh_inst+' Cumulative Solar Exposure - SB - JH',$
        xtitle='Mission Days',ytitle='Cumulative Solar Exposure (days)',psym=-3
    lineplot,sb_sd,sb_exp, title=sb_inst+' - SB',psym=-3
    lineplot,sb_sd[0:n_elements(diff_exp)-1],diff_exp, title=sb_inst+' - Difference (SB-JH)',psym=-3

    return
end
