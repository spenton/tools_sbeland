;+
; NAME:   INSERT_SIM_PLANS
;
; PURPOSE: 
;    This routine either reads in a list of IDL savefiles or uses an input 
;    data structure to insert the Solar Exposure Times in the database
;    table dbo.SimSolarExposureData.
;    The savefiles or data structure expected comes from cal_sim_exp.pro
;    and is re-formatted here to match the database table definition.
;
; CALLING SEQUENCE:
;    insert_exp_sav, indata=sim_exp, /dbinsert
;
; OPTIONAL INPUT PARAMETERS:
;    fnames - list of IDL savefiles to process
;
;    indata - the data structure to process
;
; OPTIONAL INPUT KEYWORDS:
;    insert - if set will insert the re-formatted data in the database
;
; OUTPUT PARAMETERS:
;   none
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
;    The table SimSolarExposureData has a UNIQUE constraint on:
;    instrument, orbitStartTimeGps and fieldOfViewWidth
;    The routine will fail if it tries to insert a row that already exists.
;
; REVISION HISTORY:
;   2012.03.12  SBeland
;   Revision: $Id: insert_exp_sav.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;-
;*****************************************************************************
PRO INSERT_EXP_SAV, fnames=fnames, indata=indata, dbinsert=dbinsert
    ; insert the list of SIM Solar Exposure times in the
    ; database table SORCE.dbo.SimSolarExposureData

    if size(indata,/tname) eq 'STRUCT' then begin
        data=indata
    endif else begin
        if n_elements(fnames) eq 0 then begin
            fnames =['/Users/sbeland/SORCE/data/test17/sima_hrtout_10arc.sav', $
                     '/Users/sbeland/SORCE/data/test17/sima_hrtin_10arc.sav', $
                     '/Users/sbeland/SORCE/data/test17/sima_hrtout_FOV.sav', $
                     '/Users/sbeland/SORCE/data/test17/sima_hrtin_FOV.sav', $
                     '/Users/sbeland/SORCE/data/test17/simb_hrtout_10arc.sav', $
                     '/Users/sbeland/SORCE/data/test17/simb_hrtin_10arc.sav', $
                     '/Users/sbeland/SORCE/data/test17/simb_hrtout_FOV.sav', $
                     '/Users/sbeland/SORCE/data/test17/simb_hrtin_FOV.sav']
        endif
        ; combine all the data
        data = []
        for i=0,n_elements(fnames)-1 do begin
            if file_test(fnames[i],/read) then begin
                restore,fnames[i]
                ; we expect the stucture name to be sim_exp
                data=[data,sim_exp]
            endif
        endfor
        if size(data,/tname) ne 'STRUCT' then begin
            print,'No data found in the list of files: ',fnames
            return
        endif
    endelse

    ; sort by starttime
    s=sort(data.starttime)
    data=data[s]
    instrument=strupcase(data.instrument)
    p=where(instrument eq 'SIMA', comp=cp) 
    instrumentModeId = replicate('',n_elements(data))
    instrumentModeId[p]='54'    ; SIMA
    instrumentModeId[cp]='55'   ; SIMB

    if keyword_set(dbinsert) then insert=1 else insert=0

    ; get the 1-AU correction
    ;time54 = (data[p].starttime + data[p].stoptime)/2d
    ;time55 = (data[cp].starttime + data[cp].stoptime)/2d
    ;print,'  getting SUNOBSERVERDISTANCECORRECTION from database ...'
    ;q1="SELECT microsecondsSinceGpsEpoch, sunObserverDistanceCorrection FROM SolarDistAndDopplerFixedStep2 "+$
    ;   "WHERE instrumentModeId=50 and version=7"
    ;query_database, q1, solardist, info
    ;oneau54 = interpol(solardist.SUNOBSERVERDISTANCECORRECTION, solardist.MICROSECONDSSINCEGPSEPOCH, time54)
    ;oneau55 = interpol(solardist.SUNOBSERVERDISTANCECORRECTION, solardist.MICROSECONDSSINCEGPSEPOCH, time55)


    ; prepare the database connection
    jstmt = fjava_get_jdbc_statement(user='sorce', password='sorcedb', server='sorce-db', database='SORCE')
    sorceDbExchange = OBJ_NEW("oJava_DbExchange", jstmt)
    q1 = "INSERT INTO SORCE.dbo.SimSolarExposureData VALUES("

    ; each structure contains every orbit - sort things out to fit our database table

    ; START with the SIMA - FOV data
    k = where(instrument eq 'SIMA' and data.hrtpos eq 'OUT' and data.boxwidth eq 0.0,norbits)
    if norbits gt 0 then begin
        starttime=strtrim(string(ulong64(data[k].starttime)),2)
        endtime=strtrim(string(ulong64(data[k].endtime)),2)
        insun_time=strtrim(string(data[k].insun_time),2)
        gap_sci=strtrim(string(data[k].gap_sci),2)
        gap_sc=strtrim(string(data[k].gap_sc),2)
        gap_hrt=strtrim(string(data[k].gap_hrt),2)
        same_hrt=strtrim(string(data[k].same_hrt),2)
        power_off=strtrim(string(data[k].power_off),2)
        hrt_from_plan=strtrim(string(data[k].hrt_from_plan),2)
        safehold=strtrim(string(data[k].safehold),2)
        boxwidth=strtrim(string(data[k].boxwidth),2)
        time_from_plan=strtrim(string(data[k].time_from_plan),2)
        time_std_from_plan=strtrim(string(data[k].time_std_from_plan),2)
        abs_gap=strtrim(string(data[k].abs_gap),2)
        cumul_gaps = strtrim(string(total(data[k].abs_gap, /cumul)),2)
        solar_exp_out  = data[k].solar_exp

        solar_exp_in  = dblarr(norbits)
        p_in  = where(instrument eq 'SIMA' and data.hrtpos eq 'IN' and data.boxwidth eq 0.0,count_in)
        if count_in eq norbits then begin
            solar_exp_in = data[p_in].solar_exp
        endif 
        solar_cumul_in=strtrim(string(total(solar_exp_in, /cumul)),2)
        solar_cumul_out=strtrim(string(total(solar_exp_out, /cumul)),2)
        solar_exp_in=strtrim(string(solar_exp_in),2)
        solar_exp_out=strtrim(string(solar_exp_out),2)

        for i=0L,norbits-1L do begin
            if (i mod 100) eq 0 then print,i+1,'/',norbits
            q2 = q1+instrumentModeId[k[i]]
            q2 = q2+","+starttime[i]
            q2 = q2+","+endtime[i]
            q2 = q2+","+insun_time[i]
            q2 = q2+","+gap_sci[i]
            q2 = q2+","+gap_sc[i]
            q2 = q2+","+gap_hrt[i]
            q2 = q2+","+same_hrt[i]
            q2 = q2+","+power_off[i]
            q2 = q2+","+hrt_from_plan[i]
            q2 = q2+","+safehold[i]
            q2 = q2+","+boxwidth[i]
            q2 = q2+","+time_from_plan[i]
            q2 = q2+","+time_std_from_plan[i]
            q2 = q2+","+abs_gap[i]
            q2 = q2+","+solar_exp_in[i]
            q2 = q2+","+solar_exp_out[i]
            q2 = q2+","+solar_cumul_in[i]
            q2 = q2+","+solar_cumul_out[i]
            q2 = q2+","+cumul_gaps[i]
            ; add 0.0 for the 1AU and Dose values - added later
            q2 = q2+",0.0,0.0,0.0,0.0)"
            if insert then res=sorceDbExchange->execute(q2)
        endfor
    endif
    

    
    ; SIMB - FOV data
    k = where(instrument eq 'SIMB' and data.hrtpos eq 'OUT' and data.boxwidth eq 0.0,norbits)
    if norbits gt 0 then begin
        starttime=strtrim(string(ulong64(data[k].starttime)),2)
        endtime=strtrim(string(ulong64(data[k].endtime)),2)
        insun_time=strtrim(string(data[k].insun_time),2)
        gap_sci=strtrim(string(data[k].gap_sci),2)
        gap_sc=strtrim(string(data[k].gap_sc),2)
        gap_hrt=strtrim(string(data[k].gap_hrt),2)
        same_hrt=strtrim(string(data[k].same_hrt),2)
        power_off=strtrim(string(data[k].power_off),2)
        hrt_from_plan=strtrim(string(data[k].hrt_from_plan),2)
        safehold=strtrim(string(data[k].safehold),2)
        boxwidth=strtrim(string(data[k].boxwidth),2)
        time_from_plan=strtrim(string(data[k].time_from_plan),2)
        time_std_from_plan=strtrim(string(data[k].time_std_from_plan),2)
        abs_gap=strtrim(string(data[k].abs_gap),2)
        cumul_gaps = strtrim(string(total(data[k].abs_gap, /cumul)),2)
        solar_exp_out  = data[k].solar_exp

        solar_exp_in  = dblarr(norbits)
        p_in = where(instrument eq 'SIMB' and data.hrtpos eq 'IN' and data.boxwidth eq 0.0,count_in)
        if count_in eq norbits then begin
            solar_exp_in = data[p_in].solar_exp
        endif 
        solar_cumul_in=strtrim(string(total(solar_exp_in, /cumul)),2)
        solar_cumul_out=strtrim(string(total(solar_exp_out, /cumul)),2)
        solar_exp_in=strtrim(string(solar_exp_in),2)
        solar_exp_out=strtrim(string(solar_exp_out),2)

        for i=0L,norbits-1L do begin
            if (i mod 100) eq 0 then print,i+1,'/',norbits
            q2 = q1+instrumentModeId[k[i]]
            q2 = q2+","+starttime[i]
            q2 = q2+","+endtime[i]
            q2 = q2+","+insun_time[i]
            q2 = q2+","+gap_sci[i]
            q2 = q2+","+gap_sc[i]
            q2 = q2+","+gap_hrt[i]
            q2 = q2+","+same_hrt[i]
            q2 = q2+","+power_off[i]
            q2 = q2+","+hrt_from_plan[i]
            q2 = q2+","+safehold[i]
            q2 = q2+","+boxwidth[i]
            q2 = q2+","+time_from_plan[i]
            q2 = q2+","+time_std_from_plan[i]
            q2 = q2+","+abs_gap[i]
            q2 = q2+","+solar_exp_in[i]
            q2 = q2+","+solar_exp_out[i]
            q2 = q2+","+solar_cumul_in[i]
            q2 = q2+","+solar_cumul_out[i]
            q2 = q2+","+cumul_gaps[i]
            ; add 0.0 for the 1AU and Dose values - added later
            q2 = q2+",0.0,0.0,0.0,0.0)"
            if insert then res=sorceDbExchange->execute(q2)
        endfor
    endif



    ; START with the SIMA - 10arcsec data
    k = where(instrument eq 'SIMA' and data.hrtpos eq 'OUT' and data.boxwidth eq 20.0,norbits)
    if norbits gt 0 then begin
        starttime=strtrim(string(ulong64(data[k].starttime)),2)
        endtime=strtrim(string(ulong64(data[k].endtime)),2)
        insun_time=strtrim(string(data[k].insun_time),2)
        gap_sci=strtrim(string(data[k].gap_sci),2)
        gap_sc=strtrim(string(data[k].gap_sc),2)
        gap_hrt=strtrim(string(data[k].gap_hrt),2)
        same_hrt=strtrim(string(data[k].same_hrt),2)
        power_off=strtrim(string(data[k].power_off),2)
        hrt_from_plan=strtrim(string(data[k].hrt_from_plan),2)
        safehold=strtrim(string(data[k].safehold),2)
        boxwidth=strtrim(string(data[k].boxwidth),2)
        time_from_plan=strtrim(string(data[k].time_from_plan),2)
        time_std_from_plan=strtrim(string(data[k].time_std_from_plan),2)
        abs_gap=strtrim(string(data[k].abs_gap),2)
        cumul_gaps = strtrim(string(total(data[k].abs_gap, /cumul)),2)
        solar_exp_out  = data[k].solar_exp

        solar_exp_in  = dblarr(norbits)
        p_in = where(instrument eq 'SIMA' and data.hrtpos eq 'IN' and data.boxwidth eq 20.0,count_in)
        if count_in eq norbits then begin
            solar_exp_in = data[p_in].solar_exp
        endif 
        solar_cumul_in=strtrim(string(total(solar_exp_in, /cumul)),2)
        solar_cumul_out=strtrim(string(total(solar_exp_out, /cumul)),2)
        solar_exp_in=strtrim(string(solar_exp_in),2)
        solar_exp_out=strtrim(string(solar_exp_out),2)

        for i=0L,norbits-1L do begin
            if (i mod 100) eq 0 then print,i+1,'/',norbits
            q2 = q1+instrumentModeId[k[i]]
            q2 = q2+","+starttime[i]
            q2 = q2+","+endtime[i]
            q2 = q2+","+insun_time[i]
            q2 = q2+","+gap_sci[i]
            q2 = q2+","+gap_sc[i]
            q2 = q2+","+gap_hrt[i]
            q2 = q2+","+same_hrt[i]
            q2 = q2+","+power_off[i]
            q2 = q2+","+hrt_from_plan[i]
            q2 = q2+","+safehold[i]
            q2 = q2+","+boxwidth[i]
            q2 = q2+","+time_from_plan[i]
            q2 = q2+","+time_std_from_plan[i]
            q2 = q2+","+abs_gap[i]
            q2 = q2+","+solar_exp_in[i]
            q2 = q2+","+solar_exp_out[i]
            q2 = q2+","+solar_cumul_in[i]
            q2 = q2+","+solar_cumul_out[i]
            q2 = q2+","+cumul_gaps[i]
            ; add 0.0 for the 1AU and Dose values - added later
            q2 = q2+",0.0,0.0,0.0,0.0)"
            if insert then res=sorceDbExchange->execute(q2)
        endfor
    endif

    
    ; SIMB - 10arcsec data
    k = where(instrument eq 'SIMB' and data.hrtpos eq 'OUT' and data.boxwidth eq 20.0,norbits)
    if norbits gt 0 then begin
        starttime=strtrim(string(ulong64(data[k].starttime)),2)
        endtime=strtrim(string(ulong64(data[k].endtime)),2)
        insun_time=strtrim(string(data[k].insun_time),2)
        gap_sci=strtrim(string(data[k].gap_sci),2)
        gap_sc=strtrim(string(data[k].gap_sc),2)
        gap_hrt=strtrim(string(data[k].gap_hrt),2)
        same_hrt=strtrim(string(data[k].same_hrt),2)
        power_off=strtrim(string(data[k].power_off),2)
        hrt_from_plan=strtrim(string(data[k].hrt_from_plan),2)
        safehold=strtrim(string(data[k].safehold),2)
        boxwidth=strtrim(string(data[k].boxwidth),2)
        time_from_plan=strtrim(string(data[k].time_from_plan),2)
        time_std_from_plan=strtrim(string(data[k].time_std_from_plan),2)
        abs_gap=strtrim(string(data[k].abs_gap),2)
        cumul_gaps = strtrim(string(total(data[k].abs_gap, /cumul)),2)
        solar_exp_out  = data[k].solar_exp

        solar_exp_in  = dblarr(norbits)
        p_in = where(instrument eq 'SIMB' and data.hrtpos eq 'IN' and data.boxwidth eq 0.0,count_in)
        if count_in eq norbits then begin
            solar_exp_in = data[p_in].solar_exp
        endif 
        solar_cumul_in=strtrim(string(total(solar_exp_in, /cumul)),2)
        solar_cumul_out=strtrim(string(total(solar_exp_out, /cumul)),2)
        solar_exp_in=strtrim(string(solar_exp_in),2)
        solar_exp_out=strtrim(string(solar_exp_out),2)

        for i=0L,norbits-1L do begin
            if (i mod 100) eq 0 then print,i+1,'/',norbits
            q2 = q1+instrumentModeId[k[i]]
            q2 = q2+","+starttime[i]
            q2 = q2+","+endtime[i]
            q2 = q2+","+insun_time[i]
            q2 = q2+","+gap_sci[i]
            q2 = q2+","+gap_sc[i]
            q2 = q2+","+gap_hrt[i]
            q2 = q2+","+same_hrt[i]
            q2 = q2+","+power_off[i]
            q2 = q2+","+hrt_from_plan[i]
            q2 = q2+","+safehold[i]
            q2 = q2+","+boxwidth[i]
            q2 = q2+","+time_from_plan[i]
            q2 = q2+","+time_std_from_plan[i]
            q2 = q2+","+abs_gap[i]
            q2 = q2+","+solar_exp_in[i]
            q2 = q2+","+solar_exp_out[i]
            q2 = q2+","+solar_cumul_in[i]
            q2 = q2+","+solar_cumul_out[i]
            q2 = q2+","+cumul_gaps[i]
            ; add 0.0 for the 1AU and Dose values - added later
            q2 = q2+",0.0,0.0,0.0,0.0)"
            if insert then res=sorceDbExchange->execute(q2)
        endfor
    endif



    ; START with the SIMA - all other field-of-view radius
    k = where(instrument eq 'SIMA' and data.hrtpos eq 'OUT' and data.boxwidth ne 0.0 and data.boxwidth ne 20.0,norbits)
    if norbits gt 0 then begin
        starttime=strtrim(string(ulong64(data[k].starttime)),2)
        endtime=strtrim(string(ulong64(data[k].endtime)),2)
        insun_time=strtrim(string(data[k].insun_time),2)
        gap_sci=strtrim(string(data[k].gap_sci),2)
        gap_sc=strtrim(string(data[k].gap_sc),2)
        gap_hrt=strtrim(string(data[k].gap_hrt),2)
        same_hrt=strtrim(string(data[k].same_hrt),2)
        power_off=strtrim(string(data[k].power_off),2)
        hrt_from_plan=strtrim(string(data[k].hrt_from_plan),2)
        safehold=strtrim(string(data[k].safehold),2)
        boxwidth=strtrim(string(data[k].boxwidth),2)
        time_from_plan=strtrim(string(data[k].time_from_plan),2)
        time_std_from_plan=strtrim(string(data[k].time_std_from_plan),2)
        abs_gap=strtrim(string(data[k].abs_gap),2)
        cumul_gaps = strtrim(string(total(data[k].abs_gap, /cumul)),2)
        solar_exp_out  = data[k].solar_exp

        solar_exp_in  = dblarr(norbits)
        p_in = where(instrument eq 'SIMA' and data.hrtpos eq 'IN' and data.boxwidth eq 20.0,count_in)
        if count_in eq norbits then begin
            solar_exp_in = data[p_in].solar_exp
        endif 
        solar_cumul_in=strtrim(string(total(solar_exp_in, /cumul)),2)
        solar_cumul_out=strtrim(string(total(solar_exp_out, /cumul)),2)
        solar_exp_in=strtrim(string(solar_exp_in),2)
        solar_exp_out=strtrim(string(solar_exp_out),2)

        for i=0L,norbits-1L do begin
            if (i mod 100) eq 0 then print,i+1,'/',norbits
            q2 = q1+instrumentModeId[k[i]]
            q2 = q2+","+starttime[i]
            q2 = q2+","+endtime[i]
            q2 = q2+","+insun_time[i]
            q2 = q2+","+gap_sci[i]
            q2 = q2+","+gap_sc[i]
            q2 = q2+","+gap_hrt[i]
            q2 = q2+","+same_hrt[i]
            q2 = q2+","+power_off[i]
            q2 = q2+","+hrt_from_plan[i]
            q2 = q2+","+safehold[i]
            q2 = q2+","+boxwidth[i]
            q2 = q2+","+time_from_plan[i]
            q2 = q2+","+time_std_from_plan[i]
            q2 = q2+","+abs_gap[i]
            q2 = q2+","+solar_exp_in[i]
            q2 = q2+","+solar_exp_out[i]
            q2 = q2+","+solar_cumul_in[i]
            q2 = q2+","+solar_cumul_out[i]
            q2 = q2+","+cumul_gaps[i]
            ; add 0.0 for the 1AU and Dose values - added later
            q2 = q2+",0.0,0.0,0.0,0.0)"
            if insert then res=sorceDbExchange->execute(q2)
        endfor
    endif

    
    ; SIMB - all other field-of-view radius
    k = where(instrument eq 'SIMB' and data.hrtpos eq 'OUT' and data.boxwidth ne 0.0 and data.boxwidth ne 20.0,norbits)
    if norbits gt 0 then begin
        starttime=strtrim(string(ulong64(data[k].starttime)),2)
        endtime=strtrim(string(ulong64(data[k].endtime)),2)
        insun_time=strtrim(string(data[k].insun_time),2)
        gap_sci=strtrim(string(data[k].gap_sci),2)
        gap_sc=strtrim(string(data[k].gap_sc),2)
        gap_hrt=strtrim(string(data[k].gap_hrt),2)
        same_hrt=strtrim(string(data[k].same_hrt),2)
        power_off=strtrim(string(data[k].power_off),2)
        hrt_from_plan=strtrim(string(data[k].hrt_from_plan),2)
        safehold=strtrim(string(data[k].safehold),2)
        boxwidth=strtrim(string(data[k].boxwidth),2)
        time_from_plan=strtrim(string(data[k].time_from_plan),2)
        time_std_from_plan=strtrim(string(data[k].time_std_from_plan),2)
        abs_gap=strtrim(string(data[k].abs_gap),2)
        cumul_gaps = strtrim(string(total(data[k].abs_gap, /cumul)),2)
        solar_exp_out  = data[k].solar_exp

        solar_exp_in  = dblarr(norbits)
        p_in = where(instrument eq 'SIMB' and data.hrtpos eq 'IN' and data.boxwidth eq 0.0,count_in)
        if count_in eq norbits then begin
            solar_exp_in = data[p_in].solar_exp
        endif 
        solar_cumul_in=strtrim(string(total(solar_exp_in, /cumul)),2)
        solar_cumul_out=strtrim(string(total(solar_exp_out, /cumul)),2)
        solar_exp_in=strtrim(string(solar_exp_in),2)
        solar_exp_out=strtrim(string(solar_exp_out),2)

        for i=0L,norbits-1L do begin
            if (i mod 100) eq 0 then print,i+1,'/',norbits
            q2 = q1+instrumentModeId[k[i]]
            q2 = q2+","+starttime[i]
            q2 = q2+","+endtime[i]
            q2 = q2+","+insun_time[i]
            q2 = q2+","+gap_sci[i]
            q2 = q2+","+gap_sc[i]
            q2 = q2+","+gap_hrt[i]
            q2 = q2+","+same_hrt[i]
            q2 = q2+","+power_off[i]
            q2 = q2+","+hrt_from_plan[i]
            q2 = q2+","+safehold[i]
            q2 = q2+","+boxwidth[i]
            q2 = q2+","+time_from_plan[i]
            q2 = q2+","+time_std_from_plan[i]
            q2 = q2+","+abs_gap[i]
            q2 = q2+","+solar_exp_in[i]
            q2 = q2+","+solar_exp_out[i]
            q2 = q2+","+solar_cumul_in[i]
            q2 = q2+","+solar_cumul_out[i]
            q2 = q2+","+cumul_gaps[i]
            ; add 0.0 for the 1AU and Dose values - added later
            q2 = q2+",0.0,0.0,0.0,0.0)"
            if insert then res=sorceDbExchange->execute(q2)
        endfor
    endif

END

