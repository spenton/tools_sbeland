; extract the CCD Dump data and the prism position for the corresponding time ranges
;
FUNCTION CCDDUMP, startTime, stopTime, simA=simA, simB=simB, gps=gps, missionDays=missionDays, $
             julianDays=julianDays, activity=activity



    if keyword_set(julianDays) then begin
        ; user specified time in julian days
        startTime = jd2sd(startTime)
        stopTime = jd2sd(stopTime)
    endif else if keyword_set(gps) then begin
        ; user specified timetags in gps microseconds
        startTime = gps2sd(startTime/1d6)
        stopTime = gps2sd(stopTime/1d6)
    endif 

    ; get each instances of the CCDDump during this time span
    if n_elements(activity) eq 0 then activity='CCDDumpImageLight'
    plans = get_sorce(startTime, stopTime, /missionDays, simA=simA, simB=simB, activity=activity)

    if size(plans,/tname) ne 'STRUCT' then begin
        print,'No matching planned activities for this time range'
        return,0
    endif

    if not (keyword_set(simA) or keyword_set(simB)) then begin
        ; if neither is specified default to simA
        simA=1
        simB=0
    endif else if (keyword_set(simA) and keyword_set(simB)) then begin
        ; if both are specified default to simA
        simA=1
        simB=0
    endif
    if n_elements(simB) eq 0 then simB=0

    ; the telemetry Id number corresponding to the 
    ;    prism position (CCD subpixel)
    ;    CCD integration time 
    ;    prism drive temperature
    if simB eq 1 then begin
        tlmId=[100478L, 101647L, 102212L]
    endif else begin
        tlmId=[101862L, 101395L, 101909L]
    endelse

    tmpstr = {startTime:0.0d, stopTime:0.0d, prism_position:0.0d, prism_stdev:0.0d, $
        prism_min:0.d, prism_max:0.0d, fit_center:0.0d, fit_fwhm:0.0d, fit_height:0.0d, $
        intg_period:0.0d, prism_drive_temp:0.0d, saturationWidth:0.0d, zeroCrossing:0.0d, data:ptr_new()}
    out_data=[]

    for pl=0,n_elements(plans)-1 do begin

        t0 = sd2gps(plans[pl].starttime) * 1d6
        t1 = sd2gps(plans[pl].stoptime) * 1d6
        ; get the prism position for the channel of interest
        get_sorce_telemetry, data, info, t0, t1, tlmId=tlmId, /gps
        if (*data.(0)).n_science eq 0 then continue
        ; trim 10% of times at both ends
        delta_t = ((*data.(0)).science[-1].timetag - (*data.(0)).science[0].timetag)*0.10
        t0 = (*data.(0)).science[0].timetag + delta_t
        t1 = (*data.(0)).science[-1].timetag - delta_t
        p = where((*data.(0)).science.timetag ge t0 and (*data.(0)).science.timetag le t1,count)
        if count eq 0 then continue
        position = mean((*data.(0)).science[p].dn)
        std   = sqrt(variance((*data.(0)).science[p].dn))
        mn=min((*data.(0)).science[p].dn,max=mx)
        if mn eq mx then std=0.0

        if (*data.(1)).n_housekeeping gt 0 then tmpstr.intg_period = median((*data.(1)).housekeeping.eu)
        if (*data.(2)).n_housekeeping gt 0 then tmpstr.prism_drive_temp = median((*data.(2)).housekeeping.eu)

        if strpos(activity,'CCDDump') lt 0 then begin
            ; no CCDDump data to extract - just the prism position
            center=0.0
            fwhm=0.0
            height=0.0
            tmpstr.startTime=plans[pl].starttime
            tmpstr.stopTime=plans[pl].stoptime
            tmpstr.prism_position=position
            tmpstr.prism_stdev=std
            tmpstr.prism_min=mn
            tmpstr.prism_max=mx
            tmpstr.fit_center=0.0
            tmpstr.fit_fwhm=0.0
            tmpstr.fit_height=0.0
            out_data=[out_data,tmpstr]
            print,plans[pl].starttime, plans[pl].stoptime,position,std,mn,mx,center*10.0,fwhm*10.0
            continue
        endif

        ; get the CCDdump image
        if simB eq 1 then begin
            res=get_sampled_ccd_dump(plans[pl].starttime, plans[pl].stoptime, /sim_b, type=1)
            if n_elements(res) eq 1 then $
                res=get_sampled_ccd_dump(plans[pl].starttime, plans[pl].stoptime, /sim_b, type=2)
        endif else begin
            res=get_sampled_ccd_dump(plans[pl].starttime, plans[pl].stoptime, /sim_a, type=1)
            if n_elements(res) eq 1 then $
                res=get_sampled_ccd_dump(plans[pl].starttime, plans[pl].stoptime, /sim_a, type=2)
        endelse

        ; find the center of the image
        if n_elements(res) eq 1 then begin
            center=0.0
            fwhm=0.0
            height=0.0
            print,plans[pl].starttime, plans[pl].stoptime,position,std,mn,mx,center*10.0,fwhm*10.0
            tmpstr.startTime=plans[pl].starttime
            tmpstr.stopTime=plans[pl].stoptime
            tmpstr.prism_position=position
            tmpstr.prism_stdev=std
            tmpstr.prism_min=mn
            tmpstr.prism_max=mx
            tmpstr.fit_center=center*10.0
            tmpstr.fit_fwhm=fwhm*10.0
            out_data=[out_data,tmpstr]
            continue
        endif 

        for dim=0,n_elements(res[*,0])-1 do begin
            xx=lindgen(n_elements(res[dim,*]))
            fit = gaussfit(xx,res[dim,*],coef,nterms=4)
            center = coef[1]
            height=coef[0]
            fwhm = 2.0d*SQRT(2.0d*ALOG(2.0d))*coef[2]
            ; refine the fit 
            p0 = long(center - abs(fwhm * 15.0d)/2)
            p1 = p0 + long(abs(fwhm) * 15)
            fit = gaussfit(xx[p0:p1],res[dim,p0:p1],coef,nterms=4)
            center = coef[1]
            height=coef[0]
            fwhm = 2.0d*SQRT(2.0d*ALOG(2.0d))*coef[2]
            print,plans[pl].starttime, plans[pl].stoptime,position,std,mn,mx,center*10.0,fwhm*10.0,height
            tmpstr.startTime=plans[pl].starttime
            tmpstr.stopTime=plans[pl].stoptime
            tmpstr.prism_position=position
            tmpstr.prism_stdev=std
            tmpstr.prism_min=mn
            tmpstr.prism_max=mx
            tmpstr.fit_center=center*10.0
            tmpstr.fit_fwhm=fwhm*10.0
            tmpstr.fit_height=height
            tmpstr.data = ptr_new(res[dim,*])
            p=where(res[dim,*] gt 8.0d5,count)
            if count gt 0 then tmpstr.saturationWidth = xx[p[-1]] - xx[p[0]] + 1
            x=dindgen(n_elements(res[dim,*]))
            drv = deriv(x,res[dim,*])
            ; p=where(drv[100:-2] gt drv[101:-1] and drv[101:-1] lt 0.0 and drv[100:-2] gt 5000)
            p=where(drv[100:-2] gt drv[101:-1] and drv[101:-1] lt 0.0 and res[dim,100:-2] gt 8.0d5)
            if count gt 0 then begin
                ; calculate the first zero-crossing of the derivative
                p=p[0]+100L
                tmpstr.zeroCrossing=P+1 - drv[p+1] / (drv[p+1]-drv[p])
            endif
            out_data=[out_data,tmpstr]
        endfor

    endfor

    return,out_data

end
