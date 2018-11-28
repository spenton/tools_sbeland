function marty_19_20, instrumentModeId, indata=indata, waves=waves, gettsi=gettsi, tsi19=tsi19, tsi20=tsi20, $
    mincount=mincount, smooth=smooth

    if keyword_set(smooth) then noplot=1 else noplot=0
    ; compare the data from version 19 and 20 of sim for Marty's talk oct. 2014
    if size(indata,/tname)  ne 'STRUCT' then begin
       ; xtract the data for version 19
       print,'Extracting version 19 from database ...'
       res=compare_19_20(0,4200,instrumentModeId,/mission,v19=19,v20=2021,alldata=data19,/align, noplot=noplot)
 
       ; xtract the data for version 20
       print,'Extracting version 20 from database ...'
       res=compare_19_20(0,4200,instrumentModeId,/mission,v19=20,v20=2021,alldata=data20,/align, noplot=noplot)
       indata={spect19:data19.spect19, spect20:data20.spect19}
    endif
 
    tmp=label_date(date_format='%M %Y')
    interval=730.5d
    if instrumentModeId eq 41 then inst='VIS'
    if instrumentModeId eq 43 then inst='UV'
    if instrumentModeId eq 44 then inst='IR'

    if keyword_set(gettsi) and n_elements(waves) eq 2 then begin
        tsi19=get_tsi(0,1,/sima,wrange=waves,alldata=indata.spect19,version=19,/align,/noplot,mincount=mincount)
        tsi20=get_tsi(0,1,/sima,wrange=waves,alldata=indata.spect20,version=20,/align,/noplot,mincount=mincount)
        if keyword_set(smooth) then begin
            sd19=gps2sd(tsi19.timestamp/1d6)
            c19=robust_poly_fit(sd19, tsi19.tsi, 3, yfit)
            resistant_mean, (yfit-tsi19.tsi), 5.0, mean, good=k19
            tsi19={timestamp:tsi19.timestamp[k19], tsi:tsi19.tsi[k19]}
  
            sd20=gps2sd(tsi20.timestamp/1d6)
            c20=robust_poly_fit(sd20, tsi20.tsi, 3, yfit)
            resistant_mean, (yfit-tsi20.tsi), 5.0, mean, good=k20
            tsi20={timestamp:tsi20.timestamp[k20], tsi:tsi20.tsi[k20]}
        endif
        title='SORCE SIM Integrated SSI from '+strtrim(string(waves[0],format='(f0.1)'),2)
        title=title+' to '+strtrim(string(waves[1],format='(f0.1)'),2)+' nm'
        plot_multi,gps2jd(tsi19.timestamp/1d6),smooth(tsi19.tsi,2), gps2jd(tsi20.timestamp/1d6), smooth(tsi20.tsi,2), psym=[-3,-3],$
            charsize=1.4, xtickinterval=interval, xminor=24, xtickformat='LABEL_DATE',xtitle='Date',ytitle='Irradiance',title=title,$
            label=['Version 19','Version 20'],/xst,/yst
    endif
 
    for i=0,n_elements(waves)-1 do begin
        res=compare_19_20(0,1,/mission,instrumentModeId,v19=19,v20=20,alldata=indata,plotwave=waves[i],/align, noplot=noplot)
        if n_elements(res.timestamp19) lt 10 then break
        if keyword_set(smooth) then begin
            c19=robust_poly_fit(res.timestamp19, res.v19_irrad, 3, yfit)
            resistant_mean, (yfit-res.v19_irrad), 5.0, mean, good=k19
  
            c20=robust_poly_fit(res.timestamp20, res.v20_irrad, 3, yfit)
            resistant_mean, (yfit-res.v20_irrad), 5.0, mean, good=k20
  
            res={plotwave:res.plotwave, v19_wave:res.v19_wave[k19], v20_wave:res.v20_wave[k20], $
                 v19_irrad:res.v19_irrad[k19], v20_irrad:res.v20_irrad[k20], $
                 timestamp19:res.timestamp19[k19], timestamp20:res.timestamp20[k20]}
        endif 
        title='SIM '+inst+' version=19 wl='+strtrim(string(waves[i],format='(f0.1)'),2)
        lineplot,sd2jd(res.timestamp19),res.v19_irrad,title=title,nsum=2, xtitle='Date',ytitle='Irradiance',$
            ptitle='SORCE SIM Version 19 & 20',charsize=1.4, xtickinterval=interval, xminor=24, xtickformat='LABEL_DATE'

        title='SIM '+inst+' version=20 wl='+strtrim(string(waves[i],format='(f0.1)'),2)
        lineplot,sd2jd(res.timestamp20),res.v20_irrad,title=title,nsum=2, xtitle='Date',ytitle='Irradiance',$
            ptitle='SORCE SIM Version 19 & 20',charsize=1.4, xtickinterval=interval, xminor=24, xtickformat='LABEL_DATE'
    endfor
 
    return,res

end
