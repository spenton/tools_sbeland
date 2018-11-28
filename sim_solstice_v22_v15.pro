pro sim_solstice_v22_v15, waves, sim_22, sol_15, smooth=smooth

    device,decompose=1
    tmp=label_date(date_format='%D-%M!C%Y')
    interval=730.5d
    nwv=n_elements(waves)
    !p.multi=[nwv+1,1,nwv+1]
    jd0=2452640.5d   ; 1-Jan-2003 
    jd1=2457753.5d   ; 31-Dec-2016

    colorlist = ['red', 'green', 'blue', 'orange', 'yellow', $
                'magenta', 'cyan', 'dark green', 'sienna', $
                'purple', 'dark red', 'burlywood', 'lime green', $
                'maroon', 'tomato', 'slate blue', 'dark goldenrod' ,'olive', $
                'dodger blue', 'dark grey']
    colors = cgcolor(colorlist)
    !p.background=cgcolor('white')
    window,xsize=1024,ysize=1024

    for w=0L,nwv-1L do begin

        print,' Getting data for ',waves[w],' ...'
            sim_22 = get_level3_spectrum(0,5000,43,24,version=22,min_wavelength=waves[w]-0.5,max_wavelength=waves[w]+0.5)
            wv=[sim_22.minwavelength,sim_22.maxwavelength]
            mn=min(abs(wv - waves[w]),pos)
            p=where(sim_22.minwavelength eq wv[pos])
            sim_22=sim_22[p]
            if keyword_set(smooth) then begin
                t22 = gps2sd(sim_22.NOMINALGPSTIMETAG/1d6)
                c22=robust_poly_fit(t22, sim_22.irradiance, 3, yfit)
                resistant_mean, (yfit-sim_22.irradiance), 5.0, mean, good=k22
                sim_22 = sim_22[k22]
            endif

            sol_15 = get_level3_spectrum(0,5000,9,24,version=15,min_wavelength=waves[w]-1.0,max_wavelength=waves[w]+1.0)
            wv=[sol_15.minwavelength,sol_15.maxwavelength]
            mn=min(abs(wv - waves[w]),pos)
            p=where(sol_15.minwavelength eq wv[pos])
            sol_15=sol_15[p]
            if keyword_set(smooth) then begin
                t15 = gps2sd(sol_15.NOMINALGPSTIMETAG/1d6)
                c15=robust_poly_fit(t15, sol_15.irradiance, 3, yfit)
                resistant_mean, (yfit-sol_15.irradiance), 5.0, mean, good=k15
                sol_15 = sol_15[k15]
            endif 

        labels='SIM_V22 @ '+strtrim(string(waves[w],format='(f0.1)'),2)
        labels=[labels,'SOL_V15 @ '+strtrim(string(waves[w],format='(f0.1)'),2)]
        delta = median(sim_22.irradiance) - median(sol_15.irradiance)
        sol_15_val = sol_15.irradiance+delta
        ; plot the median uncertainty from SOLSTICE
        ;md0=median(sol_15_val[-200:-1])
        ;y0=md0 - median(sol_15.irradianceuncertainty)/2d
        ;y1=md0 + median(sol_15.irradianceuncertainty)/2d
        y0=sol_15_val - sol_15.irradianceuncertainty/2d
        y1=sol_15_val + sol_15.irradianceuncertainty/2d
        mn=min([y0,y1,sim_22.irradiance],max=mx)
        dd=(mx-mn)/20d
        yrange=[mn-dd, mx+dd]

        if w eq 0 then begin
            plot,gps2jd(sim_22.NOMINALGPSTIMETAG/1d6),sim_22.irradiance, ytitle='Irradiance',/nodata, ymargin=[0,2], $
                charsize=3.0, color=cgcolor('black'),xrange=[jd0,jd1], xtickinterval=interval, xminor=24, $
                xtickname=REPLICATE(' ', 8), /xst, yrange=yrange,title='SIM_V22 vs SOLSTICE_V15'
        endif else if w lt nwv-1 then begin
            plot,gps2jd(sim_22.NOMINALGPSTIMETAG/1d6),sim_22.irradiance, ytitle='Irradiance',/nodata, ymargin=[0,0], $
                charsize=3.0, color=cgcolor('black'),xrange=[jd0,jd1], xtickinterval=interval, xminor=24, $
                xtickname=REPLICATE(' ', 8), /xst, yrange=yrange
        endif else begin
            plot,gps2jd(sim_22.NOMINALGPSTIMETAG/1d6),sim_22.irradiance, xtitle='!CDate',ytitle='Irradiance',/nodata, $
                charsize=3.0, xtickinterval=interval, xminor=24, xtickformat='LABEL_DATE',color=cgcolor('black'),xrange=[jd0,jd1],$
                /xst,yrange=yrange, ymargin=[4,0]
        endelse

        ;oplot,gps2jd(sol_15.NOMINALGPSTIMETAG/1d6),sol_15_val, psym=-3,color=colors[1]
        errplot, gps2jd(sol_15.NOMINALGPSTIMETAG/1d6),y0,y1,color=cgcolor('green'),thick=2.0
        oplot,gps2jd(sol_15.NOMINALGPSTIMETAG/1d6),sol_15_val, psym=-3,color=cgcolor('dark green')

        oplot,gps2jd(sim_22.NOMINALGPSTIMETAG/1d6),sim_22.irradiance, psym=-3,color=colors[0]

        al_legend,labels,colors=colors[0:1],psym=[-4,-4],box=0,textcolor=[cgcolor('black'),cgcolor('black')],/right,/top

    endfor

end
