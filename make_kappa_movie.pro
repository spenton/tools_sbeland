pro make_kappa_movie, outfile, all_kappa, kappa_smooth, by_days=by_days

    ;wid=10
    width=1000
    height=600
    ;window,wid,xsize=width,ysize=height

    if ~keyword_set(by_days) then begin
        ;limit the wavelengths from 210 to 306 nm
        waves = kappa_smooth.waves
        wp=where(waves ge 210.0 and waves le 307.0, nwaves)

        ; set the duration to 20 seconds
        fps = nwaves / 20

        all_timesd = dblarr(n_elements(all_kappa))
        for i=0L, n_elements(all_kappa)-1L do all_timesd[i] = gps2sd((*all_kappa[i]).time/1d6)

        kappa_fit=dblarr(n_elements(all_timesd), n_elements(kappa_smooth.waves))
        for d=0L,n_elements(all_timesd)-1L do kappa_fit[d,*]=(*all_kappa[d]).kappa_fit

        ; find matching times
        ;match, kappa_smooth.timesd, all_timesd, sa, sb

        xrange=[0d,6000]
        yrange=[0.001,0.008]

        ; the file extension determines the type of output file (mp4 is probably best)
        oVid = IDLffVideoWrite(outfile)
        vidStream = oVid.AddVideoStream(width, height, fps)

        set_plot, 'z', /copy
        device, set_resolution=[width,height], set_pixel_depth=24, decomposed=0

        for w=0L,nwaves-1L do begin
            plot,all_timesd,kappa_fit[*,wp[w]],yrange=yrange,xrange=xrange,/xst,/yst,psym=-4, xtitle='!3Mission Days!X',$
                ytitle='!3Kappa!X',title='!3SORCE SIM Kappa!X',charsize=2.0, font=-1
            oplot,kappa_smooth.timesd,kappa_smooth.kappa[*,wp[w]],color=2,thick=2.0
            ;oplot,timesd[ctd],yfit2,color=2
            xyouts,0.70,0.70,'!3'+strtrim(string(waves[wp[w]],format='(F0.2," nm")'),2)+'!X',/normal,color=2,charsize=2.0, font=-1
            ;wait, 0.05
            timestamp = oVid.put(vidStream, tvrd(true=1))
        endfor

    endif else begin
        ndays = n_elements(all_kappa)
        fps = ndays / 20

        ; the file extension determines the type of output file (mp4 is probably best)
        oVid = IDLffVideoWrite(outfile)
        vidStream = oVid.AddVideoStream(width, height, fps)

        set_plot, 'z', /copy
        device, set_resolution=[width,height], set_pixel_depth=24, decomposed=0
        xrange=[210d, 307d]
        yrange=[0.001,0.008]

        for d=0L,ndays-1L do begin
            plot,(*all_kappa[d]).wavelength,(*all_kappa[d]).kappa,yrange=yrange,xrange=xrange,/xst,/yst,psym=-4, $
                xtitle='!3Wavelength (nm)!X',ytitle='!3Kappa!X',charsize=2.0, font=-1
            oplot,(*all_kappa[d]).wavelength,(*all_kappa[d]).kappa_fit,color=2, thick=2.0
            xyouts,0.70,0.70,"!3"+string(gps2sd((*all_kappa[d]).time/1d6),format='("SD = ",F0.2)')+"!X",$
                /normal,color=2,charsize=2.0, font=-1
            timestamp = oVid.put(vidStream, tvrd(true=1))
       endfor

    endelse

    device,/close
    set_plot, 'X'
    oVid.cleanup
    print, 'File "' + outfile + '" written to current directory.'

end
