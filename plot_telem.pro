;+
; NAME:  PLOT_TELEM
;
; PURPOSE: 
;   This routine extract a specific set of telemetry items over the
;   specified time range and plots the data.
;
; CALLING SEQUENCE:
;   PLOT_TELEM, starttime, stoptime, /mission, /sima
;                      
; INPUT KEYWORDS:
;   startTime -
;      The lower time range for which data will be returned.
;      May be specified as julian days, gps microseconds, or
;      SORCE mission day numbers.
;
;   stopTime -
;      The upper time range for which data will be returned
;      May be specified as julian days, gps microseconds, or
;      SORCE mission day numbers.
;
; OPTIONAL INPUT KEYWORDS:
;   simA, simB - 
;      Specifies which SIM channel to retreive.  These are mutually
;      exclusive and defaults to simA.
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
;   psfile -
;      name of file to save the postscript plot. If not provided, will
;      use plot_data routine as a GUI.
;
; OPTIONAL OUTPUT PARAMETERS:
;   outdata - 
;      structure containing the resturned data
;
; EXAMPLE:  
;
; COMMON BLOCKS:
;      NONE
;
; REVISION HISTORY:
;   2013.12.17  SBeland
;   Revision: $Id
;-
;*****************************************************************************
pro GEN_PLOT_TELEM, xval, yval, xrange=xrange, colors=colors, labels=labels, psfile=psfile, resolution=resolution

    tek_color,0,25
    nlines=n_elements(xval)

    colorlist = ['red', 'green', 'blue', 'orange',  $
              'magenta', 'cyan', 'dark green', 'sienna', $
              'purple', 'dark red', 'burlywood', 'lime green', $
              'maroon', 'tomato', 'slate blue', 'dark goldenrod' ,'olive', $
              'dodger blue', 'dark grey',$
              'red', 'green', 'blue', 'orange',$
              'magenta', 'cyan', 'dark green', 'sienna', $
              'purple', 'dark red', 'burlywood', 'lime green']
    colors = cgcolor(colorlist)

    nlabels = n_elements(labels)
    delta_t= xrange[1]-xrange[0]
    if delta_t le 2.1 then begin
        res=label_date(date_format='%H:%I:%S', /round_up)
    endif else if delta_t le 7.1 then begin
        res=label_date(date_format='%M %D %H:%I', /round_up)
    endif else if delta_t le 30.1 then begin
        res=label_date(date_format='%M %D', /round_up)
    endif else begin
        res=label_date(date_format='%M %D %Y', /round_up)
    endelse

    lbl='  ('+string(xrange[0],format='(C(CMOI2.2, "/", CDI2.2, "/", CYI))')+')'

    if n_elements(psfile) eq 0 then begin
        ; use the plot function to display on the screen
        !P.multi=[0,1,nlines]
        for line=0,nlines-1 do begin
            if not ptr_valid(xval[line]) then continue
            if nlabels gt 0 then begin
                myplot=plot(*xval[line],*yval[line], symbol='+', /current, color=colorlist[line], $
                      layout=[1,nlines,line+1],/xstyle, /ystyle, sym_size=0.40, title=labels[line]+lbl, xrange=xrange, $
                      font_size=6, margin=[0.05,0.20,0.05,0.20], dimensions=[800,1200], xtickformat='LABEL_DATE')
            endif else begin
                myplot=plot(*xval[line],*yval[line], symbol='+', /current, color=colorlist[line], $
                      layout=[1,nlines,line+1],/xstyle, /ystyle, sym_size=0.40,  xrange=xrange, $
                      font_size=6, margin=[0.05,0.02,0.05,0.05], dimensions=[800,1200], xtickformat='LABEL_DATE')
            endelse
        endfor

    endif else if strpos(psfile,'.ps') gt 0 then begin
        ; use the plot routine to save a postscript file since the "save" method for IDLgr is soooo sloooow
        thisDevice = !d.name
        set_plot,'ps'
        device,/port,xoff=0.05,yoff=-5.1,xsize=8,ysize=16,bits=4, /color,/inches,file=psfile
        !p.background=255
        !p.color=0
        !P.multi=[0,1,nlines+5]  ; the +5 is to make all the plots fit as expected (with the yoff=-5)
        !p.font = 0
        for line=0,nlines-1 do begin
            if not ptr_valid(xval[line]) then continue
            if nlabels eq 0 then $
                plot,*xval[line],*yval[line], xrange=xrange, /xstyle, /ystyle, symsize=0.3, charsize=0.8, /nodata, xtickformat='LABEL_DATE' $
            else $
                plot,*xval[line],*yval[line], xrange=xrange, /xstyle, /ystyle, symsize=0.3, title=labels[line]+lbl, charsize=0.8, $
                    /nodata, xtickformat='LABEL_DATE'
            oplot,*xval[line], *yval[line], color=colors[line], psym=-1
        endfor
        device,/close_file
        set_plot, thisDevice
        !p.background=255
        !p.color=0
        !p.font = -1
        !p.multi=0

    endif else begin
        ; use the plot routine to save a postscript file since the "save" method for IDLgr is soooo sloooow
        thisDevice = !d.name
        set_plot,'Z'
        if n_elements(resolution) eq 2 then this_res=resolution else this_res=[1066,1600]
        device, set_resolution=this_res
        !p.background=255
        !p.color=0
        tvlct,r,g,b,/get
        !p.font = 0
        !P.multi=[0,1,nlines] 
        for line=0,nlines-1 do begin
            if not ptr_valid(xval[line]) then continue
            if nlabels eq 0 then $
                plot,*xval[line],*yval[line], xrange=xrange, /xstyle, /ystyle, symsize=0.3, charsize=1.4, /nodata, xtickformat='LABEL_DATE' $
            else $
                plot,*xval[line],*yval[line], xrange=xrange,/xstyle,  /ystyle, symsize=0.3, title=labels[line]+lbl, charsize=1.4, $
                    /nodata, xtickformat='LABEL_DATE'
            oplot,*xval[line], *yval[line], color=colors[line], psym=-1
        endfor
        if strpos(psfile,'.png') gt 0 then begin
            img=tvrd()
            write_png,psfile,img,r,g,b, xresolution=300, yresolution=300
        endif else begin
            img=tvrd()
            write_jpeg,psfile,img, quality=90
        endelse
        device,/close
        set_plot, thisDevice
        image=0b
        !p.multi=0
    endelse

END
;-----------------------------------------------------------------------------
PRO PLOT_TELEM,  startTime, stopTime, simA=simA, simB=simB,  $
         gps=gps, missionDays=missionDays, julianDays=julianDays, vms=vms, $
         outdata=outdata, items=items, tlmid=tlmid, psfile=psfile, $
         clean=clean, silent=silent, extelem=extelem, resolution=resolution


    if not keyword_set(silent) then $
        print, 'PLOT_TELEM ('+systime()+'):  Processing time range: ',startTime,' to ',stopTime

    ; adjust the startTime and stopTime to get full orbits ENDING within time range
    if keyword_set(missionDays) then begin
        ; user specified time in mission (sorce) days
        t0 = sd2gps(startTime)*1.d6
        t1 = sd2gps(stopTime)*1.d6
    endif else if keyword_set(julianDays) then begin
        ; user specified time in julian days
        t0 = jd2gps(startTime)*1.d6
        t1 = jd2gps(stopTime)*1.d6
    endif else if keyword_set(vms) then begin
        ; user specified time in VMS string format
        t0 = vms2gps(startTime)*1.d6
        t1 = vms2gps(stopTime)*1.d6
    endif else begin
        ; user specified timetags in gps microseconds
        t0 = startTime
        t1 = stopTime
    endelse

    ; if same day is provided go back one day
    if t0 eq t1 then t0-=86400d6

    if size(outdata,/tname) ne 'STRUCT' then begin
        if n_elements(items) eq 0 and n_elements(tlmid) eq 0 then begin
            items = ['uv_temp', 'vis_1_temp','ir_temp',   'esr_temp','case_temp','bench_1_temp',$
                 'prism_drive_temp', 'position','uv_array',  'vis1_array','ir_array','esr_array']
        endif 
        ; get the data
        outdata=show_sim_timeline(t0, t1, /gps,  sima=sima, simb=simb, items=items, tlmid=tlmid, externalElement=extelem, /noplot)
    endif
    items = tag_names(outdata)
    xval=ptrarr(n_elements(items))
    yval=ptrarr(n_elements(items))

    ; remove the outliers before plotting the data
    for i=0,n_elements(items)-1 do begin
        if (*outdata.(i)).n_science gt 0 then begin
            pos = where(strpos(tag_names((*outdata.(i)).science),'EU') eq 0)
            if pos ge 0 then ytemp= (*outdata.(i)).science.eu else ytemp=(*outdata.(i)).science.dn
            if keyword_set(clean) then begin
                resistant_mean, ytemp, 5.0, mean, good=keep
                xtemp = gps2jd((*outdata.(i)).science.timetag/1d6)
                xval[i]=ptr_new(xtemp[keep],/no_copy)
                yval[i]=ptr_new(ytemp[keep],/no_copy)
            endif else begin
                xtemp = gps2jd((*outdata.(i)).science.timetag/1d6)
                xval[i]=ptr_new(xtemp,/no_copy)
                yval[i]=ptr_new(ytemp,/no_copy)
            endelse
        endif else if (*outdata.(i)).n_housekeeping gt 0 then begin
            pos = where(strpos(tag_names((*outdata.(i)).housekeeping),'EU') eq 0)
            if pos ge 0 then ytemp= (*outdata.(i)).housekeeping.eu else ytemp=(*outdata.(i)).housekeeping.dn
            if keyword_set(clean) then begin
                resistant_mean, ytemp, 5.0, mean, good=keep
                xtemp=gps2jd((*outdata.(i)).housekeeping.timetag/1d6)
                xval[i]=ptr_new(xtemp[keep],/no_copy)
                yval[i]=ptr_new(ytemp[keep],/no_copy)
            endif else begin
                xtemp=gps2jd((*outdata.(i)).housekeeping.timetag/1d6)
                xval[i]=ptr_new(xtemp,/no_copy)
                yval[i]=ptr_new(ytemp,/no_copy)
            endelse
        endif else begin
            xval[i]=ptr_new(gps2jd([t0,t1]/1d6))
            yval[i]=ptr_new([0.0,0.0])
        endelse
        ; scale two XPS telem to seconds and kilo counts (as requested by Tom Woods)
        if strpos(items[i],'XPS$INTG_TM') ge 0 or strpos(items[i],'XPS$DIODE_1_DATA') ge 0 then *yval[i]/=1000d
        ptr_free,outdata.(i)
    endfor

    gen_plot_telem,xval,yval, labels=items, psfile=psfile, xrange=gps2jd([t0,t1]/1d6), resolution=resolution
    
    for i=0,n_elements(items)-1 do ptr_free, xval[i], yval[i]

END
