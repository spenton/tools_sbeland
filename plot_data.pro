; $Id: plot_data.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;+
; NAME:
;   PLOT_DATA
; PURPOSE:
;   Plot the data in a scrollable window
;
; CALLING SEQUENCE:
;   PLOT_DATA, cmd, data1 [,data2]
;
; OPTIONAL INPUTS:
;
; OPTIONAL INPUT KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; USAGE EXAMPLE:
;  IDL> x=findgen(100)
;  IDL> y=sin(x)
;  IDL> plot_data,'plot',x,y,xtitle='X',ytitle='Y',title='A Test', psym=-4
;
; NOTES:
;
; PROCEDURE CALLS:
;
; REVISION HISTORY:
;   Written S. Beland October 2001
;   Edited N.J. Cunningham 2007, to restore original states of color
;          table, plot settings (!p.*), and device,/decomposed setting
;          on Quit.
;   Further modifications to protect other device/color settings, and
;          to prevent plot_data from hijacking the current window for command
;          line plots.  NJC 12/2007
;   Allow keyword inheritance (using CALL_PROCEDURE instead of
;          EXECUTE); protect zoom state from interference by managing
;          !x and !y; use EMPTY calls to ensure that crosshairs/zoom
;          boxes show up; remove direct printing menu option; add
;          cleanup procedure; add menu item to undo last zoom.  NJC 2/2008
;-
;===========================================================
;
; Main Event handler
pro plot_data_event,event

  Catch, theError
  IF theError NE 0 THEN BEGIN
     Catch, /Cancel
     if theError eq -529 then begin
         wset,-1
         q.userwindow=!d.window
         widget_control, event.top, set_uvalue=q
     endif else begin
         widget_control, event.top, set_uvalue=q, /NO_COPY
         ;message, /reissue_last
         RETURN
     endelse
  ENDIF

  widget_control, event.top, get_uvalue=q, /NO_COPY

  tvlct,rsave,gsave,bsave,/get
  q.rsave = rsave
  q.gsave = gsave
  q.bsave = bsave
  q.puser = !p
  q.xuser = !x
  q.yuser = !y
  device,get_decomposed = decomp
  q.decomp = decomp
  userwindow = !d.window
  wset,q.plot_id
  tvlct,q.red,q.green,q.blue
  !p = q.psave
  !x = q.xsave
  !y = q.ysave
  device,decomposed = 0
  wset,q.plot_id

  if TAG_NAMES(event, /STRUCTURE_NAME) eq 'WIDGET_BASE' AND $
     event.id EQ q.main then begin ; The window has been sized

     ;; Get the new size of the window
     WIDGET_CONTROL, q.main, TLB_GET_SIZE=windowSize

     ;; Determine the change in the window size
     deltaX = windowSize[0] - q.windowSize[0]
     deltaY = windowSize[1] - q.windowSize[1]

     IF not (deltaX EQ 0.0 AND deltaY EQ 0.0) THEN begin

        widget_control, q.main, map=0

        ;; Get the pixel size of the draw widget
        plotGeometry = WIDGET_INFO(q.wplot, /GEOMETRY)

        ;; Determine the new size based on the amount the window grew
        newXSize = plotGeometry.scr_xsize + deltaX
        newYSize = plotGeometry.scr_ysize + deltaY

        ;; delete scroll bars (or window doesn't size properly)
        WIDGET_CONTROL, q.vscroll, get_value=vpos
        WIDGET_CONTROL, q.hscroll, get_value=hpos
        WIDGET_CONTROL, q.vscroll, /DESTROY
        WIDGET_CONTROL, q.hscroll, /DESTROY

        ;; delete the pixmap and create a new one
        wdelete, q.pixmapid
        window, xs=newXSize, ys=newYSize, /pixmap, /free
        q.pixmapid = !d.window
        q.wXsize = newXsize
        q.wYsize = newYsize

        widget_control, q.wplot, xsize=newXsize, ysize=newYsize

        ;; create new draw and sliders
        if q.xmin EQ q.xmin_abs AND q.xmax EQ q.xmax_abs THEN $
           sensitive=0 $
        else $
           sensitive=1
        q.hscroll = cw_fslider(q.base2, /DRAG, minimum=q.hminmax[0], maximum=q.hminmax[1], $
                       uvalue='HSCROLL', value=hpos, /SUPPRESS_VALUE, XSIZE=q.wXsize)

        widget_control, q.hscroll, sensitive=sensitive

        q.vscroll = cw_fslider(q.base2, /DRAG, /VERTICAL, $
                               minimum=q.vminmax[0], maximum=q.vminmax[1], $
                               value=vpos, uvalue='VSCROLL', $
                               /SUPPRESS_VALUE, YSIZE=q.wYsize)
        widget_control, q.main, map=1

        ;; redraw the plot
        plot_data_genplot,q

        WIDGET_CONTROL, q.main, TLB_GET_SIZE=windowSize
        q.windowSize = windowSize
     endif
  ENDIF else begin

     widget_control, event.id, get_uvalue=uvalue

     case uvalue of
        'top_menu': BEGIN       ; an option from the menubar was selected

           case event.value of
              'File.Quit': BEGIN
                 WIDGET_CONTROL, SET_UVALUE=q, event.top, /NO_COPY
                 WIDGET_CONTROL, event.top, /DESTROY
                 RETURN
              END

              'File.Create PS': BEGIN
                 keywords=q.pslocal
                 WIDGET_CONTROL, event.top, SET_UVALUE=q, /NO_COPY
                 fdecomp, keywords.filename,disk,dir,fname
                 keywords.filename=disk+dir+fname+'.ps'
                 keywords = CMPS_FORM(defaults=keywords, cancel=canceled, create=create)
                 IF canceled THEN return
                 WIDGET_CONTROL, event.top, GET_UVALUE=q, /NO_COPY
                 q.pslocal=keywords
                 if create then begin
                    thisDevice = !D.NAME
                    TVLCT, r, g, b, /GET
                    SET_PLOT, 'PS'
                    DEVICE, _EXTRA=keywords
                                ;tvlct,q.red,q.green,q.blue
                    !p.background=255
                    !p.color=0
                    plot_data_genplot,q
                    DEVICE, /CLOSE_FILE
                    SET_PLOT, thisDevice
                    wset,q.plot_id
                                ;tvlct,q.red,q.green,q.blue
                    !p.background=255
                    !p.color=0
                 endif
              END

;;               'File.Print': BEGIN
;;                  keywords = q.pslocal
;;                  thisDevice = !D.NAME
;;                  fdecomp, keywords.filename,disk,dir,fname
;;                  keywords.filename=disk+dir+fname+'.ps'
;;                  TVLCT, r, g, b, /GET
;;                  SET_PLOT, 'PS'
;;                  DEVICE, _EXTRA=keywords
;;                  tvlct,q.red,q.green,q.blue
;;                  !p.background=0
;;                  !p.color=1
;;                  res=EXECUTE(q.plot_proc)
;;                  DEVICE, /CLOSE_FILE
;;                  SET_PLOT, thisDevice
;;                  wset,q.plot_id
;;                  tvlct,q.red,q.green,q.blue
;;                  !p.background=0
;;                  !p.color=1
;;                  status=send2printer(keywords.filename)
;;               END

              'File.Create JPG': BEGIN
                 fdecomp, q.pslocal.filename,disk,dir,fname
                 filename=disk+dir+fname+'.jpg'
                 wset,q.pixmapid
                 img=tvrd(true=3)
                 wset,q.plot_id
                 file=dialog_pickfile(title='Select Output Filename (JPG)', $
                                      filter='*.jpg',file=filename)
                 if file ne '' then begin
                    write_jpeg,file,img,true=3
                    img=0b
                 endif
              END

              'File.Create PNG': BEGIN
                 fdecomp, q.pslocal.filename,disk,dir,fname
                 filename=disk+dir+fname+'.png'
                 wset,q.pixmapid
                 img=tvrd(true=1)
                 wset,q.plot_id
                 file=dialog_pickfile(title='Select Output Filename (PNG)', $
                                      filter='*.png',file=filename)
                 if file ne '' then begin
                    ;; could have chosen indexed color, but have to
                    ;; generate the indices instead of using our
                    ;; current color table!
                                ;indexed_img = color_quan(img,1,red,green,blue)
                                ;write_png,file,indexed_img,red,green,blue
                    write_png,file,img
                    img=0b
                 endif
              END

              'Display.Hide X,Y': begin
                 widget_control,q.showxy_id,get_value=str1
                 if strpos(str1,'Show') ge 0 then begin
                    ;; show the X,Y cursor coordinates
                    widget_control, q.showxy_id,set_value='Hide X,Y'
                    q.showxy=1
                 endif else begin
                    ;; Hide the X,Y cursor coordinates
                    widget_control, q.showxy_lbl,set_value=string(replicate(32b,40))
                    widget_control, q.showxy_id,set_value='Show X,Y'
                    q.showxy=0
                 endelse
              end

              'Display.Zoom': begin
                 widget_control,q.message,set_value= $
                                'Place Cursor on first corner'+ $
                                ' and push left mouse button'
                 q.state = 'ZOOM1'
              end

              'Display.Undo Last Zoom': begin
                 q.xmin = q.xmin_undo
                 q.xmax = q.xmax_undo
                 q.ymin = q.ymin_undo
                 q.ymax = q.ymax_undo
                 wset,q.pixmapid
                 plot_data_genplot,q
                 wset,q.plot_id
                 device,copy=[0,0,q.wXsize,q.wYsize,0,0,q.pixmapid]
                 q.state = 'X/Y'
                 widget_control, q.undozoom_id,sensitive=0
              end

              'Display.UnZoom X': begin
                 q.xmin = q.xmin_abs
                 q.xmax = q.xmax_abs
                 widget_control, q.hscroll, set_value=q.xmin_abs
                 widget_control, q.hscroll, sensitive=0
                 wset,q.pixmapid
                 plot_data_genplot,q
                 wset,q.plot_id
                 device,copy=[0,0,q.wXsize,q.wYsize,0,0,q.pixmapid]
              end

              'Display.UnZoom Y': begin
                 q.ymin = q.ymin_abs
                 q.ymax = q.ymax_abs
                 widget_control, q.vscroll, set_slider_max=q.ymax_abs
                 widget_control, q.vscroll, set_slider_min=q.ymin_abs
                 widget_control, q.vscroll, set_value=q.ymin_abs
                 widget_control, q.vscroll, sensitive=0
                 wset,q.pixmapid
                 plot_data_genplot,q
                 wset,q.plot_id
                 device,copy=[0,0,q.wXsize,q.wYsize,0,0,q.pixmapid]
              end

              'Display.UnZoom All': begin
                 q.xmin = q.xmin_abs
                 q.xmax = q.xmax_abs
                 q.ymin = q.ymin_abs
                 q.ymax = q.ymax_abs
                 widget_control, q.vscroll, set_slider_max=q.ymax_abs
                 widget_control, q.vscroll, set_slider_min=q.ymin_abs
                 widget_control, q.vscroll, set_value=q.ymin_abs
                 widget_control, q.hscroll, set_value=q.xmin_abs
                 widget_control, q.hscroll, sensitive=0
                 widget_control, q.vscroll, sensitive=0
                 wset,q.pixmapid
                 plot_data_genplot,q
                 wset,q.plot_id
                 device,copy=[0,0,q.wXsize,q.wYsize,0,0,q.pixmapid]
              end

              else:
           endcase
        END


        'VSCROLL': begin
           ymax_current = q.ymax
           ymin_current = q.ymin
           delta = ymax_current - ymin_current
           widget_control, q.vscroll, get_value=vscroll
           ymin_new = (vscroll-delta/2.0d) > q.ymin_abs
           ymax_new = (ymin_new + delta) < q.ymax_abs
           if ymax_new eq q.ymax_abs then ymin_new = (ymax_new - delta) > q.ymin_abs
           q.ymax = ymax_new
           q.ymin = ymin_new
           wset,q.pixmapid
           plot_data_genplot,q
           wset,q.plot_id
           device,copy=[0,0,q.wXsize,q.wYsize,0,0,q.pixmapid]
        end

        'HSCROLL': begin
           xmax_current = q.xmax
           xmin_current = q.xmin
           delta = xmax_current - xmin_current
           widget_control, q.hscroll, get_value=hscroll
           xmin_new = (hscroll-delta/2.0) > q.xmin_abs
           xmax_new = (xmin_new + delta) < q.xmax_abs
           if xmax_new eq q.xmax_abs then xmin_new = (xmax_new - delta) > q.xmin_abs
           q.xmax = xmax_new
           q.xmin = xmin_new
           wset,q.pixmapid
           plot_data_genplot,q
           wset,q.plot_id
           device,copy=[0,0,q.wXsize,q.wYsize,0,0,q.pixmapid]
        end

        'PLOT1': begin
           xd = event.x         ;device coordinates
           yd = event.y
           if q.showxy then begin
              res=convert_coord(xd,yd,/dev,/to_data,/double)
              if q.xdate then begin
                 ;; convert x value to date/time
                 caldat, res[0], month, day, year, hour, minute, sec
                 year=fix(year) - fix(year/100)*100
                 sdate=string(format="(I02,'/',I02,'/',I02,' ',I02,':',I02,':',I02)",$
                              month,day,year,hour,minute,sec)
                 widget_control, q.showxy_lbl, set_value=sdate+'   '+string(res[1])
              endif else begin
                 widget_control, q.showxy_lbl, $
                                 set_value=string(res[0])+'   '+string(res[1])
              endelse
           endif
           wset, q.plot_id
           device, copy=[0,0,q.wXsize,q.wYsize,0,0,q.pixmapid]
           if (q.state eq 'X/Y') then begin
              plots, [xd,xd], [0,q.wYsize], color=0, /dev
              plots, [0,q.wXsize], [yd,yd], color=0, /dev
              empty            ;ensure that buffered data is written to graphics device
              ;; (otherwise sometimes these lines don't appear!)
              if event.press eq 2 then begin
                 q.x1 = xd
                 q.y1 = yd
                 q.state = 'ZOOM2'
              endif
           endif

           if q.state eq 'ZOOM1' and (event.press eq 1) then begin
              plots,[xd,xd],[0,q.wYsize],color=0,/dev
              plots,[0,q.wXsize],[yd,yd],color=0,/dev
              empty            ;ensure that buffered data is written to graphics device
              q.x1 = xd
              q.y1 = yd
              widget_control,q.message,set_value= $
                             'Position at second corner and'+' push left mouse button'
              q.state = 'ZOOM2'
           endif

           if (q.state eq 'ZOOM2') then begin
              x = [q.x1,xd]
              y = [q.y1,yd]
              plots,[x[0],x[1],x[1],x[0],x[0]],[y[0],y[0],y[1],y[1],y[0]],/dev,color=0
              empty            ;ensure that buffered data is written to graphics device
              if (event.release eq 2) or (event.press eq 1) then begin
                 if abs(x[1]-x[0]) lt 2 then begin
                    widget_control,q.message,set_value=' '
                    wset,q.pixmapid
                    plot_data_genplot,q
                    wset,q.plot_id
                    device,copy=[0,0,q.wXsize,q.wYsize,0,0,q.pixmapid]
                 endif else begin
                    v = convert_coord(x,y,/dev,/to_data,/double)
                    x = v[0,*]
                    y = v[1,*]

                    q.xmin_undo = q.xmin
                    q.xmax_undo = q.xmax
                    q.ymin_undo = q.ymin
                    q.ymax_undo = q.ymax
                    widget_control, q.undozoom_id,/sensitive
                    q.xmin = min(x)
                    q.xmax = max(x)
                    q.ymin = min(y)
                    q.ymax = max(y)

                    if (q.ymin NE q.ymin_abs) OR (q.ymax NE q.ymax_abs) THEN BEGIN
                       widget_control, q.vscroll, set_value=(q.ymax+q.ymin)/2.0
                       widget_control, q.vscroll, sensitive=1
                    endif

                    if (q.xmin NE q.xmin_abs) OR (q.xmax NE q.xmax_abs) THEN BEGIN
                       widget_control, q.hscroll, set_value=(q.xmax+q.xmin)/2.0
                       widget_control, q.hscroll, sensitive=1
                    endif

                    widget_control,q.message,set_value=' '
                    wset,q.pixmapid
                    plot_data_genplot,q
                    wset,q.plot_id
                    device,copy=[0,0,q.wXsize,q.wYsize,0,0,q.pixmapid]
                    q.state = 'X/Y'
                 endelse
              endif
           endif

        end

        else:
     endcase

  endelse

  q.xsave = !x
  q.ysave = !y
  tvlct,q.rsave,q.gsave,q.bsave
  wset,q.userwindow
  !p = q.puser
  !x = q.xuser
  !y = q.yuser
  device,decomposed = q.decomp
  ;; Reset the windows user value to the updated state structure
  WIDGET_CONTROL, event.top, SET_UVALUE=q, /NO_COPY

end

;===========================================================
PRO plot_data_cleanup, tlb
  widget_control, tlb, get_uvalue = q, /no_copy

  ptr_free, q.plotkw
  wdelete, q.pixmapid

  tvlct,q.rsave,q.gsave,q.bsave
  ;wset,q.userwindow
  !p = q.puser
  !x = q.xuser
  !y = q.yuser
  device,decomposed = q.decomp
END


;===========================================================
pro plot_data_genplot,q,norange=norange
;; procedure to call the plot routine with appropriate keywords.
;; NORANGE keyword is called initially, so that any unspecified range
;; remains unset, to be determined automatically by the plotting
;; routine.

;; If NORANGE is not set, range needs to be set by program, not
;; according to users's initial keywords.

  if n_elements(*q.plotkw) gt 0 then extra = *q.plotkw
  if not keyword_set(norange) then begin
     xr = [q.xmin,q.xmax]
     yr = [q.ymin,q.ymax]
     ;; ensure that XRANGE and YRANGE tag names exist in the EXTRA
     ;; structure, with appropriate values.
     if n_elements(extra) eq 0 then extra = {xrange: xr, yrange:yr} else begin
        if total(tag_names(extra) eq 'XRANGE') gt 0 then extra.xrange = xr else $
           extra = create_struct(temporary(extra),'xrange',xr)
        if total(tag_names(extra) eq 'YRANGE') gt 0 then extra.yrange = yr else $
           extra = create_struct(temporary(extra),'yrange',yr)
     endelse
  endif

  case q.npars of
     1: call_procedure,q.plot_proc,q.d1,_strict_extra=extra
     2: call_procedure,q.plot_proc,q.d1,q.d2,_strict_extra=extra
     3: call_procedure,q.plot_proc,q.d1,q.d2,q.d3,_strict_extra=extra
     4: call_procedure,q.plot_proc,q.d1,q.d2,q.d3,q.d4,_strict_extra=extra
     5: call_procedure,q.plot_proc,q.d1,q.d2,q.d3,q.d4,q.d5,_strict_extra=extra
     6: call_procedure,q.plot_proc,q.d1,q.d2,q.d3,q.d4,q.d5,q.d6, $
                       _strict_extra=extra
     7: call_procedure,q.plot_proc,q.d1,q.d2,q.d3,q.d4,q.d5,q.d6,q.d7, $
                       _strict_extra=extra
     8: call_procedure,q.plot_proc,q.d1,q.d2,q.d3,q.d4,q.d5,q.d6,q.d7,q.d8, $
                       _strict_extra=extra
     9: call_procedure,q.plot_proc,q.d1,q.d2,q.d3,q.d4,q.d5,q.d6,q.d7,q.d8, $
                       q.d9,_strict_extra=extra
     10: call_procedure,q.plot_proc,q.d1,q.d2,q.d3,q.d4,q.d5,q.d6,q.d7,q.d8, $
                        q.d9,q.d10,_strict_extra=extra
     11: call_procedure,q.plot_proc,q.d1,q.d2,q.d3,q.d4,q.d5,q.d6,q.d7,q.d8, $
                        q.d9,q.d10,q.d11,_strict_extra=extra
     12: call_procedure,q.plot_proc,q.d1,q.d2,q.d3,q.d4,q.d5,q.d6,q.d7,q.d8, $
                        q.d9,q.d10,q.d11,q.d12,_strict_extra=extra
     13: call_procedure,q.plot_proc,q.d1,q.d2,q.d3,q.d4,q.d5,q.d6,q.d7,q.d8, $
                        q.d9,q.d10,q.d11,q.d12,q.d13,_strict_extra=extra
     14: call_procedure,q.plot_proc,q.d1,q.d2,q.d3,q.d4,q.d5,q.d6,q.d7,q.d8, $
                        q.d9,q.d10,q.d11,q.d12,q.d13,q.d14,_strict_extra=extra
     15: call_procedure,q.plot_proc,q.d1,q.d2,q.d3,q.d4,q.d5,q.d6,q.d7,q.d8, $
                        q.d9,q.d10,q.d11,q.d12,q.d13,q.d14,q.d15,_strict_extra=extra
     16: call_procedure,q.plot_proc,q.d1,q.d2,q.d3,q.d4,q.d5,q.d6,q.d7,q.d8, $
                        q.d9,q.d10,q.d11,q.d12,q.d13,q.d14,q.d15,q.d16, $
                        _strict_extra=extra
     17: call_procedure,q.plot_proc,q.d1,q.d2,q.d3,q.d4,q.d5,q.d6,q.d7,q.d8, $
                        q.d9,q.d10,q.d11,q.d12,q.d13,q.d14,q.d15,q.d16, $
                        q.d17,_strict_extra=extra
     18: call_procedure,q.plot_proc,q.d1,q.d2,q.d3,q.d4,q.d5,q.d6,q.d7,q.d8, $
                        q.d9,q.d10,q.d11,q.d12,q.d13,q.d14,q.d15,q.d16, $
                        q.d17,q.d18,_strict_extra=extra
     19: call_procedure,q.plot_proc,q.d1,q.d2,q.d3,q.d4,q.d5,q.d6,q.d7,q.d8, $
                        q.d9,q.d10,q.d11,q.d12,q.d13,q.d14,q.d15,q.d16, $
                        q.d17,q.d18,q.d19,_strict_extra=extra
     20: call_procedure,q.plot_proc,q.d1,q.d2,q.d3,q.d4,q.d5,q.d6,q.d7,q.d8, $
                        q.d9,q.d10,q.d11,q.d12,q.d13,q.d14,q.d15,q.d16, $
                        q.d17,q.d18,q.d19,q.d20,_strict_extra=extra
     21: call_procedure,q.plot_proc,q.d1,q.d2,q.d3,q.d4,q.d5,q.d6,q.d7,q.d8, $
                        q.d9,q.d10,q.d11,q.d12,q.d13,q.d14,q.d15,q.d16, $
                        q.d17,q.d18,q.d19,q.d20,q.d21,_strict_extra=extra
     22: call_procedure,q.plot_proc,q.d1,q.d2,q.d3,q.d4,q.d5,q.d6,q.d7,q.d8, $
                        q.d9,q.d10,q.d11,q.d12,q.d13,q.d14,q.d15,q.d16, $
                        q.d17,q.d18,q.d19,q.d20,q.d21,q.d22,_strict_extra=extra
     23: call_procedure,q.plot_proc,q.d1,q.d2,q.d3,q.d4,q.d5,q.d6,q.d7,q.d8, $
                        q.d9,q.d10,q.d11,q.d12,q.d13,q.d14,q.d15,q.d16, $
                        q.d17,q.d18,q.d19,q.d20,q.d21,q.d22,q.d23,_strict_extra=extra
     24: call_procedure,q.plot_proc,q.d1,q.d2,q.d3,q.d4,q.d5,q.d6,q.d7,q.d8, $
                        q.d9,q.d10,q.d11,q.d12,q.d13,q.d14,q.d15,q.d16, $
                        q.d17,q.d18,q.d19,q.d20,q.d21,q.d22,q.d23,q.d24, $
                        _strict_extra=extra
     25: call_procedure,q.plot_proc,q.d1,q.d2,q.d3,q.d4,q.d5,q.d6,q.d7,q.d8, $
                        q.d9,q.d10,q.d11,q.d12,q.d13,q.d14,q.d15,q.d16, $
                        q.d17,q.d18,q.d19,q.d20,q.d21,q.d22,q.d23,q.d24, $
         q.d25, _strict_extra=extra
     26: call_procedure,q.plot_proc,q.d1,q.d2,q.d3,q.d4,q.d5,q.d6,q.d7,q.d8, $
                        q.d9,q.d10,q.d11,q.d12,q.d13,q.d14,q.d15,q.d16, $
                        q.d17,q.d18,q.d19,q.d20,q.d21,q.d22,q.d23,q.d24, $
         q.d25,q.d26, _strict_extra=extra
     27: call_procedure,q.plot_proc,q.d1,q.d2,q.d3,q.d4,q.d5,q.d6,q.d7,q.d8, $
                        q.d9,q.d10,q.d11,q.d12,q.d13,q.d14,q.d15,q.d16, $
                        q.d17,q.d18,q.d19,q.d20,q.d21,q.d22,q.d23,q.d24, $
         q.d25,q.d26,q.d27, _strict_extra=extra
     28: call_procedure,q.plot_proc,q.d1,q.d2,q.d3,q.d4,q.d5,q.d6,q.d7,q.d8, $
                        q.d9,q.d10,q.d11,q.d12,q.d13,q.d14,q.d15,q.d16, $
                        q.d17,q.d18,q.d19,q.d20,q.d21,q.d22,q.d23,q.d24, $
         q.d25,q.d26,q.d27,q.d28, _strict_extra=extra
     29: call_procedure,q.plot_proc,q.d1,q.d2,q.d3,q.d4,q.d5,q.d6,q.d7,q.d8, $
                        q.d9,q.d10,q.d11,q.d12,q.d13,q.d14,q.d15,q.d16, $
                        q.d17,q.d18,q.d19,q.d20,q.d21,q.d22,q.d23,q.d24, $
         q.d25,q.d26,q.d27,q.d28,q.d29, _strict_extra=extra
     30: call_procedure,q.plot_proc,q.d1,q.d2,q.d3,q.d4,q.d5,q.d6,q.d7,q.d8, $
                        q.d9,q.d10,q.d11,q.d12,q.d13,q.d14,q.d15,q.d16, $
                        q.d17,q.d18,q.d19,q.d20,q.d21,q.d22,q.d23,q.d24, $
         q.d25,q.d26,q.d27,q.d28,q.d29,q.d30, _strict_extra=extra
     31: call_procedure,q.plot_proc,q.d1,q.d2,q.d3,q.d4,q.d5,q.d6,q.d7,q.d8, $
                        q.d9,q.d10,q.d11,q.d12,q.d13,q.d14,q.d15,q.d16, $
                        q.d17,q.d18,q.d19,q.d20,q.d21,q.d22,q.d23,q.d24, $
         q.d25,q.d26,q.d27,q.d28,q.d29,q.d30,q.d31, _strict_extra=extra
     32: call_procedure,q.plot_proc,q.d1,q.d2,q.d3,q.d4,q.d5,q.d6,q.d7,q.d8, $
                        q.d9,q.d10,q.d11,q.d12,q.d13,q.d14,q.d15,q.d16, $
                        q.d17,q.d18,q.d19,q.d20,q.d21,q.d22,q.d23,q.d24, $
         q.d25,q.d26,q.d27,q.d28,q.d29,q.d30,q.d31,q.d32, _strict_extra=extra
     33: call_procedure,q.plot_proc,q.d1,q.d2,q.d3,q.d4,q.d5,q.d6,q.d7,q.d8, $
                        q.d9,q.d10,q.d11,q.d12,q.d13,q.d14,q.d15,q.d16, $
                        q.d17,q.d18,q.d19,q.d20,q.d21,q.d22,q.d23,q.d24, $
          q.d25,q.d26,q.d27,q.d28,q.d29,q.d30,q.d31,q.d32,q.d33, _strict_extra=extra
     34: call_procedure,q.plot_proc,q.d1,q.d2,q.d3,q.d4,q.d5,q.d6,q.d7,q.d8, $
                        q.d9,q.d10,q.d11,q.d12,q.d13,q.d14,q.d15,q.d16, $
                        q.d17,q.d18,q.d19,q.d20,q.d21,q.d22,q.d23,q.d24, $
         q.d25,q.d26,q.d27,q.d28,q.d29,q.d30,q.d31,q.d32,q.d33,q.d34, _strict_extra=extra
     35: call_procedure,q.plot_proc,q.d1,q.d2,q.d3,q.d4,q.d5,q.d6,q.d7,q.d8, $
                        q.d9,q.d10,q.d11,q.d12,q.d13,q.d14,q.d15,q.d16, $
                        q.d17,q.d18,q.d19,q.d20,q.d21,q.d22,q.d23,q.d24, $
         q.d25,q.d26,q.d27,q.d28,q.d29,q.d30,q.d31,q.d32,q.d33,q.d34, $
         q.d35, _strict_extra=extra
     36: call_procedure,q.plot_proc,q.d1,q.d2,q.d3,q.d4,q.d5,q.d6,q.d7,q.d8, $
                        q.d9,q.d10,q.d11,q.d12,q.d13,q.d14,q.d15,q.d16, $
                        q.d17,q.d18,q.d19,q.d20,q.d21,q.d22,q.d23,q.d24, $
         q.d25,q.d26,q.d27,q.d28,q.d29,q.d30,q.d31,q.d32,q.d33,q.d34, $
         q.d35,q.d36, _strict_extra=extra
     37: call_procedure,q.plot_proc,q.d1,q.d2,q.d3,q.d4,q.d5,q.d6,q.d7,q.d8, $
                        q.d9,q.d10,q.d11,q.d12,q.d13,q.d14,q.d15,q.d16, $
                        q.d17,q.d18,q.d19,q.d20,q.d21,q.d22,q.d23,q.d24, $
         q.d25,q.d26,q.d27,q.d28,q.d29,q.d30,q.d31,q.d32,q.d33,q.d34, $
         q.d35,q.d36,q.d37, _strict_extra=extra
     38: call_procedure,q.plot_proc,q.d1,q.d2,q.d3,q.d4,q.d5,q.d6,q.d7,q.d8, $
                        q.d9,q.d10,q.d11,q.d12,q.d13,q.d14,q.d15,q.d16, $
                        q.d17,q.d18,q.d19,q.d20,q.d21,q.d22,q.d23,q.d24, $
         q.d25,q.d26,q.d27,q.d28,q.d29,q.d30,q.d31,q.d32,q.d33,q.d34, $
         q.d35,q.d36,q.d37,q.d38, _strict_extra=extra
     39: call_procedure,q.plot_proc,q.d1,q.d2,q.d3,q.d4,q.d5,q.d6,q.d7,q.d8, $
                        q.d9,q.d10,q.d11,q.d12,q.d13,q.d14,q.d15,q.d16, $
                        q.d17,q.d18,q.d19,q.d20,q.d21,q.d22,q.d23,q.d24, $
         q.d25,q.d26,q.d27,q.d28,q.d29,q.d30,q.d31,q.d32,q.d33,q.d34, $
         q.d35,q.d36,q.d37,q.d38,q.d39, _strict_extra=extra
     40: call_procedure,q.plot_proc,q.d1,q.d2,q.d3,q.d4,q.d5,q.d6,q.d7,q.d8, $
                        q.d9,q.d10,q.d11,q.d12,q.d13,q.d14,q.d15,q.d16, $
                        q.d17,q.d18,q.d19,q.d20,q.d21,q.d22,q.d23,q.d24, $
         q.d25,q.d26,q.d27,q.d28,q.d29,q.d30,q.d31,q.d32,q.d33,q.d34, $
         q.d35,q.d36,q.d37,q.d38,q.d39,q.d40, _strict_extra=extra
  endcase

end


;===========================================================
;
; Main widget driver
;
pro plot_data, plot_proc, d1, d2, d3, d4, d5, d6, d7, d8, $
               d9, d10, d11, d12, d13, d14, d15, d16, d17, d18, $
               d19, d20, d21, d22, d23, d24, d25, d26, d27, d28, d29, d30, $
               d31, d32, d33, d34, d35, d36, d37, d38, d39, d40, $
               WTITLE=wtitle, WXSIZE=wxsize, WYSIZE=wysize, GROUP=group, $
               OUTNAME=outname, border=border, xdate=xdate,_extra=extra

  ; Error handling.
  On_Error, 2
  Catch, theError
  IF theError NE 0 THEN BEGIN
     Catch, /Cancel
     message, /reissue_last
     RETURN
  ENDIF

;
; initialization
;

  npars = N_PARAMS()
  IF npars LT 2 then begin
     print,'Usage: plot_data, plot_proc, d1 [,d2, d3, d4, d5, d6,d7,d8]'
     return
  ENDIF

  ;; make sure we're using a proper display device!
  ;; (e.g. Can't use plot_data on PS device)
  devname = !d.name
  if devname ne 'X' and devname ne 'WIN' then begin
     message,'Device must be set to X or WIN.  (Currently set to '+devname+')'
     return
  endif

  IF size(plot_proc,/TNAME) NE 'STRING' THEN BEGIN
     print,'Usage: plot_data, plot_proc, d1 [,d2, d3, d4, d5, d6,d7,d8]'
     print,'       plot_proc has to be the name of a ploting procedure,'
     print,'       which must accept the XRANGE and YRANGE keywords'
     return
  ENDIF
  if not keyword_set(wXsize) then wXsize = 850
  if not keyword_set(wYsize) then wYsize = 550
  if not keyword_set(wtitle) then wtitle=''
  if not keyword_set(xdate) then xdate=0
  border = keyword_set(border)

  device,get_decomposed = decomp
  tvlct,red,green,blue,/get
  rsave = red
  gsave = green
  bsave = blue
  puser = !p
  userwindow = !d.window

  device,decomposed=0
;                0   1   2   3   4   5   6   7   8   9   10  11
  ; red(0) =   [255,  0,255,  0,  0,255,255,  0,240, 75,110,150]
  ; green(0) = [255,  0,  0,255,  0,255,  0,255,140,180,110,120]
  ; blue(0) =  [255,  0,  0,  0,255,  0,255,255,100, 40,185,150]
  loadct,0,/silent
  tek_color,0,12
  ;tvlct,red,green,blue
  !p.background=255
  !p.color=0

  ;; If we have an error after this point, we need to reset the color
  ;; plot and device settings to their previous values
  Catch, theError
  IF theError NE 0 THEN BEGIN
     Catch, /Cancel
     tvlct,rsave,gsave,bsave
     !p = puser
     device,decomposed = decomp
     wset,userwindow
     message, /reissue_last
     RETURN
  ENDIF

;  !p.region=[0.08,0.08,0.92,0.92]

  if n_elements(d1)  eq 0 then  d1=0b else npars=1 
  if n_elements(d2)  eq 0 then  d2=0b else npars=2 
  if n_elements(d3)  eq 0 then  d3=0b else npars=3 
  if n_elements(d4)  eq 0 then  d4=0b else npars=4 
  if n_elements(d5)  eq 0 then  d5=0b else npars=5 
  if n_elements(d6)  eq 0 then  d6=0b else npars=6 
  if n_elements(d7)  eq 0 then  d7=0b else npars=7 
  if n_elements(d8)  eq 0 then  d8=0b else npars=8 
  if n_elements(d9)  eq 0 then  d9=0b else npars=9 
  if n_elements(d10) eq 0 then  d10=0b else npars=10
  if n_elements(d11) eq 0 then  d11=0b else npars=11 
  if n_elements(d12) eq 0 then  d12=0b else npars=12 
  if n_elements(d13) eq 0 then  d13=0b else npars=13 
  if n_elements(d14) eq 0 then  d14=0b else npars=14 
  if n_elements(d15) eq 0 then  d15=0b else npars=15 
  if n_elements(d16) eq 0 then  d16=0b else npars=16 
  if n_elements(d17) eq 0 then  d17=0b else npars=17 
  if n_elements(d18) eq 0 then  d18=0b else npars=18 
  if n_elements(d19) eq 0 then  d19=0b else npars=19 
  if n_elements(d20) eq 0 then  d20=0b else npars=20 
  if n_elements(d21) eq 0 then  d21=0b else npars=21 
  if n_elements(d22) eq 0 then  d22=0b else npars=22 
  if n_elements(d23) eq 0 then  d23=0b else npars=23 
  if n_elements(d24) eq 0 then  d24=0b else npars=24 
  if n_elements(d25) eq 0 then  d25=0b else npars=25 
  if n_elements(d26) eq 0 then  d26=0b else npars=26 
  if n_elements(d27) eq 0 then  d27=0b else npars=27 
  if n_elements(d28) eq 0 then  d28=0b else npars=28 
  if n_elements(d29) eq 0 then  d29=0b else npars=29 
  if n_elements(d30) eq 0 then  d30=0b else npars=30 
  if n_elements(d31) eq 0 then  d31=0b else npars=31 
  if n_elements(d32) eq 0 then  d32=0b else npars=32 
  if n_elements(d33) eq 0 then  d33=0b else npars=33 
  if n_elements(d34) eq 0 then  d34=0b else npars=34 
  if n_elements(d35) eq 0 then  d35=0b else npars=35 
  if n_elements(d36) eq 0 then  d36=0b else npars=36 
  if n_elements(d37) eq 0 then  d37=0b else npars=37 
  if n_elements(d38) eq 0 then  d38=0b else npars=38 
  if n_elements(d39) eq 0 then  d39=0b else npars=39 
  if n_elements(d40) eq 0 then  d40=0b else npars=40 

  if not keyword_set(group) then $
     main = widget_base(title=wtitle, /col, mbar=menubar, uvalue='MAIN', /TLB_SIZE_EVENTS) $
  else $
     main = widget_base(title=wtitle, group=group, /col, mbar=menubar, uvalue='MAIN', /TLB_SIZE_EVENTS)

  tmp_struct = {cw_pdmenu_s, flags:0, name:''}
  top_menu_desc = [ $
                  {cw_pdmenu_s, 1, 'File'}, $ ; File menu
                  {cw_pdmenu_s, 0, 'Create PS'}, $
                  {cw_pdmenu_s, 0, 'Create JPG'}, $
                  {cw_pdmenu_s, 0, 'Create PNG'}, $
;;                   {cw_pdmenu_s, 0, 'Print'}, $
                  {cw_pdmenu_s, 2, 'Quit'}, $
                  {cw_pdmenu_s, 1, 'Display'}, $ ; Display menu
                  {cw_pdmenu_s, 0, 'Hide X,Y'}, $
                  {cw_pdmenu_s, 0, 'Zoom'}, $
                  {cw_pdmenu_s, 0, 'Undo Last Zoom'}, $
                  {cw_pdmenu_s, 0, 'UnZoom X'}, $
                  {cw_pdmenu_s, 0, 'UnZoom Y'}, $
                  {cw_pdmenu_s, 2, 'UnZoom All'} $
                  ]
  top_menu = cw_pdmenu(menubar, top_menu_desc, /mbar, /help, $
                       IDS=mbar_ids, /return_full_name, uvalue='top_menu')
  showxy_lbl = widget_label(main,value='       ',/DYNAMIC_RESIZE,/ALIGN_LEFT,Font='8x13')
  base2 = widget_base(main,col=2)
  showxy_id=mbar_ids[6]
  undozoom_id=mbar_ids[8]
  widget_control, undozoom_id, sensitive=0

;
; plot window
;
  hmin=0.0d & hmax=1.0d
  vmin=0.0d & vmax=1.0d
  wplot = widget_draw(base2,uvalue='PLOT1',retain=2, $
                      xsize=wXsize, ysize=wYsize, /button_events,/motion)
  hscroll = cw_fslider(base2, /DRAG, minimum=hmin, maximum=hmax, $
                       uvalue='HSCROLL', value=hmin, /SUPPRESS_VALUE, XSIZE=wXsize)
  vscroll = cw_fslider(base2, /DRAG, minimum=vmin, maximum=vmax, $
                       /VERTICAL, uvalue='VSCROLL', value=vmax, $
                       /SUPPRESS_VALUE, YSIZE=wYsize)
  message = widget_label(main,value='       ',/DYNAMIC_RESIZE, $
                         /ALIGN_LEFT, Font='8x13')
;
; save bases
;
  widget_control, main, /realize
  widget_control, wplot, get_value=plot_id

  WIDGET_CONTROL, main, TLB_GET_SIZE=windowSize

  window, xs=wXsize, ys=wYsize, /pixmap, /free ;create pixmap
  pixmapid = !d.window

  pslocal = CMPS_FORM(/init,/landscape,FILENAME=outname,xoffset=0.7,yoffset=0.8,$
                      xsize=9.70,ysize=7.0,/color,/encapsulated)

  q = {main:main, $
       base2:base2, $
       wplot:wplot, $
       plot_id:plot_id, $
       mbar_ids:mbar_ids, $
       message:message, $
       windowSize:windowSize, $
       wXsize:wXsize, $
       wYsize:wYsize, $
       xmin:hmin, $
       xmax:hmax, $
       ymin:vmin, $
       ymax:vmax, $
       xmin_undo: 0., $
       xmax_undo: 0., $
       ymin_undo: 0., $
       ymax_undo: 0., $
       xmax_abs:hmax, $
       xmin_abs:hmin, $
       ymax_abs:vmax, $
       ymin_abs:vmin, $
       vscroll:vscroll, $
       hscroll:hscroll, $
       rsave:rsave, $
       bsave:bsave, $
       gsave:gsave, $
       puser:puser, $
       xuser:!x, $
       yuser:!y, $
       userwindow: userwindow, $
       state:'X/Y', $
       x1:0L, $
       y1:0L, $
       psave:!p, $
       xsave:!x, $
       ysave:!y, $
       red:red, $
       blue:blue, $
       green:green, $
       decomp:decomp, $
       pixmapid:pixmapid, $
       pslocal:pslocal, $
       plot_proc: plot_proc, $
       d1: d1, $
       d2: d2, $
       d3: d3, $
       d4: d4, $
       d5: d5, $
       d6: d6, $
       d7: d7, $
       d8: d8, $
       d9: d9, $
       d10: d10, $
       d11: d11, $
       d12: d12, $
       d13: d13, $
       d14: d14, $
       d15: d15, $
       d16: d16, $
       d17: d17, $
       d18: d18, $
       d19: d19, $
       d20: d20, $
       d21: d21, $
       d22: d22, $
       d23: d23, $
       d24: d24, $
       d25: d25, $
       d26: d26, $
       d27: d27, $
       d28: d28, $
       d29: d29, $
       d30: d30, $
       d31: d31, $
       d32: d32, $
       d33: d33, $
       d34: d34, $
       d35: d35, $
       d36: d36, $
       d37: d37, $
       d38: d38, $
       d39: d39, $
       d40: d40, $
       npars:npars, $                  ;don't include plot_proc string
       plotkw:ptr_new(/allocate_heap), $ ;keywords to be passed to plot_proc
       showxy_id:showxy_id, $
       showxy_lbl:showxy_lbl, $
       showxy:1, $
       undozoom_id:undozoom_id, $
       xdate:xdate, $
       vminmax:[0.0d,0.0d], $
       hminmax:[0.0d,0.0d] }

  if n_elements(extra) gt 0 then *q.plotkw=extra

  ;; get the range of the plot by ploting to pixmap
  wset,q.pixmapid
  ;; let the plot_proc automatically set ranges, unless user supplies any
  plot_data_genplot,q,/norange
  q.xmin = !x.crange[0]
  q.xmax = !x.crange[1]
  q.ymin = !y.crange[0]
  q.ymax = !y.crange[1]
  if border then begin
     ;;add space around plot
     xdelta = DOUBLE(q.xmax-q.xmin)*0.05d
     q.xmin = DOUBLE(q.xmin) - xdelta
     q.xmax = DOUBLE(q.xmax) + xdelta
     ydelta = DOUBLE(q.ymax-q.ymin)*0.10d
     q.ymin = DOUBLE(q.ymin) - ydelta
     q.ymax = DOUBLE(q.ymax) + ydelta
  endif
  q.xmin_abs = q.xmin
  q.xmax_abs = q.xmax
  q.ymin_abs = q.ymin
  q.ymax_abs = q.ymax
  widget_control, q.hscroll, set_value=[q.xmin_abs,q.xmin_abs,q.xmax_abs]
  widget_control, q.vscroll, set_value=[q.ymin_abs,q.ymin_abs,q.ymax_abs]
  widget_control, q.hscroll, sensitive=0
  q.hminmax=[q.xmin_abs,q.xmax_abs]
  q.vminmax=[q.ymin_abs,q.ymax_abs]
  plot_data_genplot,q
  wset,q.plot_id
  plot_data_genplot,q

  q.xsave=!x
  q.ysave=!y

  WIDGET_CONTROL, main, SET_UVALUE=q
  xmanager,'plot_data', main, /no_block, cleanup = 'plot_data_cleanup'

  ;; have to reset to user values, as we've already changed
  ;; these to plot_data values, and the next event will interpret the
  ;; current settings as the user settings
  tvlct,q.rsave,q.gsave,q.bsave
  !p = puser
  !x = q.xuser
  !y = q.yuser
  device,decomposed = q.decomp
  wset,userwindow

END
