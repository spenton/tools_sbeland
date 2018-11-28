; $Id: plot_multi.pro,v 1.2 2018/11/26 16:16:02 spenton Exp $
; Routine to plot up to 20 lines on one plot.
;
; Assumes the inputs plot_multi,x0,y0,x1,y1,x2,y2,x3,y3
;
;  Written by S. Beland
;  Modified to use pointer arrays instead of structures for assembling data;
;  passes keywords to PLOT, OPLOT
;  uses FSC_COLOR routine from
;  http://dfanning.com for display-independent color.  N.J. Cunningham, Feb 2008
; 
;  Modified to wrap plot_data so as to be called by itself instead as the 
;  string parameter to plot_data.
; -----------------------------------------------------------
pro my_plot_multi, x0,y0, x1,y1, x2,y2, x3,y3, x4,y4, x5,y5, x6,y6, $
                x7,y7,x8,y8,x9,y9,x10,y10,x11,y11,x12,y12,x13,y13,x14,y14, $
                x15,y15,x16,y16,x17,y17,x18,y18,x19,y19, $
                xrange=xrange,yrange=yrange,xstyle=xstyle, $
                ystyle=ystyle,colors=colors, labels=labels, psym=psym, symsizes=symsizes, $
                linestyle=linestyle, thick=thick, nsum=nsum, _extra=extra

   ; Error handling.

On_Error, 2
Catch, theError
IF theError NE 0 THEN BEGIN
   Catch, /Cancel
   message, /reissue_last
   RETURN
ENDIF

  np = n_params()
  if np eq 0 then return
  if np eq 1 then begin
     nlines = 1
     xptr = (yptr = ptrarr(1))
     yptr[0] = ptr_new(x0)
     xptr[0] = ptr_new(indgen(n_elements(x0)))
     minx=min(*xptr[0],max=maxx)
     miny=min(*yptr[0],max=maxy)
  endif else begin
     nlines = np/2
     xptr = (yptr = ptrarr(nlines))

     minx = min(x0,max=maxx)
     miny = min(y0,max=maxy)

     switch nlines-1 of
        19: begin
           minx = min([x19,minx,maxx],max=maxx)
           miny = min([y19,miny,maxy],max=maxy)
           xptr[19] = ptr_new(x19)
           yptr[19] = ptr_new(y19)
        end
        18: begin
           minx = min([x18,minx,maxx],max=maxx)
           miny = min([y18,miny,maxy],max=maxy)
           xptr[18] = ptr_new(x18)
           yptr[18] = ptr_new(y18)
        end
        17: begin
           minx = min([x17,minx,maxx],max=maxx)
           miny = min([y17,miny,maxy],max=maxy)
           xptr[17] = ptr_new(x17)
           yptr[17] = ptr_new(y17)
        end
        16: begin
           minx = min([x16,minx,maxx],max=maxx)
           miny = min([y16,miny,maxy],max=maxy)
           xptr[16] = ptr_new(x16)
           yptr[16] = ptr_new(y16)
        end
        15: begin
           minx = min([x15,minx,maxx],max=maxx)
           miny = min([y15,miny,maxy],max=maxy)
           xptr[15] = ptr_new(x15)
           yptr[15] = ptr_new(y15)
        end
        14: begin
           minx = min([x14,minx,maxx],max=maxx)
           miny = min([y14,miny,maxy],max=maxy)
           xptr[14] = ptr_new(x14)
           yptr[14] = ptr_new(y14)
        end
        13: begin
           minx = min([x13,minx,maxx],max=maxx)
           miny = min([y13,miny,maxy],max=maxy)
           xptr[13] = ptr_new(x13)
           yptr[13] = ptr_new(y13)
        end
        12: begin
           minx = min([x12,minx,maxx],max=maxx)
           miny = min([y12,miny,maxy],max=maxy)
           xptr[12] = ptr_new(x12)
           yptr[12] = ptr_new(y12)
        end
        11: begin
           minx = min([x11,minx,maxx],max=maxx)
           miny = min([y11,miny,maxy],max=maxy)
           xptr[11] = ptr_new(x11)
           yptr[11] = ptr_new(y11)
        end
        10: begin
           minx = min([x10,minx,maxx],max=maxx)
           miny = min([y10,miny,maxy],max=maxy)
           xptr[10] = ptr_new(x10)
           yptr[10] = ptr_new(y10)
        end
        9: begin
           minx = min([x9,minx,maxx],max=maxx)
           miny = min([y9,miny,maxy],max=maxy)
           xptr[9] = ptr_new(x9)
           yptr[9] = ptr_new(y9)
        end
        8: begin
           minx = min([x8,minx,maxx],max=maxx)
           miny = min([y8,miny,maxy],max=maxy)
           xptr[8] = ptr_new(x8)
           yptr[8] = ptr_new(y8)
        end
        7: begin
           minx = min([x7,minx,maxx],max=maxx)
           miny = min([y7,miny,maxy],max=maxy)
           xptr[7] = ptr_new(x7)
           yptr[7] = ptr_new(y7)
        end
        6: begin
           minx = min([x6,minx,maxx],max=maxx)
           miny = min([y6,miny,maxy],max=maxy)
           xptr[6] = ptr_new(x6)
           yptr[6] = ptr_new(y6)
        end
        5: begin
           minx = min([x5,minx,maxx],max=maxx)
           miny = min([y5,miny,maxy],max=maxy)
           xptr[5] = ptr_new(x5)
           yptr[5] = ptr_new(y5)
        end
        4: begin
           minx = min([x4,minx,maxx],max=maxx)
           miny = min([y4,miny,maxy],max=maxy)
           xptr[4] = ptr_new(x4)
           yptr[4] = ptr_new(y4)
        end
        3: begin
           minx = min([x3,minx,maxx],max=maxx)
           miny = min([y3,miny,maxy],max=maxy)
           xptr[3] = ptr_new(x3)
           yptr[3] = ptr_new(y3)
        end
        2: begin
           minx = min([x2,minx,maxx],max=maxx)
           miny = min([y2,miny,maxy],max=maxy)
           xptr[2] = ptr_new(x2)
           yptr[2] = ptr_new(y2)
        end
        1: begin
           minx = min([x1,minx,maxx],max=maxx)
           miny = min([y1,miny,maxy],max=maxy)
           xptr[1] = ptr_new(x1)
           yptr[1] = ptr_new(y1)
        end
        0: begin
           xptr[0] = ptr_new(x0)
           yptr[0] = ptr_new(y0)
        end
     endswitch
  endelse

  colorlist = ['red', 'green', 'blue', 'orange', 'yellow', $
                'magenta', 'cyan', 'dark green', 'sienna', $
                'purple', 'dark red', 'burlywood', 'lime green', $
                'maroon', 'tomato', 'slate blue', 'dark goldenrod' ,'olive', $
                'dodger blue', 'dark grey']
  if n_elements(colors) eq 0 then begin
     ;colors = fsc_color(colorlist)
     colors = cgcolor(colorlist)
  endif else if n_elements(colors) ne nlines then begin
     ;colors = [colors,fsc_color(colorlist)]
     colors = [colors,cgcolor(colorlist)]
  endif
;;   if n_elements(colors) eq 0 then colors=indgen(nlines)+2

  if n_elements(xrange) eq 0 then xrange=[minx,maxx]
  if n_elements(yrange) eq 0 then yrange=[miny,maxy]
  if n_elements(xstyle) eq 0 then xstyle=0
  if n_elements(ystyle) eq 0 then ystyle=0
  if n_elements(linestyle) gt 1 then begin
     setlines = 1
     linest = linestyle
  endif else setlines=0
  if n_elements(thick) gt 1 then begin
     linethick = thick
  endif else linethick=replicate(0,nlines)
  case n_elements(psym) of
     0: symarr=replicate(-4,nlines)
     1: symarr = replicate(psym,nlines)
     else: symarr = psym
  endcase
  ;p=where(symarr eq 10,count)
  ;if count gt 0 then symarr[p]=-3

  case n_elements(nsum) of
     0: n_sum = replicate(1,nlines)
     1: n_sum = replicate(nsum,nlines)
     else: n_sum = nsum
  endcase
  
  plot,*xptr[0],*yptr[0], /nodata, xrange=xrange, $ 
       yrange=yrange, xstyle = xstyle, ystyle = ystyle, _extra = extra

  if not keyword_set(symsizes) then symsizes=replicate(1.0,nlines)

  for line=0,nlines-1 do begin
     if setlines eq 1 then begin 
        oplot,*xptr[line],*yptr[line], color=colors[line],psym=symarr[line], thick=linethick[line], $
              linestyle = linest[line], nsum = n_sum[line], symsize=symsizes[line], _extra = extra 
     endif else begin
        oplot,*xptr[line],*yptr[line], color=colors[line],psym=symarr[line], thick=linethick[line], $
              nsum = n_sum[line], symsize=symsizes[line], _extra = extra
     endelse
  endfor

  ; will print the labels in the order provided.
  ; if the number of labels is different that the number of datasets, it is assumed
  ; that the missing labels are for the last datasets
  nlabels = n_elements(labels)
  if nlabels gt 0 and nlabels le nlines then begin
     ; labels were provided
     xpos=xrange[0]+(xrange[1]-xrange[0])*0.8
     ypos=yrange[0]+(yrange[1]-yrange[0])*0.1
     label_psym= symarr[0:nlabels-1]
     p=where(label_psym eq 10, count)
     if count gt 0 then label_psym[p]=-3
     if setlines eq 1 then $
        al_legend,labels,psym=label_psym,colors=byte(colors[0:nlabels-1]),/box,$
                linestyle=linest[0:nlabels-1], symsize=symsizes[0:nlabels-1],_extra=extra $
     else $
        al_legend,labels,psym=label_psym,colors=byte(colors[0:nlabels-1]),/box,$
                symsize=symsizes[0:nlabels-1], _extra=extra
  endif

  ptr_free, xptr, yptr

end

; -----------------------------------------------------------
;
; This is just a warpper routine for plot_data to display multiple plots
pro plot_multi, d1, d2, d3, d4, d5, d6, d7, d8, $
               d9, d10, d11, d12, d13, d14, d15, d16, d17, d18, $
               d19, d20, d21, d22, d23, d24, d25, d26, d27, d28, d29, d30, $
               d31, d32, d33, d34, d35, d36, d37, d38, d39, d40, $
               WTITLE=wtitle, WXSIZE=wxsize, WYSIZE=wysize, GROUP=group, $
               OUTNAME=outname, border=border, xdate=xdate,_extra=extra

    resolve_routine, 'plot_data', /no_recompile

    plot_data, 'my_plot_multi',  d1, d2, d3, d4, d5, d6, d7, d8, $
               d9, d10, d11, d12, d13, d14, d15, d16, d17, d18, $
               d19, d20, d21, d22, d23, d24, d25, d26, d27, d28, d29, d30, $
               d31, d32, d33, d34, d35, d36, d37, d38, d39, d40, $
               WTITLE=wtitle, WXSIZE=wxsize, WYSIZE=wysize, GROUP=group, $
               OUTNAME=outname, border=border, xdate=xdate,_extra=extra


end

