;+
; NAME:   SHOW_SIM_TIMELINE
;
; PURPOSE: Extract SORCE SIM data from the database, scale and plot it.
;
; CALLING SEQUENCE:
;   out_data = SHOW_SIM_TIMELINE(startTime, stopTime, 
;                  items=items, simA=simA, simB=simB, 
;                  noplot=noplot, noscale=noscale, ...)
;
; INPUT PARAMETERS:
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
;   indata -
;      The data structure from a previous call to show_sim_timeline
;      or from get_sorce_telemetry can be used as input to plot
;      the data (avoids extracting from database a second time)
;      If a structure is passed, the startTime, stopTime and are
;      ignored as well as the items specification.  Only the plotting
;      keywords are used.
;
;   items -
;      The telemetry item NAMES requested. Must be 
;      used with the simA or simB keyword. 
;      By default the following are extracted:
;          esr_array, uv_array, vis1_array, vis2_array, ir_array,
;          acs_sun_presence, rad_trap_in, shutter_pos, position,
;          fssunangle0, fssunangle1
;
;   simA, simB - 
;      Specifies which SIM channel to retreive.  These are mutually
;      exclusive and defaults to simA.
;
;   noplot -
;      If specified, skips plotting the data on the screen.
;
;   xscale -
;      If specified, will subtract the starting GPS time from the
;      X axis and scale the time in seconds. xscale and xsd are
;      mutually exclusive.
;
;   xsd -
;      If specified, will display the X-axis as mission days instead
;      of GPS micro-seconds. xscale and xsd are mutually exclusive.
;
;   Any additional input keywords will be passed on to get_sorce_telemetry.
;   This includes (see the header of get_sorce_telemetry for details):
;      tlmId
;      gps
;      missionDays
;      julianDays
;      scienceOnly 
;      noScience 
;      verbose 
;      server 
;
;
; OUTPUT PARAMETERS:
;   out_data -
;      The telemetry data retrieved.  This is the structure returned
;      by get_sorce_telemetry.
;
; EXAMPLE:  
;
; COMMON BLOCKS:
;      NONE
;
; NOTE:
;
; REVISION HISTORY:
;   $Id: show_sim_timeline.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
; 
;-
;*****************************************************************************
FUNCTION SHOW_SIM_TIMELINE, startTime, stopTime, indata=indata, items=items, $
         simA=simA, simB=simB, noplot=noplot, xscale=xscale, yscale=yscale, $
         xsd=xsd, _extra=extra

   if n_elements(indata) eq 0 then begin
      ;def_itemName = ['esr_array', 'uv_array', 'vis1_array', 'vis2_array', $
      ;       'ir_array', 'mu_sun_presence', 'rad_trap_in', 'shutter_pos', $
      ;       'position', 'fssunangle0', 'fssunangle1', 'uv_temp', 'prism_drive_temp']
      ;def_tlmId_A = [100058L, 100046L, 100043L, 100055L, $
      ;      100052L, 8025L, 100331L, 102121L, $
      ;      101862L, 7871L, 7872L, 101069L, 101909L]
      ;def_tlmId_B = [100048L, 100057L, 100044L, 100056L, $
      ;      100047L, 8025L, 100582L, 100197L, $
      ;      100478L, 7871L, 7872L, 101580L, 102212L]

      def_itemName = ['uv_array', $
             'ac_sun_present', 'fssunangle0', 'fssunangle1', 'fssunanglevld', $
             'shutter_pos', 'rad_trap_in', 'rad_trap_out',  $
             'position', 'primary', 'uv_temp', 'prism_drive_temp']
      def_tlmId_A = [100046L, $
            7464L, 7871L, 7872L, 7749L, $
            102121L, 100331L, 100038L, $
            101862L, 100328L, 101069L, 101909L]
      def_tlmId_B = [100057L, $
            7464L, 7871L, 7872L, 7749L, $
            100197L, 100582L, 101579L, $
            100478L, 100329L, 101580L, 102212L]

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

      if n_elements(items) eq 0 then begin
          set_default=1
          ; check if tlmId was passed instead
          if size(extra, /tname) eq 'STRUCT' then begin
              if (where(tag_names(extra) eq 'TLMID'))[0] ge 0 then begin
                  set_default=0
              endif
          endif
          if set_default then begin
              if size(extra, /tname) eq 'STRUCT' then begin
                  ; add tlmId to the existing extra structure
                  if simB then $
                     extra = create_struct('TLMID',def_tlmId_B, extra) $
                  else $
                     extra = create_struct('TLMID',def_tlmId_A, extra) 
              endif else begin
                  ; create a new strucrure to contain the telemetry items
                  if simB then $
                     extra = create_struct('TLMID',def_tlmId_B) $
                  else $
                     extra = create_struct('TLMID',def_tlmId_A)
              endelse
          endif 
          get_sorce_telemetry, data, info, startTime, stopTime, _extra=extra

      endif else begin
          if simB then extElem='sim_b' else if simA then extElem='sim_a'
          ; check if externalElement was provided with items
          if size(extra, /tname) eq 'STRUCT' then begin
              pos = where(tag_names(extra) eq 'EXTERNALELEMENT')
              if pos ge 0 then extElem=extra.externalelement
          endif
          get_sorce_telemetry, data, info, startTime, stopTime, externalElement=extElem, $
              itemName=items, _extra=extra
      endelse

      if keyword_set(noplot) then return, data

   endif else begin
       data=indata
   endelse

   ;find the smallest start time
   if size(data,/tname) eq "UNDEFINED" then return,0

   min_time=!VALUES.F_INFINITY
   for i=0,n_tags(data)-1L DO BEGIN
       if (*data.(i)).n_science gt 0 then begin
           min_time = min([min_time,(*data.(i)).science.timetag])
       endif else if (*data.(i)).N_HOUSEKEEPING gt 0 then begin
           min_time = min([min_time,(*data.(i)).housekeeping.timetag])
       endif
   endfor

   ; plot each telemetry items (in seconds)
   tnames = tag_names(data)
   for i=0,n_tags(data)-1L DO BEGIN
      color=0
      psym=0
      if strpos(tnames[i],'RAD_TRAP') ge 0 then psym=-7
      if strpos(tnames[i],'POWER') ge 0 then psym=-7
      if strpos(tnames[i],'UV_TEMP') ge 0 then psym=-7
      if strpos(tnames[i],'SUNPRESENT') ge 0 then begin
          psym=-7
          color=[0,180,0]
      endif
      if strpos(tnames[i],'SHUTTER') ge 0 then begin
          psym=-7
          color=[160,160,160]
      endif
      if (*data.(i)).n_science gt 0 then begin
          ; normalize and offset the data
          xdata = (*data.(i)).science.timetag
          ydata = (*data.(i)).science.dn
          if keyword_set(xscale) then begin
              xdata = (xdata-min_time)/1.0d6 
          endif else if keyword_set(xsd) then begin
              xdata = gps2sd(xdata/1.0d6)
          endif
          if keyword_set(yscale) then begin
              miny=min(ydata, max=maxy)
              if miny ne maxy then begin
                  ydata = 1.0d - (double(maxy)-ydata)/(double(maxy)-double(miny)) + double(i/2.0)
              endif else begin
                  ydata += double(i/2.0)
              endelse
          endif
          lineplot, xdata, ydata, title=tnames[i],psym=psym, color=color, charsize=1.5
      endif else if (*data.(i)).N_HOUSEKEEPING gt 0 then begin
          tnn = tag_names((*data.(i)).housekeeping)
          tnn_pos = where(strpos(tnn,'EU') eq 0)
          xdata = (*data.(i)).housekeeping.timetag
          if keyword_set(xscale) then begin
              xdata = (xdata-min_time)/1.0d6 
          endif else if keyword_set(xsd) then begin
              xdata = gps2sd(xdata/1.0d6)
          endif
          if tnn_pos[0] ge 0 then $
             ydata = (*data.(i)).housekeeping.eu $
          else $
             ydata = (*data.(i)).housekeeping.dn
          if keyword_set(yscale) then begin
              miny=min(ydata, max=maxy)
              if miny ne maxy then begin
                  ydata = 1.0d - (double(maxy)-ydata)/(double(maxy)-double(miny)) + double(i/2.0)
              endif else begin
                  ydata += double(i/2.0)
              endelse
          endif 
          lineplot, xdata, ydata, title=tnames[i], psym=psym, color=color, charsize=1.5
      endif

   endfor

   return,data

END 
