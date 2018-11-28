;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Align the irradiance at each of the OBC segments for the provided
;   time series Irradiance versus mission day.
;   The irradiance is adjusted to the levels for the time range from 453 to 1570.
;   The regions to measure the median irradiance were visually identified (by Jerry H.)
;   close to the OBC where we have low solar activity.
;
; CALLING SEQUENCE:
;   ALIGN_IRRAD, timestamps, irradiance
;
; INPUT PARAMETERS:
;   timestamps -
;      Array of times.
;
;   irradiance -
;      Array of irradiance to adjust (same dimension as timestamps).
;      The input array is modifed inplace to reflect the alignment.
;
; OPTIONAL INPUT PARAMETERS:
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
;
; OPTIONAL OUTPUT PARAMETERS:
;   offsets -
;      Array of structure containing the list of indices to which each 
;      offset value was applied.
;
;   trend -
;      Perorm a polynomial fit to the data before and after the OBC and
;      align the 2 curves in irradiance thus aligning the trends in the 
;      data over the time spand between each OBC.
;
; COMMON BLOCKS:
;      NONE
;
; NOTES:
;
; REVISION HISTORY:
;   Revision: $Id: align_irrad.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;------------------------------------------------------------------
;
pro align_irrad, timestamps, irradiances, $
         gps=gps, missionDays=missionDays, julianDays=julianDays, $
         offsets=offsets, trend=trend


    ; Prestart  prestop   poststart   poststop   jump   correctDirection
    values=[$
    [166.2d, 176.4, 181.7d, 199.6, 178.0], $
    [181.7d, 199.6, 203.7, 221.6, 201.5], $
    ;[261.4, 266.9, 290.1, 295.5, 277.6], $
    [426.4, 444.1, 464.2, 482.5, 449.5], $
    ;[1540.6, 1567.0, 1581.1, 1628.6, 1569.8], $
    ;[1540.6, 1567.0, 1595.0, 1628.6, 1569.8], $
    ;[1542.6, 1551.0, 1578.8, 1592.0, 1570.0], $
    ;[1542.6, 1551.0, 1601.0, 1613.0, 1570.0], $
    ;[1532.5, 1539.3, 1605.5, 1612.0, 1570.0], $
    ;[1510.0, 1540.0, 1586.0, 1605.0, 1570.0], $
    [1533.0, 1548.0, 1579.0, 1587.0, 1570.0], $
    ;[2147.1, 2170.0, 2191.2, 2226.1, 2173.4], $
    ;[2130.0, 2147.0, 2191.2, 2208.1, 2173.4], $
        ;[2130.0, 2142.0, 2283.8, 2294.3, 2173.4], $
        ;[2128.0, 2146.5, 2183.6, 2199.2, 2173.4], $
        ;[2133.0, 2144.0, 2216.0, 2225.0, 2173.4], $
        [2133.0, 2144.0, 2257.0, 2273.0, 2173.4], $
        ;[2133.0, 2144.0, 2191.0, 2200.0, 2173.4], $
    [2401.2, 2448.8, 2471.1, 2501.2, 2478.0], $
    [2776.0, 2799.5, 2804.0, 2838.0, 2803.0], $
    [2804.0, 2838.0, 2839.8, 2893.8, 2838.0], $
    [2839.8, 2893.8, 2899.0, 2926.9, 2894.0], $
    [2907.0, 2921.0, 2974.0, 2986.0, 2927.0], $
    ;[3019.7, 3027.2, 3031.8, 3035.7, 3029.1], $
    ;[3012.7, 3025.0, 3093.2, 3104.5, 3029.1], $
    [3019.0, 3030.0, 3092.0, 3102.0, 3029.0], $
                            ;[3015.0, 3027.0, 3096.0, 3105.0, 3030.0], $
                            ;[3096.0, 3105.0, 3240.0, 3253.0, 3150.0], $
                            ;[3114.0, 3136.0, 3169.0, 3186.0, 3150.0], $
                            [3120.0, 3141.0, 3173.0, 3183.0, 3150.0], $   ; commented out SBeland 2017-01-14
    ;[3120.0, 3133.0, 3158.0, 3164.8, 3150.0], $
                            ;[3480.0, 3493.0, 3707.0, 3720.0, 3577.0], $
                            ;[3373.0, 3380.0, 3611.0, 3616.0, 3577.0], $
                            ;[3528.0, 3536.0, 3622.0, 3628.0, 3577.0], $
                            ;[3614.0, 3637.0, 4193.0, 4220.0, 3900.0], $
                            ;[3790.0, 3807.0, 4135.0, 4164.0, 3900.0], $
                            ;[3752.0, 3788.0, 4070.0, 4107.0, 3900.0], $
                            [3506.0, 3514.0, 3667.0, 3675.0, 3595.0], $     ; commented out SBeland 2017-01-14
                            [3656.0, 3664.0, 3709.0, 3718.0, 3679.0], $     ; commented out SBeland 2017-01-14
                            [3752.0, 3788.0, 4076.0, 4109.0, 3900.0]]     ; commented out SBeland 2017-01-14
    ;[3131.6, 3145.3, 3157.3, 3177.4, 3150.0]]

    values_jerry = [$
    [166.2d, 176.4, 181.7d, 199.6, 178.0], $
    [181.7d, 199.6, 203.7, 221.6, 201.5], $
    [426.4, 444.1, 464.2, 482.5, 449.5], $
    [1532.0,1550.0,1604.0,1621.0,1571.0],$
    [2133.0, 2144.0, 2257.0, 2273.0, 2173.4], $
    [2401.2, 2448.8, 2471.1, 2501.2, 2478.0], $
    [2630.0,2663.0,2818.17,2877.39,2802.], $
    [2804.0, 2838.0, 2839.8, 2893.8, 2838.0], $
    [2847.47,2877.75,2903.60,2914.91,2893.], $
    [2905.56, 2915.06, 2948.95, 2957.82, 2926.], $
    [3019.0, 3030.0, 3092.0, 3102.0, 3029.0], $
    [3120.0, 3141.0, 3173.0, 3183.0, 3150.0], $   
    [3506.0, 3514.0, 3667.0, 3675.0, 3595.0], $   
    [3645.7,3675.7, 3706.25, 3716., 3684], $
    [3759.54, 3817.82, 4141.91, 4194.50, 3845.], $
    [3752.0, 3788.0, 4076.0, 4109.0, 3900.0], $
    [4193.33, 4202.33, 4270., 4285., 4217.85]]

    ; for testing purposes try out Jerry's values
    values=values_jerry

    pos2803  = reform(where(values[4,*] eq 2803.0d))

    if keyword_set(missionDays) then begin
        ; do nothing here since already in mission days
        day453 = 449.5
    endif else if keyword_set(julianDays) then begin
        ; user specified time in julian days
        values = sd2jd(values)
        day453 = sd2jd(449.5d)
    endif else begin
        ; user specified time in gps days (default)
        values = sd2gps(values)*1.d6
        day453 = sd2gps(449.5d)*1d6
    endelse

    left0    = reform(values[0,*])
    left1    = reform(values[1,*])
    right0   = reform(values[2,*])
    right1   = reform(values[3,*])
    obctimes = reform(values[4,*])
    pos453   = reform(where(obctimes eq day453))

    obc_before = where(obctimes le day453,complement=obc_after,count0)
    s = reverse(sort(obctimes[obc_before]))

    all_pos=[]
    all_delta=[]
    for i=0,count0-1 do begin
        pos = where(timestamps le obctimes[obc_before[s[i]]],comp=cp,count)
        if count eq 0 then break

        if NOT keyword_set(trend) then begin
            medleft  = where(timestamps ge left0[obc_before[s[i]]]  and timestamps le left1[obc_before[s[i]]],countl)
            if countl eq 0 then begin
                ; double the width of the search area
                new_left0 = left0[obc_before[s[i]]] - (left1[obc_before[s[i]]] - left0[obc_before[s[i]]])
                medleft  = where(timestamps ge new_left0  and timestamps le left1[obc_before[s[i]]],countl)
            endif
            medright = where(timestamps ge right0[obc_before[s[i]]] and timestamps le right1[obc_before[s[i]]],countr)
            if countr eq 0 then begin
                ; double the width of the search area
                new_right1 = right1[obc_before[s[i]]] + (right1[obc_before[s[i]]] - right0[obc_before[s[i]]])
                medright = where(timestamps ge right0[obc_before[s[i]]] and timestamps le new_right1,countr)
            endif
            if countl gt 0 and countr gt 0 then begin
                delta = median(irradiances[medleft]) - median(irradiances[medright])
            endif else delta=0d
        endif else begin
            ;medright = where(timestamps ge obctimes[obc_before[s[i]]] and timestamps le obctimes[obc_before[s[i]]+1],countr)
            ;medleft  = where(timestamps le obctimes[obc_before[s[i]]] and timestamps ge obctimes[obc_before[s[i]]-1],countl)
            medright = where(timestamps ge right0[obc_before[s[i]]] and timestamps le right1[obc_before[s[i]]],countr)
            medleft  = where(timestamps le left0[obc_before[s[i]]] and timestamps ge left1[obc_before[s[i]]],countl)
            if countl gt 1 and countr gt 1 then begin
                cright = robust_poly_fit(timestamps[medright], irradiances[medright], 2)
                cleft  = robust_poly_fit(timestamps[medleft], irradiances[medleft], 2)
                delta = poly(obctimes[obc_before[s[i]]],cleft) - poly(obctimes[obc_before[s[i]]],cright)
            endif else delta=0d
        endelse
        irradiances[pos] -= delta
        all_pos = [all_pos, pos]
        all_delta = [all_delta, replicate(-delta,n_elements(pos))]
    endfor
   

    s = sort(obctimes[obc_after])
    for i=0,n_elements(obc_after)-1 do begin
        pos = where(timestamps gt obctimes[obc_after[s[i]]],count)
        if count eq 0 then break

        if keyword_set(trend) and obctimes[obc_after[s[i]]] ge obctimes[pos2803] then begin
            ; the trends after day 2800 are not clear so ignore
            delta = 0d
        endif else if NOT keyword_set(trend) or i eq n_elements(obc_after)-1 then begin
            medleft  = where(timestamps ge left0[obc_after[s[i]]]  and timestamps le left1[obc_after[s[i]]],countl)
            if countl eq 0 then begin
                ; double the width of the search area
                new_left0 = left0[obc_after[s[i]]] - (left1[obc_after[s[i]]] - left0[obc_after[s[i]]])
                medleft  = where(timestamps ge new_left0  and timestamps le left1[obc_after[s[i]]],countl)
            endif
            medright = where(timestamps ge right0[obc_after[s[i]]] and timestamps le right1[obc_after[s[i]]],countr)
            if countr eq 0 then begin
                ; double the width of the search area
                new_right1 = right1[obc_after[s[i]]] + (right1[obc_after[s[i]]] - right0[obc_after[s[i]]])
                medright = where(timestamps ge right0[obc_after[s[i]]] and timestamps le new_right1,countr)
            endif
            if countl gt 0 and countr gt 0 then begin
                delta = median(irradiances[medleft]) - median(irradiances[medright])
            endif else delta=0d
        endif else begin
            ;medright = where(timestamps ge obctimes[obc_after[s[i]]] and timestamps le obctimes[obc_after[s[i]]+1],countr)
            ;medleft  = where(timestamps le obctimes[obc_after[s[i]]] and timestamps ge obctimes[obc_after[s[i]]-1],countl)
            medright = where(timestamps ge right0[obc_after[s[i]]] and timestamps le right1[obc_after[s[i]]],countr)
            medleft  = where(timestamps le left0[obc_after[s[i]]] and timestamps ge left1[obc_after[s[i]]],countl)
            if countl gt 1 and countr gt 1 then begin
                cright = robust_poly_fit(timestamps[medright], irradiances[medright], 2)
                cleft  = robust_poly_fit(timestamps[medleft], irradiances[medleft], 2)
                delta = poly(obctimes[obc_after[s[i]]],cleft) - poly(obctimes[obc_after[s[i]]],cright)
            endif else delta=0d
        endelse
        irradiances[pos] += delta
        all_pos = [all_pos, pos]
        all_delta = [all_delta, replicate(-delta,n_elements(pos))]
    endfor

    if n_elements(all_pos) eq 0 then offsets=-1 else begin
        s = sort(all_pos)
        offsets={pos:all_pos[s], offset:all_delta[s]}
    endelse
 
    return

end
