; ==================================================================================
; This procedure reads in a SUSIM ASCII data file
; ==================================================================================
  function read_susim_ascii,file, julian=julian
; ----------------------------------------------------------------------------------

    if n_elements(file) eq 0 then begin
        file=dialog_pickfile(filter='susim*.ascii')
        if strlen(file) eq 0 then return,-1
    endif

    openr, lun, file, /get_lun
    ; skip the first 19 lines
    line = ''
    FOR i=1,19 DO readf, lun, line

    ; our array of wavelengths
    waves=findgen(296) + 115.5d
    uars_day=[]
    MSD=[]  ; Mean Solar Distance
    irrad=[]

    array = ''
    WHILE NOT EOF(lun) DO BEGIN 
        READF, lun, line 
        values=strsplit(line,' ',/extract)
        if strlen(line) eq 0 then continue
        if n_elements(values) eq 3 then begin
            uars_day=[uars_day,double(values[0])+double(values[1])/86400d]
            MSD = [MSD, double(values[2])]
        endif else if n_elements(values) eq 5 or n_elements(values) eq 1 then begin
            irrad = [irrad, double(values)]
        endif
    ENDWHILE
    free_lun,lun

    ; create a 2D array of irradiances [num_days x num_waves]
    ;   [0,*] is the irradiance for all wavelengths on the first day
    ;   [*,0] is the irradiance for every day for the first wavelength

    irradiance=reform(irrad, n_elements(uars_day), n_elements(waves))
    if keyword_set(julian) then begin
        return, {wavelengths:waves, DATE:ud2jd(uars_day), MEAN_SOLAR_DISTANCE:MSD, irradiance:irradiance}
    endif else begin
        return, {wavelengths:waves, DATE:uars_day, MEAN_SOLAR_DISTANCE:MSD, irradiance:irradiance}
    endelse

end
