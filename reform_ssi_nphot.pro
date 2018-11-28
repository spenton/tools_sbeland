function reform_ssi_nphot, file=file, data=data
    ; this routine simply reads the SORCE SSI combined Irradiance
    ; daily data and reformats it to a daily spectra in 
    ; number of photons / seconds / m^2 / nm
    ; The data is interpolated to a fixed grid of 1nm (to match TIMED-SEE)

    if n_elements(file) eq 0 then begin
        file='/Users/sbeland/SORCE/data/sorce_ssi_L3.sav'
        ; file=dialog_pickfile()
        ; if strlen(file) eq 0 then return,-1
    endif

    restore, file

    ; get the days covered
    pos = uniq(sorce_ssi_L3.jdn)
    days = sorce_ssi_L3.jdn[pos]
    ndays = n_elements(days)
    ; convert the dates from yyyymmdd to yyyyddd
    yy=double(strmid(strtrim(string(long(sorce_ssi_L3.ymd[pos])),2),0,4))
    mm=double(strmid(strtrim(string(long(sorce_ssi_L3.ymd[pos])),2),4,2))
    dd=double(strmid(strtrim(string(long(sorce_ssi_L3.ymd[pos])),2),6,2))
    new_date=ymd2yd(yy,mm,dd)

    ; now add the daily data from TIMED-SEE from 0 to 155 nm
    restore,'/Users/sbeland/SORCE/data/timed_see_nphotons_flux.sav'
    tsee_jdn = yd2jd(double(timed_see.date)) - 0.5d

    ; wavelengths for SORCE will go from 155 to 600 nm
    ; define wavelengths to be from 0 to 600 nm
    wave = dindgen(600)+1.0d
    nwave = n_elements(wave)


    sorce_ssi = {date:new_date, jdn:sorce_ssi_L3.jdn[pos], wave:wave, $
                 irradiance:fltarr(nwave,ndays), nphotons:fltarr(nwave,ndays)}

    ; store the SORCE data
    for daynum=0L,ndays-1L do begin
        p=where(sorce_ssi_L3.jdn eq days[daynum],count)
        if count eq 0 then continue
        new_w = (sorce_ssi_L3.w1[p]+sorce_ssi_L3.w2[p])/2.0
        wmin = min(new_w,max=wmax)
        wpos=where(wave ge wmin and wave le wmax,complement=cwpos)
        flux = interpol(sorce_ssi_L3.irradiance[p],new_w, wave, /nan)
        flux[cwpos] = 0.0
        ; convert from Watts/m^2/nm to number_of_photons/sec/m^2/nm
        nphotons = flux * 5.03d15 * wave

        ; get the timed_see data for this day
        tsee_flux=[]
        p=where(tsee_jdn eq days[daynum],count)
        if count gt 0 then begin
            tsee_flux = interpol(timed_see.sp_flux[*,p],timed_see.sp_wave, wave[0:154], /nan)
            tsee_nphotons = tsee_flux * 5.03d15 * wave[0:154]
            flux = [tsee_flux[0:154], flux[155:-1]]
            nphotons = [tsee_nphotons[0:154], nphotons[155:-1]]
        endif
        sorce_ssi.irradiance[*,daynum] = flux
        sorce_ssi.nphotons[*,daynum] = nphotons
    endfor

    return, sorce_ssi
end
