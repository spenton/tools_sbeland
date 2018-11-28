function reform_level3, indata

    ; get the list of unique wavelengths form the structure returned by get_level3_spectrum
     sw=sort(indata.minwavelength)
     wv=sw[uniq(indata[sw].minwavelength)]

     ; get the total number of unique days
     sd=sort(indata.NOMINALGPSTIMETAG)
     dd=sd[uniq(indata[sd].NOMINALGPSTIMETAG)]

     irradiance = dblarr(n_elements(dd), n_elements(wv))

     ; populate each wavelength column
     for i=0L,n_elements(wv)-1L do begin
         p=where(indata.minwavelength eq indata[wv[i]].minwavelength,count)
         match, indata[dd].NOMINALGPSTIMETAG, indata[p].NOMINALGPSTIMETAG, sa, sb
         irradiance[sa, i] = indata[p[sb]].irradiance
     endfor

     return, {INSTRUMENTMODEID: indata[0].INSTRUMENTMODEID, $
              VERSION:indata[0].version, $
              TIMESPANINHOURS:indata[0].TIMESPANINHOURS, $
              NOMINALGPSTIMETAG:indata[dd].NOMINALGPSTIMETAG, $
              wavelength:indata[wv].minwavelength, $
              irradiance:temporary(irradiance)}
end
