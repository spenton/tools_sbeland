;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Compare the specified SimCalibratedIrradiance with our reference  
;   spectra for the requested instrumentModeId and estimate a new aleph
;   to match the absolute irradiance from our reference day (453.67)
;
; CALLING SEQUENCE:
;
; INPUT PARAMETERS:
;
; OPTIONAL INPUT PARAMETERS:
;   version -
;      Requests a specific version of the SimCalibratedIrradiance, otherwise
;      pick the latest version valid for the specified time range.
;
; RETURNED PARAMETERS:
;
; OPTIONAL OUTPUT PARAMETERS:
;
; EXAMPLE:  
;
; COMMON BLOCKS:
;      NONE
;
; NOTES:
;
; REVISION HISTORY:
;   Revision: $Id: get_aleph.pro,v 1.2 2018/11/26 16:16:02 spenton Exp $
;------------------------------------------------------------------
function get_aleph, instrumentModeId, spect=spect, version=version, noplot=noplot, $
    dburl=dburl, dbdriver=dbdriver, user=user, password=password

    if instrumentModeId ge 41 then begin
        t0=453.67840d
        t1=453.69524d
    endif else begin
        ; for ESR we're only interested in the SolarIRScan 
        t0=454.21794d
        t1=454.25990d
    endelse

    if n_elements(user) eq 0 then     user='sbeland'
    if n_elements(password) eq 0 then password='sbeland1'
    if n_elements(dburl) eq 0 then    dbUrl='jdbc:sybase:Tds:lasp-db-dev:4100/SORCE_SIM_DEV'
    if n_elements(dbdriver) eq 0 then dbDriver='com.sybase.jdbc3.jdbc.SybDriver'

    ;jstmt = fjava_get_jdbc_statement(user=user, password=password, dbdriver=dbdriver, dburl=dburl) 
    ;oJavaDbExchange = OBJ_NEW("oJava_DbExchange", jstmt)

    if n_elements(version) eq 0 then version=2541
    version22=24

    irradColumns22=['Wavelength','SimCalibratedIrradiance']
    irradColumns=['Wavelength','SimCorrectedIrradiance']

    ; our reference SimCalibratedIrradiance comes from SIMA ONLY
    if instrumentModeID lt 45 then begin
        spect22_db=get_science_product(irradColumns22, t0, t1, instrumentModeId, /mission,version=version22)
    endif else begin
        ; the CalibratedIrradiances for SimA and SimB are to adjusted in version 22
        ; our reference is SimA and we'll need to interpolate SimA data to SimB wavelength grid
        spect22_dba=get_science_product(irradColumns22, t0, t1, instrumentModeId-4, /mission,version=version22)
        spect22_db=get_science_product(irradColumns22, t0, t1, instrumentModeId, /mission,version=version22)
        spect22_db.irradiance = interpol(spect22_dba.irradiance, spect22_dba.wavelengthref, spect22_db.wavelengthref)
    endelse

    ; get the data structure from get_sim_spectra which align_spectra2 expects
    spect22 = get_sim_spectra(t0, t1, instrumentModeId, /missionDays, version=version22, /noCcdCorr, /noTempCorr, /dn_only)

    ; we need to re-align version 22 with the better Kuruckz alignment of version 19
    ; fill in the get_sim_spectra structure with the corresponding wavelength and irradiance from DB
    match,spect22_db.MICROSECONDSSINCEGPSEPOCH, spect22.timetag, suba, subb
    spect22 = spect22[subb]
    spect22.wavelength = spect22_db[suba].wavelengthref
    spect22.irradiance = spect22_db[suba].irradiance

    ; align with our standard reference spectra
    spect22 = align_spectra2(t0,t1,instrumentModeId, /mission, spect=spect22, /no_plot)

    ; our SimCorrectedIrradiance we want to match to our reference
    if size(spect,/tname) ne 'STRUCT' then begin
        my_spect=get_science_product(irradColumns, t0, t1, instrumentModeId, /mission,version=version,$
            user=user, password=password,dburl=dburl, dbdriver=dbdriver)
        my_wavelength = my_spect.wavelengthref
        my_irradiance = my_spect.irradiance
    endif else begin
        ; check if the input structure comes from compare_19_20 and process_uncorr
        p=where(strpos(tag_names(spect),'SPECT20') ge 0,getit)
        if getit gt 0 then begin
           p=where(spect.SPECT20[*].timestamp[20] gt sd2gps(t0)*1d6,count)
           if count eq 0 then begin
               print,'Error:  input spectra does not cover expected time range'
               return,-1
           endif
            my_wavelength = spect.SPECT20[p[0]].wavelength
            my_irradiance = spect.SPECT20[p[0]].irradiance
        endif else begin
            ; make a copy of the input spectra
            my_wavelength = spect.wavelength
            my_irradiance = spect.irradiance
        endelse
    endelse

    ; interpolate to the new wavelength grid
    p = where(my_wavelength gt 0d)
    irrad22=interpol(spect22.irradiance, spect22.wavelength, my_wavelength[p])
    aleph = irrad22 / my_irradiance[p]

    if NOT keyword_set(noplot) then begin
       lineplot,spect22.wavelength, spect22.irradiance,title='Reference Spectrum @453 version='+strtrim(string(version22),2), $
           xtitle='Wavelength (nm)',ytitle='Irradiance',charsize=1.4
       lineplot,my_wavelength[p], my_irradiance[p],title='Spectrum @453 version='+strtrim(string(version),2)
       plot_multi,my_wavelength[p], aleph, my_wavelength[p], (gauss_smooth(aleph,3.0,/edge_trunc))>0d, $
           xtitle='Wavelength (nm)', ytitle='Aleph 20',/xst,/yst,charsize=1.4, psym=[-4,-3], thick=[1.0,3.0], $
           title='Aleph for Mode='+strtrim(string(InstrumentModeId),2)+' Version='+strtrim(string(version),2)
    endif

    s=sort(my_wavelength[p])
    if instrumentModeId ne 44 and instrumentModeId ne 48 then aleph=(gauss_smooth(aleph[s],3.0,/edge_trunc))>0d else aleph=aleph[s]>0d
    return, {instrumentModeId:instrumentModeId, wavelength:my_wavelength[p[s]], aleph:aleph}

end
