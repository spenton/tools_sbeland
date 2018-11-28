;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Extract the degradation information from the ESR table scans matching
;   the data for the requested wavelength taken at the "same" time for 
;   ESRA and ESRB.  This is done for two consecutive table scans of ESRB
;   to get a running slope of the degradation with mission time.
;
; CALLING SEQUENCE:
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
; OPTIONAL INPUT PARAMETERS:
;   version -
;      Requests a specific version of the SimProfileIntegralCal, otherwise
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
;    The table used for the irradiance is expected to have the final degradation
;    correction applied to the ESR data but no raypath or diode degradation applied to
;    the photodiode data (as in version 15 of our development database for version 19).
;
;
; REVISION HISTORY:
;   Revision: $Id: adjust_esrab_453.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;------------------------------------------------------------------
function adjust_esrab_453, version=version, $
    dburl=dburl, dbdriver=dbdriver, user=user, password=password, $
    gps=gps, missionDays=missionDays, julianDays=julianDays, $
    esrdata=esrdata


    if n_elements(version) eq 0 then version=21
    if n_elements(user) eq 0 then     user='sbeland'
    if n_elements(password) eq 0 then password='sbeland1'
    if n_elements(dburl) eq 0 then    dbUrl='jdbc:sybase:Tds:lasp-db-dev:4100/SORCE_SIM_V19'
    if n_elements(dbdriver) eq 0 then dbDriver='com.sybase.jdbc3.jdbc.SybDriver'

    ; average the delta from day 453 and 491 (to cleanup the noise at wave<500nm)
    t0=[453.0d, 491.0, 673.0, 852.0, 1034.0, 1223.0, 1419.0, 1607.5, 1790.1, 1972.0, 2160.5, 2342.1, 2525.0, 2728.2]
    t1=[454.0d, 492.1, 674.1, 853.1, 1035.1, 1224.0, 1421.2, 1609.6, 1792.3, 1974.2, 2162.5, 2344.3, 2527.1, 2730.4]

    ; limit the fit from 265 to 2950 nm
    waves=dindgen(2950.0-265.0)+265d
    ; but extend the wavelength coverage from 200 to 3000 nm
    all_wave=[dindgen(65)+200d,waves,dindgen(51)+2950d]
    delta_fit=dblarr(n_elements(t0),n_elements(waves))

    ; compare every ESRFullScan
    for tt=0,n_elements(t0)-1 do begin

        esra=get_science_product(['Wavelength','SimCorrectedIrradiance'],t0[tt],t1[tt],31,/mission,$
            version=version,user=user,password=password,dburl=dburl,dbdriver=dbdriver)
        sa=sort(esra.wavelengthref)
        esra=esra[sa]
        esrb=get_science_product(['Wavelength','SimCorrectedIrradiance'],t0[tt],t1[tt],32,/mission,$
            version=version,user=user,password=password,dburl=dburl,dbdriver=dbdriver)
        sb=sort(esrb.wavelengthref)
        esrb=esrb[sb]

        yfitb=smooth(esrb.irradiance,2)
        yfita=interpol(smooth(esra.irradiance,2),esra.wavelengthref,esrb.wavelengthref,/spline)
        md = strtrim(string(t0[tt],format='(I)'),2)
        ;lineplot,esrb.wavelengthref, yfita-yfitb, xtitle='Wavelength (nm)',ytitle='Delta Irradiance esra-esrb',$
        ;    title='Smooth(esra-esrb,2) @ '+md, ptitle='Difference in SimCorrectedIrradiance between ESRA & ESRB', charsize=1.4

        p=where(esrb.wavelengthref lt 865.0d,comp=cp)
        resistant_mean,yfita[p]-yfitb[p],3.0,avg0,goodvec=k0
        resistant_mean,yfita[cp]-yfitb[cp],4.0,avg0,goodvec=k1
        k=[p[k0],cp[k1]]

        delta_irrad = yfita[k]-yfitb[k]
        ;lineplot,esrb[k].wavelengthref, delta_irrad, title='esra-esrb Sigma>3.0 reject @ '+md

        ssetb=bspline_iterfit(esrb[k].wavelengthref, smooth(delta_irrad,10),maxiter=10,requiren=10,bkspace=5)
        delta_fit[tt,*]=bspline_valu(waves,ssetb)
        lineplot,waves, delta_fit[tt,*], title='esra-esrb @ '+md, xtitle='Wavelength (nm)',ytitle='Delta Irradiance esra-esrb',$
            ptitle='Difference in SimCorrectedIrradiance between ESRA & ESRB', charsize=1.4

    endfor

    ; pad the start and end wavelength range with constant values
    delta_irrad = total(delta_fit,1) / n_elements(t0)
    delta_irrad=[replicate(delta_irrad[0],65), delta_irrad, replicate(delta_irrad[-1],51)]
    lineplot,all_wave, delta_irrad, title='Average(esra1,esrb1)'

    return,{wavelength:all_wave, delta_irrad:delta_irrad}

end
