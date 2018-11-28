;+
; NAME:   GET_PDTAU2
;
; PURPOSE: 
;    This routine uses the photodiode spectrum from SIMA and SIMB for the
;    first part of the mission to determine F-Function covering the specific
;    spectral region.
;
; CALLING SEQUENCE:
;    result=get_pdtau2(instrumetModeId, inspectrum=inspectrum, solar54=solar54, solar55=solar55)
;
; OPTIONAL INPUT PARAMETERS:
;    none
;
; OPTIONAL INPUT KEYWORDS:
;    none 
;
; OUTPUT PARAMETERS:
;    none
;
; OPTIONAL OUTPUT PARAMETERS:
;    NONE
;
; EXAMPLE:  
;
; COMMON BLOCKS:
;    NONE
;
; NOTES:
;
; REVISION HISTORY:
;   Revision: $Id: get_pdtau2.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;-
;***************************************************************************************************

FUNCTION GET_PDTAU2, t0, t1, visdiode=visdiode, uvdiode=uvdiode, stepsize=stepsize, $
    solar54=solar54, solar55=solar55, inspectrum=inspectrum, noplot=noplot, coeffs=coeffs, $
    afact=afact, kappa=kappa, waverange=waverange, surfacefit=surfacefit, wavezero=wavezero, $
    getdegcol=getdegcol, degcol_c54=degcol_c54, degcol_c55=degcol_c55, gettau=gettau

    if n_elements(stepsize) eq 0 then stepsize=3d
    if n_elements(wavezero) eq 0 then wavezero=2500d

    if keyword_set(visdiode) then begin
        if n_elements(solar54) eq 0 or n_elements(solar55) eq 0 then begin
            ;restore,'~/SORCE/data/solarexp_54_55.sav'
            ;restore,'~/SORCE/data/solarexp_54_55_modified_v20.sav'
            ;restore,'~/SORCE/data/solarexp_54_55_mod_407nm_4371.sav'
            restore,'~/SORCE/data/solarexp_54_55_mod_407nm.sav'
        endif

        if size(afact,/tname) ne 'STRUCT' then begin
            ; get the vis diode raypath from the raytrace data
            ;readcol,'~/SORCE/data/overlap_vis1.txt',pd_aw,pd_apath,format='(d,d)',/silent
            ;readcol,'~/SORCE/data/overlap_vis1a_v20_reduced55.txt',pd_aw,pd_apath,format='(d,d)',/silent
            ;just as a test multiply the raytrace raypath by 4.8
            ;pd_apath *= 4.8d
            ;readcol,'~/SORCE/data/overlap_vis1_smooth_mod_v20.txt',pd_aw,pd_apath,format='(d,d)'
            readcol,'~/SORCE/data/overlap_vis1_mod_407nm.txt',pd_aw,pd_apath,format='(d,d)'
            ;readcol,'~/SORCE/data/overlap_vis1_mod2_uv_aligned.txt',pd_aw,pd_apath,format='(d,d)'
            ;readcol,'~/SORCE/data/overlap_vis1_smooth_mod2.txt',pd_aw,pd_apath,format='(d,d)'
            afact={wavelength:pd_aw, a:pd_apath}
        endif

        if size(kappa,/tname) ne 'STRUCT' then begin
            ; for the VIS diode, we use the same Kappa as the ESR
            ;readcol,'~/SORCE/data/ESR_kappa_453_1570_1004.txt',pdkappa_w, pdkappa_k,format='(d,d)'
            ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2001.txt',pdkappa_w, pdkappa_k,format='(d,d)'
            ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2002.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
            ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2011_reduced55.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
            ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2011_mod_v20.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
            readcol,'~/SORCE/data/ESR_kappa_453_1570_2011_mod_407nm.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
            ;readcol,'~/SORCE/data/VISAB_453-1570_kappa_2011_s5455mod2.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
            kappa={x:pdkappa_W, y:pdkappa_k}
        endif

        if size(inspectrum,/tname) ne 'STRUCT' then begin
            ; use the UNCORRECTED IRRADIANCE from version 2011 (which has the 1AU factor applied)
            ; restore SimA and SimB data
            ;restore,'~/SORCE/data/sima_vis_uncorr_2011_ccdshift.sav'
            ;restore,'~/SORCE/data/simb_vis_uncorr_2011_ccdshift.sav'
            ;inspectrum = {simA:visa_uncorr_2011_ccdshift.spect20, simB:visb_uncorr_2011_ccdshift.spect20}
            ;restore,'~/SORCE/data/sima_vis_uncorr_22_4.sav'
            ;restore,'~/SORCE/data/simb_vis_uncorr_22_4.sav'
            ;inspectrum = {simA:VISA_UNCORR_22_4.spect20, simB:VISB_UNCORR_22_4.spect20}
            restore,'~/SORCE/data/sima_vis_uncorr_2381.sav'
            restore,'~/SORCE/data/simb_vis_uncorr_2380.sav'
            inspectrum = {simA:VISA_UNCORR_2381.spect20, simB:VISB_UNCORR_2380.spect20}
        endif
        instrumentModeId=41
        yrange=[310d,950d]
        no_detdeg=0
        epsilon=10d
        order=3
        sdrange=[min([t0,t1],max=m),m]
    endif else if keyword_set(uvdiode) then begin
        if n_elements(solar54) eq 0 or n_elements(solar55) eq 0 then begin
            ;restore,'~/SORCE/data/solarexp_uv_54_55_adjusted_sol260nm.sav'
            ;restore,'~/SORCE/data/solarexp_uv_54_55_mod2.sav'
            restore,'~/SORCE/data/solarexp_uv_54_55_mod3_5100.sav'
        endif

        if size(afact,/tname) ne 'STRUCT' then begin
            ; get the uv diode raypath from the raytrace data
            ;readcol,'~/SORCE/data/overlap_uv_smooth_pos.txt',pd_aw,pd_apath,format='(d,d)'
            ;readcol,'~/SORCE/data/overlap_uv_smooth_pos_15.txt',pd_aw,pd_apath,format='(d,d)'
            ;readcol,'~/SORCE/data/overlap_uv_smooth_pos_07.txt',pd_aw,pd_apath,format='(d,d)'
            readcol,'~/SORCE/data/overlap_uv_smooth_pos_160.txt',pd_aw,pd_apath,format='(d,d)'
            afact={wavelength:pd_aw, a:pd_apath}
        endif

        if size(kappa,/tname) ne 'STRUCT' then begin
            ; for the UV diode, we use a different Kappa from the ESR
            ;readcol,'~/SORCE/data/UVAB_453-1570_kappa_2011.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
            ;readcol,'~/SORCE/data/UVAB_453-1570_kappa_2011_15.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
            ;readcol,'~/SORCE/data/UVAB_453-1570_kappa_2011_mod_15.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
            ;readcol,'~/SORCE/data/UVAB_453-1570_kappa_2011_07.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
            ;readcol,'~/SORCE/data/UVAB_453-1570_kappa_2011_160.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
            ;readcol,'~/SORCE/data/UVAB_453-1200_kappa_2011_afact160_mod407.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
            ;readcol,'~/SORCE/data/UVAB_453-1200_kappa_2011_afact160.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
            ;readcol,'~/SORCE/data/UVAB_453-1570_kappa_2011_160_s5455mod2.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
            ;readcol,'~/SORCE/data/UVAB_453-1570_kappa_2011_160_s5455mod3.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
            readcol,'~/SORCE/data/UVAB_453-1570_kappa_2011_160_s5455mod4.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
            ;readcol,'~/SORCE/data/UVAB_453-1570_kappa_2011_160_s5455_sol260nm.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
            kappa={x:pdkappa_W, y:pdkappa_k}
        endif

        if size(inspectrum,/tname) ne 'STRUCT' then begin
            ; use the UNCORRECTED IRRADIANCE from version 2011 (which has the 1AU factor applied)
            ; restore SimA and SimB data
            ;restore,'~/SORCE/data/sima_uv_uncorr_2011.sav'
            ;restore,'~/SORCE/data/simb_uv_uncorr_2011.sav'
            ;inspectrum = {simA:uva_uncorr_2011.spect20, simB:uvb_uncorr_2011.spect20}
            ;restore,'~/SORCE/data/sima_uv_uncorr_2011_ccdshift.sav'
            ;restore,'~/SORCE/data/simb_uv_uncorr_2011_ccdshift.sav'
            ;restore,'~/SORCE/data/sima_uv_uncorr_22_1.sav'
            ;restore,'~/SORCE/data/simb_uv_uncorr_22_1.sav'
            ;inspectrum = {simA:UVA_UNCORR_22_1.SPECT20, simB:UVB_UNCORR_22_1.SPECT20}
            restore,'~/SORCE/data/sima_uv_uncorr_2375.sav'
            restore,'~/SORCE/data/simb_uv_uncorr_2375.sav'
            inspectrum = {simA:UVA_UNCORR_2375.SPECT20, simB:UVB_UNCORR_2375.SPECT20}
        endif
        instrumentModeId=43
        yrange=[200d,310d]
        no_detdeg=1
        order=3
        epsilon=0.7d
        ;sdrange=[max([t0,450d]), max([t0,t1])]
        sdrange=[max([t0,100d]), max([t0,t1])]
        if n_elements(waverange) eq 0 then waverange=[200.0d,310d]
    endif

    if n_elements(waverange) eq 0 then begin
        na =LONG(n_elements(inspectrum.sima) / 2.0)
        nb =LONG(n_elements(inspectrum.simb) / 2.0)
        if ptr_valid(inspectrum.sima[0]) then begin
            mna=min((*inspectrum.sima[na]).wavelength,max=mxa)
            mnb=min((*inspectrum.simb[nb]).wavelength,max=mxb)
            waverange=[max([mna,mnb]),min([mxa,mxb])]
        endif else begin
            pa=where(inspectrum.sima[na].wavelength gt 0d)
            mna=min(inspectrum.sima[na].wavelength[pa],max=mxa)
            pb=where(inspectrum.simb[nb].wavelength gt 0d)
            mnb=min(inspectrum.simb[nb].wavelength[pb],max=mxb)
            waverange=[max([mna,mnb]),min([mxa,mxb])]
        endelse
    endif

    if ptr_valid(inspectrum.sima[0]) then begin
        sima_timestamp=dblarr(n_elements(inspectrum.sima))
        for i=0L,n_elements(sima_timestamp)-1L do sima_timestamp[i]=(*inspectrum.sima[i]).timestamp[0]

        simb_timestamp=dblarr(n_elements(inspectrum.simb))
        for i=0L,n_elements(simb_timestamp)-1L do simb_timestamp[i]=(*inspectrum.simb[i]).timestamp[0]

        k1=where(sima_timestamp gt sd2gps(453.67)*1d6)
        p=where((*inspectrum.sima[k1[0]]).wavelength ge waverange[0] and (*inspectrum.sima[k1[0]]).wavelength le waverange[1], count)
        if count eq 0 then begin
            print,'Error: no wavelengths within specified wavelength range'
            return,-1
        endif
        waves = (*inspectrum.sima[k1[0]]).wavelength[p]
    endif else begin
        k0=n_elements(inspectrum.sima[0].wavelength) / 2
        sima_timestamp=inspectrum.sima[*].timestamp[k0]
        simb_timestamp=inspectrum.simb[*].timestamp[k0]
        k1=where(inspectrum.sima[*].timestamp[k0] gt sd2gps(453.67)*1d6)
        p=where(inspectrum.sima[k1[0]].wavelength ge waverange[0] and inspectrum.sima[k1[0]].wavelength le waverange[1], count)
        if count eq 0 then begin
            print,'Error: no wavelengths within specified wavelength range'
            return,-1
        endif
        waves = inspectrum.sima[k1[0]].wavelength[p]
    endelse

    ; process every other wavelength
    if n_elements(waverange) eq 2 then begin
        p=where(waves ge min(waverange) and waves le max(waverange),count)
        if count eq 0 then begin
            print,'Error: no data within specified waverange'
            return,-1
        endif
        waves=waves[p]
        ;if instrumentModeId eq 41 then tmp=2d else tmp=4d
        ;if count gt 20 then begin
        ;    p=lindgen(count/tmp)*tmp+1
        ;    waves=waves[p]
        ;endif
    endif else begin
        ;p=lindgen(n_elements(waves)/4)*4+1
        ;waves=waves[p]
    endelse

    pa=where(sima_timestamp ge sd2gps(t0)*1d6 and sima_timestamp le sd2gps(t1)*1d6,count)
    if count lt 5 then begin
        print,'Error: not enough SIMA spectrum with timerange to fit'
        return,-1
    endif
    pb=where(simb_timestamp ge sd2gps(t0)*1d6 and simb_timestamp le sd2gps(t1)*1d6,count)
    if count lt 5 then begin
        print,'Error: not enough SIMB spectrum with timerange to fit'
        return,-1
    endif
    inspectrum={sima:inspectrum.sima[pa], simb:inspectrum.simb[pb]}

    nspctb = n_elements(inspectrum.simb)
    nspcta = n_elements(inspectrum.sima)

    ; loop over every wavelengths and find a tau for that wavelength (this will take a while)
    print,' processing '+strtrim(string(n_elements(waves)),2)+' wavelengths ...'
    nwaves=n_elements(waves)
    outdata=ptrarr(nwaves)
    for i=0L,nwaves-1 do begin 
        print,'Processing wavelength: '+strtrim(string(waves[i],format='(F0.2)'),2),format='("     ",a,"   ",$,%"\r")' 
        ;if waves[i] lt wavezero[0] then begin
            res=tnmin_pdtau2(waves[i], inspectrum=inspectrum, coeffs=coeffs, solar54=solar54, solar55=solar55, order=order, $
                step=stepsize, afact=afact, kappa=kappa, /noplot, no_detdeg=no_detdeg, sdrange=sdrange, gettau=gettau)
            if size(res,/tname) ne 'STRUCT' then continue
            outdata[i]=ptr_new({wavelength:res.wavelength, times:res.fdeg_sd, fdeg:res.fdeg, coeffs:coeffs, $
                sd54:res.sd54, degcol54:res.degcol54, sd55:res.sd55, degcol55:res.degcol55})
        ;endif else begin
        ;    npts = n_elements((*outdata[i-1]).times)
        ;    coeffs=(*outdata[i-1]).coeffs * 0d
        ;    coeffs[0]=1d
        ;    outdata[i]=ptr_new({wavelength:waves[i], times:(*outdata[i-1]).times, $
        ;        fdeg:replicate(1d,npts), coeffs:coeffs})
        ;endelse

        if i mod 4 eq 0 and NOT keyword_set(noplot) then begin
            sd=findgen(res.fdeg_sd[-1])
            plot_multi, res.fdeg_sd, res.fdeg, sd,poly(sd,coeffs),/xst,/yst,xtitle='Mission Days',$
                psym=[-4,-3], ytitle='Tau (Kappa*F_Func)', title='VisA & VisB @ '+strtrim(string(waves[i],format='(f8.2)'),2),$
                thick=[1.0,2.0];,yrange=[0.0,2.0]
        endif
    endfor

    p=where(ptr_valid(outdata) eq 1,count)
    if count eq 0 then begin
        print,'No valid fit found'
        return,-1
    endif else if count ne n_elements(outdata) then begin
        outdata=outdata[p]
    endif

    if keyword_set(surfacefit) or keyword_set(getdegcol) then begin
        ;create the arrays from the coefficients obtained above
        wv=[]
        times=[]
        degfit=[]
        degall=[]
        for i=0,n_elements(outdata)-1 do begin
            times=[times,(*outdata[i]).times]
            wv=[wv,replicate((*outdata[i]).wavelength[0], n_elements((*outdata[i]).times))]
            degfit=[degfit,poly((*outdata[i]).times, (*outdata[i]).coeffs)]
            degall=[degall,(*outdata[i]).fdeg]
        endfor
        p=where(finite(degall) eq 1, comp=cp)
        if cp[0] ne -1 then begin
            times=times[p]
            wv=wv[p]
            degfit=degfit[p]
            degall=degall[p]
        endif

        title='FFunc Degradation Mode='+strtrim(string(instrumentModeId),2)
        ; remove points too close for griddata (otherwise surface throws errors)
        grid_input,times,wv,degfit,times,wv,degfit,epsilon=epsilon

        ; get the polynomial fit
        data=dblarr(3,n_elements(times))
        ; change the axis to mission day from microseconds since sfit doesn't handle very large numbers well
        ;
        ; WE'LL NEED TO MODIFY THE COEFFICIENTS TO WORK WITH MICROSECONDS IN JAVA CODE !!!
        ;
        data[0,*]=times & data[1,*]=wv & data[2,*]=degfit
        res=sfit(data, 3, /irregular, kx=sfit_coeffs)

        ; force the FFUNC to 1.0 for wavelength above: 700nm
        k=where(data[1,*] ge wavezero[0],count)
        if count gt 0 then data[2,k]=1d

        ; do a 3-sigma clean up
        sigma3=stddev(degfit-res) * 2d
        k=where(abs(degfit-res) gt sigma3,count,complement=cp)
        if count gt 0 then res=sfit(data[*,cp], 3, /irregular, kx=sfit_coeffs)

        if keyword_set(surfacefit) then begin
            srf=surface(data[2,cp],data[0,cp],data[1,cp],irregular=grid,xtitle='Mission Day',ytitle='Wavelength',$
               ztitle='F-Function',title=title,xrange=[0,5200],yrange=yrange)
            if keyword_set(uvdiode) then begin
                ; generate 3 2D arrays with the fitted data (much cleaner)
                nwaves=n_elements(waves)/2
                waves_sub=waves[indgen(nwaves)*2]
                ndays=2600d
                days=dindgen(ndays)*2d
                wvgrid=dblarr(ndays,nwaves) 
                daygrid=dblarr(ndays,nwaves) 
                for i=0L,ndays-1 do wvgrid[i,*]=waves_sub
                for i=0,nwaves-1 do daygrid[*,i]=days
                ffgrid=sfit_poly(days,waves_sub,sfit_coeffs)
                srf=surface(ffgrid,daygrid,wvgrid,irregular=0,/overplot,color=[120,120,200])

                ct = COLORTABLE(72, /reverse)
                cntr=contour(ffgrid,days,waves_sub,/fill,rgb_table=ct,c_value=cval,xrange=[min(days,max=mx),mx],$
                    yrange=yrange,xtitle='Mission Day',ytitle='Wavelength (nm)',title=title)
                cb=colorbar()      
            endif else begin
                nwaves=n_elements(waves)/2
                waves_sub=waves[indgen(nwaves)*2]
                ndays=2600d
                days=dindgen(ndays)*2d
                wvgrid=dblarr(ndays,nwaves) 
                daygrid=dblarr(ndays,nwaves) 
                for i=0L,ndays-1 do wvgrid[i,*]=waves_sub
                for i=0,nwaves-1 do daygrid[*,i]=days
                ffgrid=sfit_poly(days,waves_sub,sfit_coeffs)
                srf=surface(ffgrid,daygrid,wvgrid,irregular=0,/overplot,color=[120,120,200])
                ;srf=surface(res,times[cp],wv[cp],irregular=0,/overplot,color=[120,120,200])
            endelse
        endif

        if keyword_set(getdegcol) then begin
            ; we'll generate a surface of degCol vs wavelength and mission day (from the solar exposure)
            nwaves=n_elements(waves)
            tmp=dblarr(nwaves)
            ; number of orbits with solarexposure > 0
            norbits=n_elements((*outdata[0]).sd54)
            ;if keyword_set(uvdiode) then begin
            ;    ; use the coefficients from sfit for a smooth function
            ;    ; NO negative values allowed (due to 3rd order polyfit extrapolation)
            ;    ffgrid=sfit_poly((*outdata[0]).sd54,waves,sfit_coeffs) > 0d   
            ;    p=where(solar54.solar_exp_orbit gt 0d)
            ;    for i=0L,n_elements(waves)-1L do ffgrid[*,i] = total(solar54[p].solar_exp_orbit * ffgrid[*,i],/cum)/86400d
            ;    degcol_c54 = {wavelength:waves, t1:sd2gps((*outdata[0]).sd54)*1d6, degcol:ffgrid, sfit_coeffs:sfit_coeffs} 
            ;endif else begin
                order=3
                if keyword_set(uvdiode) then $
                   degcol54=dblarr(norbits,nwaves) $
                else degcol54=dblarr(norbits,order+1)
                for i=0L,norbits-1 do begin
                    tmp *= 0d
                    for j=0L,nwaves-1L do tmp[j] = (*outdata[j]).degcol54[i]
                    if keyword_set(uvdiode) then $
                       degcol54[i,*]=gauss_smooth(tmp, 2.0, /edge_mirror) $
                    else degcol54[i,0:order] = robust_poly_fit(waves,tmp,order)
                endfor
                if keyword_set(uvdiode) then $
                    degcol_c54 = {wavelength:waves, t1:sd2gps((*outdata[0]).sd54)*1d6, degcol:degcol54, sfit_coeffs:sfit_coeffs}  $
                else $
                    degcol_c54 = {wavelength:waves, t1:sd2gps((*outdata[0]).sd54)*1d6, coeffs:degcol54, sfit_coeffs:sfit_coeffs}
            ;endelse

            norbits=n_elements((*outdata[0]).sd55)
            ;if keyword_set(uvdiode) then begin
            ;    ; use the coefficients from sfit for a smooth function
            ;    ; NO negative values allowed (due to 3rd order polyfit extrapolation)
            ;    ffgrid=sfit_poly((*outdata[0]).sd55,waves,sfit_coeffs) > 0d   
            ;    p=where(solar55.solar_exp_orbit gt 0d)
            ;    for i=0L,n_elements(waves)-1L do ffgrid[*,i] = total(solar55[p].solar_exp_orbit * ffgrid[*,i],/cum)/86400d
            ;    degcol_c55 = {wavelength:waves, t1:sd2gps((*outdata[0]).sd55)*1d6, degcol:ffgrid, sfit_coeffs:sfit_coeffs} 
            ;endif else begin
                if keyword_set(uvdiode) then $
                   degcol55=dblarr(norbits,nwaves) $
                else degcol55=dblarr(norbits,order+1)
                for i=0L,norbits-1 do begin
                    tmp *= 0d
                    for j=0L,nwaves-1L do tmp[j] = (*outdata[j]).degcol55[i]
                    if keyword_set(uvdiode) then $
                        degcol55[i,*]=gauss_smooth(tmp, 2.0, /edge_mirror) $
                    else $
                        degcol55[i,0:order] = robust_poly_fit(waves,tmp,order)
                endfor
               if keyword_set(uvdiode) then $
                    degcol_c55 = {wavelength:waves, t1:sd2gps((*outdata[0]).sd55)*1d6, degcol:degcol55, sfit_coeffs:sfit_coeffs}  $
                else $
                    degcol_c55 = {wavelength:waves, t1:sd2gps((*outdata[0]).sd55)*1d6, coeffs:degcol55, sfit_coeffs:sfit_coeffs}
            ;endelse


            ;nwaves=n_elements(waves)/2
            ;waves=waves[indgen(nwaves)*2]
            ;p54=where(solar54.solar_exp_orbit gt 0d,c54)
            ;p55=where(solar55.solar_exp_orbit gt 0d,c55)
            ;sd54=gps2sd(solar54[p54].t1/1d6)
            ;sd55=gps2sd(solar55[p55].t1/1d6)
            ;degcol54 = dblarr(c54,nwaves)
            ;degcol55 = dblarr(c55,nwaves)
            ;for i=0L,nwaves-1 do degcol54[*,i]=total(solar54[p54].solar_exp_orbit * sfit_poly(sd54,waves[i],sfit_coeffs),/cum)
            ;degcol54 /= 86400d
            ;for i=0L,nwaves-1 do degcol55[*,i]=total(solar55[p55].solar_exp_orbit * sfit_poly(sd55,waves[i],sfit_coeffs),/cum)
            ;degcol55 /= 86400d
            ;; for each day, get the 3rd order polynomial fit of degcol vs wavelength
            ;order=3
            ;temp54=dblarr(c54,order+1)
            ;for i=0L,c54-1 do temp54[i,*] = poly_fit(waves,degcol54[i,*],3)
            ;temp55=dblarr(c55,order+1)
            ;for i=0L,c55-1 do temp55[i,*] = poly_fit(waves,degcol55[i,*],3)
            ;degcol_c54 = {wavelength:waves, t1:solar54[p54].t1, coeffs:temp54, sfit_coeffs:sfit_coeffs}
            ;degcol_c55 = {wavelength:waves, t1:solar55[p55].t1, coeffs:temp55, sfit_coeffs:sfit_coeffs}

        endif

        coeffs=sfit_coeffs
    endif

    ;if modes[0] eq 41 then kappa_fit=[pfit[k],pfit_1[ck]]
    return, outdata

END

