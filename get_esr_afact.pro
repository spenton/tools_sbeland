;+
; NAME:   GET_ESR_AFACT
;
; PURPOSE: 
;    This routine various blocks of ESR data from SimA and SimB to
;    calculate a Kappa for each block and then adjusting the raypath
;    so that the Kappa from each block overlaps.
;    It allows us to estimate a new value of the ESR raypath instead 
;    of relying on the raytrace values.
;
; CALLING SEQUENCE:
;    result=get_esr_afact(inspectrum=inspectrum, solar54=solar54, solar55=solar55)
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
;   Revision: $Id: get_esr_afact.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;-
;*******************************************************************
function min_get_esr_afact, coeffs

    ;common esr_afact, times, solar54, solar55, esrspect, afactw, niter

    new_afact={wavelength:afactw, a:poly(afactw,coeffs)}

    res0=get_esrkappa(solar54=solar54,solar55=solar55,afact=new_afact,inspect=esrspect,timerange=[times[0],times[1]],/noplot)
    res1=get_esrkappa(solar54=solar54,solar55=solar55,afact=new_afact,inspect=esrspect,timerange=[times[1],times[2]],/noplot)
    res2=get_esrkappa(solar54=solar54,solar55=solar55,afact=new_afact,inspect=esrspect,timerange=[times[2],times[3]],/noplot)

    ; use the BSpline since it fits the data much better but limit the wavelength range
    p=where(res0.wave_fit ge 290 and res0.wave_fit le 1000)

    delta = [res1.kappa_fit[p]-res0.kappa_fit[p], res1.kappa_fit[p]-res2.kappa_fit[p]]
    out_value=stddev(delta)
    format = '("   ",i0,",  ",G15.6,"  (",'+strtrim(string(n_elements(coeffs)-1),2)+'(E,", "),E," )",$,%"\r")'
    print,niter, out_value, coeffs, format=format
    plot_multi,new_afact.wavelength,new_afact.a,/xst,/yst,title='AFACT='+string(coeffs,format='('+$
        string(n_elements(coeffs))+'F)'),label=string(out_value)
    niter+=1L
    return, out_value
end
;*******************************************************************
function min_get_esr_afact2, coeffs, df, times=times, solar54=solar54, solar55=solar55,  $
    esrspect=esrspect, afactw=afactw, niter=niter, verbose=verbose

    new_afact={wavelength:afactw, a:poly(afactw,coeffs)}

    format = '("   ",i0,",  (",'+strtrim(string(n_elements(coeffs)-1),2)+'(E,", "),E," )",$)' ;,%"\r")'
    print,niter, coeffs, format=format

    res0=get_esrkappa(solar54=solar54,solar55=solar55,afact=new_afact,inspect=esrspect,timerange=[times[0],times[1]],/noplot,verbose=verbose)
    res1=get_esrkappa(solar54=solar54,solar55=solar55,afact=new_afact,inspect=esrspect,timerange=[times[1],times[2]],/noplot,verbose=verbose)
    res2=get_esrkappa(solar54=solar54,solar55=solar55,afact=new_afact,inspect=esrspect,timerange=[times[2],times[3]],/noplot,verbose=verbose)

    ; use the BSpline since it fits the data much better but limit the wavelength range
    p=where(res0.wave_fit ge 290 and res0.wave_fit le 1000)

    delta = [res1.kappa_fit[p]-res0.kappa_fit[p], res1.kappa_fit[p]-res2.kappa_fit[p]]
    ;out_value=stddev(delta)
    out_value=median(delta)
    format = '("   ",E15.8)' ;,$,%"\r")'
    print,out_value, format=format
    ;plot_multi,new_afact.wavelength,new_afact.a,/xst,/yst,title='AFACT='+string(coeffs,format='('+$
    ;    string(n_elements(coeffs))+'F)'),label=string(out_value)
    niter+=1
    return, out_value
end
;*****************************************************************************
FUNCTION GET_ESR_AFACT, inspectrum=inspectrum, solar54=sol54, solar55=sol55, noplot=noplot, $
    ffunct=ffunct, afact=afact, verbose=verbose, version=version, order=order, maxiter=maxiter

    ;common esr_afact

    ; version 1003 has the 1AU applied to SimUncorrectedIrradiance already
    if n_elements(version) eq 0 then version=1003
    if n_elements(user) eq 0 then     user='sbeland'
    if n_elements(password) eq 0 then password='sbeland1'
    if n_elements(dburl) eq 0 then    dbUrl='jdbc:sybase:Tds:lasp-db-dev:4100/SORCE_SIM_V20'
    if n_elements(dbdriver) eq 0 then dbDriver='com.sybase.jdbc3.jdbc.SybDriver'
    if n_elements(verbose) eq 0 then verbose=0

    times=[453d,825d, 1200d, 1570d]

    if size(inspectrum,/tname) ne 'STRUCT' then begin
        ; get all of the ESR data within the specified timerange and only keep around the requested wavelength
        if keyword_set(verbose) then print,'   getting esrA ...'
        esrA=get_science_product(['Wavelength','SimUncorrectedIrradiance','SimConvertedDataNumbers'],times[0],times[-1],31,/mission,$
            version=version,user=user,password=password,dburl=dburl,dbdriver=dbdriver)
        if keyword_set(verbose) then print,'   getting esrB ...'
        esrB=get_science_product(['Wavelength','SimUncorrectedIrradiance','SimConvertedDataNumbers'],times[0],times[-1],32,/mission,$
            version=version,user=user,password=password,dburl=dburl,dbdriver=dbdriver)
        esrspect={sp31:esrA, sp32:esrB}
    endif else esrspect=inspectrum

    if n_elements(sol54) eq 0 then begin
        if keyword_set(verbose) then print,'   getting solar54 ...'
        q1 = 'SELECT orbitStartTimeGps t0, orbitEndTimeGps t1, cumulativeSolarExpHrtOut solar_exp, '+$
             'solarExposureHrtOut as solar_exp_orbit from SimSolarExposureData where instrumentModeId=54 and '+ $
             'orbitStartTimeGps>='+strtrim(string(ulong64(sd2gps(times[0]-1)*1d6)),2)+' and '+ $
             'orbitStartTimeGps<='+strtrim(string(ulong64(sd2gps(times[-1]+1)*1d6)),2)
        query_database, q1, sol54, info
    endif 
    solar54=sol54

    if n_elements(sol55) eq 0 then begin
        if keyword_set(verbose) then print,'   getting solar55 ...'
        q1 = 'SELECT orbitStartTimeGps t0, orbitEndTimeGps t1, cumulativeSolarExpHrtOut solar_exp, '+$
             'solarExposureHrtOut as solar_exp_orbit from SimSolarExposureData where instrumentModeId=55 and '+ $
             'orbitStartTimeGps>='+strtrim(string(ulong64(sd2gps(times[0]-1)*1d6)),2)+' and '+ $
             'orbitStartTimeGps<='+strtrim(string(ulong64(sd2gps(times[-1]+1)*1d6)),2)
        query_database, q1, sol55, info
    endif 
    solar55=sol55
 
    ; read the raytrace raypath as a starting point
    if n_elements(order) eq 0 then order=8
    if n_elements(afact) eq 0 then begin
        readcol,'~/SORCE/data/overlap_esr_fit.txt',esrw,esrray1,esrray2,format='(d,d,d)'
        coeff=poly_fit(esrw,esrray2,order)
        afact={wavelength:esrw, a:poly(esrw,coeff)}
    endif else begin
        coeff=poly_fit(afact.wavelength,afact.a,order)
    endelse

    afactw= afact.wavelength
    ;init_coeff=reform(coeff)
    ;scale=init_coeff

    niter=0L
    ;coeffs = AMOEBA(1.0d-6, P0=init_coeff, scale=scale, FUNCTION_VALUE=fval, FUNCTION_NAME='min_get_esr_afact', nmax=maxiter)

    functargs = {times:times, solar54:solar54, solar55:solar55, esrspect:esrspect, afactw:afactw, niter:niter, verbose:verbose}
    parinfo=replicate({value:0.0d, fixed:0, limited:[0,0], limits:[0d,0d], step:0d, tnside:2},order+1)
    parinfo.value=reform(coeff)
    ;parinfo[0].step=1d-4
    ;parinfo[1].step=1d-6
    ;parinfo[1].limited=[1,0]
    ;parinfo[1].limits=[0d,0.01d]

    coeffs = tnmin('min_get_esr_afact2', functargs=functargs, bestmin=f0, status=status, epsrel=1d-4, $
        nfev=nfev, autoderivative=1, errmsg=errmsg, parinfo=parinfo, maxiter=maxiter, /quiet)

    print,''
    print,'coeffs=',coeffs, format='(G)'
    if n_elements(coeffs) eq 1 then begin
        print,'AMOEBA failed to converge'
        return,-1
    endif

    ; now display the best fit
    new_afact={wavelength:afact.wavelength, a:poly(afact.wavelength,coeffs)}
    res0=get_esrkappa(solar54=solar54,solar55=solar55,afact=new_afact,inspect=esrspect,timerange=[times[0],times[-1]],/noplot,verbose=verbose)
    res1=get_esrkappa(solar54=solar54,solar55=solar55,afact=new_afact,inspect=esrspect,timerange=[times[0],times[1]],/noplot,verbose=verbose)
    res2=get_esrkappa(solar54=solar54,solar55=solar55,afact=new_afact,inspect=esrspect,timerange=[times[1],times[2]],/noplot,verbose=verbose)
    res3=get_esrkappa(solar54=solar54,solar55=solar55,afact=new_afact,inspect=esrspect,timerange=[times[2],times[3]],/noplot,verbose=verbose)

    if NOT keyword_set(noplot) then begin
        yfit0 = (res0.exp_coeff[0]*exp(res0.exp_coeff[1]*res0.wave_fit)+res0.exp_coeff[2])
        yfit1 = (res1.exp_coeff[0]*exp(res1.exp_coeff[1]*res1.wave_fit)+res1.exp_coeff[2])
        yfit2 = (res2.exp_coeff[0]*exp(res2.exp_coeff[1]*res2.wave_fit)+res2.exp_coeff[2])
        yfit3 = (res3.exp_coeff[0]*exp(res3.exp_coeff[1]*res3.wave_fit)+res3.exp_coeff[2])
        plot_multi,res0.wavelength,res0.kappa,res0.wave_fit,yfit0,res1.wavelength,res1.kappa,res1.wave_fit,yfit1,$
            res2.wavelength,res2.kappa,res2.wave_fit,yfit2,res3.wavelength,res3.kappa,res3.wave_fit,yfit3, $
            label=['Kappa 453->1570','ExpFit','Kappa 453->825','ExpFit','Kappa 825->1200','ExpFit','Kappa 1200->1570','ExpFit'],$
            xtitle='Wavelength (nm)', ytitle='ESR Kappa',charsize=1.4, $
            /xst,/yst,psym=[-4,-3,-4,-3,-4,-3,-4,-3], color=[0,0,2,2,3,3,4,4]

        plot_multi,new_afact.wavelength,new_afact.a, /xst,/yst,xtitle='Wavelength',ytitle='ESR Raypath',$
            title='ESR Raypath from 453->825->1200->1570',charsize=1.4

    endif


    return, new_afact

END

