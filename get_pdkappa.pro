;+
; NAME:   GET_PDKAPPA
;
; PURPOSE: 
;    This routine uses the photodiode spectrum from SIMA and SIMB for the
;    first part of the mission to determine a Kappa covering the specific
;    spectral region (to compare with the kappa obtained from the ESR).
;
; CALLING SEQUENCE:
;    result=get_pdkappa(instrumetModeId, inspectrum=inspectrum, solar54=solar54, solar55=solar55)
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
;   Revision: $Id: get_pdkappa.pro,v 1.2 2018/11/26 16:16:02 spenton Exp $
;-
;*****************************************************************************
FUNCTION GET_PDKAPPA, instrumentModeId, inspectrum=inspectpd, solar54=solar54, solar55=solar55, noplot=noplot, $
    ffunct=ffunct, afact=afact, t0=t0, t1=t1, waverange=waverange, version=version, bversion=bversion,irrad_table=irrad_table

    modes = [41,43,44]
    if n_elements(t0) eq 0 then t0=453d
    if n_elements(t1) eq 0 then t1=1570d
    ;t0=1570d
    ;t1=3200d

    match,modes,instrumentModeId,suba,subb,count=count
    if count eq 0 then begin
        print,'Error: wrong instrumentModeId (expecting '+modes,+')'
        return,-1
    endif
    ; define the modes for SimA and SimB
    modes=[modes[suba],modes[suba]+4]

    ; version 16 of development database contains SimCalibratedIrradiance without degradation
    ; and fixes for pd tempCorr, pre 453 ccdshift and mode48 dn adjustments
    ;if n_elements(version) eq 0 then version=16
    if n_elements(version) eq 0 then version=2011

    if n_elements(irrad_table) eq 0 then irrad_table='SimUncorrectedIrradiance'

    ; get the raypath
    if n_elements(afact) eq 0 then begin
        query_database, /reset
        q2 = "SELECT calibrationSetId from CalibrationMetadata where calibrationTableName='SimRayPathDegradationParams' "+$
            "and version=6 and instrumentModeId="+strtrim(string(instrumentModeId),2)
        query_database, q2, res, info
        q3='select * from SimRayPathDegradationParams where calibrationSetId='+strtrim(string(res.(0)),2)
        query_database, q3, raypath, info
        afact={wavelength:raypath.wavelength, a:raypath.SINGLEPASSAREAFRACTION}
    endif
        
    ; get the spectrum and solar exposure data for a single wavelength 
    ; and use the first spectra to determine the wavelengths array
    if size(inspectpd,/tname) eq 'STRUCT' then begin
        ; get the wavelengths from the spectra in the middle of the array
        pos=n_elements(inspectpd.sima.diode)/2
        waves = (*inspectpd.sima.diode[pos]).wavelength
    endif else begin
        res=tnmin_pddeg(t0,t1,modes,400d,coeffs=c,inspectrum=inspectpd, solar54=solar54, solar55=solar55, version=version, $
            bversion=bversion,irrad_table=irrad_table)
        waves = (*inspectpd.sima.diode[0]).wavelength
    endelse

    ; process every other wavelength
    if n_elements(waverange) eq 2 then begin
        p=where(waves ge min(waverange) and waves le max(waverange),count)
        if count eq 0 then begin
            print,'Error: no data within specified waverange'
            return,-1
        endif
        waves=waves[p]
        if count gt 20 then begin
            p=lindgen(count/2)*2+1
            waves=waves[p]
        endif
    endif else begin
        p=lindgen(n_elements(waves)/2)*2+1
        waves=waves[p]
    endelse

    pdkappa=dblarr(n_elements(waves))
    ; loop over every wavelengths and find a kappa for that wavelength (this will take a while)
    print,' processing '+strtrim(string(n_elements(waves)),2)+' wavelengths ...'
    nwaves=n_elements(waves)
    ;nskip = fix(nwaves) / 20   ; limit the number of plots to 30 by skipping the right amount of wavelengths
    nskip=1
    title0=strtrim(string(modes[0]),2)+", "+strtrim(string(modes[0]),2)
    if instrumentModeId eq 41 then ymax=2.2 else ymax=1.0
    window,0
    for i=0,nwaves-1 do begin 
        res=tnmin_pddeg(t0,t1,modes,waves[i],coeffs=coeff0,inspectrum=inspectpd, ffunct=ffunct, $
            solar54=solar54, solar55=solar55, afact=afact, version=version, bversion=bversion,irrad_table=irrad_table)
        pdkappa[i]=coeff0 
        print,i+1,nwaves,waves[i],format='($,I0," of ",I0,"  (",F0.2,")   ")'
        if (i MOD nskip) eq 0 and NOT keyword_set(noplot) then begin
            ; plot the newly processed A & B data to verify alignment
            ; align the 2 timeseries in irradiance
            meda=median(res.sima.new_irrad)
            medb=median(res.simb.new_irrad)
            delta=meda - medb
            yval = [res.sima.new_irrad, res.simb.new_irrad+delta]
            p=where(yval gt 0.d and yval lt ymax,count)
            if count gt 1 then yrange=[min(yval[p],max=m),m] else yrange=[0d,1d]
            title=title0+" Kappa="+strtrim(string(pdkappa[i],format='(d0.6)'),2)+" @ "+strtrim(string(waves[i],format='(d0.2)'),2)+" nm"
            ;plot_multi,gps2sd(res.sima.TIMESTAMP/1d6), res.sima.new_irrad, $
            ;    gps2sd(res.simb.TIMESTAMP/1d6), res.simb.new_irrad+delta,psym=[-4,-4], yrange=yrange, $
            ;    xtitle='Mission Day',ytitle='Corrected Irradiance',title=title
            plot,gps2sd(res.sima.TIMESTAMP/1d6), res.sima.new_irrad,yrange=yrange,title=title,/nodata,color=0
            oplot,gps2sd(res.sima.TIMESTAMP/1d6), res.sima.new_irrad,psym=-4, color=2
            oplot,gps2sd(res.simb.TIMESTAMP/1d6), res.simb.new_irrad+delta, psym=-4, color=3
        endif
    endfor
    print,''

    ; cleanup the kappa
    pdkappa=abs(pdkappa)

    ; clean up the data by instrumentModeId
    if modes[0] eq 41 then begin
        ww=dindgen(665*5+1)/5d +300d
        guess=[0.01d, -0.007d, 1.0d-5]
        p=where(waves gt 327d and waves lt 800d)
        dd=exponential_fit(waves[p],pdkappa[p],guess=guess, fita=[1,1,1])
        exp_fit=dd[0]*exp(ww*dd[1]) + dd[2]
        ; exponential doesn't fit very well but a 4th order polynomial does
        coeff=robust_poly_fit(waves[p],pdkappa[p],4,/double) 
        pfit=poly(ww,coeff)
        p1=where(waves gt 625)
        coeff_1=robust_poly_fit(waves[p1],pdkappa[p1],4,/double) 
        pfit_1=poly(ww,coeff_1)
        k=where(ww lt 610d,comp=ck)
        pfit=[pfit[k],pfit_1[ck]]
        ; extend the bspline to 300nm using the 4th order polynomial fit
        w0=dindgen(24)+300
        y0=poly(w0,coeff)
        p=where(waves gt 324.5d)
        ; first remove outlyers
        coeff=robust_poly_fit(waves,pdkappa,4,/double) 
        this_wv=[w0,waves[p]]
        values=poly(this_wv,coeff)
        resistant_mean, abs(values - [y0,pdkappa[p]]), 3.0, mean,good=keep
        sset=bspline_iterfit(this_wv[keep],smooth(([y0,pdkappa[p]])[keep],3),maxiter=0,requiren=10,bkspace=5,nord=4)
        kappa_fit=bspline_valu(ww,sset)
    endif else if modes[0] eq 43 then begin
        ww=dindgen(110*20+1)/20d +200d
        ; clean the peaks
        p4=where(waves gt 270d, comp=cp4)
        sset=bspline_iterfit(waves[cp4],pdkappa[cp4],requiren=10,bkspace=5,nord=3)
        fit=bspline_valu(waves[cp4],sset)
        resistant_mean,abs(pdkappa[cp4]-fit),4.0,mean,good=k0

        sset=bspline_iterfit(waves[p4],pdkappa[p4],requiren=1.25,bkspace=1.25,nord=3)
        fit=bspline_valu(waves[p4],sset)
        resistant_mean,abs(pdkappa[p4]-fit),4.0,mean,good=k1

        new_waves=[waves[cp4[k0]], waves[p4[k1]]]
        uvkappa=[pdkappa[cp4[k0]], pdkappa[p4[k1]]]

        ;if n_elements(k0) gt 0 and n_elements(k0) lt n_elements(waves) then begin
        ;    pos=intarr(n_elements(waves))
        ;    pos[k0]=1
        ;    k1=where(pos eq 0)
        ;    uvkappa[k1]=-1
        ;endif
        ; for the UV kappa, remove the solar signature by masking out some areas
        ;p=where(waves gt 307d,count)
        ;if count gt 0 then uvkappa[p]=-1
        ;p=where(waves gt 271.5d and waves lt 290.5d, count)
        ;if count gt 0 then uvkappa[p]=-1
        ;p=where(waves gt 243.7d and waves lt 264.0d, count)
        ;if count gt 0 then uvkappa[p]=-1
        ;p=where(waves gt 233.8d and waves lt 241.8d, count) 
        ;if count gt 0 then uvkappa[p]=-1
        ;p=where(waves gt 214.5d and waves lt 218.2d, count) 
        ;if count gt 0 then uvkappa[p]=-1
        ;p=where(waves gt 218.2d and waves lt 228.0d, count) 
        ;if count gt 0 then uvkappa[p]+=8.25d-5
        ;p=where(waves gt 218.2d and waves lt 223.2d, count) 
        ;if count gt 0 then uvkappa[p]+=8.25d-5
        ;p=where(waves gt 214.5d and waves lt 223.2d, count) 
        ;if count gt 0 then uvkappa[p]=-1
        ;k=where(uvkappa gt 0d)
        ;sset=bspline_iterfit(waves[k],uvkappa[k],maxiter=0,requiren=10,bkspace=5,nord=3)
        ;kappa_fit=bspline_valu(ww,sset)

        ;new_waves=waves[k0]
        ;uvkappa=pdkappa[k0]

        ; fit tthe last 3 points at 307.4, 307.8 and 308
        p0=where(new_waves gt 305d,count0)
        p1=where(new_waves ge 299d and new_waves le 307d,count1)
        if count0 gt 0 then begin
            cc=ladfit(new_waves[p1],uvkappa[p1])
            uvkappa[p0]=poly(new_waves[p0],cc)
            new_waves=[new_waves,310d]
            uvkappa=[uvkappa,poly(new_waves[-1],cc)]
        endif
        p0=where(new_waves le 267d,comp=p1)
        p2=where(new_waves ge 292d)
        p3=where(new_waves le 310d)
        uvkappa=[gauss_smooth(uvkappa[p0],10,width=20,/edge_mirror), $
                 gauss_smooth(uvkappa[p1],2,width=3,/edge_mirror)] 
        coeff0=robust_poly_fit(new_waves[p2], uvkappa[p2], 3,/double)
        uvkappa[p2]=poly(new_waves[p2],coeff0)
        ;uvkappa[p2]=smooth(uvkappa[p2],3)
        sset0=bspline_iterfit(new_waves[p3],uvkappa[p3],requiren=0.5,bkspace=0.5,nord=3)
        kappa_fit=bspline_valu(ww,sset0)
        ;kappa_fit0=bspline_valu(ww,sset0)
        ;sset1=bspline_iterfit(new_waves[p2],uvkappa[p2],requiren=10,bkspace=5,nord=3)
        ;kappa_fit1=bspline_valu(ww,sset1)
        ;p=where(ww ge 292d,comp=k)
        ;kappa_fit=[kappa_fit0[k],kappa_fit1[p]]

    endif else if modes[0] eq 44 then begin
        ww=dindgen(860*5+1)/5d +840d
        sset=bspline_iterfit(waves,smooth(pdkappa,3),maxiter=0,requiren=10,bkspace=5,nord=4)
        kappa_fit=bspline_valu(ww,sset)
    endif


    ; plot the data
    if NOT keyword_set(noplot) then begin
        title="Degradation model Kappa (Modes "+strtrim(string(modes[0]),2)+', '+strtrim(string(modes[1]),2)+')'
        title=string(title,t0,t1,format='(A,"  SD=",I0," to ",I0)')
        if n_elements(pfit) gt 0 then begin
            fit_label='Poly Fit '
            for i=0,n_elements(coeff)-1 do fit_label=fit_label+strtrim(string(coeff[i]),2)+', '
            plot_multi,waves,pdkappa,ww,kappa_fit,ww,pfit,/xst,/yst,title=title,yrange=[0d,0.004d],$
                xtitle="Wavelength (nm)",ytitle="Kappa",charsize=1.4,label=["Kappa","BSpline Fit",fit_label],psym=[-4,-3,-3]
        endif else begin
            plot_multi,waves,pdkappa,ww,kappa_fit,/xst,/yst,title=title,yrange=[0d,0.004d],$
                xtitle="Wavelength (nm)",ytitle="Kappa",charsize=1.4,label=["Kappa","BSpline Fit"],psym=[-4,-3], $
                thick=[1.0,3.0]
        endelse
    endif

    ;if modes[0] eq 41 then kappa_fit=[pfit[k],pfit_1[ck]]
    return, {wavelength:waves, kappa:pdkappa, x:ww, y:kappa_fit}

END

