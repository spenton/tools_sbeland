;+
; NAME:   GET_PDTAU
;
; PURPOSE: 
;    This routine uses the photodiode spectrum from SIMA and SIMB for the
;    first part of the mission to determine a TAU covering the specific
;    spectral region.
;
; CALLING SEQUENCE:
;    result=get_pdtau(instrumetModeId, inspectrum=inspectrum, solar54=solar54, solar55=solar55)
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
;   Revision: $Id: get_pdtau.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;-
;*****************************************************************************
FUNCTION GET_PDTAU, t0, t1, visdiode=visdiode, uvdiode=uvdiode, stepsize=stepsize, solar54=solar54, solar55=solar55, $
    inspectrum=inspectrum, noplot=noplot, ffunct=ffunct, afact=afact, kappa=kappa, waverange=waverange

    if n_elements(stepsize) eq 0 then stepsize=90d

    if keyword_set(visdiode) then begin
        ; get the vis diode raypath from the raytrace data
        readcol,'~/SORCE/data/overlap_vis1.txt',pd_aw,pd_apath,format='(d,d)',/silent
        ;just as a test multiply the raytrace raypath by 4.8
        pd_apath *= 4.8d
        afact={wavelength:pd_aw, a:pd_apath}

        if size(kappa,/tname) ne 'STRUCT' then begin
            ; for the VIS diode, we use the same Kappa as the ESR
            ;readcol,'~/SORCE/data/ESR_kappa_453_1570_1004.txt',pdkappa_w, pdkappa_k,format='(d,d)'
            ;readcol,'~/SORCE/data/ESR_kappa_453_1570_2001.txt',pdkappa_w, pdkappa_k,format='(d,d)'
            readcol,'~/SORCE/data/ESR_kappa_453_1570_2002.txt',pdkappa_w, pdkappa_k,format='(d,d)',/silent
            kappa={x:pdkappa_W, y:pdkappa_k}
        endif

        if size(inspectrum,/tname) ne 'STRUCT' then begin
            ; use the UNCORRECTED IRRADIANCE from version 2011 (which has the 1AU factor applied)
            ; restore SimA and SimB data
            restore,'~/SORCE/data/sima_vis_uncorr_2011.sav'
            restore,'~/SORCE/data/simb_vis_uncorr_2011.sav'
            inspectrum = {simA:visa_uncorr_2011.spect20, simB:visb_uncorr_2011.spect20}
        endif
    endif

    if n_elements(waverange) eq 0 then begin
        pa=where(inspectrum.sima[0].wavelength gt 0d)
        pb=where(inspectrum.simb[0].wavelength gt 0d)
        mna=min(inspectrum.sima[0].wavelength[pa],max=mxa)
        mnb=min(inspectrum.simb[0].wavelength[pb],max=mxb)
        waverange=[max([mna,mnb]),min([mxa,mxb])]
    endif

    p=where(inspectrum.sima[0].wavelength ge waverange[0] and inspectrum.sima[0].wavelength le waverange[1], count)
    if count eq 0 then begin
        print,'Error: no wavelengths within specified wavelength range'
        return,-1
    endif
    waves = inspectrum.sima[0].wavelength[p]
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

    pa=where(inspectrum.sima[*].timestamp[0] ge sd2gps(t0)*1d6 and inspectrum.sima[*].timestamp[0] le sd2gps(t1)*1d6,count)
    if count lt 5 then begin
        print,'Error: not enough SIMA spectrum with timerange to fit'
        return,-1
    endif
    pb=where(inspectrum.simb[*].timestamp[0] ge sd2gps(t0)*1d6 and inspectrum.simb[*].timestamp[0] le sd2gps(t1)*1d6,count)
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
    alldata=ptrarr(nwaves)
    for i=0L,nwaves-1 do begin 
        allwaves=[]
        alltimes=[]
        alltau=[]
        delta=(inspectrum.simb.timestamp[0] - inspectrum.simb[0].timestamp[0])
        posB=where(delta ge 0d and delta le stepsize*86400d6, countB)
        posA=where(inspectrum.sima.timestamp[0] ge inspectrum.simb[posB[0]].timestamp[0] and $
                   inspectrum.sima.timestamp[0] le inspectrum.simb[posB[-1]].timestamp[0], countA)
        print,i+1,nwaves,waves[i],format='(I0," of ",I0,"  (",F0.2,")")'
        while posB[-1] lt nspctB-1 do begin
            if n_elements(posA) ge 3 and n_elements(posB) ge 3 then begin
                inspectpd={sima:inspectrum.sima[posA], simb:inspectrum.simb[posB]}
                res=tnmin_pdtau(waves[i], inspectpd, coeffs=coeff, ffunct=ffunct, $
                    solar54=solar54, solar55=solar55, afact=afact, kappa=kappa, /noplot)
                if size(res,/tname) eq 'STRUCT' then begin
                    allwaves=[allwaves,waves[i]]
                    alltau=[alltau,coeff[0]]
                    tmp=(inspectrum.simb[posB].timestamp)[*]
                    p=where(tmp gt 0d)
                    mn=min(tmp[p],max=mx)
                    alltimes=[alltimes,gps2sd((mn+mx)/2d6)]
                endif
            endif
            mx=max(inspectrum.simb[posB].timestamp)
            print,gps2sd(mx/1d6), format='(F15.2,$,%"\r")'
            delta=(inspectrum.simb.timestamp[0] - inspectrum.simb[posB[0]+1].timestamp[0])
            posB=where(delta ge 0d and delta le stepsize*86400d6, count)
            posA=where(inspectrum.sima.timestamp[0] ge inspectrum.simb[posB[0]].timestamp[0] and $
                   inspectrum.sima.timestamp[0] le inspectrum.simb[posB[-1]].timestamp[0], count)
        endwhile
        p=where(finite(alltau) eq 1, count)
        if count gt 5 then cc=robust_poly_fit(alltimes[p],alltau[p],4) else cc=ladfit(alltimes[p],alltau[p])
        alldata[i]=ptr_new({wavelength:allwaves, times:alltimes, tau:alltau, coeffs:cc})
        plot_multi,alltimes[p],alltau[p], alltimes[p],poly(alltimes[p],cc),/xst,/yst,xtitle='Mission Days',psym=[-4,-3],$
            ytitle='Tau (Kappa*F_Func)', title='VisA & VisB @ '+strtrim(string(waves[i],format='(f8.2)'),2),thick=[1.0,2.0]
    endfor

    ; cleanup the values of Tau
    alltau=abs(alltau)

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
        sset=bspline_iterfit([w0,waves[p]],smooth([y0,pdkappa[p]],3),maxiter=0,requiren=10,bkspace=5,nord=4)
        kappa_fit=bspline_valu(ww,sset)
    endif else if modes[0] eq 43 then begin
        ww=dindgen(110*20+1)/20d +200d
        ; clean th epeaks
        uvkappa=smooth(pdkappa,3)
        sset=bspline_iterfit(waves,uvkappa,maxiter=0,requiren=10,bkspace=5,nord=3)
        kappa_fit=bspline_valu(ww,sset)
        fit=interpol(kappa_fit,ww,waves)
        resistant_mean,uvkappa-fit,3.0,mean,good=k0
        if n_elements(k0) gt 0 and n_elements(k0) lt n_elements(waves) then begin
            pos=intarr(n_elements(waves))
            pos[k0]=1
            k1=where(pos eq 0)
            uvkappa[k1]=-1
        endif
        ; for the UV kappa, remove the solar signature by masking out some areas
        p=where(waves gt 307d,count)
        if count gt 0 then uvkappa[p]=-1
        p=where(waves gt 271.5d and waves lt 290.5d, count)
        if count gt 0 then uvkappa[p]=-1
        p=where(waves gt 243.7d and waves lt 264.0d, count)
        if count gt 0 then uvkappa[p]=-1
        p=where(waves gt 233.8d and waves lt 241.8d, count) 
        if count gt 0 then uvkappa[p]=-1
        p=where(waves gt 214.5d and waves lt 218.2d, count) 
        if count gt 0 then uvkappa[p]=-1
        p=where(waves gt 218.2d and waves lt 228.0d, count) 
        if count gt 0 then uvkappa[p]+=8.25d-5
        p=where(waves gt 218.2d and waves lt 223.2d, count) 
        if count gt 0 then uvkappa[p]+=8.25d-5
        p=where(waves gt 214.5d and waves lt 223.2d, count) 
        if count gt 0 then uvkappa[p]=-1
        k=where(uvkappa gt 0d)
        sset=bspline_iterfit(waves[k],uvkappa[k],maxiter=0,requiren=10,bkspace=5,nord=3)
        kappa_fit=bspline_valu(ww,sset)
    endif else if modes[0] eq 44 then begin
        ww=dindgen(860*5+1)/5d +840d
        sset=bspline_iterfit(waves,smooth(pdkappa,3),maxiter=0,requiren=10,bkspace=5,nord=4)
        kappa_fit=bspline_valu(ww,sset)
    endif


    ; plot the data
    if NOT keyword_set(noplot) then begin
        title="Degradation model Kappa (Modes "+strtrim(string(modes[0]),2)+', '+strtrim(string(modes[1]),2)+')'
        if n_elements(coeff) gt 0 then begin
            fit_label='Poly Fit '
            for i=0,n_elements(coeff)-1 do fit_label=fit_label+strtrim(string(coeff[i]),2)+', '
            plot_multi,waves,pdkappa,ww,kappa_fit,ww,pfit,/xst,/yst,title=title,yrange=[0d,0.004d],$
                xtitle="Wavelength (nm)",ytitle="Kappa",charsize=1.4,label=["Kappa","BSpline Fit",fit_label],psym=[-4,-3,-3,-3]
        endif else begin
            plot_multi,waves,pdkappa,ww,kappa_fit,/xst,/yst,title=title,yrange=[0d,0.004d],$
                xtitle="Wavelength (nm)",ytitle="Kappa",charsize=1.4,label=["Kappa","BSpline Fit"],psym=[-4,-3,-3,-3]
        endelse
    endif

    ;if modes[0] eq 41 then kappa_fit=[pfit[k],pfit_1[ck]]
    return, {wavelength:waves, kappa:pdkappa, x:ww, y:kappa_fit}

END

