function get_profile_b, imode, fitcoeffa=fitcoeffa, fitcoeffb=fitcoeffb, profa=profa, profb=profb, solexpb=solexpb, version=version

if n_elements(imode) eq 0 then imode=[41,45]
if n_elements(version) eq 0 then version=22

if n_elements(fitcoeffa) eq 0 then begin
    if imode[0] eq 41 then begin
        readcol,'~/SORCE/data/sima_vis_ccdshift_33.txt',modeid,st0,int_period,coeff0,coeff1,coeff2,coeff3,format='(D)'
        ;readcol,'~/SORCE/data/sima_vis_ccdshift_33_2.txt',modeid,st0,int_period,coeff0,coeff1,coeff2,coeff3,format='(D)'
        fitcoeffa=[coeff0,coeff1,coeff2,coeff3]
    endif else if imode[0] eq 43 then begin
        readcol,'~/SORCE/data/sima_uv_ccdshift_33.txt',modeid,st0,int_period,coeff0,coeff1,coeff2,coeff3,format='(D)'
        fitcoeffa=[coeff0,coeff1,coeff2,coeff3]
    endif
endif 
if n_elements(fitcoeffb) eq 0 then begin
    if imode[0] eq 41 then begin
        ;readcol,'~/SORCE/data/simb_vis_ccdshift_33.txt',modeid,st0,int_period,coeff0,coeff1,coeff2,coeff3,format='(D)'
        ;readcol,'~/SORCE/data/simb_vis_ccdshift_33_2.txt',modeid,st0,int_period,coeff0,coeff1,coeff2,coeff3,format='(D)'
        readcol,'~/SORCE/data/simb_vis_ccdshift_33_3.txt',modeid,st0,int_period,coeff0,coeff1,coeff2,coeff3,format='(D)'
        fitcoeffb=[coeff0,coeff1,coeff2,coeff3]
    endif else if imode[0] eq 43 then begin
        readcol,'~/SORCE/data/simb_uv_ccdshift_33.txt',modeid,st0,int_period,coeff0,coeff1,coeff2,coeff3,format='(D)'
        fitcoeffb=[coeff0,coeff1,coeff2,coeff3]
    endif
endif

pl=get_sorce_plan(32.8,35.0,/sima,act='SolarQuickScan')
scanid =0
specta=get_sim_spectra(pl[scanid].starttime,pl[scanid].stoptime,imode[0],/mission,version=version,fitcoeff=fitcoeffa,/no1au,/verbose,/uncorr,profile_data=profa)
spectb=get_sim_spectra(pl[scanid].starttime,pl[scanid].stoptime,imode[1],/mission,version=version,fitcoeff=fitcoeffb,/no1au,/verbose,/uncorr,profile_data=profb)

if imode[0] eq 41 then begin
    ; remove glint
    p=where(specta.ccdpos_uncorr gt 12180d and specta.ccdpos_uncorr lt 12470d,count,comp=cp)
    specta=specta[cp]
    p=where(spectb.ccdpos_uncorr gt 50200d and spectb.ccdpos_uncorr lt 50500d,count,comp=cp)
    spectb=spectb[cp]
endif else if imode[0] eq 43 then begin
    ; remove glint
    p=where(specta.ccdpos_uncorr gt 34650d and specta.ccdpos_uncorr lt 34900d,count,comp=cp)
    specta=specta[cp]
    p=where(spectb.ccdpos_uncorr gt 27530d and spectb.ccdpos_uncorr lt 27780d,count,comp=cp)
    spectb=spectb[cp]
endif

newvala=interpol(specta.irradiance,specta.wavelength, spectb.wavelength,/spl)

if n_elements(solexpb) gt 0 then begin
    ; provided with an extra solar exposure for SimB
    ; calculate the prism-degradation for SimB with Kappa measure on day 371 (best fit)
    ;restore,'~/SORCE/data/kappa_371_spline.sav'
    restore,'~/SORCE/data/kappa3_371_logfit.sav'
    ; get the raypath for VIS
    readcol,'~/SORCE/data/overlap_vis1_smooth.txt', rayw, raya, format='(d,d)'
    kval = interpol(kappa3.kappa, kappa3.wavelength, spectb.wavelength, /spl)
    aval = interpol(raya, rayw, spectb.wavelength, /spl)
    pdb = (1d - aval) * exp(-kval*solexpb) + aval*exp(-kval*solexpb/2d)
    spectb.irradiance /= pdb
endif 

delta = spectb.irradiance / newvala

; sset2=bspline_iterfit(spectb[cp].ccdpos, delta[cp], maxiter=0,requiren=3,bkspace=3,nord=8)
; afit2=bspline_valu(spectb.ccdpos,sset2)

; dd=interpol(afit2,spectb.wavelength,profb.y4,/spl)
; new_profb = profb
; new_profb.y8 = profb.y8 * dd

pb=interpol(profb.y8,profb.y4,spectb.wavelength,/quad)
npb=spectb.irradiance * pb / newvala
if imode[0] eq 41 then begin
    p=where(spectb.wavelength gt 381d and spectb.wavelength lt 382d,comp=cp) 
    new_profb = profb
    new_profb.y8=interpol(npb[cp],spectb[cp].wavelength,new_profb.y4,/quad)
    ; we scale the y8 to make sure SimB is just about the same brightness as SimA again on day 118
    new_profb.y8 *= 0.995d
    sset=bspline_iterfit(new_profb.y4, new_profb.y8, maxiter=0,requiren=10,bkspace=5)
    new_profb.y8=bspline_valu(new_profb.y4, sset)
endif else if imode[0] eq 43 then begin
    p=where(spectb.wavelength gt 200.7d and spectb.wavelength lt 304.0)
    sset2=bspline_iterfit(spectb[p].ccdpos,npb[p], maxiter=0,requiren=30,bkspace=30,nord=8)
    afit2=bspline_valu(spectb.ccdpos,sset2)
    new_profb = profb
    new_profb.y8=interpol(afit2,spectb.wavelength,new_profb.y4,/quad)
endif

stop

return,new_profb

end
