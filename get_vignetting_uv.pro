;+
; GET_VIGNETTING_UV
;  
; Obtain the fraction of light vignetted by the calibration mirror
; on the UV diode beam.
; We do this by compairing consecutive spectrum with and without the 
; calibration mirror in the beam.
; It turns out this is a function of CCD position (and not a constant
; as we expected).
;
; As soon as this problem was discovered, the activity templates were
; modified to make sure the calibration mirror was always stowed after
; being used.  The situation was fixed around mission day 192
;-
function get_vignetting_uv, sima=sima, simb=simb, noplot=noplot, ratios=ratios

    if not keyword_set(sima) and not keyword_set(simb) then sima=1

    if keyword_set(simb) then begin
        ; define the times when the calibration mirror was in and out
        ; note that we only have 4 scans early in the mission with cal_in=1
        ; and these are before the first cold event on day 39 which 
        ; might have affected the CCD shift (which means we might be looking
        ; at different wavelengths when comparing the two scans)
        ; let's see what we get ...
        modeId=47
        ptitle='SimB UV Vignetting'
        ; for SimB we found, because of a cold event, a shift between the 2 
        ; datasets
        cal_out_times = [$
            [38.507630d, 38.524162d], $
            [38.524520d, 38.540850d], $
            [38.642550d, 38.658870d], $
            [38.659230d, 38.675550d], $

            [38.507630d, 38.524162d], $
            [38.524520d, 38.540850d], $
            [38.642550d, 38.658870d], $
            [38.659230d, 38.675550d]]
        cal_in_times = [$
            [40.272000d, 40.292000d], $
            [40.272000d, 40.292000d], $
            [40.272000d, 40.292000d], $
            [40.272000d, 40.292000d], $

            [40.340270d, 40.355000d], $
            [40.340270d, 40.355000d], $
            [40.340270d, 40.355000d], $
            [40.340270d, 40.355000d]]
    endif

    if keyword_set(sima) then begin
        ; define the times when the calibration mirror was in and out
        modeId=43
        ptitle='SimA UV Vignetting'
        cal_out_times = [$
            [118.15938d, 118.17622d], $
            [128.62713d, 128.64397d], $
            [132.67811d, 132.69495d], $
            [140.44567d, 140.46251d], $
            [142.94902d, 142.96587d], $
            [150.18330d, 150.20014d], $
            [164.09587d, 164.11271d], $
            [171.18819d, 171.20503d], $
            [191.10953d, 191.12637d]]
        cal_in_times = [$
            [118.36208d, 118.37887d], $
            [128.42456d, 128.44140d], $
            [132.88066d, 132.89750d], $
            [140.64855d, 140.66539d], $
            [142.67833d, 142.69517d], $
            [150.38597d, 150.40281d], $
            [164.16341d, 164.18025d], $
            [171.39083d, 171.40767d], $
            [190.63677d, 190.65361d]]
    endif


    if n_elements(cal_out_times) ne n_elements(cal_in_times) then begin
        print,'Error: not the same number of spectrum with and without vignetting'
        return,-1
    endif

    nspect = n_elements(cal_out_times[0,*])
    ratios = ptrarr(nspect)

    ;get the "standard" UV spectra and get the corresponding CCDPos
    plot_simspect,outdata=spdata,/noplot

    for i=0,nspect-1 do begin
        print,i
        spect_out = get_sim_spectra(cal_out_times[0,i],cal_out_times[1,i],modeId,/mission,/dn_only)
        if size(spect_out,/tname) ne 'STRUCT' then begin
            print,'Spect_Out from ',cal_out_times[0,i],cal_out_times[1,i],' is invalid'
            continue
        endif
        s=sort(spect_out.ccdpos)
        spect_out=spect_out[s]

        spect_in = get_sim_spectra(cal_in_times[0,i],cal_in_times[1,i],modeId,/mission,/dn_only)
        if size(spect_in,/tname) ne 'STRUCT' then begin
            print,'Spect_in from ',cal_in_times[0,i],cal_in_times[1,i],' is invalid'
            continue
        endif
        s=sort(spect_in.ccdpos)
        spect_in=spect_in[s]

        ; make sure we're comparing the same positions
        if n_elements(spect_out) ne n_elements(spect_in) then begin
            print,'Spect_out and Spect_in do not have the same number of elements ',cal_out_times[0,i]
            s=SORT(spect_out.ccdpos)
            spect_out=spect_out[s]
            p=UNIQ(spect_out.ccdpos)
            spect_out=spect_out[p]
            s=SORT(spect_in.ccdpos)
            spect_in=spect_in[s]
            p=UNIQ(spect_in.ccdpos)
            spect_in=spect_in[p]
            match,spect_out.ccdpos,spect_in.ccdpos,subo,subi
            if n_elements(subo) lt 10 then match,spect_out.ccdpos,spect_in.ccdpos,subo,subi,epsilon=8.0
            spect_in=spect_in[subi]
            spect_out=spect_out[subo]
        endif else if max(spect_out.ccdpos - spect_in.ccdpos) gt 1.0 then begin
            print,'Spect_out and Spect_in are not at the same CCD POS ',cal_out_times[0,i]
            continue
        endif

        ratios[i]=ptr_new({ccdpos:spect_out.ccdpos, ratio:spect_in.dn /spect_out.dn})

        if NOT keyword_set(noplot) then $
            lineplot,(*ratios[i]).ccdpos, (*ratios[i]).ratio,psym=-3,title=strtrim(string(cal_in_times[0,i])),charsize=1.5

    endfor

    ; smooth out the noisy regions
    if modeId eq 43 then begin
        meanr=(*ratios[0]).ratio * 0.0d
        count=0.0d
        for i=0,nspect-1 do begin
            if ptr_valid(ratios[i]) then begin
                meanr+=(*ratios[i]).ratio
                count+=1.0d
            endif
        endfor
        meanr /= count
        p0=where((*ratios[0]).ccdpos le 20500.0d)
        meanr[p0]=smooth(meanr[p0],11)
        p1=where((*ratios[0]).ccdpos gt 20500.0d and (*ratios[0]).ccdpos le 30000.0d)
        meanr[p1]=smooth(meanr[p1],7)
        p2=where((*ratios[0]).ccdpos gt 30000.0d and (*ratios[0]).ccdpos le 49670.0d)
        meanr[p2]=smooth(meanr[p2],3)
        ; leave the last few points as the average to keep the transition intact when
        ; the mirror stops vignetting at the largest prism angles (CCDPOS)

        ; vignetting should have a maximum value of 1.0
        p=where(meanr[0:-20] ge 1.0d,count)
        if count gt 0 then begin
            p=max(p)
            meanr[0:p] = 1.0d
        endif
        p=where(meanr gt 1.0d,count)
        if count gt 0 then meanr[p]=1.0d
        outdata={ccdpos:(*ratios[0]).ccdpos, vignette_corr:meanr}

    endif else if modeId eq 47 then begin
        ; for SimB, the ccdpos has a reversed relation to wavelength
        ; smooth out the regions in wavelength with low SNR
        newccdpos=lambda2ccd(ModeId, spdata.uvspect.wavelength, 20d)
        s=sort(newccdpos)
        newccdpos=newccdpos[s]
        meanr=[]
        ccdpos=[]
        for i=0,nspect-1 do begin
            if ptr_valid(ratios[i]) then begin
                meanr=[meanr,(*ratios[i]).ratio]
                ccdpos=[ccdpos,(*ratios[i]).ccdpos]
            endif
        endfor
        p=where(ccdpos ge 27340d and ccdpos le 27640d, complement=cp)
        ccdpos=ccdpos[cp]
        meanr=meanr[cp]
        s=sort(ccdpos)
        ccdpos=ccdpos[s]
        meanr=meanr[s]
        resistant_mean,meanr,5.0,mymean,good=keep
        ccdpos=ccdpos[keep]
        meanr=meanr[keep]
        sset=bspline_iterfit(ccdpos,meanr,maxiter=0,requiren=10,bkspace=5,nord=4)
        meanr=bspline_valu(newccdpos,sset)

        ; vignetting should have a maximum value of 1.0
        p=where(newccdpos gt 29000d and meanr gt 1.0d,count)
        if count gt 0 then begin
            p=min(p)
            meanr[p:-1] = 1.0d
        endif
        p=where(meanr gt 1.0d,count)
        if count gt 0 then meanr[p]=1.0d
        outdata={ccdpos:newccdpos, vignette_corr:meanr}
    endif


    if NOT keyword_set(noplot) then $
            lineplot,outdata.ccdpos, outdata.vignette_corr,psym=-5,title='Smoothed Mean Ratio',xtitle='CCDPOS',$
            ytitle='DN Ratio (Vignetted / Un-Vignetted)', ptitle=ptitle

    return,outdata

end
