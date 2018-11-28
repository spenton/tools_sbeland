; routine to get the ccdoffset data for a series of days and instrumentmodeID
function test_ccdoffset, starttimes, stoptimes, instrumentModes, version=version, $
                         noplot=noplot,plotfit=plotfit, noTempCorr=noTempCorr, noCcdCorr=noCcdCorr, $
                         refT0=refT0, refT1=refT1, order=order, avgfit=avgfit, dowave=dowave

    if n_elements(starttimes) ne n_elements(stoptimes) then begin
        print,'Error: starttimes and stoptimes don"t have the same number of elements'
        return,-1
    endif

    if keyword_set(dowave) then noTempCorr=0
    ntimes=n_elements(starttimes)
    nmodes=n_elements(instrumentModes)
    nreftimes = n_elements(refT0)
    if nreftimes gt 1 and (nreftimes ne ntimes) then begin
        print,'Error: starttimes and refT0 are not of the same size'
        return,-1
    endif
    outdata = ptrarr(ntimes * nmodes)
    if n_elements(order) eq 0 then order=3 else order=order < 5
    goodfit = 0.95
    if order lt 1 then coeffs = dblarr(ntimes,2) else coeffs = dblarr(ntimes,order+1)
    tvlct,r,g,b,/get
    ; keep the reference spectra for each mode in memory 
    refSpect = ptrarr(n_elements(instrumentModes))

    for tt=0,ntimes-1 do begin
        print,'Processing data from ',starttimes[tt]
        ccdpos = []
        offsets=[]
        for mode=0,nmodes-1 do begin
            print,'    mode=',instrumentModes[mode]
            if ptr_valid(refSpect[mode]) and nreftimes le 1 then begin
                ; same reference spectra 
                result = get_ccdoffset(starttimes[tt], stoptimes[tt], refSpect=(*refSpect[mode]),$
                    instrumentModes[mode],/mission, version=version, noTempCorr=noTempCorr, $
                    noCcdCorr=noCcdCorr, /quiet, dowave=dowave)
            endif else if nreftimes ge 1 then begin
                ; reference spectra will be different every time
                result = get_ccdoffset(starttimes[tt], stoptimes[tt], refSpect=thisref, $
                    instrumentModes[mode],/mission, version=version, noTempCorr=noTempCorr, $
                    noCcdCorr=noCcdCorr, /quiet, refT0=refT0[tt], refT1=refT1[tt], dowave=dowave)
                if ptr_valid(refSpect[mode]) then ptr_free,refSpect[mode]
                refSpect[mode] = PTR_NEW(thisref,/no_copy)
            endif else begin
                ; use the default reference spectra
                result = get_ccdoffset(starttimes[tt], stoptimes[tt], refSpect=thisref,$
                    instrumentModes[mode],/mission, version=version, noTempCorr=noTempCorr, $
                    noCcdCorr=noCcdCorr, /quiet, dowave=dowave)
                if ptr_valid(refSpect[mode]) then ptr_free,refSpect[mode]
                refSpect[mode] = PTR_NEW(thisref,/no_copy)
            endelse
            if size(result,/tname) ne 'STRUCT' then begin
                print,'NO valid fit found for ',starttimes[tt],' ',instrumentModes[mode]
                continue
            endif
            p=where(result.deriv_corr ge goodfit,count)
            if count eq 0 then begin
                print,'NO valid fit found for ',starttimes[tt],' ',instrumentModes[mode]
                continue
            endif else result=result[p]
            ccdpos = [ccdpos,result.meanccdpos]
            offsets = [offsets,result.deriv_offset]
            outdata[tt*nmodes + mode] = PTR_NEW(result,/no_copy)
        endfor

        ; get the 3rd order fit for these modes combined
        if n_elements(ccdpos) eq 0 then continue
        s=sort(ccdpos)
        ccdpos=ccdpos[s]
        offsets=offsets[s]
        if order lt 1 then begin
           fit_coef = [median(offsets), 0.0d] 
        endif else if n_elements(ccdpos) lt order+1 then begin
           fit_coef = robust_poly_fit(ccdpos, offsets, order-n_elements(offsets), /double)
           fit_coef=[fit_coef, replicate(0d, order-n_elements(offsets)+1)]
        endif else begin
            fit_coef = robust_poly_fit(ccdpos, offsets, order, /double)
        endelse
        coeffs[tt,*] = fit_coef
        if NOT keyword_set(noplot) then begin
            yfit=poly(ccdpos,fit_coef)
            title=strtrim(string(starttimes[tt],format='(F7.2)'),2)
            if tt eq 0 then begin
                lineplot,ccdpos,offsets,title=title,xtitle='CCD Position (Subpixels)',$
                ytitle='DERIV Cross-Correlation Offset (Subpixels)',psym=-3
                if keyword_set(plotfit) then lineplot,ccdpos,yfit,title=title+' FIT',psym=-3
            endif else begin
                color=[r[byte(tt+2)],g[byte(tt+2)],b[byte(tt+2)]]
                lineplot,ccdpos,offsets,title=title,xtitle='CCD Position (Subpixels)',$
                ytitle='DERIV Cross-Correlation Offset (Subpixels)',psym=-3,color=color
                if keyword_set(plotfit) then lineplot,ccdpos,yfit,title=title+' FIT',psym=-3,color=color
            endelse
        endif
    endfor

    ; calculate the fit by combining all the data
    p=where(ptr_valid(outdata) eq 1)
    data=outdata[p]

    deriv_ccdpos = []
    deriv_offsets=[]
    for i=0,n_elements(data)-1 do begin
        deriv_ccdpos=[deriv_ccdpos,(*data[i]).meanccdpos]
        deriv_offsets=[deriv_offsets,(*data[i]).deriv_offset]
    endfor
    avgfit=robust_poly_fit(deriv_ccdpos, deriv_offsets, order, /double)
    if NOT keyword_set(noplot) then lineplot,ccdpos,poly(ccdpos,avgfit),title='Fit of all Data',psym=-4,color=3

    return,{fitCoeff:coeffs[p,*], starttimes:starttimes[p], stoptimes:stoptimes[p], data:outdata[p]}

end
