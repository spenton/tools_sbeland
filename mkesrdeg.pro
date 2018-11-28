pro mkesrdeg, wavelength, solar54=solar54, solar55=solar55, data1=data1,data2=data2,$
    data3=data3,data4=data4,data5=data5,recovery=recovery, usecoeff=usecoeff

    ; get tnmin_esrdeg_all for each segment at specified wavelength
    t0=[453d, 1570, 2173, 2455, 2803]
    t1=[1570d, 2173, 2455, 2803, 3600]

    res1 = tnmin_esrdeg_all(t0[0],t1[0], wavelength, /mission, coeffs=coeff1, inspectrum=data1, $
        solar54=solar54, solar55=solar55, fit_goodness=good1, /submean, recovery=recovery)
    if keyword_set(usecoeff) then usecoeff=coeff1
    res2 = tnmin_esrdeg_all(t0[1],t1[1], wavelength, /mission, coeffs=coeff2, inspectrum=data2, $
        solar54=solar54, solar55=solar55, fit_goodness=good2, /submean, recovery=recovery, usecoeff=usecoeff)
    res3 = tnmin_esrdeg_all(t0[2],t1[2], wavelength, /mission, coeffs=coeff3, inspectrum=data3, $
        solar54=solar54, solar55=solar55, fit_goodness=good3, /submean, recovery=recovery, usecoeff=usecoeff)
    res4 = tnmin_esrdeg_all(t0[3],t1[3], wavelength, /mission, coeffs=coeff4, inspectrum=data4, $
        solar54=solar54, solar55=solar55, fit_goodness=good4, /submean, recovery=recovery, usecoeff=usecoeff)
    res5 = tnmin_esrdeg_all(t0[4],t1[4], wavelength, /mission, coeffs=coeff5, inspectrum=data5, $
        solar54=solar54, solar55=solar55, fit_goodness=good5, /submean, recovery=recovery, usecoeff=usecoeff)

    timeA = gps2sd([res1.sp31.(0), res2.sp31.(0), res3.sp31.(0), res4.sp31.(0), res5.sp31.(0)]/1d6)
    timeB = gps2sd([res1.sp32.(0), res2.sp32.(0), res3.sp32.(0), res4.sp32.(0), res5.sp32.(0)]/1d6)

    irdA = [res1.sp31.irradiance, res2.sp31.irradiance, res3.sp31.irradiance, res4.sp31.irradiance, res5.sp31.irradiance]
    new_irdA = [res1.sp31.new_irrad, res2.sp31.new_irrad, res3.sp31.new_irrad, res4.sp31.new_irrad, res5.sp31.new_irrad]
   
    irdB = [res1.sp32.irradiance, res2.sp32.irradiance, res3.sp32.irradiance, res4.sp32.irradiance, res5.sp32.irradiance]
    new_irdB = [res1.sp32.new_irrad, res2.sp32.new_irrad, res3.sp32.new_irrad, res4.sp32.new_irrad, res5.sp32.new_irrad]

    plot_multi, timeA, new_irdA, timeB, new_irdB, timeA, irdA, timeB, irdB, /xst,/yst,$
        label=['ESRA Calib','ESRB Calib','ESRA','ESRB'],xtitle='Mission Days',ytitle='Irradiance',$
        title='ESRA,B PD @ '+strtrim(string(wavelength,format='(F0.1)'))+'nm', psym=[-4,-4,-7,-7]
 
    ; re-align the various segments
    delta1=mean(res2.sp31.new_irrad[0:8]) - mean(res1.sp31.new_irrad[-8:-1])
    delta2=mean(res3.sp31.new_irrad[0:8]) - mean(res2.sp31.new_irrad[-8:-1])
    delta3=mean(res4.sp31.new_irrad[0:8]) - mean(res3.sp31.new_irrad[-8:-1])
    delta4=mean(res5.sp31.new_irrad[0:8]) - mean(res4.sp31.new_irrad[-8:-1])

    new_irdA = [res1.sp31.new_irrad, res2.sp31.new_irrad-delta1, res3.sp31.new_irrad-delta1-delta2, $
        res4.sp31.new_irrad-delta1-delta2-delta3, res5.sp31.new_irrad-delta1-delta2-delta3-delta4]

    ;delta1=mean(res2.sp32.new_irrad[0:8]) - mean(res1.sp32.new_irrad[-8:-1])
    ;delta2=mean(res3.sp32.new_irrad[0:8]) - mean(res2.sp32.new_irrad[-8:-1])
    ;delta3=mean(res4.sp32.new_irrad[0:8]) - mean(res3.sp32.new_irrad[-8:-1])
    ;delta4=mean(res5.sp32.new_irrad[0:8]) - mean(res4.sp32.new_irrad[-8:-1])

    new_irdB = [res1.sp32.new_irrad, res2.sp32.new_irrad-delta1, res3.sp32.new_irrad-delta1-delta2, $
        res4.sp32.new_irrad-delta1-delta2-delta3, res5.sp32.new_irrad-delta1-delta2-delta3-delta4]

    labels = 'ESRA   kappa=['+strtrim(string([coeff1,coeff2,coeff3,coeff4,coeff5],format='(5(E10.3,","))'),2)+']'
    labels=[labels,'ESRB   goodness=['+strtrim(string([good1,good2,good3,good4,good5],format='(5(E10.3,","))'),2)+']']
    plot_multi, timeA, new_irdA, timeB, new_irdB, /xst,/yst,$
        label=labels,xtitle='Mission Days',ytitle='Irradiance',$
        title='Aligned ESRA,B @ '+strtrim(string(wavelength,format='(F0.1)'))+'nm'


end
