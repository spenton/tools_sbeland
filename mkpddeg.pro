pro mkpddeg, modes, wavelength, solar54=solar54, solar55=solar55, data1=data1,data2=data2,data3=data3,data4=data4,data5=data5

    ; get tnmin_pddeg for each segment at specified wavelength
    t0=[453d, 1570, 2173, 2455, 2803]
    t1=[1570d, 2173, 2455, 2803, 3600]

    res1 = tnmin_pddeg(t0[0],t1[0], modes, wavelength, coeffs=coeff1, inspectrum=data1, $
        solar54=solar54, solar55=solar55, fit_goodness=good1, /submean)
    res2 = tnmin_pddeg(t0[1],t1[1], modes, wavelength, coeffs=coeff2, inspectrum=data2, $
        solar54=solar54, solar55=solar55, fit_goodness=good2, /submean)
    res3 = tnmin_pddeg(t0[2],t1[2], modes, wavelength, coeffs=coeff3, inspectrum=data3, $
        solar54=solar54, solar55=solar55, fit_goodness=good3, /submean)
    res4 = tnmin_pddeg(t0[3],t1[3], modes, wavelength, coeffs=coeff4, inspectrum=data4, $
        solar54=solar54, solar55=solar55, fit_goodness=good4, /submean)
    res5 = tnmin_pddeg(t0[4],t1[4], modes, wavelength, coeffs=coeff5, inspectrum=data5, $
        solar54=solar54, solar55=solar55, fit_goodness=good5, /submean)

    timeA = gps2sd([res1.sima.(0), res2.sima.(0), res3.sima.(0), res4.sima.(0), res5.sima.(0)]/1d6)
    timeB = gps2sd([res1.simb.(0), res2.simb.(0), res3.simb.(0), res4.simb.(0), res5.simb.(0)]/1d6)

    irdA = [res1.sima.irradiance, res2.sima.irradiance, res3.sima.irradiance, res4.sima.irradiance, res5.sima.irradiance]
    new_irdA = [res1.sima.new_irrad, res2.sima.new_irrad, res3.sima.new_irrad, res4.sima.new_irrad, res5.sima.new_irrad]
   
    irdB = [res1.simb.irradiance, res2.simb.irradiance, res3.simb.irradiance, res4.simb.irradiance, res5.simb.irradiance]
    new_irdB = [res1.simb.new_irrad, res2.simb.new_irrad, res3.simb.new_irrad, res4.simb.new_irrad, res5.simb.new_irrad]

    plot_multi, timeA, new_irdA, timeB, new_irdB, timeA, irdA, timeB, irdB, /xst,/yst,$
        label=['SimA Calib','SimB Calib','SimA','SimB'],xtitle='Mission Days',ytitle='Irradiance',$
        title='SimA,B PD @ '+strtrim(string(wavelength,format='(F0.1)'))+'nm'
 
    ; re-align the various segments
    delta1=mean(res2.sima.new_irrad[0:8]) - mean(res1.sima.new_irrad[-8:-1])
    delta2=mean(res3.sima.new_irrad[0:8]) - mean(res2.sima.new_irrad[-8:-1])
    delta3=mean(res4.sima.new_irrad[0:8]) - mean(res3.sima.new_irrad[-8:-1])
    delta4=mean(res5.sima.new_irrad[0:8]) - mean(res4.sima.new_irrad[-8:-1])

    new_irdA = [res1.sima.new_irrad, res2.sima.new_irrad-delta1, res3.sima.new_irrad-delta1-delta2, $
        res4.sima.new_irrad-delta1-delta2-delta3, res5.sima.new_irrad-delta1-delta2-delta3-delta4]

    ;delta1=mean(res2.simb.new_irrad[0:8]) - mean(res1.simb.new_irrad[-8:-1])
    ;delta2=mean(res3.simb.new_irrad[0:8]) - mean(res2.simb.new_irrad[-8:-1])
    ;delta3=mean(res4.simb.new_irrad[0:8]) - mean(res3.simb.new_irrad[-8:-1])
    ;delta4=mean(res5.simb.new_irrad[0:8]) - mean(res4.simb.new_irrad[-8:-1])

    new_irdB = [res1.simb.new_irrad, res2.simb.new_irrad-delta1, res3.simb.new_irrad-delta1-delta2, $
        res4.simb.new_irrad-delta1-delta2-delta3, res5.simb.new_irrad-delta1-delta2-delta3-delta4]

    labels = 'k=['+strtrim(string([coeff1,coeff2,coeff3,coeff4,coeff5],format='(5(E10.3,","))'),2)+']'
    labels=[labels,'g=['+strtrim(string([good1,good2,good3,good4,good5],format='(5(E10.3,","))'),2)+']']
    plot_multi, timeA, new_irdA, timeB, new_irdB, /xst,/yst,$
        label=labels,xtitle='Mission Days',ytitle='Irradiance',$
        title='Aligned SimA,B PD @ '+strtrim(string(wavelength,format='(F0.1)'))+'nm'


end
