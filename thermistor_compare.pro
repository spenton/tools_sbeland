; comparing the results from the 3 methods to transform the telemetry DN
; to degrees Celsius:
;   1) get_sorce_telemetry,
;   2) Java code
;   3) thermistor_convert from Jerry
;
pro thermistor_compare, dn, prism=prism, diode=diode, outdata=outdata, noplot=noplot

    if keyword_set(prism) then begin
        ;define the full range of DN
        if n_elements(dn) gt 0 then alldn=dn else alldn=dindgen((19500d -5000d)/10d)*10d - 19500d

        scaling = -20d / 2d^16
        vout=scaling * alldn

        ; Jerry's code
        tmpdn=alldn
        convertData = Thermistor_convert(tmpdn)
        
        ; get_sorce_telemetry
        query = "select c0,c1,c2,c3,c4,c5,c6,c7 from SORCE_CT.dbo.TelemetryAnalogConversions where tlmId=101909"
        query_database,query, dbcoeffs
        if size(dbcoeffs,/tname) ne 'STRUCT' then begin
            coeffs=[101.27743d, 0.021881277d, 2.8301190d-06, 2.2539677d-10, $
                    9.1379575d-15, 1.4924627d-19, 0d, 0d]
        endif else begin
            coeffs = [dbCoeffs.c0, dbCoeffs.c1, dbCoeffs.c2, dbCoeffs.c3, dbCoeffs.c4, dbCoeffs.c5, dbCoeffs.c6, dbCoeffs.c7]
        endelse
        get_sorceData=poly(alldn, coeffs)

        ;Java code
        query="SELECT * FROM dbo.TimSimThermistorCalibration where thermistorItemName='prism_drive_temp'"
        query+=" and instrumentName='sim_a' and version=1"
        query_database, query, javaCal
        if size(javaCal,/tname) ne 'STRUCT' then begin
           ; no database access - force values
           javaCal={resistance:19600d, referenceVoltage:7.17d, $
                    a:0.0011292410d, b:0.00023410770, c:8.7754680d-08, zeroCelsius:273.15d}
        endif
        rt=vout * javaCal.resistance / (javaCal.referenceVoltage - vout)
        javaData = 1d / (javaCal.a + javaCal.b * alog(rt) + javaCal.c * alog(rt)^3d) - javaCal.zerocelsius
        
    endif else begin   ; processing photo-diode thermistors
        ;define the full range of DN
        if n_elements(dn) gt 0 then alldn=dn else alldn=dindgen(2d^16 / 10d)*10d - 2d^15d

        scaling = -20d / 2d^16
        vout=scaling * alldn

        ; Jerry's code
        tmpdn=alldn
        convertData = minTherm_convert(tmpdn)
        
        ; get_sorce_telemetry
        query = "select c0,c1,c2,c3,c4,c5,c6,c7 from SORCE_CT.dbo.TelemetryAnalogConversions where tlmId=100778"
        query_database,query, dbcoeffs
        if size(dbcoeffs,/tname) ne 'STRUCT' then begin
            coeffs=[6.8650700d, 0.00077011700d, 1.9561100d-09, 8.7549100d-14, $
                    8.8249700d-19, 4.0919100d-23, 0d, 0d]
        endif else begin
            coeffs = [dbCoeffs.c0, dbCoeffs.c1, dbCoeffs.c2, dbCoeffs.c3, dbCoeffs.c4, dbCoeffs.c5, dbCoeffs.c6, dbCoeffs.c7]
        endelse
        get_sorceData=poly(alldn, coeffs)

        ;Java code
        query="SELECT * FROM dbo.TimSimThermistorCalibration where thermistorItemName='vis_1_temp'"
        query+=" and instrumentName='sim_a' and version=1"
        query_database, query, javaCal
        if size(javaCal,/tname) ne 'STRUCT' then begin
           ; no database access - force values
           javaCal={resistance:23340d, referenceVoltage:2.0d, gain:7.98d, $
                    a:0.0011292410d, b:0.00023410770, c:8.7754680d-08, zeroCelsius:273.15d}
        endif
        rt = javaCal.resistance * (javaCal.referenceVoltage * javaCal.gain + vout) 
        rt /= (javaCal.referenceVoltage * javaCal.gain - vout)
        javaData = 1d / (javaCal.a + javaCal.b * alog(rt) + javaCal.c * alog(rt)^3d) - javaCal.zerocelsius
     endelse

    outdata={convertData:convertData, get_sorceData:get_sorceData, javaData:javaData}

    if NOT keyword_set(noplot)  and n_elements(alldn) gt 1 then begin
        title='Comparing Conversion from DN to Degrees C for the PHOTODIODE'
        if keyword_set(prism) then title='Comparing Conversion from DN to Degrees C for the PRISM'
        plot_multi,alldn,javadata,alldn,convertData,alldn,get_sorceData,/xst,/yst,xtitle='Thermistor DN',$
            ytitle='Temperature (C)',title=title, charsize=1.4, $
            label=['Java Code', 'Thermistor_Convert','get_sorce_telemetry'], psym=[-3,-3,-3]

        ; plot the differences compared to the JavaData
        title='Differences in converted temperature with Java code for the PHOTODIODE'
        if keyword_set(prism) then title='Differences in converted temperature with Java code for the PRISM'
        plot_multi,javaData,javadata-convertData,javaData,javaData-get_sorceData,/xst,/yst,xtitle='Temperature (C)',$
            ytitle='Differences in Temperature (C)',title=title, charsize=1.4, $
            label=['Java - Thermistor_Convert','Java - get_sorce_telemetry'], psym=[1,-3],colors=[3,4]
    endif
end
