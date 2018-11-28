;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Starting with a known diode temperature correction, tweak the correction 
;   to flatten the irradiance over the specified time range at each wavelength
;   when plotting the irradiance as a function of temperature.
; TODO test with uncorrected and corrected irradiance
;
; CALLING SEQUENCE:
;   result = GET_DIODE_TEMPCORR_FROM_IRRAD(instrumentModeId, starttime=starttime, stoptime=stoptime)
;
; INPUT PARAMETERS:
;   instrumentModeId -
;      41 for VIS1, 43 for UV or 44 for IR (SimA)
;
; OPTIONAL INPUT PARAMETERS:
;   startTime -
;      The lower time range for which data will be returned.
;      May be specified as julian days, gps microseconds, or
;      SORCE mission day numbers.
;
;   stopTime -
;      The upper time range for which data will be returned
;      May be specified as julian days, gps microseconds, or
;      SORCE mission day numbers.
;   gps - 
;      If set, indicates that input times are to be interpreted as
;      microseconds since Jan 6, 1980 midnight UT.  (default)
;   missionDays - 
;      If set, indicates that input times are to be interpreted as
;      days elapsed since launch, i.e. mission days.  If only the startTime
;      argument contains a non-zero value (stopTime is not supplied), then
;      all data for the specified mission day number will be retrieved.
;   julianDays - 
;      If set, indicates that the input times are to be interpreted as
;      julian day numbers.
;
; RETURNED PARAMETERS:
;   
;
; OPTIONAL OUTPUT PARAMETERS:
;      NONE
;
; EXAMPLE:  
;
; COMMON BLOCKS:
;      NONE
;
;-
;*******************************************************************
function GET_DIODE_TEMPCORR_FROM_IRRAD, instrumentModeId, starttime=startTime, stoptime=stopTime, $
         gps=gps, missionDays=missionDays, julianDays=julianDays, version=version, dburl=dburl, diodetemp=diodetemp, $
         corrfile=corrfile, alldata=alldata, refwave=refwave, profile_data=profile_data, noplot=noplot, verbose=verbose

    ; convert the time range in mission days
    if n_elements(starttime) eq 0 or n_elements(stoptime) eq 0 then begin
        ;starttime = 1577.7970d
        ;stoptime =  1578.4390d
        ;stoptime =  1578.0842d
        ;starttime=3646d
        ;stoptime=3668d
        ;starttime=4271d
        ;stoptime=4470d
        ;starttime=5553d
        ;stoptime=5733d
        starttime=4500d
        stoptime=4800d
        missionDays=1
    endif

    if keyword_set(missionDays) then begin
       ; do nothing here since already in mission days
    endif else if keyword_set(julianDays) then begin
       ; user specified time in julian days
       starttime = jd2sd(starttime)
       stoptime = jd2sd(stoptime)
    endif else begin
       ; user specified time in gps micro seconds
       starttime = gps2sd(starttime/1.d6)
       stoptime = gps2sd(starttime/1.d6)
    endelse

    if instrumentModeId lt 41 or instrumentModeId gt 44 then begin
        print,'Error: invalid instrumentModeId (should be between 41 and 45)'
        return,-1
    endif

    if n_elements(user) eq 0 then     user='sbeland'
    if n_elements(password) eq 0 then password='sbeland1'
    if n_elements(dbdriver) eq 0 then dbDriver='com.sybase.jdbc3.jdbc.SybDriver'
    ;if n_elements(version) eq 0 then version=1005
    ;if n_elements(dburl) eq 0 then    dbUrl='jdbc:sybase:Tds:lasp-db-dev:4100/SORCE_SIM_V20'
    ;if n_elements(dburl) eq 0 then    dbUrl='jdbc:sybase:Tds:lasp-db-dev:4100/SORCE_SIM_DEV'
    ;if n_elements(corrfile) eq 0 then corrfile='/Users/sbeland/SORCE/data/sima_vis_tempcorr_20.txt'
    if n_elements(version) eq 0 then version=24
    if n_elements(dburl) eq 0 then   dburl='jdbc:sybase:Tds:lasp-db1:4100/SORCE_SIM_V24'

    ; get the data
    if n_elements(alldata) eq 0 then begin
        ; get the irradiance
        irrad_table='SimCorrectedIrradiance'
        ; do NOT get the avg_doop for this as it re-arranges the scans
        res=compare_19_20(starttime, stoptime,instrumentModeId,/mission,/skip19,v20=version,alldata=alldata, $
            /noplot,irrad_table=irrad_table, /dettemp)

        p=where(ptr_valid(alldata.(0)) eq 1,count)
        if count eq 0 then begin
            print,'Error: no valid data within specified time range'
            return,-1
        endif 
    endif

    ; get the diode temperature from the data structure (compare_19_20 is expected to have been run with /dettemp)
    if n_elements(diodetemp) eq 0 && max((tag_names((*alldata.spect20[0]))).Contains('DETTEMP')) eq 0 then begin
        table2=['ModeledDetectorTemperature']
        diodetemp = get_science_product(table2, starttime, stoptime, instrumentModeId, /mission, version=24, $ ;version, $
            dburl='jdbc:sybase:Tds:lasp-db1:4100/SORCE_SIM_V24') 
    endif else begin
        tt=[]
        ttemp=[]
        for i=0L,n_elements(alldata.(0))-1L do begin
            tt=[tt,(*alldata.spect20[i]).TIMESTAMP]
            ttemp=[ttemp,(*alldata.spect20[i]).DETTEMP]
        endfor
        s=sort(tt)
        tt=tt[s]
        ttemp=ttemp[s]
        diodetemp={MICROSECONDSSINCEGPSEPOCH:temporary(tt), TEMPERATURE:temporary(ttemp)}
    endelse

    ; get the applied diode temperature corrections
    if n_elements(corrfile) gt 0 then readcol,corrfile,wref,drdt,tref,format='(d,d,d)'

   ; get the profile integral for this mode and version
    if size(profile_data,/tname) ne 'STRUCT' then begin
        query1 = "SELECT calibrationSetId,version FROM CalibrationMetadata WHERE "
        query1 = query1+"calibrationTableName='SimProfileIntegralCal' AND instrumentModeId="+strtrim(string(instrumentModeId),2)
        if keyword_set(version) then query1=query1+" AND version<="+strtrim(string(version),2)
        query1 = query1+" ORDER BY version DESC"
        if keyword_set(verbose) then print,'query="'+query1+'"'
        query_database, query1[0], res, nrows, database='SORCE'
        if nrows eq 0 then begin
            print,'Error: no CalibrationMetadata found for the requested date/version'
            return,-1
        endif
        calibrationSetId=res[0].calibrationSetId
        if keyword_set(verbose) then $
            print,'  using SimProfileIntegralCal iMode='+strtrim(string(instrumentModeId),2)+$
            ' version='+strtrim(string(res[0].version),2)+', calibrationSetId='+strtrim(string(calibrationSetId),2)

        ; now get the SimProfileIntegralCal entry for this calibrationSetId
        query2="SELECT * FROM SimProfileIntegralCal WHERE calibrationSetId="+strtrim(string(calibrationSetId),2)
        query_database, query2, profile_data, nrows, database='SORCE'

        if nrows eq 0 then begin
            if keyword_set(verbose) then $
                print,'NO SimProfileIntegralCal entries found for the calibrationSetId='+strtrim(string(calibrationSetId),2)
            return,outdata
        endif
    endif
    
    ; get the reference wavelengths if not provided
    ;if n_elements(refwave) eq 0 then $
    ;    query_database,'SELECT stdWavelength FROM dbo.SimStandardWavelengths where version=3 and instrumentModeId='+string(instrumentModeId),$
    ;    refwave,dburl='jdbc:sybase:Tds:lasp-db1:4100/SORCE'

    ; get the profile integral
    refwave=replicate({stdwavelength:0d}, n_elements(profile_data.y4))
    refwave.stdwavelength=profile_data.y4

    ; get the diode temperature correction from the profile integral table
    if n_elements(wref) eq 0 then begin
        wref=profile_data.y4
        drdt=profile_data.y12
    endif

    if size(alldata,/tname) eq 'STRUCT' then begin
       ; subset the alldata structure to the time range specidfied
       sd=dblarr(n_elements(alldata.spect20))
       for i=0L,n_elements(alldata.spect20)-1L do sd[i]=(*alldata.spect20[i]).timestamp[0]
       pos=where(sd ge sd2gps(starttime)*1d6 and sd le sd2gps(stoptime)*1d6,count)
       if count eq 0 then return,!NULL
       this_data={spect20:alldata.spect20[pos]}
    endif

    ; loop over each wavelength
    nwaves=n_elements(refwave)
    new_drdt=dblarr(nwaves)
    nscans=n_elements(this_data.spect20)

    ; allocate enough space for the arrays
    this_irrad=dblarr(nscans*4)
    this_temp=this_irrad

    for i=0L,nwaves-1L do begin
        ; loop over every scan to get the corresponding times for this wavelength
        pos=-1L
        print,'Processing '+string(refwave[i].STDWAVELENGTH,format='(F0.2)')+' nm'
        for j=0L, nscans-1L do begin
            ss=sort((*this_data.spect20[j]).timestamp)
            ; find every scan crossing this wavelength
            p=where((*this_data.spect20[j]).wavelength[ss[1:-1]] gt refwave[i].STDWAVELENGTH and $
                (*this_data.spect20[j]).wavelength[ss[0:-2]] lt refwave[i].STDWAVELENGTH, count)
            if count eq 0 then continue
            for k=0L,count-1L do begin
                pos+=1L
                this_irrad[pos]=interpol((*this_data.spect20[j]).IRRADIANCE[ss[p[k]:p[k]+1]], (*this_data.spect20[j]).WAVELENGTH[ss[p[k]:p[k]+1]], refwave[i].STDWAVELENGTH) 
                this_time=interpol((*this_data.spect20[j]).TIMESTAMP[ss[p[k]:p[k]+1]], (*this_data.spect20[j]).WAVELENGTH[ss[p[k]:p[k]+1]], refwave[i].STDWAVELENGTH)
                this_temp[pos]=interpol(diodetemp.temperature, diodetemp.MICROSECONDSSINCEGPSEPOCH, this_time)
            endfor
        endfor

        cc=ladfit(this_temp[0:pos],this_irrad[0:pos])
        plot,this_temp[0:pos],this_irrad[0:pos],/xst,/yst,color=0,title=string(refwave[i].stdwavelength,format='(F0.2)'),/nodata
        oplot,this_temp[0:pos],this_irrad[0:pos],psym=5,color=2
        mn=min(this_temp[0:pos],max=mx)
        oplot,[mn,mx], poly([mn,mx],cc),psym=-3,color=3
        new_drdt[i] = interpol(drdt,wref,refwave[i].stdwavelength) + cc[1]
    endfor

    if n_elements(alldata) eq 0 then alldata=temporary(this_data)

    ; smooth out the data
    outdata = {wavelength:refwave.stdwavelength, drdt:new_drdt}
    p=where(new_drdt ne 0d,count)
    if instrumentModeId eq 44 then begin
        yfit=new_drdt
    endif else begin
        sset=bspline_iterfit(outdata.wavelength[p],outdata.drdt[p],requiren=5,bkspace=5, nord=3)
        yfit=bspline_valu(outdata.wavelength,sset)
        resistant_mean,outdata.drdt[p]-yfit[p], 5.0,mean,good=keep
        sset=bspline_iterfit(outdata.wavelength[p[keep]],outdata.drdt[p[keep]],requiren=10,bkspace=10, nord=3)
        yfit=bspline_valu(outdata.wavelength,sset)
    endelse
    outdata=create_struct(outdata,'drdt_fit',yfit)

    title=string(instrumentModeId,starttime,stoptime,format='("SORCE SIM mode=",I0," from ",F0.1," to ",F0.1)')
    plot_multi,wref,drdt,outdata.wavelength,outdata.drdt,outdata.wavelength,outdata.drdt_fit,$
        xtitle='Wavelength (nm)',ytitle='Diode Temperature Correction',title=title, symsize=[1.0,1.5,1.0], $
        psym=[-3,4,-3],charsize=1.4,/xst,/yst, label=['Original drdt','Optimized drdt','BSpline drdt'],thick=[3.0,1.0,3.0]

    return, outdata
   
end
