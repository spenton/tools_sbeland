;+
; NAME:   INSERT_DEGCOL_WAVE2D
;
; PURPOSE: 
;    This routine populates the SORCE/SimPrismDegColCalWave2D database table with the 
;    new detector degradation column values obtained from get_pdtau2.pro and tnmin_pdtau2.pro.
;
; CALLING SEQUENCE:
;    insert_degcol, 41, /dbinsert
;
; OPTIONAL INPUT PARAMETERS:
;    none
;
; OPTIONAL INPUT KEYWORDS:
;    dbinsert - if set will insert the re-formatted data in the database
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
;    The degradation column values are enetered for each wavelength covered by the instrumentModeId
;    for the corresponding orbit/scan.
;
;    For large files, the dbinsert can be very slow. SYSBASE has a utility called "bcp" installed
;    on the Linux bird boxes that allows a VERY MUCH faster ingest.  Use the outfile parameter
;    for saving an ASCII file in a format that bcp will be able to use.
;
;    The follwoing command can be issued from the unix prompt to ingest the ASCII file:
;       bcp SORCE..SimPrismDegColCalWave2D in outfile.bcp -U sorce -P sorcedb -S sorce-db -c -t"," -r"\n"
;  
; REVISION HISTORY:
;   Revision: $Id: insert_degcol_wave2D.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;-
;*****************************************************************************
PRO update_calibMetadata_degcol_wave2D, instMode, calId, effDate, version=version, $
    ncoeffs=ncoeffs, verbose=verbose, dbinsert=insert, filename=filename

    ; make sure we have an entry in the CalibrationMetadata table for this calibrationSetId
    q1="SELECT * from CalibrationMetadata where calibrationSetId="+strtrim(string(CalId),2)
    if insert then begin
        query_database,q1,data,nrows
        date_str=string(systime(/jul),format='(C(CYI,"-",CMOI2.2,"-",CDI2.2," ",CHI2.2,":",CMI2.2,":",CSI2.2))')
        if nrows eq 0 then begin
            ; insert a new row in the CalibrationMetadata table
            q1="INSERT INTO dbo.CalibrationMetadata(calibrationSetId, calibrationTableName, instrumentModeId, "
            q1+="version, releaseDate, effectiveDate, rationaleForNewVersion, responsibleIndividual, "
            q1+="dataDescription, sourceOfData, applicableReferences, additionalInformation) VALUES("
            q1=q1+strtrim(string(CalId),2)+", 'SimPrismDegColCalWave2D', "
            q1=q1+strtrim(string(instMode),2)+", "+strtrim(string(version),2)
            q1=q1+", '"+date_str+"', '"+effDate+"'"
            q1=q1+", 'degradation column versus wavelength and mission day from solar exposure and F-Function'"
            q1=q1+", 'Stephane Beland'"
            q1=q1+", 'Stores the degradation column representing the surface over the MISSION DAYS and WAVELENGTH'"
            q1=q1+", 'IDL save file "+filename+"'"
            q1=q1+", 'none'"
            q1=q1+", 'none')"
            if keyword_set(verbose) then print,q1
            if insert then query_database, q1, result

        endif else begin
            q1="UPDATE CalibrationMetadata SET releaseDate='"+date_str+"'"
            q1=q1+", effectiveDate='"+effDate+"'"
            q1=q1+", version="+strtrim(string(version),2)
            q1=q1+" WHERE calibrationSetId="+strtrim(string(CalId),2)
            if keyword_set(verbose) then print,q1
            if insert then query_database, q1, result

        endelse
    endif

END

;-----------------------------------------------------------------------------
PRO INSERT_DEGCOL_WAVE2D, instrumentModeId, dbinsert=dbinsert, verbose=verbose, version=version, outfile=outfile

    pos = where(instrumentModeId, count)
    if instrumentModeId lt 31 or instrumentModeId gt 48 then begin
        print,'Error: invalid InstrumentModeId ',instrumentModeId
        return
    endif

    if n_elements(version) eq 0 then version=3

    ; read the save file and re-arrange the coefficients in a 1D array
    if instrumentModeid eq 41 then begin
        ;filename= '~/SORCE/data/sima_vis_degcol_v20_mod_407.sav'
        ;filename= '~/SORCE/data/sima_vis_degcol_v22_mod_407.sav'
        filename= '~/SORCE/data/sima_vis_degcol_v23_mod_407.sav'
        restore, filename
        ;degcol=temporary(degcol_c54)
        ;degcol=temporary(degcol_c54_22)
        degcol=temporary(degcol_c54_23)
    endif else if instrumentModeId eq 45 then begin
        ;filename = '~/SORCE/data/simb_vis_degcol_v20_mod_407.sav'
        ;filename = '~/SORCE/data/simb_vis_degcol_v22_mod_407.sav'
        filename= '~/SORCE/data/simb_vis_degcol_v23_mod_407.sav'
        restore, filename
        ;degcol=temporary(degcol_c55)
        ;degcol=temporary(degcol_c55_22)
        degcol=temporary(degcol_c55_23)
    endif else if instrumentModeid eq 43 then begin
        ;filename='~/SORCE/data/sima_uv_degcol_v20_afact160_s5455mod2_2D_smooth.sav'
        ;filename='~/SORCE/data/sima_uv_degcol_v20_afact160_s5455mod2_2D_allw.sav'
        ;filename='~/SORCE/data/sima_uv_degcol_v22_afact160_s5455mod2_2D_allw.sav'
        ;filename='~/SORCE/data/sima_uv_degcol_v23_afact160_s5455mod2_2D_allw.sav'
        filename='~/SORCE/data/sima_uv_degcol_v20_afact160_s5455mod4_2D_allw.sav'
        restore, filename
        degcol_c54=temporary(degcol_c54_new4)
        ; level out everything below 210nm
        p=where(degcol_c54.wavelength le 210d)
        tmp = degcol_c54.degcol[*,p[-1]]
        for i=0,n_elements(p)-1 do degcol_c54.degcol[*,p[i]] = tmp
        ; reduce the number of wavelengths by 1/3
        k=lindgen(n_elements(degcol_c54.wavelength)/3d)*3
        s=sort(degcol_c54.t1)
        q=uniq(degcol_c54.t1[s])
        if n_elements(q) ne n_elements(s) then begin
            tmp = dblarr(n_elements(q), n_elements(k))
            for i=0L,n_elements(q) -1 do tmp[i,*] = degcol_c54.degcol[s[q[i]],k]
            degcol={wavelength:degcol_c54.wavelength[k], t1:degcol_c54.t1[s[q]], degcol:tmp}
        endif else degcol={wavelength:degcol_c54.wavelength[k], t1:degcol_c54.t1, degcol:degcol_c54.degcol[*,k]}
    endif else if instrumentModeid eq 47 then begin
        ;filename= '~/SORCE/data/simb_uv_degcol_v20_afact160_s5455mod2_2D_smooth.sav'
        ;filename='~/SORCE/data/simb_uv_degcol_v20_afact160_s5455mod2_2D_allw.sav'
        ;filename='~/SORCE/data/simb_uv_degcol_v22_afact160_s5455mod2_2D_allw.sav'
        ;filename='~/SORCE/data/simb_uv_degcol_v23_afact160_s5455mod2_2D_allw.sav'
        filename='~/SORCE/data/sima_uv_degcol_v20_afact160_s5455mod4_2D_allw.sav'
        restore, filename
        degcol_c55=temporary(degcol_c55_new4)
        ; level out everything below 210nm
        p=where(degcol_c55.wavelength le 210d)
        tmp = degcol_c55.degcol[*,p[-1]]
        for i=0,n_elements(p)-1 do degcol_c55.degcol[*,p[i]] = tmp
        ; reduce the number of wavelengths by 1/3
        k=lindgen(n_elements(degcol_c55.wavelength)/3d)*3
        s=sort(degcol_c55.t1)
        q=uniq(degcol_c55.t1[s])
        if n_elements(q) ne n_elements(s) then begin
            tmp = dblarr(n_elements(q), n_elements(k))
            for i=0L,n_elements(q) -1 do tmp[i,*] = degcol_c55.degcol[s[q[i]],k]
            degcol={wavelength:degcol_c55.wavelength[k], t1:degcol_c55.t1[s[q]], degcol:tmp}
        endif else degcol={wavelength:degcol_c55.wavelength[k], t1:degcol_c55.t1, degcol:degcol_c55.degcol[*,k]}
    endif else if instrumentModeid eq 44 then begin
        ;filename= '~/SORCE/data/sima_uv_degcol_v20_afact160_s5455mod2_2D_smooth.sav'
        ;degcol=temporary(degcol_c54)
        ;filename='~/SORCE/data/solarexp_54_55.sav'
        filename='~/SORCE/data/solarexp_54_55_raw.sav'
        restore, filename
        solar_exp=total(solar54.solar_exp_orbit,/cum)/86400d
        p=where(solar54.solar_exp_orbit gt 0d, count)
        q=uniq(solar54[p].t1)
        count=n_elements(q)
        coeffs=dblarr(count,4)
        coeffs[*,0] = solar_exp[p[q]]
        degcol={wavelength:[820d, 1000d, 1400d, 1800d], t1:solar54[p[q]].t1, coeffs:coeffs}
    endif else if instrumentModeid eq 48 then begin
        ;restore, '~/SORCE/data/sima_uv_degcol_v20_afact160_s5455mod2_2D_smooth.sav'
        ;degcol=temporary(degcol_c55)
        ;filename='~/SORCE/data/solarexp_54_55.sav'
        filename='~/SORCE/data/solarexp_54_55_raw.sav'
        restore, filename
        solar_exp=total(solar55.solar_exp_orbit,/cum)/86400d
        p=where(solar55.solar_exp_orbit gt 0d, count)
        q=uniq(solar55[p].t1)
        count=n_elements(q)
        coeffs=dblarr(count,4)
        coeffs[*,0] = solar_exp[p[q]]
        degcol={wavelength:[820d, 1000d, 1400d, 1800d], t1:solar55[p[q]].t1, coeffs:coeffs}
    endif else if instrumentModeid eq 31 then begin
        ; use the unmodified solar exposure record
        ;restore,'~/SORCE/data/solarexp_54_55.sav'
        filename='~/SORCE/data/solarexp_54_55_mod_407nm.sav'
        restore, filename
        solar_exp=total(solar54.solar_exp_orbit,/cum)/86400d
        p=where(solar54.solar_exp_orbit gt 0d, count)
        q=uniq(solar54[p].t1)
        count=n_elements(q)
        coeffs=dblarr(count,4)
        ; for version 21 set the degcol to 0.0 for all SolarIRScan (since the uncorr looks a lot flatter than corrected)
        ;coeffs[*,0] = solar_exp[p[q]]
        degcol={wavelength:[200d, 820d, 1400d, 1800d, 3000d], t1:solar54[p[q]].t1, coeffs:coeffs}
    endif else if instrumentModeid eq 32 then begin
        ; use the unmodified solar exposure record
        ;restore,'~/SORCE/data/solarexp_54_55.sav'
        filename='~/SORCE/data/solarexp_54_55_mod_407nm.sav'
        restore, filename
        solar_exp=total(solar55.solar_exp_orbit,/cum)/86400d
        p=where(solar55.solar_exp_orbit gt 0d, count)
        q=uniq(solar55[p].t1)
        count=n_elements(q)
        coeffs=dblarr(count,4)
        ; for version 21 set the degcol to 0.0 for all SolarIRScan (since the uncorr looks a lot flatter than corrected)
        ;coeffs[*,0] = solar_exp[p[q]]
        degcol={wavelength:[200d, 820d, 1400d, 1800d, 3000d], t1:solar55[p[q]].t1, coeffs:coeffs}
    endif else begin
        print, 'No data available yet for this instrumentModeId'
        return
    endelse

    ; we have 2 formats for the degradation column save file: one with polynomial fit and the other with smoothing only
    p=where(strpos(tag_names(degcol),'DEGCOL') ge 0,degcol_2D)

    query="select * from CalibrationMetadata where calibrationTableName='SimPrismDegColCalWave2D' and "
    query+=" version="+strtrim(string(version),2)+" and instrumentModeId="+strtrim(string(instrumentModeId),2)

    query_database,query, res, nrows

    if nrows gt 0 then begin
        calId = res[0].CALIBRATIONSETID
    endif else begin
        ; get a new calibrationSetId since the version to update doesn't exist
        query_database,"select max(calibrationSetId) from CalibrationMetadata", maxcal
        calId = maxcal.(0)+1
    endelse
    print,'  inserting data for calibrationSetId=',calId

    if keyword_set(dbinsert) then insert=1 else insert=0

    ; first delete all entries in the table for this particular version (version 1 only for now)
    q1 = "DELETE FROM SimPrismDegColCalWave2D WHERE calibrationSetId="+strtrim(string(calId),2)
    if keyword_set(verbose) then print,q1
    if insert then query_database,q1,result,user='sorce', password='sorcedb', server='sorce-db', database="SORCE"

    ; make sure we have an entry in the CalibrationMetadata table for this calibrationSetId
    effectiveDate='2003-01-01 00:00:00.0'
    update_calibMetadata_degcol_wave2D, instrumentModeId, calId, effectiveDate, version=version, $
        dbinsert=insert, verbose=verbose, filename=file_basename(filename)

    ; insert one row at a time
    q1="INSERT INTO dbo.SimPrismDegColCalWave2D(calibrationSetId, orbitEndGpsMicroseconds, wavelength, degradationColumn) VALUES("

    if n_elements(outfile) gt 0 then begin
        openw,unit,outfile,/get_lun
    endif else unit=0

    norbits=n_elements(degcol.T1)
    nwaves=n_elements(degcol.WAVELENGTH)
    for i=0L,norbits-1 do begin
        print,i+1,double(i+1)/double(norbits+1)*100d,format='("   ",i0,"  (",i0," % )",$,%"\r")'
        q2 = strtrim(string(calId,format='(I)'),2)
        q2 = q2+','+strtrim(string(ulong64(degcol.t1[i]),format='(I)'),2)
        if degcol_2d eq 0 then values = poly(degcol.wavelength,degcol.coeffs[i,*]) > 0d
        for j=0L,nwaves-1 do begin
            q3=q2+','+strtrim(string(degcol.wavelength[j],format='(d0.6)'),2)
            if degcol_2d eq 0 then begin
                q3=q3+','+strtrim(string(values[j],format='(d)'),2)
            endif else begin
                q3=q3+','+strtrim(string(degcol.degcol[i,j],format='(d)'),2)
            endelse
            if keyword_set(verbose) then print,q1+q3+')'
            if unit gt 0 then printf,unit,q3
            ; only write the bcp file from now on
            ;if insert then result=conn17->execute(q1+q3+')')
        endfor

    endfor
    print,''
    if unit then begin
        flush,unit
        close,unit
    endif

    return

END





