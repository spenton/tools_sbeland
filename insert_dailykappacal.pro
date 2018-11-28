;+
; NAME:   INSERT_DAILYKAPPACAL
;
; PURPOSE: 
;    This routine populates the SORCE/SimPrismDegDailyKappaCal database table with the 
;    new Kappa values obtained from comparing the degradation of ESRA and B and the 
;    UVA and UVB.

; CALLING SEQUENCE:
;    insert_dailykappacal, version=23, /dbinsert
;
; OPTIONAL INPUT PARAMETERS:
;    none
;
; OPTIONAL INPUT KEYWORDS:
;    dbinsert - if set will insert the re-formatted data in the database
;
; OUTPUT PARAMETERS:
;   none
;
; OPTIONAL OUTPUT PARAMETERS:
;      NONE
;
; EXAMPLE:  
;
; COMMON BLOCKS:
;      NONE
;
; NOTES:
;
;    For large files, the dbinsert can be very slow. SYSBASE has a utility called "bcp" installed
;    on the Linux bird boxes that allows a VERY MUCH faster ingest.  Use the outfile parameter
;    for saving an ASCII file in a format that bcp will be able to use.
;
;    The follwoing command can be issued from the unix prompt to ingest the ASCII file:
;       bcp SORCE..SimPrismDegDailyKappaCal in outfile.bcp -U sorce -P sorcedb -S sorce-db -c -t"," -r"\n"
;
; REVISION HISTORY:
;   Revision: $Id: insert_dailykappacal.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;-
;*****************************************************************************
PRO update_calibMetadata_dailykappacal, instMode, calId, effDate, version=version, $
    verbose=verbose, dbinsert=insert, filename=filename

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
            q1=q1+strtrim(string(CalId),2)+", 'SimPrismDegDailyKappaCal', "
            q1=q1+strtrim(string(instMode),2)+", "+strtrim(string(version),2)
            q1=q1+", '"+date_str+"', '"+effDate+"'"
            ;q1=q1+", 'release for version 23'"
            q1=q1+", 'test release for version 23'"
            q1=q1+", 'Stephane Beland'"
            q1=q1+", 'Daily variable Kappa calculated on days with coincident SIMA and SIMB measurements'"
            q1=q1+", '"+filename+"'"
            q1=q1+", 'none'"
            q1=q1+", 'changed to daily Kappa for release version 23')"
            if keyword_set(verbose) then print,q1
            if insert then query_database,q1,result
        endif else begin
            q1="UPDATE CalibrationMetadata SET releaseDate='"+date_str+"'"
            q1=q1+", effectiveDate='"+effDate+"'"
            q1=q1+", version="+strtrim(string(version),2)
            q1=q1+", sourceOfData='"+filename+"'"
            q1=q1+" WHERE calibrationSetId="+strtrim(string(CalId),2)
            if keyword_set(verbose) then print,q1
            if insert then query_database,q1,result
        endelse
    endif

END

;*****************************************************************************
PRO INSERT_DAILYKAPPACAL, instrumentModeId, dbinsert=dbinsert, verbose=verbose, version=version, outfile=outfile
    ; insert the list of SIM CCD Shifts in the
    ; database table SORCE.dbo.SimCcdShiftCal

    modes = [31, 32, 41, 43, 44, 45, 47, 48]
    pos = where(modes eq instrumentModeId, count)
    if count eq 0 then begin
        print,'Error: invalid InstrumentModeId ',instrumentModeId,' [',mode,']'
        return
    endif

    if n_elements(version) eq 0 then version=23
    user='sorce'
    password='sorcedb'
    server='sorce-db'
    database='SORCE'

    ; figure out the calibrationSetId (create one if it doesn't exist)
    q1="SELECT * from CalibrationMetadata where calibrationTableName='SimPrismDegDailyKappaCal' "
    q1 += " and instrumentModeId="+strtrim(string(instrumentModeId,format='(I)'),2)
    q1 += " and version="+strtrim(string(version,format='(I)'),2)
    query_database,q1,result,nrows,user=user,password=password,server=server,database=database
    if nrows ge 1 then begin
        newCalId=result.calibrationSetId
        print,'Will update calibrationSetId = ',newcalid
    endif else begin
        ; insert a new row in the CalibrationMetaData table
        query_database,'select max(calibrationSetId) as mx from CalibrationMetadata',result
        newCalId = result.mx + 1
        print,'Will create a new calibrationSetId = ',newcalid
    endelse

    case instrumentModeId of
        31: begin
            if version eq 23 then begin
                infile=[]
            endif
            end
        32: begin
            if version eq 23 then begin
                infile=[]
            endif
            end
        41: begin
            if version eq 23 then begin
                infile='~/SORCE/data/sima_vis_dailyKappa.txt' 
                readcol,infile,mode,ondate,waves,kappa,kappaUnc,delim=",",format='(I,d,d,d,d)'
            endif
            end
        43: begin
            if version eq 23 then begin
                ;infile='~/SORCE/data/simab_uv_kappa_2375.sav' 
                ;infile='~/SORCE/data/simab_uv_kappa_2375_bspline.sav' 
                infile='~/SORCE/data/simab_uv_kappa_2375_poly.sav' 
                restore,infile
                waves=(dblarr(n_elements(kappa_smooth.timesd))+1d) # kappa_smooth.waves
                ondate=(sd2gps(kappa_smooth.timesd)*1d6) # (dblarr(n_elements(kappa_smooth.waves))+1d)
                kappa=kappa_smooth.kappa
                kappaUnc = kappa * 0d
            endif
            end
        44: begin
            if version eq 23 then begin
                infile=[]
            endif
            end
        45: begin
            if version eq 23 then begin
                infile='~/SORCE/data/sima_vis_dailyKappa.txt' 
                readcol,infile,mode,ondate,waves,kappa,kappaUnc,delim=",",format='(I,d,d,d,d)'
            endif
            end
        47: begin
            if version eq 23 then begin
                ;infile='~/SORCE/data/simab_uv_kappa_2375.sav' 
                ;infile='~/SORCE/data/simab_uv_kappa_2375_bspline.sav' 
                infile='~/SORCE/data/simab_uv_kappa_2375_poly.sav' 
                restore,infile
                waves=(dblarr(n_elements(kappa_smooth.timesd))+1d) # kappa_smooth.waves
                ondate=(sd2gps(kappa_smooth.timesd)*1d6) # (dblarr(n_elements(kappa_smooth.waves))+1d)
                kappa=kappa_smooth.kappa
                kappaUnc = kappa * 0d
            endif
            end
        48: begin
            if version eq 23 then begin
                infile=[]
            endif
            end
    endcase

    if n_elements(infile) eq 0 then begin
        print,'Error: no matching input file for mode ',instrumentModeId
        return
    endif

    if keyword_set(dbinsert) then insert=1 else insert=0

    ; first delete all entries in the table for this particular version
    q1 = "DELETE FROM SimPrismDegDailyKappaCal WHERE calibrationSetId="+strtrim(string(newcalId),2)
    if keyword_set(verbose) then print,q1
    if insert then query_database,q1

    ; make sure we have an entry in the CalibrationMetadata table for this calibrationSetId
    effectiveDate='2003-01-24 00:00:01.0'
    update_calibMetadata_dailykappacal, instrumentModeId, newCalId, effectiveDate, $
        version=version, dbinsert=insert, verbose=verbose, filename=file_basename(infile)

    ; insert one row at a time
    q1 = "INSERT INTO SimPrismDegDailyKappaCal(calibrationSetId, effectiveDateGps, wavelength, kappa, kappaUncertainty ) VALUES("

    print,strtrim(string(n_elements(ondate)),2)+' rows to insert as calibrationSetId='+strtrim(string(newcalId),2)+' ...'

    if n_elements(outfile) gt 0 then begin
        openw,unit,outfile,/get_lun
    endif else unit=0

    for i=0L,n_elements(ondate)-1 do begin
        str = strtrim(string(newCalId,format='(I)'),2)
        str = str+', '+strtrim(string(ulong64(ondate[i])),2)
        str = str+', '+strtrim(string(waves[i],format='(G)'),2)
        str = str+', '+strtrim(string(kappa[i],format='(G)'),2)
        str = str+', '+strtrim(string(kappaUnc[i],format='(G)'),2)
        q2 = q1+str+')'
        if keyword_set(verbose) then print,string(i+1)+'  '+q2
        if unit gt 0 then printf,unit,str
        ;if insert then query_database,q2,res,nrows
    endfor

    print,''
    if unit gt 0 then begin
        flush,unit
        close,unit
    endif

    return

END

