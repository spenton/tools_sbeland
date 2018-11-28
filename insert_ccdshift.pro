;+
; NAME:   INSERT_CCDSHIFT
;
; PURPOSE: 
;    This routine populates the SORCE/SimCcdShiftCal database table with the 
;    measured CCD shifts measured after each OBC anomaly event.
;    The data corresponds to instrumentModeID 54 (simA) or 55 (simB).
;
; CALLING SEQUENCE:
;    insert_ccdshift, /dbinsert
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
; REVISION HISTORY:
;   Revision: $Id: insert_ccdshift.pro,v 1.2 2018/11/26 16:16:02 spenton Exp $
;-
;*****************************************************************************
PRO update_calibMetadata_ccdshift, instMode, calId, effDate, version=version, $
    verbose=verbose, dbinsert=insert, filename=filename,comment=comment

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
            q1=q1+strtrim(string(CalId),2)+", 'SimCcdShiftCal', "
            q1=q1+strtrim(string(instMode),2)+", "+strtrim(string(version),2)
            q1=q1+", '"+date_str+"', '"+effDate+"'"
            ;q1=q1+", 'release for version 23'"
            q1=q1+", 'test release for version 24'"
            q1=q1+", 'Stephane Beland'"
            if n_elements(comment) gt 0 then q1=q1+", '"+comment+"'" else $
                q1=q1+", 'Ccd shift correction using get_ccdoffset, get_ccdshift and get_ccdshift_hybrid for IR and updated version of SimProfileIntegralCal 2138'"
            q1=q1+", '"+filename+"'"
            q1=q1+", 'none'"
            q1=q1+", 'the integration time is ignored for version 24')"
            if keyword_set(verbose) then print,q1
            if insert then query_database,q1,result
        endif else begin
            q1="UPDATE CalibrationMetadata SET releaseDate='"+date_str+"'"
            q1=q1+", effectiveDate='"+effDate+"'"
            q1=q1+", version="+strtrim(string(version),2)
            if n_elements(comment) gt 0 then q1=q1+", dataDescription='"+comment+"'" 
            q1=q1+", sourceOfData='"+filename+"'"
            q1=q1+" WHERE calibrationSetId="+strtrim(string(CalId),2)
            if keyword_set(verbose) then print,q1
            if insert then query_database,q1,result
        endelse
    endif

END

;*****************************************************************************
PRO INSERT_CCDSHIFT, instrumentModeId, dbinsert=dbinsert, verbose=verbose, version=version,comment=comment
    ; insert the list of SIM CCD Shifts in the
    ; database table SORCE.dbo.SimCcdShiftCal

    modes = [31, 32, 41, 43, 44, 45, 47, 48]
    pos = where(modes eq instrumentModeId, count)
    if count eq 0 then begin
        print,'Error: invalid InstrumentModeId ',instrumentModeId,' [',mode,']'
        return
    endif

    if n_elements(version) eq 0 then version=2
    user='sorce'
    password='sorcedb'
    server='sorce-db'
    database='SORCE'

    ; figure out the calibrationSetId (create one if it doesn't exist)
    q1="SELECT * from CalibrationMetadata where calibrationTableName='SimCcdShiftCal' "
    q1 += " and instrumentModeId="+strtrim(string(instrumentModeId,format='(I)'),2)
    q1 += " and version="+strtrim(string(version,format='(I)'),2)
    query_database,q1,result,nrows,user=user,password=password,server=server,database=database
    if nrows ge 1 then begin
        newCalId=result.calibrationSetId
    endif else begin
        ; insert a new row in the CalibrationMetaData table
        query_database,'select max(calibrationSetId) as mx from CalibrationMetadata',result
        newCalId = result.mx + 1
    endelse

    case instrumentModeId of
        31: begin
            if version eq 0 then begin
                ;newCalId=1024
                infile='~/SORCE/data/sima_esr_ccdshift.txt' 
                readcol,infile,mode,str0,str1,ondate,intg, c0,c1,c2,c3,delim=",",format='(I,A,A,d,I,d,d,d,d)'
            endif else if version eq 1 then begin
                ;newCalId=1157
                ; we use the VIS offsets which match well with what we see with ESRFullScan
                ;infile='~/SORCE/data/sima_vis_dettemp_clean_v20_ccdshift.txt' 
                ; testing a new file with different offsets fit for each scan (not cleaned) from SD=0 to 450
                ;infile='~/SORCE/data/sima_vis_dettemp_v20_ccdshift.txt' 
                ; latest file obtained with align_spectra2, spl_interp  (unfortunately, the sign got reversed)
                ;infile='~/SORCE/data/sima_vis_dettemp_v21_ccdshift.txt' 
                infile='~/SORCE/data/sima_vis_intg_period_clean_ccdshift.txt' 
                readcol,infile,mode,ondate,intg, c0,c1,c2,c3,format='(I,d,f,d,d,d,d)'
            endif else if version eq 2 then begin
                ;infile='~/SORCE/data/sima_vis_ccdshift_smooth_0_4400.txt' 
                infile='~/SORCE/data/sima_vis_ccdshift_smooth_intg_period_0_4400.txt' 
                readcol,infile,mode,ondate,intg, c0,c1,c2,c3,format='(I,d,f,d,d,d,d)'
            endif else if version eq 3 then begin
                infile='~/SORCE/data/sima_esr_ccdshift_v23.txt' 
                readcol,infile,mode,ondate,intg, c0,c1,c2,c3,format='(I,d,f,d,d,d,d)'
            endif else if version eq 4 then begin
                infile='~/SORCE/data/sima_esr_ccdoffset_2053.txt' 
                readcol,infile,mode,ondate,intg, c0,c1,c2,c3,format='(I,d,f,d,d,d,d)'
            endif else if version eq 25 then begin
                infile='~/SORCE/data/sima_vis_ccdshift_20181014.txt' 
                readcol,infile,mode,ondate,intg, c0,c1,c2,c3,format='(I,d,f,d,d,d,d)'
            endif
            end
        32: begin
            if version eq 0 then begin
                ;newCalId=1025
                infile='~/SORCE/data/simb_esr_ccdshift.txt' 
                readcol,infile,mode,str0,str1,ondate,intg, c0,c1,c2,c3,delim=",",format='(I,A,A,d,I,d,d,d,d)'
            endif else if version eq 1 then begin
                ;newCalId=1158
                ; we use the VIS offsets which match well with what we see with ESRFullScan
                ;infile='~/SORCE/data/simb_vis_dettemp_clean_v20_ccdshift.txt' 
                ; latest file obtained with align_spectra2, spl_interp  (unfortunately, the sign got reversed)
                ;infile='~/SORCE/data/simb_vis_dettemp_v21_ccdshift.txt' 
                infile='~/SORCE/data/simb_vis_v20_2011_clean_ccdshift.txt' 
                readcol,infile,mode,ondate,intg, c0,c1,c2,c3,format='(I,d,f,d,d,d,d)'
            endif else if version eq 2 then begin
                infile='~/SORCE/data/simb_vis_ccdshift_smooth_0_4400.txt' 
                readcol,infile,mode,ondate,intg, c0,c1,c2,c3,format='(I,d,f,d,d,d,d)'
            endif else if version eq 3 then begin
                infile='~/SORCE/data/simb_vis_ccdshift_v23.txt' 
                readcol,infile,mode,ondate,intg, c0,c1,c2,c3,format='(I,d,f,d,d,d,d)'
            endif else if version eq 25 then begin
                infile='~/SORCE/data/sima_vis_ccdshift_20181014.txt' 
                readcol,infile,mode,ondate,intg, c0,c1,c2,c3,format='(I,d,f,d,d,d,d)'
            endif
            end
        41: begin
            if version eq 0 then begin
                ;newCalId=1026
                infile='~/SORCE/data/sima_vis_ccdshift.txt' 
                readcol,infile,mode,str0,str1,ondate,intg, c0,c1,c2,c3,delim=",",format='(I,A,A,d,I,d,d,d,d)'
            endif else if version eq 1 then begin
                ;newCalId=1153
                ;infile='~/SORCE/data/sima_vis_dettemp_clean_v20_ccdshift.txt' 
                ; we tested varying the prism temp for possible better fit and trend over large temp changes
                ;infile='~/SORCE/data/sima_vis_dettemp_prismtemp_v20_ccdshift.txt' 
                ; testing a new file with different offsets fit for each scan (not cleaned) from SD=0 to 450
                ;infile='~/SORCE/data/sima_vis_dettemp_v20_ccdshift.txt' 
                ; latest file obtained with align_spectra2, spl_interp  (unfortunately, the sign got reversed)
                ;infile='~/SORCE/data/sima_vis_dettemp_v21_ccdshift.txt' 
                ;infile='~/SORCE/data/sima_vis_v20_2011_clean_ccdshift.txt' 
                infile='~/SORCE/data/sima_vis_intg_period_clean_ccdshift.txt' 
                readcol,infile,mode,ondate,intg, c0,c1,c2,c3,format='(I,d,f,d,d,d,d)'
            endif else if version eq 2 then begin
                ;infile='~/SORCE/data/sima_vis_ccdshift_smooth_0_4400.txt' 
                ;infile='~/SORCE/data/sima_vis_ccdshift_smooth_intg_period_0_4400.txt' 
                infile='~/SORCE/data/sima_vis_ccdshift_smooth_intg_period_0_4400_fixed.txt' 
                readcol,infile,mode,ondate,intg, c0,c1,c2,c3,format='(I,d,f,d,d,d,d)'
            endif else if version eq 3 then begin
                infile='~/SORCE/data/sima_vis_ccdshift_v23.txt' 
                readcol,infile,mode,ondate,intg, c0,c1,c2,c3,format='(I,d,f,d,d,d,d)'
            endif else if version eq 4 then begin
                ;infile='~/SORCE/data/sima_vis_ccdoffset_20150905.txt' 
                ;infile='~/SORCE/data/sima_vis_ccdoffset_20150908.txt' 
                infile='~/SORCE/data/sima_vis_ccdoffset_2053.txt' 
                readcol,infile,mode,ondate,intg, c0,c1,c2,c3,format='(I,d,f,d,d,d,d)'
            endif else if version eq 25 then begin
                infile='~/SORCE/data/sima_vis_ccdshift_20181014.txt' 
                readcol,infile,mode,ondate,intg, c0,c1,c2,c3,format='(I,d,f,d,d,d,d)'
            endif
            end
        43: begin
            if version eq 0 then begin
                ;newCalId=1027
                infile='~/SORCE/data/sima_uv_ccdshift.txt' 
                readcol,infile,mode,str0,str1,ondate,intg, c0,c1,c2,c3,delim=",",format='(I,A,A,d,I,d,d,d,d)'
            endif else if version eq 1 then begin
                ;newCalId=1151
                ;infile='~/SORCE/data/sima_uv_v20_clean_ccdshift.txt' 
                ;infile='~/SORCE/data/sima_uv_v21_ccdshift.txt' 
                ;infile='~/SORCE/data/sima_uv_v20_2011_clean_ccdshift.txt' 
                infile='~/SORCE/data/sima_uv_intg_period_clean_453_ccdshift.txt' 
                readcol,infile,mode,ondate,intg, c0,c1,c2,c3,format='(I,d,f,d,d,d,d)'
            endif else if version eq 2 then begin
                ;infile='~/SORCE/data/sima_uv_ccdshift_smooth_0_4400.txt' 
                infile='~/SORCE/data/sima_uv_ccdshift_smooth_intg_period_0_4400.txt' 
                readcol,infile,mode,ondate,intg, c0,c1,c2,c3,format='(I,d,f,d,d,d,d)'
            endif else if version eq 3 then begin
                infile='~/SORCE/data/sima_uv_ccdshift_v23.txt' 
                readcol,infile,mode,ondate,intg, c0,c1,c2,c3,format='(I,d,f,d,d,d,d)'
            endif else if version eq 4 then begin
                ;infile='~/SORCE/data/sima_uv_ccdshift_20151015.txt' 
                infile='~/SORCE/data/sima_uv_ccdshift_20151020.txt' 
                readcol,infile,mode,ondate,intg, c0,c1,c2,c3,format='(I,d,f,d,d,d,d)'
            endif else if version eq 25 then begin
                infile='~/SORCE/data/sima_uv_ccdshift_20181014.txt' 
                readcol,infile,mode,ondate,intg, c0,c1,c2,c3,format='(I,d,f,d,d,d,d)'
            endif
            end
        44: begin
            if version eq 0 then begin
                ;newCalId=1028
                infile='~/SORCE/data/sima_ir_ccdshift.txt' 
                readcol,infile,mode,str0,str1,ondate,intg, c0,c1,c2,c3,delim=",",format='(I,A,A,d,I,d,d,d,d)'
            endif else if version eq 1 then begin
                ;newCalId=1155
                ;infile='~/SORCE/data/sima_ir_dettemp_clean_v20_ccdshift.txt'
                ;infile='~/SORCE/data/sima_ir_v21_ccdshift.txt'
                ;infile='~/SORCE/data/sima_ir_v20_2011_clean_ccdshift.txt'
                ;infile='~/SORCE/data/sima_ir_intg_period_clean_453_ccdshift.txt'
                infile='~/SORCE/data/sima_ir_intg_period_all_453_ccdshift.txt'
                readcol,infile,mode,ondate,intg, c0,c1,c2,c3,format='(I,d,f,d,d,d,d)'
            endif else if version eq 2 then begin
                ;infile='~/SORCE/data/sima_ir_ccdshift_smooth_0_4400.txt' 
                ;infile='~/SORCE/data/sima_ir_ccdshift_smooth_intg_period_0_4400.txt' 
                infile='~/SORCE/data/sima_ir_ccdshift_smooth_v22_0_4400.txt' 
                readcol,infile,mode,ondate,intg, c0,c1,c2,c3,format='(I,d,f,d,d,d,d)'
            endif else if version eq 3 then begin
                infile='~/SORCE/data/sima_ir_ccdshift_v23.txt' 
                readcol,infile,mode,ondate,intg, c0,c1,c2,c3,format='(I,d,f,d,d,d,d)'
            endif else if version eq 4 then begin
                infile='~/SORCE/data/sima_ir_ccdshift_20151031.txt' 
                readcol,infile,mode,ondate,intg, c0,c1,c2,c3,format='(I,d,f,d,d,d,d)'
            endif else if version eq 11 then begin
                infile='~/SORCE/data/sima_ir_ccdshift_v24_20180131.txt' 
                readcol,infile,mode,ondate,intg, c0,c1,c2,c3,format='(I,d,f,d,d,d,d)'
            endif else if version eq 12 or version eq 7 then begin
                infile='~/SORCE/data/sima_ir_ccdshift_20180203.txt' 
                readcol,infile,mode,ondate,intg, c0,c1,c2,c3,format='(I,d,f,d,d,d,d)'
            endif else if version eq 25 then begin
                ;infile='~/SORCE/data/sima_ir_ccdshift_20181014.txt' 
                infile='~/SORCE/data/sima_ir_ccdshift_20181107.txt' 
                readcol,infile,mode,ondate,intg, c0,c1,c2,c3,format='(I,d,f,d,d,d,d)'
            endif else if version eq 125 then begin
                infile='~/SORCE/data/sima_ir_ccdshift_20181120.txt' 
                readcol,infile,mode,ondate,intg, c0,c1,c2,c3,format='(I,d,f,d,d,d,d)'
            endif
            end
        45: begin
            if version eq 0 then begin
                ;newCalId=1029
                infile='~/SORCE/data/simb_vis_ccdshift.txt' 
                readcol,infile,mode,str0,str1,ondate,intg, c0,c1,c2,c3,delim=",",format='(I,A,A,d,I,d,d,d,d)'
            endif else if version eq 1 then begin
                ;newCalId=1154
                ;infile='~/SORCE/data/simb_vis_dettemp_clean_v20_ccdshift.txt'
                ; testing a new file with different offsets fit for each scan (not cleaned) from SD=0 to 450
                ;infile='~/SORCE/data/simb_vis_dettemp_v20_ccdshift.txt' 
                ; latest file obtained with align_spectra2, spl_interp  (unfortunately, the sign got reversed)
                ;infile='~/SORCE/data/simb_vis_dettemp_v21_ccdshift.txt' 
                infile='~/SORCE/data/simb_vis_v20_2011_clean_ccdshift.txt' 
                readcol,infile,mode,ondate,intg, c0,c1,c2,c3,format='(I,d,f,d,d,d,d)'
            endif else if version eq 2 then begin
                infile='~/SORCE/data/simb_vis_ccdshift_smooth_0_4400.txt' 
                readcol,infile,mode,ondate,intg, c0,c1,c2,c3,format='(I,d,f,d,d,d,d)'
            endif else if version eq 3 then begin
                infile='~/SORCE/data/simb_vis_ccdshift_v23.txt' 
                readcol,infile,mode,ondate,intg, c0,c1,c2,c3,format='(I,d,f,d,d,d,d)'
            endif else if version eq 4 then begin
                infile='~/SORCE/data/simb_vis_ccdoffset_20150916.txt' 
                readcol,infile,mode,ondate,intg, c0,c1,c2,c3,format='(I,d,f,d,d,d,d)'
            endif else if version eq 25 then begin
                infile='~/SORCE/data/simb_vis_ccdshift_20181014.txt' 
                readcol,infile,mode,ondate,intg, c0,c1,c2,c3,format='(I,d,f,d,d,d,d)'
            endif
            end
        47: begin
            if version eq 0 then begin
                ;newCalId=1030
                infile='~/SORCE/data/simb_uv_ccdshift.txt' 
                readcol,infile,mode,str0,str1,ondate,intg, c0,c1,c2,c3,delim=",",format='(I,A,A,d,I,d,d,d,d)'
            endif else if version eq 1 then begin
                ;newCalId=1152
                ;infile='~/SORCE/data/simb_uv_v20_clean_ccdshift.txt'
                ;infile='~/SORCE/data/simb_uv_v21_ccdshift.txt' 
                infile='~/SORCE/data/simb_uv_v20_2011_clean_ccdshift.txt' 
                readcol,infile,mode,ondate,intg, c0,c1,c2,c3,format='(I,d,f,d,d,d,d)'
            endif else if version eq 2 then begin
                infile='~/SORCE/data/simb_uv_ccdshift_smooth_0_4400.txt' 
                readcol,infile,mode,ondate,intg, c0,c1,c2,c3,format='(I,d,f,d,d,d,d)'
            endif else if version eq 3 then begin
                infile='~/SORCE/data/simb_uv_ccdshift_v23.txt' 
                readcol,infile,mode,ondate,intg, c0,c1,c2,c3,format='(I,d,f,d,d,d,d)'
            endif else if version eq 4 then begin
                ;infile='~/SORCE/data/simb_uv_ccdshift_20151015.txt' 
                ;infile='~/SORCE/data/simb_uv_ccdshift_20151019.txt' 
                ;infile='~/SORCE/data/simb_uv_ccdshift_20151027.txt' 
                infile='~/SORCE/data/simb_uv_ccdshift_20151102.txt' 
                readcol,infile,mode,ondate,intg, c0,c1,c2,c3,format='(I,d,f,d,d,d,d)'
            endif else if version eq 25 then begin
                infile='~/SORCE/data/simb_uv_ccdshift_20181014.txt' 
                readcol,infile,mode,ondate,intg, c0,c1,c2,c3,format='(I,d,f,d,d,d,d)'
            endif
            end
        48: begin
            if version eq 0 then begin
                ;newCalId=1082
                infile='~/SORCE/data/simb_ir_ccdshift.txt' 
                readcol,infile,mode,str0,str1,ondate,intg, c0,c1,c2,c3,delim=",",format='(I,A,A,d,I,d,d,d,d)'
            endif else if version eq 1 then begin
                ;newCalId=1156
                ;infile='~/SORCE/data/simb_ir_dettemp_clean_v20_ccdshift.txt'
                ;infile='~/SORCE/data/simb_ir_v21_ccdshift.txt'
                infile='~/SORCE/data/simb_ir_v20_2011_clean_ccdshift.txt'
                readcol,infile,mode,ondate,intg, c0,c1,c2,c3,format='(I,d,f,d,d,d,d)'
            endif else if version eq 2 then begin
                infile='~/SORCE/data/simb_ir_ccdshift_smooth_0_4400.txt' 
                readcol,infile,mode,ondate,intg, c0,c1,c2,c3,format='(I,d,f,d,d,d,d)'
            endif else if version eq 3 then begin
                infile='~/SORCE/data/simb_ir_ccdshift_v23.txt' 
                readcol,infile,mode,ondate,intg, c0,c1,c2,c3,format='(I,d,f,d,d,d,d)'
            endif else if version eq 4 then begin
                infile='~/SORCE/data/simb_ir_ccdshift_20151101.txt' 
                readcol,infile,mode,ondate,intg, c0,c1,c2,c3,format='(I,d,f,d,d,d,d)'
            endif else if version eq 25 then begin
                infile='~/SORCE/data/simb_ir_ccdshift_20181014.txt' 
                readcol,infile,mode,ondate,intg, c0,c1,c2,c3,format='(I,d,f,d,d,d,d)'
            endif
            end
    endcase

    if keyword_set(dbinsert) then insert=1 else insert=0

    ; first delete all entries in the table for this particular version
    q1 = "DELETE FROM SimCcdShiftCal WHERE calibrationSetId="+strtrim(string(newcalId),2)
    if keyword_set(verbose) then print,q1
    if insert then query_database,q1

    ; make sure we have an entry in the CalibrationMetadata table for this calibrationSetId
    effectiveDate='2003-01-24 00:00:01.0'
    update_calibMetadata_ccdshift, instrumentModeId, newCalId, effectiveDate, comment=comment, $
        version=version, dbinsert=insert, verbose=verbose, filename=file_basename(infile)

    ; insert one row at a time
    q1 = "INSERT INTO SimCcdShiftCal(calibrationSetId, effectiveDateGps, ccdIntgPeriod, c0, c1, c2, c3) VALUES("

    print,strtrim(string(n_elements(ondate)),2)+' rows to insert as calibrationSetId='+strtrim(string(newcalId),2)+' ...'

    for i=0L,n_elements(ondate)-1 do begin
        str = strtrim(string(newCalId,format='(I)'),2)
        str = str+', '+strtrim(string(ulong64(ondate[i])),2)
        str = str+', '+strtrim(string(intg[i],format='(I)'),2)
        str = str+', '+strtrim(string(c0[i],format='(G)'),2)
        str = str+', '+strtrim(string(c1[i],format='(G)'),2)
        str = str+', '+strtrim(string(c2[i],format='(G)'),2)
        str = str+', '+strtrim(string(c3[i],format='(G)'),2)
        q2 = q1+str+')'
        if keyword_set(verbose) then print,string(i+1)+'  '+q2
        if insert then query_database,q2,res,nrows
    endfor

    return

END

