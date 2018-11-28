;+
; NAME:   INSERT_PROFILEINTEGRAL
;
; PURPOSE: 
;    This routine populates the SORCE/SimProfileIntegralCal database table with the 
;    changes in wavelength after aligning with the Chance-Kurucz solar spectra.
;
; CALLING SEQUENCE:
;    insert_profileintegral, 41, /dbinsert
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
;   Revision: $Id: insert_profileintegral.pro,v 1.2 2018/11/26 16:16:02 spenton Exp $
;-
;*****************************************************************************
PRO update_calibMetadata_profile, instMode, calId, effDate, version=version, $
    verbose=verbose, dbinsert=insert, comment=comment

    if n_elements(comment) eq 0 then comment='No comment provided'
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
            q1=q1+strtrim(string(CalId),2)+", 'SimProfileIntegralCal', "
            q1=q1+strtrim(string(instMode),2)+", "+strtrim(string(version),2)
            q1=q1+", '"+date_str+"', '"+effDate+"'"
            q1=q1+", '"+comment+"'"
            q1=q1+", 'Stephane Beland'"
            q1=q1+", 'The calibrations used are based on the effectiveDate'"
            q1=q1+", 'Modified from version 4'"
            q1=q1+", 'none'"
            q1=q1+", 'none')"
            if keyword_set(verbose) then print,q1
            if insert then query_database,q1,result
        endif else begin
            q1="UPDATE CalibrationMetadata SET releaseDate='"+date_str+"'"
            q1=q1+", effectiveDate='"+effDate+"'"
            q1=q1+", rationaleForNewVersion='"+comment+"'"
            q1=q1+", version="+strtrim(string(version),2)
            q1=q1+" WHERE calibrationSetId="+strtrim(string(CalId),2)
            if keyword_set(verbose) then print,q1
            if insert then query_database,q1,result
        endelse
    endif

END

;*****************************************************************************
PRO INSERT_PROFILEINTEGRAL, instrumentModeId, dbinsert=dbinsert, verbose=verbose, version=version, comment=comment
    ; insert the Profile Integral with the new wavelengths (y4) and y8 and y13
    ; database table SORCE.dbo.SimProfileIntegralCal

    ; INSERT INTO dbo.SimProfileIntegralCal(calibrationSetId, x, y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, 
    ;    y10, y11, y12, y13, y14, y15, y16, y17, y18, y19, y20) 
    ;    VALUES(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

    modes = [31, 32, 41, 43, 44, 45, 47, 48]
    if n_elements(instrumentModeId) eq 0 then begin
        print,'Error: no InstrumentModeId provided'
        return
    endif

    pos = where(modes eq instrumentModeId, count)
    if count eq 0 then begin
        print,'Error: invalid InstrumentModeId ',instrumentModeId,' [',mode,']'
        return
    endif

    if n_elements(version) eq 0 then version=22
    user='sorce'
    password='sorcedb'
    server='sorce-db'
    database='SORCE'

    ; figure out the calibrationSetId (create one if it doesn't exist)
    q1="SELECT * from CalibrationMetadata where calibrationTableName='SimProfileIntegralCal' "
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
            calId=997
            wlim=[0.0d,3000.0d]
            coeffs=[0.0d, 1.0d, 0.0d, 0.0d]
            ; no diode temperature correction for ESR
            temp_wave=[200d, 3000d]
            temp_corr=[0d, 0d]
            refTemp=24.0d 
            end
        32: begin
            calId=998
            wlim=[0.0d, 530.0, 530.0, 1200.0, 1200.0, 1400.0, 1400.0, 1800.0, 1800.0, 3000.0]
            coeffs=[[-0.12971627d, 0.99944166d,  1.4409120d-06,  6.7865605d-09], $
                    [ 2.05750170d, 1.00001400d, -7.7423468d-06,  8.1799025d-09], $
                    [ 4.88330000d, 1.0,  0.0d,  0.0d], $
                    [-0.56023039d, 1.00005130d,  8.6114986d-07,  3.9931204d-10], $
                    [ 0.0d,        1.0d,         0.0d,           0.0d]]
            ; no diode temperature correction for ESR
            temp_wave=[200d, 3000d]
            temp_corr=[0d, 0d]
            refTemp=24.0d 
            end
        41: begin
            calId=999
            wlim=[0.0d,3000.0d]
            coeffs=[0.0d, 1.0d, 0.0d, 0.0d]
            if version eq 19 then begin
                tempcorr_file='sima_vis_tempcorr_6.txt' 
                readcol,tempcorr_file,temp_wave,temp_corr,refTemp,format='(d,d,d)'
            endif else if version eq 20 then begin
                ; ...optimized is from Alin's process
                ;tempcorr_file='VIS1_idlspect_diodtempcorr_optimized.txt'
                tempcorr_file='~/SORCE/data/sima_vis_tempcorr_20.txt'
                ;tempcorr_file='~/SORCE/data/SORCE_Si_drdt_fit.txt'    ; obtained from lab measurements Feb/2014
                readcol,tempcorr_file,temp_wave,temp_corr,refTemp,format='(d,d,d)'
                ; TEST THE SIGN OF TEMPCORR
                temp_corr*=-1d
                ; the corrections provided by Alin Tolea expect to be applied as:
                ; irrad=irrad_measured / (1 + f*(T-Tref))
                ; The Java code applies the correction as:
                ; irrad=irrad_measured * (1 + f*(Tref-T))
                ; Modify the correction so we can use the Java code with many versions of the SimProfileIntegralCal table
                ;temp_corr = (1d/(1d + temp_corr) -1d)
                
                ; test the TSIS delta_lambda/delta_temp  (new Y18)
                ;y18_file='~/SORCE/data/TSIS_dlambda_dt.txt'
                ;readcol,y18_file,y18_wave, y18_dldt, format='(d,d)'
            endif else if version eq 22 then begin
                tempcorr_file='sima_vis_tempcorr_6.txt' 
                readcol,tempcorr_file,temp_wave,temp_corr,refTemp,format='(d,d,d)'
            endif
            refTemp=median(refTemp)
            end
        43: begin
            calId=1000
            wlim=[0.0d,3000.0d]
            coeffs=[-0.20540544d, 1.0003801d, 2.0269636d-06, 2.1389755d-09]
            ; no diode temperature correction for UV diode
            temp_wave=[200d, 3000d]
            temp_corr=[0d, 0d]
            refTemp=24.0d
            end
        44: begin
            calId=1001
            wlim=[0.0d,3000.0d]
            coeffs=[0.0d, 1.0d, 0.0d, 0.0d]
            if version eq 19 then begin
                tempcorr_file='sima_ir_tempcorr_6.txt' 
                readcol,tempcorr_file,temp_wave,temp_corr,refTemp,format='(d,d,d)'
                refTemp=median(refTemp)
            endif else if version eq 20 then begin
                ; ...optimized is from Alin's process
                ;tempcorr_file='IR_idlspect_diodtempcorr_optimized.txt'
                tempcorr_file='~/SORCE/data/sima_ir_tempcorr_21.txt'
                ;tempcorr_file='~/SORCE/data/SORCE_InGaAs_drdt_fit.txt'    ; obtained from lab measurements Feb/2014
                readcol,tempcorr_file,temp_wave,temp_corr,refTemp,format='(d,d,d)'
                ; the corrections provided by Alin expect to be applied as:
                ; irrad=irrad_measured / (1 + f*(T-Tref))
                ; The Java code applies the correction as:
                ; irrad=irrad_measured * (1 + f*(Tref-T))
                ; Modify the correction so we can use the Java code with many versions of the SimProfileIntegralCal table
                ;temp_corr = (1d/(1d + temp_corr) -1d)
                refTemp=median(refTemp)
            endif else if version eq 22 then begin
                tempcorr_file='sima_ir_tempcorr_6.txt' 
                readcol,tempcorr_file,temp_wave,temp_corr,refTemp,format='(d,d,d)'
                refTemp=median(refTemp)
            endif else if version eq 24 then begin
                calId=2043
            endif
            end
        45: begin
            ; aligned with 41
            calId=1002
            wlim=[0.0d,3000.0d]
            coeffs=[-0.20124423d, 0.99988954d, 4.1470278d-07,  1.9260959d-09]
            if version eq 19 then begin
                tempcorr_file='sima_vis_tempcorr_6.txt' 
                readcol,tempcorr_file,temp_wave,temp_corr,refTemp,format='(d,d,d)'
                refTemp=median(refTemp)
            endif else if version eq 20 then begin
                ; ...optimized is from Alin's process
                ;tempcorr_file='VIS1_idlspect_diodtempcorr_optimized.txt'
                tempcorr_file='~/SORCE/data/sima_vis_tempcorr_20.txt'
                ;tempcorr_file='~/SORCE/data/SORCE_Si_drdt_fit.txt'    ; obtained from lab measurements Feb/2014
                readcol,tempcorr_file,temp_wave,temp_corr,refTemp,format='(d,d,d)'
                ; the corrections provided by Alin Tolea expect to be applied as:
                ; irrad=irrad_measured / (1 + f*(T-Tref))
                ; The Java code applies the correction as:
                ; irrad=irrad_measured * (1 + f*(Tref-T))
                ; Modify the correction so we can use the Java code with many versions of the SimProfileIntegralCal table
                ;temp_corr = (1d/(1d + temp_corr) -1d)
                refTemp=median(refTemp)
            endif else if version eq 22 then begin
                tempcorr_file='~/SORCE/data/sima_vis_tempcorr_6.txt' 
                readcol,tempcorr_file,temp_wave,temp_corr,refTemp,format='(d,d,d)'
                refTemp=median(refTemp)
            endif else if version eq 23 then begin
                delta_wv_file='~/SORCE/data/simb_vis_delta_wavelength.txt' 
                readcol,delta_wv_file, dwavelength, ccdpos, deltawv,format='(d,d,d)'
            endif
            end
        47: begin
            calId=1003
            wlim=[0.0d,3000.0d]
            coeffs=[-0.35628367d, 1.0007715d, 1.9804726d-06, 3.0432922d-09]
            ; no diode temperature correction for UV diode
            temp_wave=[200d, 3000d]
            temp_corr=[0d, 0d]
            refTemp=24.0d
            end
        48: begin
            calId=1081
            wlim=[0.0d,3000.0d]
            coeffs=[0.0d, 1.0d, 0.0d, 0.0d]
            if version eq 19 then begin
                tempcorr_file='sima_ir_tempcorr_6.txt' 
                readcol,tempcorr_file,temp_wave,temp_corr,refTemp,format='(d,d,d)'
                refTemp=median(refTemp)
            endif else if version eq 20 then begin
                ; ...optimized is from Alin's process
                ;tempcorr_file='IR_idlspect_diodtempcorr_optimized.txt'
                tempcorr_file='~/SORCE/data/sima_ir_tempcorr_21.txt'
                ;tempcorr_file='~/SORCE/data/SORCE_InGaAs_drdt_fit.txt'    ; obtained from lab measurements Feb/2014
                readcol,tempcorr_file,temp_wave,temp_corr,refTemp,format='(d,d,d)'
                ; the corrections provided by Alin expect to be applied as:
                ; irrad=irrad_measured / (1 + f*(T-Tref))
                ; The Java code applies the correction as:
                ; irrad=irrad_measured * (1 + f*(Tref-T))
                ; Modify the correction so we can use the Java code with many versions of the SimProfileIntegralCal table
                ;temp_corr = (1d/(1d + temp_corr) -1d)
                refTemp=median(refTemp)
            endif else if version eq 22 then begin
                tempcorr_file='sima_ir_tempcorr_6.txt' 
                readcol,tempcorr_file,temp_wave,temp_corr,refTemp,format='(d,d,d)'
                refTemp=median(refTemp)
            endif
            end
    endcase

    query = 'SELECT * from SimProfileIntegralCal where calibrationSetId='+strtrim(string(calId),2)
    query_database,query, res,nrows

    if instrumentModeId eq 44 and calId eq 1001 then begin
        ; for the IR diode, the last point of y10 is bad in version 1001
        ; use the ts_fcast function to extrapolate additional data points
        ; to get a better fit
        npts=10
        newres=res[-npts:-1]
        newres.x=ts_fcast(res[0:-2].x,200,npts)
        newres.y0=res[0].y0
        newres.y1=res[0].y1
        newres.y2=res[0].y2
        newres.y3=ts_fcast(res[0:-2].y3,200,npts)
        newres.y4=ts_fcast(res[0:-2].y4,200,npts)
        newres.y5=ts_fcast(res[0:-2].y5,200,npts)
        newres.y6=ts_fcast(res[0:-2].y6,200,npts)
        newres.y7=ts_fcast(res[0:-2].y7,200,npts)
        newres.y8=ts_fcast(res[0:-2].y8,200,npts)
        newres.y9=ts_fcast(res[0:-2].y9,200,npts)
        newres.y10=ts_fcast(res[0:-2].y10,200,npts)
        newres.y11=ts_fcast(res[0:-2].y11,200,npts)
        newres.y12=ts_fcast(res[0:-2].y12,200,npts)
        newres.y13=ts_fcast(res[0:-2].y13,200,npts)
        newres.y14=ts_fcast(res[0:-2].y14,200,npts)
        newres.y15=ts_fcast(res[0:-2].y15,200,npts)
        newres.y16=res[0].y16
        newres.y17=res[0].y17
        newres.y18=ts_fcast(res[0:-2].y18,200,npts)
        newres.y19=ts_fcast(res[0:-2].y19,200,npts)
        newres.y20=ts_fcast(res[0:-2].y20,200,npts)
        res=[res[0:-2],newres]
    endif

    ; correct the wavelengths
    if version eq 19 then begin
        ny4=res.y4
        for i =0,n_elements(wlim)/2 -1 do begin
            p=where(res.y4 ge wlim[i*2] and res.y4 le wlim[i*2+1], count)
            if count eq 0 then begin
                print,'Error: no wavelength found between ',wlim[i*2:i*2+1]
                return
            endif
            ny4[p] = poly(res[p].y4,coeffs[*,i])
        endfor
    endif else begin
        ; use the updated instrument model parameters (aligned to Kurucz) to get the wavelengths
        ny4 = ccd2lambda(instrumentModeId, res.x, res.y0) 
    endelse

    if version eq 23 and instrumentModeId eq 45 then begin
        ; for version=23 and mode 45, we simply update the corresponding wavelength at each ccdpos
        ;restore,'~/SORCE/data/simb_vis_profileIntegral_id2037_no_deltaw.sav'
        ;res=prof23b
        ; to force alignment with the spectra of SimA on day 453.67
        ; the original correction turned out to be too much after processing
        ; try half
        ;res.y4=res.y4 + interpol(deltawv, dwavelength, res.y4)
        ; SBeland 20160527 -> to get Kappa to stop increasing after 700nm, we arbitrarily adjust the
        ; refTemp of SIMB for the diodeTempCorr (3.0C seems like the right value)
        ;res.y2 -= 2.5d
        ;
        ; this latest file (SBeland 20160531 contains a better alignment on day 32, refTemp-=2.5
        ; and a better y8 so VISA and VISB are the same on day 32
        restore,'~/SORCE/data/simb_vis_profileIntegral_id2037_refTemp_newWave.sav'
        res=prof23b
    endif else if version eq 23 and instrumentModeId eq 47 then begin
        ;restore,'~/SORCE/data/simb_uv_profile_20151013.sav'
        ;restore,'~/SORCE/data/simb_uv_profile_20151021.sav'
        ;restore,'~/SORCE/data/simb_uv_profile_20151026.sav'
        ;restore,'~/SORCE/data/simb_uv_profile_20151102.sav'
        restore,'~/SORCE/data/simb_uv_profile_20160606.sav'
        res=new_profuvb
    endif else if version eq 23 and instrumentModeId eq 48 then begin
        ;restore,'~/SORCE/data/simb_ir_profile_20151030.sav'
        ;restore,'~/SORCE/data/simb_ir_profile_20151101.sav'
        restore,'~/SORCE/data/simb_ir_profile_20151102.sav'
        res=new_profirb
    endif else if version eq 24 and instrumentModeId eq 44 then begin
        ;restore,'~/SORCE/data/sima_ir_profile_20180110.sav'
        restore,'~/SORCE/data/sima_ir_profile_20180124.sav'
        res=new_prof44
        ; update the detector temperature correction
        ;tempcorr_file='~/SORCE/data/sima_ir_tempcorr_v24_7.txt' 
        ;tempcorr_file='~/SORCE/data/sima_ir_tempcorr_24.txt' 
        tempcorr_file='~/SORCE/data/sima_ir_tempcorr_24_8.txt' 
        readcol,tempcorr_file,temp_wave,temp_corr,refTemp,format='(d,d,d)'
        res.y12 = interpol(temp_corr, temp_wave,res.y4,/lsq)
        res.y2 = median(reftemp)
    endif else begin
        ; now correct each column by interpolating
        ; we're correcting: y1, y3, y5, y6, y7, y8, y9, y10, y11, y12, y13, y14, y15, y18, y19, y20
        temp = interpol(res.y1, res.y4, ny4, /lsq)
        res.y1 = temp
        temp = interpol(res.y3, res.y4, ny4, /lsq)
        res.y3 = temp
        temp = interpol(res.y5, res.y4, ny4, /lsq)
        res.y5 = temp
        temp = interpol(res.y6, res.y4, ny4, /lsq)
        res.y6 = temp
        temp = interpol(res.y7, res.y4, ny4, /lsq)
        res.y7 = temp
        temp = interpol(res.y8, res.y4, ny4, /lsq)
        res.y8 = temp
        temp = interpol(res.y9, res.y4, ny4, /lsq)
        res.y9 = temp
        temp = interpol(res.y10, res.y4, ny4, /lsq)
        res.y10 = temp
        temp = interpol(res.y11, res.y4, ny4, /lsq)
        res.y11 = temp
        ; y12  profIntegralDetTCoef 
        if n_elements(temp_corr) eq 0 then begin
            temp = interpol(res.y12, res.y4, ny4, /lsq)
            res.y12 = temp
        endif else begin
            temp = interpol(temp_corr, temp_wave, ny4, /lsq)
            res.y12 = temp
            res.y2 = res.y2 * 0d + refTemp
        endelse
        temp = interpol(res.y13, res.y4, ny4, /lsq)
        res.y13 = temp
        temp = interpol(res.y14, res.y4, ny4, /lsq)
        res.y14 = temp
        temp = interpol(res.y15, res.y4, ny4, /lsq)
        res.y15 = temp
        ; y18 wavelengthRefTCoef
        if instrumentModeId eq 32 then begin
            ; mode 32 has the wavelengthRefTCoef slightly off -> recalculate
            wave1 = ccd2lambda(32, res.x, res.y0-1d) 
            wave0 = ccd2lambda(32, res.x, res.y0) 
            dwdt = (wave1 - wave0)/ (-1d)
            temp = interpol(dwdt, wave0, ny4, /lsq)
            res.y18 = temp
        endif else begin
            temp = interpol(res.y18, res.y4, ny4, /lsq)
            res.y18 = temp
        endelse
        if n_elements(y18_wave) gt 0 then begin
            temp = interpol(y18_dldt, y18_wave, ny4, /lsq)
            res.y18 = temp
        endif
        temp = interpol(res.y19, res.y4, ny4, /lsq)
        res.y19 = temp
        temp = interpol(res.y20, res.y4, ny4, /lsq)
        res.y20 = temp

        res.y4=ny4
    
    endelse


    if keyword_set(dbinsert) then insert=1 else insert=0

    ; first delete all entries in the table for this particular version (version 1 only for now)
    q1 = "DELETE FROM SimProfileIntegralCal WHERE calibrationSetId="+strtrim(string(newcalId),2)
    if keyword_set(verbose) then print,q1
    if insert then query_database,q1

    ; make sure we have an entry in the CalibrationMetadata table for this calibrationSetId
    effectiveDate='2003-02-01 00:00:01.0'
    update_calibMetadata_profile, instrumentModeId, newCalId, effectiveDate, $
        version=version, dbinsert=insert, verbose=verbose, comment=comment

    ; insert one row at a time
    q1 = "INSERT INTO SimProfileIntegralCal(calibrationSetId, x, y0, y1, y2, y3, y4, y5, "
    q1 += "y6, y7, y8, y9, y10, y11, y12, y13, y14, y15, y16, y17, y18, y19, y20) VALUES("

    print,strtrim(string(n_elements(res)),2)+' rows to insert as calibrationSetId='+strtrim(string(newcalId),2)+' ...'

    for i=0L,n_elements(res)-1 do begin
        str = strtrim(string(newCalId,format='(I)'),2)
        str = str+', '+strtrim(string(res[i].x,format='(I)'),2)
        str = str+', '+strtrim(string(res[i].y0,format='(D0.4)'),2)
        str = str+', '+strtrim(string(res[i].y1),2)
        str = str+', '+strtrim(string(res[i].y2,format='(D0.4)'),2)
        str = str+', '+strtrim(string(res[i].y3),2)
        str = str+', '+strtrim(string(res[i].y4),2)
        str = str+', '+strtrim(string(res[i].y5),2)
        str = str+', '+strtrim(string(res[i].y6),2)
        str = str+', '+strtrim(string(res[i].y7),2)
        str = str+', '+strtrim(string(res[i].y8),2)
        str = str+', '+strtrim(string(res[i].y9),2)
        str = str+', '+strtrim(string(res[i].y10),2)
        str = str+', '+strtrim(string(res[i].y11),2)
        str = str+', '+strtrim(string(res[i].y12),2)
        str = str+', '+strtrim(string(res[i].y13),2)
        str = str+', '+strtrim(string(res[i].y14),2)
        str = str+', '+strtrim(string(res[i].y15),2)
        str = str+', '+strtrim(string(res[i].y16),2)
        str = str+', '+strtrim(string(res[i].y17),2)
        str = str+', '+strtrim(string(res[i].y18),2)
        str = str+', '+strtrim(string(res[i].y19),2)
        str = str+', '+strtrim(string(res[i].y20),2)
        q2 = q1+str+')'
        if keyword_set(verbose) then print,string(i+1)+'  '+q2
        if insert then query_database,q2,result,nrows
    endfor

    return

END

