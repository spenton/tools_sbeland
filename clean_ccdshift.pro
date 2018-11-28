;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Clean a series of parameters from align_spectra2 over the time range
;   to remove clear spikes in the results and to smooth out the fit.
;
; CALLING SEQUENCE:
;   clean_ccdshift, instrumentModeId, infile=infile, outfile=outfile, polyfit=polyfit
;
; INPUT PARAMETERS:
;   instrumentModeId -
;      The instrument mode of interest are currently limited to:
;      41	SIM_A	VIS1
;      43	SIM_A	UV
;      44	SIM_A	IR
;      45	SIM_B	VIS1
;      47	SIM_B	UV
;      31	SIM_A	ESR
;      32	SIM_B	ESR
;
; OPTIONAL INPUT PARAMETERS:
;   infile - 
;      Name of IDL save file containing the list of plans and coefficients
;      from align_spectra2 for the whole mission.  If not specified, will use
;      a default file name.
;   outfile - 
;      Name of IDL save file to write out with the new list of coefficients
;      after the cleanup.  If not specified, will use a default file name.
;   polyfit -
;      If this flag is set, new CCDFIT values will be generated for the list
;      of plans from infile. This will require to get the prism temperature
;      at each of the prism position for each scan (very slow).
;
; RETURNED PARAMETERS:
;   None
;
; OPTIONAL OUTPUT PARAMETERS:
;   None
;
; EXAMPLE:  
;
; COMMON BLOCKS:
;      NONE
;
;
; REVISION HISTORY:
;   Revision: $Id: clean_ccdshift.pro,v 1.2 2018/11/26 16:16:02 spenton Exp $
;-
;*******************************************************************
pro clean_ccdshift, instrumentModeId, infile=infile, outfile=outfile, $
    polyfit=polyfit, coeffs=coeffs, ccdfit=ccdfit

    if n_elements(infile) eq 0 then begin
        case instrumentModeId of
            41: infile='~/SORCE/data/align_spectra2_41_dettemp.sav'
            43: infile='~/SORCE/data/align_spectra2_43all.sav'
            ;44: infile='~/SORCE/data/align_spectra2_44.sav'
            44: infile='~/SORCE/data/align_spectra2_44_dettemp.sav'
            45: infile='~/SORCE/data/align_spectra2_45_dettemp.sav'
            47: infile='~/SORCE/data/align_spectra2_47all.sav'
            48: infile='~/SORCE/data/align_spectra2_48_dettemp.sav'
            else: begin
                print,'Error: invalid instrumentModeId ',instrumentModeId
                return
            end
        endcase
    endif

    if n_elements(outfile) eq 0 then begin
        outfile=strcompress(infile,/remove_all)
        p=strpos(outfile,'.sav')
        if p lt 0 then p=strlen(outfile)
        outfile=strmid(infile,0,p)+'_clean.sav'
    endif

    ; get the data from the save file
    sObj = OBJ_NEW('IDL_Savefile', infile)
    sNames = sObj->Names()

    k=where(strpos(sNames,'COEF') ge 0,count)
    if count eq 0 then begin
        print,'Error: no COEF variable in savefile'
        return
    endif
    sObj->Restore, sNames[k[0]]
    coef = scope_varfetch(sNames[k[0]], level=0)

    k=where(strpos(sNames,'PL') ge 0,count)
    if count eq 0 then begin
        print,'Error: no PL variable in savefile'
        return
    endif
    sObj->Restore, sNames[k[0]]
    plans = scope_varfetch(sNames[k[0]], level=0)

    k=where(strpos(sNames,'GOOD') ge 0,count)
    if count eq 0 then begin
        print,'Error: no GOOD variable in savefile'
        return
    endif
    sObj->Restore, sNames[k[0]]
    goodness = scope_varfetch(sNames[k[0]], level=0)

    k=where(strpos(sNames,'INTG_PERIOD') ge 0,count)
    if count eq 0 then begin
        print,'Error: no INTG_PERIOD variable in savefile'
        return
    endif
    sObj->Restore, sNames[k[0]]
    intg_period = scope_varfetch(sNames[k[0]], level=0)

    good=where(goodness ge 0.97,count)
    if count lt 10 then good=where(goodness ge 0.90,count)
    if count lt 10 then begin
        print,'Error: not enough GOOD data to process'
        return
    endif

    new_coef=coef

    ; the cleaning process involves a broad gauss_smoothing
    ; then a resistant_mean to get rid of outliers of more than 3 sigma
    ; then a new gauss_smooth with only remaining valid points

    ; We'll do one coefficient at a time
    ; for each segment between OBC events

    obctimes = [0.0d, 449.5, 1570.0, 2173, 2456, 2806, 2838, 2894, 2929, 3028, 3035, 3152, 3569, 3826, 5000]
    for cpos=0, n_elements(coef[0,*])-1 do begin
        if cpos eq 0 then sigma1=4.0 else sigma1=6.0
        for obcpos=0,n_elements(obctimes)-2 do begin
           tpos = where(plans[good].starttime ge obctimes[obcpos] and plans[good].stoptime le obctimes[obcpos+1], count)
           if count lt 10 then continue
           resistant_mean,coef[good[tpos],cpos],sigma1,mean,good=keep
           keep=(good[tpos])[keep]
           ; do not filter the first part since we have clear intg_period differences
           if obcpos gt 0 then begin
               if n_elements(keep) gt 10 then sigma2=1.0 else sigma2=0.1
               temp = gauss_smooth(coef[keep,cpos], sigma2, /edge_truncate)
               resistant_mean,temp - coef[keep,cpos],5.0,mean,good=keep1
               keep=keep[keep1]
               temp = gauss_smooth(coef[keep,cpos], sigma2, /edge_truncate)
           endif else begin
               temp = coef[keep,cpos]
           endelse
           p=where(temp ne 0d)
           tpos_all = where(plans.starttime ge obctimes[obcpos] and plans.stoptime le obctimes[obcpos+1], count)
           new_coef[tpos_all,cpos] = interpol(temp[p], plans[keep[p]].starttime, plans[tpos_all].starttime)
        endfor
        k=where(plans.starttime ge obctimes[1])
        k0=where(new_coef[*,cpos] ne 0d and plans.starttime ge obctimes[1])
        resistant_mean, new_coef[k0,cpos],5.0,mean,good=k1
        new_coef[k,cpos] = interpol(new_coef[k0[k1],cpos], plans[k0[k1]].starttime, plans[k].starttime)
    endfor

    coeffs=temporary(new_coef)
    ccdfit=dblarr(n_elements(plans),4)

    if keyword_set(polyfit) then begin
        ; we need to get the CCDpos and prismtemp for each of the corresponding plan
        jstmt = fjava_get_jdbc_statement(user='sorce', password='sorcedb', server='sorce-db', database='SORCE_SIM_V19')
        dbConn = OBJ_NEW("oJava_DbExchange", jstmt)
        q1a="SELECT microsecondsSinceGpsEpoch, temperature FROM dbo.PrismDriveTemperature where "
        q1a=q1a+" instrumentModeId="+strtrim(string(instrumentModeId),2)+" AND "
        q1a=q1a+" version=19 AND "

        q2a="SELECT microsecondsSinceGpsEpoch, prismPosition FROM SimConvertedDataNumbers WHERE "
        q2a=q2a+" instrumentModeId="+strtrim(string(instrumentModeId),2)+" AND "
        q2a=q2a+" version=19 AND "

        for i=0L, n_elements(plans)-1L do begin
            print,plans[i].starttime,format='($,F10.4)'
            t0=sd2gps(plans[i].starttime)*1d6
            ; add 25 minutes to start of scan
            t1=t0+25.d*60d6
            q3=" microsecondsSinceGpsEpoch>="+strtrim(ulong64(t0),2)+" AND microsecondsSinceGpsEpoch<="+strtrim(ulong64(t1),2)
            q3=q3+" ORDER by microsecondsSinceGpsEpoch ASC"
            temp_data=dbConn->getAllValues(q1a+q3)
            ccd_data=dbConn->getAllValues(q2a+q3)
            if size(temp_data,/tname) ne 'STRUCT' or size(ccd_data,/tname) ne 'STRUCT' then continue
            mn1=min(temp_data.(0),max=mx1)
            mn2=min(ccd_data.(0),max=mx2)
            if mn2 lt mn1 or mx2 gt mx1 then $
                prismtemp=interpol(temp_data.temperature,temp_data.microsecondsSinceGpsEpoch, ccd_data.microsecondsSinceGpsEpoch) $
            else $
                prismtemp=interpol(temp_data.temperature,temp_data.microsecondsSinceGpsEpoch, ccd_data.microsecondsSinceGpsEpoch,/spline) 

            ; get the new wavelengths with updated GammaZ and PixSize
            new_wave = ccd2lambda(instrumentModeId, ccd_data.prismposition, prismtemp, gamz=coeffs[i,0], pixsize=coeffs[i,1]) 
            ; get the corresponding prismposition with the default GammaZ and PixSize for the new wavelength
            new_ccdpos = lambda2ccd(instrumentModeId, new_wave, prismtemp)
            delta_ccdpos = new_ccdpos - ccd_data.prismposition
            ccdfit[i,*] = robust_poly_fit(ccd_data.prismposition, delta_ccdpos, 3, /double)
        endfor
        save,plans,ccdfit,coeffs,goodness, intg_period, file=outfile

    endif else begin

        save,plans,coeffs,goodness, intg_period, file=outfile

    endelse

    return

 end

