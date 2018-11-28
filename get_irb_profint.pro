;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Generate a profile integral table for the SIMB IR detector based on
;   the SimA IR data.  The wavelength scale is obtained from scans of IRA  
;   and IRB on day 453. The profileIntegral is obtained by matching the 
;   irradiance from IRA and IRB.
;
; CALLING SEQUENCE:
;
; INPUT PARAMETERS:
;
; OPTIONAL INPUT PARAMETERS:
;
; RETURNED PARAMETERS:
;
; OPTIONAL OUTPUT PARAMETERS:
;
; EXAMPLE:  
;
; COMMON BLOCKS:
;      NONE
;
; NOTES:
;
; REVISION HISTORY:
;   Revision: $Id: get_irb_profint.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;------------------------------------------------------------------
function get_irb_profint, noplot=noplot, dbinsert=dbinsert, verbose=verbose


    refT0=453.67840d
    refT1=453.69524d
    newCalId = 1081
    version = 19
    instrumentModeId = 48

    ; we'll ignore the wavelength shift with temperature because both SimA and SimB are very stable
    ; for this scan and the difference between the 2 prisms is smaller than 0.04 degrees.

    ; get the spectra for mode 44 (SimA IR)
    spa=get_sim_spectra(reft0,reft1,/mission,44,/notempcorr,/noccdcorr,/no1au,/dn_only,/verbose)

    ; get the spectra for mode 48 (SimB IR)
    spb=get_sim_spectra(reft0,reft1,/mission,48,/notempcorr,/noccdcorr,/no1au,/dn_only,/verbose)

    ; re-arrange the ccd position of SimB to roughlu align with SimA
    spb.ccdpos_uncorr = max(spb.ccdpos)-spb.ccdpos+16805d

    ; align IRB to IRA spectrum
    res_align=align_spectra(0,1,44,refspect=spa,spect=spb,coeffs=coeffs,order=3)

    if NOT keyword_set(noplot) then begin
        plot_multi,spa.wavelength,spa.dn,res_align.wavelength,res_align.dn,xtitle='Wavelength (nm)',$
            ytitle='DN',label=['SPA','SPB Aligned'],/xst,/yst
        plot_multi,spa.ccdpos_uncorr,spa.wavelength,res_align.ccdpos_uncorr,res_align.wavelength,xtitle='CCD Position',$
            ytitle='Wavelength (nm)',label=['SPA','SPB Aligned'],/xst,/yst 
    endif

    ; fit the original CCD position vs wavelengths with a bspline
    sset=bspline_iterfit(res_align.ccdpos,res_align.wavelength,maxiter=0,requiren=10,bkspace=5)
    mn=min(spb.ccdpos,max=mx)
    mn = floor(mn / 10.) * 10.
    mx = ceil(mx / 10.) * 10.
    ccdpos = findgen((mx - mn)/10. +1)*10. + mn
    wavelength = bspline_valu(ccdpos,sset)

    ; now scale the profile integral from SimA to match the DNs of SimB
    query_database,/reset
    query_database, 'SELECT * FROM dbo.SimProfileIntegralCal where calibrationSetId=1013', profInt_A

    profInt_B = replicate(profInt_a[0],n_elements(ccdpos))
    profInt_B.calibrationSetId=newCalId
    profInt_B.x = ccdpos
    profInt_B.y3 = interpol(profInt_A.y3, profInt_A.y4, wavelength)
    profInt_B.y4 = wavelength
    profInt_B.y5 = interpol(profInt_A.y5, profInt_A.y4, wavelength)
    profInt_B.y6 = interpol(profInt_A.y6, profInt_A.y4, wavelength)
    profInt_B.y7 = interpol(profInt_A.y7, profInt_A.y4, wavelength)

    ; we need to scale the profileIntegral of IRB to match the irradiance calculated with IRA
    dna = interpol(spa.dn, spa.wavelength, wavelength, /spline)
    dnb = interpol(res_align.dn, res_align.wavelength, wavelength, /spline)
    y8 = interpol(profInt_A.y8, profInt_A.y4, wavelength, /spline) 
    profInt_B.y8 = y8 * dnb / dna

    profInt_B.y9 = interpol(profInt_A.y9, profInt_A.y4, wavelength)
    profInt_B.y10 = interpol(profInt_A.y10, profInt_A.y4, wavelength)
    profInt_B.y11 = interpol(profInt_A.y11, profInt_A.y4, wavelength)
    profInt_B.y12 = interpol(profInt_A.y12, profInt_A.y4, wavelength)
    profInt_B.y13 = interpol(profInt_A.y13, profInt_A.y4, wavelength)
    profInt_B.y14 = interpol(profInt_A.y14, profInt_A.y4, wavelength)
    profInt_B.y15 = interpol(profInt_A.y15, profInt_A.y4, wavelength)
    profInt_B.y16 = interpol(profInt_A.y16, profInt_A.y4, wavelength)
    profInt_B.y17 = interpol(profInt_A.y17, profInt_A.y4, wavelength)
    profInt_B.y18 = interpol(profInt_A.y18, profInt_A.y4, wavelength)
    profInt_B.y19 = interpol(profInt_A.y19, profInt_A.y4, wavelength)
    profInt_B.y20 = interpol(profInt_A.y20, profInt_A.y4, wavelength)

    stmt17 = fjava_get_jdbc_statement(user='sorce', password='sorcedb', server='sorce-db', database="SORCE")
    conn17 = OBJ_NEW("oJava_DbExchange", stmt17)

    ; make sure we have an entry in the CalibrationMetadata table for this calibrationSetId
    q1="SELECT * from CalibrationMetadata where calibrationSetId="+strtrim(string(newCalId),2)
    if keyword_set(dbinsert) then begin
        date_str=string(systime(/jul),format='(C(CYI,"-",CMOI2.2,"-",CDI2.2," ",CHI2.2,":",CMI2.2,":",CSI2.2))')
        query_database,q1,data,nrows
        if nrows eq 0 then begin
            q1="INSERT INTO dbo.CalibrationMetadata(calibrationSetId, calibrationTableName, instrumentModeId, "
            q1+="version, releaseDate, effectiveDate, rationaleForNewVersion, responsibleIndividual, "
            q1+="dataDescription, sourceOfData, applicableReferences, additionalInformation) VALUES("
            q1=q1+strtrim(string(newCalId),2)+", 'SimProfileIntegralCal', "
            q1=q1+strtrim(string(instrumentModeId),2)+", "+strtrim(string(version),2)
            q1=q1+", '"+date_str+"', '2003-01-01 00:00:00.0'"
            q1=q1+", 'New version 19 by comparing data from mode 44 and mode 48 on day 453.67 to determine wavelength and y8'"
            q1=q1+", 'Stephane Beland'"
            q1=q1+", 'ProfileIntegral and wavelength as a function of ccd spot position'"
            q1=q1+", 'From mode 44 version 19 (calibrationSetId=1013)'"
            q1=q1+", 'none'"
            q1=q1+", 'none')"
            if keyword_set(verbose) then print,q1
            if keyword_set(dbinsert) then result=conn17->execute(q1)
        endif else begin
            ; update the releaseDate only
            q1="UPDATE CalibrationMetadata SET releaseDate='"+date_str+"'"
            q1=q1+" WHERE calibrationSetId="+strtrim(string(newCalId),2)
            if keyword_set(verbose) then print,q1
            if keyword_set(dbinsert) then result=conn17->execute(q1)
        endelse
    endif

    q1 = "INSERT INTO SimProfileIntegralCal(calibrationSetId, x, y0, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13, y14, y15, y16, y17, y18, y19, y20) VALUES("
    for i=0L,n_elements(profInt_B)-1 do begin
        str = strtrim(string(profInt_B[i].calibrationSetId,format='(I)'),2)
        str = str+', '+strtrim(string(profInt_B[i].x,format='(I)'),2)
        str = str+', '+strtrim(string(profInt_B[i].y0,format='(D0.2)'),2)
        str = str+', '+strtrim(string(profInt_B[i].y1,format='(D0.2)'),2)
        str = str+', '+strtrim(string(profInt_B[i].y2,format='(D0.2)'),2)
        str = str+', '+strtrim(string(profInt_B[i].y3),2)
        str = str+', '+strtrim(string(profInt_B[i].y4),2)
        str = str+', '+strtrim(string(profInt_B[i].y5),2)
        str = str+', '+strtrim(string(profInt_B[i].y6),2)
        str = str+', '+strtrim(string(profInt_B[i].y7),2)
        str = str+', '+strtrim(string(profInt_B[i].y8),2)
        str = str+', '+strtrim(string(profInt_B[i].y9),2)
        str = str+', '+strtrim(string(profInt_B[i].y10),2)
        str = str+', '+strtrim(string(profInt_B[i].y11),2)
        str = str+', '+strtrim(string(profInt_B[i].y12),2)
        str = str+', '+strtrim(string(profInt_B[i].y13),2)
        str = str+', '+strtrim(string(profInt_B[i].y14),2)
        str = str+', '+strtrim(string(profInt_B[i].y15),2)
        str = str+', '+strtrim(string(profInt_B[i].y16),2)
        str = str+', '+strtrim(string(profInt_B[i].y17),2)
        str = str+', '+strtrim(string(profInt_B[i].y18),2)
        str = str+', '+strtrim(string(profInt_B[i].y19),2)
        str = str+', '+strtrim(string(profInt_B[i].y20),2)
        q2 = q1+str+')'
        if keyword_set(verbose) then print,q2
        if keyword_set(dbinsert) then profInt_Bult=conn17->execute(q2)
    endfor

    return,profInt_b

end
