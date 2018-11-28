;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Converts wavelengths to corresponding CCD subpixel positions for the specified
;   instrument mode.
;
; CALLING SEQUENCE:
;   ccdpos = Lambda2CCD(instrumentModeId, wavelength, prismTemp)
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
;   wavelength -
;      wavelength to convert (nm)
;   prismTemp - 
;      Prism temperature at the time of the measurement
;
; OPTIONAL INPUT PARAMETERS:
;   gamz - 
;      New value of gamma-Z to use (Radians).
;   pixsize - 
;      CCD pixel size to use (mm).
;   wedge -
;      New value for the prism wedge angle (Radians).
;
; RETURNED PARAMETERS:
;   An array of the same size as wavelength with the corresponding ccd position.
;
; EXAMPLE:  
;
; COMMON BLOCKS:
;      SIM_DEF
;
; NOTES:
;
; REVISION HISTORY:
;   Revision: $Id: lambda2ccd.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;------------------------------------------------------------------
;
function lambda2ccd, instrumentModeId, lambda, tfp, gamz=gamz, pixsize=pixsize, $
         focal=focal, fref=fref, yab=yab, cz=cz, wedge=wedge, tref=tref, tsis=tsis
    common SIM_DEF, sim_model

    ; check if the common block has been defined
    if n_elements(sim_model) eq 0 then get_sim_model

    pos=where(sim_model.imode eq instrumentModeId, count)
    if count eq 0 then return,-1 else pos=pos[0]

    if n_elements(gamz) eq 0 then gamz = sim_model.gamz[pos]
    if n_elements(pixsize) eq 0 then pixsize = sim_model.pix[pos]
    if n_elements(wedge) eq 0 then wedge = sim_model.thp[pos]
    if n_elements(focal) eq 0 then focal = sim_model.foc[pos]
    if n_elements(fref) eq 0 then fref = sim_model.fref[pos]
    if n_elements(yab) eq 0 then yab = sim_model.yab[pos]
    if n_elements(cz) eq 0 then cz = sim_model.cz[pos]
    if n_elements(tref) eq 0 then tref = sim_model.TREF_SILICA

    THP2 = 2d*wedge     ; Twice prism angle
    STHP2 = SIN(THP2)   ; sin(twice prism angle)
    CTHP2 = COS(THP2)   ; twice cos(twice prism angle)
    PHI = ATAN(yab / focal)  ; DEVIATION ANGLE PHI AT SLIT JSLIT
    SP = SIN(PHI)
    CP = COS(PHI)

    ; get the corresponding index of refraction (convert from nm to microns)
    LAM = LAMBDA / 1d3
    LSQ = LAM * LAM
    LAM1 = 1d/LAM
    if keyword_set(TSIS) then begin
        ; from Dave Harber's code as of 08/08/2016
        p = [-3.6323d, 9.474d, 133.788d, 0.19303d, 2861.201d]
        tref = 20d
        
        ; from Dave Harber's 8/17/2017 Fitting
        KS1 = 0.64693485d
        KS2 = 0.45821956d
        KS3 = 0.92665664d
        LS1 = 0.0042958495d
        LS2 = 0.013094520d
        LS3 = 101.53727d
        Y = KS1*LSQ/(LSQ-LS1) + $
            KS2*LSQ/(LSQ-LS2) + $
            KS3*LSQ/(LSQ-LS3) + 1D

        dn_dt = p[0] + p[1]*(LSQ*1d6)/abs((LSQ*1d6) - p[2]^2d) + p[3]*(LSQ*1d6)/abs((LSQ*1d6) - p[4]^2d)
        N = SQRT(Y) * (1d + 1d-6*dn_dt*(TFP-TREF))
    endif else begin
        Y = sim_model.KS1*LSQ/(LSQ-sim_model.LS1) + $
            sim_model.KS2*LSQ/(LSQ-sim_model.LS2) + $
            sim_model.KS3*LSQ/(LSQ-sim_model.LS3) + 1D
        ; apply the temperature correction
        DELTA_N = POLY(LAM1, sim_model.c_tcn) * (TFP - TREF) * 1D-6
        N = SQRT(Y) + DELTA_N
    endelse

    X = 1D/N^2  ; 1/N^2
    SGSQ = CP*STHP2^2 + X*CTHP2*SP^2 + SP*STHP2*SQRT( (CTHP2+X*CP)^2 - (1D -X)^2 )
    SGSQ = SGSQ / (2D*X*(CTHP2+CP))
    NEW_GAMMA = ASIN(SQRT(SGSQ))  ; ARCSIN OF THE SQUARE ROOT OF SIN_GAMMA^2
    YY = fref * TAN(2D * (NEW_GAMMA-GAMZ) ) ; CCD SPOT POSITION
    CCD = cz + YY/pixsize

    RETURN,CCD
END
