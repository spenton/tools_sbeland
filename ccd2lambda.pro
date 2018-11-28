;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Converts CCD subpixel positions to wavelengths for the specified
;   instrument mode.
;
; CALLING SEQUENCE:
;   wavelength = CCD2Lambda(instrumentModeId, ccdPos, prismTemp)
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
;   ccdPos -
;      subpixel position 
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
;   tsis -
;      Uses the Sellmeier coefficients and prism temperature correction
;      obtained from TSIS ground calibration (from Dave Harber)
;
; RETURNED PARAMETERS:
;   An array of the same size as ccdPos with the corresponding wavelength (in nm).
;
; EXAMPLE:  
;
; COMMON BLOCKS:
;      SIM_DEF
;
; NOTES:
;
; REVISION HISTORY:
;   Revision: $Id: ccd2lambda.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;------------------------------------------------------------------
;
function ccd2lambda, instrumentModeId, ccd, tfp, gamz=gamz, pixsize=pixsize, $
             wedge=wedge, focal=focal, cz=cz, fref=fref, yab=yab, tref=tref, tsis=tsis
    common SIM_DEF, sim_model

    ; check if the common block has been defined
    if n_elements(sim_model) eq 0 then get_sim_model

    pos=where(sim_model.imode eq instrumentModeId, count)
    if count eq 0 then return,-1 else pos=pos[0]

    if n_elements(wedge) eq 0 then wedge = sim_model.thp[pos]
    if n_elements(pixsize) eq 0 then pixsize = sim_model.pix[pos]
    if n_elements(cz) eq 0 then cz = sim_model.CZ[pos]
    if n_elements(gamz) eq 0 then gamz = sim_model.gamz[pos]
    if n_elements(fref) eq 0 then fref = sim_model.fref[pos]
    if n_elements(focal) eq 0 then focal = sim_model.foc[pos]
    if n_elements(yab) eq 0 then yab = sim_model.yab[pos]
    if n_elements(tref) eq 0 then tref = sim_model.TREF_SILICA

    THP2 = 2d*wedge  ; Twice prism angle
    STHP2 = SIN(THP2)      ; sin(twice prism angle)
    CTHP22 = 2d*COS(THP2)   ; twice cos(twice prism angle)

    YY = pixsize*(CCD-CZ)    ;
    NEW_GAMMA = gamz + 0.5D * ATAN(YY/FREF )  ; PRISM INPUT ANGLE GAMMA

    PHI = ATAN(YAB / focal )  ; DEVIATION ANGLE PHI AT Y
    SGAM = SIN(NEW_GAMMA)
    SGPHI = SIN(NEW_GAMMA - PHI)
    INDEX_REFRAC = SQRT(SGAM^2 + CTHP22 * SGAM * SGPHI + SGPHI^2) / STHP2

    LSQ = 1D + INDEX_REFRAC*0D ; INITIAL WILD GUESS FOR LAMBDA SQUARED
    NSQ_IN = INDEX_REFRAC*INDEX_REFRAC

    ; ITERATE NSQ VS LAMSQ
    TOL = 1D-10
    F = 1D12
    NITER = 0  ; INIT ITERATION COUNTER
    ITERMAX = 30

    if keyword_set(TSIS) then begin
        ; define the TSIS dn_dt coefficients
        ;p = [-3.6323d, 11.706d, 133.788d, 0.25868d, 2861.201d]
        p = [-3.6323d, 9.474d, 133.788d, 0.19303d, 2861.201d]
        tref = 20d
        
        ; from Dave Harber's 8/17/2017 Fitting
        KS1 = 0.64693485d
        KS2 = 0.45821956d
        KS3 = 0.92665664d
        LS1 = 0.0042958495d
        LS2 = 0.013094520d
        LS3 = 101.53727d
    endif else begin
        C_TCN=sim_model.C_TCN
        KS1=sim_model.KS1
        KS2=sim_model.KS2
        KS3=sim_model.KS3
        LS1=sim_model.LS1
        LS2=sim_model.LS2
        LS3=sim_model.LS3
    endelse


    WHILE MAX(ABS(F)) GT TOL DO BEGIN
        NITER = NITER + 1 ; ITERATION COUNTER
        ; INDEX FUNCTION
        NSQ = KS1*LSQ/(LSQ-LS1) + $
              KS2*LSQ/(LSQ-LS2) + $
              KS3*LSQ/(LSQ-LS3) +1.D0 
        F = NSQ-NSQ_IN
        DF_DLSQ = KS1*LS1/(LSQ-LS1)^2 + $
                  KS2*LS2/(LSQ-LS2)^2 + $
                  KS3*LS3/(LSQ-LS3)^2 ; ACTUALLY NEGATIVE OF THIS
        LSQ = (LSQ + F/DF_DLSQ) >0.02d ; CLIP UNDERSHOOT (USED TO HAVE 0.04 HERE, CAUSED TROUBLE AT LAM = 0.2)

        ; HAVE TO ITERATE THE TEMPERATURE CORRECTION ALSO, AS THE WAVELENGTH BECOMES TRUE
        ; SKIP THE TEMPERATURE CORRECTION FOR NOW
        IF NITER GE 2 THEN BEGIN
            if NOT keyword_set(TSIS) then begin
                X = 1D/SQRT(LSQ)   ; POLYNOMINAL X-AXIS.  See FIT_SILICA_TEMP_COEFS2.PRO
                DN = POLY(X, C_TCN) * (TFP-TREF) * 1D-6 
                NSQ_IN = (INDEX_REFRAC-DN)^2D

            endif else begin
                ; use TSIS dn_dt for testing (the coefficients are for lambda in nm instead of microns)
                ; dn_dt is in ppm
                dn_dt = p[0] + p[1]*(LSQ*1d6)/abs((LSQ*1d6) - p[2]^2d) + p[3]*(LSQ*1d6)/abs((LSQ*1d6) - p[4]^2d)
                NSQ_IN = (INDEX_REFRAC*(1d - 1d-6*dn_dt*(TFP-TREF))) ^ 2d
            endelse
        ENDIF

        IF NITER GE ITERMAX THEN BREAK
    ENDWHILE

    ; return the wavelengths in nm
    RETURN, SQRT(LSQ) *1d3

end

