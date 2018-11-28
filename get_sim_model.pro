;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Generate and populate the COMMON block containing the instrument model parameters.
;
; CALLING SEQUENCE:
;   get_sim_model
;
; INPUT PARAMETERS:
;   NONE
;
; OPTIONAL INPUT PARAMETERS:
;   NONE
;
; RETURNED PARAMETERS:
;   NONE
;
; OPTIONAL OUTPUT PARAMETERS:
;   NONE
;
; EXAMPLE:  
;
; COMMON BLOCKS:
;   SIM_DEF
;
;
; REVISION HISTORY:
;   Revision: $Id: get_sim_model.pro,v 1.1 2018/10/01 16:57:13 spenton Exp $
;------------------------------------------------------------------
;
pro get_sim_model
    common SIM_DEF, sim_model

    ; change in index of refraction with temperature modeled by a 5th order polynomial
    c_tcn = [10.91141546d, 4.130107857d, -7.400374025d, 4.110286051d, -0.937939315d, 0.078332038d]
    
    ; these new coefficients were obtained from TSIS
    ;c_tcn = [14.672115325927734d, -37.142608642578125d,  95.550659179687500d, -132.09014892578125d, 103.68243408203125d, -43.525695800781250d, 7.6180877685546875D]

    TREF_SILICA = 20.0D    ; ZEMAX AND MALITSON SELLMEIER TEMPERATURE
    
    ; SELLMEIER COEFS
    KS1 = 0.6961663D
    KS2 = 0.4079426D  
    KS3 = 0.8974794D
    LS1 = 0.004679148D
    LS2 = 0.01351206D
    LS3 = 97.934D

    iMode = [31,32,41,43,44,45,47,48]

    thp = dblarr(n_elements(iMode))
    pix = thp
    cz = thp
    gamz = thp
    fref = thp
    foc = thp
    yab = thp

    ; 31: ESR A
        ;thp[0] = 34.516068d * !DPI / 180d   ; prism wedge angle in radians 0: sima, 1:simb
        thp[0] = 34.497529d * !DPI / 180d   ; updated values from Juan Fontenla May 2006
        ;pix[0] = -0.0012886090d            ; pixel size in mm
        ;pix[0] = -0.001288718933d            ; updated to match Kurutz on day 453
        pix[0] = -0.0012886101216771186d      ; updated to match Kurutz on day 453 (SBeland 2014-05-28)
        cz[0] =  45384d                     ; pixel position when beam is folded back on entrance slit
        ;gamz[0]= 57.040113d * !DPI / 180d  ; GammaZ in radians
        ;gamz[0]= 57.03981216240179d * !DPI / 180d  ; updated to match Kurutz on day 453
        gamz[0]= 57.040119483936316d * !DPI / 180d  ; updated to match Kurutz on day 453 (SBeland 2014-05-28)
        fref[0] = 391.15d                   ; Reference mirror to CCD
        foc[0] =   403.15d                   ; Entrance slit to prism rotation axis (mm)
        yab[0]= 34.976d
    ;32: ESR B
        ;thp[1] = 34.53d * !DPI / 180d   ; prism wedge angle in radians 0: sima, 1:simb
        thp[1] = 34.5078d * !DPI / 180d   ; updated values from Juan Fontenla May 2006
        ;pix[1] = 0.0013d            ; pixel size in mm
        ;pix[1] = 0.001295295d            ; updated values from Juan Fontenla May 2006
        pix[1] = 0.001296133116004082d            ; updated to match Kurutz on day 453
        cz[1] =  14615d                     ; pixel position when beam is folded back on entrance slit
        ;gamz[1]= 56.84556d * !DPI / 180d  ; GammaZ in radians
        ;gamz[1]= 56.843113d * !DPI / 180d  ; updated values from Juan Fontenla May 2006
        gamz[1]= 56.84169886261460d * !DPI / 180d  ; updated to match Kurutz on day 453
        ;fref[1] = 391.7d                   ; Reference mirror to CCD
        fref[1] = 390.08d                   ; updated values from Juan Fontenla May 2006
        ;foc[1] =   408.43d                   ; Entrance slit to prism rotation axis (mm)
        foc[1] =   402.098d                   ; updated values from Juan Fontenla May 2006
        yab[1] = 34.959d
    ;41: VIS A
        ;thp[2] = 34.516068d * !DPI / 180d   ; prism wedge angle in radians 0: sima, 1:simb
        thp[2] = 34.490364d * !DPI / 180d   ; updated values from Juan Fontenla May 2006
        ;pix[2] = -0.0012886090d            ; pixel size in mm
        pix[2] = -0.001288616302358402d            ; pixel size in mm aligned with Kurucz on day 453.67
        cz[2] =  45384d                     ; pixel position when beam is folded back on entrance slit
        ;gamz[2]= 57.040113d * !DPI / 180d  ; GammaZ in radians
        gamz[2]= 57.04009570428246d * !DPI / 180d  ; GammaZ in radians aligned with Kurucz on day 453.67
        fref[2] = 391.15d                   ; Reference mirror to CCD
        foc[2] =   403.15d                   ; Entrance slit to prism rotation axis (mm)
        ;yab[2]= 54.908d
        yab[2]= 49.924d                       ; updated values from Juan Fontenla May 2006
    ;43: UV A
        thp[3] = 34.516068d * !DPI / 180d   ; prism wedge angle in radians 0: sima, 1:simb
        ;pix[3] = -0.001288609d            ; pixel size in mm
        pix[3] = -0.0012901790d            ; updated to match Kurutz on day 453.67
        cz[3] =  45384d                     ; pixel position when beam is folded back on entrance slit
        ;gamz[3]= 57.040113d * !DPI / 180d  ; GammaZ in radians
        gamz[3]= 57.03773690d * !DPI / 180d  ; updated to match Kurutz on day 453.67
        fref[3] = 391.15d                   ; Reference mirror to CCD
        foc[3] =   403.15d                   ; Entrance slit to prism rotation axis (mm)
        yab[3]= -10.106d
    ;44: IR A
        ;thp[4] = 34.516068d * !DPI / 180d   ; prism wedge angle in radians 0: sima, 1:simb
        thp[4] = 34.48554d * !DPI / 180d   ; updated values from Juan Fontenla May 2006
        ;pix[4] = -0.0012886090d            ; pixel size in mm
        pix[4] = -0.001288621730d            ; updated to match Kurutz on day 453.67
        ;pix[4] = -0.0012882619643413598d     ; updated to match Kurutz with TSIS Y18
        cz[4] =  45384d                     ; pixel position when beam is folded back on entrance slit
        ;gamz[4]= 57.040113d * !DPI / 180d  ; GammaZ in radians
        gamz[4]= 57.040087017140d * !DPI / 180d  ; updated to match Kurutz on day 453.67
        ;gamz[4]= 57.041129785346733d * !DPI / 180d  ; updated to match Kurutz with TSIS Y18
        fref[4] = 391.15d                   ; Reference mirror to CCD
        foc[4] =   403.15d                   ; Entrance slit to prism rotation axis (mm)
        yab[4]= 59.988d
    ;45: VIS B
        ;thp[5] = 34.53d * !DPI / 180d   ; prism wedge angle in radians 0: sima, 1:simb
        thp[5] = 34.504933d * !DPI / 180d   ; updated values from Juan Fontenla May 2006
        ;pix[5] = 0.0013d            ; pixel size in mm
        ;pix[5] = 0.001295295d            ; updated values from Juan Fontenla May 2006
        ;pix[5] = 0.001297917878d            ; updated to match Kurutz on day 453.67
        ;pix[5] = 0.0012970639604886053d      ; updated to match SIMA on day 453.67
        pix[5] = 0.0012979179116849309d      ; updated to match SIMA on day 453.67  V23 SBeland 6/3/2016
        cz[5] =  14615d                     ; pixel position when beam is folded back on entrance slit
        ;gamz[5]= 56.84556d * !DPI / 180d  ; GammaZ in radians
        ;gamz[5]= 56.843113d * !DPI / 180d  ; updated values from Juan Fontenla May 2006
        ;gamz[5]= 56.837154719256d * !DPI / 180d  ; updated to match Kurutz on day 453.67
        ;gamz[5]= 56.839101355178485d * !DPI / 180d  ; updated to match SIMA on day 453.67
        gamz[5]= 56.837154658725716d * !DPI / 180d  ; updated to match SIMA on day 453.67 V23 SBeland 6/3/2016
        ;fref[5] = 391.7d                   ; Reference mirror to CCD
        fref[5] = 390.08d                   ; updated values from Juan Fontenla May 2006
        ;foc[5] =   408.43d                   ; Entrance slit to prism rotation axis (mm)
        foc[5] =   402.098d                   ; updated values from Juan Fontenla May 2006
        yab[5] = 50.084d
    ;47: UV B
        ;thp[6] = 34.53d * !DPI / 180d   ; prism wedge angle in radians 0: sima, 1:simb
        thp[6] = 34.519877d * !DPI / 180d   ; updated values from Juan Fontenla May 2006
        ;pix[6] = 0.0013d            ; pixel size in mm
        ;pix[6] = 0.001295295d            ; updated values from Juan Fontenla May 2006
        ;pix[6] = 0.0012990106d            ; updated values from align_spectra2 Nov. 2013
        ;pix[6] = 0.0012984572225789473d    ; updated values from align_spectra2 Nov. 2015
        pix[6] = 0.0012988281844549651d    ; updated values from align_spectra2 to SimA June 2016
        cz[6] =  14615d                     ; pixel position when beam is folded back on entrance slit
        ;gamz[6]= 56.84556d * !DPI / 180d  ; GammaZ in radians
        ;gamz[6]= 56.843113d * !DPI / 180d  ; updated values from Juan Fontenla May 2006
        ;gamz[6]= 56.840331d * !DPI / 180d  ; updated values from align_spectra2 Nov. 2013 aligned with mode 43 on day 453.67
        ;gamz[6]= 56.840292432565164d * !DPI / 180d  ; updated values from align_spectra2 Nov. 2015 
        gamz[6]= 56.840230377229680d * !DPI / 180d  ; updated values from align_spectra2 to SimA June 2016
        ;fref[6] = 391.7d                   ; Reference mirror to CCD
        fref[6] = 390.08d                   ; updated values from Juan Fontenla May 2006
        ;foc[6] =   408.43d                   ; Entrance slit to prism rotation axis (mm)
        foc[6] =   402.098d                   ; updated values from Juan Fontenla May 2006
        yab[6] = -10.040d
    ;48: IR B
        ;thp[7] = 34.53d * !DPI / 180d   ; prism wedge angle in radians 0: sima, 1:simb
        thp[7] = 34.505201d * !DPI / 180d   ; updated values from Juan Fontenla May 2006
        ;pix[7] = 0.0013d            ; pixel size in mm
        ;pix[7] = 0.001295295d            ; updated values from Juan Fontenla May 2006
        ;pix[7] = 0.001295331966d            ; updated to match Kurutz on day 453.67
        pix[7] = 0.0012907336566113222d            ; updated to match Kurutz on day 453.67
        cz[7] =  14615d                     ; pixel position when beam is folded back on entrance slit
        ;gamz[7]= 56.84556d * !DPI / 180d  ; GammaZ in radians
        ;gamz[7]= 56.843113d * !DPI / 180d  ; updated values from Juan Fontenla May 2006
        ;gamz[7]= 56.843030393245d * !DPI / 180d  ; updated to match Kurutz on day 453.67
        gamz[7]= 56.863954671044318d * !DPI / 180d  ; updated to match Kurutz on day 453.67
        ;fref[7] = 391.7d                   ; Reference mirror to CCD
        fref[7] = 390.08d                   ; updated values from Juan Fontenla May 2006
        ;foc[7] =   408.43d                   ; Entrance slit to prism rotation axis (mm)
        foc[7] =   402.098d                   ; updated values from Juan Fontenla May 2006
        yab[7] = 60.079d

    sim_model = {KS1:KS1, KS2:KS2, KS3:KS3, LS1:LS1, LS2:LS2, LS3:LS3, C_TCN:C_TCN, TREF_SILICA:TREF_SILICA, $
         iMode:iMode, thp:thp, pix:pix, cz:cz, gamz:gamz, fref:fref, foc:foc, yab:yab}


end

