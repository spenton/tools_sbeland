

FUNCTION REF_INDEX,LAM,TFP

; usage:    N = INDEX(LAMBDA,TFP)  ; LAMBDA IN MICRONS

; EVALUATE the Sellmeier1 equation for fused silica

;  INCLUDES TEMPERATURE EFFECT

; CALLING PROG IS RESPONSIBLE TO FIND PROPER TEMPERATURE.



    common SIM_DEF, sim_model



    ; check if the common block has been defined

    if n_elements(sim_model) eq 0 then get_sim_model



    LSQ = LAM*LAM  ; WAVELENGTHS IN MICRONS

    Y = sim_model.KS1*LSQ/(LSQ-sim_model.LS1)+sim_model.KS2*LSQ/(LSQ-sim_model.LS2)+sim_model.KS3*LSQ/(LSQ-sim_model.LS3)  ; Y = N^2-1

    X = 1.D0/LAM   ; POLYNOMIAL X-AXIS.  See FIT_SILICA_TEMP_COEFS2.PRO

    DELTA_N = 1.D-6*(sim_model.C_TCN[0]+X*(sim_model.C_TCN[1]+X*(sim_model.C_TCN[2]+X*(sim_model.C_TCN[3]+X*(  $

    sim_model.C_TCN[4]+X*sim_model.C_TCN[5] )))))*(TFP-sim_model.TREF_SILICA)

    N = SQRT(Y+1.D0)+ DELTA_N

    RETURN,N   ; INDEX

END



