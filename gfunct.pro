;+
;
; From the IDL online help example.
;
; Fit a function of the form F(x) = a * exp(b*x) + c to sample pairs contained in arrays X and Y. 
; The partial derivatives are easily computed symbolically:
;
; df/da = EXP(b*x)
; df/db = a * x * EXP(b*x)
; df/dc = 1.0
;-
PRO gfunct, X, A, F, pder
    bx = EXP(A[1] * X)
    F = A[0] * bx + A[2]
 
    ;If the procedure is called with four parameters, calculate the
    ;partial derivatives.
    IF N_PARAMS() GE 4 THEN $
        pder = [[bx], [A[0] * X * bx], [replicate(1.0, N_ELEMENTS(X))]]
END

