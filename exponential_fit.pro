function my_expo, x, a
; EXPONENTIAL.PRO
; y=a0*exp(a1*x)+a2 function in format for LMFIT
  y = a[0]*exp(a[1]*x)+a[2]
  return, [[y], [exp(a[1]*x)], [a[0]*x*exp(a[1]*x)], replicate(1.0,n_elements(x))]
end

function my_expo_oneminus, x, a
; EXPONENTIAL.PRO
; y=1 - (a0*exp(a1*x)+a2) function in format for LMFIT
  y = 1d - (a[0]*exp(a[1]*x)+a[2])
  return, [[y], [-exp(a[1]*x)], [-a[0]*x*exp(a[1]*x)], replicate(-1.0,n_elements(x))]
end

function exponential_fit, x, y, guess=guess, oneminus=oneminus, $
         sigma=sigma, errors=errors, yfit=yfit, fita=fita
;+
; NAME:
;   EXPONENTIAL_FIT
; PURPOSE:
;   Fit an exponential to the data of the form y=a*exp(b*x).  It's
;   just a wrapper for LMFIT, as always.
;
; CALLING SEQUENCE:
;   fit=exponential_fit(x, y[, errors=errors, guess=guess, sigma=sigma])
;
; INPUTS:
;   X, Y -- X and Y values to be fit.
;
; KEYWORD PARAMETERS:
;   GUESS -- vector containing the guess values [a,b] in y=a*exp(b*x)
;   SIGMA -- Named variable to contain the errors in the fit parameters.
;   ERRORS -- error in the y_values
; 
; REQUIRES:
;  EXPONENTIAL.pro -- See FITFCNS
;
; OUTPUTS:
;   COEFFS -- Coefficients of the powerlaw fit in vector form [a,b];
;
; MODIFICATION HISTORY:
;       Documented.
;       Wed Nov 21 11:49:45 2001, Erik Rosolowsky <eros@cosmic>
;-



  if not keyword_set(guess) then begin
    goodind = where(y gt 0)
    ;coeffs = linear_fit(x[goodind], alog(y[goodind]))
    coeffs = robust_poly_fit(x[goodind], alog(y[goodind]),1,/double)
    guess = [exp(coeffs[0]), coeffs[1], 1.0d-5]
  endif
  a = guess
  if keyword_set(oneminus) then fname="my_expo_oneminus" else fname="my_expo"
  if n_elements(fita) eq 0 then fita=intarr(n_elements(a))+1

  for i = 0, 3 do begin
    fit = lmfit(x, y, a, measure_errors=errors, /double, sigma=sigma, $
                function_name=fname, conv=conv, fita=fita)
    if conv eq 1 then break
  endfor
  res=my_expo(x,a)
  yfit=res[0:n_elements(x)-1]

  return, a
end
