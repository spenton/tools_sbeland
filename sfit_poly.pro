;+
; NAME:   SFIT_POLY
;
; PURPOSE: 
;    This routine uses the coefficients returned from SFIT routine and returns
;    an 2D array with the value evaluated at the specified X and Y locations.
;
; CALLING SEQUENCE:
;    result=sfit_poly(xarray, yarray, coefficients)
;
; INPUT PARAMETERS:
;    xarray - one dimensional array with the x locations
;    yarray - one dimensional array with the y locations
;    coefficients - coefficients returned by the routine SFIT which fits a 
;        polynomial surface to a set of X,Y and Z values
;
; OPTIONAL INPUT PARAMETERS:
;    none
;
; OPTIONAL INPUT KEYWORDS:
;    none 
;
; OUTPUT VALUES:
;    result - resturns n x m array where n and m are the number of elements in 
;        xarray and yarray respectively
;
; OUTPUT PARAMETERS:
;    none
;
; OPTIONAL OUTPUT PARAMETERS:
;    NONE
;
; EXAMPLE:  
;
; COMMON BLOCKS:
;    NONE
;
; NOTES:
;    The coefficients are expected to have the format when SFIT is called 
;    WITHOUT the MAX_DEGREE set.
;    For example, if SFIT is called with Degree=2 and MAX_DEGREE is not set, 
;    then the terms returned are [[k, y, y2], [x, xy, xy2], [x2, x2y, x2y2]]
;
; REVISION HISTORY:
;   Revision: $Id: sfit_poly.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;-
;*****************************************************************************
function sfit_poly, x, y, coeff
    
    ndim=size(coeff,/n_dim)
    nx=n_elements(x)
    ny=n_elements(y)
    if nx eq 0 or ny eq 0 then begin
        print,'Error: no input data provided'
        return,-1
    endif
    
    ; keep the same datatype as input
    coeff=double(coeff)
    xx=dblarr(nx,ny)
    yy=xx
    zz=xx
    for i=0,ny-1 do xx[*,i]=x
    for i=0,nx-1 do yy[i,*]=y

    sz=size(coeff,/dim)
    for i=0,sz[0]-1 do begin
        for j=0,sz[1]-1 do begin
            zz += coeff[j,i] * xx^i * yy^j
        endfor
    endfor

    return,zz
end
