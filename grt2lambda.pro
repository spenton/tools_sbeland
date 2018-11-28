;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Converts SOLSTICE grating position to "initial" wavelengths.
;   These are approximate as the zero point of the grating position changes.
;
; CALLING SEQUENCE:
;   wavelength = GRT2LAMBDA(grtPos, prismTemp)
;
; INPUT PARAMETERS:
;   grtPos -
;      grating position 
;
; OPTIONAL INPUT PARAMETERS:
;   None
;
; RETURNED PARAMETERS:
;   An array of the same size as grtPos with the corresponding wavelength (in nm).
;
; EXAMPLE:  
;
; COMMON BLOCKS:
;   None
;
; NOTES:
;
; REVISION HISTORY:
;   Revision: $Id: grt2lambda.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;------------------------------------------------------------------
;
function grt2lambda, grtPos, instrumentId=instrumentId

    ; TODO
    ; i'm guessing the constants are slightly different for each instrument mode

    offset = 239532.38d
    stepSize = 2.4237772022101214D-6
    d = 277.77777777777777D
    phiGInRads = 0.08503244115716374D

    ang1 = (offset - double(grtPos)) * stepSize
    wavelength = 2d * d * sin(ang1) * cos(phiGInRads / 2d)

    return, wavelength

end

