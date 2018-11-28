;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Converts SOLSTICE wavelengths to grating position.
;
; CALLING SEQUENCE:
;   wavelength = LAMBDA2GRT(wavelengths)
;
; INPUT PARAMETERS:
;   wavelengths -
;      wavelengths in nm
;
; OPTIONAL INPUT PARAMETERS:
;   None
;
; RETURNED PARAMETERS:
;   An array of the same size as wavelengths with the corresponding grating position.
;
; EXAMPLE:  
;
; COMMON BLOCKS:
;   None
;
; NOTES:
;
; REVISION HISTORY:
;   Revision: $Id: lambda2grt.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;------------------------------------------------------------------
;
function lambda2grt, wavelength, instrumentId=instrumentId

    ; TODO
    ; i'm guessing the constants are slightly different for each instrument mode

    offset = 239532.38d
    stepSize = 2.4237772022101214D-6
    d = 277.77777777777777D
    phiGInRads = 0.08503244115716374D
    cosphi = cos(phiGInRads / 2d)

    ang1 = asin(wavelength / 2d / d / cosphi)
    grtPos = offset - ang1 / stepSize

    return, grtPos

end

