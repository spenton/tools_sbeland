  function irradcnv, UnitsID, Lambda
;
; Return the irradiance conversion factor, to convert _from_ watt/m^3 _to_
; the specified irradiance unit.  Wavelength(s) (Lambda, nm) are required
; for conversions to photon fluxes, but this argument is ignored in all
; other cases.
;
; B. Knapp, 1998-08-18
;
  UnitsValid = ['W/M^3', 'STANDARD', 'W/CM^3', 'MW/M^2/NM',$
     'ERGS/S/CM^2/A', 'PHOTONS/S/CM^2/NM', 'PHOTONS/S/CM^2/A' ]
  UnitsConv  = [ 1.0,     1.0,        1.E-06,   1.E-06,    $
      1.E-07,          503.402,             50.3402           ]
;
; Show usage?
  if n_params() eq 0 then begin
     print,' '
     print,' IrradCnv returns the conversion factor _from_ Watt/m^3'
     print,' _to_ the specified irradiance units.  Wavelength(s) in'
     print,' nm are required to convert to Photon fluxes.'
     print,' '
     print,'     c = IrradCnv( UnitsID [, Lambda] )'
     return,' '
  endif
;
; Check input
  UnitsID_loc = strupcase( strtrim( UnitsID, 2 ) )
  pick = (where( UnitsValid eq UnitsID_loc ))[0]
  if pick lt 0 then begin
     message, 'Invalid irradiance unit ID', /info
     return, 1.
  endif
;
; Return the appropriate conversion factor
  if strpos( UnitsID_loc, 'PHOTON' ) ge 0 then begin
;
;    Lambda must be defined
     if n_elements( Lambda ) eq 0 then begin
        message,'Lambda not defined', /info
        return, 1.
     endif
     return, UnitsConv[pick]*Lambda
  endif else begin
     return, UnitsConv[pick]
  endelse
;
  end
