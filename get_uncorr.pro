;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Extract either the time series or a specific scan from the IDL data
;   data structure in the save file containing the uncorrected irradiance
;   form version 20 of the SORCE SIM instrument.
;
; CALLING SEQUENCE:
;   values = get_uncorr(indata, timerange, wavelength=wavelength)
;
; INPUT PARAMETERS:
;   indata -
;      The input data structure (array of pointers, one for each scan)
;   timerange -
;      The timerange to extract data for (in Julian days).
;      If timerange contains a single value, the spectra closest in time
;      to the requested value will be returned (irradiance vs wavelength).
;      If the timeranege contains 2 values, it will be assumed the user is 
;      requesting a time series and the wavelength is required.
;   wavelength -
;      In the case where timerange has 2 values, the program will interpolate
;      the data for each spectra within those dates at the requested wavelength.
;
; OPTIONAL INPUT PARAMETERS:
;   LSQUADRATIC -
;       flag to the routine interpol that defines a least squares quadratic fit
;   QUADRATIC
;       flag to the routine interpol that defines a quadratic fit
;   SPLINE
;       flag to the routine interpol that defines a cubic spline  fit
;   CLOSEST
;       flag to skip the interpolation and return the irradiance from the 
;       closest wavelength relative to the requested wavelength
;
; RETURNED PARAMETERS:
;   A structure with either the single spectra or the time series.
;
; OPTIONAL OUTPUT PARAMETERS:
;      NONE
;
; EXAMPLE:  
;   Get the time series at 289.5nm from the UV data
;      IDL> res=get_uncorr(uva_uncorr_v20,sd2jd([0d,4200d]),wavelength=289.5,/spline)
;      IDL> help,res
;      ** Structure <9f144b8>, 3 tags, length=134864, data length=134864, refs=1:
;         TIMEJD          DOUBLE    Array[8429]
;         WAVELENGTH      FLOAT     Array[8429]
;         IRRADIANCE      FLOAT     Array[8429]
;
;   Get the spectra for mission day 453.67 (2004/04/21.67)
;      IDL> uv_spect =  get_uncorr(uva_uncorr_v20,sd2jd(453.67d))
;      IDL> help,uv_spect
;      ** Structure <9a635c8>, 3 tags, length=18224, data length=18224, refs=2:
;         TIMEJD          DOUBLE    Array[1139]
;         WAVELENGTH      FLOAT     Array[1139]
;         IRRADIANCE      FLOAT     Array[1139]
;      IDL> vis_spect = get_uncorr(visa_uncorr_v20,sd2jd(453.67d))
;      IDL> ir_spect =  get_uncorr(ira_uncorr_v20,sd2jd(453.67d))
;
; COMMON BLOCKS:
;      NONE
;
; REVISION HISTORY:
;   Revision: $Id: get_uncorr.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;-
;
FUNCTION GET_UNCORR, indata, timerange, wavelength=wavelength, help=help, $
    LSQUADRATIC=LSQUADRATIC, QUADRATIC=QUADRATIC, SPLINE=SPLINE, CLOSEST=CLOSEST

    if (n_params() lt 2) or keyword_set(help) then begin
        doc_library, 'GET_UNCORR'
        return, -1
    endif

    if size(indata,/tname) ne 'POINTER' then begin
        print,'Error: the input data is expected to be an array of pointers'
        return,-1
    endif

    ntimes=n_elements(timerange)
    if ntimes eq 0 then begin
        print,'Error: no timerange was provided'
        return,-1
    endif

    nscans = n_elements(indata)
    if ntimes eq 1 then begin
        ; return the data for a single scan closest in time to the requested timerange
        temp=dblarr(nscans)
        for i=0L,nscans-1 do temp[i]=(*indata[i]).timejd[0]
        value=min(abs(temp-timerange[0]), pos)
        outdata=(*indata[pos])
    endif else if ntimes ge 2 then begin
        ; return the data for the time series at the requested wavelength
        if n_elements(wavelength) eq 0 then begin
            print,'Error: a time range was specified but the requested wavelength is missing'
            return,-1
        endif
        wavelength=wavelength[0]
        waves=[]
        irradiance=[]
        timejd=[]
        for i=0L,nscans-1 do begin
            if (*indata[i]).timejd[0] lt timerange[0] or  (*indata[i]).timejd[0] gt timerange[1] then continue
            if keyword_set(CLOSEST) then begin
                mn=min(abs((*indata[i]).wavelength - wavelength),pos)
                waves = [waves,(*indata[i]).wavelength[pos]]
                irrd = (*indata[i]).irradiance[pos]
            endif else begin
                irrd=interpol((*indata[i]).irradiance, (*indata[i]).wavelength, wavelength, LSQUADRATIC=LSQUADRATIC, QUADRATIC=QUADRATIC, SPLINE=SPLINE)
                waves=[waves,wavelength]
            endelse
            irradiance=[irradiance,irrd]
            val=min(abs((*indata[i]).wavelength - wavelength),pos)
            timejd=[timejd,  (*indata[i]).timejd[pos]]
        endfor
        outdata={timejd:timejd, wavelength:waves, irradiance:irradiance}
    endif 


    return, outdata
END

