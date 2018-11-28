;+
; Author: Stephane Beland
;
; PURPOSE: 
;
; CALLING SEQUENCE:
;
; INPUT PARAMETERS:
;   startTime -
;      The lower time range for which data will be returned.
;      May be specified as julian days, gps microseconds, or
;      SORCE mission day numbers.
;
;   stopTime -
;      The upper time range for which data will be returned
;      May be specified as julian days, gps microseconds, or
;      SORCE mission day numbers.
;
; OPTIONAL INPUT PARAMETERS:
;   version -
;      Requests a specific version of the SimProfileIntegralCal, otherwise
;      pick the latest version valid for the specified time range.
;   closest -
;      Skips the interpolation to get the irradiance at requested wavelength
;      but simply uses the value from the closest wavelenegth
;
; RETURNED PARAMETERS:
;
; OPTIONAL OUTPUT PARAMETERS:
;   coeffs -
;      The coefficients of the best polynomial fit.
;   status -
;      Returned STATUS of the mpfitfun routine.
;
; EXAMPLE:  
;
; COMMON BLOCKS:
;      NONE
;
; NOTES:

;
; REVISION HISTORY:
;   Revision: $Id: tnmin_esrdeg.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;------------------------------------------------------------------
;
;  Function to calculate the difference between Calibrated irradiance for SimA and SimB
;  using the same degradation model
function mytnmin_func, x, df, solexp31=solexp31, ird31=ird31, solexp32=solexp32, ird32=ird32, $
    new_ird31=new_ird31, new_ird32=new_ird32
    ; measured_irrad = SolarIrrad * (1 + exp(-abs(x[0]*solexp + x[1])) - exp(-abs(x[1])) )
    ; SolarIrrad = measured_irrad / (1 + exp(-abs(x[0]*solexp + x[1])) - exp(-abs(x[1])) )
    new_ird31 = ird31 / (1d + exp(-abs(x[0]*solexp31 + x[1])) - exp(-abs(x[1])))
    new_ird32 = ird32 / (1d + exp(-abs(x[0]*solexp32 + x[1])) - exp(-abs(x[1])))
    stats = moment(new_ird32 - new_ird31,mdev=out_value)
    ; evaluate the derivative 
    dfd0  = (exp(-abs(x[0]*solexp31 + x[1])) * solexp31) / ird31
    dfd0 -= (exp(-abs(x[0]*solexp32 + x[1])) * solexp32) / ird32
    dfd1  = (exp(-abs(x[1])) - exp(-abs(x[0]*solexp32 + x[1]))) / ird32
    dfd1 -= (exp(-abs(x[1])) - exp(-abs(x[0]*solexp31 + x[1]))) / ird31
    df=moment(df0, mdev=out0)
    df=moment(df1, mdev=out1)
    df = [out0,out1]
    print,x[0],x[1],out_value
    return, out_value
end

;*******************************************************************
;  Function to calculate the difference between Calibrated irradiance for SimA and SimB
;  using the same degradation model
function mytnmin_func2, x, df, solexp31=solexp31, ird31=ird31, solexp32=solexp32, ird32=ird32, $
    new_ird31=new_ird31, new_ird32=new_ird32
    ; measured_irrad = SolarIrrad * exp(-abs(x[0]*solexp + x[1]))
    ; SolarIrrad = measured_irrad / exp(-abs(x[0]*solexp + x[1]))
    ; SolarIrrad = measured_irrad * exp(abs(x[0]*solexp + x[1]))
    new_ird31 = ird31 * exp(abs(x[0]*solexp31 + x[1]))
    new_ird32 = ird32 * exp(abs(x[0]*solexp32 + x[1]))
    stats = moment(new_ird32 - new_ird31,mdev=out_value)
    out_value=abs(stats[0])
    ; evaluate the derivative 
    ;temp31 = ird31 * exp(abs(x[0]*solexp31 + x[1])) * (x[0]*solexp31 + x[1]) / abs(x[0]*solexp31 + x[1])
    ;temp32 = ird32 * exp(abs(x[0]*solexp32 + x[1])) * (x[0]*solexp32 + x[1]) / abs(x[0]*solexp32 + x[1])
    ;dfd0 = temp31 * solexp31 - temp32 * solexp32
    ;dfd1 = temp31 - temp32
    ;df=moment(df0, mdev=out0)
    ;df=moment(df1, mdev=out1)
    ;df = [out0,out1]
    print,x[0],x[1],out_value
    return, out_value
end

;*******************************************************************
;  Function to calculate the difference between Calibrated irradiance for SimA and SimB
;  using the same degradation model
function mytnmin_func3, x, new_ird31=new_ird31, new_ird32=new_ird32
    common toto, solexp31, ird31, solexp32, ird32
    ; measured_irrad = SolarIrrad * exp(-abs(x[0]*solexp + x[1]))
    ; SolarIrrad = measured_irrad / exp(-abs(x[0]*solexp + x[1]))
    ; SolarIrrad = measured_irrad * exp(abs(x[0]*solexp + x[1]))
    new_ird31 = ird31 * exp(abs(x[0]*solexp31 + x[1]))
    new_ird32 = ird32 * exp(abs(x[0]*solexp32 + x[1]))
    stats = moment(new_ird32 - new_ird31,sdev=out_value)
    return, out_value
end

function mytnmin_dfunc3, x
    common toto, solexp31, ird31, solexp32, ird32
    ; evaluate the derivative 
    temp31 = ird31 * exp(abs(x[0]*solexp31 + x[1])) * (x[0]*solexp31 + x[1]) / abs(x[0]*solexp31 + x[1])
    temp32 = ird32 * exp(abs(x[0]*solexp32 + x[1])) * (x[0]*solexp32 + x[1]) / abs(x[0]*solexp32 + x[1])
    dfd0 = moment(temp31 * solexp31 - temp32 * solexp32, sdev=out0)
    dfd1 = moment(temp31 - temp32, sdev=out1)
    return, [out0,out1]
end

;*******************************************************************
;  Function to calculate the difference between Calibrated irradiance for SimA and SimB
;  using the same degradation model
function mytnmin_func4, x, df, solexp31=solexp31, ird31=ird31, solexp32=solexp32, ird32=ird32, $
    new_ird31=new_ird31, new_ird32=new_ird32
    ; measured_irrad = SolarIrrad * (1 - (x[1]*solexp + x[0]))
    ; SolarIrrad = measured_irrad / (1 - (x[1]*solexp + x[0])) 
    new_ird31 = ird31 / (1d - (x[0] + x[1]*solexp31))
    new_ird32 = ird32 / (1d - (x[0] + x[1]*solexp32))
    stats = moment(new_ird32 - new_ird31,sdev=out_value)
    out_value=abs(stats[0])
    ; evaluate the derivative 
    df=[1.0,1.0]
    print,x[0],x[1],out_value
    return, out_value
end

;*******************************************************************
;  Function to calculate the difference between Calibrated irradiance for SimA and SimB
;  using the same degradation model
function mytnmin_func5, x, df, solexp31=solexp31, ird31=ird31, solexp32=solexp32, ird32=ird32, $
    new_ird31=new_ird31, new_ird32=new_ird32
    ; measured_irrad = SolarIrrad * (1 - (x[0]*solexp + x[1]))
    ; SolarIrrad = measured_irrad / (1 - (x[0]*solexp + x[1])) 
    new_ird31 = (x[0] + x[1]*solexp31) * exp(-abs(x[2]*solexp31)) - x[1]*solexp31
    new_ird32 = (x[0] + x[1]*solexp32) * exp(-abs(x[2]*solexp32)) - x[1]*solexp32
    stats = moment(new_ird32 - new_ird31,sdev=out_value)
    out_value=abs(stats[0])
    df=[1.0,1.0,1.0]
    print,x[0],x[1],x[2],out_value
    return, out_value
end

;*******************************************************************
function mytnmin_func6, x, df, solexp31=solexp31, ird31=ird31, solexp32=solexp32, ird32=ird32, $
    new_ird31=new_ird31, new_ird32=new_ird32
    ; measured_irrad = SolarIrrad * (1 - x[0]) * exp(-x[1]*solexp)
    ; SolarIrrad = measured_irrad / (1 - x[0]) * exp(-x[1]*solexp)
    new_ird31 = ird31 / ((1d - x[0]) * exp(-x[1]*solexp31)) + ird31[0]
    new_ird32 = ird32 / ((1d - x[0]) * exp(-x[1]*solexp32)) + ird31[0]
    stats = moment((new_ird32 - new_ird31),mdev=out_value)
    ; evaluate the derivative 
    temp0=ird32*exp(x[1]*solexp32) - ird31*exp(x[1]*solexp31)
    temp1=temp0 /(1d - x[0]) - stats[0]
    temp2=ird32*solexp32*exp(x[1]*solexp32) - ird31*solexp31*exp(x[1]*solexp31)
    dfdx0 = TOTAL( temp0 * temp1/abs(temp1) ) / n_elements(ird31) / (1d - x[0])^2.0d
    dfdx1 = TOTAL( temp2 * temp1/abs(temp1) ) / n_elements(ird31) / (1d - x[0])
    df=[dfdx0,dfdx1]
    print,x[0],x[1],out_value
    return, out_value
end

;*******************************************************************
function mytnmin_func7, x, df, solexp31=solexp31, ird31=ird31, solexp32=solexp32, ird32=ird32, $
    new_ird31=new_ird31, new_ird32=new_ird32
    ; measured_irrad = SolarIrrad * (1 - x[0]*solexp)
    ; SolarIrrad = measured_irrad / (1 - x[0]*solexp)
    new_ird31 = ird31 / (1d - x[0]*solexp31)
    new_ird32 = ird32 / (1d - x[0]*solexp32)
    F = new_ird32 - new_ird31
    ; stats[0] is the average
    stats = moment(F,mdev=out_value)
    ; evaluate the derivative 
    signs = dblarr(n_elements(F)) + 1d
    p=where((F-stats[0]) lt 0d,count)
    if count gt 0 then signs[p]=-1d
    temp0=ird32*solexp32 / (1d - x[0]*solexp32)^2d
    temp1=ird31*solexp31 / (1d - x[0]*solexp31)^2d
    df = TOTAL(signs * (temp0 - temp1)) / double(n_elements(ird31))
    print,x[0],df,out_value
    return, out_value
end

;*******************************************************************
function mytnmin_func8, x, df, solexp31=solexp31, ird31=ird31, solexp32=solexp32, ird32=ird32, $
    new_ird31=new_ird31, new_ird32=new_ird32
    ; measured_irrad = SolarIrrad * (1 - x[0]*solexp)
    ; SolarIrrad = measured_irrad / (1 - x[0]*solexp)
    new_ird32 = ird32 / (1d - x[0]*solexp32)
    new_ird31 = ird31 / (1d - x[1]*solexp31)
    F = new_ird32 - new_ird31
    ; stats[0] is the average
    stats = moment(F,mdev=out_value)
    ; evaluate the derivative 
    signs = dblarr(n_elements(F)) + 1d
    p=where((F-stats[0]) lt 0d,count)
    if count gt 0 then signs[p]=-1d
    temp0=ird32*solexp32 / (1d - x[0]*solexp32)^2d
    df0 = TOTAL(signs * temp0) / double(n_elements(ird32))
    temp1=ird31*solexp31 / (1d - x[1]*solexp31)^2d
    df1 = TOTAL(signs * temp1) / double(n_elements(ird31))
    df=[df0,df1]
    print,x,df,out_value
    return, out_value
end

;*******************************************************************
function mytnmin_func9, x, df, solexp31=solexp31, ird31=ird31, solexp32=solexp32, ird32=ird32, $
    new_ird31=new_ird31, new_ird32=new_ird32
    ; SolarIrrad = measured_irrad  - x[0]*solexp)
    new_ird31 = ird31 - x[0]*solexp31
    new_ird32 = ird32 - x[0]*solexp32
    F = new_ird32 - new_ird31
    ; stats[0] is the average
    stats = moment(F,mdev=out_value)
    ; evaluate the derivative 
    signs = dblarr(n_elements(F)) + 1d
    p=where((F-stats[0]) lt 0d,count)
    if count gt 0 then signs[p]=-1d
    df = TOTAL(signs * (solexp31 - solexp32)) / double(n_elements(ird31))
    print,x[0],df,out_value
    return, out_value
end

;*******************************************************************
function mytnmin_func10, x, new_ird31=new_ird31, new_ird32=new_ird32
    common toto, solexp31, ird31, solexp32, ird32
    ; SolarIrrad = measured_irrad  - x[0]*solexp)
    new_ird31 = ird31 - x[0]*solexp31
    new_ird32 = ird32 - x[0]*solexp32
    F = new_ird32 - new_ird31
    ; stats[0] is the average
    stats = moment(F,mdev=out_value)
    return, out_value
end


function mytnmin_dfunc10, x
    common toto, solexp31, ird31, solexp32, ird32
    ; evaluate the derivative 
    new_ird31 = ird31 - x[0]*solexp31
    new_ird32 = ird32 - x[0]*solexp32
    F = new_ird32 - new_ird31
    stats = moment(F,mdev=out_value)
    signs = dblarr(n_elements(F)) + 1d
    p=where((F-stats[0]) lt 0d,count)
    if count gt 0 then signs[p]=-1d
    df = TOTAL(signs * (solexp31 - solexp32)) / double(n_elements(ird31))
    return, df
end

;*******************************************************************
function mytnmin_func11, x, new_ird31=new_ird31, new_ird32=new_ird32
    common toto, solexp31, ird31, solexp32, ird32
    ; SolarIrrad = measured_irrad  - x[0]*solexp)
    new_ird31 = ird31 - x[0]*solexp31
    new_ird32 = ird32 - x[0]*solexp32
    F = new_ird32 - new_ird31
    return, abs(mean(F))
end


function mytnmin_dfunc11, x
    common toto, solexp31, ird31, solexp32, ird32
    ; evaluate the derivative 
    new_ird31 = ird31 - x[0]*solexp31
    new_ird32 = ird32 - x[0]*solexp32
    F = new_ird32 - new_ird31
    val=mean(F)
    if val gt 0 then signs=1d else signs=-1d
    df = signs * TOTAL(solexp31 - solexp32) / double(n_elements(ird31))
    return, df
end

;*******************************************************************
function mytnmin_func12, x, df, solexp31=solexp31, ird31=ird31, solexp32=solexp32, ird32=ird32, $
    new_ird31=new_ird31, new_ird32=new_ird32
    ; SolarIrrad = measured_irrad  - x[0]*solexp)
    new_ird31 = ird31 - x[0]*solexp31
    new_ird32 = ird32 - x[0]*solexp32
    F = new_ird32 - new_ird31
    val=mean(F)
    if val gt 0 then signs=1d else signs=-1d
    df = signs * TOTAL(solexp31 - solexp32) / double(n_elements(ird31))
    return, abs(val)
end


;*******************************************************************
function mytnmin_func13, x, df, solexp31=solexp31, ird31=ird31, solexp32=solexp32, ird32=ird32, $
    new_ird31=new_ird31, new_ird32=new_ird32
    ; measured_irrad = SolarIrrad * exp(-abs(x[0]*solexp + x[1]))
    ; SolarIrrad = measured_irrad / exp(-abs(x[0]*solexp + x[1]))
    ; SolarIrrad = measured_irrad * exp(abs(x[0]*solexp + x[1]))
    new_ird31 = ird31 * exp(abs(x[0]*solexp31 + x[1]))
    new_ird32 = ird32 * exp(abs(x[0]*solexp32 + x[1]))
    F = new_ird32 - new_ird31
    out_value = abs(mean(F))
    ; evaluate the derivative 
    ;temp31 = ird31 * exp(abs(x[0]*solexp31 + x[1])) * (x[0]*solexp31 + x[1]) / abs(x[0]*solexp31 + x[1])
    ;temp32 = ird32 * exp(abs(x[0]*solexp32 + x[1])) * (x[0]*solexp32 + x[1]) / abs(x[0]*solexp32 + x[1])
    ;dfd0 = temp31 * solexp31 - temp32 * solexp32
    ;dfd1 = temp31 - temp32
    ;df=moment(df0, mdev=out0)
    ;df=moment(df1, mdev=out1)
    ;df = [out0,out1]
    print,x[0],x[1],out_value
    return, out_value
end

;*******************************************************************
function mytnmin_func14, x, df, solexp31=solexp31, ird31=ird31, solexp32=solexp32, ird32=ird32, $
    new_ird31=new_ird31, new_ird32=new_ird32
    ; measured_irrad = SolarIrrad * (1 - x[0]) * exp(-abs(x[1]*solexp))
    ; SolarIrrad = measured_irrad / ((1 - x[0]) * exp(-abs(x[0]*solexp)))
    new_ird31 = ird31 / ((1d - x[0]) * exp(abs(x[1]*solexp31)))
    new_ird32 = ird32 / ((1d - x[0]) * exp(abs(x[1]*solexp32)))
    F = new_ird32 - new_ird31
    out_value = abs(mean(F))
    ; evaluate the derivative 
    ;temp31 = ird31 * exp(abs(x[0]*solexp31 + x[1])) * (x[0]*solexp31 + x[1]) / abs(x[0]*solexp31 + x[1])
    ;temp32 = ird32 * exp(abs(x[0]*solexp32 + x[1])) * (x[0]*solexp32 + x[1]) / abs(x[0]*solexp32 + x[1])
    ;dfd0 = temp31 * solexp31 - temp32 * solexp32
    ;dfd1 = temp31 - temp32
    ;df=moment(df0, mdev=out0)
    ;df=moment(df1, mdev=out1)
    ;df = [out0,out1]
    print,x[0],x[1],out_value
    return, out_value
end

;*******************************************************************
function mytnmin_func15, x, df, solexp31=solexp31, ird31=ird31, solexp32=solexp32, ird32=ird32, $
    new_ird31=new_ird31, new_ird32=new_ird32
    ; measured_irrad = SolarIrrad * (exp(-abs(x[0]*solexp)) + exp(-abs(x[1]*solexp))) / 2
    ; SolarIrrad = measured_irrad / ((exp(-abs(x[0]*solexp)) + exp(-abs(x[1]*solexp))) / 2)
    new_ird31 = ird31 / ((exp(-abs(x[0]*solexp31)) + exp(-abs(x[1]*solexp31))) / 2d)
    new_ird32 = ird32 / ((exp(-abs(x[0]*solexp32)) + exp(-abs(x[1]*solexp32))) / 2d)
    F = new_ird32 - new_ird31
    out_value = abs(mean(F))
    ; evaluate the derivative 
    ;temp31 = ird31 * exp(abs(x[0]*solexp31 + x[1])) * (x[0]*solexp31 + x[1]) / abs(x[0]*solexp31 + x[1])
    ;temp32 = ird32 * exp(abs(x[0]*solexp32 + x[1])) * (x[0]*solexp32 + x[1]) / abs(x[0]*solexp32 + x[1])
    ;dfd0 = temp31 * solexp31 - temp32 * solexp32
    ;dfd1 = temp31 - temp32
    ;df=moment(df0, mdev=out0)
    ;df=moment(df1, mdev=out1)
    ;df = [out0,out1]
    print,x[0],x[1],out_value
    return, out_value
end

;*******************************************************************


function tnmin_esrdeg, wavelength, version=version, coeffs=coeffs, status=status, $
    dburl=dburl, dbdriver=dbdriver, user=user, password=password, nfev=nfev, niter=niter, $
    errmsg=errmsg, inspectrum=inspectrum, closest=closest, solar54=solar54, solar55=solar55

    common toto, solarexp31, irrad31, solarexp32, irrad32

    ; look at the ESRFullScans for SimA and SimB that overlap in time
    ; We found 8 periods from start of mission to day 1422 - start with these to test

    ;t0=[278.52850d, 453.02056, 491.08897, 673.03946, 852.04106, 1034.0675, 1223.0246, 1419.0515, 1607.5097, 1790.1583, 1972.0486, 2160.5427]
    ;t1=[279.50538d, 453.99000, 492.07001, 674.01414, 853.01159, 1035.0364, 1223.9937, 1421.1005, 1609.5584, 1792.2063, 1974.1072, 2162.4683]
    t0=[453.02056d, 491.08897d, 673.03946d, 852.04106d, 1034.0675d, 1223.0246d, 1419.0515d]
    t1=[453.99000d, 492.07001d, 674.01414d, 853.01159d, 1035.0364d, 1223.9937d, 1421.1005d]

    ; check if spectrum were provided (prevent from querying db every time)
    if size(inspectrum,/tname) ne 'STRUCT' then $
        inspectrum = {sp31:ptrarr(n_elements(t0)), sp32:ptrarr(n_elements(t0))}

    if n_elements(user) eq 0 then     user='sbeland'
    if n_elements(password) eq 0 then password='sbeland1'
    if n_elements(dburl) eq 0 then    dbUrl='jdbc:sybase:Tds:lasp-db-dev:4100/SORCE_SIM_V19'
    if n_elements(dbdriver) eq 0 then dbDriver='com.sybase.jdbc3.jdbc.SybDriver'

    ; version 15 of development database contains SimCalibratedIrradiance without degradation
    if n_elements(version) eq 0 then version=15

    ; for now only process one wavelength at a time
    if n_elements(wavelength) gt 1 then wavelength=wavelength[0]

    ; get the SimSolarExposureData for modes 31 and 32
    if size(solar54,/tname) ne "STRUCT" then begin
        q1 = 'SELECT orbitStartTimeGps t0, orbitEndTimeGps t1, cumulativeSolarExpHrtOut solar_exp from '+$
             'SimSolarExposureData where instrumentModeId=54'
        query_database, q1, solar54, info
    endif

    if size(solar55,/tname) ne "STRUCT" then begin
        q1 = 'SELECT orbitStartTimeGps t0, orbitEndTimeGps t1, cumulativeSolarExpHrtOut solar_exp from '+$
             'SimSolarExposureData where instrumentModeId=55'
        query_database, q1, solar55, info
    endif


    ; get the corresponding solar_exp and irradiance at the specified
    ; wavelength for all of the spectrum for the listed days
    irrad31=[]
    solarexp31=[]
    time31=[]
    irrad32=[]
    solarexp32=[]
    time32=[]
    query_database, /reset
    for i=0,n_elements(inspectrum.sp31)-1 do begin
        if NOT PTR_VALID(inspectrum.sp31[i]) then begin
            print,'  extracting spectra mode 31 from ',t0[i]
            sp = get_science_product(['Wavelength','SimCalibratedIrradiance'],t0[i],t1[i],/mission,31,$
                version=version,dburl=dburl,dbdriver=dbdriver,user=user,password=password)
            s=sort(sp.wavelengthref)
            sp=sp[s]
            inspectrum.sp31[i] = ptr_new(sp)
        endif
        pmin=min(abs((*inspectrum.sp31[i]).wavelengthref-wavelength),pos)
        if not keyword_set(closest) then begin
            irrad = interpol((*inspectrum.sp31[i]).irradiance, (*inspectrum.sp31[i]).wavelengthref, wavelength,/lsq)
        endif else begin
            closest_wavelength = min(abs((*inspectrum.sp31[i]).wavelengthref - wavelength),pos)
            irrad = ((*inspectrum.sp31[i]).wavelengthref)[pos]
        endelse
        irrad31=[irrad31,irrad]
        ; use the cumulative solar_exposure at the end of the previous orbit
        p=where(solar54.t1 le (*inspectrum.sp31[i])[pos].MICROSECONDSSINCEGPSEPOCH,count)
        solarexp31=[solarexp31, solar54[p[-1]].solar_exp/86400d]
        time31=[time31,(*inspectrum.sp31[i])[pos].MICROSECONDSSINCEGPSEPOCH]

        if NOT PTR_VALID(inspectrum.sp32[i]) then begin
            print,'  extracting spectra mode 32 from ',t0[i]
            sp = get_science_product(['Wavelength','SimCalibratedIrradiance'],t0[i],t1[i],/mission,32,$
                version=version,dburl=dburl,dbdriver=dbdriver,user=user,password=password)
            s=sort(sp.wavelengthref)
            sp=sp[s]
            inspectrum.sp32[i] = ptr_new(sp)
        endif
        pmin=min(abs((*inspectrum.sp32[i]).wavelengthref-wavelength),pos)
        if not keyword_set(closest) then begin
            irrad = interpol((*inspectrum.sp32[i]).irradiance, (*inspectrum.sp32[i]).wavelengthref, wavelength,/lsq)
        endif else begin
            irrad = ((*inspectrum.sp32[i]).wavelengthref)[pos]
        endelse
        irrad32=[irrad32,irrad]
        ; use the cumulative solar_exposure at the end of the previous orbit
        p=where(solar55.t1 le (*inspectrum.sp32[i])[pos].MICROSECONDSSINCEGPSEPOCH,count)
        solarexp32=[solarexp32, solar55[p[-1]].solar_exp/86400d]
        time32=[time32,(*inspectrum.sp32[i])[pos].MICROSECONDSSINCEGPSEPOCH]
    endfor

    ; now we want to minimize the difference in irradiance between 31 and 32
    ; when fitting the same exponential degradation model
    ; convert the cumulative exposure in days
    functargs = {solexp31:solarexp31, ird31:irrad31, solexp32:solarexp32, ird32:irrad32}
    parinfo=replicate({value:0d, fixed:0, step:0.1d, tnside:2}, 2)
    parinfo[*].value=-0.0003d

    ;res=ladfit([solarexp32,solarexp31],[irrad32,irrad31])
    ;print,res
    ;coeff0=[1.8d, -3d-4, 0.1d]
    coeff0=[0.01d, 0.001d]
    print,'  initial guess: ',coeff0
    coeffs = tnmin('mytnmin_func15', coeff0, functargs=functargs, bestmin=f0, status=status, $
             nfev=nfev, niter=niter, autoderivative=1, errmsg=errmsg)
    res=mytnmin_func15(coeffs, df, solexp31=solarexp31, ird31=irrad31, solexp32=solarexp32, ird32=irrad32, $
        new_ird31=new_ird31, new_ird32=new_ird32)

    ;DFPMIN, coeff0, 1.0d-16, fmin, 'mytnmin_func11', 'mytnmin_dfunc10', iter=niter,/double
    ;coeff0=[-3d-3,-3d-3]
    ;powell,coeff0, [[1.0d,1.0d],[1.0d,1.0d]], 1d-16, fmin, 'mytnmin_func11'
    ;coeffs=coeff0
    ;res=amoeba(1d-16,function_name='mytnmin_func10',function_value=coeffs, ncalls=ncalls, nmax=10000, p0=[-3d-4], scale=2d)
    ;res=mytnmin_func11(coeffs, new_ird31=new_ird31, new_ird32=new_ird32)
    print,'  final value: ',coeffs

    ; apply the resulting coefficients and calculate the new CalibratedIrradiance
    ;new_ird31 = irrad31 / (1d + exp(-abs(coeffs[0]*solarexp31 + coeffs[1])) - exp(-abs(coeffs[1])))
    ;new_ird32 = irrad32 / (1d + exp(-abs(coeffs[0]*solarexp32 + coeffs[1])) - exp(-abs(coeffs[1])))
    ;new_ird31 = irrad31 / exp(-abs(coeffs[0]*solarexp31 + coeffs[1]))
    ;new_ird32 = irrad32 / exp(-abs(coeffs[0]*solarexp32 + coeffs[1]))
    n31=n_elements(new_ird31)
    n32=n_elements(new_ird32)
    sp31={MICROSECONDSSINCEGPSEPOCH:time31, wavelength:replicate(wavelength,n31), solarexp:solarexp31, $
        irradiance:irrad31, new_irrad:new_ird31}
    sp32={MICROSECONDSSINCEGPSEPOCH:time32, wavelength:replicate(wavelength,n32), solarexp:solarexp32, $
        irradiance:irrad32, new_irrad:new_ird32}
    return, {sp31:sp31, sp32:sp32}

end
