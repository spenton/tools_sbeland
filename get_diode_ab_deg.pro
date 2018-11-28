;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Extract the degradation information from the photo-diodes comparing
;   SimA to SimB matching the data for the requested wavelength 
;   taken at the "same" time for SimA and SimB diode.
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
;
; RETURNED PARAMETERS:
;
; OPTIONAL OUTPUT PARAMETERS:
;
; EXAMPLE:  
;
; COMMON BLOCKS:
;      NONE
;
; NOTES:
;    The table used for the irradiance is expected to have the final degradation
;    correction applied but no raypath or diode degradation applied to
;    the photodiode data (as in version 21 of our development database for version 19).
;
;
; REVISION HISTORY:
;   Revision: $Id: get_diode_ab_deg.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;------------------------------------------------------------------
function get_diode_ab_deg, starttime, stoptime, instModeIdA, version=version, $
    gps=gps, missionDays=missionDays, julianDays=julianDays, $
    dburl=dburl, dbdriver=dbdriver, user=user, password=password, $
    inspectA=inspectA, inspectB=inspectB, percent=percent, alignobc=alignobc, $
    irradCol=irradCol, degColumns=degColumns, noplot=noplot, afactor=afactor

    obctimes = sd2gps([0.0d, 441.5, 1570.0d, 2173d, 2455d, 2803d, 5000d])*1d6

    if keyword_set(missionDays) then begin
       ; user specified time in mission (sorce) days
       t0 = sd2gps(startTime)*1.d6
       t1 = sd2gps(stopTime)*1.d6
    endif else if keyword_set(julianDays) then begin
       ; user specified time in julian days
       t0 = jd2gps(startTime)*1.d6
       t1 = jd2gps(stopTime)*1.d6
    endif else begin
       ; user specified timetags in gps microseconds
       t0 = startTime
       t1 = stopTime
    endelse

    if n_elements(version) eq 0 then version=21
    if n_elements(user) eq 0 then     user='sbeland'
    if n_elements(password) eq 0 then password='sbeland1'
    if n_elements(dburl) eq 0 then    dbUrl='jdbc:sybase:Tds:lasp-db-dev:4100/SORCE_SIM_V19'
    if n_elements(dbdriver) eq 0 then dbDriver='com.sybase.jdbc3.jdbc.SybDriver'

    ; limit the wavelength range according to the specified instrumentModeId
    if instModeIdA eq 41 then begin
        instrument=54
        instModeIdB=45
        nsum=3
    endif else if instModeIdA eq 43 then begin
        instrument=54
        instModeIdB=47
        nsum=3
    endif else if instModeIdA eq 44 then begin
        instrument=54
        instModeIdB=48
        nsum=3
    endif else begin
        print,'Error: wrong instrumentModeId (should be between 41 and 48)'
        return,-1
    endelse

    ; the order here is important for the query_database call later
    if size(irradCol,/tname) ne 'STRING' then begin
        irradColumns=['Wavelength','SimCorrectedIrradiance','SimConvertedDataNumbers']
    endif else begin
        irradColumns=['Wavelength',irradCol[0],'SimConvertedDataNumbers']
    endelse
    if keyword_set(degColumns) then irradColumns=[irradColumns,'SimDegradationColumn','SimPrismTransDegradation']

    query_database,/reset

    if size(inspectB,/tname) ne 'STRUCT' then begin
        ; extract the spectrum from SolarQuickScan24 during this time range 
        ; which are closest to the ESR data
        plansB = get_sorce_plan(t0, t1, /gps,instrument='SIMB',activity='SolarQuickScan24')
        nspectB = n_elements(plansB)
        print,'   '+strtrim(string(nspectB),2)+' spectrum to process ...'
        diodeb_data=ptrarr(nspectB)
        diodeb_starttime=dblarr(nspectB)

        ; get all of the SolarQuickScan24 data
        ; using get_science_spectra is very slow here because of the amount of pd data
        ; do a direct query instead
        for pl=0, nspectB-1L do begin
            gps0 = plansb[pl].starttime*1d6
            gps1 = plansb[pl].stoptime*1d6
            if not keyword_set(degColumns) then begin
                q1='SELECT w.microsecondsSinceGpsEpoch,w.instrumentModeId,w.version,w.wavelengthRef,ird.irradiance '+$
                   'FROM '+irradColumns[0]+' w, '+irradColumns[1]+' ird where '+$
                   'w.version='+strtrim(string(version),2)+' and w.version=ird.version and '+$
                   'w.instrumentModeId='+strtrim(string(instModeIdB),2)+' and '+$ 
                   'w.instrumentModeId=ird.instrumentModeId and '+$
                   'w.microsecondsSinceGpsEpoch>='+strtrim(string(ulong64(gps0)),2)+' and '+$
                   'w.microsecondsSinceGpsEpoch<='+strtrim(string(ulong64(gps1)),2)+' and '+$
                   'w.microsecondsSinceGpsEpoch=ird.microsecondsSinceGpsEpoch'
            endif else begin
                q1='SELECT w.microsecondsSinceGpsEpoch,w.instrumentModeId,w.version,w.wavelengthRef,ird.irradiance,'+$
                   'dc.degradationColumn, pd.prismTransDegradation FROM '+irradColumns[0]+' w, '+$
                   irradColumns[1]+' ird, '+irradColumns[3]+' dc, '+irradColumns[4]+' pd where '+$
                   'w.version='+strtrim(string(version),2)+' and w.version=ird.version and '+$
                   'dc.version=ird.version and pd.version=ird.version and '+$
                   'w.instrumentModeId='+strtrim(string(instModeIdB),2)+' and '+$ 
                   'w.instrumentModeId=ird.instrumentModeId and '+$
                   'dc.instrumentModeId=ird.instrumentModeId and pd.instrumentModeId=ird.instrumentModeId and '+$
                   'w.microsecondsSinceGpsEpoch>='+strtrim(string(ulong64(gps0)),2)+' and '+$
                   'w.microsecondsSinceGpsEpoch<='+strtrim(string(ulong64(gps1)),2)+' and '+$
                   'w.microsecondsSinceGpsEpoch=ird.microsecondsSinceGpsEpoch and '+$
                   'dc.microsecondsSinceGpsEpoch=ird.microsecondsSinceGpsEpoch and '+$
                   'pd.microsecondsSinceGpsEpoch=ird.microsecondsSinceGpsEpoch'
            endelse
            query_database, q1, res, user=user,password=password,dburl=dburl,dbdriver=dbdriver
            if size(res,/tname) ne 'STRUCT' then begin
                print,'Error: no spectra found for this plan ',plansb[pl]
                continue
            endif
            s=sort(res.wavelengthref)
            res= res[s]
            if size(afactor,/tname) eq 'STRUCT' and keyword_set(degColumns) then begin
                ; if afactor is provided we expect it to be a structure with wavelength and a_factor tags
                afact=interpol(afactor.a_factor, afactor.wavelength, res.wavelengthref,/lsq)

                ; get the kappa
                query_database, /reset
                q2 = "SELECT calibrationSetId from CalibrationMetadata where calibrationTableName='SimPrismDegKappaCal' "+$
                    "and version=19 and instrumentModeId="+strtrim(string(instrument),2)
                query_database, q2, result, info 
                q3='select x,y from SimPrismDegKappaCal where calibrationSetId='+strtrim(string(result.(0)),2)
                query_database, q3, kappa, info
                kappa.x = 10d ^ kappa.x

                ; calculate the prism degradation
                kappadeg=interpol(kappa.y,kappa.x,res.wavelengthref,/spline)
                degcol=res.DEGRADATIONCOLUMN
                new_prismTransDeg = exp(-kappadeg * degcol) - afact * exp(-kappadeg * degcol) + afact * exp(-kappadeg * degcol / 2d)
                res.irradiance = res.irradiance * res.prismTransDegradation / new_prismTransDeg
            endif
            diodeb_data[pl]= ptr_new(res)
            diodeb_starttime[pl]= res[0].(0)
        endfor
        keep = where(ptr_valid(diodeb_data),count)
        if count eq 0 then begin
            print,'Error:  no valid diode data found'
            return,-1
        endif
        inspectB = {starttime:diodeb_starttime[keep], diode:diodeb_data[keep]}
    endif

    if size(inspectA,/tname) ne 'STRUCT' then begin
        ; extract the spectrum from SolarQuickScan24 during this time range 
        ; which are closest to the ESR data
        plansA = get_sorce_plan(t0, t1, /gps,instrument='SIMA',activity='SolarQuickScan24')
        planA_pos=[]
        for i=0L,n_elements(inspectB.starttime)-1L do begin
            mn=min(abs(inspectB.starttime[i]/1d6 - plansA.starttime),pos)
            ; only keep if diodeA scan is within 2 days of diodeB scan
            if mn lt 2d*86400d then planA_pos=[planA_pos,pos]
        endfor
        ; only keep the unique plans
        s=sort(planA_pos)
        q=uniq(planA_pos[s])
        plansA=plansA[planA_pos[s[q]]]
        nspectA = n_elements(plansA)
        print,'   '+strtrim(string(nspectA),2)+' spectrum to process ...'
        diodeA_data=ptrarr(nspectA)
        diodeA_starttime=dblarr(nspectA)

        ; get all of the SolarQuickScan24 data
        ; using get_science_spectra is very slow here because of the amount of pd data
        ; do a direct query instead
        for pl=0, nspectA-1L do begin
            gps0 = plansA[pl].starttime*1d6
            gps1 = plansA[pl].stoptime*1d6
            if not keyword_set(degColumns) then begin
                q1='SELECT w.microsecondsSinceGpsEpoch,w.instrumentModeId,w.version,w.wavelengthRef,ird.irradiance '+$
                   'FROM '+irradColumns[0]+' w, '+irradColumns[1]+' ird where '+$
                   'w.version='+strtrim(string(version),2)+' and w.version=ird.version and '+$
                   'w.instrumentModeId='+strtrim(string(instModeIdA),2)+' and '+$ 
                   'w.instrumentModeId=ird.instrumentModeId and '+$
                   'w.microsecondsSinceGpsEpoch>='+strtrim(string(ulong64(gps0)),2)+' and '+$
                   'w.microsecondsSinceGpsEpoch<='+strtrim(string(ulong64(gps1)),2)+' and '+$
                   'w.microsecondsSinceGpsEpoch=ird.microsecondsSinceGpsEpoch'
            endif else begin
                q1='SELECT w.microsecondsSinceGpsEpoch,w.instrumentModeId,w.version,w.wavelengthRef,ird.irradiance,'+$
                   'dc.degradationColumn, pd.prismTransDegradation FROM '+irradColumns[0]+' w, '+$
                   irradColumns[1]+' ird, '+irradColumns[3]+' dc, '+irradColumns[4]+' pd where '+$
                   'w.version='+strtrim(string(version),2)+' and w.version=ird.version and '+$
                   'dc.version=ird.version and pd.version=ird.version and '+$
                   'w.instrumentModeId='+strtrim(string(instModeIdA),2)+' and '+$ 
                   'w.instrumentModeId=ird.instrumentModeId and '+$
                   'dc.instrumentModeId=ird.instrumentModeId and pd.instrumentModeId=ird.instrumentModeId and '+$
                   'w.microsecondsSinceGpsEpoch>='+strtrim(string(ulong64(gps0)),2)+' and '+$
                   'w.microsecondsSinceGpsEpoch<='+strtrim(string(ulong64(gps1)),2)+' and '+$
                   'w.microsecondsSinceGpsEpoch=ird.microsecondsSinceGpsEpoch and '+$
                   'dc.microsecondsSinceGpsEpoch=ird.microsecondsSinceGpsEpoch and '+$
                   'pd.microsecondsSinceGpsEpoch=ird.microsecondsSinceGpsEpoch'
            endelse
            query_database, q1, res, user=user,password=password,dburl=dburl,dbdriver=dbdriver
            if size(res,/tname) ne 'STRUCT' then begin
                print,'Error: no spectra found for this plan ',plansA[pl]
                continue
            endif
            s=sort(res.wavelengthref)
            diodeA_data[pl]= ptr_new(res[s])
            if size(afactor,/tname) eq 'STRUCT' and keyword_set(degColumns) then begin
                ; if afactor is provided we expect it to be a structure with wavelength and a_factor tags
                afact=interpol(afactor.a_factor, afactor.wavelength, res.wavelengthref,/lsq)

                ; get the kappa
                if size(kappa,/tname) ne 'STRUCT' then begin
                    query_database, /reset
                    q2 = "SELECT calibrationSetId from CalibrationMetadata where calibrationTableName='SimPrismDegKappaCal' "+$
                        "and version=19 and instrumentModeId="+strtrim(string(instrument),2)
                    query_database, q2, result, info
                    q3='select x,y from SimPrismDegKappaCal where calibrationSetId='+strtrim(string(result.(0)),2)
                    query_database, q3, kappa, info
                    kappa.x = 10d ^ kappa.x
                endif

                ; calculate the prism degradation
                kappadeg=interpol(kappa.y,kappa.x,res.wavelengthref,/spline)
                degcol=res.DEGRADATIONCOLUMN
                new_prismTransDeg = exp(-kappadeg * degcol) - afact * exp(-kappadeg * degcol) + afact * exp(-kappadeg * degcol / 2d)
                res.irradiance = res.irradiance * res.prismTransDegradation / new_prismTransDeg
            endif

            diodeA_starttime[pl]= res[0].(0)
        endfor
        inspectA = {starttime:diodeA_starttime, diode:diodeA_data}
    endif

    nspect = n_elements(inspectA.diode)
    nwaves=24.0
    wavestep = ceil(n_elements((*inspectA.diode[0]).wavelengthref) / nwaves)
    wavepos = (lindgen(nwaves)*wavestep+1L) < (n_elements((*inspectA.diode[0]).wavelengthref)-1L)
    waves = (*inspectA.diode[0])[wavepos].wavelengthref
    outdata = ptrarr(nwaves)

    ; for each wavelength
    for w=0,n_elements(waves)-1 do begin
        print,'processing '+strtrim(string(waves[w],format='(F10.1)'),2)+' nm'

        ; loop through every esr value for this prismpos
        pda_timetag=[]
        pdb_timetag=[]
        pda_irrad=[]
        pdb_irrad=[]
        delta_irrad=[]
        pda_degCol=[]
        pdb_degCol=[]
        pda_prismDeg=[]
        pdb_prismDeg=[]
        for i=0L,nspect-1 do begin
            if NOT (ptr_valid(inspectA.diode[i]) and ptr_valid(inspectA.diode[i])) then begin
                print,'   invalid pointer at '+strtrim(string(waves[w]),2)
                continue
            endif
            irrad=interpol((*inspectA.diode[i]).irradiance, (*inspectA.diode[i]).wavelengthref, waves[w], /spline)
            pda_irrad = [pda_irrad, irrad]
            irrad=interpol((*inspectB.diode[i]).irradiance, (*inspectB.diode[i]).wavelengthref, waves[w], /spline)
            pdb_irrad = [pdb_irrad, irrad]
            delta_irrad = [delta_irrad, pda_irrad[-1]-pdb_irrad[-1]]
            mn=min(abs((*inspectA.diode[i]).wavelengthref - waves[w]),posa)
            mn=min(abs((*inspectB.diode[i]).wavelengthref - waves[w]),posb)
            pda_timetag=[pda_timetag,(*inspectA.diode[i])[posa].(0)]
            pdb_timetag=[pdb_timetag,(*inspectB.diode[i])[posb].(0)]
            if keyword_set(degColumns) then begin
                pda_degCol=[pda_degCol, (*inspectA.diode[i])[posa].degradationColumn]
                pdb_degCol=[pdb_degCol, (*inspectB.diode[i])[posb].degradationColumn]
                pda_prismDeg=[pda_prismDeg, (*inspectA.diode[i])[posa].prismTransDegradation]
                pdb_prismDeg=[pdb_prismDeg, (*inspectB.diode[i])[posb].prismTransDegradation]
            endif
        endfor
        if keyword_set(degColumns) then begin
            outdata[w] = ptr_new({pda_timetag:pda_timetag, pdb_timetag:pdb_timetag, $
                wavelength:waves[w], delta_irrad:delta_irrad, $
                pda_irrad:pda_irrad, pdb_irrad:pdb_irrad, $
                pda_degcol:pda_degcol, pdb_degcol:pdb_degcol, $
                pda_prismDeg:pda_prismDeg, pdb_prismDeg:pdb_prismDeg})
        endif else begin
            outdata[w] = ptr_new({pda_timetag:pda_timetag, pdb_timetag:pdb_timetag, $
                wavelength:waves[w], delta_irrad:delta_irrad, $
                pda_irrad:pda_irrad, pdb_irrad:pdb_irrad})
        endelse

        if keyword_set(alignobc) then begin
            for i=0,n_elements(obctimes)-3 do begin
                pleft = where(pda_timetag ge obctimes[i] and pda_timetag lt obctimes[i+1],count)
                if count eq 0 then continue
                pright = where(pda_timetag ge obctimes[i+1] and pda_timetag lt obctimes[i+2],count)
                if count eq 0 then continue
                p0=n_elements(pleft)-21 > 0
                cleft=robust_poly_fit(pda_timetag[pleft[p0:-2]], delta_irrad[pleft[p0:-2]], 1, /double)
                p1=n_elements(pright)-1 <20 
                cright=robust_poly_fit(pda_timetag[pright[1:p1]], delta_irrad[pright[1:p1]], 1, /double)
                v0=poly(pda_timetag[pright[0]],cleft)
                v1=poly(pda_timetag[pright[0]],cright)
                delta_irrad[pright] += (v0-v1)
            endfor
            (*outdata[w]).delta_irrad = delta_irrad
        endif

        if not keyword_set(noplot) then begin
            title='ModeId '+strtrim(string(instModeIdA,format='(I0)'),2)+' - '+strtrim(string(instModeIdB,format='(I0)'),2)+$
                ' @ '+strtrim(string(waves[w],format='(F10.1)'),2)
            resistant_mean,delta_irrad,5.0,mean,goodvec=k
            if keyword_set(percent) then begin
                lineplot,gps2sd(pda_timetag[k]/1d6),delta_irrad[k]/pda_irrad[k]*100d,title=title,xtitle='Mission Day',$
                    ytitle='Delta Irradiance (Percent)', ptitle='Differences in Irradiance between '+$
                    strtrim(string(instModeIdA,format='(I0)'),2)+' and '+strtrim(string(instModeIdB,format='(I0)'),2), nsum=nsum
            endif else begin
                lineplot,gps2sd(pda_timetag[k]/1d6),delta_irrad[k],title=title,xtitle='Mission Day',ytitle='Delta Irradiance',$
                    ptitle='Differences in Irradiance between '+strtrim(string(instModeIdA,format='(I0)'),2)+' and '+$
                    strtrim(string(instModeIdB,format='(I0)'),2),nsum=nsum
            endelse

        endif

    endfor  ; loop for each ESR wavelength

    if not keyword_set(noplot) then begin
        ; plot the corresponding surface in percent changes in irradiance
        xout=gps2sd((*outdata[0]).(0)/1d6)
        yout=[]
        x=[]
        y=[]
        z=[]
        for i=0,nwaves-1 do begin
            resistant_mean,(*outdata[i]).delta_irrad,3.0,mean,goodvec=k
            x=[x,gps2sd((*outdata[i]).(0)[k]/1d6)]
            y=[y,(*outdata[i]).wavelength[k]]
            if keyword_set(percent) then $
                z=[z,((*outdata[i]).delta_irrad[k])/((*outdata[i]).pda_irrad[k])*100d] $
            else $
                z=[z,(*outdata[i]).delta_irrad[k]] 
            yout=[yout,median((*outdata[i]).wavelength[k])]
        endfor
        grid=griddata(x,y,z,smooth=3,xout=xout,yout=yout,/grid)
        if keyword_set(percent) then ztitle='Percent Delta Irradiance' else ztitle='Delta Irradiance'
        title=strtrim(string(instModeIdA,format='(I0)'),2)+' vs '+strtrim(string(instModeIdB,format='(I0)'),2)
        isurface,grid,xout,yout,xtitle='Mission Day',ytitle='Wavelength (nm)',ztitle=ztitle, view_title=title 
    endif



    return,outdata

end
