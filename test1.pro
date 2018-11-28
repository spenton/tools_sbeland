restore,'~/SORCE/data/sima_vis_uncorr_2011.sav'
restore,'~/SORCE/data/simb_vis_uncorr_2011.sav'
restore,'~/SORCE/data/everyOrbitExpos_AB.sav'
readcol,'~/SORCE/data/overlap_vis1_smooth.txt',afact_a,afact_w,format='(d,d)'
afact={wavelength:afact_w, a:afact_a}

ww=488d

nspecta=n_elements(visa_uncorr_2011.spect20)
nspectb=n_elements(visb_uncorr_2011.spect20)

visa488=dblarr(nspecta)
sda488=dblarr(nspecta)

visb488=dblarr(nspectb)
sdb488=dblarr(nspectb)

for i=0,nspecta-1 do begin 
    s=sort(VISA_UNCORR_2011.spect20[i].wavelength) 
    visa488[i]=interpol(VISA_UNCORR_2011.spect20[i].irradiance[s], VISA_UNCORR_2011.spect20[i].wavelength[s], ww, /spline) 
    mn=min(abs(VISA_UNCORR_2011.spect20[i].wavelength-ww),pos) 
    sda488[i]=VISA_UNCORR_2011.spect20[i].(0)[pos] 
endfor

for i=0,nspectb-1 do begin 
    s=sort(VISB_UNCORR_2011.spect20[i].wavelength) 
    visb488[i]=interpol(VISB_UNCORR_2011.spect20[i].irradiance[s], VISB_UNCORR_2011.spect20[i].wavelength[s], ww, /spline) 
    mn=min(abs(VISB_UNCORR_2011.spect20[i].wavelength-ww),pos) 
    sdb488[i]=VISB_UNCORR_2011.spect20[i].(0)[pos] 
endfor

; average the data in 3 days increment
irradA=visa488[0]
sda=sda488[0]
irradB=visb488[0]
sdb=sdb488[0]
avgdays=3d
for i=1L,n_elements(sdb488)-1L do begin
    p=where(abs(sdb488[i] - sdb488) le avgdays*86400d6,count)
    if count eq 1 then begin
        tmpird=visb488[p]
        tmpsd=sdb488[p]
    endif else if count gt 1 then begin
        resistant_mean, visb488[p], 3.0, tmpird
        tmpsd=mean(sdb488[p])
    endif
    ; avoid duplicate
    if tmpsd eq sdb[-1] then continue
    sdb=[sdb,tmpsd]
    irradB=[irradB,tmpird]
endfor

for i=1L,n_elements(sda488)-1L do begin
    p=where(abs(sda488[i] - sda488) le avgdays*86400d6,count)
    if count eq 1 then begin
        tmpird=visa488[p]
        tmpsd=sda488[p]
    endif else if count gt 1 then begin
        resistant_mean, visa488[p], 3.0, tmpird
        tmpsd=mean(sda488[p])
    endif
    ; avoid duplicate
    if tmpsd eq sdb[-1] then continue
    sda=[sda,tmpsd]
    irradA=[irradA,tmpird]
endfor

sda=gps2sd(sda/1d6)
sdb=gps2sd(sdb/1d6)

stop

;interpolate the irradiance of SimA to match times of SimB
irradA = interpol(irradA, sda, sdb)
sda=sdb

; now get the Kappa*F (Tau) over this time period with a span of ~90days
step=90d


end
