function simb_ccd, outstr=outstr

    restore,'~/SORCE/data/SIMB_UV_OFFSET.sav'

    ; group data by bins of 10 days and average the offset within each group
    p=where(uvboffset ne -1.0)
    bs1=bs1[p]
    bs2=bs1[p]
    uvboffset=uvboffset[p]
    p0=0L
    p1=n_elements(bs1)-1L
    val_uv=[]
    uvbs1=[]
    uvbs2=[]
    
    REPEAT BEGIN
        p=where(abs(bs1[p0]-bs1) lt 10.0,comp=cp,count)
        if count gt 0 then begin
           val_uv=[val_uv,mean(uvboffset[p])]
           uvbs1=[uvbs1,min([bs1[p],bs2[p]])]
           uvbs2=[uvbs2,max([bs1[p],bs2[p]])]
        endif

        if cp[0] ge 0 then begin
            bs1=bs1[cp]
            bs2=bs2[cp]
            uvboffset=uvboffset[cp]
        endif

    ENDREP UNTIL cp[0] lt 0

    restore,'~/SORCE/data/SIMB_VIS_OFFSET.sav'

    ; group data by bins of 10 days and average the offset within each group
    p=where(visboffset ne -1.0)
    bs1=bs1[p]
    bs2=bs1[p]
    visboffset=visboffset[p]
    p0=0L
    p1=n_elements(bs1)-1L
    val_vis=[]
    visbs1=[]
    visbs2=[]

    REPEAT BEGIN
        p=where(abs(bs1[p0]-bs1) lt 10.0,comp=cp,count)
        if count gt 0 then begin
           val_vis=[val_vis,mean(visboffset[p])]
           visbs1=[visbs1,min([bs1[p],bs2[p]])]
           visbs2=[visbs2,max([bs1[p],bs2[p]])]
        endif

        if cp[0] ge 0 then begin
            bs1=bs1[cp]
            bs2=bs2[cp]
            visboffset=visboffset[cp]
        endif

    ENDREP UNTIL cp[0] lt 0

    uvstr=[]
    for i=0,n_elements(uvbs1)-1 do begin
        str1 = " '2, 47, "+'"Oct 12 2012 00:00", "'+jd2syb(sd2jd(uvbs1[i]))+'", '+$
               strtrim(string(ulong64(sd2gps(uvbs1[i])*1d6)),2)+', 550, '+$
               strtrim(string(val_uv[i],format='(F0.2)'),2)+", 0.0, 0.0, 0.0', $"
        uvstr=[uvstr,str1]
        str1 = " '2, 47, "+'"Oct 12 2012 00:00", "'+jd2syb(sd2jd(uvbs2[i]))+'", '+$
               strtrim(string(ulong64(sd2gps(uvbs2[i])*1d6)),2)+', 550, '+$
               strtrim(string(val_uv[i],format='(F0.2)'),2)+", 0.0, 0.0, 0.0', $"
        uvstr=[uvstr,str1]
    endfor

    visstr=[]
    for i=0,n_elements(visbs1)-1 do begin
        str1 = " '2, 45, "+'"Oct 12 2012 00:00", "'+jd2syb(sd2jd(visbs1[i]))+'", '+$
               strtrim(string(ulong64(sd2gps(visbs1[i])*1d6)),2)+', 550, '+$
               strtrim(string(val_vis[i],format='(F0.2)'),2)+", 0.0, 0.0, 0.0', $"
        visstr=[visstr,str1]
        str1 = " '2, 45, "+'"Oct 12 2012 00:00", "'+jd2syb(sd2jd(visbs2[i]))+'", '+$
               strtrim(string(ulong64(sd2gps(visbs2[i])*1d6)),2)+', 550, '+$
               strtrim(string(val_vis[i],format='(F0.2)'),2)+", 0.0, 0.0, 0.0', $"
        visstr=[visstr,str1]
    endfor

    outdata={visbs1:visbs1, visbs2:visbs2, val_vis:val_vis, uvbs1:uvbs1, uvbs2:uvbs2, val_uv:val_uv}
    outstr = {uvstr:uvstr,visstr:visstr}

    return,outdata

end
