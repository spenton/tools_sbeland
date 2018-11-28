pro solstice_level3_to_pointers, out9=out9, out11=out11

dbDriver='com.sybase.jdbc3.jdbc.SybDriver'                                                                           
dbUrl='jdbc:sybase:Tds:lasp-db-dev:4100/vanier_db'                                                                   
user='vanier'                                                                                                        
password='vanier1' 


q9="select nominalGpsTimetag,instrumentModeId,version,minWavelength,maxWavelength,irradiance from Level3SolarSpectra s, Level3SolarSpectraData d where s.spectraId=d.spectraId  and s.version=113 and s.instrumentModeId=9 and s.timeSpanInHours=24.0 and s.nominalGpsTimetag>727358413000000. and s.nominalGpsTimetag<1090324816000000."

q11="select nominalGpsTimetag,instrumentModeId,version,minWavelength,maxWavelength,irradiance from Level3SolarSpectra s, Level3SolarSpectraData d where s.spectraId=d.spectraId  and s.version=113 and s.instrumentModeId=11 and s.timeSpanInHours=24.0 and s.nominalGpsTimetag>727358413000000. and s.nominalGpsTimetag<1090324816000000."


query_database,q9,res9, nrows,dbdriver=dbdriver,dburl=dburl,user=user,password=password
query_database,q11,res11, nrows,dbdriver=dbdriver,dburl=dburl,user=user,password=password


;arrange the data by day
hist9 = histogram(res9.NOMINALGPSTIMETAG, bin=86400d6, reverse_indices=rev9)
npts9=n_elements(hist9)
out9=[]
for i=0L,npts9-2L do begin
    if rev9[i] eq rev9[i+1] then continue
    ind = rev9[rev9[i] : rev9[i+1]-1]
    tmp={timejd:gps2jd(res9[ind].NOMINALGPSTIMETAG/1d6), wavelength:(res9[ind].MINWAVELENGTH + res9[ind].MAXWAVELENGTH)/2d, irradiance:res9[ind].IRRADIANCE}
    out9=[out9,ptr_new(tmp)]
endfor


;arrange the data by day
hist11 = histogram(res11.NOMINALGPSTIMETAG, bin=86400d6, reverse_indices=rev11)
npts11=n_elements(hist11)
out11=[]
for i=0L,npts11-2L do begin
    if rev11[i] eq rev11[i+1] then continue
    ind = rev11[rev11[i] : rev11[i+1]-1]
    tmp={timejd:gps2jd(res11[ind].NOMINALGPSTIMETAG/1d6), wavelength:(res11[ind].MINWAVELENGTH + res11[ind].MAXWAVELENGTH)/2d, irradiance:res11[ind].IRRADIANCE}
    out11=[out11,ptr_new(tmp)]
endfor




end
