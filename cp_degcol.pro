;+
;   CP_DEGCOL
;
; A quick hack to copy the content of a version of SimPrismDegColCalWave2D
; to another.
; This was necessary for version 1 and version 2 for all modes to leave version 1
; intact for possible re-processing of data release version 21. The version 2 in
; SimPrismDegColCalWave2D were fisrt insert as a new calibration in preparation 
; data release 22 but we decided to simply extend the record for version 1 instead.
;-
pro cp_degcol, instmode, fromver, tover, dbinsert=dbinsert

    ; get the calibrationSetId corresponding to instmode and fromver
    q1='select * from CalibrationMetadata where calibrationTableName="SimPrismDegColCalWave2D" and '
    q1+='instrumentModeId='+strtrim(string(instmode),2)+' and version='+strtrim(string(fromver),2)
    query_database,q1,from_calib
    if size(from_calib,/tname) ne 'STRUCT' then return
    from_calib=from_calib.calibrationsetid

    ; get the calibrationSetId corresponding to instmode and tover
    q2='select * from CalibrationMetadata where calibrationTableName="SimPrismDegColCalWave2D" and '
    q2+='instrumentModeId='+strtrim(string(instmode),2)+' and version='+strtrim(string(tover),2)
    query_database,q2,to_calib
    if size(to_calib,/tname) ne 'STRUCT' then return
    to_calib=to_calib.calibrationsetid

    print,'copying from ',from_calib,' to ',to_calib

    ; extract the data in the fromver
    q3='select * from SimPrismDegColCalWave2D where calibrationSetId='+strtrim(string(from_calib),2)
    query_database,q3,from_data

    ; delete the data in the to_calib
    q4 = "DELETE FROM SimPrismDegColCalWave2D WHERE calibrationSetId="+strtrim(string(to_calib),2)
    if keyword_set(dbinsert) then begin
        print,'deleting entries in ',to_calib
        query_database,q4,res
    endif

    ; create the BCP file to ingest 
    nrows=n_elements(from_data)
    filename='~/SORCE/data/sim_degcol_m'+strtrim(string(instmode),2)+'_'+strtrim(string(from_calib),2)+$
        '_to_'+strtrim(string(to_calib),2)+'.bcp'
    print,'Writing ',nrows,' in file ',filename

    fmt='(I0,",",I0,",",G0.9,",",G0.15)'
    forprint,replicate(to_calib,nrows), ulong64(from_data.orbitEndGpsMicroseconds), from_data.wavelength, $
        from_data.degradationColumn, format=fmt, text=filename

    print,''
    print,'Done'
end

