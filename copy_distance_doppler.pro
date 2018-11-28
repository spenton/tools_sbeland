;+
; NAME:   COPY_DISTANCE_DOPPLER
;
; PURPOSE: 
;   This routine copies the data from the SolarDistAndDopplerFixedStep2
;   in the SORCE_NEW databasere to the SORCE database (since the SORCE
;   database table is not being populated (as of today) but some tools
;   still rely on the data being avaialable there.
;
; CALLING SEQUENCE:
;   COPY_DISTANCE_DOPPLER, version=version, database=database
;                      
; OPTIONAL INPUT KEYWORDS:
;   version -
;      Version number to extract from the SORCE_NEW database.
;   database - 
;      Name of database to get the updated data (default SORCE_NEW).
;
; OPTIONAL OUTPUT PARAMETERS:
;      NONE
;
; OPTIONAL OUTPUT PARAMETERS:
;      NONE
;
; EXAMPLE:  
;
; COMMON BLOCKS:
;      NONE
;
; REVISION HISTORY:
;   2013.12.17  SBeland
;   Revision: $Id
;-
;*****************************************************************************
PRO COPY_DISTANCE_DOPPLER, version=version, database=database, verbose=verbose, dbinsert=dbinsert


    if n_elements(version) eq 0 then version=107
    if n_elements(database) eq 0 then database='SORCE_NEW'

    tables=['SolarDistAndDopplerFixedStep2', 'SpacecraftGeoreferenceData', 'CurrentTarget', 'TargetViewParameters']

    for i=0,n_elements(tables)-1 do begin
        ; get the latest entry in the SORCE database
        q1='select top 1 * from '+tables[i]+' where instrumentModeId=50 and version=6 order by microsecondsSinceGpsEpoch desc'
        query_database,/reset
        query_database,q1,v6_data
        if n_elements(v6_data) eq 0 then begin
           print,'No data found in '+tables[i]+' for verion=6'
           continue
        endif
        t0=v6_data[0].microsecondsSinceGpsEpoch

        ; get the missing data from SORCE_NEW
        q2='select * from '+tables[i]+' where instrumentModeId=50 and version='+strtrim(string(version),2)
        q2=q2+' and microsecondsSinceGpsEpoch>='+strtrim(string(ulong64(t0)),2)+' order by microsecondsSinceGpsEpoch'
        query_database,q2,v107_data, database=database
        nrows=n_elements(v107_data)
        if nrows eq 0 then begin
            print,'No rows to copy for '+tables[i]
            continue
        endif
        ; update the version
        v107_data.version=6

        ; get the list of column names
        q3="select b.name from sysobjects a, syscolumns b where a.id = b.id and a.name='"+tables[i]+"'"
        query_database, q3, names
        ncols=n_elements(names)
        
        ; insert one row at a time
        q4 = "INSERT INTO "+tables[i] +'('
        for j=0,ncols-2 do q4=q4+names[j].name+', '
        q4=q4+names[-1].name+') VALUES('
        
        print,strtrim(string(n_elements(v107_data)),2)+' rows to insert ...'
        stmt17 = fjava_get_jdbc_statement(user='sorce', password='sorcedb', server='sorce-db', database="SORCE")
        conn17 = OBJ_NEW("oJava_DbExchange", stmt17)

        for j=0L,nrows-1 do begin
            ; we assume all the columns are numbers and let the insert handle the type conversion
            str=''
            for k=0,ncols-2 do begin
                str = str+strtrim(string(v107_data[j].(k),format='(G)'),2)+', '
            endfor
            str = str+strtrim(string(v107_data[j].(ncols-1),format='(G)'),2)+') '
            if keyword_set(verbose) then print,string(j+1)+'  '+q4+str
            if keyword_set(dbinsert) then result=conn17->execute(q4+str)
        endfor

    endfor

END
