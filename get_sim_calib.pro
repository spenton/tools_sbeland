;+
; Author: Stephane Beland
;
; PURPOSE: 
;   Return a structure containing information about the calibration
;   information in the database used for the calculation of SIM science
;   data numbers.
;
; CALLING SEQUENCE:
;   result = GET_SIM_CALIB(instrumentModeId, dbTable=dbTable, version=version)
;
; INPUT PARAMETERS:
;   instrumentModeId -
;      The instrument mode of interest:
;      31	SIM_A	ESR
;      32	SIM_B	ESR
;      41	SIM_A	VIS1
;      42	SIM_A	VIS2
;      43	SIM_A	UV
;      44	SIM_A	IR
;      45	SIM_B	VIS1
;      46	SIM_B	VIS2
;      47	SIM_B	UV
;      48	SIM_B	IR
;      54	SIM_A	GLOBAL
;      55	SIM_B	GLOBAL
;
; OPTIONAL INPUT PARAMETERS:
;   dbTable - 
;      The name of the database tabel to get the information from.
;      If not specified, the program will look at all entries in the
;      table CalibrationMetadata for the specified instrumentModId and
;      offer the list to the user to select from.
;   version -
;      Specifies the exact version of the calibration data to extract. 
;      If not specified, will offer a list to the user obtained from 
;      the table CalibrationMetadata.
;   latest - 
;      If specified, will automatically extract the calibration data 
;      from the latest version found in the table CalibrationMetadata.
;
; OPTIONAL OUTPUT PARAMETERS:
;   nrows - 
;      Returns the number of rows found.
;   dbTable -
;      Returns the database table selected if an undefined variable was provided.
;   version -
;      Returns the version selected if an undefined variable was provided.
;
; RETURNED PARAMETERS:
;   A structure whose fields correspond to the returned
;   database column names. 
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
;   Revision: $Id: get_sim_calib.pro,v 1.1 2018/11/26 16:13:48 spenton Exp $
;-
;
function get_sim_calib, instrumentModeId, dbTable=dbTable, version=version, $
                        latest=latest, nrows=nrows

    nrows=0
    ; validate the instrumentModeId
    if n_elements(instrumentModeId) eq 0 then begin
        doc_library,'get_sim_calib'
        return,-1
    endif else begin
        ; verify we got the right instrument
        valid_inst = [31,32,41,42,43,44,45,46,47,48,54,55]
        instrumentModeId = fix(instrumentModeId)
        if (where(valid_inst eq instrumentModeId))[0] lt 0 then  begin
            doc_library,'get_sim_calib'
            return,-1
        endif
    endelse

    ; prepare the database connection
    jstmt = fjava_get_jdbc_statement(user=name, password=password, dburl=dbURL, $
        dbdriver=dbdriver, server=server, database="SORCE")
    oJavaDbExchange = OBJ_NEW("oJava_DbExchange", jstmt)

    if n_elements(dbTable) eq 0 then begin
        ; get list of tables from CalibrationMetadata for this instrumentModeId
        query1 = 'SELECT distinct calibrationTableName FROM SORCE..CalibrationMetadata '+$
                 'WHERE instrumentModeId='+strtrim(instrumentModeId,2)
        res=oJavaDbExchange->getAllValues(query1)
        ; prompt the user for the table to use
        print,'List of database Tables :'
        for i=0,n_elements(res)-1 do print,res[i].calibrationTableName,' ['+strtrim(i+1,2)+']'
        dbnum=''
        read,dbnum,prompt='Type the number corresponding to the desired Table (0 to exit): '
        if strlen(dbnum) eq 0 then return,-1
        dbnum=fix(dbnum)
        if dbnum le 0 or dbnum gt n_elements(res) then return,-1
        dbTable = res[dbnum-1].calibrationTableName
    endif

    if n_elements(version) eq 0 then begin
        ; get the list of available versions from CalibrationMetaData for instModeId and dbTable
        query1 = "SELECT version FROM SORCE..CalibrationMetadata WHERE "+$
                 "instrumentModeId="+strtrim(instrumentModeId,2)+$
                 " AND calibrationTableName='"+dbTable+"' "+$
                 "ORDER by version ASC"
        res=oJavaDbExchange->getAllValues(query1)
        if keyword_set(latest) or n_elements(res) eq 1 then begin
            version = max(res.version)
        endif else begin
            print,'List of available Versions:'
            print,strtrim(res.version,2)
            vnum=''
            read,vnum,prompt='Type the desired Version: '
            if strlen(vnum) eq 0 then return,-1
            version=fix(vnum)
            p=where(version eq res.version,count)
            if count eq 0 then begin
                print,''
                print,'you typed an invalid version number - quitting'
                return,-1
            endif
        endelse
    endif

    ; extract the information for the specified instrumentModeId and version
    query1 = 'SELECT t.* FROM SORCE..'+dbTable+' t, SORCE..CalibrationMetadata m '+$
             'WHERE t.calibrationSetId = m.calibrationSetId '+$
             'AND m.instrumentModeId='+strtrim(instrumentModeId,2)+' '+$
             'AND m.version='+strtrim(version,2)

    res=oJavaDbExchange->getAllValues(query1)
    if size(res,/tname) ne 'STRUCT' then return,-1

    nrows=n_elements(res)
    return, res

end
