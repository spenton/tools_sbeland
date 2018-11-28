;+
; NAME: 
;     READ_TSIS_SIM_SRF()
; PURPOSE:               
;     Read in an ascii file with the TSIS SIM SRF log data.
; EXPLANATION:              
;     The routine is somewhat general in that it will create a structure
;     with tagnames found in the file.  It expects to have in the first line
;     the list of global settings, and their corresponding values in the 
;     second line.  The third line contains the variable names each column
;     of data that follows (expecting the data to start on line 4).
; CALLING SEQUENCE: 
;     IDL> data = read_tsis_sim_srf([filename])
;    
; OPTIONAL INPUTS:
;     filename  -  the file to read. If empty, user gets prompted for one
;               
; RETURNS
;     data  - structure containing the global parameters and the data arrays.

; EXAMPLES:
;     IDL> data=read_tsis_sim_srf('TSIS_SIM_SRF_03-15-12_07-31.csv',header=hdr, header_only=header_only)
;     IDL> help,data,/st
;
; REVISION HISTORY:
;          Version 1, SBeland, April 26, 2012
;
; DEPENDENCIES:
;    Uses the ASTRON library routines readcol.pro, valid_num.pro and gettok.pro
;    The ASTRON library can be found at: http://idlastro.gsfc.nasa.gov/
;
;-            

function read_tsis_sim_srf, filename, hdr=header, header_only=header_only

    if n_params() eq 0 then begin
        filename = dialog_pickfile(title='Select a TSIS csv file to read',/must_exist,/read,filter='*.csv')
        ; check if the dialog was cancelled
        if strlen(filename) eq 0 then return,-1
    endif

    ; make sure the filname exists
    result = FILE_SEARCH(filename,/test_read)
    if strlen(result) eq 0 then begin
        print,'Error: '+filename+' does not exist'
        return,-1
    endif

    ; read the header
    hdr1=''
    hdr2='' 
    openr,unit,filename,/get_lun
    readf,unit,hdr1
    readf,unit,hdr2
    header=[hdr1,hdr2]

    ; determine if we have a Global_log file
    if strpos(hdr1,'/') le 5 then begin
        ; we have a Gloabl_Log file starting with the comments on the first line
        ; and the column lables on the second line
        tags=[]
        hdr_struct=[]
        repeat tags=[tags,gettok(hdr2,',')] until strlen(hdr2) eq 0
        for i=0,n_elements(tags)-1 do begin
            pos=strpos(tags[i],'(')
            if pos ge 0 then tags[i] = strtrim(strmid(tags[i],0,pos),2)
            repeat begin
                pos=strpos(tags[i],' ')
                if pos ge 0 then tags[i] = strmid(tags[i],0,pos)+'_'+strmid(tags[i],pos+1)
            endrep until pos eq -1
            hdr_struct=create_struct(hdr_struct, tags[i], 0.0d)
        endfor
        nlines = file_lines(filename) -2L

    endif else begin
        tags=strsplit(hdr1,',',/extract,count=count)
        values=[]
        repeat values=[values,gettok(hdr2,',')] until strlen(hdr2) eq 0
        valid=valid_num(values,results)

        ; create the output structure
        hdr_struct=[]
        for i=0,n_elements(tags)-1 do begin
            ; insert the double result if if was a valid number
            if valid[i] then $
                hdr_struct=create_struct(hdr_struct, tags[i], results[i]) $
            else $
                hdr_struct=create_struct(hdr_struct, tags[i], values[i])
        endfor
        hdr3='' 
        readf,unit,hdr3
        tg=strsplit(hdr3,',',/extract,count=count)
        nlines = file_lines(filename) -3L
    endelse
    outdata = replicate(hdr_struct,nlines)
    indata=dblarr(n_elements(tags),nlines)
    readf,unit,indata
    free_lun,unit

    for i=0,n_elements(tags)-1 do begin
        outdata[*].(i) = reform(indata[i,*])
    endfor
    return, outdata

end
