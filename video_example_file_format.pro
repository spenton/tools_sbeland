PRO VIDEO_EXAMPLE_FILE_FORMAT

width = 500

height = 500

fps = 20

 

surf = SURFACE(/TEST, /BUFFER, DIMENSIONS=[width,height])

 

; Each of the following lines produces a file in

; a different format.

;oVid = IDLffVideoWrite('video_example_file_format.webm')

;oVid = IDLffVideoWrite('video_example_file_format.gif')

;oVid = IDLffVideoWrite('video_example_file_format.swf')

oVid = IDLffVideoWrite('video_example_file_format.bin', FORMAT='mp4')

; Prints out a list of supported file formats

PRINT, "Supported file formats: ", oVid.GetFormats()

vidStream = oVid.AddVideoStream(width, height, fps)

 

FOR i = 0, 90 do begin

  surf.Rotate, 4, /YAXIS

  frame = surf.CopyWindow()

  !NULL = oVid.Put(vidStream, frame)

ENDFOR

 

oVid = 0

END

