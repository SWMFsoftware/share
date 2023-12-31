; User only needs to set the directory name
; either in SWQU format, or the
; full location
; calling - make_WL_compare,CommonDir='${DIRPATH}/'
; For eg, .r make_WL_compare
; make_WL_compare,CommonDir='/Users/name/Desktop/Results/'
; Last Update:
; NishthaS Dec 29,2023

function set_inst,file
  tmp=strsplit(file,'/',/extract,count=count)
  tmp1 = strsplit(tmp[count-1],'_',/extract) ;last one
  SatName = tmp1[1]
  InstName = tmp1[2]            ; INSTRUMENT NAME
  ; outfile,colortable,occulter,size,range
  if SatName eq 'soho' then begin
     if InstName eq 'c2' then begin
        inst_file = ['LASCO_C2_SIM_','C2','2.2','12','6']
     endif
     if InstName eq 'c3' then begin
        inst_file = ['LASCO_C3_SIM_','C3','4.0','64.','20.']
     endif
  endif
  if SatName eq 'sta' then begin
     if InstName eq 'cor1' then begin
        inst_file = ['STEREOA_COR1_SIM_','COR1A',$
                     '1.6','8','4']
     endif
     if InstName eq 'cor2' then begin
        inst_file = ['STEREOA_COR2_SIM_','COR2A',$,
                     '4.0','30','15']
     endif
  endif
  if SatName eq 'stb' then begin
     if InstName eq 'cor1' then begin
        inst_file = ['STEREOB_COR1_SIM_','COR1B',$
                     '1.6','8','4']
     endif
     if InstName eq 'cor2' then begin
        inst_file = ['STEREOB_COR2_SIM_','COR2B',$
                     '4.0','30','15']
     endif
  endif
  return,inst_file
end
; Returns all details pertaining to the instrument

function read_swmf_wl,file      ; for all files except background
  print,'Reading File = ',file
  restore,'los_template.sav'
  data = read_ascii(file,template=template)
  xx=data.field1
  yy=data.field2
  wl=data.field3
  wl_image=fltarr(300,300)
  nn=long(0)
  for i=0, 299 do begin
     for j=0, 299 do begin
        wl_image[i,j]=wl[nn]
        nn=nn+1
     endfor
  endfor
  return,wl_image
end

function fk_corona,file  ; For background only
  print,'FK corona for background = ',file
  restore,'los_template.sav'
; Read background & calculate f-k corona contribution
  data = read_ascii(file,template=template)
  xx0=data.field1
  yy0=data.field2
  wl0=data.field3
  rr0=sqrt(xx0^2+yy0^2)
  cos2t0=2.0*(xx0/rr0)^2-1.0
  bk0=rr0^(-3.07)+82.0*rr0^(-6.95)
  bf0_temp=rr0^(0.22*cos2t0-2.47)
  index0=where(abs(rr0-2.3) lt 0.01)
  bk0_avg=average(bk0[index0])
  bf0_avg=average(bf0_temp[index0])
  cc=bk0_avg/bf0_avg
  bf0=bf0_temp*cc
  ratio_fk0=bf0/bk0
  fk_image0=fltarr(300,300)
  nn=long(0)
  for i=0, 299 do begin
     for j=0, 299 do begin
        fk_image0[i,j]=ratio_fk0[nn]
        nn=nn+1
     endfor
  endfor
  return,fk_image0
end

function extract_time,file  ; Extract time of image
  tmp1=strsplit(file,'/',/extract,count=count)
  tmp = strsplit(tmp1[count-1],'_',/extract)
  time = strsplit(tmp[4],'t',/extract)
  timing=FIX(time)/10000        ; Hrs
  tmp1=FIX(time) MOD 10000
  mins=tmp1/100                 ;Mins
;  print,STRING(mins,format='%03d')
  secs=tmp1 MOD 100
  plot_time = $
     String(timing,Format='%02d')+'Hr'+STRING(mins,FORMAT='%02d')+'Mins'
  return,plot_time
end

pro make_plot,dimage,inst,time,outfolder
  prange=[min(dimage>0),max(dimage)]
  print,'Prange WL= ',prange
  prange=[0.96,1.25]
  wdef,3,800,800
  case inst[1] of
     'C2': range=FLOAT(inst[4])
     'C3': range=FLOAT(inst[4])
     'COR1A': range=FLOAT(inst[4])
     'COR2A': range=FLOAT(inst[4])
     'COR1B': range=FLOAT(inst[4])
     'COR2B': range=FLOAT(inst[4])
  endcase
  dxy = inst[3]
  MAP=MAKE_MAP(DIMAGE,XC=0,YC=0,DX=dxy/300.,DY=dxy/300.)
  PLOT_MAP,MAP,TITLE='AWSoM CME SIM '+inst[1]+' T= '+trim(time), $
           CHARSIZE=2.0,XTITLE='SOLAR RADII',YTITLE='SOLAR RADII',$
           dmin=prange[0],dmax=prange[1], xthick=2,charthick=2,$
           ythick=2,xrange=[-range,range],yrange=[-range,range]
  plot_map_colorbar,prange,charsize=2.5,charthick=2
  tvcircle,inst[2],0,0,color=0,/data,/fill
  IMAGER=TVRD(TRUE=1)
  WRITE_PNG,outfolder+inst[0]+TRIM(time)+'.png',IMAGER,R,G,B
end

pro sim_movie,instrument,outfolder
  framerate=5
  xsize=800
  ysize=800
  video_file=outfolder+instrument+'movie.mp4'
  video_img = IDLffVideoWrite(video_file)
  stream_img = video_img.AddVideoStream(xsize,ysize,framerate)
  files_img = File_Search(outfolder+instrument+'*.png',Count=filecnt)
  for k=0,filecnt-1 do begin
     image_img = read_png(files_img[k])
     if (n_elements(image_img) eq 0) then break
     voidimg = video_img.Put(stream_img,image_img)
  endfor
  video_img.cleanup
  video_img=0
  image_img=0
end

pro final,file,outfolder,nlos
  if (nlos ne 0) then begin
     file_info = set_inst(file[0])
     print,'No. of LOS files found for '+file_info[1]+'= ',nlos
; Background
     wl_image0 = read_swmf_wl(file[0])
     fk_image0 = fk_corona(file[0])
; Next set of files
     for p=1,nlos-1 do begin
        wl_image=read_swmf_wl(file[p])
        dimage=wl_image/wl_image0
        dimage_corwl= (dimage+fk_image0)/(1.0+fk_image0)
        bad=WHERE(FINITE(dimage_corwl, /NAN))
        dimage_corwl[bad]=-999
        time_tmp = extract_time(file[p])
        make_plot,dimage_corwl,file_info,time_tmp,outfolder
     endfor
     sim_movie,file_info[0],outfolder
  endif else begin
     print,'No files found'
  endelse
end

pro make_WL_compare,CommonDir=CommonDir
  if (not keyword_set(CommonDir)) then CommonDir = 'Results/'
  
  ;CommonDir = '/Users/nishthas/Desktop/Results/run_test/'
  print,CommonDir
  nfolders = file_search(CommonDir+'run*/run*/',count=nfiles)
  if (nfiles eq 0) then begin
     nfolders = file_search(CommonDir+'run*/',count=nfilestmp)
     nfiles = nfilestmp
  endif
  if (nfiles eq 0) then begin
     nfolders = file_search(CommonDir+'SC/',count=nfilestmp)
     nfiles = nfilestmp
     nfolders = CommonDir
  endif
;;;;;
  for m=0,nfiles-1 do begin
     simdir = nfolders[m]+'/SC/'
     if (file_search(simdir+'/output/') ne simdir+'/output/') then $
        file_mkdir,simdir+'/output/'
     outdir = simdir+'/output/'
     print,'Simulation Folder',simdir
; C2
     loadct,3
     filename= file_search(simdir+'los_soho_c2*.dat',count=nlos)
     final,filename,outdir,nlos
; C3
     loadct,1
     filename= file_search(simdir+'los_soho_c3*.dat',count=nlos)
     final,filename,outdir,nlos
; NEXT COR1A
     secchi_colors,'cor1',/load
     filename = file_search(simdir+'los_sta_cor1*.dat',count=nlos)
     final,filename,outdir,nlos
; NEXT COR2A
     secchi_colors,'cor2',/load
     filename = file_search(simdir+'los_sta_cor2*.dat',count=nlos)
     final,filename,outdir,nlos
; NEXT COR1B
     secchi_colors,'cor1',/load
     filename = file_search(simdir+'los_stb_cor1*.dat',count=nlos)
     final,filename,outdir,nlos
; NEXT COR2B
     secchi_colors,'cor2',/load
     filename = file_search(simdir+'los_stb_cor2*.dat',count=nlos)
     final,filename,outdir,nlos
  endfor
end

