pro compare_AIA, TimeEvent=TimeEvent, varnames=varnames, nvars=nvars,        $
                 dataSim=dataSim, nx=nx, ny=ny,                              $
                 dir_obs=dir_obs, dir_plot=dir_plot,                         $
                 extra_plt_info=extra_plt_info, TypePlotFile=TypePlotFile,   $
                 UseTimePlotName=UseTimePlotName,                            $
                 CharSizeLocal = CharSizeLocal, unitlog=unitlog,             $
                 NameSat = NameSat, xs_map=xs_map,ys_map=ys_map

  NameSub = 'compare_AIA'

  if (not keyword_set(UseTimePlotName)) then UseTimePlotName = 0
  if (not keyword_set(TypePlotFile))    then TypePlotFile    = 'png'

  if (not keyword_set(dir_obs))  then dir_obs  ='./obsdata/'
  if (not keyword_set(dir_plot)) then dir_plot ='./output/'
  
  if (not keyword_set(xs_map)) then xs_map = 512
  if (not keyword_set(ys_map)) then ys_map = 512

  ;; params for plotting
  imagefov = [-1200,1200]
  plotmin  = 1.e-2

  InstPlt = 'aia'

  set_FilePlotName, TimeEvent, UseTimePlotName = UseTimePlotName,    $
                    TimeStrFile = TimeStrFile, fileplot = fileplot,  $
                    dir_plot = dir_plot, InstPlt = InstPlt

  ;;-------------------------------------------------------------------------
  ;; Download the data

  download_images, TimeEvent = TimeEvent, CaseInst = 'aia',       $
                   InstPlt = InstPlt, filename_I = filename_I,    $
                   dir_obs = dir_obs, TimeStrFile = TimeStrFile,  $
                   unitlog = unitlog

  ;;-------------------------------------------------------------------------
  ;; prepare the image for plotting

  process_aia, filename_I(0), aia_map  = aia94_map,  $
               xy_map = xs_map, ys_map = ys_map

  process_aia, filename_I(1), aia_map  = aia131_map, $
               xy_map = xs_map, ys_map = ys_map

  process_aia, filename_I(2), aia_map  = aia171_map, $
               xy_map = xs_map, ys_map = ys_map

  process_aia, filename_I(3), aia_map  = aia193_map, $
               xy_map = xs_map, ys_map = ys_map

  process_aia, filename_I(4), aia_map  = aia211_map, $
               xy_map = xs_map, ys_map = ys_map

  process_aia, filename_I(5), aia_map  = aia304_map, $
               xy_map = xs_map, ys_map = ys_map

  process_aia, filename_I(6), aia_map  = aia335_map, $
               xy_map = xs_map, ys_map = ys_map

  ;;-------------------------------------------------------------------------
  ;; processing the swmf data
  rSizeImage=1.98
  dxy = rSizeImage*2000./nx

  make_map_swmf_data,dataSim=dataSim, varnames=varnames, namevar='aia:94',  $
                     nx=nx, ny = ny, dxy=dxy, TimeEvent=TimeEvent,          $
                     map_out =aia94_image
  make_map_swmf_data,dataSim=dataSim, varnames=varnames, namevar='aia:131', $
                     nx=nx, ny = ny, dxy=dxy, TimeEvent=TimeEvent,          $
                     map_out =aia131_image
  make_map_swmf_data,dataSim=dataSim, varnames=varnames, namevar='aia:171', $
                     nx=nx, ny = ny, dxy=dxy, TimeEvent=TimeEvent,          $
                     map_out =aia171_image
  make_map_swmf_data,dataSim=dataSim, varnames=varnames, namevar='aia:193', $
                     nx=nx, ny = ny, dxy=dxy, TimeEvent=TimeEvent,          $
                     map_out =aia193_image
  make_map_swmf_data,dataSim=dataSim, varnames=varnames, namevar='aia:211', $
                     nx=nx, ny = ny, dxy=dxy, TimeEvent=TimeEvent,          $
                     map_out =aia211_image
  make_map_swmf_data,dataSim=dataSim, varnames=varnames, namevar='aia:304', $
                     nx=nx, ny = ny, dxy=dxy, TimeEvent=TimeEvent,          $
                     map_out =aia304_image
  make_map_swmf_data,dataSim=dataSim, varnames=varnames, namevar='aia:335', $
                     nx=nx, ny = ny, dxy=dxy, TimeEvent=TimeEvent,          $
                     map_out =aia335_image

  ;;-------------------------------------------------------------------------
  ;; plotting

  case TypePlotFile of
     'png': begin
        set_plot,'X'

        nx_plot=3
        ny_plot=2
        !p.multi=[0,3,2]
        cal_pos,nx_plot,ny_plot,y_x_ratio=y_x,x1=x1,y1=y1,x2=x2,y2=y2,$
                tm=0.1,dy=0,dx=0.5,lm=0.5
        
        wdef,2,800,600

        ny_offset = 0
     end
     'eps': begin
        set_plot,'PS'
        device,/encapsulated
        device,filename=fileplot + extra_plt_info+'.eps', $
               /color,bits_per_pixel=8
        device,xsize=20,ysize=26

        nx_plot=3
        ny_plot=4
        !p.multi=[0,3,4]
        cal_pos,nx_plot,ny_plot,y_x_ratio=y_x,x1=x1,y1=y1,x2=x2,y2=y2,$
                tm=-0.5,dy=-0.5,dx=0.5,lm=0.8,bm=0.3

        ny_offset = 2
     end
  endcase

  ;;-----------------------------  94 -----------------------------------------
  aia_lct,wavelnth=94,/load
  pos=[x1[0],y1[1+ny_offset],x2[0],y2[1+ny_offset]]
  plotmax=max([max(aia94_map.data),max(aia94_image.data)])
  plot_map_local, aia94_image, position=pos, charsize=CharSizeLocal,   $
                  dmin=plotmin, dmax=plotmax,                          $
                  xrange=imagefov,yrange=imagefov,                     $
                  xtitle='', title='Model AIA 94'

  pos=[x1[0],y1[0+ny_offset],x2[0],y2[0+ny_offset]]
  plot_map_local, aia94_map, position=pos, charsize=CharSizeLocal,     $
                  dmin=plotmin, dmax=plotmax,                          $
                  xrange=imagefov,yrange=imagefov,                     $
                  xtitle='', title='SDO AIA 94'

  ;;----------------------------- 171 -----------------------------------------
  aia_lct,wavelnth=171,/load
  pos=[x1[1],y1[1+ny_offset],x2[1],y2[1+ny_offset]]
  plotmax=max([max(aia171_map.data),max(aia171_image.data)])
  plot_map_local, aia171_image, position=pos, charsize=CharSizeLocal,  $
                  dmin=plotmin, dmax=plotmax,                          $
                  xrange=imagefov,yrange=imagefov,                     $
                  xtitle='', ytitle='', title='Model AIA 171'

  pos=[x1[1],y1[0+ny_offset],x2[1],y2[0+ny_offset]]
  plot_map_local, aia171_map,position=pos, charsize=CharSizeLocal,     $
                  dmin=plotmin, dmax=plotmax,                          $
                  xrange=imagefov,yrange=imagefov,                     $
                  xtitle='', ytitle='',title='SDO AIA 171'

  ;;----------------------------- 193 -----------------------------------------
  aia_lct,wavelnth=193,/load
  pos=[x1[2],y1[1+ny_offset],x2[2],y2[1+ny_offset]]
  plotmax=max([max(aia193_map.data),max(aia193_image.data)])
  plot_map_local, aia193_image, position=pos, charsize=CharSizeLocal,   $
                  dmin=plotmin, dmax=plotmax,                           $
                  xrange=imagefov, yrange=imagefov,                     $
                  xtitle='',ytitle='',title='Model AIA 193'

  pos=[x1[2],y1[0+ny_offset],x2[2],y2[0+ny_offset]]
  plot_map_local, aia193_map, position=pos, charsize=CharSizeLocal,     $
                  dmin=plotmin, dmax=plotmax,                           $
                  xrange=imagefov,yrange=imagefov,                      $
                  xtitle='',ytitle='',title='SDO AIA 193'

  ;;--------------------------- save png part 1 ------------------------------
  case TypePlotFile of
     'png': begin
        imager=tvrd(true=1)
        write_png,fileplot+'_part1'+extra_plt_info+'.png',imager,r,g,b
        ; wdef,3,1200,900
        wdef,3,800,600
     end
     'eps':
  endcase

  ;;----------------------------- 131 -----------------------------------------
  aia_lct,wavelnth=131,/load
  pos=[x1[0],y1[1],x2[0],y2[1]]
  plotmax=max([max(aia131_map.data),max(aia131_image.data)])
  plot_map_local, aia131_image, position=pos, charsize=CharSizeLocal,  $
                  dmin=plotmin, dmax=plotmax,                          $
                  xrange=imagefov,yrange=imagefov,                     $
                  xtitle='',title='Model AIA 131'
  pos=[x1[0],y1[0],x2[0],y2[0]]
  plot_map_local, aia131_map, position=pos, charsize=CharSizeLocal,    $
                  dmin=plotmin, dmax=plotmax,                          $
                  xrange=imagefov,yrange=imagefov,                     $
                  title='SDO AIA 131'

  ;;----------------------------- 211 -----------------------------------------
  aia_lct,wavelnth=211,/load
  pos=[x1[1],y1[1],x2[1],y2[1]]
  plotmax=max([max(aia211_map.data),max(aia211_image.data)])
  plot_map_local, aia211_image, position=pos, charsize=CharSizeLocal,  $
                  dmin=plotmin, dmax=plotmax,                          $
                  xrange=imagefov,yrange=imagefov,                     $
                  xtitle='',ytitle='',title='Model AIA 211'
  pos=[x1[1],y1[0],x2[1],y2[0]]
  plot_map_local, aia211_map, position=pos, charsize=CharSizeLocal,    $
                  dmin=plotmin, dmax=plotmax,                          $
                  xrange=imagefov,yrange=imagefov,                     $
                  ytitle='',title='SDO AIA 211'

  ;;----------------------------- 331 -----------------------------------------
  aia_lct,wavelnth=335,/load
  pos=[x1[2],y1[1],x2[2],y2[1]]
  plotmax=max([max(aia335_map.data),max(aia335_image.data)])
  plot_map_local, aia335_image, position=pos, charsize=CharSizeLocal,  $
                  dmin=plotmin, dmax=plotmax,                          $
                  xrange=imagefov,yrange=imagefov,                     $
                  xtitle='',ytitle='',title='Model AIA 335'
  pos=[x1[2],y1[0],x2[2],y2[0]]
  plot_map_local, aia335_map, position=pos, charsize=CharSizeLocal,    $
                  dmin=plotmin, dmax=plotmax,                          $
                  xrange=imagefov,yrange=imagefov,                     $
                  ytitle='',title='SDO AIA 335'

  ;;----------------------------- save ---------------------------------------
  case TypePlotFile of
     'png': begin
        imager=tvrd(true=1)
        write_png,fileplot+'_part2'+extra_plt_info+'.png',imager,r,g,b

        ;; reset
        wdelete,2,3
     end
     'eps': begin
        device,/close
        
        ;; reset window
        set_plot,'X'
     end
  endcase

  !p.multi=0
  loadct,0

  out: printf, unitlog, NameSub, ' is finished.'
  printf, unitlog, ''
end
