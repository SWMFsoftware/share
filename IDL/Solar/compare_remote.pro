pro compare_remote, dir_sim=dir_sim, dir_plot=dir_plot,     $
                    extra_plt_info=extra_plt_info,          $
                    UseTimePlotName=UseTimePlotName,        $
                    CharSizeLocal=CharSizeLocal,            $
                    TypePlotFile=TypePlotFile,              $
                    xs_map=xs_map, ys_map=ys_map

  if (not keyword_set(dir_sim))  then dir_sim  = './simdata/'
  if (not keyword_set(dir_plot)) then dir_plot = './output/'
  if (not keyword_set(extra_plt_info)) then begin
     extra_plt_info = ''
  endif else begin
     if (strmid(extra_plt_info,0,1) ne '_') then $
        extra_plt_info = '_' + extra_plt_info
  endelse

  if (not keyword_set(UseTimePlotName)) then UseTimePlotName = 0
  if (not keyword_set(TypePlotFile))    then TypePlotFile  = 'eps'
  if (not keyword_set(CharSizeLocal))   then CharSizeLocal = 1.5

  if (not keyword_set(xs_map)) then xs_map = 512
  if (not keyword_set(ys_map)) then ys_map = 512
  
  get_lun, unitlog
  openw, unitlog, 'log_compare_remote.log'

  files_sim = file_search(dir_sim+'/los*dat', count = nSimFile)

  if (nSimFile eq 0) then begin
     printf, unitlog, ' no simulation data'
     close, unitlog
     free_lun, unitlog
     return
  endif

  for i = 0, n_elements(files_sim) -1 do begin
     file_sim = files_sim(i)

     read_swmf_remote_tec, file_sim, TimeEvent=TimeEvent, ObsXyz=ObsXyz,    $
                           varnames=varnames, nvars=nvars, data=dataSim,    $
                           nx = nx, ny = ny, rMaxSim = rMaxSim, DoWl =DoWl, $
                           DoAIA=DoAIA, DoXRT=DoXRT, DoEUV=DoEUV

     determine_NameSat, TimeEvent, ObsXyz, NameSat, file_sim=file_sim

     printf, unitlog, '------------------------------------------------------'
     printf, unitlog, 'For ', file_sim
     printf, unitlog, 'NameSat = ', NameSat
     printf, unitlog, 'DoWl    = ', DoWl,  format='(a,i2)'
     printf, unitlog, 'DoAIA   = ', DoAIA, format='(a,i2)'
     printf, unitlog, 'DoXRT   = ', DoXRT, format='(a,i2)'
     printf, unitlog, 'DoEUV   = ', DoEUV, format='(a,i2)'
     printf, unitlog, ''

     compare_remote_one, TimeEvent=TimeEvent, varnames=varnames,              $
                         nvars=nvars, dataSim=dataSim, nx=nx, ny=ny,          $
                         CharSizeLocal=CharSizeLocal, dir_plot=dir_plot,      $
                         extra_plt_info=extra_plt_info,                       $
                         UseTimePlotName=UseTimePlotName, unitlog=unitlog,    $
                         NameSat=NameSat, rMaxSim=rMaxSim,                    $
                         DoWl =DoWl, DoAIA=DoAIA, DoXRT=DoXRT, DoEUV=DoEUV,   $
                         TypePlotFile=TypePlotFile, xs_map=xs_map,            $
                         ys_map=ys_map
  endfor

  printf, unitlog, '------------------------------------------------------'

  close, unitlog
  free_lun, unitlog
end
