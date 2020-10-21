;set up the start and end time for obtaining the in-situ observations
pro compare_insitu, dir_sim=dir_sim, dir_plot=dir_plot,     $
                    extra_plt_info=extra_plt_info,          $
                    UseTimePlotName=UseTimePlotName,        $
                    CharSizeLocal=CharSizeLocal,            $
                    DoPlotTe=DoPlotTe, Model=Model,         $
                    dir_obs=dir_obs

  if (not keyword_set(dir_sim)) then begin
     if (file_test('./simdata', /directory)) then begin
        dir_sim  = './simdata/'
        print, ' Uses the default dir_sim = ./simdata'
     endif else begin
        print, ' Please specify the directory containing simulation results'
        return
     endelse
  endif

  if (not keyword_set(dir_plot)) then begin
     if (file_test('./output', /directory) eq 0) then file_mkdir, './output'
     dir_plot = './output'
     print, ' Saves into the default dir_plot = ./output'
  endif

  if (not keyword_set(dir_obs)) then begin
     if (file_test('./obsdata', /directory) eq 0) then file_mkdir, './output'
     dir_obs = './obsdata'
     print, ' Saves into the default dir_obs = ./obsdata'
  endif

  if (not keyword_set(extra_plt_info)) then begin
     extra_plt_info = ''
  endif else begin
     if (strmid(extra_plt_info,0,1) ne '_') then $
        extra_plt_info = '_' + extra_plt_info
  endelse
  if (not keyword_set(UseTimePlotName)) then UseTimePlotName = 0
  if (not keyword_set(CharSizeLocal))   then CharSizeLocal = 2.5
  if (not keyword_set(DoPlotTe))        then DoPlotTe = 0

  if (not keyword_set(Model)) then Model = 'AWSoM'

  files_sim = file_search(dir_sim+'/*sat', count = nSimFile)

  if (nSimFile eq 0) then begin
     print, ' no simulation data'
     return
  endif

  for i = 0, n_elements(files_sim) -1 do begin
     file_sim     = files_sim(i)

     read_swmf_sat, file_sim, time_swmf, n_swmf, ux_swmf, uy_swmf,        $
                    uz_swmf, bx_swmf, by_swmf, bz_swmf, ti_swmf, te_swmf, $
                    ut_swmf, ur_swmf, B_swmf, DoContainData=DoContainData,$
                    TypeData=TypeData, TypePlot=TypePlot,                 $
                    start_time=start_time, end_time=end_time
     
     if DoContainData ne 1 then begin
        print, " Error: filename=", file_sim, " does not contain any data"
        continue
     endif

     get_insitu_data, start_time, end_time, TypeData, u_obs, n_obs, tem_obs,  $
                      mag_obs, time_obs, DoContainData=DoContainData

     if DoContainData ne 1 then begin
        print, " Error: no observational data are found."
        continue
     endif

     ;; save the observation for QU.
     EventTime = strmid(time_swmf(n_elements(time_swmf)/2),0,19)
     EventTime=repstr(EventTime,'-','_')
     EventTime=repstr(EventTime,':','_')

     ObsFileName = dir_obs + '/' + strmid(TypePlot,1) + '_' +EventTime  + '.out'
     w = fltarr(n_elements(u_obs),5)
     w = [[time_obs], [n_obs], [u_obs], [tem_obs], [mag_obs]]
     x = indgen(n_elements(u_obs))
     varname = ['count', 'time', 'density', 'velocity', 'temperature', 'magnetic_field']
     save_pict,ObsFileName,TypeData+' Observational data',varname,w,x
     print, ' saving the observation file to ObsFileName: ', ObsFileName

     if (UseTimePlotName) then begin
        fileplot = dir_plot + '/' + start_time + TypePlot                   $
                   + extra_plt_info + '.eps'
     endif else begin
        CR_number = fix(tim2carr(time_swmf(n_elements(time_swmf)/2),/dc))
        fileplot = dir_plot + '/CR'+string(CR_number,format='(i4)') +TypePlot  $
                   + extra_plt_info + '.eps'
     endelse

     set_plot,'PS'
     device,/encapsulated
     device,filename=fileplot,/color,bits_per_pixel=8
     device,xsize=20,ysize=20

     plot_insitu, time_obs, u_obs,  n_obs,  tem_obs, mag_obs,             $
                  time_swmf, ur_swmf, n_swmf,  ti_swmf,  te_swmf, B_swmf, $
                  start_time, end_time, typeData=typeData,                $
                  charsize=CharSizeLocal, DoPlotTe = DoPlotTe,            $
                  legendNames=Model

     device,/close_file

     print, ' finished plotting for TypeData: ', TypeData
  endfor
end

