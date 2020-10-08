;set up the start and end time for obtaining the in-situ observations
pro compare_insitu, dir_sim=dir_sim, dir_plot=dir_plot,     $
                    extra_plt_info=extra_plt_info,          $
                    UseTimePlotName=UseTimePlotName,        $
                    CharSizeLocal=CharSizeLocal,            $
                    DoPlotTe=DoPlotTe, Model=Model

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
     file_lowcase = strlowcase(file_sim)

     type = 'none'

     print, ' file_lowcase=', file_lowcase

     if (strpos(file_lowcase, 'earth') ge 0) then begin
        type      = 'OMNI'
        TypePlot  = '_omni'
     endif else if (strpos(file_lowcase, 'sta')     ge 0 or $
                    strpos(file_lowcase, 'stereoa') ge 0) then begin
        type      = 'Stereo A'
        TypePlot  = '_sta'
     endif else if (strpos(file_lowcase, 'stb')     ge 0 or $
                    strpos(file_lowcase, 'stereob') ge 0) then begin
        type      = 'Stereo B'
        TypePlot  = '_stb'
     endif

     read_swmf_sat, file_sim, time_swmf, n_swmf, ux_swmf, uy_swmf,        $
                    uz_swmf, bx_swmf, by_swmf, bz_swmf, ti_swmf, te_swmf, $
                    ut_swmf, ur_swmf, B_swmf, DoContainData=DoContainData
     
     if DoContainData ne 1 then begin
        print, " Error: filename=", file_sim, " does not contain any data"
        continue
     endif

     start_time = time_swmf(0)
     end_time   = time_swmf(n_elements(time_swmf)-1)

     get_insitu_data, start_time, end_time, type, u_obs, n_obs, tem_obs,  $
                      mag_obs, time_obs, DoContainData=DoContainData

     if DoContainData ne 1 then begin
        print, " Error: no observational data are found."
        continue
     endif

     if (UseTimePlotName) then begin
        fileplot = dir_plot + start_time + TypePlot                   $
                   + extra_plt_info + '.eps'
     endif else begin
        CR_number = fix(tim2carr(time_swmf(n_elements(time_swmf)/2),/dc))
        fileplot = dir_plot + '/CR'+string(CR_number,format='(i4)') +TypePlot  $
                   + extra_plt_info + '.eps'
     endelse

     plot_insitu, time_obs, u_obs,  n_obs,  tem_obs, mag_obs,             $
                  time_swmf, ur_swmf, n_swmf,  ti_swmf,  te_swmf, B_swmf, $
                  start_time, end_time, fileplot=fileplot, type=type,     $
                  charsize=CharSizeLocal, DoPlotTe = DoPlotTe,            $
                  legendNames=Model
     print, ' finished plotting for type: ', type
  endfor
end

