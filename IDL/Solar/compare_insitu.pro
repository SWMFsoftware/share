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

     compare_insitu_one, file_sim=file_sim, extra_plt_info=extra_plt_info, $
                         UseTimePlotName=UseTimePlotName,                  $
                         CharSizeLocal=CharSizeLocal, DoPlotTe=DoPlotTe,   $
                         Model=Model, dir_obs=dir_obs, dir_plot=dir_plot,  $
                         DoSaveObs=1

  endfor
end

