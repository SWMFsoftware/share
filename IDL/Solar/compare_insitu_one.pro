;set up the start and end time for obtaining the in-situ observations
pro compare_insitu_one, file_sim=file_sim,                      $
                        extra_plt_info=extra_plt_info,          $
                        UseTimePlotName=UseTimePlotName,        $
                        CharSizeLocal=CharSizeLocal,            $
                        DoPlotTe=DoPlotTe, Model=Model,         $
                        dir_obs=dir_obs, dir_plot=dir_plot,     $
                        DoSaveObs=DoSaveObs,DoLogT=DoLogT

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

  if (not isa(DoSaveObs)) then DoSaveObs = 1

  print, "compare_insitu_one: DoSaveObs=", DoSaveObs

  read_swmf_sat, file_sim, time_swmf, n_swmf, ux_swmf, uy_swmf,        $
                 uz_swmf, bx_swmf, by_swmf, bz_swmf, ti_swmf, te_swmf, $
                 ut_swmf, ur_swmf, B_swmf, DoContainData=DoContainData,$
                 TypeData=TypeData, TypePlot=TypePlot,                 $
                 start_time=start_time, end_time=end_time
  
  if DoContainData ne 1 then begin
     print, "compare_insitu_one error: filename=", $
            file_sim, " does not contain any data"
     return
  endif

  get_insitu_data, start_time, end_time, TypeData, u_obs, n_obs, tem_obs,  $
                   mag_obs, time_obs, DoContainData=DoContainData

  if DoContainData ne 1 then begin
     print, "compare_insitu_one: error: no observational data are found."
     return
  endif

  ;; save the observation for QU.
  if (DoSaveObs eq 1) then begin
     EventTime = strmid(time_swmf(n_elements(time_swmf)/2),0,19)
     EventTime=repstr(EventTime,'-','_')
     EventTime=repstr(EventTime,':','_')

     print, "compare_insitu_one: dir_obs =", dir_obs

     utc_obs = anytim2utc(cdf2utc(time_obs),/external)
     yy = utc_obs.YEAR
     mo = utc_obs.MONTH
     dy = utc_obs.DAY
     hh = utc_obs.HOUR
     mm = utc_obs.MINUTE
     ss = utc_obs.SECOND

     ObsFileName = dir_obs + '/' + strmid(TypePlot,1) + '_' +EventTime + '.out'
     w = fltarr(n_elements(u_obs),5)
     w = [[yy], [mo], [dy], [hh], [mm], [ss], [n_obs], [u_obs], [tem_obs], [mag_obs]]
     x = indgen(n_elements(u_obs))
     varname = ['count', 'year', 'mo', 'dy', 'hr', 'mn', 'sc', 'Rho', $
                'V_tot', 'Temperature', 'B_tot']
     save_pict,ObsFileName,TypeData+' Observational data',varname,w,x
     print, 'compare_insitu_one: saving the observation file to ObsFileName: ', ObsFileName
  endif

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
               legendNames=Model, DoLogT=DoLogT,                       $
               file_dist=fileplot.replace('.eps','.txt')

  device,/close_file
end

