;set up the start and end time for obtaining the in-situ observations
pro compare_diff_models, dir_sim=dir_sim, dir_plot=dir_plot,                  $
                         file1=file1, file2=file2, file3=file3, file4=file4,  $
                         legendNames=legendNames, CharSizeLocal=CharSizeLocal,$
                         DoPlotTe=DoPlotTe

  if (not keyword_set(dir_sim))   then dir_sim  = '../simdata/'
  if (not keyword_set(dir_plot))  then dir_plot = '../output/'
  if (not keyword_set(CharSizeLocal))   then CharSizeLocal = 2.5

  file1 = dir_sim+file1
  file2 = dir_sim+file2
  if (keyword_set(file3)) then file3 = dir_sim+file3
  if (keyword_set(file4)) then file4 = dir_sim+file4

  file_lowcase = strlowcase(file1)

  type = 'none'

  if (strpos(file_lowcase, 'earth') ge 0) then begin
     type = 'OMNI'
     TypePlot  = '_omni'
  endif else if (strpos(file_lowcase, 'sta')     ge 0 or $
                 strpos(file_lowcase, 'stereoa') ge 0) then begin
     type = 'Stereo A'
     TypePlot  = '_sta'
  endif else if (strpos(file_lowcase, 'stb')     ge 0 or $
                 strpos(file_lowcase, 'stereob') ge 0) then begin
     type = 'Stereo B'
     TypePlot  = '_stb'
  endif

  read_swmf_sat, file1, time_swmf1, n_swmf1, ux_swmf1, uy_swmf1,       $
                 uz_swmf1, bx_swmf1, by_swmf1, bz_swmf1, ti_swmf1, te_swmf1, $
                 ut_swmf1, ur_swmf1, B_swmf1

  read_swmf_sat, file2, time_swmf2, n_swmf2, ux_swmf2, uy_swmf2,       $
                 uz_swmf2, bx_swmf2, by_swmf2, bz_swmf2, ti_swmf2, te_swmf2, $
                 ut_swmf2, ur_swmf2, B_swmf2

  if (keyword_set(file3)) then $
     read_swmf_sat, file3, time_swmf3, n_swmf3, ux_swmf3, uy_swmf3,       $
                    uz_swmf3, bx_swmf3, by_swmf3, bz_swmf3, ti_swmf3,     $
                    te_swmf3, ut_swmf3, ur_swmf3, B_swmf3

  if (keyword_set(file4)) then $
     read_swmf_sat, file4, time_swmf4, n_swmf4, ux_swmf4, uy_swmf4,       $
                    uz_swmf4, bx_swmf4, by_swmf4, bz_swmf4, ti_swmf4,     $
                    te_swmf4, ut_swmf4, ur_swmf4, B_swmf4
  
  start_time = time_swmf1(0)
  end_time   = time_swmf1(n_elements(time_swmf1)-1)

  CR_number = fix(tim2carr(time_swmf1(n_elements(time_swmf1)/2),/dc))
  fileplot = dir_plot + '/CR'+string(CR_number,format='(i4)') +TypePlot $
             + '_tmp' + '.eps'

  case type of
     'OMNI': begin
        obs_data=spdfgetdata('OMNI_COHO1HR_MERGED_MAG_PLASMA', $
                             ['ABS_B' , 'V', 'N','T'],         $
                             [start_time, end_time])

        if(not isa(obs_data)) then begin
           print, ' Could not obtain OMNI data for the simulation file: ',$
                  file1
           goto, out
        endif

        u_obs    = obs_data.v.dat
        n_obs    = obs_data.n.dat
        tem_obs  = obs_data.t.dat
        mag_obs  = obs_data.abs_b.dat
        time_obs = obs_data.epoch.dat
     end

     'Stereo A': begin
        obs_data=spdfgetdata('STA_COHO1HR_MERGED_MAG_PLASMA',       $
                             ['B' , 'plasmaSpeed', 'plasmaDensity', $
                              'plasmaTemp'],[start_time, end_time])

        if(not isa(obs_data)) then begin
           print, ' Could not obtain Stereo A data for the simulation file: ',$
                  file1
           goto, out
        endif

        u_obs    = obs_data.plasmaSpeed.dat
        n_obs    = obs_data.plasmaDensity.dat
        tem_obs  = obs_data.plasmaTemp.dat
        mag_obs  = obs_data.B.dat
        time_obs = obs_data.epoch.dat
     end

     'Stereo B': begin
        obs_data=spdfgetdata('STB_COHO1HR_MERGED_MAG_PLASMA',       $
                             ['B' , 'plasmaSpeed', 'plasmaDensity', $
                              'plasmaTemp'],[start_time, end_time])

        if(not isa(obs_data)) then begin
           print, ' Could not obtain Stereo B data for the simulation file: ',$
                  file1
           goto, out
        endif

        u_obs    = obs_data.plasmaSpeed.dat
        n_obs    = obs_data.plasmaDensity.dat
        tem_obs  = obs_data.plasmaTemp.dat
        mag_obs  = obs_data.B.dat
        time_obs = obs_data.epoch.dat
     end

     'none': begin
        print, ' case: unknown type, skip file: ', file1
        goto, out
     end
  endcase

  plot_insitu, time_obs, u_obs,  n_obs,  tem_obs, mag_obs,                   $
               time_swmf1, ur_swmf1, n_swmf1,  ti_swmf1,  te_swmf1, B_swmf1, $
               start_time, end_time, fileplot=fileplot, type=type,           $
               charsize=CharSizeLocal, legendNames = legendNames,            $
               time_simu2=time_swmf2, u_simu2=ur_swmf2, n_simu2=n_swmf2,     $
               ti_simu2=ti_swmf2, te_simu2=te_swmf2, b_simu2=b_swmf2,        $
               time_simu3=time_swmf3, u_simu3=ur_swmf3, n_simu3=n_swmf3,     $
               ti_simu3=ti_swmf3, te_simu3=te_swmf3, b_simu3=b_swmf3,        $
               time_simu4=time_swmf4, u_simu4=ur_swmf4, n_simu4=n_swmf4,     $
               ti_simu4=ti_swmf4, te_simu4=te_swmf4, b_simu4=b_swmf4,        $
               DoPlotTe=DoPlotTe

  print, ' finished plotting for type: ', type

  out: 
end
