;set up the start and end time for obtaining the in-situ observations
pro compare_insitu, dir_sim=dir_sim, dir_plot=dir_plot,     $
                    extra_plt_info=extra_plt_info,          $
                    UseTimePlotName=UseTimePlotName,        $
                    CharSizeLocal=CharSizeLocal

  if (not keyword_set(dir_sim))   then dir_sim  = './simdata/'
  if (not keyword_set(dir_plot))  then dir_plot = './output/'
  if (not keyword_set(extra_plt_info)) then begin
     extra_plt_info = ''
  endif else begin
     if (strmid(extra_plt_info,0,1) ne '_') then $
        extra_plt_info = '_' + extra_plt_info
  endelse
  if (not keyword_set(UseTimePlotName)) then UseTimePlotName = 0
  if (not keyword_set(CharSizeLocal))   then CharSizeLocal = 2.5

  files_sim = file_search(dir_sim+'/*sat', count = nSimFile)

  if (nSimFile eq 0) then begin
     print, ' no simulation data'
     return
  endif

  for i = 0, n_elements(files_sim) -1 do begin
     file_sim    = files_sim(i)
     file_lowcase = strlowcase(file_sim)

     type = 'none'

     adapt_realization = ''

     print, 'file_lowcase=', file_lowcase

     if (strpos(file_lowcase, 'adapt')   ge 0) then begin
        index_adapt=strpos(file_lowcase, 'adapt')
        if (strmid(file_lowcase, index_adapt+5,1) eq '_') then begin
           adapt_realization = '_adapt_'+strmid(file_lowcase, index_adapt+6, 2)
        endif else begin
           adapt_realization = '_adapt_'+strmid(file_lowcase, index_adapt+5, 2)
        endelse
     endif

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

     read_swmf_sat, file_sim, time_swmf, n_swmf, ux_swmf, uy_swmf,       $
                    uz_swmf, bx_swmf, by_swmf, bz_swmf, ti_swmf, te_swmf, $
                    ut_swmf, ur_swmf, B_swmf
     
     start_time = time_swmf(0)
     end_time   = time_swmf(n_elements(time_swmf)-1)

     if (UseTimePlotName) then begin
        fileplot = dir_plot + start_time + TypePlot +           $
                   adapt_realization + extra_plt_info + '.eps'
     endif else begin
        CR_number = fix(tim2carr(time_swmf(n_elements(time_swmf)/2),/dc))
        fileplot = dir_plot + '/CR'+string(CR_number,format='(i4)') +TypePlot $
                   + adapt_realization + extra_plt_info + '.eps'
     endelse

     case type of
        'OMNI': begin
           obs_data=spdfgetdata('OMNI_COHO1HR_MERGED_MAG_PLASMA', $
                                ['ABS_B' , 'V', 'N','T'],         $
                                [start_time, end_time])

           if(not isa(obs_data)) then begin
              print, ' Could not obtain OMNI data for the simulation file: ',$
                 file_sim
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
              print, ' Could not obtain Stereo A data for the simulation file: ', file_sim
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
              print, ' Could not obtain Stereo B data for the simulation file: ', file_sim
              goto, out
           endif

           u_obs    = obs_data.plasmaSpeed.dat
           n_obs    = obs_data.plasmaDensity.dat
           tem_obs  = obs_data.plasmaTemp.dat
           mag_obs  = obs_data.B.dat
           time_obs = obs_data.epoch.dat
        end

        'none': begin
           print, ' case: unknown type, skip file: ', file_sim
           goto, out
        end
     endcase
     
     ;Set DoPlotTe=0 to exclude Te plot
     ;DoPlotTe = 0   
     
     plot_insitu, time_obs, u_obs,  n_obs,  tem_obs, mag_obs,             $
                  time_swmf, ur_swmf, n_swmf,  ti_swmf,  te_swmf, B_swmf, $
                  start_time, end_time, fileplot=fileplot, type=type,     $
                  charsize=CharSizeLocal, DoPlotTe=DoPlotTe

     print, ' finished plotting for type: ', type
     out:

     print,'Calculating integrated curve distance between obs. & SWMF output at 1 AU'
     ;Normalization for OMNI data
     t_norm=10                  ;10 days
     u_norm = max(u_obs) - min(u_obs(where(u_obs gt 0)))
     n_norm = max(n_obs) - min(n_obs(where(n_obs gt 0)))
     tem_norm = max(tem_obs) - min(tem_obs(where( tem_obs gt 0)))
     mag_norm = max(mag_obs) - min(mag_obs(where( mag_obs gt 0)))

     print,'Normalizations:'
     help,t_norm,u_norm,n_norm,tem_norm,mag_norm
     
     ;time in units of t_norm days
     t_obs=time_obs/(24.*60.*60.*1.e3)/t_norm
     ; Converting swmf time (YYYY-MO-DDTHH:MM:SS) to epoch time in sec
     TIMESTAMPTOVALUES,time_swmf+'Z',year=yy,month=mo,day=dy,hour=hh,min=mm,sec=ss
     cdf_epoch,t_swmf,yy,mo,dy,hh,mm,ss,/compute_epoch
     t_swmf=t_swmf/(24.*60.*60.*1.e3)/t_norm

     index_u = where(u_obs gt 0)
     index_n = where(n_obs gt 0)
     index_T = where(tem_obs gt 0)
     index_B = where(mag_obs gt 0)

     dist_int_u=curve_int_distance(t_obs(index_u),u_obs(index_u),t_swmf,ur_swmf)
     dist_int_t=curve_int_distance(t_obs(index_T),tem_obs(index_T)/tem_norm,$
                                   t_swmf,ti_swmf/tem_norm)
     dist_int_n=curve_int_distance(t_obs(index_n),n_obs(index_n)/n_norm,$
                                   t_swmf,n_swmf/n_norm)
     dist_int_b=curve_int_distance(t_obs(index_B),mag_obs(index_B)/mag_norm,$
                                   t_swmf,B_swmf/mag_norm)
     print,'Int Curve distance (ur,n,ti,b) for file',file_lowcase,' is: '$
           ,dist_int_u,dist_int_n,dist_int_t,dist_int_b

  endfor

end
