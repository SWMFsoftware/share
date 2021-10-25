pro compare_insitu_from_list, filename_list=filenmae_list, dir_plot=dir_plot,   $
                              DoPlotTe=DoPlotTe, CharSizeLocal=CharSizeLocal,   $
                              nyMaxLegend=nyMaxLegend

  if (not keyword_set(filenmae_list))  then filename_list  = './filename_list.txt'
  if (not keyword_set(dir_plot))       then dir_plot = './output/'
  if (not keyword_set(DoPlotTe))       then DoPlotTe = 0
  if (not keyword_set(CharSizeLocal))  then CharSizeLocal = 2.5
  if (not keyword_set(nyMaxLegend))    then nyMaxLegend=3

  umax_plot = 900
  nmax_plot = 35
  tmax_plot = 1e6
  bmax_plot = 25

  umax = 0
  nmax = 0
  tmax = 0
  bmax = 0

  ;; first iteration to get the CR number and the range of the plot...
  openr, lun, filename_list, /get_lun
  line = ''
  ;; everything before #START is header...
  while not eof(lun) do begin
     readf,lun,line
     if line ne '#START' then continue
     break
  end
  while not eof(lun) do begin
     readf,lun,line
     strings_local = strsplit(line,',', /EXTRACT)
     filename_sat = strings_local[0]

     read_swmf_sat, filename_sat, time_swmf, n_swmf, ux_swmf, uy_swmf,       $
                    uz_swmf, bx_swmf, by_swmf, bz_swmf, ti_swmf, te_swmf,    $
                    ut_swmf, ur_swmf, B_swmf, DoContainData=DoContainData,   $
                    TypeData=TypeData, TypePlot=TypePlot,                    $
                    start_time=start_time, end_time=end_time

     if DoContainData ne 1 then begin
        print, " Error: filename=", file_sim_adapt, " does not contain any data"
        continue
     endif

     CR_number = fix(tim2carr(time_swmf(n_elements(time_swmf)/2),/dc))

     if (not keyword_set(CR_number_plot)) then $
        CR_number_plot = CR_number

     if CR_number_plot ne CR_number then begin
        print, 'CR_number_plot /= CR_number, the comparison must be done with the same CR'
        print, 'correct the list to proceed. '
        print, 'CR_number_plot =', CR_number_plot
        print, 'CR_number      =', CR_number
        return
     endif

     umax = max([umax,max(ut_swmf)*1.3])
     nmax = max([nmax,max(n_swmf )*1.3])
     tmax = max([tmax,max(ti_swmf)*1.3])
     bmax = max([bmax,max(B_swmf )*1.3]*1e5)
     if DoPlotTe then tmax = max([tmax,max(te_swmf)*1.3])
  end

  get_insitu_data, start_time, end_time, TypeData, u_obs, n_obs, tem_obs,  $
                   mag_obs, time_obs, DoContainData=DoContainData

  if DoContainData ne 1 then begin
     print, " Error: no observational data are found."
     return
  endif
  
  ;; adjust again with the observation
  umax_plot = min([umax_plot,max([umax,max(u_obs)*1.3])])
  nmax_plot = min([nmax_plot,max([nmax,max(n_obs)*1.3])])
  tmax_plot = min([tmax_plot,max([tmax,max(tem_obs)*1.3])])
  bmax_plot = min([bmax_plot,max([bmax,max(mag_obs)*1.3])])

  ;; set the filename for the output figure...
  fileplot = dir_plot + '/CR'+string(CR_number_plot,format='(i4)') +TypePlot $
             + '_compare_list' + '.eps'

  set_plot,'PS'
  device,/encapsulated
  device,filename=fileplot,/color,bits_per_pixel=8
  device,xsize=20,ysize=20

  IsOverPlot       = 0
  nLegendPlotTotal = 0
  
  openr, lun, filename_list, /get_lun
  while not eof(lun) do begin
     readf,lun,line
     if line ne '#START' then continue
     break
  end
  while not eof(lun) do begin
     readf,lun,line
     strings_local = strsplit(line,',', /EXTRACT)
     filename_sat = strings_local[0]
     
     read_swmf_sat, filename_sat, time_swmf, n_swmf, ux_swmf, uy_swmf,       $
                    uz_swmf, bx_swmf, by_swmf, bz_swmf, ti_swmf, te_swmf,    $
                    ut_swmf, ur_swmf, B_swmf, DoContainData=DoContainData,   $
                    TypeData=TypeData, TypePlot=TypePlot,                    $
                    start_time=start_time, end_time=end_time

     DoLegendIn = 1
     ;; set the legend and DoLegendIn
     if (strings_local[1] eq 'none') then begin
        ;; none means don't write the legend for this line
        DoLegendIn = 0
     endif else if (strings_local[1] eq 'default') then begin
        ;; not so sure whether a user should use default, but just in case...
        if (strpos(filename_sat, 'AWSoM2T') ge 0) then begin
           legendNameIn = 'AWSoM-2T'
        endif else if (strpos(filename_sat, 'AWSoMR') ge 0) then begin
           legendNameIn = 'AWSoM-R'
        endif else begin
           legendNameIn = 'AWSoM'
        endelse
     endif else begin
        legendNameIn = strings_local[1]
     endelse

     ;; set the line thickness and the color
     linethickIn = float(strings_local[2])
     colorIn     = fix(strings_local[3])

     ;; calculate the position for the legend, 0.02112 is from the
     ;; print out value when calling legend in procedures_local...
     nxLegend = nLegendPlotTotal/nyMaxLegend
     nyLegend = nLegendPlotTotal mod nyMaxLegend
     legendPosLIn=[0.15+0.2*nxLegend,0.95-0.02112*nyLegend]
           
     plot_insitu, time_obs, u_obs,  n_obs,  tem_obs, mag_obs,                 $
                  time_swmf, ur_swmf, n_swmf,  ti_swmf,  te_swmf, B_swmf,     $
                  start_time, end_time, typeData=typeData,                    $
                  charsize=CharSizeLocal, DoPlotTe = DoPlotTe,                $
                  legendNames=legendNameIn, DoShowDist=0,                     $
                  ymax_I=[umax_plot,nmax_plot,tmax_plot,bmax_plot],           $
                  IsOverPlot=IsOverPlot, DoLogT=1, linethick=linethickIn,     $
                  colorLocal=colorIn, DoLegend=DoLegendIn,                    $
                  nLegendPlot=nLegendPlot,legendPosL=legendPosLIn

     ;; nLegendPlot is the number of legend printed in plot_insitu
     nLegendPlotTotal = nLegendPlotTotal + nLegendPlot
     IsOverPlot = 1
  end

  device,/close_file
end
