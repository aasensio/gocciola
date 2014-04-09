function T_raylj, T, nu
	 f = 6.62606876d-27 / 1.3806503d-16 * nu
	 TRJ =  f / (exp(f/T)-1.d0)	 
	 return, TRJ
end 
function set_source_units, stash, onoff
    widget_control, stash, GET_UVALUE=state

    state.unit = onoff
    widget_control, stash, SET_UVALUE=state
    return, state   
end

function set_logx, stash, onoff
	 widget_control, stash, GET_UVALUE=state

 	 state.logx = onoff
	 widget_control, stash, SET_UVALUE=state
	 return, state	 
end

function set_logy, stash, onoff
	 widget_control, stash, GET_UVALUE=state

 	 state.logy = onoff
	 widget_control, stash, SET_UVALUE=state
	 return, state	 
end

function set_sourtexc, stash, onoff
	 widget_control, stash, GET_UVALUE=state

 	 state.sour_texc = onoff
	 
	 widget_control, stash, SET_UVALUE=state
	 return, state	 
end

function set_normalizey, stash, onoff
	 widget_control, stash, GET_UVALUE=state

 	 state.normalizey = onoff
	 
	 widget_control, stash, SET_UVALUE=state
	 return, state	 
end

function set_normalizex, stash, onoff
	 widget_control, stash, GET_UVALUE=state

 	 state.normalizex = onoff
	 
	 widget_control, stash, SET_UVALUE=state
	 return, state	 
end

function set_overplot, stash, onoff, CLEAN=clean
	 widget_control, stash, GET_UVALUE=state

 	 if (not keyword_set(CLEAN)) then state.overplot = onoff
	 if (onoff eq 0) then state.drawn = 0
	 
	 widget_control, stash, SET_UVALUE=state
	 return, state	 
end

function set_active, stash, line
	 widget_control, stash, GET_UVALUE=state

 	 state.drawn(line) = 1
	 
	 widget_control, stash, SET_UVALUE=state
	 return, state	 
end

pro plot_source, state
@vars.common
	 k = state.current
	 if (state.overplot eq 0) then begin
	 	  state.drawn = 0
		  state.drawn(k) = 1
	 endif
	 ngraph = n_elements(where(state.drawn eq 1))
	 xaxis = dblarr(geometry.Nradius)
	 yaxis = dblarr(ngraph,geometry.Nradius)
	 
	 j = 0
	 
	 if (state.sour_texc eq 0) then begin
	 	  for i = 0, transition.Ntran-1 do begin
		   	if (state.drawn(i) eq 1) then begin
					 yaxis(j,*) = transition.S(i,*) 
					 if (state.normalizey eq 0) then yaxis(j,*) = yaxis(j,*) * transition.B(i,*)
					 j = j + 1
				endif
		  endfor
	 endif else begin
	 	  for i = 0, transition.Ntran-1 do begin
	 		 	if (state.drawn(i) eq 1) then begin
					 yaxis(j,*) = transition.texc(i,*) 
					 j = j + 1
				endif
		  endfor
	 endelse
	 
	 xaxis = geometry.r
	 if (state.normalizex eq 1) then xaxis = xaxis / min(geometry.r)
	 
	 minim = min(yaxis)
	 maxim = max(yaxis)
	 
	 if (state.logx eq 1) then titlex = 'log r'
	 if (state.logx eq 0) then titlex = 'r'
	 if (state.logy eq 1 and state.sour_texc eq 1) then titley = 'log T!Iexc!N'
	 if (state.logy eq 1 and state.sour_texc eq 0) then titley = 'log S'
	 if (state.logy eq 0 and state.sour_texc eq 1) then titley = 'T!Iexc!N'
	 if (state.logy eq 0 and state.sour_texc eq 0) then titley = 'S'
	 
	 plot, xaxis, yaxis(0,*), XLOG=state.logx, YLOG=state.logy, $
	 	  yrange=[minim,maxim], xtitle=titlex, ytitle=titley
	 for i = 1, ngraph-1 do begin
	 	  oplot, xaxis, yaxis(i,*), linestyle=i
	 endfor
	 if (state.normalizey eq 0) then begin
		 if (state.sour_texc eq 0) then begin
	 		  oplot, xaxis, transition.B(k,*), linestyle=2
		 endif else begin
	 		  temp = 1.d0 + 1.4744991554d-47 * transition.freq(k)^3 / transition.B(k,*)
			  temp = 4.7992375477d-11 * transition.freq(k) / alog(temp)		  
	 		  oplot, xaxis, temp, linestyle=2
		 endelse
	 endif
	 
	 tau1 = interpol(xaxis,transition.tau(k,*),1.d0)
	 arrowsize = (!Y.CRANGE[1]-!Y.CRANGE[0]) / 25.

	 ind = cerca(transition.tau(k,*),tau1(0))

	 arrow, tau1, yaxis(ind)+arrowsize, tau1, yaxis(ind), /data, thick=2.
	 
	 total_opacity = transition.tau(k,0);total(transition.tau(k,*),2)
	 texto_data = ' Total line opacity : '+strtrim(string(total_opacity),2)+string(10B)+string(13B)
	 
; Calculate Rayleigh-Jeans brightness temperature
; Tb=Tbg(RJ)*exp(-tau)+Texc(RJ)*(1-exp(-tau))-Tbg(RJ)
; where T(RJ)=c^2/(2*k*nu^2) * B_nu(T)
	 tbackground = 2.7296d0
	 texc = mean(transition.texc(k,0:geometry.Nradius-2))
	 b1 = t_raylj(texc, transition.freq(k)) - t_raylj(tbackground, transition.freq(k))
	 t_brightness = (1.d0 - exp(-total_opacity)) * b1
	 texto_data = texto_data + 'Rayl-J brightness temperature : '+$
	 	  strtrim(string(t_brightness),2)+string(10B)+string(13B)
	 
	 dvdr = mean(deriv(atmos.r,atmos.v_mac)) * 3.085678d18
	 zco = mean(atmos.n)
	 X = 3.08567d18 * mean(atmos.nh2) / (dvdr * t_brightness)
	 
	 n = n_elements(atmos.r)
	 vturb = sqrt(2.d0 * 1.381d-16 * mean(atmos.T) / (28.*1.66d-24) + 1.d10) / 1.d5
	 X = (atmos.r(n-1)-atmos.r(0)) * mean(atmos.nh2) / (t_brightness * vturb) 
	 X = 0.94d20 * sqrt(mean(atmos.nh2)) / t_brightness
	 texto_data = texto_data + 'CO-to-H2 conversion : '+strtrim(string(X),2)
	 
	 widget_control, state.texto, set_value=texto_data
	 

end


pro fill_lines_source, state
@vars.common
;---------------------------------	
; Fill the level list
;---------------------------------	
	 lines = strarr(transition.Ntran)
	 case state.unit of
   	 0: lines = strtrim(string(transition.up),2)+' ->  '+strtrim(string(transition.low),2)+$
		  '   :      '+strtrim(string(transition.freq),2)
   	 1: lines = strtrim(string(transition.up),2)+' ->  '+strtrim(string(transition.low),2)+$
		  '   :      '+strtrim(string(transition.freq / 2.99792458d10),2) ; E(cm^-1) = 1/c * E(Hz)
   	 2: lines = strtrim(string(transition.up),2)+' ->  '+strtrim(string(transition.low),2)+$
	 	  '   :      '+strtrim(string(6.62606876d-27 / 1.3806503d-16 * transition.freq),2) ; E(K)=h/k * E(Hz)
   	 3: lines = strtrim(string(transition.up),2)+' ->  '+strtrim(string(transition.low),2)+$
	 	  '   :      '+strtrim(string(2.99792458d14 / transition.freq),2) ; E(microns)= 1.d6 * c / E(Hz)		  
   	 4: lines = strtrim(string(transition.up),2)+' ->  '+strtrim(string(transition.low),2)+$
	 	  '   :      '+strtrim(string(2.99792458d18 / transition.freq),2) ; E(A) = 1.d10 * c / E(Hz) 		  
    endcase
	 widget_control, state.lineList, SET_VALUE=lines
end

pro Viewsource_event, Event
	 
	 stash = widget_info( Event.handler, /CHILD)
	 widget_control, stash, GET_UVALUE=state

	 widget_control, Event.id, GET_UVALUE=Action
	 	 
	 case Action of
	 	  'QUIT': begin
		   	widget_control, Event.top, /DESTROY
		  end
		  
		  'POSTCRIPT': begin	
		   	name = dialog_pickfile()
				if (name ne '') then begin
					 abre_ps, name, /landscape
					 plot_source, state
					 cierra_ps
				endif
		  end
		  
		  'LEVELLIST': begin
				state.current = Event.index
		   	if (state.overplot eq 1) then state = set_active(stash,state.current)				
		   	plot_source, state
		  end
		  
		  'LOGX_OFF': begin	
				plot_source, set_logx(stash,0)
		  end

		  'LOGX_ON': begin
				plot_source, set_logx(stash,1)
		  end	
		  
		  'LOGY_OFF': begin	
				plot_source, set_logy(stash,0)
		  end

		  'LOGY_ON': begin
				plot_source, set_logy(stash,1)
		  end		  
		  
		  'SOURCE_ON': begin				
				plot_source, set_sourtexc(stash,0)
		  end		  
		  
		  'SOURCE_OFF': begin
				plot_source, set_sourtexc(stash,1)
		  end	
		  
		  'OPLOT_ON': begin				
	 	   	state = set_overplot(stash,1)
		  end		  
		  
		  'OPLOT_OFF': begin
				state = set_overplot(stash,0)
		  end		  
		  
		  'NORMALIZ_ON': begin				
	 	   	plot_source, set_normalizey(stash,1)
		  end		  
		  
		  'NORMALIZ_OFF': begin
				plot_source, set_normalizey(stash,0)
		  end		
		  
		  'RADIUSN_ON': begin				
	 	   	plot_source, set_normalizex(stash,1)
		  end		  
		  
		  'RADIUSN_OFF': begin
				plot_source, set_normalizex(stash,0)
		  end		  
		  
		  'CLEAN': begin
				state = set_overplot(stash,0, /CLEAN)
				state.current = 0
				plot_source, state
		  end		  
	 	  'HZ': begin
              fill_lines_source, set_source_units(stash, 0)
        end             

        'CM': begin
              fill_lines_source, set_source_units(stash, 1)
        end             

        'KELVIN': begin
              fill_lines_source, set_source_units(stash, 2)
        end         		  
		  
	 	  'MICRONS': begin
              fill_lines_source, set_source_units(stash, 3)
        end         		  
		  
	 	  'ANGSTROMS': begin
              fill_lines_source, set_source_units(stash, 4)
        end         		  		  
	 
	 endcase
end

function viewsource_init
@vars.common
	 state = {baseWidget: 0L, drawWidget: 0L, lineList: 0L, logx: 0, logy: 0, sour_texc: 0, $ 
	 	  current: 0, drawn: intarr(transition.Ntran), overplot: 0, normalizey: 1, normalizex: 1, $
		  unit: 0, texto: 0L}
	 
;---------------------------------	
; Main window
;---------------------------------	
	 state.baseWidget = widget_base(TITLE='Source function / excitation temperature', /ROW, $
	 	  MBAR=menuBar, RESOURCE_NAME='Viewsource')
	 populBase = widget_base(state.baseWidget, /COLUMN)

;---------------------------------	
; File menu
;---------------------------------	
	 fileMenu  = widget_button(menuBar, VALUE='File', /MENU)

;---------------------------------	
; Quit button
;---------------------------------	
	 quitButton = widget_button(fileMenu, VALUE='Quit', UVALUE='QUIT', $
                             RESOURCE_NAME='quitbutton')

;---------------------------------	
; Postcript button
;---------------------------------	
	 postcriptButton = widget_button(fileMenu, VALUE='Postcript', UVALUE='POSTCRIPT', $
                             RESOURCE_NAME='postcript')									  
	 
;---------------------------------	
; Drawing window
;---------------------------------	
	 bigFrame = widget_base(populBase, /ROW)
	 drawFrame = widget_base(bigFrame, /ROW) ;/FRAME)
    state.drawWidget = widget_draw(drawFrame, XSIZE=600, YSIZE=450, $
                                 /BUTTON_EVENTS, UVALUE='SET_SCALE')
;---------------------------------	
; Text window
;---------------------------------												
	 msgWidget = widget_base(bigFrame, /ROW)
	 state.texto = widget_text(msgWidget, value='', XSIZE=50, YSIZE=10)
											
	 botones1 = widget_base(populBase, /FRAME, /ROW)

;---------------------------------	
; Linear/logaritmic axis for X	 
;---------------------------------	
	 linlogxFrame = widget_base(botones1, /COLUMN, /EXCLUSIVE)
	 linxButton = widget_button(linlogxFrame, VALUE='Linear X', $
                               UVALUE='LOGX_OFF')
	 logxButton = widget_button(linlogxFrame, VALUE='Log X', $
                               UVALUE='LOGX_ON')	
	 widget_control, (state.logx) ? logxButton : linxButton, /SET_BUTTON										 									 
	 
;---------------------------------	
; Linear/logaritmic axis for X	 
;---------------------------------	
	 linlogyFrame = widget_base(botones1, /COLUMN, /EXCLUSIVE)
	 linyButton = widget_button(linlogyFrame, VALUE='Linear Y', $
                               UVALUE='LOGY_OFF')
	 logyButton = widget_button(linlogyFrame, VALUE='Log Y', $
                               UVALUE='LOGY_ON')										 										 
	 widget_control, (state.logy) ? logyButton : linyButton, /SET_BUTTON										 
	 
;---------------------------------	
; Source function/excitation temperature
;---------------------------------	
	 sour_texcFrame = widget_base(botones1, /COLUMN, /EXCLUSIVE)
	 sourButton = widget_button(sour_texcFrame, VALUE='Source', $
                               UVALUE='SOURCE_ON')
	 texcButton = widget_button(sour_texcFrame, VALUE='Excitation temperature', $
                               UVALUE='SOURCE_OFF')										 										 
	 widget_control, (state.sour_texc) ? texcButton : sourButton, /SET_BUTTON		

;---------------------------------	
; Normalize to Planck's function
;---------------------------------	
	 normalizFrame = widget_base(botones1, /COLUMN, /EXCLUSIVE)
	 normalizButton = widget_button(normalizFrame, VALUE='Normalize (Planck`s)', $
                               UVALUE='NORMALIZ_ON')
	 notnormalizButton = widget_button(normalizFrame, VALUE='Not normalize', $
                               UVALUE='NORMALIZ_OFF')										 										 
	 widget_control, (state.normalizey) ? normalizButton : notnormalizButton, /SET_BUTTON		
	 
;---------------------------------	
; Plot/overplot	 
;---------------------------------	
	 plotoplotFrame = widget_base(botones1, /COLUMN, /EXCLUSIVE)
	 plotButton = widget_button(plotoplotFrame, VALUE='Plot', $
                               UVALUE='OPLOT_OFF')
	 oplotButton = widget_button(plotoplotFrame, VALUE='Overplot', $
                               UVALUE='OPLOT_ON')										 										 
	 widget_control, (state.overplot) ? oplotButton : plotButton, /SET_BUTTON										 	 

;---------------------------------	
; Radius normalization
;---------------------------------	
	 radiusFrame = widget_base(botones1, /COLUMN, /EXCLUSIVE)
	 radnorButton = widget_button(radiusFrame, VALUE='Radius normalization', $
                               UVALUE='RADIUSN_ON')
	 radnotnorButton = widget_button(radiusFrame, VALUE='Radius', $
                               UVALUE='RADIUSN_OFF')										 										 
	 widget_control, (state.normalizex) ? radnorButton : radnotnorButton, /SET_BUTTON										 	 	 

;---------------------------------	
; Cleaning
;---------------------------------		 
	 cleanFrame = widget_base(botones1, /ROW)
	 cleanButton = widget_button(cleanFrame, VALUE='Clean', $
                               UVALUE='CLEAN')	 
;---------------------------------	
; Level list
;---------------------------------	
	 transFrame = widget_base( populBase, /ROW, /FRAME )
	 lineFormat = 'Transition (!3Å)'
	 lineFrame  = widget_base(transFrame, /COLUMN)
	 lineLabel  = widget_label(lineFrame, VALUE='Line list        ')
	 state.lineList = widget_list(lineFrame, UVALUE='LEVELLIST', $
                                   RESOURCE_NAME='list', $
                                   YSIZE=10, XSIZE=18, VALUE=lineFormat)		
	 
	 unitsFrame = widget_base(transFrame, /COLUMN, /EXCLUSIVE)
    hzButton = widget_button(unitsFrame, VALUE='Hz', UVALUE='HZ')
    cmButton = widget_button(unitsFrame, VALUE='cm^-1', UVALUE='CM')
    kelvinButton = widget_button(unitsFrame, VALUE='Kelvin', UVALUE='KELVIN')
	 micronsButton = widget_button(unitsFrame, VALUE='Micron', UVALUE='MICRONS')
	 angsButton = widget_button(unitsFrame, VALUE='Angstrom', UVALUE='ANGSTROMS')
    widget_control, hzButton, /SET_BUTTON 
														  
	 fill_lines_source, state						 

	 widget_control, widget_info(state.baseWidget, /CHILD), SET_UVALUE=state
	 
	 return, state
end


pro viewsource, LEADER=leader
	 if (NOT keyword_set(LEADER)) then leader=0
	 
	 state = viewsource_init()
	 
	 widget_control, state.baseWidget, /REALIZE, GROUP_LEADER=leader
	 
	 xmanager, 'Viewsource', state.baseWidget, $
    	  EVENT_HANDLER='Viewsource_event', GROUP_LEADER=leader	 
end
