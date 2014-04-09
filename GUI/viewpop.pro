function set_pop_units, stash, onoff
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

function set_normalizex, stash, onoff
	 widget_control, stash, GET_UVALUE=state

 	 state.normalizex = onoff
	 widget_control, stash, SET_UVALUE=state
	 return, state	 
end

function set_logy, stash, onoff
	 widget_control, stash, GET_UVALUE=state

 	 state.logy = onoff
	 widget_control, stash, SET_UVALUE=state
	 return, state	 
end

function set_popdepart, stash, onoff
	 widget_control, stash, GET_UVALUE=state

 	 state.pop_depart = onoff
	 
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

function set_popactive, stash, line
	 widget_control, stash, GET_UVALUE=state

 	 state.drawn(line) = 1
	 state.current = line
	 
	 widget_control, stash, SET_UVALUE=state
	 return, state	 
end

pro plot_popul, state
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
	 
	 case state.pop_depart of 
	 	  0: begin
				  for i = 0, atom.Nlev-1 do begin
		   			if (state.drawn(i) eq 1) then begin
							 yaxis(j,*) = atom.pop(i,*) 
							 j = j + 1
						endif
				  endfor
		     end
	 	  1: begin
				  for i = 0, atom.Nlev-1 do begin
	 		 			if (state.drawn(i) eq 1) then begin
							 yaxis(j,*) = atom.b(i,*) 
							 j = j + 1
						endif
				  endfor
		     end
		  2: begin
				  for i = 0, atom.Nlev-1 do begin
	 		 			if (state.drawn(i) eq 1) then begin
							 yaxis(j,*) = atom.popi(i,*) 
							 j = j + 1
						endif
				  endfor
		     end
		  3: begin
				  for i = 0, atom.Nlev-1 do begin
	 		 			if (state.drawn(i) eq 1) then begin
							 yaxis(j,*) = atom.pop(i,*) / atom.b(i,*) 
							 j = j + 1
						endif
				  endfor
		     end	  
	 endcase	 	  
	 
	 xaxis = geometry.r
	 if (state.normalizex eq 1) then xaxis = xaxis / min(geometry.r)
	 
	 minim = min(yaxis)
	 maxim = max(yaxis)
	 
	 if (state.logx eq 1) then titlex = 'log r'
	 if (state.logx eq 0) then titlex = 'r'
	 if (state.logy eq 1 and state.pop_depart eq 1) then titley = 'log n/n!ILTE!N'
	 if (state.logy eq 1 and state.pop_depart eq 0) then titley = 'log n'
	 if (state.logy eq 0 and state.pop_depart eq 1) then titley = 'n/n!ILTE!N'
	 if (state.logy eq 0 and state.pop_depart eq 0) then titley = 'n'
	 
	 plot, xaxis, yaxis(0,*), XLOG=state.logx, YLOG=state.logy, $
	 	  yrange=[minim,maxim], xtitle=titlex, ytitle=titley
	 for i = 1, ngraph-1 do begin
	 	  oplot, xaxis, yaxis(i,*), linestyle=i
	 endfor
end


pro fill_levels_pop, state
@vars.common
;---------------------------------	
; Fill the level list
;---------------------------------	
	 lines = strarr(atom.Nlev)
	 lines = string(atom.levels)
	 case state.unit of
	 	  0: lines = string(atom.levels)
		  1: lines = string(atom.levels / 2.99792458d10) ; E(cm^-1) = 1/c * E(Hz)
		  2: lines = string(6.62606876d-27 / 1.3806503d-16 * atom.levels) ; E(K)=h/k * E(Hz)
	 endcase
	 
	 widget_control, state.lineList, SET_VALUE=lines
end

pro Viewpops_event, Event
	 
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
					 plot_popul, state
					 cierra_ps
				endif
		  end
		  
		  'LEVELLIST': begin
				k = Event.index				
;		   	if (state.overplot eq 1) then 
				state = set_popactive(stash,k)				
		   	plot_popul, state
		  end
		  
		  'LOGX_OFF': begin	
				plot_popul, set_logx(stash,0)
		  end

		  'LOGX_ON': begin
				plot_popul, set_logx(stash,1)
		  end	
		  
		  'LOGY_OFF': begin	
				plot_popul, set_logy(stash,0)
		  end

		  'LOGY_ON': begin
				plot_popul, set_logy(stash,1)
		  end		  
		  
		  'POPUL_ON': begin				
				plot_popul, set_popdepart(stash,0)
		  end		  
		  
		  'POPUL_OFF': begin
				plot_popul, set_popdepart(stash,1)
		  end	
		  
		  'INITIAL_ON': begin
				plot_popul, set_popdepart(stash,2)
		  end	

		  'LTE_ON': begin
				plot_popul, set_popdepart(stash,3)
		  end			  
		  
		  'OPLOT_ON': begin				
	 	   	state = set_overplot(stash,1)
		  end		  
		  
		  'OPLOT_OFF': begin
				state = set_overplot(stash,0)
		  end		  
		  
		  'RADIUSN_ON': begin				
	 	   	plot_popul, set_normalizex(stash,1)
		  end		  
		  
		  'RADIUSN_OFF': begin
				plot_popul, set_normalizex(stash,0)
		  end		  
		  
		  'CLEAN': begin
				state = set_overplot(stash,0, /CLEAN)
				state.current = 0
				plot_popul, state
		  end		  
		  
		  'HZ': begin
	 	   	fill_levels_pop, set_pop_units(stash, 0)
		  end		  
		  
		  'CM': begin
	 	   	fill_levels_pop, set_pop_units(stash, 1)				
		  end		  
		  
		  'KELVIN': begin
	 	   	fill_levels_pop, set_pop_units(stash, 2)				
		  end		  		  		  
	 
	 endcase
end

function viewpop_init
@vars.common
	 state = {baseWidget: 0L, drawWidget: 0L, lineList: 0L, logx: 0, logy: 0, pop_depart: 0, $ 
	 	  current: 0, drawn: intarr(atom.Nlev), overplot: 0, normalizex: 1, unit: 0, pop_lvg: 0}
	 
;---------------------------------	
; Main window
;---------------------------------	
	 state.baseWidget = widget_base(TITLE='View populations', /ROW, MBAR=menuBar, $
                                  RESOURCE_NAME='Viewpop')
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
	 drawFrame = widget_base(populBase, /FRAME)
    state.drawWidget = widget_draw(drawFrame, XSIZE=600, YSIZE=450, $
                                 /BUTTON_EVENTS, UVALUE='SET_SCALE')
											
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
; Populations/departure coefficients
;---------------------------------	
	 pop_departFrame = widget_base(botones1, /COLUMN, /EXCLUSIVE)
	 populButton = widget_button(pop_departFrame, VALUE='Population', $
                               UVALUE='POPUL_ON')
	 departButton = widget_button(pop_departFrame, VALUE='Departure', $
                               UVALUE='POPUL_OFF')										 										 
	 initialButton = widget_button(pop_departFrame, VALUE='Initial', $
                               UVALUE='INITIAL_ON')	
	 lteButton = widget_button(pop_departFrame, VALUE='LTE', $
                               UVALUE='LTE_ON')											 
										 										 
	 widget_control, (state.pop_depart) ? departButton : populButton, /SET_BUTTON		

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
	 cleanFrame = widget_base(botones1, /COLUMN)
	 cleanButton = widget_button(cleanFrame, VALUE='Clean', $
                               UVALUE='CLEAN')	 
;---------------------------------	
; Level list
;---------------------------------	
	 transFrame = widget_base( populBase, /ROW, /FRAME )
	 lineFormat = 'Level (Hz)'
	 lineFrame  = widget_base(transFrame, /COLUMN)
	 lineLabel  = widget_label(lineFrame, VALUE='        Energy levels        ')
	 state.lineList = widget_list(lineFrame, UVALUE='LEVELLIST', $
                                   RESOURCE_NAME='list', $
                                   YSIZE=10, VALUE=lineFormat)		
											  
	 unitsFrame = widget_base(transFrame, /COLUMN, /EXCLUSIVE)
	 hzButton = widget_button(unitsFrame, VALUE='Hz', UVALUE='HZ')
	 cmButton = widget_button(unitsFrame, VALUE='cm^-1', UVALUE='CM')
	 kelvinButton = widget_button(unitsFrame, VALUE='Kelvin', UVALUE='KELVIN')
	 widget_control, hzButton, /SET_BUTTON	 
											  
	 fill_levels_pop, state						 

	 widget_control, widget_info(state.baseWidget, /CHILD), SET_UVALUE=state
	 
	 return, state
end


pro viewpop, LEADER=leader
	 if (NOT keyword_set(LEADER)) then leader=0
	 
	 state = viewpop_init()
	 
	 widget_control, state.baseWidget, /REALIZE, GROUP_LEADER=leader
	 
	 xmanager, 'Viewpops', state.baseWidget, $
    	  EVENT_HANDLER='Viewpops_event', GROUP_LEADER=leader	 
end
