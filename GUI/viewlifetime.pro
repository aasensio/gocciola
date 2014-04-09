function set_lifetime_logx, stash, onoff
	 widget_control, stash, GET_UVALUE=state

 	 state.logx = onoff
	 widget_control, stash, SET_UVALUE=state
	 return, state	 
end

function set_lifetime_normalizex, stash, onoff
	 widget_control, stash, GET_UVALUE=state

 	 state.normalizex = onoff
	 widget_control, stash, SET_UVALUE=state
	 return, state	 
end

function set_lifetime_logy, stash, onoff
	 widget_control, stash, GET_UVALUE=state

 	 state.logy = onoff
	 widget_control, stash, SET_UVALUE=state
	 return, state	 
end

pro plot_lifetime, state
@vars.common
	 k = state.current

	 xaxis = dblarr(geometry.Nradius)
	 yaxis = dblarr(geometry.Nradius)
	 emision = fltarr(geometry.Nradius)
	 absorcion = fltarr(geometry.Nradius)
	 
	 for i = 0, transition.Ntran-1 do begin
	 	  l = transition.low(i) - 1
		  u = transition.up(i) - 1

		  if (k eq u) then emision = emision + transition.Rul(i,*)
		  if (k eq l) then absorcion = absorcion + transition.Rlu(i,*)
	 endfor
	 
	 yaxis = 1.d0 / (emision + absorcion)
	 
	 xaxis = geometry.r
	 if (state.normalizex eq 1) then xaxis = xaxis / min(geometry.r)
	 
	 minim = min(yaxis)
	 maxim = max(yaxis)
	 
	 if (state.logx eq 1) then titlex = 'log r'
	 if (state.logx eq 0) then titlex = 'r'
	 if (state.logy eq 1) then titley = 'log t'
	 if (state.logy eq 0) then titley = 't (s)'
	 
	 plot, xaxis, yaxis, XLOG=state.logx, YLOG=state.logy, $
	 	  yrange=[minim,maxim], xtitle=titlex, ytitle=titley	 
end


pro fill_levels_lifetime, state
@vars.common
;---------------------------------	
; Fill the level list
;---------------------------------	
	 lines = strarr(atom.Nlev)
	 lines = string(atom.levels)
	 widget_control, state.lineList, SET_VALUE=lines
end

pro Viewlifetime_event, Event
	 
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
					 plot_lifetime, state
					 cierra_ps
				endif
		  end
		  
		  'LEVELLIST': begin
				state.current = Event.index		   
		   	plot_lifetime, state
		  end
		  
		  'LOGX_OFF': begin	
				plot_lifetime, set_lifetime_logx(stash,0)
		  end

		  'LOGX_ON': begin
				plot_lifetime, set_lifetime_logx(stash,1)
		  end	
		  
		  'LOGY_OFF': begin	
				plot_lifetime, set_lifetime_logy(stash,0)
		  end

		  'LOGY_ON': begin
				plot_lifetime, set_lifetime_logy(stash,1)
		  end		  
		  		  
		  'RADIUSN_ON': begin				
	 	   	plot_lifetime, set_lifetime_normalizex(stash,1)
		  end		  
		  
		  'RADIUSN_OFF': begin
				plot_lifetime, set_lifetime_normalizex(stash,0)
		  end		  
		  	 
	 endcase
end

function viewlifetime_init
@vars.common
	 state = {baseWidget: 0L, drawWidget: 0L, lineList: 0L, logx: 0, logy: 0, current: 0, $
	 normalizex: 1}
	 
;---------------------------------	
; Main window
;---------------------------------	
	 state.baseWidget = widget_base(TITLE='View lifetimes', /ROW, MBAR=menuBar, $
                                  RESOURCE_NAME='Viewlifetime')
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
; Linear/logaritmic axis for Y
;---------------------------------	
	 linlogyFrame = widget_base(botones1, /COLUMN, /EXCLUSIVE)
	 linyButton = widget_button(linlogyFrame, VALUE='Linear Y', $
                               UVALUE='LOGY_OFF')
	 logyButton = widget_button(linlogyFrame, VALUE='Log Y', $
                               UVALUE='LOGY_ON')										 										 
	 widget_control, (state.logy) ? logyButton : linyButton, /SET_BUTTON										 
	 
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
; Level list
;---------------------------------	
	 transFrame = widget_base( populBase, /ROW, /FRAME )
	 lineFormat = 'Level (Hz)'
	 lineFrame  = widget_base(transFrame, /COLUMN)
	 lineLabel  = widget_label(lineFrame, VALUE='        Energy levels        ')
	 state.lineList = widget_list(lineFrame, UVALUE='LEVELLIST', $
                                   RESOURCE_NAME='list', $
                                   YSIZE=10, VALUE=lineFormat)		
											  
	 fill_levels_lifetime, state						 

	 widget_control, widget_info(state.baseWidget, /CHILD), SET_UVALUE=state
	 
	 return, state
end


pro viewlifetime, LEADER=leader
	 if (NOT keyword_set(LEADER)) then leader=0
	 
	 state = viewlifetime_init()
	 
	 widget_control, state.baseWidget, /REALIZE, GROUP_LEADER=leader
	 
	 xmanager, 'Viewlifetime', state.baseWidget, $
    	  EVENT_HANDLER='Viewlifetime_event', GROUP_LEADER=leader	 
end
