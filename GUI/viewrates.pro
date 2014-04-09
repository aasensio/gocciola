;---------------------------------
; This function sets the Y axis with logarithm scale
;---------------------------------
function set_logy_rates, stash, onoff
	 widget_control, stash, GET_UVALUE=state

 	 state.logy = onoff
	 widget_control, stash, SET_UVALUE=state
	 return, state	 
end

;---------------------------------
; This function sets the X axis with logarithm scale
;---------------------------------
function set_logx_rates, stash, onoff
	 widget_control, stash, GET_UVALUE=state

 	 state.logx = onoff
	 widget_control, stash, SET_UVALUE=state
	 return, state	 
end

;---------------------------------
; This function sets the upward or downward rates
;---------------------------------
function set_updown_rates, stash, onoff
	 widget_control, stash, GET_UVALUE=state

 	 state.updown = onoff
	 widget_control, stash, SET_UVALUE=state
	 return, state	 
end

;---------------------------------
; This function sets the current line being plotted
;---------------------------------
function set_line_rates, stash, onoff
	 widget_control, stash, GET_UVALUE=state

 	 state.current = onoff
	 widget_control, stash, SET_UVALUE=state
	 return, state	 
end

;---------------------------------
; This function plots the radiative rates
;---------------------------------
pro plot_rates, state
@vars.common
	 k = state.current	 	
	 	 	 
	 xaxis = geometry.r
	 if (state.updown eq 1) then begin
	 	  titley = 'R!Ilu!N (s!E-1!N)'
	 	  yaxis = transition.Rlu(k,*)		  
	 endif else begin
	     titley = 'R!Iul!N (s!E-1!N)'		  		  
	 	  yaxis = transition.Rul(k,*)
	 endelse
	 titlex = 'Radius (km)'	 
	 	 	 	 
	 plot, xaxis, yaxis, YLOG=state.logy, xtitle=titlex, ytitle=titley, XLOG=state.logx
end


;---------------------------------	
; This function fills the line list
;---------------------------------	
pro fill_lines, state
@vars.common

	 lines = strarr(transition.Ntran)
	 lines = string(transition.up)+'    '+string(transition.low)
	 widget_control, state.lineList, SET_VALUE=lines
end

pro Viewrates_event, Event
	 
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
	 				 plot_rates, state
		   		 cierra_ps
				endif
		  end
		  
		  'LEVELLIST': begin
		   	plot_rates, set_line_rates(stash, Event.index)
		  end
		  		  
		  'LOGY_OFF': begin					
		   	plot_rates, set_logy_rates(stash,0)
		  end

		  'LOGY_ON': begin		  		
				plot_rates, set_logy_rates(stash,1)			
		  end		
		  
		  'LOGX_OFF': begin					
		   	plot_rates, set_logx_rates(stash,0)
		  end

		  'LOGX_ON': begin		  		
				plot_rates, set_logx_rates(stash,1)			
		  end		
		  		  
		  'UP': begin			   	
				plot_rates, set_updown_rates(stash,1)
		  end

		  'DOWN': begin		   	
				plot_rates, set_updown_rates(stash,0)
		  end		  		  
		  	 
	 endcase
end

function viewrates_init
@vars.common
	 state = {baseWidget: 0L, drawWidget: 0L, lineList: 0L, logy: 0, updown: 0, $ 
	 	  current: 0, logx: 0}
	 
;---------------------------------	
; Main window
;---------------------------------	
	 state.baseWidget = widget_base(TITLE='Radiative rates', /ROW, $
	 	  MBAR=menuBar, RESOURCE_NAME='ViewRates')
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
; Linear/logaritmic axis for Y
;---------------------------------	
	 linlogyFrame = widget_base(botones1, /COLUMN, /EXCLUSIVE)
	 linyButton = widget_button(linlogyFrame, VALUE='Linear Y', $
                               UVALUE='LOGY_OFF')
	 logyButton = widget_button(linlogyFrame, VALUE='Log Y', $
                               UVALUE='LOGY_ON')										 										 
	 widget_control, (state.logy) ? logyButton : linyButton, /SET_BUTTON										 

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
; Upward/downward rates
;---------------------------------	
	 up_downFrame = widget_base(botones1, /COLUMN, /EXCLUSIVE)
	 upButton = widget_button(up_downFrame, VALUE='Upward', $
                               UVALUE='UP')
	 downButton = widget_button(up_downFrame, VALUE='Downward', $
                               UVALUE='DOWN')										 										 
	 widget_control, (state.updown) ? upButton : downButton, /SET_BUTTON		

;---------------------------------	
; Level list
;---------------------------------	
	 transFrame = widget_base( populBase, /ROW, /FRAME )
	 lineFormat = 'Transition (!3Å)'
	 lineFrame  = widget_base(transFrame, /COLUMN)
	 lineLabel  = widget_label(lineFrame, VALUE='Line list (!3Å)       ')
	 state.lineList = widget_list(lineFrame, UVALUE='LEVELLIST', $
                                   RESOURCE_NAME='list', $
                                   YSIZE=10, VALUE=lineFormat)		
	 											  
	 fill_lines, state						 

	 widget_control, widget_info(state.baseWidget, /CHILD), SET_UVALUE=state
	 
	 return, state
end


pro viewradrates, LEADER=leader
	 if (NOT keyword_set(LEADER)) then leader=0
	 
	 state = viewrates_init()
	 
	 widget_control, state.baseWidget, /REALIZE, GROUP_LEADER=leader
	 
	 xmanager, 'ViewRates', state.baseWidget, $
    	  EVENT_HANDLER='Viewrates_event', GROUP_LEADER=leader	 
end
