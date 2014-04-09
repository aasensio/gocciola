function set_logy, stash, onoff
	 widget_control, stash, GET_UVALUE=state

 	 state.logy = onoff
	 widget_control, stash, SET_UVALUE=state
	 return, state	 
end

function set_tran, stash, onoff
	 widget_control, stash, GET_UVALUE=state

 	 state.k = onoff
	 widget_control, stash, SET_UVALUE=state
	 return, state	 
end

pro plot_back, state
@vars.common

	 xaxis = geometry.r / 1.d5
	 yaxis = atom.backopa(state.k,*)
	 
	 titley = '!7v!3!Ic!N (cm!E-1!N)'
	 titlex = 'r [km]'	 
	 
	 plot, xaxis, yaxis, YLOG=state.logy, xtitle=titlex, ytitle=titley
	 
end

pro Viewback_event, Event
	 
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
					 plot_iter, state
					 cierra_ps
				endif
		  end
		  		  
		  'LOGY_OFF': begin	
				plot_back, set_logy(stash,0)
		  end

		  'LOGY_ON': begin
				plot_back, set_logy(stash,1)
		  end
		  
		  'LEVELLIST': begin	
				plot_back, set_tran(stash,Event.index)
		  end		  		  
	 
	 endcase
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


function viewback_init
@vars.common
	 state = {baseWidget: 0L, drawWidget: 0L, logy: 0, k:0, lineList:0L}
	 
;---------------------------------	
; Main window
;---------------------------------	
	 state.baseWidget = widget_base(TITLE='View background opacity', /ROW, MBAR=menuBar, $
                                  RESOURCE_NAME='Viewback')
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
	 linlogyFrame = widget_base(botones1, /COLUMN, /EXCLUSIVE)
	 linyButton = widget_button(linlogyFrame, VALUE='Linear Y', $
                               UVALUE='LOGY_OFF')
	 logyButton = widget_button(linlogyFrame, VALUE='Log Y', $
                               UVALUE='LOGY_ON')										 										 
	 widget_control, (state.logy) ? logyButton : linyButton, /SET_BUTTON										 
	 
;---------------------------------	
; Transition list 
;---------------------------------	
	 transFrame = widget_base( populBase, /ROW, /FRAME )
    lineFormat = 'Transition (!3Å)'
    lineFrame  = widget_base(transFrame, /COLUMN)
    lineLabel  = widget_label(lineFrame, VALUE='Line list        ')
    state.lineList = widget_list(lineFrame, UVALUE='LEVELLIST', $
                              RESOURCE_NAME='list', $
                              YSIZE=10, VALUE=lineFormat)	 
	 
	 widget_control, widget_info(state.baseWidget, /CHILD), SET_UVALUE=state
	 
	 fill_lines, state
	 	 
	 return, state
end


pro viewback, LEADER=leader
	 if (NOT keyword_set(LEADER)) then leader=0
	 
	 state = viewback_init()
	 
	 widget_control, state.baseWidget, /REALIZE, GROUP_LEADER=leader
	 
	 plot_back, state
	 
	 xmanager, 'Viewbackopacity', state.baseWidget, $
    	  EVENT_HANDLER='Viewback_event', GROUP_LEADER=leader	 
end
