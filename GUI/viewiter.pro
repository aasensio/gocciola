function set_logy, stash, onoff
	 widget_control, stash, GET_UVALUE=state

 	 state.logy = onoff
	 widget_control, stash, SET_UVALUE=state
	 return, state	 
end

pro plot_iter, state
@vars.common
	 xaxis = indgen(iterat.Niter)
	 yaxis = iterat.rel_change
	 
	 titley = 'R!Ic!N'
	 titlex = 'Iteration'
	 if (state.logy eq 1) then titley = 'log R!Ic!N'
	 
	 plot, xaxis, yaxis, YLOG=state.logy, xtitle=titlex, ytitle=titley
	 
end

pro Viewiter_event, Event
	 
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
				plot_iter, set_logy(stash,0)
		  end

		  'LOGY_ON': begin
				plot_iter, set_logy(stash,1)
		  end		  		  
	 
	 endcase
end

function viewiter_init
@vars.common
	 state = {baseWidget: 0L, drawWidget: 0L, logy: 0}
	 
;---------------------------------	
; Main window
;---------------------------------	
	 state.baseWidget = widget_base(TITLE='View iterations', /ROW, MBAR=menuBar, $
                                  RESOURCE_NAME='Viewiter')
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
	 
	 widget_control, widget_info(state.baseWidget, /CHILD), SET_UVALUE=state
	 	 
	 return, state
end


pro viewiter, LEADER=leader
	 if (NOT keyword_set(LEADER)) then leader=0
	 
	 state = viewiter_init()
	 
	 widget_control, state.baseWidget, /REALIZE, GROUP_LEADER=leader
	 
	 plot_iter, state
	 
	 xmanager, 'Viewiter', state.baseWidget, $
    	  EVENT_HANDLER='Viewiter_event', GROUP_LEADER=leader	 
end
