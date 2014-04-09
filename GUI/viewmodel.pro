function set_norm_model, stash, onoff
	 widget_control, stash, GET_UVALUE=state

 	 state.normalize = onoff
	 widget_control, stash, SET_UVALUE=state
	 return, state	 
end

pro plotmodel, state
@vars.common
	 openr, 2, '../'+strtrim(atom.model,2)
	 lb, 2, 22
	 readf, 2, n_levels
	 E = fltarr(n_levels)
	 g = fltarr(n_levels)
	 orig = fltarr(n_levels)
	 lb, 2, 5
	 readf, 2, n_transitions
	 low = fltarr(n_transitions)
	 up = low
	 A = low
	 lb, 2, 5
	 temp1 = 0.
	 temp2 = 0
	 temp3 = 0.
	 temp4 = ''
	 xoff = 1
	 yoff = 0

	 for i = 0, n_levels-1 do begin
	 	  readf, 2, temp1, temp2, temp3, temp4
		  E(i) = temp1 - yoff
		  orig(i) = xoff
		  g(i) = temp2  
		  if ( (i ne 0) and (g(i) eq 1)) then begin
				xoff = xoff + 1
				orig(i) = xoff				
				if (state.normalize eq 1) then begin
					 yoff = yoff + E(i) - E(i-1) - 20
					 E(i) = E(i) - (E(i) - E(i-1)) + 20
				endif
		  endif
		  
	 endfor
	 
	 lb, 2, 2
	 temp1 = 0
	 temp2 = 0
	 temp3 = 0.
	 for i = 0, n_transitions-1 do begin
	 	  readf, 2, temp1, temp2, temp3
		  low(i) = temp2-1
		  up(i) = temp1-1
		  A(i) = temp3
	 endfor	 	 
	 close, 2
	 
	 ymax = max(E)+20
	 ymin = min(E)-20
	 xmin = 0
	 xmax = max(orig) + 2
	 ylabel = ''
	 if (state.normalize eq 0) then ylabel = 'E [cm!E-1!N]'
	 plot, [0,0],[0,0],/nodata,xrange=[xmin,xmax],yrange=[ymin,ymax], ytitle=ylabel
	 for i = 0, n_levels-1 do begin
	 	  oplot, [orig(i)-0.5,orig(i)+0.5],[E(i),E(i)]
		  xyouts, orig(i)+0.7,E(i), strtrim(string(i),2)
	 endfor
end

pro Viewmodel_event, Event
	 
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
					 plotmodel, state
					 cierra_ps
				endif
		  end
		  
		  'NORM_ON': begin				
		   	plotmodel, set_norm_model(stash,1)
		  end
		  
		  'NORM_OFF': begin	
				plotmodel, set_norm_model(stash,0)
		  end
	 
	 endcase
end

function viewmodel_init
@vars.common
	 state = {baseWidget: 0L, drawWidget: 0L, normalize:0}
	 
;---------------------------------	
; Main window
;---------------------------------	
	 state.baseWidget = widget_base(TITLE='Atomic / molecular model', /ROW, $
	 	  MBAR=menuBar, RESOURCE_NAME='Viewmodel')
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
	 normFrame = widget_base(botones1, /COLUMN, /EXCLUSIVE)
	 normalButton = widget_button(normFrame, VALUE='Fictitial axis', $
                               UVALUE='NORM_ON')
	 notnormalButton = widget_button(normFrame, VALUE='Energy axis', $
                               UVALUE='NORM_OFF')	
	 widget_control, notnormalButton, /SET_BUTTON										 									 
	 
	 widget_control, widget_info(state.baseWidget, /CHILD), SET_UVALUE=state
	 return, state
end

pro viewmodel, LEADER=leader
	 if (not keyword_set(LEADER)) then leader=0
	 state = viewmodel_init()
	 
	 widget_control, state.baseWidget, /REALIZE, GROUP_LEADER=leader
	 plotmodel, state
	 
	 xmanager, 'Viewmodel', state.baseWidget, $
    	  EVENT_HANDLER='Viewmodel_event', GROUP_LEADER=leader
end
