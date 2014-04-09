function set_rad_atmos, stash, onoff
	 widget_control, stash, GET_UVALUE=state

 	 state.radial = onoff
	 widget_control, stash, SET_UVALUE=state
	 return, state	 
end

function set_logx_atmos, stash, onoff
	 widget_control, stash, GET_UVALUE=state
 	 state.logx = onoff
	 widget_control, stash, SET_UVALUE=state	
	 return, state
end

function set_logy_atmos, stash, onoff
	 widget_control, stash, GET_UVALUE=state
 	 state.logy = onoff
	 widget_control, stash, SET_UVALUE=state	
	 return, state
end

pro plot_atmos, state
@vars.common
	 
	 if (state.radial eq 1) then begin
	 	 !p.multi = [0,3,2]
		 plot, atmos.r, atmos.T, xtitle='Radius [cm]', ytitle='T [K]'
		 plot, atmos.r, atmos.NH2, xtitle='Radius [cm]', ytitle='N(H!I2!N) [cm!E-3!N]'
		 plot, atmos.r, atmos.N, xtitle='Radius [cm]', ytitle='N / N(H!I2!N) [cm!E-3!N]'
		 plot, atmos.r, atmos.N_e, xtitle='Radius [cm]', ytitle='N!Ie!N [cm!E-3!N]'
		 plot, atmos.r, atmos.T_dust, xtitle='Radius [cm]', ytitle='T!Idust!N [K]'
		 plot, atmos.r, atmos.V_mac, xtitle='Radius [cm]', ytitle='V!Imac!N [km/s]'
	 endif
	 
	 if (state.radial eq 0) then begin
	 	 !p.multi = [0,2,1]
		 plot, atmos.p, xtitle='Number', ytitle = 'Characteristics [cm]',XLOG=state.logx, YLOG=state.logy
		 plot, atmos.mu_exit, xtitle='Number', ytitle = '!7l!3!Iexit!N',XLOG=state.logx, YLOG=state.logy
	 endif
	 
	 if (state.radial eq 2) then begin
	 	 !p.multi = 0
		 theta = dindgen(200)/200. * 2.d0 * !DPI
		 plot, [0.0,0.0],[0.0,0.0],/nodata,xrange=[0,max(atmos.r)], $ 
		   	yrange=[0,max(atmos.r)], xstyle=1, ystyle=1,XLOG=state.logx, YLOG=state.logy
		 for i = 0, atmos.N_radius-1 do begin
				oplot, atmos.r(i)*cos(theta), atmos.r(i)*sin(theta)				
		 endfor
	 endif
	 !p.multi=0
	 if (state.radial eq 3) then begin
		  plot, atmos.r, atmos.T, xtitle='Radius [cm]', ytitle='T [K]',XLOG=state.logx, YLOG=state.logy
	 endif
	 
	 if (state.radial eq 4) then begin
	 	  plot, atmos.r, atmos.NH2, xtitle='Radius [cm]', ytitle='N(H!I2!N) [cm!E-3!N]',XLOG=state.logx, YLOG=state.logy
	 endif
	 
	 if (state.radial eq 5) then begin
	 	  plot, atmos.r, atmos.N, xtitle='Radius [cm]', ytitle='N / N(H!I2!N) [cm!E-3!N]',XLOG=state.logx, YLOG=state.logy
	 endif
	 
	 if (state.radial eq 6) then begin
	 	  plot, atmos.r, atmos.N_e, xtitle='Radius [cm]', ytitle='N!Ie!N [cm!E-3!N]',XLOG=state.logx, YLOG=state.logy
	 endif
	 
	 if (state.radial eq 7) then begin
	 	  plot, atmos.r, atmos.T_dust, xtitle='Radius [cm]', ytitle='T!Idust!N [K]',XLOG=state.logx, YLOG=state.logy
	 endif
	 
	 if (state.radial eq 8) then begin
	 	  plot, atmos.r, atmos.v_mac, xtitle='Radius [cm]', ytitle='V!Imac!N [km/s]',XLOG=state.logx, YLOG=state.logy
	 endif
	 	 	 
	 !p.multi = 0
end

pro Viewatmos_event, Event
	 
	 stash = widget_info( Event.handler, /CHILD)
	 widget_control, stash, GET_UVALUE=state

	 widget_control, Event.id, GET_UVALUE=Action
	 	 
	 case Action of
	 	  'QUIT': begin
		   	widget_control, Event.top, /DESTROY
		  end
		  
		  'LOGX_ON': begin
				plot_atmos, set_logx_atmos(stash,1)
		  end
		  
		  'LOGX_OFF': begin
				plot_atmos, set_logx_atmos(stash,0)
		  end
		  
		  'LOGY_ON': begin
				plot_atmos, set_logy_atmos(stash, 1)
		  end
		  
		  'LOGY_OFF': begin
				plot_atmos, set_logy_atmos(stash, 0)
		  end
		  
		  'POSTCRIPT': begin	
		   	name = dialog_pickfile()
				if (name ne '') then begin
					 abre_ps, name, /landscape
					 plot_atmos, state
					 cierra_ps
				endif
		  end
		  
		  'ALL': begin			   	
				plot_atmos, set_rad_atmos(stash,1)
		  end

		  'CHARACT': begin			   	
				plot_atmos, set_rad_atmos(stash,0)
		  end		  
		  
		  'SCHEME': begin
		   	plot_atmos, set_rad_atmos(stash,2)
		  end
		  
		  'TEMP' : begin
		   	plot_atmos, set_rad_atmos(stash,3)
		  end
		  
		  'MOLHYD' : begin
		   	plot_atmos, set_rad_atmos(stash,4)
		  end
		  
		  'ABUNDANCE' : begin
		   	plot_atmos, set_rad_atmos(stash,5)
		  end
		  
		  'ELECTRON' : begin
		   	plot_atmos, set_rad_atmos(stash,6)
		  end
		  
		  'DUST' : begin
		   	plot_atmos, set_rad_atmos(stash,7)
		  end
		  
		  'MACVEL' : begin
		   	plot_atmos, set_rad_atmos(stash,8)
		  end		  		  		  		  		  
		  		  	 
	 endcase
end

function viewatmos_init
@vars.common
	 state = {baseWidget: 0L, drawWidget: 0L, lineList: 0L, radial: 2, logx: 0, logy: 0}
	 
;---------------------------------	
; Main window
;---------------------------------	
	 state.baseWidget = widget_base(TITLE='Atmosphere model', /ROW, $
	 	  MBAR=menuBar, RESOURCE_NAME='Viewatmos')
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
; Radial variables button
;---------------------------------	

	 radialMenu = widget_button(fileMenu, VALUE='Radial variables',/menu)									  
	 allButton = widget_button(radialMenu, VALUE='All radial variables', UVALUE='ALL', $
                             RESOURCE_NAME='all')									  	 
	 TButton = widget_button(radialMenu, VALUE='Temperature', UVALUE='TEMP', $
		   		 	   	RESOURCE_NAME='temperature')		  
	 molHButton = widget_button(radialMenu, VALUE='Molecular hydrogen abundance', UVALUE='MOLHYD', $
		   		 	   	RESOURCE_NAME='molhyd')		  								
	 nButton = widget_button(radialMenu, VALUE='Abundance', UVALUE='ABUNDANCE', $
		   		 	   	RESOURCE_NAME='abundance')		  								
	 neButton = widget_button(radialMenu, VALUE='Electron density', UVALUE='ELECTRON', $
		   		 	   	RESOURCE_NAME='electron')		  								
	 TdustButton = widget_button(radialMenu, VALUE='Dust temperature', UVALUE='DUST', $
		   		 	   	RESOURCE_NAME='dust')		  								
	 vButton = widget_button(radialMenu, VALUE='Macroscopic velocity', UVALUE='MACVEL', $
		   		 	   	RESOURCE_NAME='macvel')		  																								
									  
;---------------------------------	
; Characteristics button
;---------------------------------	
	 characButton = widget_button(fileMenu, VALUE='Characteristics', UVALUE='CHARACT', $
                             RESOURCE_NAME='charact')									  									  									  

;---------------------------------	
; Drawing button
;---------------------------------	
	 drawButton = widget_button(fileMenu, VALUE='Scheme', UVALUE='SCHEME', $
                             RESOURCE_NAME='scheme')									  									  									  
									  	 
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
	 
											
	 widget_control, widget_info(state.baseWidget, /CHILD), SET_UVALUE=state
	 
	 return, state
end


pro viewatmos, LEADER=leader
	 if (NOT keyword_set(LEADER)) then leader=0
	 
	 state = viewatmos_init()
	 
	 widget_control, state.baseWidget, /REALIZE, GROUP_LEADER=leader
	 
	 plot_atmos, state
	 
	 xmanager, 'Viewatmos', state.baseWidget, $
    	  EVENT_HANDLER='Viewatmos_event', GROUP_LEADER=leader	 
end
