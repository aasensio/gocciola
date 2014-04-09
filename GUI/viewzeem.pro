;---------------------------------
; This function sets the lower range in x axis
;---------------------------------
function set_spectrum_stokes, stash, onoff
    widget_control, stash, GET_UVALUE=state
	 
	 state.stokes = onoff	 
    widget_control, stash, SET_UVALUE=state
	 return, state
end

;---------------------------------
; This function sets the y-axis type
;---------------------------------
function set_spectrum_yaxis, stash, onoff
    widget_control, stash, GET_UVALUE=state

    state.inten_temp = onoff
    widget_control, stash, SET_UVALUE=state
    return, state   
end


;---------------------------------
; This function sets the flag to use a wavelength axis or a frequency axis
;---------------------------------
function set_freqwave, stash, onoff
	 widget_control, stash, GET_UVALUE=state

 	 state.freqwave = onoff
	 widget_control, stash, SET_UVALUE=state
	 return, state	 
end

;---------------------------------
; This function sets the current line being plotted
;---------------------------------
function set_line, stash, onoff
	 widget_control, stash, GET_UVALUE=state

 	 state.current = onoff
	 widget_control, stash, SET_UVALUE=state
	 return, state	 
end

;---------------------------------
; This function plots a line
;---------------------------------
pro plot_zeeman, state
@vars.common
	 
	 !p.multi = 0
	 
	 if (state.stokes eq 4) then begin
	 	  !p.multi = [0,2,2]		  
	 endif
	 
	 lamb = zeeman.lambda(state.current,*)
	 case state.freqwave of
; Wavelength		  
		  0 : begin
					 titlex = 'Wavelength (!3Å)'
		   	end
; Frequency	 	  
		  1 : begin
		   		 lamb = 2.99792458d10 / (lamb * 1.d-8)
					 titlex = 'Frequency (Hz)'		   		 
		   	end
; Velocity		  
		  2 : begin
		   		 lamb = 2.99792458d5 * (1.d0 - zeeman.lambda0(state.current) / lamb)
		   		 titlex = 'Velocity (km/s)'
		   	end
; Wavenumber		  
		  3 : begin
		   		 lamb = 1.d0 / (lamb * 1.d-8)
		   		 titlex = 'Wavenumber (cm!E-1!N)'
		   	end
	 endcase
	 
	 case state.stokes of
	 	  0 : begin					 
					 S = zeeman.SI(state.current,*) 
	 	   		 plot, lamb, S, xtitle=titlex, ytitle='I'
				end
				
	 	  1 : begin
					 S = zeeman.SQ(state.current,*) / zeeman.SI(state.current,*)
	 	   		 plot, lamb, S, xtitle=titlex, ytitle='Q/I'		  
				end
				
	 	  2 : begin
					 S = zeeman.SU(state.current,*) / zeeman.SI(state.current,*)
	 	   		 plot, lamb, S, xtitle=titlex, ytitle='U/I'		  
				end
				
	 	  3 : begin
					 S = zeeman.SV(state.current,*) / zeeman.SI(state.current,*)
	 	   		 plot, lamb, S, xtitle=titlex, ytitle='V/I'		  
	 	   	end
				
		  4 : begin
					 u = zeeman.SI(state.current,*) 
	 	   		 plot, lamb, u, xtitle=titlex, ytitle='I'
					 
					 S = zeeman.SQ(state.current,*) / zeeman.SI(state.current,*)
	 	   		 plot, lamb, S, xtitle=titlex, ytitle='Q/I'
					 
					 S = zeeman.SU(state.current,*) / zeeman.SI(state.current,*)
	 	   		 plot, lamb, S, xtitle=titlex, ytitle='U/I'
					 
					 S = zeeman.SV(state.current,*) / zeeman.SI(state.current,*)
	 	   		 plot, lamb, S, xtitle=titlex, ytitle='V/I'		  
		   	end
	 endcase	 	  
	 		 
end


;---------------------------------	
; This function fills the line list
;---------------------------------	
pro fill_lines_zeeman, state
@vars.common

	 lines = strarr(zeeman.n_lines)
	 lines = string(zeeman.lambda0)
	 widget_control, state.lineList, SET_VALUE=lines
end

pro Viewzeem_event, Event
	 
	 stash = widget_info( Event.handler, /CHILD)
	 widget_control, stash, GET_UVALUE=state

	 widget_control, Event.id, GET_UVALUE=Action
	 	 
	 case Action of
	 	  'QUIT': begin
		   	!x.range = 0
				!y.range = 0
				widget_control, Event.top, /DESTROY
		  end
		  
		  'POSTCRIPT': begin	
		   	name = dialog_pickfile()
				if (name ne '') then begin
					 abre_ps, name, /landscape					 
	 				 plot_zeeman, state
		   		 cierra_ps
				endif
		  end
		  
		  'LEVELLIST': begin
		   	plot_zeeman, set_line(stash, Event.index)
		  end
		  		  		  
		  'FREQ_ON': begin			   
				plot_zeeman, set_freqwave(stash,1)
		  end

		  'FREQ_OFF': begin
				plot_zeeman, set_freqwave(stash,0)
		  end		  
		  
		  'FREQ_VEL': begin
				plot_zeeman, set_freqwave(stash,2)
		  end	
		  
		  'FREQ_NUMBER': begin		   					
				plot_zeeman, set_freqwave(stash,3)
		  end		  
		  		  		  
		  'INTENSITY': begin
				  plot_zeeman, set_spectrum_yaxis(stash, 0)
        end
		  
		  'TEMPERATURE': begin	 
              plot_zeeman, set_spectrum_yaxis(stash, 1)
        end	  
		  
		  'JANSKY': begin	 
              plot_zeeman, set_spectrum_yaxis(stash, 2)
        end	  
		  
		  'WCM2': begin	 
              plot_zeeman, set_spectrum_yaxis(stash, 3)
        end	  
		  
		  'I_ON': begin	 
              plot_zeeman, set_spectrum_stokes(stash, 0)
        end	  
		  
		  'Q_ON': begin	 
              plot_zeeman, set_spectrum_stokes(stash, 1)
        end	  
		  
		  'U_ON': begin	 
              plot_zeeman, set_spectrum_stokes(stash, 2)
        end	  
		  
		  'V_ON': begin	 
              plot_zeeman, set_spectrum_stokes(stash, 3)
        end	 
		   
		  'IQUV_ON': begin	 
              plot_zeeman, set_spectrum_stokes(stash, 4)
        end	  		  		  		  
		  	 
	 endcase
end

function viewzeem_init
@vars.common
	 state = {baseWidget: 0L, drawWidget: 0L, labelMu: 0L, lineList: 0L, stokes: 4L,$
	 	  freqwave: 0L, inten_temp: 0L, current: 0L}
	 
;---------------------------------	
; Main window
;---------------------------------	
	 state.baseWidget = widget_base(TITLE='Zeeman', /ROW, $
	 	  MBAR=menuBar, RESOURCE_NAME='Viewzeeman')
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
; I, Q, U, V, or IQUV
;---------------------------------	
	 stokesFrame = widget_base(botones1, /COLUMN, /EXCLUSIVE)
	 IButton = widget_button(stokesFrame, VALUE='I', $
                               UVALUE='I_ON')
	 QButton = widget_button(stokesFrame, VALUE='Q', $
                               UVALUE='Q_ON')
	 UButton = widget_button(stokesFrame, VALUE='U', $
                               UVALUE='U_ON')
	 VButton = widget_button(stokesFrame, VALUE='V', $
                               UVALUE='V_ON')
	 IQUVButton = widget_button(stokesFrame, VALUE='IQUV', $
                               UVALUE='IQUV_ON')										 										 										 										 
	 widget_control, IQUVButton, /SET_BUTTON		
	 
;---------------------------------	
; Frequency/wavelength
;---------------------------------	
	 freqwaveFrame = widget_base(botones1, /COLUMN, /EXCLUSIVE)
	 freqButton = widget_button(freqwaveFrame, VALUE='Frequency axis', $
                               UVALUE='FREQ_ON')
	 waveButton = widget_button(freqwaveFrame, VALUE='Wavelength axis', $
                               UVALUE='FREQ_OFF')										 										 
	 velocButton = widget_button(freqwaveFrame, VALUE='Velocity axis', $
                               UVALUE='FREQ_VEL')										 										 										 
	 wavnumButton = widget_button(freqwaveFrame, VALUE='Wavenumber axis', $
                               UVALUE='FREQ_NUMBER')										 										 										 										 
	 widget_control, waveButton, /SET_BUTTON	
	 
;---------------------------------	
; Intensity/Temperature/mJ
;---------------------------------	 
	 intenFrame = widget_base(botones1, /COLUMN, /EXCLUSIVE)										 
	 intenButton = widget_button(intenFrame, VALUE='Flux (cgs)', $
                               UVALUE='INTENSITY')										 										 										 
	 temperButton = widget_button(intenFrame, VALUE='Temperature', $
                               UVALUE='TEMPERATURE')
	 janskyButton = widget_button(intenFrame, VALUE='Flux (Jy)', $
                               UVALUE='JANSKY')										 
	 wcm2Button = widget_button(intenFrame, VALUE='Flux (W*cm^-2)', $
                               UVALUE='WCM2')										 										 
										 
	 widget_control, (state.inten_temp) ? temperButton : intenButton, /SET_BUTTON	 		 
										 										 
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
	 fill_lines_zeeman, state
	 	 	 
	 widget_control, widget_info(state.baseWidget, /CHILD), SET_UVALUE=state
	 
	 return, state
end


pro viewzeem, LEADER=leader
	 if (NOT keyword_set(LEADER)) then leader=0
	 
	 state = viewzeem_init()
	 
	 widget_control, state.baseWidget, /REALIZE, GROUP_LEADER=leader
	 
	 xmanager, 'Viewzeeman', state.baseWidget, $
    	  EVENT_HANDLER='Viewzeem_event', GROUP_LEADER=leader	 
end
