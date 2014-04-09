pro analyze_event, event
	 widget_control, Event.id, GET_UVALUE=Action
	 
	 stash = widget_info( Event.handler, /CHILD)
	 widget_control, stash, GET_UVALUE=state
	 	 
	 case Action of
    	  'QUIT':  begin      
         	widget_control, Event.top, /DESTROY
    	  end
		  
		  'ATMOS': begin
		   	viewatmos, LEADER=event.top
		  end
		  
		  'MODEL': begin
		   	viewmodel, LEADER=event.top
		  end
		  
		  'BACKOPA': begin
		   	viewback, LEADER=event.top
		  end
		  
		  'POPULATION': begin
		   	viewpop, LEADER=event.top
		  end		  
		  
		  'SOURCEFUN': begin
		   	viewsource, LEADER=event.top
		  end		  		  
		  
		  'SPECTRUM': begin
		   	viewspectrum, LEADER=event.top
		  end
		  
		  'FLUX': begin
		   	viewflux, LEADER=event.top
		  end
		  		  		  
		  'RADRATES': begin
		   	viewradrates, LEADER=event.top
		  end
		  
		  'ZEEMAN': begin
		   	viewzeem, LEADER=event.top
		  end
		  
		  'ITERATION': begin
		   	viewiter, LEADER=event.top
		  end
		  
		  'LIFETIMES': begin
		   	viewlifetime, LEADER=event.top
		  end
	 endcase
	 


end


function analyze_init
	 state = {baseWidget: 0L, logoWidget: 0L, directory: '', specButton: 0L, fluxButton: 0L, $
	 	  backButton: 0L, radratesButton: 0L, zeemanButton: 0L}
	 state.baseWidget = widget_base(TITLE='Analyzer', MBAR=menuBar)
	 
	 analyzeBase = widget_base(state.baseWidget, /COLUMN)
	 state.logoWidget = widget_draw(analyzeBase, XSIZE=450, YSIZE=100, /FRAME)	 	 	 
	  
	 fileMenu = widget_button(menuBar, VALUE='File', /MENU)
	 quitButton = widget_button(fileMenu, VALUE='Quit', UVALUE='QUIT')
	 
	 atmosMenu = widget_button(menuBar, VALUE='Models', /MENU)
	 sgridButton = widget_button(atmosMenu, VALUE='Atmosphere', UVALUE='ATMOS')
	 atomicButton = widget_button(atmosMenu, VALUE='Atomic/molecular', UVALUE='MODEL')
	 state.backButton = widget_button(atmosMenu, VALUE='Background opacity', UVALUE='BACKOPA')
	 
	 levelMenu = widget_button(menuBar, VALUE='Levels', /MENU)
	 populButton = widget_button(levelMenu, VALUE='Population', UVALUE='POPULATION')	 
	 lifetButton = widget_button(levelMenu, VALUE='Lifetimes', UVALUE='LIFETIMES')
	 
	 transMenu = widget_button(menuBar, VALUE='Transitions', /MENU)
	 sourceButton = widget_button(transMenu, VALUE='Source function', UVALUE='SOURCEFUN')
	 state.specButton = widget_button(transMenu, VALUE='Spectrum', UVALUE='SPECTRUM')
	 state.fluxButton = widget_button(transMenu, VALUE='Flux', UVALUE='FLUX')
	 state.radratesButton = widget_button(transMenu, VALUE='Radiative rates', UVALUE='RADRATES')
	 state.zeemanButton = widget_button(transMenu, VALUE='Zeeman synthesis', UVALUE='ZEEMAN')
	 
	 iterMenu = widget_button(menuBar, VALUE='Iterations', /MENU)
	 iterButton = widget_button(iterMenu, VALUE='Iterations', UVALUE='ITERATION')
	 
	 widget_control, widget_info(state.baseWidget, /CHILD), SET_UVALUE=state
	 
	 return, state
end

; STOP : stops after reading all data
; OUT : puts in the variable all the data

pro analyze, RESULTSFILE=resultsfile, ITERATIONFILE=iterationfile, SPECTRUMFILE=spectrumfile, $
	 ATMOSFILE=atmosfile, FLUXFILE=fluxfile, BACKFILE=backfile, RADRATESFILE=radratesfile, $
	 ZEEMANFILE=zeemanfile, STOP=stop, OUT=out
@vars.common
@files.common

	 !p.multi = 0
	 !p.charsize = 1.0
;	 !p.region  = [0,0,1,1]
;	 !p.noerase = 0

	 state = analyze_init()
;---------------------------------
; Read results file
;---------------------------------
	 if (not keyword_set(RESULTSFILE)) then begin
	 	  read_result, '../Results/output.dat'
	 endif else begin
	 	  read_result, resultsfile
	 endelse

	 if (transition.Spectrum eq 0) then begin
	 	 widget_control, state.specButton, sensitive=0
		 widget_control, state.fluxButton, sensitive=0
	 endif

;---------------------------------
; Read background opacity file
;---------------------------------
	 if (not keyword_set(BACKFILE)) then begin
	 	  read_backopa, '../Results/background.dat'
	 endif else begin
	 	  read_backopa, backfile
	 endelse
	 
	 if (atom.BackOpacity eq 0) then begin
	 	  widget_control, state.backButton, sensitive = 0
	 endif
	 	 
;---------------------------------
; Read iteration file
;---------------------------------
	 if (not keyword_set(ITERATIONFILE)) then begin
	 	  read_iter, '../Results/iteration.dat'
	 endif else begin
	 	  read_iter, iterationfile
	 endelse
	 
;---------------------------------
; Read radiative rates file
;---------------------------------
	 if (not keyword_set(RADRATESFILE)) then begin
	 	  read_radrates, '../Results/radrates.dat'
	 endif else begin
	 	  read_radrates, radratesfile
	 endelse	 
	 
	 if (transition.radrates eq 0) then begin
	 	  widget_control, state.radratesButton, sensitive = 0
	 endif
	 
;---------------------------------
; Read spectrum file
;---------------------------------

	 if (not keyword_set(SPECTRUMFILE)) then begin
	 	  spectrfile = '../Results/spectrum.dat'	
	 endif else begin
	 	  spectrfile = spectrumfile
	 endelse
	 	 
;---------------------------------
; Read flux file
;---------------------------------
;	 if (not keyword_set(NOSPECTRUM)) then begin
		 if (not keyword_set(FLUXFILE)) then begin
	 		  read_flux, '../Results/flux.dat'
		 endif else begin
	 		  read_flux, fluxfile
		 endelse
;	 endif	 

;---------------------------------
; Read atmosphere file
;---------------------------------
	 if (not keyword_set(ATMOSFILE)) then begin
	 	  read_atmos, '../Results/atmosphere.dat'
	 endif else begin
	 	  read_atmos, atmosfile
	 endelse	 
	 
;---------------------------------
; Read Zeeman synthesis
;---------------------------------
	 salid = 0
	 if (not keyword_set(ZEEMANFILE)) then begin
;	 	  read_zeeman, '../Results/zeeman.dat', salid
	 endif else begin
;	 	  read_zeeman, radratesfile, salid
	 endelse	 
	 
	 if (salid eq 0) then begin
	 	  widget_control, state.radratesButton, sensitive = 0
	 endif
	 
	 
	 if (transition.radrates eq 0) then begin
	 	  widget_control, state.radratesButton, sensitive = 0
	 endif	  
	 
	 if (keyword_set(STOP)) then stop
	 
	 if (arg_present(OUT)) then begin
	 	  out = {atom:atom, transition:transition, geometry:geometry, spectrum:spectrum,$
		   	iterat:iterat, atmos:atmos}
		  out.atom = atom
		  out.transition = transition
		  out.geometry = geometry
		  out.spectrum = spectrum
		  out.iterat = iterat
		  out.atmos = atmos	 
	 endif else begin
	 
		 widget_control, state.baseWidget, /REALIZE

;		 read_gif, 'logo.gif', logo, r, g, b
	;	 tvlct, r, g, b
;		 widget_control, state.logoWidget, GET_VALUE=WindowNo
;   	 wset, WindowNo
;		 tv, logo 

		 xmanager, 'Analyze', state.baseWidget, EVENT_HANDLER='Analyze_Event'
	 endelse
end
