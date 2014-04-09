; This procedure reads the output file from molALI
pro read_result, fichero
@vars.common

;---------------------------------
; Constants
;---------------------------------
	PH = 6.626176d-27
	PK = 1.380662d-16
	PC = 2.99792458d10

	PHK = PH / PK
	PHC = 2.d0 * PH / PC^2


	get_lun, u

	openr, u, fichero

	lb, u, 4
;---------------------------------
; Read active transitions
;---------------------------------
	readf, u, n_active

	lb, u, 1
;---------------------------------
; Read number of atomic/molecular levels
;---------------------------------
	readf, u, n_levels

	lb, u, 1
;---------------------------------
; Read number of grid points
;---------------------------------
	readf, u, n_radios

	lb, u, 1
;---------------------------------
; Read number of emergent mu
;---------------------------------
	readf, u, n_caract

	lb, u, 1
;---------------------------------
; Read number of frequency points for each line
;---------------------------------
	readf, u, n_frq

	 lb, u, 1
;---------------------------------
; Read if spectrum is present
;---------------------------------
	 readf, u, spectrum_present

	 lb, u, 1
;---------------------------------
; Read if background opacity has been used
;---------------------------------
	readf, u, back_opa_present

	 lb, u, 1
;---------------------------------
; Read if radiative rates are present
;---------------------------------
	readf, u, radrates_present

	lb, u, 1
	nombre_modelo = ''
;---------------------------------
; Read name of atomic/molecular model
;---------------------------------
	readf, u, nombre_modelo

;---------------------------------
; Geometry structure
;---------------------------------
	geometry = { Nradius: n_radios, Nmu: n_caract, r: dblarr(n_radios), mu: dblarr(n_caract) }

;---------------------------------
; Atomic structure
;---------------------------------
	atom = { Nlev: n_levels, pop: dblarr(n_levels,n_radios), b: dblarr(n_levels,n_radios), $
		levels: dblarr(n_levels), popi: dblarr(n_levels,n_radios), model:nombre_modelo, $
		BackOpacity: back_opa_present, backopa: dblarr(n_active,n_radios)}

;---------------------------------
; Transition structure
;---------------------------------
	transition = { Spectrum: spectrum_present, Nfrq: n_frq, Ntran: n_active, Lstar: dblarr(n_active,n_radios), $
	 	  Jbar: dblarr(n_active,n_radios), Jbar20: dblarr(n_active,n_radios), $
		  tau: dblarr(n_active,n_radios), S: dblarr(n_active,n_radios), B: dblarr(n_active,n_radios),$
		  Texc: dblarr(n_active,n_radios), opacity: dblarr(n_active,n_radios), $
		  low: intarr(n_active), up: intarr(n_active), freq: dblarr(n_active), $
		  Rul: dblarr(n_active,n_radios), Rlu: dblarr(n_active,n_radios), radrates: radrates_present }

;---------------------------------
; Spectrum structure
;---------------------------------
	 spectrum = { Nlambda: transition.Nfrq*transition.Ntran, $
	 	  lambda: dblarr(transition.Ntran,transition.Nfrq), $
	 	  I: dblarr(geometry.nmu,transition.Ntran,transition.Nfrq), $
		  ILTE: dblarr(geometry.nmu,transition.Ntran,transition.Nfrq), $
		  flux: dblarr(transition.Ntran,transition.Nfrq), $
		  fluxLTE: dblarr(transition.Ntran,transition.Nfrq) }

	lb, u, 1
	levels = intarr(2,n_active)
	frec = dblarr(n_active)
	temp = dblarr(3)
;---------------------------------
; Read transitions
;---------------------------------
	for i = 0, n_active-1 do begin
		readf, u, temp
		levels(0,i) = fix(temp(0))
		levels(1,i) = fix(temp(1))
		frec(i) = temp(2)
	endfor

	transition.low = levels(1,*)
	transition.up = levels(0,*)
	transition.freq = frec

	lb, u, 3
;---------------------------------
; Read shell grid
;---------------------------------
	temp = dblarr(n_radios)
	readf, u, temp
	geometry.r = temp

	lb, u, 3
;---------------------------------
; Read angular grid
;---------------------------------
	temp = dblarr(n_caract)
	readf, u, temp
	geometry.mu = temp

	lb, u, 3
;---------------------------------
; Read energy levels
;---------------------------------
	temp = dblarr(n_levels)
	readf, u, temp
	atom.levels = temp

	lb, u, 3
;---------------------------------
; Read final populations
;---------------------------------
	temp = dblarr(n_levels+1,n_radios)
	readf, u, temp
	atom.pop = temp(1:n_levels,*)

	lb, u, 6
;---------------------------------
; Read departure coefficients
;---------------------------------
	readf, u, temp
	atom.b = temp(1:n_levels,*)

	lb, u, 6
;---------------------------------
; Read initial population
;---------------------------------
	readf, u, temp
	atom.popi = temp(1:n_levels,*)

	lb, u, 6
;---------------------------------
; Read lstar
;---------------------------------
	temp = dblarr(n_active+1,n_radios)
	readf, u, temp
	transition.Lstar = temp(1:n_active,*)

	lb, u, 6
;---------------------------------
; Read Jbar
;---------------------------------
	readf, u, temp
	transition.Jbar = temp(1:n_active,*)

	lb, u, 6
;---------------------------------
; Read Jbar20
;---------------------------------
;	lb, u, n_radios
	readf, u, temp
	transition.Jbar20 = temp(1:n_active,*)


	lb, u, 14
;---------------------------------
; Read source functions
;---------------------------------
	temp = dblarr(5,n_radios)
	for i = 1, n_active do begin
		readf, u, temp

		transition.tau(i-1,*) = temp(0,*)
		transition.S(i-1,*) = temp(1,*)
		transition.B(i-1,*) = temp(2,*)
		transition.Texc(i-1,*) = temp(3,*)
		transition.opacity(i-1,*) = temp(4,*)

		if (i ne n_active) then lb, u, 9
	endfor


	close, u
	free_lun, u

end
