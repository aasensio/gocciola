; This function reads n lines of rubbish
pro lb, unit, n
	a = ''
	for i = 1, n do begin
		readf, unit, a
	endfor
end 

; This procedure reads the output file of jac_molec.f90
pro leer, fichero

; Constants	
	PH = 6.626176d-27
	PK = 1.380662d-16
	PC = 2.99792458d10
	
	PHK = PH / PK
	PHC = 2.d0 * PH / PC^2
	
	u = 2
	openr, u, fichero

	lb, u, 4
; Read active transitions	
	readf, u, n_active
	
	lb, u, 1
; Read number of atomic/molecular levels
	readf, u, n_levels
	
	lb, u, 1
; Read number of grid points
	readf, u, n_radios
	
	lb, u, 1
	levels = intarr(2,n_active)
	frec = dblarr(n_active)
	temp = dblarr(3)
; Read transitions
	for i = 0, n_active-1 do begin
		readf, u, temp
		levels(0,i) = fix(temp(0))
		levels(1,i) = fix(temp(1))
		frec(i) = temp(2)
	endfor
	
	lb, u, 3
; Read shell grid
	r = dblarr(n_radios)
	readf, u, r
	
	lb, u, 3
; Read final populations
	pop = dblarr(n_levels+1,n_radios)
	readf, u, pop
	
	lb, u, 6
; Read departure coefficients
	depart = dblarr(n_levels+1,n_radios)
	readf, u, depart
	
	lb, u, 6
; Read lstar
	lstar = dblarr(n_active+1,n_radios)
	readf, u, lstar
	
	lb, u, 6
; Read Jbar
	Jbar = dblarr(n_active+1,n_radios)
	readf, u, Jbar
	
	lb, u, 11
; Read source functions
	tau = dblarr(n_active,n_radios)
	S = dblarr(n_active,n_radios)
	B = dblarr(n_active,n_radios)
	temp = dblarr(4,n_radios)
	for i = 1, n_active do begin
		readf, u, temp
		
		tau(i-1,*) = temp(0,*)
		S(i-1,*) = temp(1,*)
		B(i-1,*) = temp(2,*)
		
		if (i ne n_active) then lb, u, 6
	endfor

; Equivalent temperatures	
	T = dblarr(n_active,n_radios)
	for i = 0, n_active-1 do begin
		T(i,*) = PHK * frec(i) / alog(1.d0 + (PHC * frec(i)^3 / (S(i,*) * B(i,*))))
	endfor
	
	close,2

	stop
end
