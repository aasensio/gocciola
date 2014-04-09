pro read_one_spectrum, fichero, tran, mu
@vars.common
;	 if (transition.Spectrum eq 1) then begin
	 	  get_lun, u
		  openr, u, fichero
		  a=''
		  readf, u, a
		  tam = strlen(a) + 2
		  point_lun, u, 0

; NLTE line
		  offset = mu * spectrum.Nlambda + tran * transition.Nfrq
		  point_lun, u, long(offset * tam * 1L)
		  temp = dblarr(2,transition.Nfrq)

		  readf, u, temp

		  spectrum.lambda(tran,*) = temp(0,*)

	 	  spectrum.I(mu,tran,*) = temp(1,*)

; LTE line
		  orig = ((geometry.Nmu-1) * spectrum.Nlambda + transition.Ntran * transition.Nfrq ) * tam
		  offset = mu * spectrum.Nlambda + tran * transition.Nfrq
		  point_lun, u, long(orig + offset * tam * 1L)
		  temp = dblarr(2,transition.Nfrq)

		  readf, u, temp

	 	  spectrum.ILTE(mu,tran,*) = temp(1,*)

		  close, u
		  free_lun, u
;	 endif
end

; Reads a line
; It returns the spectrum of the tran line at angle mu in the structure estruc
pro read_line, fichero, tran, mu, estruc
@vars.common
;	 if (transition.Spectrum eq 1) then begin
	 	  get_lun, u
		  openr, u, fichero
		  a=''
		  ;readf, u, a
		  tam = strlen(a) + 1

		  point_lun, u, 0

; NLTE line
		  offset = mu * spectrum.Nlambda + tran * transition.Nfrq
		  point_lun, u, long(offset * tam * 1L)
		  temp = dblarr(2,transition.Nfrq)

		  readf, u, temp

		  estruc.spectrum.lambda(tran,*) = temp(0,*)

	 	  estruc.spectrum.I(mu,tran,*) = temp(1,*)

		  close, u
		  free_lun, u
;	 endif
end

pro read_spectrum, fichero
@vars.common


; If an emergent spectrum is present
; if (transition.Spectrum eq 1) then begin
	 get_lun, u
	 openr, u, fichero

;	 lb, u, 1
;---------------------------------
; Read NLTE spectrum
;---------------------------------
	 temp = dblarr(2,geometry.Nmu*spectrum.Nlambda)
	 readf, u, temp
	 for i = 1, geometry.Nmu do begin
		 for j = 1, transition.Ntran do begin
		   	start = (j-1) * transition.Nfrq + (i-1) * spectrum.Nlambda
				ending = j * transition.Nfrq + (i-1) * spectrum.Nlambda
				spectrum.lambda(j-1,*) = temp(0,start:ending-1)
		   	spectrum.I(i-1, j-1, *) = temp(1,start:ending-1)
		 endfor
	 endfor

;	 lb, u, 1
;---------------------------------
; Read LTE spectrum
;---------------------------------
	 temp = dblarr(2,geometry.Nmu*spectrum.Nlambda)
	 readf, u, temp
	 for i = 1, geometry.Nmu do begin
		 for j = 1, transition.Ntran do begin
		   	start = (j-1) * transition.Nfrq + (i-1) * spectrum.Nlambda
				ending = j * transition.Nfrq + (i-1) * spectrum.Nlambda
		   	spectrum.ILTE(i-1, j-1, *) = temp(1,start:ending-1)
		 endfor
	 endfor

	 close, u
	 free_lun, u

; endif
end
