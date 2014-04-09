pro read_flux, fichero
@vars.common


; If an emergent spectrum is present
 if (transition.Spectrum eq 1) then begin
	 get_lun, u
	 openr, u, fichero
	 lb, u, 1
;---------------------------------
; Read NLTE flux
;---------------------------------
	 temp = dblarr(2,spectrum.Nlambda)
	 readf, u, temp
	 for j = 1, transition.Ntran do begin
	 	  start = (j-1) * transition.Nfrq
		  ending = j * transition.Nfrq
		  spectrum.lambda(j-1,*) = temp(0,start:ending-1)
	 	  spectrum.flux(j-1, *) = temp(1,start:ending-1)

; Little trick to avoid extremal strange values of the flux
;		  media = mean(spectrum.flux(j-1,*))
;		  std = stddev(spectrum.flux(j-1,*))
;		  ind_no = where(abs(spectrum.flux(j-1,*)-media) gt 3.*std)
;		  ind_yes = where(abs(spectrum.flux(j-1,*)-media) le 3.*std)
;		  if (ind_no(0) ne -1) then begin
;		   	spectrum.flux(j-1,ind_no) = mean(spectrum.flux(j-1,5:10))
;		  endif
	 endfor

;	 lb, u, 1
;---------------------------------
; Read LTE flux
;---------------------------------
	 temp = dblarr(2,spectrum.Nlambda)
	 readf, u, temp
	 for j = 1, transition.Ntran do begin
	 	  start = (j-1) * transition.Nfrq
	 	  ending = j * transition.Nfrq
	 	  spectrum.fluxLTE(j-1, *) = temp(1,start:ending-1)

; Little trick to avoid extremal strange values of the flux
		  media = mean(spectrum.fluxLTE(j-1,*))
		  std = stddev(spectrum.fluxLTE(j-1,*))
		  ind_no = where(abs(spectrum.fluxLTE(j-1,*)-media) gt 3.*std)
		  ind_yes = where(abs(spectrum.fluxLTE(j-1,*)-media) le 3.*std)
		  if (ind_no(0) ne -1) then begin
	 		  spectrum.fluxLTE(j-1,ind_no) = mean(spectrum.fluxLTE(j-1,5:10))
		  endif
	 endfor

	 close, u
	 free_lun, u

 endif
end
