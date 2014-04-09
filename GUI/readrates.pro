pro read_radrates, fichero
@vars.common

	 
; If an emergent spectrum is present
 if (transition.radrates eq 1) then begin
	 get_lun, u
	 openr, u, fichero
	
	 temp = dblarr(2,geometry.Nradius)
	 for i = 1, transition.Ntran do begin
	 	  readf, u, temp
		  transition.Rul(i-1,*) = temp(0,*)
		  transition.Rlu(i-1,*) = temp(1,*)
	 endfor

	 close, u
	 free_lun, u

 endif
end
