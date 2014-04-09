pro read_backopa, fichero
@vars.common
	 
	 if (atom.BackOpacity eq 1) then begin
		 get_lun, u
		 openr, u, fichero

		 temp = dblarr(geometry.Nradius)

		 for i = 0, transition.Ntran-1 do begin
	 		  readf, u, temp
			  atom.backopa(i,*) = temp
		 endfor

		 close, u
		 free_lun, u 
	 endif
end
