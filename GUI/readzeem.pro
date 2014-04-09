pro read_zeeman, fichero, salida
@vars.common
	 get_lun, u
	 salida = 0
	 z = findfile(fichero) eq ''

	 if (z(0) eq 0) then begin
		 openr, u, fichero	 

	; Read the number of synthesized lines	
		 readf, u, n_zeem_lines

	; Read the first wavelength	 
		 readf, u, wave

	; Read the number of frequencies of the line	 
		 readf, u, n_zeem_frqs

		 zeeman = {n_lines: n_zeem_lines, n_frq: n_zeem_frqs, lambda: dblarr(n_zeem_lines,n_zeem_frqs), $
	 		  SI: dblarr(n_zeem_lines,n_zeem_frqs), SQ: dblarr(n_zeem_lines,n_zeem_frqs), $
			  SU: dblarr(n_zeem_lines,n_zeem_frqs), SV: dblarr(n_zeem_lines,n_zeem_frqs), $
			  lambda0 : dblarr(n_zeem_lines) }

		 temp2 = dblarr(5,n_zeem_frqs)
		 for i = 0, n_zeem_lines-1 do begin
	 		  readf, u, temp2
			  zeeman.lambda0(i) = wave
			  zeeman.lambda(i,*) = wave + temp2(0,*) / 1.d3
			  zeeman.SI(i,*) = temp2(1,*)
			  zeeman.SQ(i,*) = temp2(2,*)
			  zeeman.SU(i,*) = temp2(3,*)
			  zeeman.SV(i,*) = temp2(4,*)
			  if (i ne n_zeem_lines-1) then begin
		   		readf, u, wave, n_zeem_frqs 
			  endif		  
		 endfor
		 close, u
		 salida = 1
	 endif
end
