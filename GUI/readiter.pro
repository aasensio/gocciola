pro read_iter, fichero
@vars.common
	 n_iter = n_lineas(fichero)
	 
	 iterat = {Niter: n_iter, rel_change: dblarr(n_iter) }
	 
	 get_lun, u
	 openr, u, fichero
	 temp = dblarr(2,n_iter)
	 readf, u, temp
	 
	 iterat.rel_change = temp(1,*)
	 
	 close, u
	 free_lun, u 
end
