pro read_atmos, fichero
@vars.common
	 get_lun, u
	 openr, u, fichero
	 lb, u, 4
; Read number of radial shells
	 readf, u, nr

	 lb, u, 1
; Read number of characteristics
	 readf, u, np
	 
	 atmos = {N_radius: nr, N_charac: np, r: dblarr(nr), T: dblarr(nr), NH2: dblarr(nr), $
	 	  N: dblarr(nr), N_e: dblarr(nr), V_mic: dblarr(nr), T_dust: dblarr(nr), V_mac: dblarr(nr), $
		  A_v: dblarr(nr), p: dblarr(np), mu_exit: dblarr(np) }
		  
; Read radial variables
	 temp = dblarr(9,nr)	 
	 lb, u, 5
	 readf, u, temp
	 
	 atmos.r = temp(0,*)
	 atmos.T = temp(1,*)
	 atmos.NH2 = temp(2,*)
	 atmos.N = temp(3,*)
	 atmos.N_e = temp(4,*)
	 atmos.V_mic = temp(5,*)
	 atmos.T_dust = temp(6,*)
	 atmos.V_mac = temp(7,*)
	 atmos.A_v = temp(8,*)
	 
	 lb, u, 4
; Read characteristics and external mu
	 temp = dblarr(2,np)
	 readf, u, temp
	 
	 atmos.p = temp(0,*)
	 atmos.mu_exit = temp(1,*)
	 
	 close, u
	 free_lun, u
end
