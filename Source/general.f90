module general
use variables
use maths
implicit none
contains
! *********************************************************
! *********************************************************
! GENERAL ROUTINES
! *********************************************************
! *********************************************************

! ---------------------------------------------------------
! Rotates a vector to put the index 1 into N and N into 1
! Vectors are real(kind=8)
! ---------------------------------------------------------
	function rotate(a)
	real(kind=8), INTENT(IN) :: a(:)
	real(kind=8) :: rotate(size(a))
	integer :: i, n
		n = size(a)
		do i = 1, n
			rotate(i) = a(n-i+1)
		enddo
	end function rotate

! ---------------------------------------------------------
! Read lines that don't mind in a file
! INPUT:
!     - unit : unit to whom the file is associated
!     - n : number of lines to read
! ---------------------------------------------------------
	subroutine lb(unit, n)
	integer, INTENT(IN) :: unit, n
	integer :: i
		do i = 1, n
			read(unit,*)
		enddo
	end subroutine lb

! ---------------------------------------------------------
! Calculates the partition function for a temperature using the
! given coefficients or using the energy levels
! INPUT:
!     - T : temperature
! OUTPUT:
!     Partition function like Sauval & Tatum (1984) or calculated using the
!     energy levels
! ---------------------------------------------------------
	function part_func(T)
	real(kind=8) :: part_func
	real(kind=8), INTENT(IN) :: T
	real(kind=8) :: theta, result
	integer :: i

		part_func = 0.d0
		if (n_partition /= 0) then
			theta = log10(5040.d0 / T)
			result = 0.d0

			do i = 1, n_partition
				result = result + partition(i) * theta**(i-1.d0)
			enddo

			part_func = 10**(result)
		else
			do i = 1, nl
				part_func = part_func + dlevel(2,i) * dexp(-dlevel(1,i) * PHK / T)  !g*exp(-h*nu/(k*T))
			enddo			
		endif
		
	end function part_func

! ---------------------------------------------------------
! Calculates LTE populations
! ---------------------------------------------------------	
	subroutine poplte
	real(kind=8), allocatable :: u(:), fi(:), fl(:)
	real(kind=8) :: kte, tp, fac, tot
	integer :: ip, i, ipl
	
		allocate(u(ni))
		allocate(fi(ni))
		allocate(fl(nl))

		do ip = 1, n_radios
			ipl = nl * (ip-1)
			tp = datph(1, ip)
			kte = PHK / tp			
			
			u = 0.d0

! Calculate the atomic partition function
			if (index(tipo_modelo,'ATOM') /= 0) then

				do i = 1, nl
					fl(i) = dlevel(2,i) * dexp(-dlevel(1,i) * kte)  !g*exp(-h*nu/(k*T))
					u(nli(i)) = u(nli(i)) + fl(i)
				enddo
				
			else if (index(tipo_modelo,'MOLECULE') /= 0) then

! Calculate the molecular partition function
				u(1) = part_func(tp)

				do i = 1, nl
					fl(i) = dlevel(2,i) * dexp(-dlevel(1,i) * kte)  !g*exp(-h*nu/(k*T))				
				enddo
				
			endif
			
! If the atom is multiionic
			if (ni > 1) then
				
				do i = 1, nl
					fl(i) = fl(i) / u(nli(i))
				enddo
				
				fac = datph(2,ip) * CI / (tp * dsqrt(tp))
				
				do i = 1, ni-1
					fi(i) = fac * u(i) / u(i+1) * dexp(dion(1,i+1)*kte)
				enddo
				
				fi(ni) = 1.d0
				tot = fi(ni)
				
				do i = ni - 1, 1, -1
					fi(i) = fi(i) * fi(i+1)
					tot = tot + fi(i)
				enddo
				tot = datph(3,ip) *  nh(ip) / tot		
				
				do i = 1, nl
					popl(i+ipl) = fl(i) * fi(nli(i)) * tot
				enddo
			else

! Molecules are not multiionic, so the program flow always takes this way
				tot = datph(3,ip) * nh(ip) / u(1)   ! abundance*nH/part_func

				do i = 1, nl
					popl(i + ipl) = fl(i) * tot      ! n_tot * degen / U * exp(-h*nu/k*T)
				enddo
				
			endif

		enddo 

		deallocate(u)
		deallocate(fi)
		deallocate(fl)
		
	end subroutine poplte

! ---------------------------------------------------------
! Calculates background radiation for every frequency (fuera(:))
! INPUT:
!     - transicion : transition number
! OUTPUT:
!     Puts background radiation in fuera(:)
! ---------------------------------------------------------
	subroutine calcula_fondo(transicion)
	integer, INTENT(IN) :: transicion
	integer :: i
	real(kind=8) :: frecu0, frecu

! Obtain central frequency of the transition
		frecu0 = dtran(2,transicion)
		fuera = 0.d0
		
		if (background_flag == 1) then
			do i = 1, nfrq
				frecu = 1.d0*(i - nfrq / 2 + 1) / n_freq_ud_doppler
! Obtain frequency axis
				frecu = frecu0 + doppler_width(n_radios,transicion) * frecu
! Planck function at temperature tback 
				fuera(i) = planck(frecu,tback) !2.d0 * PK * frecu**2 / PC**2 * tback
			enddo					
		endif

	end subroutine calcula_fondo
	
! ---------------------------------------------------------
! Calculates incident radiation in the central core due to a source
! It can be a blackbody or a greybody
! INPUT:
!     - transicion : transition number
! ---------------------------------------------------------
	subroutine calcula_centro(transicion)
	integer, INTENT(IN) :: transicion
	integer :: i
	real(kind=8) :: frecu0, frecu, opac, lam

		frecu0 = dtran(2,transicion)
		dentro = 0.d0
		
		if (index(central_source_type,'BLACKB') /= 0) then		
			if (central_flag == 2 .or. central_flag == 3) then
				do i = 1, nfrq
					frecu = 1.d0*(i - nfrq / 2 + 1) / n_freq_ud_doppler
! Obtain central frequency of the transition				
					frecu = frecu0 + doppler_width(1,transicion) * frecu
! Planck function at temperature tcentral				
					dentro(i) = planck(frecu,tcentral)
				enddo		
			endif
		endif
		
		if (index(central_source_type,'GREYB') /= 0) then		
			if (central_flag == 2 .or. central_flag == 3) then
				do i = 1, nfrq
					frecu = 1.d0*(i - nfrq / 2 + 1) / n_freq_ud_doppler
! Obtain central frequency of the transition				
					frecu = frecu0 + doppler_width(1,transicion) * frecu

					lam = 1.d4 * PC / frecu
					opac = ref_opacity * (ref_wavelength / lam) ** (spectral_index)

! Planck function at temperature tcentral	multiplied by the opacity factor
					dentro(i) = planck(frecu,tcentral) * (1.d0 - dexp(-opac))
				enddo		
			endif
		endif

	end subroutine calcula_centro
	
! ---------------------------------------------------------
! Calculates background radiation for every frequency (fuera(:))
! INPUT:
!     - transicion : transition number
! OUTPUT:
!     Puts background radiation in fuera(:)
! ---------------------------------------------------------
	subroutine calcula_fondo_bf(transicion)
	integer, INTENT(IN) :: transicion
	integer :: i
	real(kind=8) :: frecu

		fuera = 0.d0
		
		if (background_flag == 1) then
			do i = 1, nfrq_bf
				frecu = bf_cross_section(1,transicion,i)

! Planck function at temperature tback 
				fuera_bf(i) = planck(frecu,tback) !2.d0 * PK * frecu**2 / PC**2 * tback
			enddo					
		endif

	end subroutine calcula_fondo_bf
	
! ---------------------------------------------------------
! Calculates incident radiation in the central core due to a source
! It can be a blackbody or a greybody
! INPUT:
!     - transicion : transition number
! ---------------------------------------------------------
	subroutine calcula_centro_bf(transicion)
	integer, INTENT(IN) :: transicion
	integer :: i
	real(kind=8) :: frecu, opac, lam

		dentro = 0.d0
		
		if (index(central_source_type,'BLACKB') /= 0) then		
			if (central_flag == 2 .or. central_flag == 3) then
				do i = 1, nfrq_bf
					frecu = bf_cross_section(1,transicion,i)
					
! Planck function at temperature tcentral				
					dentro_bf(i) = planck(frecu,tcentral)
				enddo		
			endif
		endif
		
		if (index(central_source_type,'GREYB') /= 0) then		
			if (central_flag == 2 .or. central_flag == 3) then
				do i = 1, nfrq_bf
					frecu = bf_cross_section(1,transicion,i)

					lam = 1.d4 * PC / frecu
					opac = ref_opacity * (ref_wavelength / lam) ** (spectral_index)

! Planck function at temperature tcentral	multiplied by the opacity factor
					dentro_bf(i) = planck(frecu,tcentral) * (1.d0 - dexp(-opac))
				enddo		
			endif
		endif

	end subroutine calcula_centro_bf	

!----------------------------------------------------------
! Generates line profile
!----------------------------------------------------------
	subroutine genera_perfiles
	integer :: ir, ip, ifrq, inout, i
	real(kind=8) :: dir, angmu
	real(kind=8), allocatable :: x(:)

		allocate(x(nfrq))

		profile = 0.d0
		
		do i = 1, nfrq
			x(i) = 1.d0*(i - nfrq / 2 - 1 ) / n_freq_ud_doppler
			perfil(i) = 1.d0 / SQRTPI * voigt(0.d0,x(i),0)
		enddo		
		
! Going in and out
		do inout = 1, 2
		
			dir = (-1)**inout
! Going through all the shells
			do ir = 1, n_radios
			
! Going through all the intersections in each shell
				do ip = 1, n_cortes(ir)
				
					angmu = mu(ir,ip)

! All the frequency points
					do ifrq = 1, nfrq 
						x(ifrq) = 1.d0*(ifrq - nfrq / 2 - 1 ) / n_freq_ud_doppler - dir * angmu * datph(7,ir)
						profile(ifrq,ip,ir,inout) = 1.d0 / SQRTPI * voigt(0.d0,x(ifrq),0)

					enddo !ifrq
					
				enddo !ip
				
			enddo !ir
			
		enddo !inout
		

		deallocate(x)
		
	end subroutine genera_perfiles

! ---------------------------------------------------------
! Calculates frequency integration weights using a Gaussian integration scheme
! ---------------------------------------------------------
	subroutine freq_weights
	real(kind=8) :: div
	integer :: ir, ip, inout

		
		call genera_perfiles
		wtnu = 0.d0
		
		wtnu(2:nfrq-1) = 1.0d0
      wtnu(1) = 0.5d0
      wtnu(nfrq) = 0.5d0
		wtnu = 0.5d0 * wtnu

		
! Normalize the profile so that int(phi(x),x)=1
      div = 0.d0
      div = sum(wtnu(:) * perfil(:))
      div = 1.d0 / div
      wtnu = wtnu * div
		
! Going in and out
		do inout = 1, 2
		
! Going through all the shells
			do ir = 1, n_radios
			
! Going through all the intersections in each shell
				do ip = 1, n_cortes(ir)
					
! Set the frequency weights
					frqwt(2:nfrq-1,ip,ir,inout) = 1.d0
					frqwt(1,ip,ir,inout) = 0.5d0
					frqwt(nfrq,ip,ir,inout) = 0.5d0					
					frqwt(:,ip,ir,inout) = 0.5d0 * frqwt(:,ip,ir,inout)
					
! Normalize the profile so that int(phi(x),x)=1
					div = 0.d0
					div = sum(frqwt(:,ip,ir,inout) * profile(:,ip,ir,inout))
					div = 1.d0 / div
					frqwt(:,ip,ir,inout) = frqwt(:,ip,ir,inout) * div
					
				enddo !ip
				
			enddo !ir
			
		enddo !inout		
		
	end subroutine freq_weights

! ---------------------------------------------------------
! Returns angular integration weights. It is obtained using trapezoidal
! integration, and returns the weight for each intersection
! INPUT:
!     - capa : shell
!     - caract : characteristic. They define an intersection point
! OUTPUT:
!     - peso : returns the integration weight
! ---------------------------------------------------------	
	subroutine ang_weights(capa, caract, peso)
	integer, INTENT(IN) :: capa, caract
	real(kind=8), INTENT(INOUT) :: peso
	integer :: cortes

		peso = 0.d0
		cortes = n_cortes(capa)
		
		if ( (caract > 1) .and. (caract < cortes) ) then
			peso = 0.5d0 * dabs( mu(capa,caract+1) - mu(capa,caract-1) )
		else
			if (caract == 1) then
				peso = 0.5d0 * dabs( mu(capa,2) - mu(capa,1) )
			endif
			if (caract == cortes) then
				peso = 0.5d0 * dabs( mu(capa,cortes) - mu(capa,cortes-1) )
			endif
		endif
		
	end subroutine ang_weights
	
! ---------------------------------------------------------
! Returns impact parameter integration weights. It is obtained using trapezoidal
! integration, and returns the weight for each intersection
! INPUT:
!     - caract : characteristic. 
! OUTPUT:
!     - peso : returns the integration weight
! ---------------------------------------------------------	
	function impact_weights(caract)
	real(kind=8) :: impact_weights
	integer, INTENT(IN) :: caract
	real(kind=8) :: peso
	integer :: cortes

		peso = 0.d0
		cortes = n_caract
		
		if ( (caract > 1) .and. (caract < cortes) ) then
			peso = 0.5d0 * dabs( p(caract+1) - p(caract-1) )
		else
			if (caract == 1) then
				peso = 0.5d0 * dabs( p(2) - p(1) )
			endif
			if (caract == cortes) then
				peso = 0.5d0 * dabs( p(cortes) - p(cortes-1) )
			endif
		endif
		impact_weights = peso
	end function impact_weights	


! ---------------------------------------------------------
! Returns the number of intersections between a shell and all the characteristics
! INPUT:
!     - capa : shell where we want to know the number of intersections
! ---------------------------------------------------------	
	function n_cortes(capa)
	integer :: n_cortes
	integer, INTENT(IN) :: capa
		
		if (index(geometry_type,'SPHERICAL') /= 0) then
			if (capa == n_radios) then
				n_cortes = caract_core + capa - 1
			else 
				n_cortes = caract_core + capa
			endif		
		else
			n_cortes = angleset_pp
		endif
		
	end function n_cortes
	
	
! ---------------------------------------------------------
! Returns the time in seconds
! ---------------------------------------------------------
	function second()
   real :: second
   integer :: time_array(8)
   real :: tiempo
   
	call date_and_time(values=time_array)
	tiempo = time_array(5)*3600.*1000.+time_array(6)*60.*1000.+time_array(7)*1000.+time_array(8)
	second = tiempo / 1000.
   	
	end function second

! ---------------------------------------------------------
! Returns Planck's function for a frequency and temperature in cgs
! INPUT:
!     - nu : frequency
!     - T : temperature
! ---------------------------------------------------------		  
	function planck(nu,T)
	real(kind=8) :: planck
	real(kind=8), INTENT(IN) :: nu, T
	real(kind=8) :: temporal, exponen
		
		planck = 0.d0

		if (T /= 0.d0) then
			exponen = PHK * nu / T
			if (exponen < 50) then
				temporal = dexp(exponen) - 1.d0	  
				planck = PHC * nu**3.d0 / temporal
			else
				planck = 0.d0
			endif
		else
			planck = 0.d0
		endif	
		
	end function planck

! ---------------------------------------------------------
! Returns the equivalent temperature of a blackbody which emits the intensity inten at frequency nu
! INPUT:
!     - nu : frequency
!     - inten : intensity in frequency units
! ---------------------------------------------------------		  
	function T_planck_equiv(nu,inten)
	real(kind=8) :: T_planck_equiv
	real(kind=8), INTENT(IN) :: nu, inten
	real :: temporal
		
		T_planck_equiv = 0.d0

		if (inten /= 0.d0) then
			temporal = 1.e0 + PHC * nu**3.e0 / inten
			if (temporal > 0.d0) then
				T_planck_equiv = PHK * nu / alog(temporal)
			else
				T_planck_equiv = 0.d0
			endif
		else
			T_planck_equiv = 0.d0
		endif	
		
	end function T_planck_equiv

! ---------------------------------------------------------
! Returns line opacity, dust opacity, continuum opacity, line source function,
! dust source function and continuum source function for a transition
! INPUT:
!     - it : transition
! OUTPUT:
!     - min_tau : minimum optical depth to be considered as thin
!     - chil, kappa, kappa_dust, Sl, B, B_dust
! ---------------------------------------------------------		  

	subroutine opacity(it, min_tau, from, to)
	integer, INTENT(IN) :: it, from, to
	integer, INTENT(INOUT) :: min_tau
	real(kind=8) :: glu, acon, sig, tau0, chim, chip
	integer :: up, low, ip1, ip, ipl, from2

! Initialize opacities and source functions to 0
		kappa = 0.d0
		kappa_dust = 0.d0
		chil = 0.d0
		B = 0.d0
		B_dust = 0.d0		
		Sl = 0.d0		
		
		up = itran(1,it)
		low = itran(2,it)
		glu = dlevel(2,low) / dlevel(2,up)
		acon = PHC * (dtran(2,it))**3.d0  !acon = (2*h*nu**3)/c**2

		tau0 = 0.d0
		chip = 0.d0
		chim = 0.d0
		ip1 = 1

		from2 = from
		if (from == 0) then 
			from2 = 1
		endif
				
		do ip = from2, to
			ipl = nl * (ip-1)
			sig = dtran(1,it) / doppler_width(ip,it) !sig=(h*nu*Blu)/(4*pi*Dnudoppler)

! Line opacity and line source function
 			chil(ip) = sig * (pop(low+ipl) - glu * pop(up+ipl))

! MASER SOLUTION???
			if (chil(ip) < 0.d0) then
!				print *, 'Negative absorption detected'
				chil(ip) = -chil(ip)
			endif
			Sl(ip) = acon * glu * pop(up+ipl) / (pop(low+ipl) - glu * pop(up+ipl))

! MASER SOLUTION???
			if (Sl(ip) < 0.d0) then
!				print *, 'Negative source function detected'
				Sl(ip) = 0.d0
			endif
			
! Continuum opacity and continuum source function (not dust)
			if (back_opac_flag == 1) then
				kappa(ip) = chic(it,ip)
!				B(ip) = acon * glu * popl(up+ipl) / (popl(low+ipl) - glu * popl(up+ipl))
				B(ip) = Sc(it,ip)
			endif
				
			
! Dust opacity and dust source function 
			if (dust_flag /= 0) then			
				call dust_opacity(ip,it)
			endif									

		enddo	!ip

! If you want to connect this feature, comment the following line
		ip1 = 1
		min_tau = n_radios + 1 - ip1
		
	end subroutine opacity
	
! ---------------------------------------------------------
! Returns line opacity and line source function for a transition and a shell
! INPUT:
!     - it : transition
! OUTPUT:
!     - chi_out : opacity
!		- S_out : source function
! ---------------------------------------------------------		  

	subroutine opacity_one(it, ip, chi_out, S_out)
	integer, INTENT(IN) :: it, ip
	real(kind=8), INTENT(INOUT) :: chi_out, S_out
	real(kind=8) :: glu, acon, sig
	integer :: up, low, ipl

! Initialize opacities and source functions to 0
		chi_out = 0.d0
		S_out = 0.d0
		
		up = itran(1,it)
		low = itran(2,it)
		glu = dlevel(2,low) / dlevel(2,up)
		acon = PHC * (dtran(2,it))**3.d0  !acon = (2*h*nu**3)/c**2
		
		ipl = nl * (ip-1)
		sig = dtran(1,it) / doppler_width(ip,it) !sig=(h*nu*Blu)/(4*pi*Dnudoppler)

! Line opacity and line source function
		chi_out = sig * (pop(low+ipl) - glu * pop(up+ipl))			
		S_out = acon * glu * pop(up+ipl) / (pop(low+ipl) - glu * pop(up+ipl))		
		
	end subroutine opacity_one
	
! ---------------------------------------------------------
! Returns dust opacity for indicated wavelength
! ---------------------------------------------------------		  

	subroutine dust_opacity(ip, it)
	integer, INTENT(IN) :: ip, it
	real(kind=8) :: lamb, lamb_ref, dielectric(2), opac_graf_par, opac_graf_per, opac_silic
 	complex(kind=8) :: eps_graf_par, eps_graf_per, eps_silic, uno, dos
	real(kind=8) :: mdust_mgas, dust_density, dust_radius, Q0, exponent
	
		uno = cmplx(1.d0,0.d0)
		dos = cmplx(2.d0,0.d0)
	
! Central wavelength in cm		
		lamb = (PC / dtran(2,it))

		if (dust_flag == 1) then
! kappa_dust = A_v * (lambda/lambda_ref)^(-expav)
			lamb_ref = (PC / lambda_ref)
			kappa_dust(ip) = 1.d-21 * nh(ip) * (lamb / lamb_ref)**(-expav) 
! Planck's function for dust source function
			B_dust(ip) = planck(dtran(2,it),datph(6,ip))			
		endif
				
		if (dust_flag == 2) then
			dielectric(1) = interpolate(graphite_lambda,eps_graphite_par(1,:),lamb*1.d4)
			dielectric(2) = interpolate(graphite_lambda,eps_graphite_par(2,:),lamb*1.d4)
			eps_graf_par = cmplx(dielectric(1),dielectric(2))
			
			dielectric(1) = interpolate(graphite_lambda,eps_graphite_per(1,:),lamb*1.d4)
			dielectric(2) = interpolate(graphite_lambda,eps_graphite_per(2,:),lamb*1.d4)
			eps_graf_per = cmplx(dielectric(1),dielectric(2))
			
			dielectric(1) = interpolate(silicate_lambda,eps_silicate(1,:),lamb*1.d4)
			dielectric(2) = interpolate(silicate_lambda,eps_silicate(2,:),lamb*1.d4)
			eps_silic = cmplx(dielectric(1),dielectric(2))


			opac_graf_par = 1.d0 / (UMA * 1.00794d0) * 10.d0**C_graf * 8.d0 * PI**2 / lamb *&
				imaginary((eps_graf_par-uno)/(eps_graf_par+dos)) * (dsqrt(0.25d-4)-dsqrt(0.005d-4))
				
			opac_graf_per = 1.d0 / (UMA * 1.00794d0) * 10.d0**C_graf * 8.d0 * PI**2 / lamb *&
				imaginary((eps_graf_per-uno)/(eps_graf_per+dos)) * (dsqrt(0.25d-4)-dsqrt(0.005d-4))
				
			opac_silic = 1.d0 / (UMA * 1.00794d0) * 10.d0**C_sil * 8.d0 * PI**2 / lamb *&
				imaginary((eps_silic-uno)/(eps_silic+dos)) * (dsqrt(0.25d-4)-dsqrt(0.005d-4))								

! Maybe it is chi=kappa*nh2*mu/m_gas_dust				
			kappa_dust(ip) = nh(ip) * 2.01588d0 * UMA * (2.d0 / 3.d0 * opac_graf_per + &
				1.d0 / 3.d0	* opac_graf_par + opac_silic)
			B_dust(ip) = planck(dtran(2,it),datph(6,ip))
				
		endif
		
! See the thesis of Fabrice Herpin
		if (dust_flag == 3) then
! kappa_dust = n_d*sigma_d*Q
! Reference wavelength = 80 microns = 80.d-4 cm
			lamb_ref = 80.d-4
!			kappa_dust(ip) = (2.d0 * 1.d-2 * UMA * nh(ip) * 3.d0) / (4.d0 * 3.d0 * 5.d-6) *&
!				2.d-3 * (lamb / lamb_ref)**(-1.1d0) 

			mdust_mgas = 1.d-2
			dust_density = 3.d0   ! rho_d = 3 g cm^-3
			dust_radius = 5.d-6   ! a = 0.05 microns
			Q0 = 2.d-3 				 ! Q0 = 2.d-3  (efficiency)
			exponent = -1.1d0
			kappa_dust(ip) = (2.d0 * mdust_mgas * UMA * nh(ip) * 3.d0) / &
				(4.d0 * dust_density * dust_radius) * Q0 * (lamb / lamb_ref)**(exponent) 
				
! Planck's function for dust source function
			B_dust(ip) = planck(dtran(2,it),datph(6,ip))			
		endif		
		
	end subroutine dust_opacity	

end module general
