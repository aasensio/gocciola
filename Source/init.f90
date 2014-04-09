module initial
use variables
use atomic
use atmosfer
use general
use background_opac
use overlap
use data_write
use maths
implicit none
contains

! *********************************************************
! *********************************************************
! INITIALIZATION ROUTINES
! *********************************************************
! *********************************************************

	
! ---------------------------------------------------------
! Initialize all variables
! ---------------------------------------------------------

   subroutine init

   integer :: i, arr_size
   real(kind=8), allocatable :: x_quadr(:), w_quadr(:)
   		
		print *, 'Opening configuration file (config.dat)'
		open(unit=2,file='config.dat',action='read',status='old')
		call lb(2,4)
		
      ! Read if verbose mode is active or not
		read(2,*) verbose_flag
      call lb(2,2)
      
		! Read if populations are going to be calculated
		read(2,*) converge_flag
		select case(converge_flag)
			case(0)
				print *, 'Calculating populations...'
			case(1)
				print *, 'Reading populations from file...'	
			case(2)
				print *, 'Synthesizing LTE spectrum'
			case(3)
				print *, 'Read population from a previous run'
			case(4)
				print *, 'Read population from a previous run and do Zeeman synthesis of a line'
		end select
		
		call lb(2,2)
! Read if the atmosphere is interpolated when the spectrum is obtained
		read(2,*) interpolate_atmosphere
		if (interpolate_atmosphere /= 0) then
			if (converge_flag == 1) then
				print *, 'Interpolating the atmosphere for the calculation of the emerging spectrum'
			else
! If interpolate_atmosphere/=0 and not calculating spectrum, let it be 0 because we use the original atmosphere
				interpolate_atmosphere = 0
			endif
		endif

		call lb(2,2)
! Read initialization method
		read(2,*) geometry_type
		if (index(geometry_type,'SPHERICAL') /= 0) then
			print *, 'Solving a SPHERICAL problem'
		endif
		if (index(geometry_type,'PLANEP') /= 0) then
			print *, 'Solving a PLANE-PARALLEL problem'
		endif
		
		call lb(2,2)		
! Read angle set for plane-parallel
		read(2,*) angleset_pp
		if (index(geometry_type,'PLANEP') /= 0) then
			print *, 'Using ',angleset_pp,'-angle quadrature'
			print *, 'Calculating the appropiate angular weights'
			allocate(mu_pp(angleset_pp))
			allocate(wtmu(angleset_pp))
			allocate(x_quadr(2*angleset_pp))
			allocate(w_quadr(2*angleset_pp))

			call gauleg(-1.d0,1.d0,x_quadr,w_quadr,2*angleset_pp)
			mu_pp = x_quadr(1:angleset_pp)
			wtmu = w_quadr(1:angleset_pp)
			
			deallocate(x_quadr)
			deallocate(w_quadr)
		endif
		
		call lb(2,2)		
! Read initialization method
		read(2,*) lvg_flag
		if (lvg_flag == 0.d0) then
			print *, 'Initializing problem with LTE'
		else
			print *, 'Initializing problem with LVG with precision : ', lvg_flag
		endif
		
		call lb(2,2)						
! Read number of characteristics going through the core
		read(2,*) caract_core
		print *, 'Number of characteristics going through core : ',caract_core
				
		call lb(2,2)
! Read the flag to activate dust 
		read(2,*) dust_flag
		select case(dust_flag)
			case(0) 
				print *, 'Dust opacity not included'
			case(1)
				print *, 'Dust opacity included in a power law model'
			case(2) 
				print *, 'Using realistic dust opacity'

! Read complex refractive index (dielectric function) for graphite
				open(unit=3,file='Dust/graphite.dat',action='read',status='old')
				print *, ' - Including graphite grains'
				call lb(3,2)
				read(3,*) arr_size
				allocate(eps_graphite_par(2,arr_size))
				allocate(eps_graphite_per(2,arr_size))
				allocate(graphite_lambda(arr_size))
				do i = 1, arr_size
					read(3,*) graphite_lambda(i), eps_graphite_par(1,i), eps_graphite_par(2,i),&
						eps_graphite_per(1,i), eps_graphite_per(2,i)
					graphite_lambda(i) = graphite_lambda(i) * 1.d-4    ! in cm
				enddo
				close(3)
				
				open(unit=3,file='Dust/silicate.dat',action='read',status='old')
				print *, ' - Including silicate grains'
				call lb(3,2)
				read(3,*) arr_size
				allocate(eps_silicate(2,arr_size))
				allocate(silicate_lambda(arr_size))
				do i = 1, arr_size
					read(3,*) silicate_lambda(i), eps_silicate(1,i), eps_silicate(2,i)
					silicate_lambda(i) = silicate_lambda(i) * 1.d-4    ! in cm
				enddo				
				close(3)				
				
		end select
		
		call lb(2,2)
! Read the flag to activate background radiation		
		read(2,*) background_flag
		if (background_flag == 1) then 
			print *, 'Including microwave background radiation' 
		else
			 print *, 'Microwave background radiation not included'
		endif
		
		call lb(2,2)
! Read the flag to activate background opacity
		read(2,*) back_opac_flag
		if (back_opac_flag == 1) then 
			print *, 'Including background opacity' 
		else
			 print *, 'Background opacity not included'
		endif
		
		call lb(2,2)
! Read the flag to activate central source
! FALTA VER QUE PASA SI central_flag==3
		read(2,*) central_flag
		select case(central_flag)
			case(3)
				print *, 'Including a central source not filling the central core'
				repeat = 1.d0		
			case(2)
				print *, 'Including a central source filling the central core'
				repeat = 1.d0
				repeat = 2.d0
			case(1)
				print *, 'Including an empty core'
				repeat = 2.d0
			case(0)
				print *, 'Including an absorbing core'
				repeat = 1.d0
		end select
					
		call lb(2,2)
! Read background temperature
		read(2,*) tback
		print *, 'Background temperature : ', tback

		call lb(2,2)
! Read the type of central source's spectrum
		read(2,*) central_source_type
		print *, 'Using a central body of type : ', central_source_type		
		
		call lb(2,2)
! Read the reference data for a greybody source
		read(2,*) ref_wavelength, ref_opacity, spectral_index
		if (index(central_source_type,'GREYB') /= 0) then
			print *, 'Ref. wavelength : ',ref_wavelength, ' - Ref. opacity : ',ref_opacity,&
				' - Spectral index : ',spectral_index
		endif
		
		call lb(2,2)
! Read central source temperature
		read(2,*) tcentral
		print *, 'Central source temperature : ', tcentral
		
		call lb(2,2)
! Read atmosphere file
		read(2,*) nombre_atmosf
		print *, 'Reading atmosphere : ', nombre_atmosf
		call read_atmosf(nombre_atmosf)

! Interpolate the atmosphere if needed (only when calculating the emerging spectrum)		
		if (interpolate_atmosphere /= 0) then
			call interpol_atmosphere
		endif
		
		call lb(2,2)
! Read distance to the source
		read(2,*) distancia
		distancia = distancia * PARSEC
		print *, 'Distance : ', distancia, ' cm = ', distancia/PARSEC, ' pc'
		
		call lb(2,2)
! Read angular size
		read(2,*) ang_size
!		ang_size = 2.d0 * datan(r(n_radios) / distancia) * 206265.d0
		ang_size = ang_size / 2.d0	
		print *, 'Angular size (radius): ', ang_size, ' arcsec'

		call lb(2,2)
! Read stellar radius
		read(2,*) star_radius
		print *, 'Star radius : ', star_radius, ' cm'
		
		do i = 1, caract_core
! Get the number of impact parameters crossing the central star in boundary type 3
			if (p(i) < star_radius) impact_star_radius = i
		enddo
		if (central_flag == 3) then
			print *, 'Maximum impact parameter which goes through the central body ',impact_star_radius, p(impact_star_radius)
		endif
		
		call lb(2,2)
! Read telescope's sigma
		read(2,*) tel_fwhm, tel_fwhm_cut
		print *, 'Telescope''s full width half maximum (FWHM) : ', tel_fwhm
		if (tel_fwhm_cut /= 0) then
			print *, 'Cutting the observing beam at ',tel_fwhm_cut,' arcsec'
		else
! If tel_fwhm_cut=0 then no cut is wanted
			tel_fwhm_cut = 1.d80
		endif
		
		call lb(2,2)
! Read observing instrument
		read(2,*) observing_instrument
		print *, 'Observing instrument : ', observing_instrument
      
		call lb(2,2)
! Read if output spectrum has to be saved
		read(2,*) output_spectrum
		
		call lb(2,2)
! Read if radiative rates are saved
		read(2,*) radrates_flag
		
		call lb(2,2)
! Read if overlapping of lines are taken into account
		read(2,*) overlap_flag

		call lb(2,2)
! Read if collisional coefficientes are going to be used
		read(2,*) include_collis
		
		call lb(2,2)
! Read output file
		read(2,*) nombre_salida						
		print *, 'Output file : ', nombre_salida
		
		call lb(2,2)
! Read iteration file
		read(2,*) nombre_iteration
		print *, 'Iteration file : ', nombre_iteration							
		
		call lb(2,2)
! Read emerging spectrum file
		read(2,*) nombre_emerging
		print *, 'Emerging spectrum output file : ', nombre_emerging		

		call lb(2,2)
! Read emerging flux output file
		read(2,*) nombre_flux
		print *, 'Emerging flux output file : ', nombre_flux				

		call lb(2,2)		
! Read iteration file
		read(2,*) nombre_atmosf_out
		print *, 'Atmosphere output file : ', nombre_atmosf_out
		
		call lb(2,2)		
! Read background opacity file
		read(2,*) nombre_background_out
		print *, 'Background opacity output file : ', nombre_background_out		

		call lb(2,2)		
! Read radiative rates file
		read(2,*) nombre_radrates
		print *, 'Radiative rates output file : ', nombre_radrates	
		
		call lb(2,2)
! Read collisional rates file
		read(2,*) nombre_collisrates
		print *, 'Collisional rates output file : ', nombre_collisrates

		call lb(2,2)		
! Read maximum number of iterations
		read(2,*) iter_max
		print *, 'Maximum number of iterations : ', iter_max

		call lb(2,2)
! Read number of iterations before any acceleration
		read(2,*) itracc
		itracc = -max(0,itracc)
		print *, 'Number of iterations before any acceleration : ', -itracc

		call lb(2,2)
! Read number of iterations before re-acceleration
		read(2,*) nintracc
		nintracc = max(0,nintracc)
		print *, 'Number of iterations before re-acceleration : ', nintracc
		
		call lb(2,2)
! Read acceleration order
		read(2,*) nord
		nord = min(nordm,nord)
		print *, 'Acceleration order : ', nord
		
		call lb(2,2)
! Type of acceleration (Ng or no acceleration)
		read(2,*) iacel
		if (iacel == 1) then 
			print *, 'Type of acceleration : NG'
		else
			print *, 'No acceleration'
		endif

		call lb(2,2)
! Value to turn on polynomic expansion of exponentials
		read(2,*) trick
		print *, 'Value to turn on polynomic expansion of exponentials : ', trick
		
		call lb(2,2)
! Precision
		read(2,*) precision
		print *, 'Solution precision (final maximum relative change) : ', precision		

		
		call lb(2,2)
! Iterative improvement of the solution
		read(2,*) iterat_improv
		if (iterat_improv == 0) then
			print *, 'Not using iterative improvement of the solution to the linearized system'
		else
			print *, 'Using iterative improvement of the solution to the linearized system'
		endif
		
		call lb(2,2)
! Linear solver algorithm
		read(2,*) linear_solver_algorithm
		if (linear_solver_algorithm == 0) then
			print *, 'Using LU solver for the linear systems'
		else
			print *, 'Using SVD solver for the linear systems'
		endif
		
		call lb(2,2)
! SNTB acceleration
		read(2,*) SNTB_flag
		if (SNTB_flag == 0) then
			print *, 'Not using SNTB acceleration'
		else
			print *, 'Using SNTB acceleration'
		endif		
		
		call lb(2,2)
! Optical depth to be considered as thin		
		read(2,*) tauip1
		print *, 'Optical depth to be considered as thin : ', tauip1
		
		call lb(2,2)
! Number of points in the gaussian quadrature used in LVG for the calculation of the escape prob.
		read(2,*) gaussq_n
		print *, 'Number of points for the LVG gaussian quadrature : ', gaussq_n

		call lb(2,2)
! File for the atomic/molecular model
		read(2,*) nombre_modelo
		print *, 'Using atomic/molecular model : ', nombre_modelo
		
		call lb(2,2)
		read(2,*) include_bound_free

		print *, 'Calculating angular network for intersections between shells and characteristics'
		if (index(geometry_type,'SPHERICAL') /= 0) then
			allocate(mu(n_radios,n_caract))
			call mugrid
		else
			allocate(mu(n_radios,angleset_pp))
			do i = 1, n_radios
				mu(i,:) = abs(mu_pp)
			enddo
		endif
			
		
! Allocate memory for all variables depending on the size of the grid
		allocate(B(n_radios))
		allocate(B_dust(n_radios))
		allocate(chil(n_radios))
		allocate(kappa(n_radios))
		allocate(kappa_dust(n_radios))
!		allocate(chie(n_radios))
		allocate(Sl(n_radios))
		allocate(Lstar(n_radios))
		allocate(Jbar(n_radios))
		
		


		call lb(2,2)
! Read number of Doppler widths for the profile
		read(2,*) n_freq_perfil
		print *, 'Number of Doppler widths to sample the line profile : ', n_freq_perfil

		call lb(2,2)
! Read number of frequency points for each Doppler width
		read(2,*) n_freq_ud_doppler
		print *, 'Number of frequency points to sample each Doppler width : ', n_freq_ud_doppler

		call lb(2,2)
! Read reference wavelength
		read(2,*) lambda_ref
		print *, 'Reference wavelength : ', lambda_ref,' microns'

! Put it in s^(-1)
		lambda_ref = PC / (lambda_ref * 1.d-4)
		
		call lb(2,2)
! Read dust opacity wavelength dependence exponent
		read(2,*) expav
		print *, 'Dust opacity wavelength dependence exponent : ', expav

		call lb(2,2)
! Read helium abundance (only for a molecular model)
		read(2,*) xhel		
		if (index(tipo_modelo,'MOLECULE') /= 0) then
			print *, 'Helium abundace : ', xhel
		endif
		
		call lb(2,2)
! Read atomic hydrogen abundance (only for a molecular model)
		read(2,*) xhid
		if (index(tipo_modelo,'MOLECULE') /= 0) then
			print *, 'Atomic hydrogen abundance : ', xhid
		endif

! Reading atomic/molecular model
		print *, 'Reading atomic/molecular model'
		call atom_model(nombre_modelo)
		
! Quantities for the bound-free case
		if (ni > 1) then
			allocate(Jbarb(n_radios))
			allocate(Lstarb(n_radios))
			allocate(Lstarc(n_radios))
			allocate(Ibf(n_radios))
			
			allocate(Jbar_total_bf_a(n_radios*nbf))
			allocate(Jbar_total_bf_b(n_radios*nbf))
			
			allocate(Lstar_total_bf_a(n_radios*nbf))
			allocate(Lstar_total_bf_b(n_radios*nbf))
			allocate(Lstar_total_bf_c(n_radios*nbf))
			allocate(I_total_bf(n_radios*nbf))
		endif			
				
		allocate(Jbar_total(n_radios*nact))
		allocate(lstar_total(n_radios*nact))
		allocate(scratch(nl*n_radios,nord+2))
		allocate(pop(nl*n_radios))
		allocate(popl(nl*n_radios))
		allocate(popi(nl*n_radios))
		allocate(tau(nact,n_radios))
		allocate(chic(nact,n_radios))
		allocate(Sc(nact,n_radios))
		
		allocate(chic_bf(nt,nfrq_bf,n_radios))
		allocate(Sc_bf(nt,nfrq_bf,n_radios))

! Maximum velocity in the atmosphere
		V_max = maxval(dabs(datph(7,:)))
		print *, 'Maximum velocity in Doppler units : ', V_max
		
! Total frequency points (taking account of the maximum Doppler displacement)
		nfrq = n_freq_perfil * n_freq_ud_doppler + 1 + 2 * V_max * n_freq_ud_doppler
		allocate(profile(nfrq,n_caract,n_radios,2))
		allocate(frqwt(nfrq,n_caract,n_radios,2))
		allocate(overl(nact,nact,n_radios))  
!		allocate(overlap_chi(nfrq,n_radios))
		allocate(In(nfrq,n_caract))
		allocate(wtnu(nfrq))
		allocate(perfil(nfrq))
		allocate(dentro(nfrq))
		allocate(fuera(nfrq))
		allocate(dentro_bf(nfrq_bf))
		allocate(fuera_bf(nfrq_bf))

! Calculating frequency integration weights
		print *, 'Calculating frequency integration weights...'
		call freq_weights
		
! LTE initialization
		print *, 'LTE initialization'
		call poplte
		pop(1:nl*n_radios) = popl(1:nl*n_radios)		
		popi(1:nl*n_radios) = popl(1:nl*n_radios)	
		
! Calculate optical depth for each line
		print *, 'Calculating optical depth for each line in LTE'
		call calc_tau
		
! Calculate background opacity
		if (back_opac_flag == 1) then
			print *, 'Calculating background opacity'
			call background
		endif		

! Read type of error to calculate
		call lb(2,2)
		read(2,*) error_flag
		select case(error_flag)
			case(0)
				print *, 'Using convergence error'
			case(1)
				print *, 'Using true error (needs a previous run with the same configuration)'
				print *, 'If not, it may crash'
				allocate(pop_final(nl*n_radios))
				call read_pops_true_error(pop_final)
		end select
! Read iterative scheme to do
		call lb(2,2)
		read(2,*) tipo_iteracion
		select case(tipo_iteracion)
			case(1)!'LAMBDA') 
				omega_sor = 1.d0
				print *, 'Starting iterative process with Lambda-iteration...'
			case(2)!'LAMBDAJAC')
				omega_sor = 1.d0
				print *, 'Starting iterative process with Lambda and then Jacobi...'
			case(3)!'JACOBI')
				omega_sor = 1.d0
				print *, 'Starting iterative process with Jacobi iteration...'
			case(4)!'GS')
				omega_sor = 1.d0
				print *, 'Starting iterative process with Gauss-Seidel iteration...'
			case(5)!'SOR')
				print *, 'Starting iterative process with Successive Overrelaxation iteration...'
				call lb(2,2)
				read(2,*) omega_sor
				print *, 'Omega parameter : ', omega_sor
			case(6)!'SEMIGS')
				omega_sor = 1.d0
				print *, 'Starting iterative process with Semi Gauss-Seidel iteration...'
			case(7)!'SEMISOR')
				print *, 'Starting iterative process with Semi Successive Overrelaxation iteration...'
				call lb(2,2)
				read(2,*) omega_sor
				print *, 'Omega parameter : ', omega_sor
			case(8)!'CEP')
				print *, 'Starting iterative process with CEP method...'
! Allocate memory for the derivatives of Jbar and the source function				
				allocate(dJbar_totaldn(n_radios*nact,n_radios*nl))
				allocate(dJbardn(n_radios,n_radios*nl))
				allocate(dSldn(n_radios,n_radios*nl))
				n_quadr_beta = 80
      		allocate(x_e3(n_quadr_beta))
				allocate(w_e3(n_quadr_beta))
      		call gauleg(-7.d0,7.d0,x_e3,w_e3,n_quadr_beta)
		end select
      
				
! Close file
		close(2)
		
! If overlapping is included, then initialize the variables
		if (overlap_flag == 1) then
			print *, 'Initializing overlapping...'
			call init_overlap
		endif
					
   end subroutine init
end module initial
