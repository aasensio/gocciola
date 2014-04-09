module lvg_func
use variables
use statis
use atmosfer
use data_write
use maths
implicit none
contains

! *********************************************************
! *********************************************************
! LVG INITIALIZATION ROUTINES
! *********************************************************
! *********************************************************

!-----------------------------------------------------------------
! Calculates Jbar using the LVG approximation
!-----------------------------------------------------------------
	subroutine LVGformal_solver(itr)
	integer, INTENT(IN) :: itr
	integer :: idir, dir, k, k0, k1, kdel, cortes, inout
	real(kind=8) :: Q, tau_lvg, integrand, beta, derv, mu_c
	real(kind=8), allocatable :: chi0(:), s0(:)
	real(kind=8) :: opac, source
	real(kind=8), allocatable :: vel(:)
	real(kind=8), allocatable :: gauss_weights(:), gauss_abs(:)

	allocate(chi0(nfrq))
	allocate(s0(nfrq))
	allocate(gauss_weights(gaussq_n))
	allocate(gauss_abs(gaussq_n))
	
	Jbar = 0.d0
	
	call gauleg(0.d0,1.d0,gauss_abs,gauss_weights,gaussq_n)		
	
	allocate(vel(n_radios))

! Calculate outer boundary condition for this transition
	call calcula_fondo(itr)
	call calcula_centro(itr)

	idir = -1 ; k0 = n_radios ; k1 = 1 ; kdel = -1 ; inout = 1
	
   vel = datph(7,:)
   
! If the velocity profile is 0. PUT THIS IN A BETTER WAY

	if (maxval(datph(7,:)) == 0.d0 .and. minval(datph(7,:)) == 0.d0) then
   	do k = 1, n_radios
      	vel(k) = (k-1.d0) / (n_radios-1.d0)
      enddo
   endif
		
! Go through all the shells inside out
	do k = 1, n_radios
	
		beta = 0.d0		
		mu_c = dsqrt(1.d0 - r(1)**2/r(k)**2)
		
! Calculate the number of intersections between the characteristics and the shell
		cortes = n_cortes(k)
		
! Calculate radial velocity field gradient.
		if (k /= 1 .and. k /= n_radios) then
!			derv = lag_deriv(r(k-1),r(k),r(k+1),vel(k-1),vel(k),vel(k+1),1)
			derv = deriv(r(k-1),r(k),r(k+1),vel(k-1),vel(k),vel(k+1),1)
		else
			if (k == 1) then
!				derv = lag_deriv(r(1),r(2),r(3),vel(1),vel(2),vel(3),0)
				derv = deriv(r(1),r(2),r(3),vel(1),vel(2),vel(3),0)				
			else
				if (k == n_radios) then
!					derv = lag_deriv(r(n_radios-2),r(n_radios-1),r(n_radios),&
!						vel(n_radios-2),vel(n_radios-1),vel(n_radios),2)
					derv = deriv(r(n_radios-2),r(n_radios-1),r(n_radios),&
						vel(n_radios-2),vel(n_radios-1),vel(n_radios),2)						
				endif
			endif
		endif		

! Calculate line integrated opacity
		chi0 = chil(k) * perfil + kappa(k) + kappa_dust(k)
!		opac = sum(wtnu * chi0) / doppler_width(k,itr)  
		opac = chil(k) / doppler_width(k,itr)

! Calculate line integrated source function		
		s0 = ( chil(k) * Sl(k) * perfil + kappa(k) * B(k) + kappa_dust(k) * B_dust(k) ) / opac
!		source = sum(wtnu * s0)
		source = Sl(k)

! Integration of beta using gaussian quadrature with 24 points
		beta = 0.d0
		do dir = 1, gaussq_n 
			Q = dabs(vel(k)/r(k) * (1.d0 - gauss_abs(dir)**2) + gauss_abs(dir)**2 * derv)

			if (Q /= 0.d0) then
				tau_lvg = opac / Q
				integrand = (1.d0-exp(-tau_lvg))/tau_lvg
			else
				integrand = 0.d0				
			endif
			
			beta = beta + gauss_weights(dir) * integrand
		enddo
		
! Average intensity in Sobolev approximation
		Jbar(k) = (1.d0 - beta) * source + 0.5d0*(1.d0-mu_c) * beta * dentro(nfrq/2) +&
			beta * fuera(nfrq/2)
	enddo	!k

	
	
	deallocate(vel)
	deallocate(chi0)
	deallocate(s0)
	deallocate(gauss_weights)
	deallocate(gauss_abs)
		
	end subroutine LVGformal_solver		

!-----------------------------------------------------------------
! Calculates Jbar using LVG approximation and fills Jbar_total matrix
!-----------------------------------------------------------------

	subroutine calcJbarLVG
	integer :: it, ip, ipt
		do it = 1, nr

! Calculate opacities and source functions
			call opacity(it,n_radios,1,n_radios)

! Do the formal solution of the RTE using LVG approximation
			call LVGformal_solver(it)	

! Build up the Lstar_total and Jbar_total matrices. These are necessary for the construction
! of the rate matrix
			do ip = 1, n_radios
				ipt = nr*(ip-1)
				Lstar_total(it+ipt) = 0.d0 
				Jbar_total(it+ipt) = Jbar(ip)					
			enddo			

		enddo
	end subroutine calcJbarLVG

!-----------------------------------------------------------------
! Calcula las poblaciones en LVG
!-----------------------------------------------------------------

	subroutine lvgpops
	integer :: i, j
	real(kind=8), allocatable :: nh_ant(:)
	real(kind=8) :: factor, omega_bak
			
		omega_bak = omega_sor
		omega_sor = 1.d0
		print *, 'LVG solution to : ', dabs(lvg_flag)
		relat_error = 1.d0
		i = 0		
		
		if (lvg_flag < 0.d0) then				
			print *, 'Running LVG in slow mode'
			allocate(nh_ant(n_radios))
			nh_ant = nh
			factor = 1.d10
			
			do j = 1, 10
				nh = nh_ant * factor
				i = 0
				relat_error = 1.d0
				print *, 'Multiplying density by factor : ', factor
								
				do while (relat_error > dabs(lvg_flag) .and. i < 1000)
					i = i + 1
					call calcJbarLVG

					call corrige_poblaciones	
               print *, 'Relative error : ', relat_error				
				enddo				
				print *, 'Relative error : ', relat_error
				
				factor = 10.d0**(9.d0-j)

			enddo
			
			nh = nh_ant
			deallocate(nh_ant)		
			
		else
			do while (relat_error > dabs(lvg_flag) .and. i < 1000)
				i = i + 1

				call calcJbarLVG
				call corrige_poblaciones            

				print *, 'Step: ',i, '  Relat. error : ', relat_error

			enddo
		endif
		
! Initial population		
		popi(1:nl*n_radios) = pop(1:nl*n_radios)
		if (i == 1000) then
			print *, 'Reached 1000 iterations in LVG initialization. Trying these values'
		endif

		omega_sor = omega_bak
		
	end subroutine lvgpops


end module lvg_func
