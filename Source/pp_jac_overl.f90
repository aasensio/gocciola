module formal_jacobi_pp_overlap
use variables
use general
use short_ch
use overlap
implicit none
contains

! *********************************************************
! *********************************************************
! JACOBI FORMAL SOLUTION ROUTINES INCLUDING OVERLAPPING
! *********************************************************
! *********************************************************

!-----------------------------------------------------------------
! Calculate Jbar_total and Lstar_total for this model and atmosphere
!-----------------------------------------------------------------
	subroutine formal_solver_jacobi_pp_overlap
	integer :: it, ip, ipt, npmip1
	
		do it = 1, nr

! Calculate opacities and source functions
			call opacity(it,npmip1,1,n_radios)
			
			print *, '       Transition : ', it

! Do the formal solution of the RTE
			call rtpp_jacobi_overlap(npmip1, it)	

! Build up the Lstar_total and Jbar_total matrices. These are necessary for the construction
! of the rate matrix
			do ip = 1, n_radios
				ipt = nr*(ip-1)
				Lstar_total(it+ipt) = lstar(ip)
				Jbar_total(it+ipt) = Jbar(ip)					
			enddo			
			
		enddo
		
	end subroutine formal_solver_jacobi_pp_overlap

	
!-----------------------------------------------------------------------
! Formal solution of the RTE for Jacobi iterative scheme in plane-parallel geometry
! INPUT:
!     - np : number of shells to be taken into account
!     - itr : transition number
!-----------------------------------------------------------------------	
	subroutine rtpp_jacobi_overlap(np, itr)
	integer, INTENT(IN) :: np, itr
	integer :: idir, dir, k, k0, k1, kdel, km, kp, frq, cortes, inout
	real(kind=8), dimension(nfrq) :: chim, chi0, chip, sm, s0, sp, dtm, dtp, exu
	real(kind=8), dimension(nfrq) :: chim_ov, emm_ov, chi0_ov, em0_ov, chip_ov, emp_ov
	real(kind=8) :: sdelta, psim, psi0, psip, short, dm, dp

	real(kind=8), allocatable :: mu_pp(:), wtmu(:)
	
	Jbar = 0.d0 ; Lstar = 0.d0 ; In = 0.d0 ; Jbar20(itr,:) = 0.d0

	chim_ov = 0.d0
	emm_ov = 0.d0
	chi0_ov = 0.d0
	em0_ov = 0.d0
	chip_ov = 0.d0
	emp_ov = 0.d0		

	allocate(mu_pp(angleset_pp))
	allocate(wtmu(angleset_pp))

! Set angle quadrature
	select case(angleset_pp)
		case(1)
			mu_pp(1) = 1.d0 / dsqrt(3.d0)
			wtmu(1) = 1.d0
		case(3)
			mu_pp(1) = 0.80847425d0
			mu_pp(2) = 0.57735027d0
			mu_pp(3) = 0.11417547d0

			wtmu(1) = 2.d0 / 6.d0
			wtmu(2) = 2.d0 / 6.d0
			wtmu(3) = 2.d0 / 6.d0
	end select
	
! Calculate outer boundary condition for this transition
	call calcula_fondo(itr)

! Calculate inner boundary condition for this transition
	call calcula_centro(itr)

!---------------------
! INCOMING SECTION
!---------------------
	idir = -1 ; k0 = np ; k1 = 1 ; kdel = -1 ; inout = 1 
!	mu_pp = -1.d0 / dsqrt(3.d0) ; wtmu = 1.d0

! Outer boundary condition
! n_caract pasa a valer 1, si solo usamos un angulo
	do k = 1, n_caract
		In(:,k) = fuera
	enddo	
	
	cortes = angleset_pp

	do dir = 1, cortes
		Jbar(k0) = Jbar(k0) + 0.5d0 * sum( wtmu(dir) * frqwt(:,dir,k0,inout) *&
			 profile(:,dir,k0,inout) * In(:,dir) )
!		Jbar20(itr,k0) = Jbar20(itr,k0) + FOURSQRTTWO * &
!				sum( (3.d0*mu(k0,dir)**2-1.d0) * wtmu * wtnu * perfil * In(:,dir) )
	enddo	
	
! Go through all the shells inside out
	do k = k0+kdel,k1,kdel
	
! Calculate the number of intersections between the characteristics and the shell
		cortes = angleset_pp
		
		km = k - kdel				
			
! Go through the intersection points in each shell
		do dir = 1, cortes

! If overlapping, then calculate the opacity of the lines which overlap the one we are in
			call overlap_opacity(itr,km,dir,inout,chim_ov,emm_ov)
			
! Opacities in point O and M of short-characteristics
			chim = chil(km) * profile(:,dir,km,inout) + chim_ov + kappa(km) + kappa_dust(km)						

			sm = ( chil(km) * profile(:,dir,km,inout) * Sl(km) + emm_ov + kappa(km) * B(km) + &
				kappa_dust(km) * B_dust(km) ) / chim
				
			
			dm = dabs(1.d0/mu_pp(dir) * (r(km) - r(k)) )			

			call overlap_opacity(itr,k,dir,inout,chi0_ov,em0_ov)
	
			chi0 = chil(k) * profile(:,dir,k,inout) + chi0_ov + kappa(k) + kappa_dust(k)			
			s0 = ( chil(k) * profile(:,dir,k,inout) * Sl(k) + em0_ov + kappa(k) * B(k) + &
				kappa_dust(k) * B_dust(k) ) / chi0			
			
			if (k /= k1) then
				kp = k + kdel
				
				call overlap_opacity(itr,kp,dir,inout,chip_ov,emp_ov)
				
				chip = chil(kp) * profile(:,dir,kp,inout) + chip_ov + kappa(kp) + kappa_dust(kp)
				sp = ( chil(kp) * profile(:,dir,kp,inout) * Sl(kp) + emp_ov + kappa(kp) * B(kp) + &
					kappa_dust(kp) * B_dust(kp) ) / chip						 
				dp = dabs(1.d0/mu_pp(dir) * (r(k) - r(kp)) )
			else
				dp = dm
				sp = sm
				chip = chim
			endif

! Calculate optical depths
			dtm = 0.5d0 * (chi0 + chim) * dm
			dtp = 0.5d0 * (chi0 + chip) * dp
			
! Calculate reduction factor
			where (dtm > trick)
				exu = dexp(-dtm)
			elsewhere
				exu = 1.d0 - dtm + dtm**2.d0 / 2.d0
			endwhere
			
			do frq = 1, nfrq
											
! Parabolic SC
				if (k /= k1) then					
					call par_sc(dtm(frq),dtp(frq),psim,psi0,psip)
					short = psim * sm(frq) + psi0 * s0(frq) + psip * sp(frq)
					In(frq,dir) = In(frq,dir) * exu(frq) + short					
					
					sdelta = chil(k) * profile(frq,dir,k,inout) / chi0(frq)
					short = psi0 * sdelta
					Lstar(k) = Lstar(k) + 0.5d0 * wtmu(dir) * frqwt(frq,dir,k,inout) *&
						 profile(frq,dir,k,inout) * short
				else
					
! Linear SC	
					call lin_sc(dtm(frq),psim,psi0)
					short = psim * sm(frq) + psi0 * s0(frq)
					In(frq,dir) = In(frq,dir) * exu(frq) + short
					
					sdelta = chil(k) * profile(frq,dir,k,inout) / chi0(frq)
					short = psi0 * sdelta
					Lstar(k) = Lstar(k) + 0.5d0 * wtmu(dir) * frqwt(frq,dir,k,inout) *&
						 profile(frq,dir,k,inout) * short
			
				endif
			
			enddo !frq
			
! For each shell, calculate the contribution to Jbar of the intensities
			Jbar(k) = Jbar(k) + 0.5d0 * sum( wtmu(dir) * frqwt(:,dir,k,inout) *&
				 profile(:,dir,k,inout) * In(:,dir) )
!			Jbar20(itr,k) = Jbar20(itr,k) + FOURSQRTTWO * &
!				sum( (3.d0*mu(k,dir)**2-1.d0) * wtmu * wtnu * perfil * In(:,dir) )
			
		enddo !dir
	
	enddo !k

!---------------------
! LOWER BOUNDARY CONDITION
!---------------------
	idir = 1 ; k0 = 1 ; k1 = np ; kdel = 1 ; inout = 2
!	mu_pp = 1.d0 / dsqrt(3.d0) ; wtmu = 1.d0
	
! Intersections with the core shell
	cortes = angleset_pp

! Core boundary condition
	if (central_flag /= 1) then
		do k = 1, cortes
			In(:,k) = dentro
		enddo
	endif

! Contribution of boundary condition to Jbar
	do dir = 1, cortes	
		Jbar(k0) = Jbar(k0) + 0.5d0 * sum( wtmu(dir) * frqwt(:,dir,k0,inout) *&
			 profile(:,dir,k0,inout) * In(:,dir) )
!		Jbar20(itr,k0) = Jbar20(itr,k0) + FOURSQRTTWO * &
!				sum( (3.d0*mu(k0,dir)**2-1.d0) * wtmu * wtnu * perfil * In(:,dir) )
	enddo

!---------------------
! OUTGOING SECTION		
!---------------------

! Go through all the shell from core+1 to the exterior one
	do k = k0+kdel,k1,kdel
	
! Calculate the number of intersections between the characteristic and the shell
		cortes = angleset_pp
		
		do dir = 1, cortes 

			call overlap_opacity(itr,k,dir,inout,chi0_ov,em0_ov)

			chi0 = chil(k) * profile(:,dir,k,inout) + chi0_ov + kappa(k) + kappa_dust(k)
			s0 = ( chil(k) * profile(:,dir,k,inout) * Sl(k) + em0_ov + kappa(k) * B(k) + &
				kappa_dust(k) * B_dust(k) ) / chi0
			
			km = k - kdel

			call overlap_opacity(itr,km,dir,inout,chim_ov,emm_ov)
	
			chim = chil(km) * profile(:,dir,km,inout) + chim_ov + kappa(km) + kappa_dust(km)
			sm = ( chil(km) * profile(:,dir,km,inout) * Sl(km) + emm_ov + kappa(km) * B(km) + &
				kappa_dust(km) * B_dust(km) ) / chim
			dm = dabs(1.d0/mu_pp(dir) * (r(km) - r(k)))

! If we are not in the last one, use parabolic SC
			if ( k /= k1 ) then
				kp = k + kdel
				
				call overlap_opacity(itr,kp,dir,inout,chip_ov,emp_ov)
				
				chip = chil(kp) * profile(:,dir,kp,inout) + chip_ov + kappa(kp) + kappa_dust(kp)
				sp = ( chil(kp) * profile(:,dir,kp,inout) * Sl(kp) + emp_ov + kappa(kp) * B(kp) + &
					kappa_dust(kp) * B_dust(kp) ) / chip
				dp = dabs(1.d0/mu_pp(dir) * (r(kp) - r(k)))
			else
				chip = chim
				sp = sm
				dp = dm
			endif
			
			dtm = 0.5d0 * (chi0 + chim) * dm
			dtp = 0.5d0 * (chi0 + chip) * dp

! Calculate reduction factor
			where (dtm > trick)
				exu = dexp(-dtm)
			elsewhere
				exu = 1.d0 - dtm + dtm**2.d0 / 2.d0
			endwhere
			
			do frq = 1, nfrq
				
! Parabolic SC
				if (k /= k1) then
	
					call par_sc(dtm(frq),dtp(frq),psim,psi0,psip)
					short = psim * sm(frq) + psi0 * s0(frq) + psip * sp(frq)
					In(frq,dir) = In(frq,dir) * exu(frq) + short
					
					sdelta = chil(k) * profile(frq,dir,k,inout) / chi0(frq)
					short = psi0 * sdelta
					Lstar(k) = Lstar(k) + 0.5d0 * wtmu(dir) * frqwt(frq,dir,k,inout) * &
						profile(frq,dir,k,inout) * short
				else

! Linear SC
					call lin_sc(dtm(frq),psim,psi0)
					short = psim * sm(frq) + psi0 * s0(frq)
					In(frq,dir) = In(frq,dir) * exu(frq) + short
					
					sdelta = chil(k) * profile(frq,dir,k,inout) / chi0(frq)
					short = psi0 * sdelta
					Lstar(k) = Lstar(k) + 0.5d0 * wtmu(dir) * frqwt(frq,dir,k,inout) * &
						profile(frq,dir,k,inout) * short
					
				endif			
			
			enddo !frq
			
! Calculate the contribution to Jbar of the intensity in each point 

			Jbar(k) = Jbar(k) + 0.5d0 * sum( wtmu(dir) * frqwt(:,dir,k,inout) * &
				profile(:,dir,k,inout) * In(:,dir) )

!			Jbar20(itr,k) = Jbar20(itr,k) + FOURSQRTTWO * &
!				sum( (3.d0*mu(k,dir)**2-1.d0) * wtmu * wtnu * perfil * In(:,dir) )			

		enddo	!dir	
		
	enddo !k

	deallocate(wtmu)
	deallocate(mu_pp)
	
	
	end subroutine rtpp_jacobi_overlap

!-----------------------------------------------------------------------
! Formal solution of the RTE for Jacobi iterative scheme in plane-parallel geometry
! INPUT:
!     - np : number of shells to be taken into account
!     - itr : transition number
! It returns only the intensity (for the calculation of the flux)
!-----------------------------------------------------------------------	
	subroutine rtpp_jacobi_overlap_only_int(np, itr)
	integer, INTENT(IN) :: np, itr
	integer :: idir, dir, k, k0, k1, kdel, km, kp, frq, cortes, inout
	real(kind=8), dimension(nfrq) :: chim, chi0, chip, sm, s0, sp, dtm, dtp, exu
	real(kind=8), dimension(nfrq) :: chim_ov, emm_ov, chi0_ov, em0_ov, chip_ov, emp_ov
	real(kind=8) :: psim, psi0, psip, short, dm, dp
	
	In = 0.d0

	chim_ov = 0.d0
	emm_ov = 0.d0
	chi0_ov = 0.d0
	em0_ov = 0.d0
	chip_ov = 0.d0
	emp_ov = 0.d0		

! Set angle quadrature
!	select case(angleset_pp)
!		case(1)
!			mu_pp(1) = 1.d0 / dsqrt(3.d0)
!			wtmu(1) = 1.d0
!		case(3)
!			mu_pp(1) = 0.80847425d0
!			mu_pp(2) = 0.57735027d0
!			mu_pp(3) = 0.11417547d0

!			wtmu(1) = 2.d0 / 6.d0
!			wtmu(2) = 2.d0 / 6.d0
!			wtmu(3) = 2.d0 / 6.d0
!	end select
	
! Calculate outer boundary condition for this transition
	call calcula_fondo(itr)

! Calculate inner boundary condition for this transition
	call calcula_centro(itr)

!---------------------
! LOWER BOUNDARY CONDITION
!---------------------
	idir = 1 ; k0 = 1 ; k1 = np ; kdel = 1 ; inout = 2
!	mu_pp = 1.d0 / dsqrt(3.d0) ; wtmu = 1.d0
	
! Intersections with the core shell
	cortes = angleset_pp

! Core boundary condition
	if (central_flag /= 1) then
		do k = 1, cortes
			In(:,k) = dentro
		enddo
	endif

!---------------------
! OUTGOING SECTION		
!---------------------

! Go through all the shell from core+1 to the exterior one
	do k = k0+kdel,k1,kdel
	
! Calculate the number of intersections between the characteristic and the shell
		cortes = angleset_pp
		
		do dir = 1, cortes 

			call overlap_opacity(itr,k,dir,inout,chi0_ov,em0_ov)

			chi0 = chil(k) * profile(:,dir,k,inout) + chi0_ov + kappa(k) + kappa_dust(k)
			s0 = ( chil(k) * profile(:,dir,k,inout) * Sl(k) + em0_ov + kappa(k) * B(k) + &
				kappa_dust(k) * B_dust(k) ) / chi0
			
			km = k - kdel

			call overlap_opacity(itr,km,dir,inout,chim_ov,emm_ov)
	
			chim = chil(km) * profile(:,dir,km,inout) + chim_ov + kappa(km) + kappa_dust(km)
			sm = ( chil(km) * profile(:,dir,km,inout) * Sl(km) + emm_ov + kappa(km) * B(km) + &
				kappa_dust(km) * B_dust(km) ) / chim
			dm = dabs(1.d0/mu_pp(dir) * (r(km) - r(k)))

! If we are not in the last one, use parabolic SC
			if ( k /= k1 ) then
				kp = k + kdel
				
				call overlap_opacity(itr,kp,dir,inout,chip_ov,emp_ov)
				
				chip = chil(kp) * profile(:,dir,kp,inout) + chip_ov + kappa(kp) + kappa_dust(kp)
				sp = ( chil(kp) * profile(:,dir,kp,inout) * Sl(kp) + emp_ov + kappa(kp) * B(kp) + &
					kappa_dust(kp) * B_dust(kp) ) / chip
				dp = dabs(1.d0/mu_pp(dir) * (r(kp) - r(k)))
			else
				chip = chim
				sp = sm
				dp = dm
			endif
			
			dtm = 0.5d0 * (chi0 + chim) * dm
			dtp = 0.5d0 * (chi0 + chip) * dp

! Calculate reduction factor
			where (dtm > trick)
				exu = dexp(-dtm)
			elsewhere
				exu = 1.d0 - dtm + dtm**2.d0 / 2.d0
			endwhere
			
			do frq = 1, nfrq
				
! Parabolic SC
				if (k /= k1) then
	
					call par_sc(dtm(frq),dtp(frq),psim,psi0,psip)
					short = psim * sm(frq) + psi0 * s0(frq) + psip * sp(frq)
					In(frq,dir) = In(frq,dir) * exu(frq) + short
				else

! Linear SC
					call lin_sc(dtm(frq),psim,psi0)
					short = psim * sm(frq) + psi0 * s0(frq)
					In(frq,dir) = In(frq,dir) * exu(frq) + short
				endif			
			
			enddo !frq
			
		enddo	!dir	
		
	enddo !k
	
	
	end subroutine rtpp_jacobi_overlap_only_int
end module formal_jacobi_pp_overlap
