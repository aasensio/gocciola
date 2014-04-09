module formal_jacobi_pp
use variables
use general
use short_ch
use maths
implicit none
contains

! *********************************************************
! *********************************************************
! JACOBI FORMAL SOLUTION ROUTINES
! *********************************************************
! *********************************************************

!-----------------------------------------------------------------
! Calculate Jbar_total and Lstar_total for this model and atmosphere
!-----------------------------------------------------------------
	subroutine formal_solver_jacobi_pp
	integer :: it, ip, ipt, npmip1
	
		do it = 1, nr

! Calculate opacities and source functions
			call opacity(it,npmip1,1,n_radios)
			
			if (verbose_flag == 1) then
				if (mod(it,1) == 0) then
					write(*,*) '       Transition (bb) : ', it
				endif
         endif

! Do the formal solution of the RTE
			call rtpp_jacobi(npmip1, it)	
				
		
! Build up the Lstar_total and Jbar_total matrices. These are necessary for the construction
! of the rate matrix			
			do ip = 1, n_radios
				ipt = nr*(ip-1)
				Lstar_total(it+ipt) = lstar(ip)
				Jbar_total(it+ipt) = Jbar(ip)
			enddo
			
		enddo
								
		if (include_bound_free == 1) then
			do it = 1, nbf

				if (verbose_flag == 1) then
					if (mod(it,1) == 0) then
						write(*,*) '       Transition (bf) : ', it+nr
					endif
         	endif

! Do the formal solution of the RTE
				call rtpp_jacobi_bf(npmip1, it+nr)

! Build up the Lstar_total and Jbar_total matrices. These are necessary for the construction
! of the rate matrix
				do ip = 1, n_radios
					ipt = nbf*(ip-1)
					Jbar_total_bf_a(it+ipt) = Jbar(ip)
					Jbar_total_bf_b(it+ipt) = Jbarb(ip)

					Lstar_total_bf_a(it+ipt) = Lstar(ip)
					Lstar_total_bf_b(it+ipt) = Lstarb(ip)
					Lstar_total_bf_c(it+ipt) = Lstarc(ip)
					I_total_bf(it+ipt) = Ibf(ip)
				enddo			

			enddo
		endif
		
	end subroutine formal_solver_jacobi_pp

	
!-----------------------------------------------------------------------
! Formal solution of the RTE for Jacobi iterative scheme in plane-parallel geometry
! INPUT:
!     - np : number of shells to be taken into account
!     - itr : transition number
!-----------------------------------------------------------------------	
	subroutine rtpp_jacobi(np, itr)
	integer, INTENT(IN) :: np, itr
	integer :: idir, dir, k, k0, k1, kdel, km, kp, frq, cortes, inout
	real(kind=8), dimension(nfrq) :: chim, chi0, chip, sm, s0, sp, dtm, dtp, exu, rlu
	real(kind=8) :: sdelta, psim, psi0, psip, short, dm, dp
	
	Jbar = 0.d0 ; Lstar = 0.d0 ; In = 0.d0 ; Jbar20(itr,:) = 0.d0


! Set angle quadrature
!	select case(angleset_pp)
!		case(1)
!			mu_pp(1) = 1.d0 / dsqrt(3.d0)
!			wtmu(1) = 1.d0
!		case(2)
!			mu_pp(1) = 0.80847425d0
!			mu_pp(2) = 0.57735027d0
!			mu_pp(3) = 0.11417547d0

!			wtmu(1) = 2.d0 / 6.d0
!			wtmu(2) = 2.d0 / 6.d0
!			wtmu(3) = 2.d0 / 6.d0
!		case(3)
!			mu_pp(1) = 0.9324695142d0
!			mu_pp(2) = 0.6612093864d0
!			mu_pp(3) = 0.2386191860d0

!			wtmu(1) = 0.1713244923d0
!			wtmu(2) = 0.3607615730d0
!			wtmu(3) = 0.4679139345d0

!		case(4)
!			mu_pp(1) = 0.9602898564d0
!			mu_pp(2) = 0.7966664774d0
!			mu_pp(3) = 0.5255324099d0
!			mu_pp(4) = 0.1834346424d0

!			wtmu(1) = 0.1012285362d0
!			wtmu(2) = 0.2223810344d0
!			wtmu(3) = 0.3137066458d0
!			wtmu(4) = 0.3626837833d0

!	end select
	
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
		Jbar20(itr,k0) = Jbar20(itr,k0) + FOURSQRTTWO * &
				sum( (3.d0*mu_pp(dir)**2-1.d0) * wtmu(dir) * frqwt(:,dir,k0,inout) * profile(:,dir,k0,inout) * In(:,dir) )
	enddo	
		
! Go through all the shells inside out
	do k = k0+kdel,k1,kdel
	
! Calculate the number of intersections between the characteristics and the shell
		cortes = angleset_pp
		
		km = k - kdel				
			
! Go through the intersection points in each shell
		do dir = 1, cortes
			
! Opacities in point O and M of short-characteristics			
			chim = chil(km) * profile(:,dir,km,inout) + kappa(km) + kappa_dust(km)						
			sm = ( chil(km) * profile(:,dir,km,inout) * Sl(km) + kappa(km) * B(km) + &
				kappa_dust(km) * B_dust(km) ) / chim
			dm = dabs(1.d0/mu_pp(dir) * (r(km) - r(k)) )			

			chi0 = chil(k) * profile(:,dir,k,inout) + kappa(k) + kappa_dust(k)
			
			s0 = ( chil(k) * profile(:,dir,k,inout) * Sl(k) + kappa(k) * B(k) + kappa_dust(k) * B_dust(k) ) / chi0
			
			rlu = chil(k) * profile(:,dir,k,inout) / chi0

			if (k /= k1) then
				kp = k + kdel
				chip = chil(kp) * profile(:,dir,kp,inout) + kappa(kp) + kappa_dust(kp)
				sp = ( chil(kp) * profile(:,dir,kp,inout) * Sl(kp) + kappa(kp) * B(kp) + kappa_dust(kp) * B_dust(kp) )&
							 / chip						 
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
						 profile(frq,dir,k,inout) * short * rlu(frq)
				else
					
! Linear SC	
					call lin_sc(dtm(frq),psim,psi0)
					short = psim * sm(frq) + psi0 * s0(frq)
					In(frq,dir) = In(frq,dir) * exu(frq) + short
					
					sdelta = chil(k) * profile(frq,dir,k,inout) / chi0(frq)
					short = psi0 * sdelta
					Lstar(k) = Lstar(k) + 0.5d0 * wtmu(dir) * frqwt(frq,dir,k,inout) *&
						 profile(frq,dir,k,inout) * short * rlu(frq)
			
				endif
			
			enddo !frq
			
! For each shell, calculate the contribution to Jbar of the intensities
			Jbar(k) = Jbar(k) + 0.5d0 * sum( wtmu(dir) * frqwt(:,dir,k,inout) *&
				 profile(:,dir,k,inout) * In(:,dir) )
			Jbar20(itr,k) = Jbar20(itr,k) + FOURSQRTTWO * &
				sum( (3.d0*mu_pp(dir)**2-1.d0) * wtmu(dir) * frqwt(:,dir,k,inout) * profile(:,dir,k,inout) * In(:,dir) )
			
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
		Jbar20(itr,k0) = Jbar20(itr,k0) + FOURSQRTTWO * &
			sum( (3.d0*mu_pp(dir)**2-1.d0) * wtmu(dir) * frqwt(:,dir,k0,inout) * profile(:,dir,k0,inout) * In(:,dir) )
	enddo


!---------------------
! OUTGOING SECTION		
!---------------------

! Go through all the shell from core+1 to the exterior one
	do k = k0+kdel,k1,kdel
	
! Calculate the number of intersections between the characteristic and the shell
		cortes = angleset_pp
		
		do dir = 1, cortes 

			chi0 = chil(k) * profile(:,dir,k,inout) + kappa(k) + kappa_dust(k)
			s0 = ( chil(k) * profile(:,dir,k,inout) * Sl(k) + kappa(k) * B(k) + kappa_dust(k) * B_dust(k) ) / chi0

			rlu = chil(k) * profile(:,dir,k,inout) / chi0
			
			km = k - kdel
			chim = chil(km) * profile(:,dir,km,inout) + kappa(km) + kappa_dust(km)
			sm = ( chil(km) * profile(:,dir,km,inout) * Sl(km) + kappa(km) * B(km) + kappa_dust(km) * B_dust(km) ) / chim
			dm = dabs(1.d0/mu_pp(dir) * (r(km) - r(k)))

! If we are not in the last one, use parabolic SC
			if ( k /= k1 ) then
				kp = k + kdel
				chip = chil(kp) * profile(:,dir,kp,inout) + kappa(kp) + kappa_dust(kp)
				sp = ( chil(kp) * profile(:,dir,kp,inout) * Sl(kp) + kappa(kp) * B(kp) + kappa_dust(kp) * B_dust(kp) ) / chip
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
					
					sdelta = 1.d0 !chil(k) * profile(frq,dir,k,inout) / chi0(frq)
					short = psi0 * sdelta
					Lstar(k) = Lstar(k) + 0.5d0 * wtmu(dir) * frqwt(frq,dir,k,inout) * &
						profile(frq,dir,k,inout) * short * rlu(frq)
				else

! Linear SC
					call lin_sc(dtm(frq),psim,psi0)
					short = psim * sm(frq) + psi0 * s0(frq)
					In(frq,dir) = In(frq,dir) * exu(frq) + short
					
					sdelta = 1.d0 !chil(k) * profile(frq,dir,k,inout) / chi0(frq)
					short = psi0 * sdelta
					Lstar(k) = Lstar(k) + 0.5d0 * wtmu(dir) * frqwt(frq,dir,k,inout) * &
						profile(frq,dir,k,inout) * short * rlu(frq)
					
				endif			
			
			enddo !frq
			
! Calculate the contribution to Jbar of the intensity in each point 

			Jbar(k) = Jbar(k) + 0.5d0 * sum( wtmu(dir) * frqwt(:,dir,k,inout) * &
				profile(:,dir,k,inout) * In(:,dir) )

			Jbar20(itr,k) = Jbar20(itr,k) + FOURSQRTTWO * &
				sum( (3.d0*mu_pp(dir)**2-1.d0) * wtmu(dir) * frqwt(:,dir,k,inout) * profile(:,dir,k,inout) * In(:,dir) )

		enddo	!dir	
		
	enddo !k

		
	end subroutine rtpp_jacobi
	
!-----------------------------------------------------------------------
! Formal solution of the RTE for Jacobi iterative scheme in plane-parallel geometry
! INPUT:
!     - np : number of shells to be taken into account
!     - itr : transition number
!-----------------------------------------------------------------------	
	subroutine rtpp_jacobi_bf(np, itr)
	integer, INTENT(IN) :: np, itr
	integer :: idir, dir, k, k0, k1, kdel, km, kp, frq, cortes, inout, up, low, ipl
	real(kind=8), dimension(nfrq_bf) :: chim, chi0, chip, sm, s0, sp, dtm, dtp, exu, rlu, bf_weight, bf_cross_sec, bf_freq
	real(kind=8), dimension(nfrq_bf) :: bf_boltzmann, bf_twohnu3, Slu
	real(kind=8) :: In_bf(nfrq,n_caract)
	real(kind=8) :: sdelta, psim, psi0, psip, short, dm, dp
	real(kind=8) :: nlow_m, nlow_0, nlow_p, nup_m, nup_0, nup_p
	real(kind=8) :: ratio_lte_m, ratio_lte_0, ratio_lte_p
	
	Jbar = 0.d0
	Lstar = 0.d0
	In_bf = 0.d0
	Jbarb = 0.d0
	Lstarb = 0.d0
	Lstarc = 0.d0
	Ibf = 0.d0
	
! Calculate outer boundary condition for this transition
	call calcula_fondo_bf(itr)

! Calculate inner boundary condition for this transition
	call calcula_centro_bf(itr)

!---------------------
! INCOMING SECTION
!---------------------
	idir = -1 ; k0 = np ; k1 = 1 ; kdel = -1 ; inout = 1 

! Outer boundary condition
! n_caract pasa a valer 1, si solo usamos un angulo
	do k = 1, n_caract
		In_bf(:,k) = fuera_bf
	enddo	
	
	cortes = angleset_pp
	
	bf_freq = bf_cross_section(1,itr,:)
	bf_cross_sec = bf_cross_section(2,itr,:)
	bf_weight = bf_cross_section(3,itr,:)
	bf_boltzmann = exp(-PH*bf_freq / (PK * datph(1,k0)))
	bf_twohnu3 = 2.d0*PH*bf_freq**3 / PC**2
		
	Ibf(k0) = 4.d0 * PI * sum(bf_weight * bf_cross_sec * bf_boltzmann * bf_twohnu3)
	
	do dir = 1, cortes
		Jbar(k0) = Jbar(k0) + 4.d0 * PI * sum( wtmu(dir) * bf_weight * bf_cross_sec * In_bf(:,dir) )
		Jbarb(k0) = Jbarb(k0) + 4.d0 * PI * sum( wtmu(dir) * bf_weight * bf_cross_sec * bf_boltzmann * In_bf(:,dir) )		
	enddo
	
	up = itran(1,itr)
	low = itran(2,itr)
	
! Go through all the shells inside out
	do k = k0+kdel,k1,kdel
	
		bf_boltzmann = exp(-PH*bf_freq / (PK * datph(1,k)))
				
! Calculate the number of intersections between the characteristics and the shell
		cortes = angleset_pp
		
		km = k - kdel
											
! Go through the intersection points in each shell
		do dir = 1, cortes
			
! Opacities in point O and M of short-characteristics

! Point M
			ipl = nl * (km-1)
			nlow_m = pop(low+ipl)
			nup_m = pop(up+ipl)
			ratio_lte_m = popl(low+ipl) / popl(up+ipl)
							
			chim = (nlow_m - nup_m * ratio_lte_m * exp(-PH*bf_freq / (PK * datph(1,km)))) * bf_cross_sec + &
				chic_bf(itr,:,km) + kappa_dust(km)
! Emissivity
			sm = bf_twohnu3 * nup_m * ratio_lte_m * exp(-PH*bf_freq / (PK * datph(1,km))) * bf_cross_sec
			sm = ( sm + chic_bf(itr,:,km) * Sc_bf(itr,:,km) + kappa_dust(km) * B_dust(km) ) / chim
			dm = dabs(1.d0/mu_pp(dir) * (r(km) - r(k)) )			
			
! Point O
			ipl = nl * (k-1)
			nlow_0 = pop(low+ipl)
			nup_0 = pop(up+ipl)
			ratio_lte_0 = popl(low+ipl) / popl(up+ipl)

			chi0 = (nlow_0 - nup_0 * ratio_lte_0 * exp(-PH*bf_freq / (PK * datph(1,k)))) * bf_cross_sec + &
				chic_bf(itr,:,k) + kappa_dust(k)
			s0 = bf_twohnu3 * nup_0 * ratio_lte_0 * exp(-PH*bf_freq / (PK * datph(1,k))) * bf_cross_sec
			s0 = ( s0 + chic_bf(itr,:,k) * Sc_bf(itr,:,k) + kappa_dust(k) * B_dust(k) ) / chi0
			
			Slu = nup_0 * ratio_lte_0 * exp(-PH*bf_freq / (PK * datph(1,k)))						
			Slu = bf_twohnu3 * Slu / ( nlow_0 - Slu)
												
			rlu = (nlow_0 - nup_0 * ratio_lte_0 * exp(-PH*bf_freq / (PK * datph(1,k)))) * bf_cross_sec / chi0

! Point P
			if (k /= k1) then
			
				kp = k + kdel
				
				ipl = nl * (kp-1)
				nlow_p = pop(low+ipl)
				nup_p = pop(up+ipl)
				ratio_lte_p = popl(low+ipl) / popl(up+ipl)
								
				chip = (nlow_p - nup_p * ratio_lte_p * exp(-PH*bf_freq / (PK * datph(1,kp)))) * bf_cross_sec + &
					chic_bf(itr,:,kp) + kappa_dust(kp)
				sp = bf_twohnu3 * nup_p * ratio_lte_p * exp(-PH*bf_freq / (PK * datph(1,kp))) * bf_cross_sec
				sp = ( sp + chic_bf(itr,:,kp) * Sc_bf(itr,:,kp) + kappa_dust(kp) * B_dust(kp) ) / chip
				
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
			
			do frq = 1, nfrq_bf
											
! Parabolic SC
				if (k /= k1) then					
					call par_sc(dtm(frq),dtp(frq),psim,psi0,psip)
					short = psim * sm(frq) + psi0 * s0(frq) + psip * sp(frq)
					In_bf(frq,dir) = In_bf(frq,dir) * exu(frq) + short					
					
					sdelta = 1.d0 !chil(k) * profile(frq,dir,k,inout) / chi0(frq)
					short = psi0 * sdelta
					Lstar(k) = Lstar(k) + wtmu(dir) * bf_weight(frq) * bf_cross_sec(frq) * short * rlu(frq) * Slu(frq)
					Lstarb(k) = Lstarb(k) + wtmu(dir) * bf_weight(frq) * bf_cross_sec(frq) * bf_boltzmann(frq) * short * &
						rlu(frq) * Slu(frq)
					Lstarc(k) = Lstarc(k) + wtmu(dir) * bf_weight(frq) * bf_cross_sec(frq) * bf_boltzmann(frq) * &
						bf_twohnu3(frq) * short * rlu(frq)
				else
					
! Linear SC	
					call lin_sc(dtm(frq),psim,psi0)
					short = psim * sm(frq) + psi0 * s0(frq)
					In_bf(frq,dir) = In_bf(frq,dir) * exu(frq) + short
					
					sdelta = 1.d0 !chil(k) * profile(frq,dir,k,inout) / chi0(frq)
					short = psi0 * sdelta
					Lstar(k) = Lstar(k) + wtmu(dir) * bf_weight(frq) * bf_cross_sec(frq) * short * rlu(frq) * Slu(frq)
					Lstarb(k) = Lstarb(k) + wtmu(dir) * bf_weight(frq) * bf_cross_sec(frq) * bf_boltzmann(frq) * short * &
						rlu(frq) * Slu(frq)
					Lstarc(k) = Lstarc(k) + wtmu(dir) * bf_weight(frq) * bf_cross_sec(frq) * bf_boltzmann(frq) * &
						bf_twohnu3(frq) * short * rlu(frq)
			
				endif
			
			enddo !frq
						
! For each shell, calculate the contribution to Jbar of the intensities
			Jbar(k) = Jbar(k) + 4.d0 * PI * sum( wtmu(dir) * bf_weight * bf_cross_sec * In_bf(:,dir) )
			Jbarb(k) = Jbarb(k) + 4.d0 * PI * sum( wtmu(dir) * bf_weight * bf_cross_sec * bf_boltzmann * In_bf(:,dir) )
			
		enddo !dir
		
		Ibf(k) = 4.d0 * PI * sum(bf_weight * bf_cross_sec * bf_boltzmann * bf_twohnu3)
	
	enddo !k
		
!---------------------
! LOWER BOUNDARY CONDITION
!---------------------
	idir = 1 ; k0 = 1 ; k1 = np ; kdel = 1 ; inout = 2
	
! Intersections with the core shell
	cortes = angleset_pp

! Core boundary condition
	if (central_flag /= 1) then
		do k = 1, cortes
			In_bf(:,k) = dentro_bf
		enddo
	endif
	
	bf_boltzmann = exp(-PH*bf_freq / (PK * datph(1,k0)))
		
	do dir = 1, cortes
		Jbar(k0) = Jbar(k0) + 4.d0 * PI * sum( wtmu(dir) * bf_weight * bf_cross_sec * In_bf(:,dir) )
		Jbarb(k0) = Jbarb(k0) + 4.d0 * PI * sum( wtmu(dir) * bf_weight * bf_cross_sec * bf_boltzmann * In_bf(:,dir) )		
	enddo
	
!---------------------
! OUTGOING SECTION		
!---------------------

! Go through all the shell from core+1 to the exterior one
	do k = k0+kdel,k1,kdel
	
		bf_boltzmann = exp(-PH*bf_freq / (PK * datph(1,k)))
		
		km = k - kdel
	
! Calculate the number of intersections between the characteristic and the shell
		cortes = angleset_pp
		
		do dir = 1, cortes 

! Point M
			ipl = nl * (km-1)
			nlow_m = pop(low+ipl)
			nup_m = pop(up+ipl)
			ratio_lte_m = popl(low+ipl) / popl(up+ipl)
				
			chim = (nlow_m - nup_m * ratio_lte_m * exp(-PH*bf_freq / (PK * datph(1,km)))) * bf_cross_sec + &
				chic_bf(itr,:,km) + kappa_dust(km)
			sm = bf_twohnu3 * nup_m * ratio_lte_m * exp(-PH*bf_freq / (PK * datph(1,km))) * bf_cross_sec
			sm = ( sm + chic_bf(itr,:,km) * Sc_bf(itr,:,km) + kappa_dust(km) * B_dust(km) ) / chim
			dm = dabs(1.d0/mu_pp(dir) * (r(km) - r(k)) )			

! Point O
			ipl = nl * (k-1)
			nlow_0 = pop(low+ipl)
			nup_0 = pop(up+ipl)
			ratio_lte_0 = popl(low+ipl) / popl(up+ipl)

			chi0 = (nlow_0 - nup_0 * ratio_lte_0 * exp(-PH*bf_freq / (PK * datph(1,k)))) * bf_cross_sec + &
				chic_bf(itr,:,k) + kappa_dust(k)
			s0 = bf_twohnu3 * nup_0 * ratio_lte_0 * exp(-PH*bf_freq / (PK * datph(1,k))) * bf_cross_sec
			s0 = ( s0 + chic_bf(itr,:,k) * Sc_bf(itr,:,k) + kappa_dust(k) * B_dust(k) ) / chi0
			
			Slu = nup_0 * ratio_lte_0 * exp(-PH*bf_freq / (PK * datph(1,k)))			
			Slu = bf_twohnu3 * Slu / ( nlow_0 - Slu)
						
			rlu = (nlow_0 - nup_0 * ratio_lte_0 * exp(-PH*bf_freq / (PK * datph(1,k)))) * bf_cross_sec / chi0

! Point P
			if (k /= k1) then
				ipl = nl * (km-1)
				nlow_p = pop(low+ipl)
				nup_p = pop(up+ipl)
				ratio_lte_p = popl(low+ipl) / popl(up+ipl)

				kp = k + kdel
				
				chip = (nlow_p - nup_p * ratio_lte_p * exp(-PH*bf_freq / (PK * datph(1,kp)))) * bf_cross_sec + &
					chic_bf(itr,:,kp) + kappa_dust(kp)
				sp = bf_twohnu3 * nup_p * ratio_lte_p * exp(-PH*bf_freq / (PK * datph(1,kp))) * bf_cross_sec
				sp = ( sp + chic_bf(itr,:,kp) * Sc_bf(itr,:,kp) + kappa_dust(kp) * B_dust(kp) ) / chip
				
				dp = dabs(1.d0/mu_pp(dir) * (r(k) - r(kp)) )
			else
				dp = dm
				sp = sm
				chip = chim
			endif
			
			dtm = 0.5d0 * (chi0 + chim) * dm
			dtp = 0.5d0 * (chi0 + chip) * dp

! Calculate reduction factor
			where (dtm > trick)
				exu = dexp(-dtm)
			elsewhere
				exu = 1.d0 - dtm + dtm**2.d0 / 2.d0
			endwhere
			
			do frq = 1, nfrq_bf
				
! Parabolic SC
				if (k /= k1) then
	
					call par_sc(dtm(frq),dtp(frq),psim,psi0,psip)
					short = psim * sm(frq) + psi0 * s0(frq) + psip * sp(frq)
					In_bf(frq,dir) = In_bf(frq,dir) * exu(frq) + short
					
					sdelta = 1.d0 !chil(k) * profile(frq,dir,k,inout) / chi0(frq)
					short = psi0 * sdelta
					Lstar(k) = Lstar(k) + wtmu(dir) * bf_weight(frq) * bf_cross_sec(frq) * short * rlu(frq) * Slu(frq)
					Lstarb(k) = Lstarb(k) + wtmu(dir) * bf_weight(frq) * bf_cross_sec(frq) * bf_boltzmann(frq) * short * &
						rlu(frq) * Slu(frq)
					Lstarc(k) = Lstarc(k) + wtmu(dir) * bf_weight(frq) * bf_cross_sec(frq) * bf_boltzmann(frq) * &
						bf_twohnu3(frq) * short * rlu(frq)
				else

! Linear SC
					call lin_sc(dtm(frq),psim,psi0)
					short = psim * sm(frq) + psi0 * s0(frq)
					In_bf(frq,dir) = In_bf(frq,dir) * exu(frq) + short
					
					sdelta = 1.d0 !chil(k) * profile(frq,dir,k,inout) / chi0(frq)
					short = psi0 * sdelta
					Lstar(k) = Lstar(k) + wtmu(dir) * bf_weight(frq) * bf_cross_sec(frq) * short * rlu(frq) * Slu(frq)
					Lstarb(k) = Lstarb(k) + wtmu(dir) * bf_weight(frq) * bf_cross_sec(frq) * bf_boltzmann(frq) * short * &
						rlu(frq) * Slu(frq)
					Lstarc(k) = Lstarc(k) + wtmu(dir) * bf_weight(frq) * bf_cross_sec(frq) * bf_boltzmann(frq) * &
						bf_twohnu3(frq) * short * rlu(frq)
					
				endif			
			
			enddo !frq
			
! Calculate the contribution to Jbar of the intensity in each point 

			Jbar(k) = Jbar(k) + 4.d0 * PI * sum( wtmu(dir) * bf_weight * bf_cross_sec * In_bf(:,dir) )
			Jbarb(k) = Jbarb(k) + 4.d0 * PI * sum( wtmu(dir) * bf_weight * bf_cross_sec * bf_boltzmann * In_bf(:,dir) )

		enddo	!dir	
		
	enddo !k

		
	end subroutine rtpp_jacobi_bf	

!-----------------------------------------------------------------------
! Formal solution of the RTE for Jacobi iterative scheme in plane-parallel geometry
! INPUT:
!     - np : number of shells to be taken into account
!     - itr : transition number
! It returns only the intensity (for the calculation of the flux)
!-----------------------------------------------------------------------	
	subroutine rtpp_jacobi_only_intensity(np, itr)
	integer, INTENT(IN) :: np, itr
	integer :: idir, dir, k, k0, k1, kdel, km, kp, frq, cortes, inout
	real(kind=8), dimension(nfrq) :: chim, chi0, chip, sm, s0, sp, dtm, dtp, exu
	real(kind=8) :: psim, psi0, psip, short, dm, dp

	
	In = 0.d0

! Set angle quadrature
!	select case(angleset_pp)
!		case(1)
!			mu_pp(1) = 1.d0 / dsqrt(3.d0)
!			wtmu(1) = 1.d0
!		case(3)
!			mu_pp(1) = 0.80847425d0
!			mu_pp(2) = 0.57735027d0
!			mu_pp(3) = 0.11417547d0
!
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

			chi0 = chil(k) * profile(:,dir,k,inout) + kappa(k) + kappa_dust(k)
			s0 = ( chil(k) * profile(:,dir,k,inout) * Sl(k) + kappa(k) * B(k) + kappa_dust(k) * B_dust(k) ) / chi0
			
			km = k - kdel
			chim = chil(km) * profile(:,dir,km,inout) + kappa(km) + kappa_dust(km)
			sm = ( chil(km) * profile(:,dir,km,inout) * Sl(km) + kappa(km) * B(km) + kappa_dust(km) * B_dust(km) ) / chim
			dm = dabs(1.d0/mu_pp(dir) * (r(km) - r(k)))

! If we are not in the last one, use parabolic SC
			if ( k /= k1 ) then
				kp = k + kdel
				chip = chil(kp) * profile(:,dir,kp,inout) + kappa(kp) + kappa_dust(kp)
				sp = ( chil(kp) * profile(:,dir,kp,inout) * Sl(kp) + kappa(kp) * B(kp) + kappa_dust(kp) * B_dust(kp) ) / chip
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
	
	
	end subroutine rtpp_jacobi_only_intensity
end module formal_jacobi_pp
