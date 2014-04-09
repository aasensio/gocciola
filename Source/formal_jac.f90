
module formal_jacobi
use variables
use general
use short_ch
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
	subroutine formal_solver_jacobi
	integer :: it, ip, ipt, npmip1
		
		do it = 1, nr

! Calculate opacities and source functions
			call opacity(it,npmip1,1,n_radios)
			
			if (verbose_flag == 1) then
         	if (mod(it,10) == 0) then
					write(*,*) '       Transition : ', it
				endif
         endif

! Do the formal solution of the RTE
			call rtspher_jacobi(npmip1, it)	

! Build up the Lstar_total and Jbar_total matrices. These are necessary for the construction
! of the rate matrix
			do ip = 1, n_radios
				ipt = nr*(ip-1)
				Lstar_total(it+ipt) = lstar(ip)
				Jbar_total(it+ipt) = Jbar(ip)					
			enddo			
			
		enddo
		
	end subroutine formal_solver_jacobi

	
!-----------------------------------------------------------------------
! Formal solution of the RTE for Jacobi iterative scheme.
! INPUT:
!     - np : number of shells to be taken into account
!     - itr : transition number
!-----------------------------------------------------------------------	
	subroutine rtspher_jacobi(np, itr)
	integer, INTENT(IN) :: np, itr
	integer :: idir, dir, k, k0, k1, kdel, km, kp, frq, cortes, inout
	real(kind=8), dimension(nfrq) :: chim, chi0, chip, sm, s0, sp, dtm, dtp, exu, rlu
	real(kind=8) :: wtmu, sdelta, psim, psi0, psip, short, dm, dp
	
	kp = 0
	Jbar = 0.d0 ; Lstar = 0.d0 ; In = 0.d0 ; Jbar20(itr,:) = 0.d0

! Calculate outer boundary condition for this transition
	call calcula_fondo(itr)

!---------------------
! INCOMING SECTION
!---------------------
	idir = -1 ; k0 = np ; k1 = 1 ; kdel = -1 ; inout = 1

! Outer boundary condition
	do k = 1, n_caract
		In(:,k) = fuera
	enddo	
	cortes = n_cortes(k0)
	do dir = 1, cortes
		call ang_weights(k0,dir,wtmu)
		Jbar(k0) = Jbar(k0) + 0.5d0 * sum( wtmu * frqwt(:,dir,k0,inout) *&
			 profile(:,dir,k0,inout) * In(:,dir) )
		Jbar20(itr,k0) = Jbar20(itr,k0) + FOURSQRTTWO * &
			sum( (3.d0*mu(k0,dir)**2-1.d0) * wtmu * frqwt(:,dir,k0,inout) * profile(:,dir,k0,inout) * In(:,dir) )
	enddo	
	
! Go through all the shells inside out
	do k = k0+kdel,k1,kdel
	
! Calculate the number of intersections between the characteristics and the shell
		cortes = n_cortes(k)
		
		km = k - kdel				
			
! Go through the intersection points in each shell
		do dir = 1, cortes
		
! Opacities in point O and M of short-characteristics
			chim = chil(km) * profile(:,dir,km,inout) + kappa(km) + kappa_dust(km)						
			sm = ( chil(km) * profile(:,dir,km,inout) * Sl(km) + kappa(km) * B(km) + kappa_dust(km) * B_dust(km) ) / chim

			chi0 = chil(k) * profile(:,dir,k,inout) + kappa(k) + kappa_dust(k)			
			s0 = ( chil(k) * profile(:,dir,k,inout) * Sl(k) + kappa(k) * B(k) + kappa_dust(k) * B_dust(k) ) / chi0

			rlu = chil(k) * profile(:,dir,k,inout) / chi0
		
! Generate the angular integration weight for this point
			call ang_weights(k,dir,wtmu)
		
			dm = dabs( dsqrt( r(km)**2 - p(dir)**2 ) - dsqrt( r(k)**2 - p(dir)**2 ) )
			
! If we are not in the tanget point nor in the last shell, we can use parabolic SC
			if (dir /= cortes .and. k /= k1) then
			 	kp = k + kdel
				chip = chil(kp) * profile(:,dir,kp,inout) + kappa(kp) + kappa_dust(kp)
				sp = ( chil(kp) * profile(:,dir,kp,inout) * Sl(kp) + kappa(kp) * B(kp) + kappa_dust(kp) * B_dust(kp) )&
						 / chip						 
				dp = dabs( dsqrt( r(kp)**2 - p(dir)**2.d0 ) - dsqrt( r(k)**2 - p(dir)**2.d0 ) )
			else
! If not, we use linear SC
				chip = chim
				sp = sm
				dp = dm
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
				if (dir /= cortes .and. k /= k1) then					
					call par_sc(dtm(frq),dtp(frq),psim,psi0,psip)
					short = psim * sm(frq) + psi0 * s0(frq) + psip * sp(frq)
					In(frq,dir) = In(frq,dir) * exu(frq) + short					
					
					sdelta = chil(k) * profile(frq,dir,k,inout) / chi0(frq)
					short = psi0 * sdelta
					Lstar(k) = Lstar(k) + 0.5d0 * wtmu * frqwt(frq,dir,k,inout) *&
						 profile(frq,dir,k,inout) * short * rlu(frq)

				else
					if (dir == cortes) then
						call par_sc(dtm(frq),dtp(frq),psim,psi0,psip)
						short = psim * sm(frq) + psi0 * s0(frq) + psip * sp(frq)
						In(frq,dir) = In(frq,dir) * exu(frq) + short
					
						sdelta = chil(k) * profile(frq,dir,k,inout) / chi0(frq)
						short = psi0 * sdelta
						Lstar(k) = Lstar(k) + repeat*0.5d0 * wtmu * frqwt(frq,dir,k,inout) *&
							 profile(frq,dir,k,inout) * short * rlu(frq)
					else
! Linear SC	
						call lin_sc(dtm(frq),psim,psi0)
						short = psim * sm(frq) + psi0 * s0(frq)
						In(frq,dir) = In(frq,dir) * exu(frq) + short
					
						sdelta = chil(k) * profile(frq,dir,k,inout) / chi0(frq)
						short = psi0 * sdelta

						Lstar(k) = Lstar(k) + 0.5d0 * wtmu * frqwt(frq,dir,k,inout) *&
							 profile(frq,dir,k,inout) * short * rlu(frq)
					endif
			
				endif
			
			enddo !frq
			
! For each shell, calculate the contribution to Jbar of the intensities
			Jbar(k) = Jbar(k) + 0.5d0 * sum( wtmu * frqwt(:,dir,k,inout) *&
				 profile(:,dir,k,inout) * In(:,dir) )
			Jbar20(itr,k) = Jbar20(itr,k) + FOURSQRTTWO * &
				sum( (3.d0*mu(k,dir)**2-1.d0) * wtmu * frqwt(:,dir,k,inout) * profile(:,dir,k,inout) * In(:,dir) )
			
		enddo !dir
	
	enddo !k

!---------------------
! LOWER BOUNDARY CONDITION
!---------------------
	idir = 1 ; k0 = 1 ; k1 = np ; kdel = 1 ; inout = 2
	
! Intersections with the core shell
	cortes = n_cortes(k0) 

! Calculate inner boundary condition for this transition
	call calcula_centro(itr)
	
! Core boundary condition (only for absorbing and a full source)
! If vacuum boundary condition, then do nothing
	if (central_flag /= 1 .and. central_flag /= 3) then		
		do k = 1, cortes
			In(:,k) = dentro			
		enddo
	endif

! If central star is present (only for sources smaller than the central core)
	if (central_flag == 3) then
		do k = 1, impact_star_radius
			In(:,k) = dentro
		enddo
	endif

! Contribution of boundary condition to Jbar
	do dir = 1, cortes
		call ang_weights(k0,dir,wtmu)
		Jbar(k0) = Jbar(k0) + 0.5d0 * sum( wtmu * frqwt(:,dir,k0,inout) *&
			 profile(:,dir,k0,inout) * In(:,dir) )
		Jbar20(itr,k0) = Jbar20(itr,k0) + FOURSQRTTWO * &
			sum( (3.d0*mu(k0,dir)**2-1.d0) * wtmu * frqwt(:,dir,k0,inout) * profile(:,dir,k0,inout) * In(:,dir) )
	enddo

!---------------------
! OUTGOING SECTION		
!---------------------		

! Go through all the shell from core+1 to the exterior one
	do k = k0+kdel,k1,kdel
	
! Calculate the number of intersections between the characteristic and the shell
		cortes = n_cortes(k)

		if (k /= k1) then
! Because we are going out, last intersection with the shell is the tangent one. We do not account
! for it and we will do it later (only if we are not in the last shell)		
			cortes = cortes - 1
		endif		
		
		do dir = 1, cortes 

! If we are not in the last one, use parabolic SC
			if ( k /= k1 ) then
				kp = k + kdel
				chip = chil(kp) * profile(:,dir,kp,inout) + kappa(kp) + kappa_dust(kp)
				sp = ( chil(kp) * profile(:,dir,kp,inout) * Sl(kp) + kappa(kp) * B(kp) + kappa_dust(kp) * B_dust(kp) ) / chip
			endif

			chi0 = chil(k) * profile(:,dir,k,inout) + kappa(k) + kappa_dust(k)
			s0 = ( chil(k) * profile(:,dir,k,inout) * Sl(k) + kappa(k) * B(k) + kappa_dust(k) * B_dust(k) ) / chi0

			rlu = chil(k) * profile(:,dir,k,inout) / chi0
			
			call ang_weights(k,dir,wtmu)

! As we have substracted 1 to the value of cortes, there will always be an M point
			km = k - kdel
			chim = chil(km) * profile(:,dir,km,inout) + kappa(km) + kappa_dust(km)
			sm = ( chil(km) * profile(:,dir,km,inout) * Sl(km) + kappa(km) * B(km) + kappa_dust(km) * B_dust(km) ) / chim
			dm = dabs( dsqrt( r(km)**2 - p(dir)**2.d0 ) - dsqrt( r(k)**2 - p(dir)**2.d0 ) )				

! If we are not in the last shell, there will always be a P point
			if (k /= k1) then
				dp = dabs( dsqrt( r(kp)**2 - p(dir)**2 ) - dsqrt( r(k)**2 - p(dir)**2 ) )
			else
				dp = 0.d0
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
					Lstar(k) = Lstar(k) + 0.5d0 * wtmu * frqwt(frq,dir,k,inout) * &
						profile(frq,dir,k,inout) * short * rlu(frq)
				else

! Linear SC
					call lin_sc(dtm(frq),psim,psi0)
					short = psim * sm(frq) + psi0 * s0(frq)
					In(frq,dir) = In(frq,dir) * exu(frq) + short
					
					sdelta = chil(k) * profile(frq,dir,k,inout) / chi0(frq)
					short = psi0 * sdelta
					Lstar(k) = Lstar(k) + 0.5d0 * wtmu * frqwt(frq,dir,k,inout) * &
						profile(frq,dir,k,inout) * short * rlu(frq)
				endif			
			
			enddo !frq
			
! Calculate the contribution to Jbar of the intensity in each point 

			Jbar(k) = Jbar(k) + 0.5d0 * sum( wtmu * frqwt(:,dir,k,inout) * &
				profile(:,dir,k,inout) * In(:,dir) )
			Jbar20(itr,k) = Jbar20(itr,k) + FOURSQRTTWO * &
				sum( (3.d0*mu(k,dir)**2-1.d0) * wtmu * frqwt(:,dir,k,inout) * profile(:,dir,k,inout) * In(:,dir) )

		enddo	!dir	

! Now take into account that the direction loop was done until cortes-1 if we are in an internal
! shell, because we don't have to propagate the intensity at the tangent point yet
		if (k /= k1) then
			call ang_weights(k,cortes+1,wtmu)
			Jbar(k) = Jbar(k) + 0.5d0 * sum( wtmu * frqwt(:,cortes+1,k,inout) * &
				profile(:,cortes+1,k,inout) * ( In(:,cortes+1) ) )
			Jbar20(itr,k) = Jbar20(itr,k) + FOURSQRTTWO * &
				sum( (3.d0*mu(k,dir)**2-1.d0) * wtmu * frqwt(:,dir,k,inout) * profile(:,dir,k,inout) * In(:,dir) )
		endif
		
	enddo !k
	
	
	end subroutine rtspher_jacobi

!-----------------------------------------------------------------------
! Formal solution of the RTE for Jacobi iterative scheme.
! INPUT:
!     - np : number of shells to be taken into account
!     - itr : transition number
! It returns only the intensity (for the calculation of the flux)
!-----------------------------------------------------------------------	
	subroutine rtspher_jacobi_only_intensity(np, itr)
	integer, INTENT(IN) :: np, itr
	integer :: idir, dir, k, k0, k1, kdel, km, kp, frq, cortes, inout
	real(kind=8), dimension(nfrq) :: chim, chi0, chip, sm, s0, sp, dtm, dtp, exu
	real(kind=8) :: psim, psi0, psip, short, dm, dp
		
	kp = 0
	In = 0.d0

! Calculate outer boundary condition for this transition
	call calcula_fondo(itr)

!---------------------
! INCOMING SECTION
!---------------------
	idir = -1 ; k0 = np ; k1 = 1 ; kdel = -1 ; inout = 1

! Outer boundary condition
	do k = 1, n_caract
		In(:,k) = fuera
	enddo	

	cortes = n_cortes(k0)
	
! Go through all the shells inside out
	do k = k0+kdel,k1,kdel
	
! Calculate the number of intersections between the characteristics and the shell
		cortes = n_cortes(k)
		
		km = k - kdel				
			
! Go through the intersection points in each shell
		do dir = 1, cortes
		
! Opacities in point O and M of short-characteristics
			chim = chil(km) * profile(:,dir,km,inout) + kappa(km) + kappa_dust(km)						
			sm = ( chil(km) * profile(:,dir,km,inout) * Sl(km) + kappa(km) * B(km) + kappa_dust(km) * B_dust(km) ) / chim

			chi0 = chil(k) * profile(:,dir,k,inout) + kappa(k) + kappa_dust(k)			
			s0 = ( chil(k) * profile(:,dir,k,inout) * Sl(k) + kappa(k) * B(k) + kappa_dust(k) * B_dust(k) ) / chi0
		
			dm = dabs( dsqrt( r(km)**2 - p(dir)**2 ) - dsqrt( r(k)**2 - p(dir)**2 ) )
			
! If we are not in the tanget point nor in the last shell, we can use parabolic SC
			if (dir /= cortes .and. k /= k1) then
			 	kp = k + kdel
				chip = chil(kp) * profile(:,dir,kp,inout) + kappa(kp) + kappa_dust(kp)
				sp = ( chil(kp) * profile(:,dir,kp,inout) * Sl(kp) + kappa(kp) * B(kp) + kappa_dust(kp) * B_dust(kp) )&
						 / chip						 
				dp = dabs( dsqrt( r(kp)**2 - p(dir)**2.d0 ) - dsqrt( r(k)**2 - p(dir)**2.d0 ) )
			else
! If not, we use linear SC
				chip = chim
				sp = sm
				dp = dm
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
				if (dir /= cortes .and. k /= k1) then					
					call par_sc(dtm(frq),dtp(frq),psim,psi0,psip)
					short = psim * sm(frq) + psi0 * s0(frq) + psip * sp(frq)
					In(frq,dir) = In(frq,dir) * exu(frq) + short					
				else
					if (dir == cortes) then
						call par_sc(dtm(frq),dtp(frq),psim,psi0,psip)
						short = psim * sm(frq) + psi0 * s0(frq) + psip * sp(frq)
						In(frq,dir) = In(frq,dir) * exu(frq) + short
					else
! Linear SC	
						call lin_sc(dtm(frq),psim,psi0)
						short = psim * sm(frq) + psi0 * s0(frq)
						In(frq,dir) = In(frq,dir) * exu(frq) + short
					endif
			
				endif
			
			enddo !frq
			
			
		enddo !dir
	
	enddo !k

!---------------------
! LOWER BOUNDARY CONDITION
!---------------------
	idir = 1 ; k0 = 1 ; k1 = np ; kdel = 1 ; inout = 2
	
! Intersections with the core shell
	cortes = n_cortes(k0) 

! Calculate inner boundary condition for this transition
	call calcula_centro(itr)
	
! Core boundary condition (only for absorbing and a full source)
! If vacuum boundary condition, then do nothing
	if (central_flag /= 1 .and. central_flag /= 3) then		
		do k = 1, cortes-1
			In(:,k) = dentro			
		enddo
	endif

! If central star is present (only for sources smaller than the central core)
	if (central_flag == 3) then
		do k = 1, impact_star_radius
			In(:,k) = dentro
		enddo
	endif

!---------------------
! OUTGOING SECTION		
!---------------------		

! Go through all the shell from core+1 to the exterior one
	do k = k0+kdel,k1,kdel
	
! Calculate the number of intersections between the characteristic and the shell
		cortes = n_cortes(k)

		if (k /= k1) then
! Because we are going out, last intersection with the shell is the tangent one. We do not account
! for it and we will do it later (only if we are not in the last shell)		
			cortes = cortes - 1
		endif		
		
		do dir = 1, cortes 

! If we are not in the last one, use parabolic SC
			if ( k /= k1 ) then
				kp = k + kdel
				chip = chil(kp) * profile(:,dir,kp,inout) + kappa(kp) + kappa_dust(kp)
				sp = ( chil(kp) * profile(:,dir,kp,inout) * Sl(kp) + kappa(kp) * B(kp) + kappa_dust(kp) * B_dust(kp) ) / chip
			endif

			chi0 = chil(k) * profile(:,dir,k,inout) + kappa(k) + kappa_dust(k)
			s0 = ( chil(k) * profile(:,dir,k,inout) * Sl(k) + kappa(k) * B(k) + kappa_dust(k) * B_dust(k) ) / chi0

! As we have substracted 1 to the value of cortes, there will always be an M point
			km = k - kdel
			chim = chil(km) * profile(:,dir,km,inout) + kappa(km) + kappa_dust(km)
			sm = ( chil(km) * profile(:,dir,km,inout) * Sl(km) + kappa(km) * B(km) + kappa_dust(km) * B_dust(km) ) / chim
			dm = dabs( dsqrt( r(km)**2 - p(dir)**2.d0 ) - dsqrt( r(k)**2 - p(dir)**2.d0 ) )				

! If we are not in the last shell, there will always be a P point
			if (k /= k1) then
				dp = dabs( dsqrt( r(kp)**2 - p(dir)**2 ) - dsqrt( r(k)**2 - p(dir)**2 ) )
			else
				dp = 0.d0
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
	
	
	end subroutine rtspher_jacobi_only_intensity
end module formal_jacobi
