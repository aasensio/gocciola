module formal_gs_pp
use variables
use general
use statis
use short_ch
implicit none
contains

! *********************************************************
! *********************************************************
! GAUSS-SEIDEL FORMAL SOLUTION ROUTINES
! *********************************************************
! *********************************************************

!-----------------------------------------------------------------
! Does a Gauss-Seidel or SOR iteration
!-----------------------------------------------------------------
	subroutine formal_solver_gs_pp(which)
	integer :: which
		
		relat_error_p = relat_error
		relat_error = 0.d0	
		true_error = 0.d0

! Make a formal solution of RTE using GS or SOR (population changes are made inside)		
		call rtpp_gs(n_radios,which)	

	end subroutine formal_solver_gs_pp

	
!-----------------------------------------------------------------------
! Formal solution of the RTE for GS or SOR iterative scheme. Population changes
! are made inside the formal solver as mentioned in Trujillo Bueno & Fabiani Bendicho (1995)
! INPUT:
!     - np : number of shells to be taken into account
! OUTPUT:
!     Populations are changed
!-----------------------------------------------------------------------	
	
	subroutine rtpp_gs(np, which)
	integer, INTENT(IN) :: np, which
	integer :: idir, dir, k, k0, k1, kdel, km, kp, frq, cortes, nip, it, ipt, inout
	real(kind=8), allocatable :: chim(:), chi0(:), chip(:), sm(:), s0(:), sp(:), dtm(:), dtp(:), exu(:), true(:)
	real(kind=8), allocatable :: Inte(:,:,:), IM(:,:,:), opam(:,:,:), opap(:,:,:), dtmold(:,:,:), smold(:,:,:), sfm(:,:,:), sfp(:,:,:)
	real(kind=8), allocatable :: dism(:), disp(:)
	real(kind=8), allocatable :: Lsta(:,:), Jba(:,:)
	real(kind=8) :: sdelta, psim, psi0, psip, short, dm, dp, bond, psim_rev, psi0_rev, psip_rev
	real(kind=8) :: short_rev
!	real(kind=8) :: repetir = 1.d0
	
	allocate(chim(nfrq))
	allocate(chi0(nfrq))
	allocate(chip(nfrq))
	allocate(sm(nfrq))
	allocate(s0(nfrq))
	allocate(sp(nfrq))
	allocate(dtm(nfrq))
	allocate(dtp(nfrq))
	allocate(exu(nfrq))
	allocate(true(nfrq))

	allocate(Inte(nfrq,n_caract,nact))
	allocate(IM(nfrq,n_caract,nact))
	allocate(opam(nfrq,n_caract,nact))
	allocate(opap(nfrq,n_caract,nact))
	allocate(dtmold(nfrq,n_caract,nact))
	allocate(smold(nfrq,n_caract,nact))
	allocate(sfm(nfrq,n_caract,nact))
	allocate(sfp(nfrq,n_caract,nact))

	allocate(dism(n_caract))
	allocate(disp(n_caract))

	allocate(Lsta(n_radios,nact))
	allocate(Jba(n_radios,nact))
	
	Jba = 0.d0 ; Lsta = 0.d0 ; Inte = 0.d0 ; IM = 0.d0; dtmold = 0.d0 ; smold = 0.d0 ; true = 0.d0

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

!---------------------
! INCOMING SECTION
!---------------------
	idir = -1 ; k0 = np ; k1 = 1 ; kdel = -1 ; inout = 1
!	mu_pp = -1.d0 / dsqrt(3.d0) ; wtmu = 1.d0

! Boundary condition in the outer shell	
	cortes = angleset_pp
	
	do it = 1, nr		
		do k = 1, cortes
			call calcula_fondo(it)
			Inte(:,k,it) = fuera
		enddo
	enddo

	
! We go through the shells outside in
	do k = k0+kdel,k1,kdel
	
		do it = 1, nr

! Calculate opacities and source functions (line, background and dust)
! CHECK THE INTERVAL !!!!! I verify inside opacity that k+kdel > 0
			call opacity(it,nip,k+kdel,k-kdel)			
			
! Calculate the number of intersections between the characteristics and the shell
			cortes = angleset_pp

			km = k - kdel

! Go through the intersection points in each shell
			do dir = 1, cortes

! Opacities in point O and M of short-characteristics
				chim = chil(km) * profile(:,dir,km,inout) + kappa(km) + kappa_dust(km)
				sm = ( chil(km) * profile(:,dir,km,inout) * Sl(km) + kappa(km) * B(km) + kappa_dust(km) * B_dust(km) ) / chim
				dm = dabs(1.d0/mu_pp(dir) * (r(k)-r(km)))
				
				chi0 = chil(k) * profile(:,dir,k,inout) + kappa(k) + kappa_dust(k)
				s0 = ( chil(k) * profile(:,dir,k,inout) * Sl(k) + kappa(k) * B(k) + kappa_dust(k) * B_dust(k) ) / chi0

! If we are not in the tanget point nor in the last shell, we can use parabolic SC
				if (k /= k1) then
			 		kp = k + kdel
					chip = chil(kp) * profile(:,dir,kp,inout) + kappa(kp) + kappa_dust(kp)
					sp = ( chil(kp) * profile(:,dir,kp,inout) * Sl(kp) + kappa(kp) * B(kp) + kappa_dust(kp) * B_dust(kp) )&
							 / chip
					dp = dabs(1.d0/mu_pp(dir) * (r(k)-r(kp)))
				else
! If not, we use linear SC
					chip = chim
					sp = sm
					dp = dm
				endif

! Calculate optical depths
				dtm = 0.5d0 * (chi0 + chim) * dm
				dtp = 0.5d0 * (chi0 + chip) * dp

				do frq = 1, nfrq

					if (dtm(frq) > trick) then
               	exu(frq) = dexp(-dtm(frq))
            	else
               	exu(frq) = 1.d0 - dtm(frq) + dtm(frq)**2 / 2.d0
            	endif				
					
					bond = Inte(frq,dir,it) * exu(frq)
! Parabolic SC
					if (k /= k1) then

						call par_sc(dtm(frq),dtp(frq),psim,psi0,psip)
						short = psim * sm(frq) + psi0 * s0(frq) + psip * sp(frq)
						Inte(frq,dir,it) = bond + short
						Jba(k,it) = Jba(k,it) + 0.5d0 * wtmu(dir) * frqwt(frq,dir,k,inout) * &
							 profile(frq,dir,k,inout) * bond

!						sdelta = chil(k) * profile(frq,dir,k,inout) / chi0(frq)
!						short = psi0 * sdelta
!						Lstar(k) = Lstar(k) + 0.5 * wtmu * wtnu(frq) * perfil(frq) * short
					else

! Linear SC
						call lin_sc(dtm(frq),psim,psi0)
						short = psim * sm(frq) + psi0 * s0(frq)
						Inte(frq,dir,it) = bond + short
						Jba(k,it) = Jba(k,it) + 0.5d0 * wtmu(dir) * frqwt(frq,dir,k,inout) * &
							profile(frq,dir,k,inout) * Inte(frq,dir,it)

						sdelta = chil(k) * profile(frq,dir,k,inout) / chi0(frq)
						short = psi0 * sdelta
						Lsta(k,it) = Lsta(k,it) + 0.5d0 * wtmu(dir) * frqwt(frq,dir,k,inout) * &
							profile(frq,dir,k,inout) * short

					endif

				enddo !frq

			enddo !dir
			
		enddo !it
	
	enddo !k

!---------------------
! LOWER BOUNDARY CONDITION
!---------------------
	idir = 1 ; k0 = 1 ; k1 = np ; kdel = 1 ; inout = 2
!	mu_pp = 1.d0 / dsqrt(3.d0) ; wtmu = 1.d0
	
! Intersection with the core shell
	cortes = angleset_pp

! Boundary condition in the core
	if (central_flag /= 1) then
		do it = 1, nr
			do k = 1, cortes
				call calcula_centro(it)
				Inte(:,k,it) = dentro
			enddo
		enddo
	endif
	
! Contribution to Jbar(core) of the boundary condition
	do it = 1, nr
		do dir = 1, cortes
			Jba(k0,it) = Jba(k0,it) + 0.5d0 * sum( wtmu(dir) * frqwt(:,dir,k0,inout) * &
				profile(:,dir,k0,inout) * Inte(:,dir,it) )
		enddo
	enddo

!---------------------
! MAKE POPULATION CORRECTION
!---------------------

! Build up Lstar_total and Jbar_total matrices
	ipt = 0
	do it = 1, nr
		Lstar_total(it+ipt) = Lsta(k0,it)
		Jbar_total(it+ipt) = Jba(k0,it)
	enddo !it

! Do the correction
	call rate_eq(k0, relat_error, true_error)
		
!---------------------
! OUTGOING SECTION		
!---------------------		

! Go through all the shells, from core+1 to the outer shell
	do k = k0+kdel,k1,kdel
		
		do it = 1, nr
			
! Calculate opacities and source functions
			call opacity(it,nip,k-kdel,k+kdel)		  ! This could be from km to kp
			
! Calculate number of intersections between the characteristic and the shell
			cortes = angleset_pp

			do dir = 1, cortes 

! If we are not in the last shell, use parabolic SC
				

				chi0 = chil(k) * profile(:,dir,k,inout) + kappa(k) + kappa_dust(k)
				s0 = ( chil(k) * profile(:,dir,k,inout) * Sl(k) + kappa(k) * B(k) + kappa_dust(k) * B_dust(k) ) / chi0

! As we have substracted 1 to the value of cortes, there will always be an M point
				km = k - kdel
				chim = chil(km) * profile(:,dir,km,inout) + kappa(km) + kappa_dust(km)
				sm = ( chil(km) * profile(:,dir,km,inout) * Sl(km) + kappa(km) * B(km) + kappa_dust(km) * B_dust(km) ) / chim
				dm = dabs(1.d0/mu_pp(dir) * (r(km) - r(k)))				


				if ( k /= k1 ) then
					kp = k + kdel
					chip = chil(kp) * profile(:,dir,kp,inout) + kappa(kp) + kappa_dust(kp)
					sp = ( chil(kp) * profile(:,dir,kp,inout) * Sl(kp) + kappa(kp) * B(kp) + kappa_dust(kp) * B_dust(kp) ) / chip
					dp = dabs(1.d0/mu_pp(dir) * (r(kp) - r(k)))
				else
					dp = dm
					sp = sm
					chip = chim
				endif
				
! Save this quantities for later usage
				dism(dir) = dm
				opam(:,dir,it) = chim
				sfm(:,dir,it) = sm
				
				disp(dir) = dp
				opap(:,dir,it) = chip
				sfp(:,dir,it) = sp
								
				dtm = 0.5d0 * (chi0 + chim) * dm
				dtp = 0.5d0 * (chi0 + chip) * dp

				do frq = 1, nfrq
					if (dtm(frq) > trick) then
               	exu(frq) = dexp(-dtm(frq))
            	else
               	exu(frq) = 1.d0 - dtm(frq) + dtm(frq)**2 / 2.d0
            	endif	
					
					IM(frq,dir,it) = Inte(frq,dir,it)
					bond = Inte(frq,dir,it) * exu(frq)

! Parabolic SC
					if (k /= k1) then

! Calculate SC in the IN direction
						call par_sc(dtp(frq),dtm(frq),psim_rev,psi0_rev,psip_rev)
						short_rev = psim_rev * sp(frq) + psi0_rev * s0(frq) + psip_rev * sm(frq)

! Calculate SC in the OUT direction
						call par_sc(dtm(frq),dtp(frq),psim,psi0,psip)
						short = psim * sm(frq) + psi0 * s0(frq) + psip * sp(frq)

						Inte(frq,dir,it) = bond + short
						true(frq) = short_rev

! Take into account that we did not calculate Lstar before
						sdelta = chil(k) * profile(frq,dir,k,inout) / chi0(frq)
						short = (psi0 + psi0_rev) * sdelta
						Lsta(k,it) = Lsta(k,it) + 0.5d0 * wtmu(dir) * frqwt(frq,dir,k,inout) * &
							profile(frq,dir,k,inout) * short
					else

! Linear SC
						call lin_sc(dtm(frq),psim,psi0)
						short = psim * sm(frq) + psi0 * s0(frq)
						Inte(frq,dir,it) = bond + short
						true(frq) = 0.d0

						sdelta = chil(k) * profile(frq,dir,k,inout) / chi0(frq)
						short = psi0 * sdelta
						Lsta(k,it) = Lsta(k,it) + 0.5d0 * wtmu(dir) * frqwt(frq,dir,k,inout) * &
							profile(frq,dir,k,inout) * short

					endif			

				enddo !frq

! Calculate the contribution to Jbar of each point (correcter with the IN intensity)
				Jba(k,it) = Jba(k,it) + 0.5d0 * sum( wtmu(dir) * frqwt(:,dir,k,inout) * &
					profile(:,dir,k,inout) * ( true + Inte(:,dir,it) ) )

			enddo	!dir	
			
		enddo !it
		
!---------------------
! MAKE POPULATION CORRECTION
!---------------------
! Build up the Lstar_total and Jbar_total matrices
		ipt = nr * (k-1)
		do it = 1, nr
			Lstar_total(it+ipt) = Lsta(k,it)
			Jbar_total(it+ipt) = Jba(k,it)
		enddo !it		

! Do the correction
		call rate_eq(k, relat_error, true_error)
		

! Now that we have changed the population, we have to correct the contribution of the boundary
! intensities to this point before going out
! Only do it if the complete GS is being used		

		if (which == 1) then
		
			cortes = angleset_pp
			do it = 1, nr

				call opacity(it, nip, k, k)

				do dir = 1, cortes

					dm = dism(dir)
					dp = disp(dir)
					chi0 = chil(k) * profile(:,dir,k,inout) + kappa(k) + kappa_dust(k)
					s0 = ( chil(k) * profile(:,dir,k,inout) * Sl(k) + kappa(k) * B(k) + kappa_dust(k) * B_dust(k) ) / chi0

					sm = sfm(:,dir,it)
					sp = sfp(:,dir,it)

					dtm = 0.5d0 * (opam(:,dir,it) + chi0) * dm
					dtp = 0.5d0 * (opap(:,dir,it) + chi0) * dp

					do frq = 1, nfrq
						if (dtm(frq) > trick) then
      	         	exu(frq) = dexp(-dtm(frq))
         	   	else
	               	exu(frq) = 1.d0 - dtm(frq) + dtm(frq)**2 / 2.d0
   	         	endif

						bond = IM(frq,dir,it) * exu(frq)

						if (k /= k1) then
							call par_sc(dtm(frq),dtp(frq),psim,psi0,psip)
							short = psim * sm(frq) + psi0 * s0(frq) + psip * sp(frq)
						else
							call lin_sc(dtm(frq),psim,psi0)
							short = psim * sm(frq) + psi0 * s0(frq)
						endif

						Inte(frq,dir,it) = bond + short

					enddo !frq

				enddo !dir

			enddo !it		
			
		endif
		
	enddo !k

	deallocate(chim)
	deallocate(chi0)
	deallocate(chip)
	deallocate(sm)
	deallocate(s0)
	deallocate(sp)
	deallocate(dtm)
	deallocate(dtp)
	deallocate(exu)
	deallocate(true)

	deallocate(Inte)
	deallocate(IM)
	deallocate(opam)
	deallocate(opap)
	deallocate(dtmold)
	deallocate(smold)
	deallocate(sfm)
	deallocate(sfp)

	deallocate(dism)
	deallocate(disp)
	
	deallocate(Lsta)
	deallocate(Jba)
		
	
	end subroutine rtpp_gs
end module formal_gs_pp
