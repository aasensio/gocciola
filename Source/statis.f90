module statis
use variables
use maths
implicit none
contains
! *********************************************************
! *********************************************************
! RATE EQUATIONS ROUTINES
! *********************************************************
! *********************************************************

!-----------------------------------------------------------------
! Puts collisional rates in the matrix C for the atmosphere point ip
! C(i,j) is the downward collision rate (i>j)
! In the following subroutine C(low,up) goes in A(up,low)
!-----------------------------------------------------------------
	subroutine collis(ip,C)
	integer, INTENT(IN) :: ip
	real(kind=8), INTENT(INOUT) :: C(:, :)
	integer :: it, up, low

	C = 0.d0	
	if (include_collis == 1) then	
		do it = 1, nt
			up = itran(1,it)
			low = itran(2,it)			
			C(up,low) = collision(ip,it)
		enddo
	endif
		
	end subroutine collis


!-----------------------------------------------------------------
! Solves statistical equilibrium equations for a point in the atmosphere (ip)
! It is based on the preconditioning of Rybicki & Hummer (1991) following the scheme presented in
! Socas-Navarro & Trujillo Bueno (ApJ 490,383 1997) for obtaning the approximate corrections
! to the populations
!-----------------------------------------------------------------
	subroutine rate_eq(ip, dxmax, true_e)
	integer, INTENT(IN) :: ip
	real(kind=8), INTENT(INOUT) :: dxmax, true_e
	integer :: low, up, ipl, ipt, i, j, il, it
	real(kind=8) :: cul, acon, blu, glu, djbar, poptot
	integer, allocatable :: eliminate(:), indx(:)
	real(kind=8), allocatable :: A(:,:), L_old(:,:), A_LU(:,:), n_old(:), B(:), delta(:)
	real(kind=8) :: dxmaxold, ratio_lte

		allocate(A(nl,nl))
		allocate(A_LU(nl,nl))
		allocate(L_old(nl,nl))
		allocate(n_old(nl))
		allocate(B(nl))
		allocate(delta(nl))
		allocate(eliminate(1))
		allocate(indx(nl))

		
! If Lambda-iteration is chosen, then use Lambda* = 0
		if (tipo_iteracion == 1 .or. tipo_iteracion == 2) Lstar_total = 0.d0
		
! Offsets for the storage
		ipl = nl*(ip-1)
		ipt = nr*(ip-1)

		A = 0.d0
		
! Upward collisional rates A(low,up). The non-preconditioned matrix L_old is the same as the preconditioned
! matrix, because the collisional rates need not preconditioning
 		call collis(ip,A)
 		
!  		Jbar_total = 0.d0

! Calculate downward rates and interchange positions
		do i = 1, nl
			do j = 1, i-1
				cul = A(i,j)
! Rate into de j-th level				
				A(j,i) = cul
! Rate into de i-th level				
! Detailed balance
				A(i,j) = popl(i+ipl) / popl(j+ipl) * cul
			enddo
		enddo
		L_old = A

		
! Active radiative rates (preconditioning)
		do it = 1, nact

		
			up = itran(1,it)
			low = itran(2,it)
			acon = PHC * (dtran(2,it))**3  !acon = (2*h*nu**3)/c**2
			blu = PI4H * dtran(1,it) / dtran(2,it)
			glu = dlevel(2,low) / dlevel(2,up)

! dJbar = Jbar - Lambda * Slu^old			
			djbar = Jbar_total(it+ipt) - Lstar_total(it+ipt) * acon * pop(up+ipl) * glu /&
					 ( pop(low+ipl) - glu * pop(up+ipl) )
			
! -n_u^new * B_ul * dJbar - n_u^new * Aul * (1-Lstar^*)
			A(low,up) = A(low,up) + glu * blu * (acon * (1.d0 - Lstar_total(it+ipt)) + djbar )

! n_l^new * B_lu * dJbar
			A(up,low) = A(up,low) + blu * djbar

			L_old(low,up) = L_old(low,up) + glu * blu * (acon + Jbar_total(it+ipt))
			L_old(up,low) = L_old(up,low) + blu * Jbar_total(it+ipt)
			
		enddo

! Active radiative bound-free rates (preconditioning)
		if (include_bound_free == 1) then
			do it = nact+1, nact+nbf

				up = itran(1,it)
				low = itran(2,it)

				ratio_lte = popl(low+ipl) / popl(up+ipl)

! Everything multiplied by -n_u^new
				A(low,up) = A(low,up) + ratio_lte * &
					(Jbar_total_bf_a(it+ipt) - Lstar_total_bf_b(it+ipt) + I_total_bf(it+ipt) - Lstar_total_bf_c(it+ipt))

! Everything multiplied by n_l^new			
				A(up,low) = A(up,low) + (Jbar_total_bf_a(it+ipt) - Lstar_total_bf_a(it+ipt))

				L_old(low,up) = L_old(low,up) + ratio_lte * (Jbar_total_bf_a(it+ipt) + I_total_bf(it+ipt))
				L_old(up,low) = L_old(up,low) + Jbar_total_bf_a(it+ipt)

			enddo
		endif

! Background radiative rates
		do it = nact + 1, nr
		
			up = itran(1,it)
			low = itran(2,it)
			acon = PHC * (dtran(2,it))**3  !acon = (2*h*nu**3)/c**2
			blu = PI4H * dtran(1,it) / dtran(2,it)
			glu = dlevel(2,low) / dlevel(2,up)
			
			A(low,up) = A(low,up) + glu * blu * (acon + Jbar_total(it+ipt))
			A(up,low) = A(up,low) + blu * Jbar_total(it+ipt)

			L_old(low,up) = L_old(low,up) + glu * blu * (acon + Jbar_total(it+ipt))
			L_old(up,low) = L_old(up,low) + blu * Jbar_total(it+ipt)
			
		enddo		

! This linear system is homogeneous with null determinant. We have to substitute one equation
! with a closure relation. Sum over one column returns the total rate out of level i
		do i = 1, nl

			do j = 1, i-1
				A(i,i) = A(i,i) - A(j,i)
				L_old(i,i) = L_old(i,i) - L_old(j,i)
			enddo
			do j = i+1, nl
				A(i,i) = A(i,i) - A(j,i)
				L_old(i,i) = L_old(i,i) - L_old(j,i)
			enddo
			
		enddo		
! Particle conservation equation. Eliminate the equation with the greatest population in the
! last iterative step
		
! 		print *, A
! 		stop
		
		eliminate = maxloc(pop(ipl+1:ipl+nl))		
		B = 0.d0
		A(eliminate(1), 1:nl) = 1.d0
		L_old(eliminate(1), 1:nl) = 1.d0

		
		poptot = datph(3,ip) * nh(ip)  ! abund(referred to H) * nh
		B(eliminate(1)) = poptot

		
! A_LU is the A matrix, not the LU decomposition
		A_LU = A
		
! 		print *, A
! 		stop
		
		do il = 1, nl
			n_old(il) = pop(il+ipl)
		enddo

		
		delta = B - matmul(L_old,n_old)
		
! If SLAP is used, use the previous populations as initialization
		call linsolve(A, delta, indx, 0)

! Maximum relative change. dxmax has to be initialized before entering		
		if (error_flag == 0) then
			do il = 1, nl
! Jacobi and GS (omega=1)
				B(il) = n_old(il) + omega_sor * delta(il)
				dxmaxold = dxmax
				if (B(il) /= 0.d0) dxmax = max(dxmax, dabs( B(il) - n_old(il) ) / dabs(B(il)) )
				if (dxmax /= dxmaxold) level_highest_change = il
! 				write(14,*) ip, il, dabs( B(il) - n_old(il) ) / dabs(B(il))
				pop(il+ipl) = B(il)
			enddo

! True error. dxmax has to be initialized before entering		
		else
			do il = 1, nl
! Jacobi and GS (omega=1)
				B(il) = n_old(il) + omega_sor * delta(il)
				if (B(il) /= 0.d0) true_e = max(true_e, dabs( B(il) - pop_final(il+ipl) ) / &
					dabs(pop_final(il+ipl)) )
! 				write(14,*) ip, il, pop_final(il+ipl), dabs( B(il) - pop_final(il+ipl) ) / &
! 					dabs(pop_final(il+ipl))
				if (B(il) /= 0.d0) dxmax = max(dxmax, dabs( B(il) - n_old(il) ) / dabs(B(il)) )
				pop(il+ipl) = B(il) 
			enddo
		endif

!		write(12,*) pop

		deallocate(A)
		deallocate(A_LU)
		deallocate(B)
		deallocate(L_old)
		deallocate(n_old)
		deallocate(delta)
		deallocate(eliminate)
		deallocate(indx)

	end subroutine rate_eq
	
!-----------------------------------------------------------------
! Makes population correction and calculates the maximum relative change in this step
!-----------------------------------------------------------------
	subroutine corrige_poblaciones		
	real(kind=8), allocatable :: pop_ant(:), delta(:)
	integer :: ip
	
		allocate(pop_ant(nl*n_radios))
		allocate(delta(nl*n_radios))

		relat_error_p = relat_error
		relat_error = 0.d0
		true_error = 0.d0
		pop_ant = pop		
		
		do ip = 1, n_radios
			call rate_eq(ip,relat_error,true_error)			
		enddo
		
! If any population is negative, use only half of the change		
!		delta = pop - pop_ant
!		minim = minval(pop)
!		if (minim < 0.d0) then
!			print *, 'Negative populations. Reducing the change'
!			pop = pop_ant + delta / 2.d0
!		endif
		
!		if (relat_error > relat_error_p) then
!			print *, 'Relat_error : ',relat_error
!
!			pop = pop_ant + delta / 2.d0
!			
!			dxmax = 0.d0
!			do ip = 1, n_radios
!				ipl = nl*(ip-1)
!				do il = 1, nl				
!					if (pop(il+ipl) /= 0.d0) then
!						dxmax = max(dxmax, dabs( pop(il+ipl) - pop_ant(il+ipl) ) / dabs(pop(il+ipl)) )
!					endif
!				enddo
!			enddo
!			relat_error = dxmax
!			
!			print *, 'Reducing change...', relat_error
!		endif

		deallocate(pop_ant)
		deallocate(delta)
		
	end subroutine corrige_poblaciones	
end module statis
