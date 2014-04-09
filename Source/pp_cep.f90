module formal_cep
use variables
use general
use short_ch
use maths
use atmosfer
use statis
implicit none
contains

! *********************************************************
! *********************************************************
! CEP METHOD SOLUTION ROUTINES
! *********************************************************
! *********************************************************
   
!---------------------------------------------------------
! This function calculates the beta function
!---------------------------------------------------------
	function beta2(tau_in)
   real(kind=8) :: beta2, tau_in, salida, tau, coef
   real(kind=8) :: d, b, q, dbdx
   integer :: k
   
! KROLIK & McKEE APPROXIMATION   
   	coef = 1.d0/(6.d0*sqrt(3.d0))
   	tau = tau_in / dsqrt(PI)
   	if (tau < 1.d-4) then 
			salida = 1.d0
		else
			if (tau >= 3.41d0) then
   			salida = 1.d0 / (dsqrt(PI)*tau) * ( dsqrt(log(tau)) + 0.25d0 / dsqrt(log(tau)) + 0.14d0)
   		else
				salida = 1.d0 - 0.415d0*tau + 0.355d0*tau*log(tau)
   			dbdx = 0.355 -0.415 + 0.355*log(tau)
   			k = 1
   			d = tau * coef
   			b = 2.d0 * d
   			q = 1.d0
          	do while (q > 1.d-3)          
            	salida = salida - d * tau
            	dbdx = dbdx - b
             	k = k + 1             	
             	d = -d * tau * sqrt((k+1.d0)/(k+2.d0))*(k-1.d0)/(k*(k+2.d0))
             	b = (k+1.d0) * d
             	q = abs(b/dbdx)
            enddo	    			
	   	endif
		endif
   	
! Calculation of the beta function by integrating the E3 exponential integral function
!   	if (tau_in == 0.d0) then
!			beta2 = 1.d0
!			return
!   	endif
		
!		salida = 0.d0
!   	if (tau_in /= 0.d0) then
!   		do i = 1, n_quadr_beta
!   			salida = salida + w_e3(i)*(0.5d0-expint(3,tau_in*exp(-x_e3(i)**2)/sqrt(PI))) / tau_in
!   		enddo
!   	endif
   	   	
   	beta2 = salida
   	   	
   end function beta2

!-------------------------------------------------------------------
! Solves a system of nonlinear equations using the Newton method
!-------------------------------------------------------------------
	subroutine mnewt(ntrial,x,n,tolx,tolf,derivatives)
	integer :: n,ntrial
	real(kind=8) :: tolf,tolx,x(n)
	logical :: derivatives
	integer :: i,k
	real(kind=8) :: errf,errx,fjac(n,n),fvec(n),p(n), tol
	integer, allocatable :: indx(:)
			
		errf = 1.d10
		do k=1,ntrial
			call usrfun(x,n,fvec,fjac, derivatives)   !User subroutine supplies function values at x in fvec
                
			do i=1,n  !Check function convergence.
				errf=errf+dabs(fvec(i))
			enddo 
			if (errf <= tolf) then
				return
			endif
			p = -fvec

			allocate(indx(n))
			tol = 1.d-10
			call ludcmp(fjac,indx,tol)
			call lubksb(fjac,indx,p)
			deallocate(indx)
		                			                
			errx=0.d0  ! Check root convergence.
			x = x + p
                
			do i=1,n   !Update solution.
				errx=errx+dabs(p(i))                    
			enddo 
			if(errx <= tolx) then				
				print *, 'Iteration ',k, '   rel_error: ',errx
				write(12,*) k, errx
				return
			endif
			
			print *, 'Iteration ',k, '   rel_error: ',errx
			write(12,*) k, errx
		enddo 		
        
	end subroutine mnewt
	
!-----------------------------------------------------------------
! Calculate Jbar_total and Lstar_total for this model and atmosphere
!-----------------------------------------------------------------
	subroutine formal_solver_jacobi_cep(x)
	real(kind=8) :: x(:)
	integer :: up, low, it, ip, ipl, ipt, ip1
	real(kind=8) :: sig, glu, acon, chim, chip, tau0
		
		do it = 1, nr

! Calculate opacities and source functions
!			call opacity(it,npmip1,1,n_radios)
			up = itran(1,it)
			low = itran(2,it)
			sig = dtran(1,it) / dtran(4,it) !sig=(h*nu*Blu)/(4*pi*Dnudoppler)
			glu = dlevel(2,low) / dlevel(2,up)
			acon = 1.4743d-2 * (1.d-15 * dtran(2,it))**3  !acon = (2*h*nu**3)/c**2

! Line opacity, line source function and Planck function. 			

			tau0 = 0.d0
			chip = 0.d0
			chim = 0.d0
			ip1 = 1
			do ip = 1, n_radios
				ipl = nl * (ip-1)
				
				chil(ip) = sig * (x(low+ipl) - glu * x(up+ipl))
				Sl(ip) = acon * glu * x(up+ipl) / (x(low+ipl) - glu * x(up+ipl))
				
				dSldn(ip,up+ipl) = acon*glu / (x(low+ipl) - glu * x(up+ipl))**2 * x(low+ipl)
				dSldn(ip,low+ipl) = -acon*glu / (x(low+ipl) - glu * x(up+ipl))**2 * x(up+ipl)
				
				B(ip) = acon * glu * popl(up+ipl) / (popl(low+ipl) - glu * popl(up+ipl))
								
			enddo
						
			if (verbose_flag == 1) then
				if (mod(it,10) == 0) then
					write(*,*) '       Transition : ', it
				endif
         endif
			
! Calculate outer boundary condition for this transition
        	call calcula_fondo(it)

! Calculate inner boundary condition for this transition
        	call calcula_centro(it)

! Do the formal solution of the RTE
			call rtpp_cep(n_radios, it)	

! Build up the Lstar_total and Jbar_total matrices. These are necessary for the construction
! of the rate matrix
			do ip = 1, n_radios
				ipt = nr*(ip-1)
				Lstar_total(it+ipt) = lstar(ip)
				Jbar_total(it+ipt) = Jbar(ip)
				dJbar_totaldn(it+ipt,:) = dJbardn(ip,:)
			enddo			
			
		enddo
		
	end subroutine formal_solver_jacobi_cep
	
!---------------------------------------------------------
! This subroutine solves the radiative transfer equation using the 
! coupled escape probability scheme developed by Moshe Elitzur
!---------------------------------------------------------
   subroutine rtpp_cep(np,it)
	integer :: np, it
   integer :: i, j
   real(kind=8) :: tau_local, escape, delta_tau(4), alpha_tau(4)

   Jbar = 0.d0
	Lstar = 0.d0
!	flux = 0.d0
	dJbardn = 0.d0

	do i = 1, np
		tau_local = dabs(tau(it,i) - tau(it,i-1))
	
		escape = beta2(tau_local)
		
		Jbar(i) = Sl(i) * (1.d0 - escape)
		dJbardn(i,:) = dSldn(i,:) * (1.d0 - escape)
		
		Lstar(i) = (1.d0 - escape)
				
! Boundary condition
! There is a problem because the CEP method is based on zones, rather than in grid points
		escape = beta2(tau(it,i))
		Jbar(i) = Jbar(i) + 0.5d0*(fuera(1)+dentro(1)) * tau(it,i) * escape / tau_local
				
		escape = beta2(tau(it,i-1))
		Jbar(i) = Jbar(i) - 0.5d0*(fuera(1)+dentro(1)) * tau(it,i-1) * escape / tau_local
								
		do j = 1, np
			if (j /= i) then
				alpha_tau = 0.d0
				
				delta_tau(1) = dabs(tau(it,i) - tau(it,j))
				escape = beta2(delta_tau(1))
				alpha_tau(1) = delta_tau(1) * escape 
												
				delta_tau(2) = dabs(tau(it,i-1) - tau(it,j))
				escape = beta2(delta_tau(2))				
				alpha_tau(2) = delta_tau(2) * escape 
												
				delta_tau(3) = dabs(tau(it,i) - tau(it,j-1))				
				escape = beta2(delta_tau(3))
				alpha_tau(3) = delta_tau(3) * escape 
								
				delta_tau(4) = dabs(tau(it,i-1) - tau(it,j-1))
				escape = beta2(delta_tau(4))
				alpha_tau(4) = delta_tau(4) * escape 
								
				Jbar(i) = Jbar(i) + 0.5d0 * Sl(j) / tau_local * &
					(alpha_tau(1) - alpha_tau(2) - alpha_tau(3) + alpha_tau(4))
					
				dJbardn(i,:) = dJbardn(i,:) + 0.5d0 * dSldn(j,:) / tau_local * &
					(alpha_tau(1) - alpha_tau(2) - alpha_tau(3) + alpha_tau(4))
				
			endif
		enddo

!		flux(np) = flux(np) + (1.d0-Jbar(i)/Sl(i)) * Sl(i) * tau_local
	enddo
	
   end subroutine rtpp_cep	
	
!-------------------------------------------------------------------
! Returns the value of the equations of the nonlinear set to be solved
!-------------------------------------------------------------------
	subroutine funcv_analytic(x,n,fvec,fjac)
	integer :: n
	real(kind=8) :: x(n), fvec(n), fjac(n,n)
	integer :: ip
	integer :: low, up, ipl, ipt, i, j, il, it, il2
	real(kind=8) :: A(nl, nl), A2(nl,nl), B(nl), populat(nl), producto(nl), cul, acon, blu, glu, djbar, poptot
	real(kind=8), allocatable :: jac(:,:,:)

		allocate(jac(nl,nl,n))
		
		call calc_tau		
		call formal_solver_jacobi_cep(x)
				
		do ip = 1, n_radios

! Offsets for storing
			ipl = nl*(ip-1)
			ipt = nr*(ip-1)
		
			A = 0.d0
			jac = 0.d0
		
! Upward collisional rates A(low,up)
			call collis(ip,A)		

! Calculate the downward rates and interchange the positions
			do i = 1, nl
				do j = 1, i-1
					cul = A(i,j)
! Rate into de j-th level				
					A(j,i) = cul
! Rate into de i-th level				
					A(i,j) = popl(i+ipl) / popl(j+ipl) * cul
				enddo
			enddo
					
! Active radiative rates (preconditioning)
			do it = 1, nact
		
				up = itran(1,it)
				low = itran(2,it)
				acon = 1.4743d-2 * (1.d-15 * dtran(2,it))**3  !acon = (2*h*nu**3)/c**2
				blu = PI4H * dtran(1,it) / dtran(2,it)
				glu = dlevel(2,low) / dlevel(2,up)
			
				djbar = Jbar_total(it+ipt)
			
				A(low,up) = A(low,up) + glu * blu * (acon + djbar )				
				A(up,low) = A(up,low) + blu * djbar
				
				jac(low,up,:) = jac(low,up,:) + glu * blu * dJbar_totaldn(it+ipt,:)
				jac(up,low,:) = jac(up,low,:) + blu * dJbar_totaldn(it+ipt,:)
			
			enddo		
		
! Background radiative rates
			do it = nact + 1, nr
		
				up = itran(1,it)
				low = itran(2,it)
				acon = 1.4743d-2 * (1.d-15 * dtran(2,it))**3  !acon = (2*h*nu**3)/c**2
				blu = PI4H * dtran(1,it) / dtran(2,it)
				glu = dlevel(2,low) / dlevel(2,up)
			
				A(low,up) = A(low,up) + glu * blu * (acon + Jbar_total(it+ipt))
				A(up,low) = A(up,low) + blu * Jbar_total(it+ipt)
				
				jac(low,up,:) = jac(low,up,:) + glu * blu * dJbar_totaldn(it+ipt,:)
				jac(up,low,:) = jac(up,low,:) + blu * dJbar_totaldn(it+ipt,:)
			
			enddo		
		
! The system is homogeneous with null determinant. We have to substitute one of the equations by a closure
! relation. The sum over one column gives the total output rate from level i
			do i = 1, nl

				do j = 1, i-1
					A(i,i) = A(i,i) - A(j,i)
					jac(i,i,:) = jac(i,i,:) - jac(j,i,:)
				enddo
				do j = i+1, nl
					A(i,i) = A(i,i) - A(j,i)
					jac(i,i,:) = jac(i,i,:) - jac(j,i,:)
				enddo
				
			enddo		
			
			A2 = A
					
! Conservation equation
			B = 0.d0
			A(nl, 1:nl) = 1.d0
			jac(nl, 1:nl, :) = 0.d0
		
			poptot = datph(3,ip) * nh(ip)
			B(nl) = poptot
			
			do il = 1, nl
				populat(il) = x(il+ipl)
			enddo
			
! Return the rate equations			
			producto = matmul(A,populat)-B
						
			do il = 1, nl				
				fvec(il+ipl) = producto(il)								
			enddo

! Analytical Jacobian. Since the equations are F_j = sum(A_jk * n_k, k) - b_j, the Jacobian J_ji is 
! J_ji = dF_j/dn_i = sum(dA_jk/dn_i * n_k,k) + sum(A_jk * dn_k/dn_i,k) = sum(dA_jk/dn_i * n_k,k) + A_ji
			do il = 1, nl
				do il2 = 1, nl
					fjac(il+ipl,il2+ipl) = A(il,il2)
				enddo
				do il2 = 1, nl*n_radios
					fjac(il+ipl,il2) = fjac(il+ipl,il2) + sum(jac(il,:,il2)*populat)
				enddo
			enddo
						
		enddo
		
		deallocate(jac)
		
	end subroutine funcv_analytic	

!-------------------------------------------------------------------
! Returns the value of the equations of the nonlinear set to be solved
!-------------------------------------------------------------------
	subroutine funcv(x,n,fvec)
	integer :: n
	real(kind=8) :: x(n), fvec(n)
	integer :: ip
	integer :: low, up, ipl, ipt, i, j, il, it
	real(kind=8) :: A(nl, nl), B(nl), populat(nl), producto(nl), cul, acon, blu, glu, djbar, poptot

		call formal_solver_jacobi_cep(x)
		
		do ip = 1, n_radios

! Offsets for storing
			ipl = nl*(ip-1)
			ipt = nr*(ip-1)
		
			A = 0.d0
		
! Upward collisional rates A(low,up)
			call collis(ip,A)		

! Calculate the downward rates and interchange the positions
			do i = 1, nl
				do j = 1, i-1
					cul = A(i,j)
! Rate into de j-th level				
					A(j,i) = cul
! Rate into de i-th level				
					A(i,j) = popl(i+ipl) / popl(j+ipl) * cul
				enddo
			enddo
		
! Active radiative rates (preconditioning)
			do it = 1, nact
		
				up = itran(1,it)
				low = itran(2,it)
				acon = 1.4743d-2 * (1.d-15 * dtran(2,it))**3  !acon = (2*h*nu**3)/c**2
				blu = PI4H * dtran(1,it) / dtran(2,it)
				glu = dlevel(2,low) / dlevel(2,up)
			
				djbar = Jbar_total(it+ipt)
			
				A(low,up) = A(low,up) + glu * blu * (acon + djbar )
				A(up,low) = A(up,low) + blu * djbar
			
			enddo		
		
! Background radiative rates
			do it = nact + 1, nr
		
				up = itran(1,it)
				low = itran(2,it)
				acon = 1.4743d-2 * (1.d-15 * dtran(2,it))**3  !acon = (2*h*nu**3)/c**2
				blu = PI4H * dtran(1,it) / dtran(2,it)
				glu = dlevel(2,low) / dlevel(2,up)
			
				A(low,up) = A(low,up) + glu * blu * (acon + Jbar_total(it+ipt))
				A(up,low) = A(up,low) + blu * Jbar_total(it+ipt)
			
			enddo		
		
! The system is homogeneous with null determinant. We have to substitute one of the equations by a closure
! relation. The sum over one column gives the total output rate from level i
			do i = 1, nl

				do j = 1, i-1
					A(i,i) = A(i,i) - A(j,i)
				enddo
				do j = i+1, nl
					A(i,i) = A(i,i) - A(j,i)
				enddo
				
			enddo		
		
! Conservation equation
			B = 0.d0
			A(nl, 1:nl) = 1.d0
		
			poptot = datph(3,ip) * nh(ip)
			B(nl) = poptot
			
			do il = 1, nl
				populat(il) = pop(il+ipl)
			enddo
			
! Return the rate equations			
			producto = matmul(A,populat)-B
						
			do il = 1, nl				
				fvec(il+ipl) = producto(il)
			enddo		
		enddo			
		
	end subroutine funcv
	
!-------------------------------------------------------------------
! Calculates the Jacobian using forward-differences
!-------------------------------------------------------------------            
	subroutine fdjac(n,x,fvec,df)
	integer :: n
	real(kind=8) :: df(n,n),fvec(n),x(n),EPS
	PARAMETER (EPS=1.d-4)
	integer :: i,j
	real(kind=8) :: h,temp,f(n)
		
		do j=1,n			
			temp=x(j)
			h=EPS*dabs(temp)
			if(h.eq.0.d0)h=EPS
			x(j)=temp+h
			h=x(j)-temp
			call funcv(x,n,f)
			x(j)=temp
			do i=1,n 
				df(i,j)=(f(i)-fvec(i))/h
			enddo 
		enddo 
        
	end subroutine fdjac	
	
!-------------------------------------------------------------------
! Returns the equations and the Jacobian of the nonlinear set to be solved at point x
!-------------------------------------------------------------------
	subroutine usrfun(x,n,fvec,fjac, derivatives)
	integer :: n
	logical :: derivatives
	real(kind=8) :: x(n), fvec(n), fjac(n,n)
	
		fjac = 0.d0
		
! Analytical derivatives
		if (derivatives) then
			call funcv_analytic(x,n,fvec,fjac)
		else
! Numerical derivatives
			call funcv(x,n,fvec)
			call fdjac(n,x,fvec,fjac)
		endif
				
	end subroutine usrfun	

end module formal_cep
