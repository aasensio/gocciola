module flux_module
use variables
implicit none
contains

! *********************************************************
! *********************************************************
! THESE ROUTINES CALCULATE THE EMERGING FLUX IN A CORRECT WAY
! RESAMPLING THE EXTERNAL MU AXIS
! *********************************************************
! *********************************************************

! ---------------------------------------------------------
! Performs the formal solution of the RT equation for a plane-parallel atmosphere
! for a given source function and opacity
! Opacity and source function have the index as
! opacity(depth,frequency)
! ---------------------------------------------------------	
	function formal_sol_vect(height,opacit, source, boundary)
	real*8 :: height(:), opacit(:,:), source(:,:), boundary(:), mu
	real*8 :: formal_sol_vect(size(boundary))
	integer :: k, km, kp, idir, k0, kf, n_depths
	real*8, allocatable :: chim(:), chi0(:), chip(:), sm(:), s0(:), sp(:), Inten(:), dtp(:), dtm(:), exu(:)
	real*8, allocatable :: psim(:), psi0(:), psip(:)
	real*8 :: dm, dp

		n_depths = size(boundary)
		allocate(chim(n_depths))
		allocate(chi0(n_depths))
		allocate(chip(n_depths))
		allocate(sm(n_depths))
		allocate(s0(n_depths))
		allocate(sp(n_depths))
		allocate(Inten(n_depths))
		allocate(dtp(n_depths))
		allocate(dtm(n_depths))
		allocate(exu(n_depths))
		allocate(psim(n_depths))
		allocate(psi0(n_depths))
		allocate(psip(n_depths))		

! Boundary condition
		Inten = boundary
		
		idir = 1
		k0 = 2
		kf = n_depths
				
		do k = k0, kf, idir
			mu = mus(k)
! Parabolic short-characteristics
			if (k /= kf) then
				km = k - idir
				kp = k + idir
				chim = opacit(km,:)
				chi0 = opacit(k,:)
				chip = opacit(kp,:)
				sm = source(km,:)
				s0 = source(k,:)
				sp = source(kp,:)
				dm = dabs(1.d0/mu*((height(k)) - (height(km))))
				dp = dabs(1.d0/mu*((height(kp)) - (height(k))))
			else
! Linear short-characteristics			
				km = k - idir
				chim = opacit(km,:)
				chi0 = opacit(k,:)
				chip = 0.d0
				sm = source(km,:)
				s0 = source(k,:)
				sp = 0.d0
				dm = dabs(1.d0/mu*((height(k)) - (height(km))))
				dp = 0.d0
			endif
			
			where(chim == 0.d0)	
				chim = 1.d-10
			endwhere
			where(chi0 == 0.d0)     
                                chi0 = 1.d-10
                        endwhere
			where(chip == 0.d0)     
                                chip = 1.d-10
                        endwhere
			dtm = 0.5d0 * (chim + chi0) * dm
			dtp = 0.5d0 * (chi0 + chip) * dp
			
			where (dtm >= 1.d-4) 
				exu = dexp(-dtm)
			elsewhere
				exu = 1.d0 - dtm + 0.5d0 * dtm**2.d0
			endwhere
			
			if (k /= kf) then
				call par_sc_vect(dtm,dtp,psim,psi0,psip)
				Inten = Inten * exu + psim*sm + psi0*s0 + psip*sp
			else
				call lin_sc_vect(dtm,psim,psi0)
				Inten = Inten * exu + psim*sm + psi0*s0
			endif
									
		enddo

		formal_sol_vect = Inten

		deallocate(chim)
		deallocate(chi0)
		deallocate(chip)
		deallocate(sm)
		deallocate(s0)
		deallocate(sp)
		deallocate(Inten)
		deallocate(dtp)
		deallocate(dtm)
		deallocate(exu)
		deallocate(psim)
		deallocate(psi0)
		deallocate(psip)
	
	end function formal_sol_vect

end module flux_module