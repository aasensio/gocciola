module accel
use variables
use statis
use maths
implicit none
contains
! *********************************************************
! *********************************************************
! ACCELERATION ROUTINES
! *********************************************************
! *********************************************************

!-----------------------------------------------------------------
! Executes acceleration if they are activated (Ng)
!-----------------------------------------------------------------
	subroutine aceleracion

		if (iacel == 1 .and. itracc > 0) then
			if (ng(pop,nl*n_radios,nord,scratch)) then
				print *, 'Ng acceleration'
				itracc = -nintracc
			endif
		endif
	
	end subroutine aceleracion
	
! ---------------------------------------------------------
! Performs an Ng acceleration of the solution
! ---------------------------------------------------------

	function ng(y,m,n,yy)
	logical :: ng, ng0
	integer, INTENT(IN) :: n, m
	real(kind=8), INTENT(INOUT) :: y(m), yy(m,*)
	real(kind=8) :: A(7,7), C(7), di, dy, wt
	integer :: k, i, j
	integer, save :: ntry
	data ntry /0/
	
	if (n < 1 .or. n > 5) return
	ntry = ntry + 1
	yy(:,ntry) = y
	ng = .false.
	if (ntry <= n+1) return
		
	A = 0.d0 ; C = 0.d0
	
	do k = 1, m
		wt = 1.d0 / (1.d0 + dabs(y(k)))
		do i = 1, n
			dy = yy(k,ntry-1) - yy(k,ntry)
			di = wt * (dy + yy(k,ntry-i) - yy(k,ntry-i-1) )
			c(i) = c(i) + di*dy
			do j = 1, n
				A(i,j) = A(i,j) + di * (dy + yy(k,ntry-j) - yy(k,ntry-j-1) )
			enddo !j
		enddo !i
	enddo !k
	
	call llslv(A,C,n,7)
	
	ng = .true.
	do i = 1, n
		y(:) = y(:) + c(i) * ( yy(:,ntry-i) - yy(:,ntry) )
	enddo !i

	entry ng0
	ntry = 0
	
	end function ng
end module accel
