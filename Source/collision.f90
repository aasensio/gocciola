module coll_data
use variables
implicit none
contains

! *********************************************************
! *********************************************************
! COLLISION ROUTINES
! *********************************************************
! *********************************************************

! ----------------------------------------------------------
! Calculates collisional rates for CO with H, H2 and He
! Ayres & Wiedemann (1989)  ApJ 338, 1033
! ----------------------------------------------------------

	subroutine co_collision(unit)
	integer, INTENT(IN) :: unit
	integer :: i, from, to
	real*8, PARAMETER :: OMEGA0 = 2103.d0
	real*8, allocatable :: beta(:), temp(:)
	
		allocate(beta(n_radios))
		allocate(temp(n_radios))

		read(unit,*) from, to
		do i = from, to
			beta = PHK * PC * OMEGA0 / datph(1,:)
			read(unit,*) itran(1,i), itran(2,i)
! Collision with H2
!			temp = 19.1d0 - 0.069d0 * 64.d0 * beta**ONETHIRD
!			temp = exp(temp) / ( beta	* (1.d0 - dexp(-beta)))
!			collision(:,i) = 4.2d-13 * temp * nh(:) 
			
			temp = PK * datph(1,:) * dexp(19.1d0 - 64.d0 * datph(1,:)**(-ONETHIRD)) 
			temp = temp / (1.d0 - dexp(-beta)) 
			collision(:,i) = temp * nh(:) / 50.d0
			
! Collision with H
!			temp = 18.1d0 - 0.069d0 * 3.d0 * beta**ONETHIRD
!			temp = exp(temp) / ( beta	* (1.d0 - dexp(-beta)))
!			collision(:,i) = collision(:,i) + 4.2d-13 * temp * nh(:) * xhid
			
			temp = PK * datph(1,:) * dexp(18.1d0 - 3.d0 * datph(1,:)**(-ONETHIRD)) 
			temp = temp / (1.d0 - dexp(-beta)) 
			collision(:,i) = collision(:,i) !+ temp * nh(:) * xhid 
			
! Collision with He
!			temp = 19.1d0 - 0.069d0 * 87.d0 * beta**ONETHIRD
!			temp = exp(temp) / ( beta	* (1.d0 - dexp(-beta)))
!			collision(:,i) = collision(:,i) + 4.2d-13 * temp * nh(:) * xhel

			temp = PK * datph(1,:) * dexp(19.1d0 - 87.d0 * datph(1,:)**(-ONETHIRD)) 
			temp = temp / (1.d0 - dexp(-beta)) 
			collision(:,i) = collision(:,i) !+ temp * nh(:) * xhel
			
! Collision with electrons (using Thompson's schematic cross section)
			temp = (1.d0 + beta) + 19.d0 * dexp(-3.22d0 * beta) * (1.d0 + 4.22*beta)
			temp = 1.4d-9 / dsqrt(beta) * temp
			collision(:,i) = collision(:,i) !+ temp * datph(2,i) 
			collision(:,i) = 0.d0
		enddo

		deallocate(beta)
		deallocate(temp)
	end subroutine co_collision


end module coll_data
