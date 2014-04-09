module background_opac
use variables
use background_opacity_module, only : background_opacity
implicit none
contains
! *********************************************************
! *********************************************************
! GLOBAL VARIABLES
! *********************************************************
! *********************************************************

! ---------------------------------------------------------
! Initialize all variables
! ---------------------------------------------------------

	subroutine background
	real(kind=8) :: chio, etao, chi_e
! 	real*4 :: freq, tel, pel, pg, rho
	real(kind=8) :: freq, lambda_in, Tel, Pel, Pg, PHminus, PHplus, PH2, PH2plus
	integer :: it, i, ihtot, ifudge, j
		
! The bound-bound transitions		
		ihtot = 2

		etao = 0.d0 !????????
		IFUDGE=1
		print *, 'Continuum opacity for bound-bound transitions...'
		do it = 1, nact
			freq = dtran(2,it)
			do i = 1, n_radios				
				Tel = datph(1,i)
				pel = datph(2,i) * 1.3806503d-16 * Tel
				Pg = datph(8,i) 				

				open(unit=5, file='ion.dat', status='old')					  
! 				call readion(1)
				close(5)
! 				rho = UMA * atw * (pg-pel) / tel / BK
! 				call sopas(1, ihtot, freq, pel, tel, pg, chio, chi_e, etao)
				lambda_in = 2.99792458d10 / freq * 1.d8
				call background_opacity(tel, pel, pg, PHminus, PHplus, PH2, PH2plus, lambda_in, chio, chi_e)
				
! Background opacities are divided by rho. But we need it in cm^-1
				chic(it,i) = (chio + chi_e) !* rho
				Sc(it,i) = etao / chio 
			
			enddo

		enddo
		
		print *, 'Continuum opacity for bound-free transitions...'
! Now the bound-free transitions		
		do it = nact+1, nact+nbf
						
			do i = 1, n_radios				
				Tel = datph(1,i)
				pel = datph(2,i) * 1.3806503d-16 * Tel
				Pg = datph(8,i) 				

				open(unit=5, file='ion.dat', status='old')					  
! 				call readion(1)
				close(5)
! 				rho = UMA * atw * (pg-pel) / tel / BK
				
				do j = 1, nfrq_bf
					freq = bf_cross_section(1,it,j)
			
! 					call sopas(1, ihtot, freq, pel, tel, pg, chio, chi_e, etao)       
					lambda_in = 2.99792458d10 / freq * 1.d8
					call background_opacity(tel, pel, pg, PHminus, PHplus, PH2, PH2plus, lambda_in, chio, chi_e)

! Background opacities are divided by rho. But we need it in cm^-1
					chic_bf(it,j,i) = (chio + chi_e) !* rho
					Sc_bf(it,j,i) = etao / chio 
				enddo

			enddo

		enddo
		
	end subroutine background
	
	
end module background_opac
