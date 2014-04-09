module overlap
use variables
use general
implicit none
contains

! *********************************************************
! *********************************************************
! OVERLAPPING ROUTINES
! *********************************************************
! *********************************************************

	
! ---------------------------------------------------------
! Fills the structure of overlapping lines
! overl(tran_i,n_overl,shell) : this array indicates several things:
!		- The first column indicates the number of lines which overlap a given one (n)
!		- For each row given by tran_i, there are n different numbers when changing the second
!			index (n_overl) which indicate the transition with which there is an overlap
! ---------------------------------------------------------

	subroutine init_overlap
	integer :: i, n_overl, j, l
	real(kind=8) :: f1, f2, df1, df2, distancia, criterium
	
! We go through every transition obtaining the overlappings with every other transition
		overl = 0
		do l = 1, n_radios
			do i = 1, nr
				f1 = dtran(2,i)
				df1 = doppler_width(l,i)
				n_overl = 0
				j = 1
				do while(j <= nr)
! If i and j are different transitions				
					if (j /= i) then
						f2 = dtran(2,j)						
						df2 = doppler_width(l,j)
						
						distancia = dabs(f1-f2)
!						print *, 'ditancia : ',distancia
						criterium = (df1+df2) * FWHM
!						print *, 'criterio : ',criterium
						
						if (distancia < criterium) then
							n_overl = n_overl + 1
							overl(n_overl+1,i,l) = j
						endif
						
					endif
					j = j + 1
				enddo
! Now store the number of overlapping transitions for transition i				
				overl(1,i,l) = n_overl
			enddo	
						
		enddo
				
	end subroutine init_overlap
	
! ---------------------------------------------------------
! Returns the opacity of the overlapping transitions
! Given the transition, the radius, the caract and the direction, it returns the
! opacity to add to the normal opacity due to the overlapping transitions and the source
! function. It returns a frequency axis.
! ---------------------------------------------------------
	subroutine overlap_opacity(transition,r,caract,inout,chi_ov,emm_ov)
	integer, INTENT(IN) :: transition, r, caract, inout
	real(kind=8), INTENT(INOUT) :: chi_ov(nfrq), emm_ov(nfrq)
	integer :: i, j, ifrq, n_overl
	real(kind=8) :: chi_in, Sl_in, emmis_in, x, x_new, angmu, dir, prof
	
	
		chi_ov = 0.d0
		emm_ov = 0.d0
		
! Get the number of overlapping transitions
		n_overl = overl(1,transition,r)
		
! Get the mu angle at this point	and the direction
		angmu = mu(r,caract)
		dir = (-1.d0)**inout
		
! Go through every frequency
		do ifrq = 1, nfrq

			x = 1.d0*(ifrq - nfrq / 2 - 1 ) / n_freq_ud_doppler - dir * angmu * datph(7,r)	
! Frequencies of the actual transition
			x = dtran(2,transition) + x * doppler_width(r,transition)

! Go through every transition which overlaps the given one
			do i = 1, n_overl
		
! Get which transitions are
				j = overl(i+1,transition,r)
			
! Get the opacity and source function (and emmisivity) of the overlapping transition			
				call opacity_one(j,r,chi_in,Sl_in)
				emmis_in = Sl_in * chi_in
				x_new = (x-dtran(2,j)) / doppler_width(r,j)							
				prof = 1.d0 / SQRTPI * voigt(0.d0,x_new,0)
				chi_ov(ifrq) = chi_ov(ifrq) + chi_in * prof
				emm_ov(ifrq) = emm_ov(ifrq) + emmis_in * prof
				
			enddo !i
			
		enddo !ifrq
		
	end subroutine overlap_opacity

end module overlap
