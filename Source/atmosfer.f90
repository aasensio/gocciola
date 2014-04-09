module atmosfer
use variables
use general
implicit none
contains

! *********************************************************
! *********************************************************
! ATMOSPHERE ROUTINES
! *********************************************************
! *********************************************************

!-----------------------------------------------------------------
! Calculates optical depth in each line
!-----------------------------------------------------------------		
	
	subroutine calc_tau
	integer :: it, i, ipl, up, low
	real(kind=8) :: sig, glu, acon, agnus, chilp, chilm
	
		do it = 1, nact
			tau(it,1) = 0.d0
			up	= itran(1,it)
			low = itran(2,it)
			glu = dlevel(2,low) / dlevel(2,up)
			acon = PHC * (dtran(2,it)) ** 3.d0   !2*h*nu^3/c^2
			agnus = acon * glu
			
			chilm = 0.d0
			do i = 1, n_radios
				ipl = nl * (i-1)
				sig = dtran(1,it) / doppler_width(i,it)
				chilp = perfil(nfrq/2-1) * sig * ( pop(low+ipl) - glu * pop(up+ipl) ) + chic(it,i)
				
! ARREGLAR ESTO				
				if (i /= 1) then
					tau(it,i) = tau(it,i-1) + 0.5d0 * (chilm + chilp) * (r(i) - r(i-1))
				endif
				chilm = chilp
			enddo
			
		enddo
		
	end subroutine calc_tau

! ---------------------------------------------------------
! Fills mu variable with the angle cosines between the ray and the normal
! to each intersection point between the shell and the characteristic
! ---------------------------------------------------------		
	subroutine mugrid
	integer :: i, cortes
		
		mu = 2.d0
		do i = 1, n_radios
			if (i == n_radios) then
				cortes = caract_core + i - 1
			else 
				cortes = caract_core + i
			endif
			
			mu(i,1:cortes) = dsqrt( r(i)**2.d0 - p(1:cortes)**2.d0 ) / r(i) 			
			
		enddo
		
	end subroutine mugrid

! ---------------------------------------------------------
! Reads atmosphere from file
! INPUT:
!     - fich : file where the atmosphere is defined
! ---------------------------------------------------------	
	subroutine read_atmosf(fich)
	character*60, INTENT(IN) :: fich
	integer :: i
	real(kind=8) :: del
		
		open(unit=3,file=fich,status='old',action='read')
		call lb(3,3)
		read(3,*) n_radios

! Allocate memory for shell, characteristic, atmosphere parameters and hydrogen abundance grid
		if (index(geometry_type,'SPHERICAL') /= 0) then
			n_caract = n_radios + caract_core - 1
		else
			n_caract = angleset_pp
		endif
		print *, 'Total characteristics number : ', n_caract
		allocate(r(n_radios))
		allocate(p(n_caract))
! If Zeeman effect is active, then include in the atmosphere model the magnetic field and its 
! orientation
		if (converge_flag == 4) then
			allocate(datph(11,n_radios))
		else
			allocate(datph(8,n_radios))
		endif
		
		allocate(thermal_v(n_radios))
		allocate(nh(n_radios))				
		
		
! If Zeeman effect is included, read the magnetic field	
		if (converge_flag == 4) then
			do i = 1, n_radios
				read(3,*) r(i), datph(1,i), nh(i), datph(3,i), datph(2,i), datph(4,i), datph(6,i),&
						datph(7,i), datph(8,i), datph(9,i), datph(10,i), datph(11,i)
			enddo
		else
		
! If Zeeman effect is included, skip the magnetic field information if available
			do i = 1, n_radios
				read(3,*) r(i), datph(1,i), nh(i), datph(3,i), datph(2,i), datph(4,i), datph(6,i),&
						datph(7,i), datph(8,i)
			enddo		
		endif
		close(3)

! Generate characteristics grid
		if (index(geometry_type,'SPHERICAL') /= 0) then
			p(caract_core+1:n_caract) = r(1:n_radios-1)
			del = r(1) / (caract_core ) 
			do i = 1, caract_core
				p(i) = del * i - del
			enddo
		endif
		
! Calculate visual absorption in each shell. Later on, we will average between two shells to
! obtain visual absorption in each shell
		datph(5,1) = 1.d-21 * nh(1) * r(1)
		do i = 2, n_radios
			datph(5,i) = 1.d-21 * nh(i) * (r(i) - r(i-1))
		enddo	
	
	end subroutine read_atmosf

! ---------------------------------------------------------
! Interpolates the atmosphere resizing all the atmospheric variables
! ---------------------------------------------------------	
	subroutine interpol_atmosphere
	real(kind=8), allocatable :: temp(:), temp2(:,:), eje_anterior(:), eje_nuevo(:)
	integer :: i
	real(kind=8) :: del
		
		n_radios_old = n_radios
		n_caract_old = n_caract

		n_radios = interpolate_atmosphere

! Allocate the two axis for the interpolation
		allocate(eje_anterior(n_radios_old))
		allocate(eje_nuevo(n_radios))
		do i = 1, n_radios_old
			eje_anterior(i) = 1.d0 * i / n_radios_old
		enddo
		do i = 1, n_radios
			eje_nuevo(i) = 1.d0 * i / n_radios
		enddo

		print *, 'Interpolating the atmosphere...'

! Allocate memory for shell, characteristic, atmosphere parameters and hydrogen abundance grid
		n_caract = n_radios + caract_core - 1

		print *, 'New total characteristics number : ', n_caract

		allocate(temp(n_radios))
		call spline(eje_anterior,r,eje_nuevo,temp,'NO  ')
		deallocate(r)
		allocate(r(n_radios))
		r = temp
		
		deallocate(p)
		allocate(p(n_caract))
		
! If Zeeman effect is active, then include in the atmosphere model the magnetic field and its 
! orientation
		if (converge_flag == 4) then
			allocate(temp2(11,n_radios))
			do i = 1, 11
				call spline(eje_anterior,datph(i,:),eje_nuevo,temp2(i,:),'NO  ')
			enddo
			deallocate(datph)
			allocate(datph(11,n_radios))
			datph = temp2
			deallocate(temp2)
		else
			allocate(temp2(8,n_radios))
			do i = 1, 8
				call spline(eje_anterior,datph(i,:),eje_nuevo,temp2(i,:),'NO  ')
			enddo
			deallocate(datph)
			allocate(datph(8,n_radios))
			datph = temp2
			deallocate(temp2)
		endif
		
		
		call spline(eje_anterior,nh,eje_nuevo,temp,'NO  ')
		deallocate(nh)
		allocate(nh(n_radios))
		nh = temp

		deallocate(thermal_v)
		allocate(thermal_v(n_radios))
				
		
! Generate characteristics grid
		if (index(geometry_type,'SPHERICAL') /= 0) then
			p(caract_core+1:n_caract) = r(1:n_radios-1)
			del = r(1) / (caract_core ) 
			do i = 1, caract_core
				p(i) = del * i - del
! Get the number of impact parameters crossing the central star in boundary type 3
				if (p(i) < star_radius) impact_star_radius = i
			enddo
		endif
		
! Calculate visual absorption in each shell. Later on, we will average between two shells to
! obtain visual absorption in each shell
		datph(5,1) = 1.d-21 * nh(1) * r(1)
		do i = 2, n_radios
			datph(5,i) = 1.d-21 * nh(i) * (r(i) - r(i-1))
		enddo	
		

		deallocate(eje_anterior)
		deallocate(eje_nuevo)
		deallocate(temp)

	end subroutine interpol_atmosphere
	
end module atmosfer
