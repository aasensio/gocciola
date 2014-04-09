module atomic
use variables
use general
use coll_data
implicit none
contains

! *********************************************************
! *********************************************************
! ATOMIC OR MOLECULAR MODEL ROUTINES
! *********************************************************
! *********************************************************

! ----------------------------------------------------------
! Calculates collisional coefficients for each transition
! INPUT : 
!		- tipo : type of collision (see README for more details)
! OUTPUT : 
!		- fills collision(frq,shell) with the correct value
! It reads the initial and final transition which accounts for these type of transition
! This way, one can use different schemes of collisions for different transitions (for example,
! use a formula for ro-vibrational CO transitions and a table for interpolation for pure rotational
! transitions)
! ----------------------------------------------------------
	subroutine colisiones(tipo)
	character(len=20), INTENT(IN) :: tipo
	character(len=4) :: option
	integer :: i, j, n_temperatures, from, to
	real(kind=8), allocatable :: densidad(:), temper(:), colis_coef(:), temporal(:)
	real(kind=8) :: C0, gu_gl
	
		allocate(densidad(n_radios))
		allocate(temporal(n_radios))

!		select case(tipo_colision)

!======================================================
! Fixed collisional rates for all the atmosphere
!======================================================
!			case('FIXED')
			if (index(tipo,'FIXED') /= 0 .and. index(tipo,'FIXED_H') == 0 .and. &
				index(tipo,'FIXED_H2') == 0 .and. index(tipo,'FIXED_HE') == 0) then
					
					print *, 'Using fixed collisional rates for all the atmosphere'
					read(3,*) from, to
					do i = from, to
						read(3,*) itran(1,i), itran(2,i), dtran(3,i)
						collision(:,i) = dtran(3,i)
					enddo
			endif

!======================================================
! Fixed collisional cross-section for collision with H
!======================================================
!			case('FIXED_H')
			if (index(tipo,'FIXED_H') /= 0 .and. index(tipo,'FIXED_H2') == 0 .and. &
				index(tipo,'FIXED_HE') == 0) then

				read(3,*) from, to
				print *, 'Using fixed collisional cross-section with H for all the atmosphere'
				
! If we are dealing with an atomic model, then nh(:) is atomic hydrogen density
				if (index(tipo_modelo,'ATOM') /= 0) then
					densidad = nh(:)	

! If we are dealing with a molecular model, then nh(:) is molecular hydrogen density
! and n(H)=xhid*n(H2)
				else
					densidad = xhid * nh(:)
				endif
				
				do i = from, to
					read(3,*) itran(1,i), itran(2,i), dtran(3,i)
! Add to the current collisional rate the calculated here					
					collision(:,i) = collision(:,i) + densidad(:) * dtran(3,i)
				enddo
				
			endif
			
!======================================================
! Fixed collisional cross-section for collision with H2
!======================================================
!			case('FIXED_H2')

			if (index(tipo,'FIXED_H2') /= 0) then	
				
				read(3,*) from, to
				print *, 'Using fixed collisional cross-section with H2 for all the atmosphere'
				
! If we are dealing with an atomic model, then nh(:) is atomic hydrogen density
				if (index(tipo_modelo,'ATOM') /= 0) then
					densidad = nh(:)	

! If we are dealing with a molecular model, then nh(:) is molecular hydrogen density					
				else
					densidad = nh(:)
				endif
				
				do i = from, to
					read(3,*) itran(1,i), itran(2,i), dtran(3,i)					
! Add to the current collisional rate the calculated here										
					collision(:,i) = collision(:,i) + densidad(:) * dtran(3,i)
				enddo
				
			endif


!======================================================
!======================================================
!			case('FIXED_HE')
			if (index(tipo,'FIXED_HE') /= 0) then
				
			endif
			
!======================================================
! Interpoled collisional cross-section for collision with H2
!======================================================
!			case('INTERPOL_H2')
			
			if (index(tipo,'INTERPOL_H2') /= 0 .and. index(tipo,'INTERPOL_H2_SQRTT') == 0) then
								
				read(3,*) from, to
				print *, 'Using interpolated collisional rates with H2'
				
! If we are dealing with an atomic model, then nh(:) is atomic hydrogen density
				if (index(tipo_modelo,'ATOM') /= 0) then
					densidad = nh(:)	

! If we are dealing with an molecular model, then nh(:) is molecular hydrogen density
				else
					densidad = nh(:)
				endif
								
! Read the number of temperatures where we give the value of the collisional rates
				read(3,*) n_temperatures
				allocate(temper(n_temperatures))
				allocate(colis_coef(n_temperatures))
				
				do i = 1, n_temperatures
					read(3,*) temper(i)
				enddo
				
				call lb(3,2)
				do i = from, to
! We read the values of the collisional rates and we interpolate using a spline interpolation				
					read(3,*) itran(1,i), itran(2,i), (colis_coef(j), j=1,n_temperatures)
					option = 'NO  '
					call spline(temper,colis_coef,datph(1,:),temporal,option)
					collision(:,i) = collision(:,i) + densidad(:) * temporal

				enddo
				
			endif
				
!======================================================
! Interpoled collisional cross-section for collision with H2
!======================================================
!			case('INTERPOL_H2_SQRTT')

			if (index(tipo,'INTERPOL_H2_SQRTT') /= 0) then
			
				read(3,*) from, to
				print *, 'Using interpolated collisional rates with H2 with sqrt(T)'
				
! If we are dealing with an atomic model, then nh(:) is atomic hydrogen density
				if (index(tipo_modelo,'ATOM') /= 0) then
					densidad = nh(:)	

! If we are dealing with an molecular model, then nh(:) is molecular hydrogen density
				else
					densidad = nh(:)
				endif
								
! Read the number of temperatures where we give the value of the collisional rates
				read(3,*) n_temperatures
				allocate(temper(n_temperatures))
				allocate(colis_coef(n_temperatures))
				
				do i = 1, n_temperatures
					read(3,*) temper(i)
				enddo
				
				call lb(3,2)
				do i = from, to
! We read the values of the collisional rates and we interpolate using a spline interpolation				
					read(3,*) itran(1,i), itran(2,i), (colis_coef(j), j=1,n_temperatures)
					option = 'NO  '
					call spline(temper,colis_coef,datph(1,:),temporal,option)
					collision(:,i) = collision(:,i) + densidad(:) * temporal
				enddo
				
			endif
								
!======================================================
! Interpoled collisional cross-section for collision with H
!======================================================
!			case('INTERPOL_H')
			if (index(tipo,'INTERPOL_H') /= 0 .and. index(tipo,'INTERPOL_H2') == 0 .and. &
				index(tipo,'INTERPOL_H2_SQRTT') == 0) then
				
			endif

!======================================================
! Interpoled collisional cross-section for collision with He
!======================================================			
!			case('INTERPOL_HE')						

			if (index(tipo,'INTERPOL_HE') /= 0) then
			
				read(3,*) from, to
				print *, 'Using interpolated collisional rates with He'
				
! If we are dealing with an atomic model, then nh(:) is atomic hydrogen density
				if (index(tipo_modelo,'ATOM') /= 0) then
					densidad = nh(:) * xhel	

! If we are dealing with an molecular model, then nh(:) is molecular hydrogen density
				else
					densidad = nh(:) * xhel
				endif
								
! Read the number of temperatures where we give the value of the collisional rates
				read(3,*) n_temperatures
				allocate(temper(n_temperatures))
				allocate(colis_coef(n_temperatures))
				
				do i = 1, n_temperatures
					read(3,*) temper(i)
				enddo
				
				call lb(3,2)
				do i = from, to
! We read the values of the collisional rates and we interpolate using a spline interpolation				
					read(3,*) itran(1,i), itran(2,i), (colis_coef(j), j=1,n_temperatures)
					option = 'NO  '
					call spline(temper,colis_coef,datph(1,:),temporal,option)
					collision(:,i) = collision(:,i) + densidad(:) * temporal			
				enddo
				
			endif
			
!======================================================
! Interpoled collisional cross-section for collision with electrons
!======================================================			
!			case('INTERPOL_EL')

			if (index(tipo,'INTERPOL_EL') /= 0) then
			
				read(3,*) from, to
				print *, 'Using interpolated collisional rates with electrons'
				
! If we are dealing with an atomic model, then nh(:) is atomic hydrogen density
				if (index(tipo_modelo,'ATOM') /= 0) then
					densidad = datph(2,:)

! If we are dealing with an molecular model, then nh(:) is molecular hydrogen density
				else
					densidad = datph(2,:)
				endif
								
! Read the number of temperatures where we give the value of the collisional rates
				read(3,*) n_temperatures
				allocate(temper(n_temperatures))
				allocate(colis_coef(n_temperatures))
				
				do i = 1, n_temperatures
					read(3,*) temper(i)
				enddo
				
				call lb(3,2)
				do i = from, to
! We read the values of the collisional rates and we interpolate using a spline interpolation				
					read(3,*) itran(1,i), itran(2,i), (colis_coef(j), j=1,n_temperatures)
					option = 'NO  '
					call spline(temper,colis_coef,datph(1,:),temporal,option)
					gu_gl = dlevel(2,itran(1,i)) / dlevel(2,itran(2,i))
					collision(:,i) = collision(:,i) + densidad(:) * temporal * gu_gl * sqrt(datph(1,:))
				enddo
				
			endif			

!======================================================			
! Collisional cross-section given by a function depending on the molecule
!======================================================						
!			case('FUNCTION')
			
			if (index(tipo,'FUNCTION') /= 0) then
! Use a function to calculate the collisional rates depending on the molecule			
				
				if (index(atom_mol_nombre,'CO') /= 0) then
! Collisional rates from Ayres & Wiedemann (1989)
						call co_collision(3)
				endif
				
			endif
			
!======================================================			
! Interpolated collisional cross-section for ions
!======================================================			
!			case('OMEGA')

			if (index(tipo,'OMEGA') /= 0) then
			
				print *, 'Using interpolated collisional rates for ions'
				
! If we are dealing with an atomic model, then nh(:) is atomic hydrogen density
				if (index(tipo_modelo,'ATOM') /= 0) then
					densidad = datph(2,:)	
				endif
								
! Read the number of temperatures where we give the value of the collisional rates
				read(3,*) n_temperatures
				allocate(temper(n_temperatures))
				allocate(colis_coef(n_temperatures))
				
				do i = 1, n_temperatures
					read(3,*) temper(i)
				enddo
				
				call lb(3,2)
				do i = 1, nt
! We read the values of the collisional rates and we interpolate using a spline interpolation				
					read(3,*) itran(1,i), itran(2,i), (colis_coef(j), j=1,n_temperatures)
					option = 'NO  '
					call spline(temper,colis_coef,datph(1,:),temporal,option)
					C0 = PION * PI * A0**2 * sqrt( 8.d0 / (PI * PK * PME) ) 
					collision(:,i) = collision(:,i) + C0 * densidad(:) * temporal / &
						( sqrt(datph(1,:)) * dlevel(2,itran(1,i)))					
				enddo
				
			endif
					
		if (allocated(densidad)) deallocate(densidad)
		if (allocated(temporal)) deallocate(temporal)
		if (allocated(temper)) deallocate(temper)
		if (allocated(colis_coef)) deallocate(colis_coef)
	
	end subroutine colisiones

! ----------------------------------------------------------
! Read an atomic or molecular model from a file
! INPUT:
!	   - fich : file where the model is
! OUTPUT:
!     A lot of variables defining the model are changed
! ----------------------------------------------------------
	subroutine atom_model(fich)
	character(len=40), INTENT(IN) :: fich
	integer :: i, j, up, low, ip
	real(kind=8) :: vtherm

		open(unit=3,file=fich,status='old',action='read')
		
		call lb(3,7)

! Read if it is an atom or a molecule
		read(3,*) tipo_modelo
		
		call lb(3,2)
! Read the name of the atom or molecule
		read(3,*) atom_mol_nombre
		if (index(tipo_modelo,'ATOM') /= 0) then
			print *, 'Reading atomic model : ', atom_mol_nombre
		else if (index(tipo_modelo,'MOLECULE') /= 0) then
			print *, 'Reading molecular model : ', atom_mol_nombre
		endif
		
		call lb(3,2)
! Read the mass of the element
		read(3,*) mass
		print *, 'Mass of the element : ', mass, ' amu'
		
		call lb(3,2)
! Read the number of ionization stages
		read(3,*) ni
		allocate(dion(1,ni))			
		print *, 'Number of ionization stages : ', ni
		
		call lb(3,2)
! Read the ionization potentials
		do i = 1, ni
			read(3,*) dion(1,i)
		enddo

		call lb(3,2)		
! Read the number of levels in the model and allocate memory for some data
		read(3,*) nl
		print *, 'Number of levels of the model : ', nl
		allocate(nli(nl))
		allocate(dlevel(2,nl))
		allocate(label(nl))
		
		call lb(3,2)
! Read the number of transitions (including non-radiative transitions) and
! allocate memory for the structures
		read(3,*) nt
		print *, 'Number of transitions : ', nt
		allocate(itran(2,nt))
		allocate(threshold_bf(nt))		
		allocate(dtran(4,nt))
		allocate(doppler_width(n_radios,nt))
		allocate(collision(n_radios,nt))
		
		call lb(3,2)
! Read the number of radiative transitions		
		read(3,*) nr, nbf
		print *, 'Number of bound-bound radiative transitions : ', nr
		print *, 'Number of bound-free radiative transitions : ', nbf
				
		call lb(3,2)
! Read the number of active transitions
		read(3,*) nact
		print *, 'Number of active transitions : ', nact
		
! Allocate memory for Jbar20
		allocate(Jbar20(nact, n_radios))
				
		call lb(3,2)
! Read the energy levels, degeneracy and the ionization stage they belong		
		print *, 'Reading energy levels...'
		do i = 1, nl
			read(3,*) dlevel(1,i), dlevel(2,i), nli(i), label(i)
		enddo
		

! If it is a molecule, energy levels are given in cm^-1, then we have to convert it to Hz
		if (index(tipo_modelo,'MOLECULE') /= 0) then
			dlevel(1,:) = PC * dlevel(1,:)
		endif
		
		call lb(3,2)
! Read radiative transitions (upper level, lower level, Aij)		
		print *, 'Reading radiative transitions and Einstein coefficients'
		do i = 1, nr
			read(3,*) itran(1,i), itran(2,i), dtran(1,i)
		enddo

! Calculate the frequency of each active transition
		do i = 1, nr
			dtran(2,i) = dabs(dlevel(1,itran(1,i)) - dlevel(1,itran(2,i)))
		enddo

! Using Aij, precalculate several things : c^2/(8*pi*nu^2)*(gu/gl)*Aul = h*nu/(4*pi)*Blu
		do i = 1, nr
			up = itran(1,i)
			low = itran(2,i)
			dtran(1,i) = dtran(1,i) / (PI8C2 * dtran(2,i)**2.d0) * (dlevel(2,up) / dlevel(2,low))
		enddo
		
		call lb(3,2)
		
! Read bound-free transitions
		print *, 'Reading bound-free transitions...'
		do i = nr+1, nr+nbf
			print *, 'hey'
			read(3,*) itran(1,i), itran(2,i), nfrq_bf, threshold_bf(i)

! We assume that all the cross-sections are given with the same number of wavelenghts
! bf_cross_section(1) is frequency
! bf_cross_section(1) is cross_section
! bf_cross_section(1) is trapezoidal weight
			if (i == nr+1) then
				allocate(bf_cross_section(3,nt,nfrq_bf))
			endif
! Read the cross-section			
			do j = 1, nfrq_bf
				read(3,*) bf_cross_section(1,i,j), bf_cross_section(2,i,j)

! Frequency
				bf_cross_section(1,i,j) = PC / (1.d-8*bf_cross_section(1,i,j))
			enddo
			
! Trapezoidal weights				
			bf_cross_section(3,i,1) = 0.5d0 * (bf_cross_section(1,i,2)-bf_cross_section(1,i,1))
			bf_cross_section(3,i,nfrq_bf) = 0.5d0 * (bf_cross_section(1,i,nfrq_bf)-bf_cross_section(1,i,nfrq_bf-1))
			do j = 2, nfrq_bf-1
				bf_cross_section(3,i,j) = 0.5d0 * (bf_cross_section(1,i,j+1)-bf_cross_section(1,i,j-1))
			enddo
			
! alpha' = alpha / (h*nu)
			bf_cross_section(2,i,:) = bf_cross_section(2,i,:) / (PH*bf_cross_section(1,i,:))
			
			do j = 1, nfrq_bf
!				print *, bf_cross_section(1,i,j), bf_cross_section(2,i,j), bf_cross_section(3,i,j)
			enddo
			
			call lb(3,1)
		enddo
		
! Read collisional rates		
		read(3,*) tipo_colision
		collision = 0.d0

! Read all the possible kinds of collisions until END is reached
		do while (index(tipo_colision,'END') /= 1)
			call colisiones(tipo_colision)		
			call lb(3,1)
			read(3,*) tipo_colision
		enddo

! Calculate Doppler width in each shell and each transition
		do ip = 1, n_radios

! Vtherm = sqrt(2*K*T/M + v_micro**2)
	      vtherm = dsqrt(2.d0 * PKUMA * datph(1,ip) / mass + (datph(4,ip) * 1.d5)**2.d0)		
			thermal_v(ip) = vtherm
			
! Transform velocities to Doppler units (first transform to cgs)
			datph(7,ip) = 1.d5 * datph(7,ip) / vtherm			
			
! Fill the Doppler width delta_nu=nu_0 * v/c 
      	do i=1,nr
      		dtran(4,i) = dtran(2,i) * vtherm / PC
				doppler_width(ip,i) = dtran(2,i) * vtherm / PC
			enddo				
			
		enddo
		
! Read partition function if it is a molecular model
		if (index(tipo_modelo,'MOLECULE') /= 0) then
			
			call lb(3,2)
			read(3,*) n_partition

! If n_partition=0, then build up partition function using energy levels

			if (n_partition /= 0) then
				print *, 'Reading coefficients for the partition function'
				
! Coefficients for the fit (Sauval & Tatum (1984))
				allocate(partition(n_partition))

				do i = 1, n_partition
					read(3,*) partition(i)
				enddo
			endif
			
		endif		
		
		close(3)	
	end subroutine atom_model

end module atomic
