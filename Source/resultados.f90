module data_write
use variables
use emerging
use general
implicit none
contains
! *********************************************************
! *********************************************************
! RESULTS OUTPUT ROUTINES
! *********************************************************
! *********************************************************

!-----------------------------------------------------------------	
! Write final results in a file
!-----------------------------------------------------------------
	subroutine write_results
	integer :: i, ip, ipl, ipt, it, up, low, min_tau
	real(kind=8) :: sig, glu, acon, agnus, chilp, chilm, snlte
	real(kind=8), allocatable :: slte(:), sb(:), tem(:), opac(:)

		allocate(slte(n_radios))
		allocate(sb(n_radios))
		allocate(tem(n_radios))
		allocate(opac(n_radios))
	
		open (UNIT=3,FILE=nombre_salida,STATUS='replace',ACTION='write')	
	
	
		write(3,*) '-----------------------------------------------------'
		write(3,*) '                GENERAL DATA                         '
		write(3,*) '-----------------------------------------------------'
		write(3,*) 'N. active transitions'
		write(3,*) nact
		write(3,*) 'N. levels'
		write(3,*) nl
		write(3,*) 'N. grid points'
		write(3,*) n_radios
		write(3,*) 'N. characteristics'
		write(3,*) n_caract
		write(3,*) 'N. frequencies'
		write(3,*) nfrq
		write(3,*) 'Output spectrum present'
		write(3,*) output_spectrum		
		write(3,*) 'Background opacity included'
		write(3,*) back_opac_flag
		write(3,*) 'Radiative rates included'
		write(3,*) radrates_flag
		write(3,*) 'Atomic/molecular model used'
		write(3,*) nombre_modelo
		write(3,*) 'Transitions (upper, lower, frec)'
		do i = 1, nact
			write(3,FMT='(I5, I5, E12.5)') itran(1,i), itran(2,i), dtran(2,i)
		enddo
		
		write(3,*) '-----------------------------------------------------'
		write(3,*) '                   RADIAL GRID                       '
		write(3,*) '-----------------------------------------------------'
		
		do ip = 1, n_radios
			write(3,FMT='(1P5E12.5)') r(ip)
		enddo
		
		write(3,*) '-----------------------------------------------------'
		write(3,*) '                ANGULAR GRID                         '
		write(3,*) '-----------------------------------------------------'
		
		if (index(geometry_type,'PLANEP') /= 0) then
			do ip = 1, n_caract
				write(3,FMT='(1P5E12.5)') mu_pp(ip), wtmu(ip)
			enddo
		else
			do ip = 1, n_caract
				write(3,FMT='(1P5E12.5)') mu(n_radios,ip)
			enddo
		endif
		
		write(3,*) '-----------------------------------------------------'
		write(3,*) '                 ENERGY LEVELS                       '
		write(3,*) '-----------------------------------------------------'		
		
		do ip = 1, nl
			write(3,FMT='(1P5E12.5)') dlevel(1,ip)
		enddo
		
		write(3,*) '-----------------------------------------------------'
		write(3,*) '             FINAL POPULATIONS                       '
		write(3,*) '-----------------------------------------------------'
		
		do ip = 1, n_radios
			ipl = nl * (ip-1)
			write(3,FMT='(I3, (1P5E13.5))') ip, (pop(i+ipl), i = 1, nl)
		enddo
		
		write(3,*)
		write(3,*) '*****************************************************'
		write(3,*)
		
		write(3,*) '-----------------------------------------------------'
		write(3,*) '        FINAL DEPARTURE COEFFICIENTS                 '
		write(3,*) '-----------------------------------------------------'
		
		do ip = 1, n_radios
			ipl = nl * (ip-1)
			write(3,FMT='(I3, (1P5E13.5))') ip, (pop(i+ipl)/popl(i+ipl), i = 1, nl)
		enddo
		
		write(3,*)
		write(3,*) '*****************************************************'
		write(3,*)
		
		write(3,*) '-----------------------------------------------------'
		write(3,*) '        INITIAL POPULATION (LTE or LVG)              '
		write(3,*) '-----------------------------------------------------'
		
		do ip = 1, n_radios
			ipl = nl * (ip-1)
			write(3,FMT='(I3, (1P5E13.5))') ip, (popi(i+ipl), i = 1, nl)
		enddo		
		
		write(3,*)
		write(3,*) '*****************************************************'
		write(3,*)
		
		write(3,*) '-----------------------------------------------------'
		write(3,*) '                    FINAL LSTAR                      '
		write(3,*) '-----------------------------------------------------'
		
		do ip = 1, n_radios
			ipt = nact * (ip-1)
			write(3,FMT='(I3, (1P5E13.5))') ip, (lstar_total(i+ipt), i = 1, nact)
		enddo
		
		write(3,*)
		write(3,*) '*****************************************************'
		write(3,*)
		
		write(3,*) '-----------------------------------------------------'
		write(3,*) '                    FINAL JBAR                       '
		write(3,*) '-----------------------------------------------------'
		
		do ip = 1, n_radios
			ipt = nact * (ip-1)
			write(3,FMT='(I3, (1P5E13.5))') ip, (Jbar_total(i+ipt), i = 1, nact)
		enddo
		
		write(3,*)
		write(3,*) '*****************************************************'
		write(3,*)
		
		write(3,*) '-----------------------------------------------------'
		write(3,*) '                  FINAL JBAR20                       '
		write(3,*) '-----------------------------------------------------'
		
		do ip = 1, n_radios
			write(3,FMT='(I3, (1P5E13.5))') ip, (Jbar20(i,ip), i = 1, nact)
		enddo
		
		write(3,*)
		write(3,*) '*****************************************************'
		write(3,*)
		
		write(3,*) '-----------------------------------------------------'
		write(3,*) '                  SOURCE FUNCTIONS                   '
		write(3,*) '-----------------------------------------------------'
		tau = 0.d0
		do it = 1, nact
			tau(it,n_radios) = 0.d0
			up	= itran(1,it)
			low = itran(2,it)
			glu = dlevel(2,low) / dlevel(2,up)
			acon = PHC * (dtran(2,it)) ** 3   !2*h*nu^3/c^2
			agnus = acon * glu
			call opacity(it,min_tau,1,n_radios)

			chilm = 0.d0			
			do i = n_radios, 1, -1
				sig = dtran(1,it) / doppler_width(i,it)
				ipl = nl * (i-1)
				chilp = chil(i) + kappa(i)
				opac(i) = chilp
				if (i /= n_radios) tau(it,i) = tau(it,i+1) + 0.5d0 * (chilm + chilp) * dabs(r(i+1) - r(i))
				chilm = chilp
				snlte = agnus * pop(up+ipl) / (pop(low+ipl) - glu*pop(up+ipl))
				slte(i) = agnus * popl(up+ipl) / (popl(low+ipl) - glu*popl(up+ipl))
				sb(i) = snlte / slte(i)
				if ((1.d0 + (PHC * dtran(2,it)**3 / snlte)) > 0.d0) then
					tem(i) = PHK * dtran(2,it) / log(1.d0 + (PHC * dtran(2,it)**3 / snlte))
				else
					tem(i) = 0.d0
					!print *, 'Hey'
				endif
				call calcula_centro(it)
			enddo
			
			write(3,*) '*****************************************************'
			write(3,*) '    TRANSITION IT=', it
			write(3,*) '         ', up, ' -> ',low
			write(3,*) ' lambda : ', PC / dtran(2,it) * 1.d8,' angstrom'
			write(3,*) ' freq : ', dtran(2,it), ' Hz'
			write(3,*) '    PLANCK FUNCTION AT SURFACE = ', slte(n_radios)
			write(3,*) ' Equivalent internal boundary condition temperature = ', T_planck_equiv(dtran(2,it),dentro(1))
			write(3,*) '     tau          S/B         B        T_exc      Opacity'
			do i = 1, n_radios
				write(3,FMT='(1P5E12.3)') tau(it,i), sb(i), slte(i), tem(i), opac(i)
			enddo
			
			write(3,*)		
		enddo
		
		write(3,*) 'Output spectrum present'
		write(3,*) output_spectrum
		
		close(3)	
		
		if (output_spectrum == 1) then 
			print *, 'Writing emerging spectrum...'
			call emerge
		endif
		
! Save another copy of the final results for later use by the true error version
! It includes only the radial grid and the populations
		if (error_flag == 0) then
			open (UNIT=3,FILE=trim(adjustl(nombre_salida))//'.final',STATUS='replace',ACTION='write')	
	
	
			write(3,*) '-----------------------------------------------------'
			write(3,*) '                GENERAL DATA                         '
			write(3,*) '-----------------------------------------------------'
			write(3,*) 'N. levels'
			write(3,*) nl
			write(3,*) 'N. grid points'
			write(3,*) n_radios
			
			write(3,*) '-----------------------------------------------------'
			write(3,*) '                   RADIAL GRID                       '
			write(3,*) '-----------------------------------------------------'

			do ip = 1, n_radios
				write(3,FMT='(1P5E21.15)') r(ip)
			enddo
		
			write(3,*) '-----------------------------------------------------'
			write(3,*) '             FINAL POPULATIONS                       '
			write(3,*) '-----------------------------------------------------'

			do ip = 1, n_radios
				ipl = nl * (ip-1)
				write(3,FMT='(I3, (1P5E13.5))') ip, (pop(i+ipl), i = 1, nl)
			enddo
			
			close(3)
		endif

		deallocate(slte)
		deallocate(sb)
		deallocate(tem)
		deallocate(opac)

	end subroutine write_results

!-----------------------------------------------------------------	
! Write atmosphere to a file
!-----------------------------------------------------------------
	subroutine write_atmos
	integer :: i, j
	
		open (UNIT=3,FILE=nombre_atmosf_out,STATUS='replace',ACTION='write')	
		
		write(3,*) '-----------------------------------------------------'
		write(3,*) '                GENERAL DATA                         '
		write(3,*) '-----------------------------------------------------'		
		write(3,*) 'N. grid points'
		write(3,*) n_radios
		write(3,*) 'N. characteristics'
		write(3,*) n_caract		
		
		write(3,*)
		write(3,*) '-----------------------------------------------------'
		write(3,*) '                 RADIAL VARIABLES                    '
		write(3,*) '-----------------------------------------------------'		
		write(3,*) '  r    T   n(H2)    n   n_e   v_mic   T_dust  v_mac  A_v'
		do i = 1, n_radios
			write(3,'(1P9E13.5)') r(i), datph(1,i), nh(i), datph(3,i), datph(2,i), datph(4,i), &
				datph(6,i), datph(7,i)*thermal_v(i)/1.d5, datph(5,i)
		enddo
		write(3,*)
		write(3,*) '-----------------------------------------------------'
		write(3,*) '                IMPACT PARAMETERS                    '
		write(3,*) '-----------------------------------------------------'		
		do i = 1, n_caract					
			write(3,'(1P2E13.5)') p(i), mu(n_radios,i)	
		enddo
		
		close(3)
		
		
		open (UNIT=3,FILE=trim(adjustl(nombre_atmosf_out))//'.complete',STATUS='replace',ACTION='write')	
		
		write(3,*) '-----------------------------------------------------'
		write(3,*) '                GENERAL DATA                         '
		write(3,*) '-----------------------------------------------------'		
		write(3,*) 'N. grid points'
		write(3,*) n_radios
		write(3,*) 'N. characteristics'
		write(3,*) n_caract		
		
		write(3,*)
		write(3,*) '-----------------------------------------------------'
		write(3,*) '                 RADIAL VARIABLES                    '
		write(3,*) '-----------------------------------------------------'		
		write(3,*) '  r    T   n(H2)    n   n_e   v_mic   T_dust  v_mac  A_v'
		do i = 1, n_radios
			write(3,'(1P9E13.5)') r(i), datph(1,i), nh(i), datph(3,i), datph(2,i), datph(4,i), &
				datph(6,i), datph(7,i)*thermal_v(i)/1.d5, datph(5,i)
		enddo
		write(3,*)
		write(3,*) '-----------------------------------------------------'
		write(3,*) '                IMPACT PARAMETERS                    '
		write(3,*) '-----------------------------------------------------'		
		do j = 1, n_radios
			write(3,FMT='(A,1X,I4,1X,A)') ' ----- Shell', j, '-----'
			do i = 1, n_caract			
				if (mu(j,i) <= 1.d0) then
					write(3,'(1P2E13.5)') p(i), mu(j,i)
				endif
			enddo			
		enddo
		
		close(3)
		
	end subroutine write_atmos

!-----------------------------------------------------------------	
! Write iteration step
!-----------------------------------------------------------------	
	subroutine write_iteration
		if (error_flag == 0) then
			write(4,*) iter, relat_error
		else
			write(4,*) iter, true_error
		endif
	end subroutine write_iteration
	
!-----------------------------------------------------------------	
! Write background opacity if needed
!-----------------------------------------------------------------	
	subroutine write_background_opa
	integer :: it, i
		if (back_opac_flag == 1) then
			open (UNIT=3,FILE=nombre_background_out,STATUS='replace',ACTION='write')	

			do it = 1, nact
				do i = 1, n_radios
					write(3,*) chic(it,i)
				enddo
			enddo
			
			close(3)
		endif
		
	end subroutine write_background_opa
	
!-----------------------------------------------------------------	
! Write radiative rates for all active transitions (upward and downward)
!-----------------------------------------------------------------	
	subroutine write_radrates
	integer :: it, i, ipt, up, low
	real(kind=8) :: Rul, Rlu, glu, blu, bul, aul, acon
	
		if (radrates_flag == 1) then
			open (UNIT=3,FILE=nombre_radrates,STATUS='replace',ACTION='write')	

			do it = 1, nact
				up = itran(1,it)
				low = itran(2,it)
				glu = dlevel(2,low) / dlevel(2,up)
				acon = PHC * (dtran(2,it))**3
				blu = PI4H * dtran(1,it) / dtran(2,it)
				bul = glu * blu
				aul = acon * bul
				do i = 1, n_radios
					ipt = nr*(i-1)					
! Emmision rate					
					Rul =  aul + bul * Jbar_total(it+ipt)
! Absorption rate
					Rlu = blu * Jbar_total(it+ipt)
					write(3,*) Rul, Rlu
				enddo
			enddo
			
			close(3)
		endif
		
	end subroutine write_radrates


!-----------------------------------------------------------------
! Write collisional rates for all transitions (only deexcitation)
!-----------------------------------------------------------------
	subroutine write_collisrates
	integer :: it, ip
		open(UNIT=3,FILE=nombre_collisrates,STATUS='replace',ACTION='write')
		write(3,*) 'First shell collisions'
		do it = 1, nt 		
			write(3,*) itran(1,it), itran(2,it), collision(1,it)/nh(1)
		enddo
		do it = 1, nt 
			do ip = 1, n_radios
				write(3,*) itran(1,it), itran(2,it), collision(ip,it)/nh(ip)
			enddo
		enddo
		close(3)
	end subroutine write_collisrates
	
!-----------------------------------------------------------------	
! Read populations from file
!-----------------------------------------------------------------	
	subroutine read_populations	
	integer :: int1, i, ip, ipl, temp, n_radios_read, n_caract_read
	real(kind=8) :: b(nl)
	character*80 :: t1
	real(kind=8), allocatable :: temp_anterior(:), temp_nuevo(:), eje_anterior(:), eje_nuevo(:), pop_new(:)

! If we are going to interpolate the atmosphere, there are n_radios_old and n_caract_old data
! After that, we have to interpolate
		print *, 'Reading populations...'
		n_radios_read = n_radios
		n_caract_read = n_caract

		if (interpolate_atmosphere /= 0) then
			n_radios_read = n_radios_old
			n_caract_read = n_caract_old
		endif

		open (UNIT=3,FILE=nombre_salida,STATUS='old',ACTION='read')
		call lb(3,22)

		call lb(3,nact)
		call lb(3,3)
		
		call lb(3,n_radios_read)
		call lb(3,3)
		
		call lb(3,n_caract_read)
		call lb(3,3)
		
		call lb(3,nl)
		call lb(3,3)
				
		do ip = 1, n_radios_read
			ipl = nl * (ip-1)
			read(3,*) int1, (pop(i+ipl), i = 1, nl)
		enddo
		
		call lb(3,6)
				
		do ip = 1, n_radios_read
			ipl = nl * (ip-1)
			read(3,FMT='(I3, (1P5E13.5))') int1, (b(i), i = 1, nl)
			do i = 1, nl
				popl(i+ipl) = pop(i+ipl) / b(i)
			enddo
		enddo
		
		close(3)		
		
! If spectrum has been calculated, then indicate it in the results file		
		open (UNIT=3,FILE=nombre_salida,STATUS='old',ACTION='read')
		open (UNIT=4,FILE='caca.temporal',STATUS='replace',ACTION='readwrite')			
		
		temp = 1
		do while(temp > 0)
			read(3,FMT='(A80)',END=12) t1
			write(4,*) trim(adjustl(t1))
		enddo

12		close(3)
		rewind(4)
		open (UNIT=3,FILE=nombre_salida,STATUS='replace',ACTION='write')			
		
		temp = 1
		do while(temp > 0)
			read(4,FMT='(A80)',END=13) t1
			if (temp == 15) t1 = '1'
			write(3,*) trim(adjustl(t1))
			temp = temp + 1
		enddo
		
13		close(4)
		close(3)

! If we are interpolating the atmosphere, do it for pop and popl
		allocate(eje_anterior(n_radios_old))
		allocate(eje_nuevo(n_radios))
		allocate(temp_anterior(n_radios_old))
		allocate(temp_nuevo(n_radios))
		allocate(pop_new(nl*n_radios))

		do i = 1, n_radios_old
			eje_anterior(i) = 1.d0 * i / n_radios_old
		enddo
		do i = 1, n_radios
			eje_nuevo(i) = 1.d0 * i / n_radios
		enddo
		
! Interpolate pop. For each level, we interpolate through the atmosphere
		do i = 1, nl
			do ip = 1, n_radios_old
				ipl = nl * (ip-1)
				temp_anterior(ip) = pop(i+ipl)
			enddo
			
!			call spline(eje_anterior,temp_anterior,eje_nuevo,temp_nuevo,'NO  ')
			call linear_interpol(eje_anterior,temp_anterior,eje_nuevo,temp_nuevo)
			do ip = 1, n_radios
				ipl = nl * (ip-1)
				pop_new(i+ipl) = temp_nuevo(ip)
			enddo
		enddo

		deallocate(pop)
		allocate(pop(nl*n_radios))
		pop = pop_new

! Interpolate popl.		
		do i = 1, nl
			do ip = 1, n_radios_old
				ipl = nl * (ip-1)
				temp_anterior(ip) = popl(i+ipl)
			enddo

!			call spline(eje_anterior,temp_anterior,eje_nuevo,temp_nuevo,'NO  ')
			call linear_interpol(eje_anterior,temp_anterior,eje_nuevo,temp_nuevo)
			
			do ip = 1, n_radios
				ipl = nl * (ip-1)
				pop_new(i+ipl) = temp_nuevo(ip)
			enddo
		enddo

		deallocate(popl)
		allocate(popl(nl*n_radios))
		popl = pop_new

		deallocate(eje_anterior)
		deallocate(eje_nuevo)
		deallocate(temp_anterior)
		deallocate(temp_nuevo)
		deallocate(pop_new)
			
	end subroutine read_populations
	
!-----------------------------------------------------------------	
! Read populations from file for using the true error
!-----------------------------------------------------------------	
	subroutine read_pops_true_error(salida)	
	real(kind=8), INTENT(INOUT) :: salida(nl*n_radios)
	integer :: int1, i, ip, ipl, n_r	
	real(kind=8), allocatable :: radial_grid(:), pop_temp(:), pop_temp2(:), pop_temp3(:)

		open (UNIT=3,FILE=trim(adjustl(nombre_salida))//'.final',STATUS='old',ACTION='read')	

		call lb(3,6)	
		read(3,*) n_r
		allocate(radial_grid(n_r))
		allocate(pop_temp(nl*n_r))
		allocate(pop_temp2(n_r))
		allocate(pop_temp3(n_radios))

		call lb(3,3)	

! Read the radial grid. It can be finer or coarser than the one we are calculating
		do ip = 1, n_r
			read(3,*) radial_grid(ip)
		enddo

		call lb(3,3)

		do ip = 1, n_r
			ipl = nl * (ip-1)
			read(3,*) int1, (pop_temp(i+ipl), i = 1, nl)
		enddo		
				
! In case the number of radial shells is different from the atmosphere we are calculating, then
! do an interpolation

		do ipl = 1, nl
! Read the population in each point of the old atmosphere (with different radial grid)
			do ip = 1, n_r
				pop_temp2(ip) = pop_temp(ipl+nl*(ip-1))
			enddo

! Interpolate to the new grid
			call spline(radial_grid,pop_temp2,r,pop_temp3,'NO  ')

! Put this interpolated values into the final population
			do ip = 1, n_radios
				salida(ipl+nl*(ip-1)) = pop_temp3(ip)
			enddo
			
		enddo
				
		close(3)
		
		deallocate(radial_grid)
		deallocate(pop_temp)
		deallocate(pop_temp2)		
		deallocate(pop_temp3)
		
	end subroutine read_pops_true_error
	
end module data_write
