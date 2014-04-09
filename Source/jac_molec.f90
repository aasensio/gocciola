! *********************************************************
! *********************************************************
! *********************************************************
! MAIN PROGRAM
! *********************************************************
! *********************************************************
! *********************************************************

! This code solves the non-LTE problem for an atomic or molecular model using Jacobi, Gauss-Seidel
! or Successive Overrelaxation iterative schemes in a spherically symmetric atmosphere
! ANDRES ASENSIO RAMOS 
! Last revision (21/02/2001)

program jacobi
use general
use initial
use formal_jacobi
use formal_jacobi_overlap
use formal_jacobi_pp
use formal_jacobi_pp_overlap
use formal_cep
use formal_gs
use formal_gs_overlap
use formal_gs_pp
use formal_gs_pp_overlap
use statis
use atmosfer
use accel
use constantes
use variables
use lvg_func
use data_write
use dealloc
use emerging
use zeeman
implicit none
	real :: anterior, total, actual
	real(kind=8) :: estimation, conv_error

   
! Initialize everything
   call init
	
! Write the model atmosphere
	call write_atmos
	

! If needed, read the population from a previous run
	if (converge_flag == 3) then
		call read_pops_true_error(pop)

! Then do normal iteration		
		converge_flag = 0
	endif
		
	if (converge_flag == 0) then

	! Open iteration file	
		open (UNIT=4,FILE=nombre_iteration,STATUS='replace',ACTION='write')	
		
	! LVG populations if flag is active
		if (lvg_flag /= 0.d0) then
			anterior = second()
			call lvgpops
			actual = second()
			print *, 'Time for LVG initialization : ', actual-anterior, ' s'
		endif

		relat_error = 1.d10
		iter = 1
		anterior = second()
		total = 0.d0

! ***********************************************************
! ***********************************************************
! JACOBI
! ***********************************************************
! ***********************************************************
!		if (tipo_iteracion == 'JACOBI' .or. tipo_iteracion =='LAMBDA' .or. tipo_iteracion =='LAMBDAJAC') then
		if (tipo_iteracion == 3 .or. tipo_iteracion == 1 .or. tipo_iteracion == 2) then

! Iterate until reaching desired precision or maximum number of iterations
			do while (relat_error > precision .and. iter <= iter_max)

! Calculate Jbar_total and Lstar_total using a formal solution of the RTE 
! based on short-characteristics
				print *, 'Formal solver'
				if (index(geometry_type,'SPHERICAL') /= 0) then
					if (overlap_flag /= 0) then
						call formal_solver_jacobi_overlap
					else
						call formal_solver_jacobi
					endif
				else
					if (overlap_flag /= 0) then				
						call formal_solver_jacobi_pp_overlap
					else
						call formal_solver_jacobi_pp
					endif
				endif
				
! Make the correction for the population using preconditioning by Rybicki & Hummer (1991)
				print *, 'Population correction'
      		call corrige_poblaciones
				print *, 'Maximum change at level ',level_highest_change

				iter = iter + 1
				itracc = itracc + 1

! Make acceleration if needed
				print *, 'Acceleration test'
				call aceleracion

				actual = second()
				total = total + actual - anterior
				
				if (error_flag == 0) then
					print *, 'Iteration ',iter-1, ' ---  Relative error : ',relat_error
				else
					print *, 'Iteration ',iter-1, ' ---  Relative error : ', relat_error									
					print *, '              ---  True error :     ', true_error
				endif
				
				estimation = relat_error / relat_error_p
				if (estimation < 1.d0) then 
					print *, 'Relative error curve slope : ', estimation
					conv_error = relat_error * estimation / (1.d0 - estimation)
					print *, 'Convergence error : ',conv_error
				endif
				print *, 'Iteration time : ', actual - anterior, ' s  / Total time : ', total, ' s'
				print *
				call write_iteration!(actual-anterior)				
				anterior = second()
				
!				if (tipo_iteracion == 'LAMBDAJAC' .and. relat_error < 0.1d0) then
				if (tipo_iteracion == 2 .and. relat_error < 0.1d0) then
!					tipo_iteracion = 'JACOBI'
					tipo_iteracion = 3
					print *, 'Nuevo ',tipo_iteracion
					print *, 'Switched to Jacobi'
				endif

! Calculate optical depth in each transition
				call calc_tau
			enddo

		endif

! ***********************************************************
! ***********************************************************
! GAUSS-SEIDEL OR SUCCESSIVE OVERRELAXATION
! ***********************************************************
! ***********************************************************

!		if (tipo_iteracion == 'GS' .or. tipo_iteracion == 'SOR') then
		if (tipo_iteracion == 4 .or. tipo_iteracion == 5) then


! Iterate until reaching desired precision or maximum number of iterations
			do while (relat_error > precision .and. iter <= iter_max)

! Make a GS or SOR iteration
				if (index(geometry_type,'SPHERICAL') /= 0) then
					if (overlap_flag /= 0) then
						call formal_solver_gs_overlap(1)
					else
						call formal_solver_gs(1)
					endif
				else
					if (overlap_flag /= 0) then				
						call formal_solver_gs_pp_overlap(1)
					else
						call formal_solver_gs_pp(1)
					endif
				endif

				iter = iter + 1
				itracc = itracc + 1

! Make acceleration if needed
				call aceleracion

				actual = second()
				total = total + actual - anterior
				
				if (error_flag == 0) then
					print *, 'Iteration ',iter-1, ' ---  Relative error : ', relat_error
				else
					print *, 'Iteration ',iter-1, ' ---  Relative error : ', relat_error					
					print *, '              ---  True error :     ', true_error
				endif
				
				estimation = relat_error / relat_error_p
				if (estimation < 1.d0) then 
					print *, 'Relative error curve slope : ', estimation
					print *, 'Omega estimation : ', 2.d0 / (1.d0 + dsqrt(1.d0 - estimation))
					conv_error = relat_error * estimation / (1.d0 - estimation)
					print *, 'Convergence error : ',conv_error
				endif
				print *, 'Iteration time : ', actual - anterior, ' s  / Total time : ', total, ' s'
				print *
				call write_iteration!(actual-anterior)				
				anterior = second()

! Calculate optical depth in each transition
				call calc_tau

			enddo

		endif	
		
! ***********************************************************
! ***********************************************************
! SEMI GAUSS-SEIDEL OR SEMI-SUCCESSIVE OVERRELAXATION
! ***********************************************************
! ***********************************************************

!		if (tipo_iteracion == 'SEMIGS' .or. tipo_iteracion == 'SEMISOR') then
		if (tipo_iteracion == 6 .or. tipo_iteracion == 7) then


! Iterate until reaching desired precision or maximum number of iterations
			do while (relat_error > precision .and. iter <= iter_max)

! Make a GS or SOR iteration
				if (index(geometry_type,'SPHERICAL') /= 0) then
					if (overlap_flag /= 0) then
						call formal_solver_gs_overlap(0)
					else
						call formal_solver_gs(0)
					endif
				else
					if (overlap_flag /= 0) then				
						call formal_solver_gs_pp_overlap(0)
					else
						call formal_solver_gs_pp(0)
					endif
				endif

				iter = iter + 1
				itracc = itracc + 1

! Make acceleration if needed
				call aceleracion

				actual = second()
				total = total + actual - anterior
				if (error_flag == 0) then
					print *, 'Iteration ',iter-1, ' ---  Relative error : ', relat_error
				else
					print *, 'Iteration ',iter-1, ' ---  Relative error : ', relat_error									
					print *, '              ---  True error :     ', true_error
				endif
				estimation = relat_error / relat_error_p
				if (estimation < 1.d0) then 
					print *, 'Relative error curve slope : ', estimation
					print *, 'Omega estimation : ', 2.d0 / (1.d0 + dsqrt(1.d0 - estimation))
					conv_error = relat_error * estimation / (1.d0 - estimation)
					print *, 'Convergence error : ',conv_error
				endif
				print *, 'Iteration time : ', actual - anterior, ' s  / Total time : ', total, ' s'
				print *
				call write_iteration!(actual-anterior)				
				anterior = second()

! Calculate optical depth in each transition
				call calc_tau

			enddo

		endif			
		
! ***********************************************************
! ***********************************************************
! NEWTON METHOD WITH THE CEP METHOD (Elitzur & Asensio Ramos MNRAS, 2006)
! ***********************************************************
! ***********************************************************

		if (tipo_iteracion == 8) then

! Precalculation of the beta function   
!   		call precalculate_beta
			call mnewt(20,pop,n_radios*nl,1.d-4,1.d-4,.TRUE.)

		endif		


! Write the final results	
		print *, 'Writing final results...'
		call write_results

! Write radiative rates
		print *, 'Writing radiative rates...'
		call write_radrates	

! Write collisional rates
		print *, 'Writing collisional rates...'
		call write_collisrates

! Write background opacity
		print *, 'Writing background opacity...'
		call write_background_opa

! Close iteration file
		close(4)	

	endif

! Read populations from previous calculation and write spectrum	
	if (converge_flag == 1) then
		call read_populations
		call write_results
		!call emerge
		print *, 'Writing atmosphere...'
		call write_atmos
		print *, 'Writing radiative rates...'
		call write_radrates
		print *, 'Writing background opacity...'
		call write_background_opa
	endif
	
! Synthesize LTE spectrum
	if (converge_flag == 2) then
		call emerge
	endif	

! Zeeman synthesis of a line
	if (converge_flag == 4) then
		call read_pops_true_error(pop)
		call zeeman_synth
	endif	

! Clean all variables	
	print *, 'Cleaning all the memory allocated...'
	call clean
	
end program jacobi
