module emerging
use variables
use short_ch
use formal_jacobi
use formal_jacobi_overlap
use formal_jacobi_pp
use formal_jacobi_pp_overlap
use general
implicit none
contains

! *********************************************************
! *********************************************************
! ROUTINES FOR EMERGING SPECTRUM
! *********************************************************
! *********************************************************

!-----------------------------------------------------------------
! This routine returns the sigma in microns of the selected instrument for a given 
! wavelength
!-----------------------------------------------------------------
	function instrum_sigma(wavelength)
   real(kind=8) :: instrum_sigma, wavelength, sigma
   
   	sigma = 0.d0
   	select case(observing_instrument)

      	case('LWS')
         	
            if (wavelength < 90.d0) then
            	sigma = 0.29d0
            else
            	sigma = 0.6d0
            endif
            
         case('SWS')
         
      end select
      
      instrum_sigma = sigma
      
   end function instrum_sigma

!-----------------------------------------------------------------	
! Convolve spectrum with the instrument's response
!-----------------------------------------------------------------	
	subroutine convolve_spectrum2(flux,output)
   real(kind=8) :: flux(nact,nfrq), output(nact,nfrq), sigma2, frecu0, fwhm
   real(kind=8) :: step
   real(kind=8), allocatable :: signal(:), filter(:), result(:), frec(:)
	integer :: f1, it, j, delta, n
   
   	output = 0.d0
		do it = 1, nr
			frecu0 = dtran(2,it)

! Calculate the FWHM of the instrument in Hz
			fwhm = frecu0**2 / PC * (instrum_sigma(PC * 1.d4 / frecu0) * 1.d-4)

! Calculate the sigma^2 of the telescope at this frequency         
         sigma2 = fwhm**2 / (4.d0 * alog(2.e0))
         step = doppler_width(n_radios,it) / n_freq_ud_doppler
         n = int(3.d0 * dsqrt(sigma2) / step)
         
         allocate(signal(2*n))
         allocate(filter(2*n))
         allocate(result(2*n))
         allocate(frec(2*n))
         delta = n - nfrq/2
         signal(n-delta:n+delta) = flux(it,:)   
         do j = 1, 2*n
         	frec(j) = frecu0 + (j-n) * step
         enddo      
         
         filter = dexp(-( frec - frecu0 )**2/sigma2)
         
! Go through all frequency points of the line
			do f1 = 1, 2*n
         	result(f1) = sum(filter(2*n-f1+1:2*n)*signal(1:f1)) / sum(filter)
         enddo
         
         output(it,:) = result(n-delta:n+delta)
		enddo

   end subroutine convolve_spectrum2


!-----------------------------------------------------------------	
! Convolve spectrum with the instrument's response
!-----------------------------------------------------------------	
	subroutine convolve_spectrum(flux,output)
   real(kind=8) :: flux(nact,nfrq), output(nact,nfrq), frecu, sigma2, t, frecu0, fwhm
   real(kind=8) :: frecu0_filter, normal
	integer :: f1, f2, it
   
   	output = 0.d0
		do it = 1, nr
			frecu0 = dtran(2,it)

! Calculate the FWHM of the instrument in Hz
			fwhm = frecu0**2 / PC * (instrum_sigma(PC * 1.d4 / frecu0) * 1.d-4)

! Calculate the sigma^2 of the telescope at this frequency         
         sigma2 = fwhm**2 / (4.d0 * alog(2.0))
         
! Go through all frequency points of the line
			do f1 = 1, nfrq
				frecu = 1.d0*(f1 - nfrq / 2 + 1) / n_freq_ud_doppler
				frecu = frecu0 + doppler_width(n_radios,it) * frecu         
				frecu0_filter = frecu
               
! For each frequency, add the contribution of all the points of the line weighted by the 
! instrument profile centered at frequency frecu            
				normal = 0.d0
            do f2 = 1, nfrq
					frecu = 1.d0*(f2 - nfrq / 2 + 1) / n_freq_ud_doppler
					frecu = frecu0 + doppler_width(n_radios,it) * frecu					
	            t = (frecu - frecu0_filter)**2 / sigma2
               
            	output(it,f1) = output(it,f1) + flux(it,f2) * dexp(-t)
               normal = normal + dexp(-t)
            enddo
            output(it,f1) = output(it,f1) / normal
			enddo
		enddo

   end subroutine convolve_spectrum


!-----------------------------------------------------------------	
! Calculates telescope's response function for a given impact parameter
!-----------------------------------------------------------------	
	function tel_response(imp)
	real(kind=8) :: tel_response
	real(kind=8), INTENT(IN) :: imp
	real(kind=8) :: coc, norm, d_fwhm

! Telescope's sigma in cm (USE APROPIATELY SIGMA AND FWHM)

		d_fwhm = distancia * dtan(tel_fwhm / 206265.d0) 		
				
		coc = (imp / d_fwhm)  

! Normalization constant to assure energy conservation
		norm = 1.d0/(PI*d_fwhm)!1.d0!2.d0 / (dsqrt(PI) * d_fwhm)
!		tel_response = norm * dexp(-coc**2)	
		tel_response = norm * dexp(-imp**2/d_fwhm)	

		
	end function tel_response

!-----------------------------------------------------------------	
! Calculates emerging spectrum
! INPUT:
!		- which : if =1 then write NLTE spectrum
!					 if =2 then write LTE spectrum
!					 if =3 then write both spectra
! OUTPUT:
!     - Puts in In(:,:) the emerging spectrum in LTE and in NLTE
!-----------------------------------------------------------------	
	subroutine emerge
	integer :: it, f, nim, k
	character*5 :: type_spec
	real(kind=8) :: frecu0, frecu
!	real(kind=8) :: pop_temp(nl*n_radios), emer(nact,nfrq,n_caract), flux(nact,nfrq)
!	real(kind=8) :: convolved(nact,nfrq)
	real(kind=8), allocatable :: pop_temp(:), emer(:,:,:), flux(:,:)
	real(kind=8), allocatable :: convolved(:,:)
	real(kind=8) :: temp, response
	real(kind=8) :: mu_pp(3), wtmu_pp(3), rpar(n_caract,2), erre, sigma
	
!		do it = 1, n_caract
!			write(13,*) p(it), tel_response(p(it))
!		enddo
		
		allocate(pop_temp(nl*n_radios))
		allocate(emer(nact,nfrq,n_caract))
		allocate(flux(nact,nfrq))
		allocate(convolved(nact,nfrq))
		
		select case (converge_flag)			
			case(0)
				type_spec = trim(adjustl(''))
			case(1)
				type_spec = trim(adjustl(''))
			case(2)
				type_spec = trim(adjustl('.LTE'))
		end select
				
		open (UNIT=9,FILE=trim(adjustl(nombre_emerging))//type_spec,STATUS='replace',ACTION='write')
		open (UNIT=8,FILE=trim(adjustl(nombre_flux))//type_spec,STATUS='replace',ACTION='write')
      write(8,*) nact, nfrq

				
		emer = 0.d0
		flux = 0.d0
		pop_temp = 0.d0

		do k = 1, n_caract				
			if (k == 1) then
            	rpar(k,1) = 0.d0
            	rpar(k,2) = 0.5d0*p(k+1)
            else if (k < n_caract) then 
            	rpar(k,1) = rpar(k-1,2)
            	rpar(k,2)= p(k)
            else
            	rpar(k,1) = rpar(k-1,2)
				rpar(k,2) = r(n_radios)
			endif

			if (k == 1) then
            	rpar(k,1) = 0.d0
            	rpar(k,2) = p(k+1)
            else if (k < n_caract) then 
            	rpar(k,1) = rpar(k-1,2)
            	rpar(k,2)= p(k)
            else
            	rpar(k,1) = rpar(k-1,2)
				rpar(k,2) = r(n_radios)
			endif
			print *, rpar(k,1), rpar(k,2)
		enddo

! Put the weights in arcsec
! Now ang_size is the radial angular size
		rpar = rpar / r(n_radios) * ang_size
		
		if (converge_flag == 1 .or. converge_flag == 0) then 
!-----------------------------------------------------------------	
! WRITE NLTE SPECTRUM
!-----------------------------------------------------------------			
		print *, 'Writing NLTE spectrum...'
		flux = 0.d0
!		write(5,*) 'NLTE spectrum'
!		write(6,*) 'NLTE flux'		

		do it = 1, nr 

! Central wavelength of the transition
			frecu0 = dtran(2,it)	

			if (verbose_flag == 1) then
         		if (mod(it,10) == 0) then
					write(*,*) 'Transition : ', it
				endif
			endif
			
			call opacity(it,nim,1,n_radios)	

! Do the formal solution and put the emerging spectrum in In(:,:) 
			if (index(geometry_type,'SPHERICAL') /= 0) then
				if (overlap_flag /= 0) then
					call rtspher_jacobi_overlap_only_int(nim,it)
				else
					call rtspher_jacobi_only_intensity(nim,it)
				endif
												
				
				do k = 1, n_caract
!					call ang_weights(n_radios,k,wtmu)
					response = tel_response(p(k))
					erre = ( rpar(k,1) + rpar(k,2) ) / 2.d0

					sigma = tel_fwhm**2/ (4.d0 * dlog(2.d0))

					emer(it,1:nfrq,k) = In(1:nfrq,k)
								 					 
! These two lines are for implementing a cut in the observing beam
					if (erre < tel_fwhm_cut) then
! Go through the frequency axis
						do f = 1, nfrq

							frecu = 1.d0*(f - nfrq / 2 + 1) / n_freq_ud_doppler
! Obtain frequency axis
							frecu = frecu0 + doppler_width(n_radios,it) * frecu

!					flux(it,f) = flux(it,f) + wtmu * mu(n_radios,k) * In(f,k) 					
!						response = 1.d0   ! Still not working properly
! I think that I have to use In(f,n_caract-k) instead of In(f,k) when doing the variable change from
! mu to p
!						temp = p(k) * In(f,n_caract-k) * response 
!						temp = 2.d0 * PI * p(k) * In(f,k) * response 
!						flux(it,f) = flux(it,f) + temp * impact_weights(k) / (2.d0 * distancia**2)
!						flux(it,f) = flux(it,f) + temp * impact_weights(k) / (r(n_radios-1)**2)
! Pepe version
							flux(it,f) = flux(it,f) + PI*(rpar(k,2)**2-rpar(k,1)**2)*In(f,k)/(PI*sigma)*&
								dexp(-erre**2/sigma)

! Pepe version (try with dividing by sqrt(PI)*sigma^(1/4)
!						flux(it,f) = flux(it,f) + factor*PI*(rpar(k,2)**2-rpar(k,1)**2)*In(f,k)/sqrt(PI*sqrt(sigma))*&
!							dexp(-erre**2/sigma)

						
! One has to be careful if the impact parameter mesh is not fine enough for mu almost equal to 1
! The problem is that there is not outward spectrum for mu=0 and there is a gap between the last
! characteristic and the final radial shell. For the moment, I use r(n_radios-1) instead of
! r(n_radios), but it should be necessary to perform a formal solution for a ray which is very
! close to the mu=0

						enddo !f

					endif
				
				enddo !k

			else
				if (overlap_flag /= 0) then
					call rtpp_jacobi_overlap_only_int(nim,it)
				else
					call rtpp_jacobi_only_intensity(nim,it)
				endif
				
! Set angle quadrature
				select case(angleset_pp)
					case(1)
						mu_pp(1) = 1.d0 / dsqrt(3.d0)
						wtmu_pp(1) = 1.d0
					case(3)
						mu_pp(1) = 0.80847425d0
						mu_pp(2) = 0.57735027d0
						mu_pp(3) = 0.11417547d0

						wtmu_pp(1) = 2.d0 / 6.d0
						wtmu_pp(2) = 2.d0 / 6.d0
						wtmu_pp(3) = 2.d0 / 6.d0
				end select
				
! Assume that we are in the source
				do k = 1, angleset_pp
					response = tel_response(p(k)) 
					emer(it,1:nfrq,k) = In(1:nfrq,k)
! Go through the frequency axis
					do f = 1, nfrq

						frecu = 1.d0*(f - nfrq / 2 + 1) / n_freq_ud_doppler
! Obtain frequency axis
						frecu = frecu0 + doppler_width(n_radios,it) * frecu

!					flux(it,f) = flux(it,f) + wtmu * mu(n_radios,k) * In(f,k) 					
						response = 1.d0   ! Still not working properly
! I think that I have to use In(f,n_caract-k) instead of In(f,k) when doing the variable change from
! mu to p
						temp = mu_pp(k) * In(f,k) * response 
						flux(it,f) = flux(it,f) + 0.5d0 * temp * wtmu_pp(k)
! REVISAR EL CALCULO DEL FLUJO
	
					enddo !f
				
				enddo !k
				
			endif
			

		enddo	! it
      
! Convolve spectrum with the instrument's spectral response
!      call convolve_spectrum(flux,convolved)
		convolved = flux
		
! Now write spectrum to file
		do k = 1, n_caract
!			write(5,*) 'Characteristic : ', k			
			do it = 1, nr
				frecu0 = dtran(2,it)	
				do f = 1, nfrq
					frecu = 1.d0*(f - nfrq / 2 + 1) / n_freq_ud_doppler
					frecu = frecu0 + doppler_width(n_radios,it) * frecu
					write(9,FMT='(E23.17,5X,E23.17)') frecu, emer(it,f,k)
					
! Flux has an integration over mu, then it does not depend on mu. We write the flux only once
					if (k == 1) then
!						write(6,FMT='(E23.17,5X,E23.17)') frecu, flux(it,f) !/ r(n_radios)**2
						write(8,FMT='(E23.17,5X,E23.17)') frecu, convolved(it,f) !/ r(n_radios)**2
						write(7,FMT='(E23.17,5X,E23.17)') frecu, flux(it,f) !/ r(n_radios)**2                  
					endif
					
				enddo
			enddo
		enddo

		endif
		
		if (converge_flag == 1 .or. converge_flag == 2 .or. converge_flag == 0) then
!-----------------------------------------------------------------	
! WRITE LTE SPECTRUM
!-----------------------------------------------------------------					
!		write(5,*) 'LTE spectrum'
!		write(6,*) 'LTE flux'
		print *, 'Writing LTE spectrum...'		
		flux = 0.d0
! Exchange LTE populations with NLTE populations and make a backup of NLTE to restart it at the end		
		pop_temp = pop
		pop = popl
		
		do it = 1, nr 

! Central wavelength of the transition
			frecu0 = dtran(2,it)	
			
			if (verbose_flag == 1) then
         		if (mod(it,10) == 0) then
					write(*,*) 'Transition : ', it
				endif
			endif

						
			call opacity(it,nim,1,n_radios)	

! Do the formal solution and put in In(:,:) the emerging spectrum
			if (index(geometry_type,'SPHERICAL') /= 0) then
				if (overlap_flag /= 0) then
					call rtspher_jacobi_overlap_only_int(nim,it)
				else
					call rtspher_jacobi_only_intensity(nim,it)
				endif
				
				do k = 1, n_caract

					response = tel_response(p(k))
					erre = ( rpar(k,1) + rpar(k,2) ) / 2.d0
					sigma = tel_fwhm**2/ (4.d0 * dlog(2.d0))

					emer(it,1:nfrq,k) = In(1:nfrq,k)
!					call ang_weights(n_radios,k,wtmu)
!					response = tel_response(p(k)) 				
! These two lines are for implementing a cut in the observing beam
					if (erre < tel_fwhm_cut) then

! Go through the frequency axis
						do f = 1, nfrq

							frecu = 1.d0*(f - nfrq / 2 + 1) / n_freq_ud_doppler
! Obtain frequency axis
							frecu = frecu0 + doppler_width(n_radios,it) * frecu

							flux(it,f) = flux(it,f) + PI*(rpar(k,2)**2-rpar(k,1)**2)*In(f,k) / (sqrt(PI)*sigma)*&
								dexp(-erre**2/sigma**2)
					
						enddo !f

					endif
				
				enddo !k

				
			else
				if (overlap_flag /= 0) then
					call rtpp_jacobi_overlap_only_int(nim,it)
				else
					call rtpp_jacobi_only_intensity(nim,it)
				endif
				
				do k = 1, angleset_pp
					response = tel_response(p(k)) 
					emer(it,1:nfrq,k) = In(1:nfrq,k)
! Go through the frequency axis
					do f = 1, nfrq

						frecu = 1.d0*(f - nfrq / 2 + 1) / n_freq_ud_doppler
! Obtain frequency axis
						frecu = frecu0 + doppler_width(n_radios,it) * frecu

!					flux(it,f) = flux(it,f) + wtmu * mu(n_radios,k) * In(f,k) 					
						response = 1.d0   ! Still not working properly
! I think that I have to use In(f,n_caract-k) instead of In(f,k) when doing the variable change from
! mu to p
						temp = mu_pp(k) * In(f,k) * response 
						flux(it,f) = flux(it,f) + 0.5d0 * temp * wtmu_pp(k)
! REVISAR EL CALCULO DEL FLUJO
	
					enddo !f
				
				enddo !k
				
			endif


		enddo	! it
		
! Now write spectrum to file		
		do k = 1, n_caract
			write(*,*) 'Characteristic : ', k, nr, nfrq
			do it = 1, nr
				frecu0 = dtran(2,it)	
				do f = 1, nfrq
					frecu = 1.d0*(f - nfrq / 2 + 1) / n_freq_ud_doppler
					frecu = frecu0 + doppler_width(n_radios,it) * frecu
					write(9,FMT='(E23.17,5X,E23.17)') frecu, emer(it,f,k)
					
! Flux has an integration over mu, then it does not depend on mu. We write the flux only once
					if (k == 1) then
						write(8,FMT='(E23.17,5X,E23.17)') frecu, flux(it,f) !/ r(n_radios)**2
					endif					
					
				enddo
			enddo
		enddo		

		endif
				
! Restart NLTE populations at the end
		pop = pop_temp	
		
		close(6)
		close(5)

		deallocate(pop_temp)
		deallocate(emer)
		deallocate(flux)
		deallocate(convolved)
	
	end subroutine emerge

end module emerging
