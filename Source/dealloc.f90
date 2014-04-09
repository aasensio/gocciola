module dealloc
use variables
implicit none
contains

! *********************************************************
! *********************************************************
! FINAL CLEANING OF MEMORY
! *********************************************************
! *********************************************************

! ----------------------------------------------------------
! Deallocates all allocated variables
! ----------------------------------------------------------

	subroutine clean
		return	
		print *, 'Cleaning allocated memory...'
		if (allocated(r)) deallocate(r)
		if (allocated(p)) deallocate(p)
		if (allocated(datph)) deallocate(datph)
		if (allocated(nh)) deallocate(nh)
		if (allocated(dion)) deallocate(dion)
		if (allocated(nli)) deallocate(nli)
		if (allocated(dlevel)) deallocate(dlevel)
		if (allocated(itran)) deallocate(itran)
		if (allocated(dtran)) deallocate(dtran)
		if (allocated(doppler_width)) deallocate(doppler_width)
		if (allocated(collision)) deallocate(collision)
		if (allocated(partition)) deallocate(partition)
		if (allocated(mu)) deallocate(mu)
		if (allocated(B)) deallocate(B)
		if (allocated(B_dust)) deallocate(B_dust)
		if (allocated(chil)) deallocate(chil)
		if (allocated(kappa)) deallocate(kappa)
		if (allocated(kappa_dust)) deallocate(kappa_dust)
		if (allocated(overl)) deallocate(overl)
		if (allocated(Sl)) deallocate(Sl)
		if (allocated(Lstar)) deallocate(Lstar)
		if (allocated(Jbar)) deallocate(Jbar)
		if (allocated(In)) deallocate(In)
		if (allocated(wtnu)) deallocate(wtnu)
		if (allocated(perfil)) deallocate(perfil)
		if (allocated(dentro)) deallocate(dentro)
		if (allocated(fuera)) deallocate(fuera)
		if (allocated(Jbar_total)) deallocate(Jbar_total)
		if (allocated(Lstar_total)) deallocate(Lstar_total)
		if (allocated(scratch)) deallocate(scratch)
! Esta falla en GS
		if (allocated(pop_final)) deallocate(pop_final)
		if (allocated(tau)) deallocate(tau)
		if (allocated(chic)) deallocate(chic)
		if (allocated(Sc)) deallocate(Sc)
		if (allocated(eps_silicate)) deallocate(eps_silicate)
		if (allocated(eps_graphite_par)) deallocate(eps_graphite_par)
		if (allocated(eps_graphite_per)) deallocate(eps_graphite_per)		
		if (allocated(graphite_lambda)) deallocate(graphite_lambda)
		if (allocated(silicate_lambda)) deallocate(silicate_lambda)
		if (allocated(label)) deallocate(label)
		if (allocated(popl)) deallocate(popl)
		if (allocated(popi)) deallocate(popi)		
		if (allocated(pop)) deallocate(pop)				
		if (allocated(wtmu)) deallocate(wtmu)
		if (allocated(mu_pp)) deallocate(mu_pp)

	
	end subroutine clean

end module dealloc
