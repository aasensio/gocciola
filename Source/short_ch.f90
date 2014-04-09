module short_ch
use variables
implicit none
contains

! *********************************************************
! *********************************************************
! SHORT-CHARACTERISTICS FORMAL SOLUTION ROUTINES
! *********************************************************
! *********************************************************

! ---------------------------------------------------------
! Solves RTE using short-characteristics. Source function is supposed to vary parabolically
! between points M, O and P and RTE is integrated analitically, resulting in
!            I(O)=I(M)*exp(-tau(MO))+psi(M)*S(M)+psi(O)*S(O)+psi(P)*S(P)
! where psi are functions of the optical distance tau(MO) and tau(OP)
! It returns the value of the three psi coefficients
! ---------------------------------------------------------
   subroutine par_sc(dtm, dtp, psim, psi0, psip)
   real*8, INTENT(IN) :: dtm, dtp
   real*8, INTENT(INOUT) :: psim, psi0, psip
   real*8 :: exu, u0, u1 ,u2, d2, d3, d4
         if (dtm.ge.trick) then
            exu=dexp(-dtm)
            u0=1.d0-exu
            u1=dtm-1.d0+exu
            u2=dtm**2-2.d0*dtm+2.d0-2.d0*exu
         else
            d2=dtm**2
            d3=dtm**3
            d4=dtm**4
            u0=dtm-(d2/2.d0)
            u1=(d2/2.d0)-(d3/6.d0)
            u2=(d3/3.d0)-(d4/12.d0)
        endif

        if (dtm * dtp /= 0.d0 .and. dtm**2 /= 0.d0 .and. dtp**2 /= 0.d0) then			  
			  psim=(u2-u1*(dtp+2.d0*dtm))/(dtm*(dtm+dtp))+u0
      	  psi0=(u1*(dtm+dtp)-u2)/(dtm*dtp)
      	  psip=(u2-dtm*u1)/(dtp*(dtm+dtp))
	 	  else
		  	  psim = 0.d0
			  psi0 = 0.d0
			  psip = 0.d0
		  endif
		  
   end subroutine par_sc

! ---------------------------------------------------------
! Solves RTE using short-characteristics. Source function is supposed to vary linearly
! between points M and O and RTE is integrated analitically, resulting in
!            I(O)=I(M)*exp(-tau(MO))+psi(M)*S(M)+psi(O)*S(O)
! where psi are functions of the optical distance tau(MO)
! It returns the value of the two psi coefficients
! ---------------------------------------------------------
   subroutine lin_sc(dtm, psim, psi0)
   real*8, INTENT(IN) :: dtm
	real*8, INTENT(INOUT) :: psim, psi0
   real*8 :: exu, u0, u1, c0, cm, d2

      if (dtm.ge.trick) then
         exu=dexp(-dtm)
         u0=1.d0-exu
         u1=dtm-1.d0+exu

         c0=u1/dtm
         cm=u0-c0
      else   
         d2=dtm**2.d0
         c0=(dtm/2.d0)-(d2/6.d0)
         cm=(dtm/2.d0)-(d2/3.d0)
      endif
		psi0 = c0
		psim = cm

   end subroutine lin_sc
	
end module short_ch
