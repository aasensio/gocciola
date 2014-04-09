!********************************************************
! Convolution module for applying the observing resolution
!********************************************************
module convol
use variables
use emerging
implicit none
contains


! ---------------------------------------------------------
! Returns the mean of a vector
! ---------------------------------------------------------
	function mean2(v)
	real(kind=8) :: mean2, v(:)
	integer :: n
		n = size(v)
		mean2 = sum(v) / n
	end function mean2





!-----------------------------------------------------------------
! This subroutine performs the FFT of a complex signal of size nn, given in 
! the pseudo-complex form, taking 2*nn arrays
!-----------------------------------------------------------------
SUBROUTINE four1(data,nn,isign)
INTEGER isign,nn
REAL(kind=8) data(2*nn)
INTEGER i,istep,j,m,mmax,n
REAL(kind=8) tempi,tempr
DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      n=2*nn
      j=1
      do 11 i=1,n,2
        if(j.gt.i)then
          tempr=data(j)
          tempi=data(j+1)
          data(j)=data(i)
          data(j+1)=data(i+1)
          data(i)=tempr
          data(i+1)=tempi
        endif
        m=n/2
1       if ((m.ge.2).and.(j.gt.m)) then
          j=j-m
          m=m/2
        goto 1
        endif
        j=j+m
11    continue
      mmax=2
2     if (n.gt.mmax) then
        istep=2*mmax
        theta=6.28318530717959d0/(isign*mmax)
        wpr=-2.d0*sin(0.5d0*theta)**2
        wpi=sin(theta)
        wr=1.d0
        wi=0.d0
        do 13 m=1,mmax,2
          do 12 i=m,n,istep
            j=i+mmax
            tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
            tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
            data(j)=data(i)-tempr
            data(j+1)=data(i+1)-tempi
            data(i)=data(i)+tempr
            data(i+1)=data(i+1)+tempi
12        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
13      continue
        mmax=istep
      goto 2
      endif
      return
      END subroutine four1
		
!-----------------------------------------------------------------
! This subroutine convolves the signal with the filter. The size of the filter has to
! be lower than the size of the signal
! The filter has to be given in wrap-around order
!-----------------------------------------------------------------
subroutine convolve(signal,filter)
real(kind=8), INTENT(INOUT) :: signal(:)
real(kind=8), INTENT(IN) :: filter(:)
real(kind=8), allocatable :: signal_complex(:), filter_complex(:), ans(:)
real(kind=8) :: filter_completo(size(signal)), a, b, c, d
integer :: m, n, n2, i

	n = size(signal)
	m = size(filter)
	n2 = 2*n
	
	allocate(signal_complex(n2))
	allocate(filter_complex(n2))
	allocate(ans(n2))
	
! Zero padding of the filter
	do i = 1, (m-1)/2
		filter_completo(i) = filter(i)
		filter_completo(n+1-i) = filter(m+1-i)
	enddo
	do i = (m+3)/2,n-(m-1)/2
		filter_completo(i) = 0.0
	enddo

! Put the filter in pseudo-complex format
	do i = 1, n
		filter_complex(2*(i-1)+1) = filter_completo(i)
	enddo
	
! Put the signal in pseudo-complex format
	do i = 1, n
		signal_complex(2*(i-1)+1) = signal(i)
	enddo	

! FFT both the filter and the signal
	call four1(filter_complex,n,1)
	call four1(signal_complex,n,1)
	
! Do the complex multiplication
	do i = 1, n
		a = signal_complex(2*(i-1)+1)
		b = signal_complex(2*i)
		c = filter_complex(2*(i-1)+1)
		d = filter_complex(2*i)
		ans(2*(i-1)+1) = a*c - b*d
		ans(2*i) = a*d + b*c
	enddo

! Inverse FFT of the product	
	call four1(ans,n,-1)
	
! In the inverse FFT, the result is multiplied by n
	ans = ans / n
	
! Get the real part
	do i = 1, n
		signal(i) = ans(2*(i-1)+1)
	enddo
		
	deallocate(signal_complex)
	deallocate(filter_complex)
	deallocate(ans)
	
end subroutine convolve

!-------------------------------------------------------------------
! Return the Stokes profiles convolved with the given filter 
!-------------------------------------------------------------------
subroutine convolve_macroturbulence(dw,emerg,emerg_conv)
real(kind=8) :: dw(:), emerg(:), emerg_conv(:)
real(kind=8), allocatable :: filter(:), ans(:), padding(:)
real(kind=8) :: macro_hz, x, wave, step, macro_A
integer :: i, m, n, n2

		n = size(dw)
      wave = mean2(dw)
		macro_hz = instrum_sigma(3.d14 / wave)
      step = dw(2) - dw(1)
		m = 3 * macro_hz / step
		m = 2 * m + 1
		allocate(filter(m))

! Build the filter. It is a gaussian of sigma=macro_A, but it has to be built in
! the correct way in order to have a working convolution.
		macro_A = 1.d0  ! DANGER!!!!!!!!!!!!
		do i = 1, m/2+1
			x = (i-1.d0)*step / macro_A
			filter(i) = exp(-x**2)
		enddo
		
		do i = m/2+2, m
			filter(i) = filter(m-i+1)
		enddo

! Normalize the filter
		filter = filter / sum(filter)
		
! Find the nearest higher power of 2		
		i = 1
     	do while (2.d0**i < n)
			i = i + 1
		enddo     
     	n2=2.d0**i
		
		allocate(padding(n2))
		allocate(ans(2*n2))
				
		padding(1:n) = emerg(:)
		padding(n:n2) = padding(n)
		call convolve(padding,filter)
		emerg_conv(:) = padding(1:n)
		
		if (allocated(filter)) deallocate(filter)
		if (allocated(ans)) deallocate(ans)

end subroutine convolve_macroturbulence

end module convol
