module maths
use variables
implicit none
contains
! *********************************************************
! *********************************************************
! MATHEMATICAL ROUTINES
! *********************************************************
! *********************************************************

!-----------------------------------------------------------------
! Returns the value of the exponential integral En(x)
!-----------------------------------------------------------------
      FUNCTION expint(n,x)
      INTEGER n,MAXIT
      REAL(kind=8) expint,x,EPS,FPMIN,EULER
      PARAMETER (MAXIT=100,EPS=1.e-7,FPMIN=1.e-30,EULER=.5772156649)
      INTEGER i,ii,nm1
      REAL(kind=8) a,b,c,d,del,fact,h,psi
      nm1=n-1
      if(n.lt.0.or.x.lt.0..or.(x.eq.0..and.(n.eq.0.or.n.eq.1)))then
        print *,'bad arguments in expint'
		  stop
      else if(n.eq.0)then
        expint=exp(-x)/x
      else if(x.eq.0.)then
        expint=1./nm1
      else if(x.gt.1.)then
        b=x+n
        c=1./FPMIN
        d=1./b
        h=d
        do 11 i=1,MAXIT
          a=-i*(nm1+i)
          b=b+2.
          d=1./(a*d+b)
          c=b+a/c
          del=c*d
          h=h*del
          if(abs(del-1.).lt.EPS)then
            expint=h*exp(-x)
            return
          endif
11      continue
        print *,'continued fraction failed in expint'
		  stop
      else
        if(nm1.ne.0)then
          expint=1./nm1
        else
          expint=-log(x)-EULER
        endif
        fact=1.
        do 13 i=1,MAXIT
          fact=-fact*x/i
          if(i.ne.nm1)then
            del=-fact/(i-nm1)
          else
            psi=-EULER
            do 12 ii=1,nm1
              psi=psi+1./ii
12          continue
            del=fact*(-log(x)+psi)
          endif
          expint=expint+del
          if(abs(del).lt.abs(expint)*EPS) return
13      continue
        print *, 'series failed in expint'
		  stop
      endif
      return
      END function expint
		
!-----------------------------------------------------------------
! Returns the weights (w) and the abscissas (x) for a Gaussian integration using the 
! Gauss-Legendre formula, using n points
!-----------------------------------------------------------------
	subroutine gauleg(x1,x2,x,w,n)
	integer, INTENT(IN) :: n
	real(kind=8), INTENT(IN) :: x1,x2
	real(kind=8), INTENT(INOUT) :: x(n),w(n)
	real(kind=8), parameter :: eps = 3.d-14
	integer :: i,j,m
	real(kind=8) :: p1,p2,p3,pp,xl,xm,z,z1
      
	m=(n+1)/2
   xm=0.5d0*(x2+x1)
   xl=0.5d0*(x2-x1)
   do i=1,m
   	z=dcos(3.141592654d0*(i-.25d0)/(n+.5d0))
1  	continue
   	p1=1.d0
   	p2=0.d0
   	do j=1,n
   		p3=p2
      	p2=p1
      	p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
		enddo
   	pp=n*(z*p1-p2)/(z*z-1.d0)
   	z1=z
   	z=z1-p1/pp
  		if(abs(z-z1).gt.EPS)goto 1
   	x(i)=xm-xl*z
   	x(n+1-i)=xm+xl*z
   	w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
   	w(n+1-i)=w(i)
	enddo
	
	end subroutine gauleg

! ---------------------------------------------------------
! Returns the mean of the vector x
! ---------------------------------------------------------
	function mean(x)
	real(kind=8) :: mean
	real(kind=8), INTENT(IN) :: x(:)
		mean = sum(x) / dble(size(x))
	end function mean
	
! ---------------------------------------------------------
! Derivative of a function using a 3-point Lagrange interpolation formulae
! ---------------------------------------------------------
	function deriv(x0, x1, x2, f0, f1, f2, which)
	real(kind=8) :: deriv
	real(kind=8), INTENT(IN) :: x0, x1, x2, f0, f1, f2
	integer, INTENT(IN) :: which
	real(kind=8) :: h1, h2, c0, c1, c2, x
	
		x = 0.d0
		if (which == 0) x = x0
		if (which == 1) x = x1
		if (which == 2) x = x2
		h1 = x1 - x0
		h2 = x2 - x1
		c0 = 1.d0 / (h1*(h1+h2))
		c1 = 1.d0 / (h1*h2)
		c2 = 1.d0 / (h2*(h1+h2))
		
		deriv = f0 * c0 * ( (x-x2) + (x-x1) ) - f1 * c1 * ( (x-x0) + (x-x2) ) + &
			f2 * c2 * ( (x-x0) + (x-x1) ) 
			
	end function deriv

! ---------------------------------------------------------
! Derivative of a function using a 3-point Lagrange interpolation formulae
! It seems to be wrong
! ---------------------------------------------------------
	function lag_deriv(x0, x1, x2, f0, f1, f2, which)
	real(kind=8) :: lag_deriv
	real(kind=8), INTENT(IN) :: x0, x1, x2, f0, f1, f2
	integer, INTENT(IN) :: which
	real(kind=8) :: h1, h2, c0, c1, c2
		h1 = x1 - x0
		h2 = x2 - x1
		c0 = 1.d0 / (h1*(h1+h2))
		c1 = 1.d0 / (h1*h2)
		c2 = 1.d0 / (h2*(h1+h2))

		lag_deriv = 0.d0
		
		if (which == 0) then
			lag_deriv = -c0 * f0 * (2d0 * h1 + h2) - c1 * f1 * (h1 + h2) - c2 * f2 * h1
		endif
		
		if (which == 1) then
			lag_deriv = -c0 * f0 * h2 + c1 * f1 * (h1 - h2) + c2 * f2 * h1
		endif
		
		if (which == 2) then			
			lag_deriv = c0 * f0 * h2 + c1 * f1 * (h1 + h2) + c2 * f2 * (2.d0 * h2 + h1)
		endif
		
	end function lag_deriv
! ---------------------------------------------------------
! LU decomposition of a matrix
!  INPUT:
!		- a is the matrix to decompose
!		
!  OUTPUT:
!		- a is the LU decomposition of a
!		- indx is a vector that records the row permutation effected by the partial pivoting
!		- d takes values +1/-1 depending on whether the number of row interchanges was odd or even
! ---------------------------------------------------------
	subroutine ludcmp(a,indx,d)
	integer, INTENT(INOUT) :: indx(:)
	real(kind=8), INTENT(INOUT) :: a(:,:), d
	real(kind=8), parameter :: TINY = 1.d-20
	integer :: i, imax, j, k, n
	real(kind=8) :: aamax, dum, sum, vv(size(a,1))
		d = 1.d0
		n = size(a,1)

		imax = 0
		
		do i = 1, n
			aamax = 0.d0	
			aamax = maxval(dabs(a(i,:)))
			if (aamax == 0.d0) print *, 'Singular matrix in LU decomposition'
			vv(i) = 1.d0 / aamax
		enddo
		
		do j = 1, n
			do i = 1, j-1
				sum = a(i,j)
				do k = 1, i-1
					sum = sum - a(i,k) * a(k,j)
				enddo
				a(i,j) = sum
			enddo
			aamax = 0.d0
			do i = j, n
				sum = a(i,j)
				do k = 1, j-1
					sum = sum - a(i,k) * a(k,j)
				enddo
				a(i,j) = sum
				dum = vv(i) * dabs(sum)
				if (dum >= aamax) then
					imax = i
					aamax = dum
				endif				
			enddo
			if (j /= imax) then
				do k = 1, n
					dum = a(imax,k)
					a(imax,k) = a(j,k)
					a(j,k) = dum
				enddo
				d = -d
				vv(imax) = vv(j)
			endif
			indx(j) = imax
			if (a(j,j) == 0.d0) a(j,j) = TINY
			if (j /= n) then
				dum = 1.d0 / a(j,j)
				do i = j+1, n
					a(i,j) = a(i,j) * dum
				enddo
			endif
		enddo
	
	end subroutine ludcmp

! ---------------------------------------------------------
! Solves the set of equations AX=b where A is the LU decomposition of a matrix
!  INPUT:
!		- a is the LU decomposition of the system matrix
!		- b is the right hand side vector of the system
! 		- indx is the vector returned by ludcmp
!  OUTPUT:
!		- b is the solution of the system
! ---------------------------------------------------------
	subroutine lubksb(a,indx,b)
	real(kind=8), INTENT(IN) :: a(:,:)
	real(kind=8), INTENT(INOUT) :: b(:)
	integer, INTENT(IN) :: indx(:)
	integer :: i, ii, n, j, ll
	real(kind=8) :: sum
		n = size(a,1)
		ii = 0
		do i = 1, n
			ll = indx(i)
			sum = b(ll)
			b(ll) = b(i)
			if (ii /= 0) then
				do j = ii, i-1
					sum = sum - a(i,j) * b(j)
				enddo
			else if (sum /= 0.d0) then
				ii = i
			endif
			b(i) = sum
		enddo
		do i = n, 1, -1
			sum = b(i)
			do j = i+1, n
				sum = sum - a(i,j) * b(j)
			enddo
			b(i) = sum / a(i,i)
		enddo
	end subroutine lubksb

! ---------------------------------------------------------
! This subroutine solves a linear system of equations using a BiCGStab iterative method
! ---------------------------------------------------------		  
	subroutine bicgstab(a,b)
	real(kind=8), INTENT(INOUT) :: a(:,:), b(:)
	real(kind=8) :: x(size(b)), r(size(b)), rhat(size(b)), p(size(b)), phat(size(b)), v(size(b))
	real(kind=8) :: s(size(b)), shat(size(b)), t(size(b)), delta(size(b))
	real(kind=8) :: rho, rho_ant, alpha, beta, omega, relative
	integer :: i
		
		alpha = 0.d0
		omega = 0.d0
! Initial solution
		x = 1.d0
	
		r = b - matmul(A,x)
		rhat = r
		rho = 1.d0
		relative = 1.d10
		i = 1
		do while (relative > 1.d-10)
			rho_ant = rho
			rho = sum(rhat*r)
			if (rho == 0) then 
				stop
			endif
			if (i == 1) then
				p = r
			else
				beta = (rho/rho_ant) * (alpha/omega)
				p = r + beta * (p - omega * v)
			endif
			phat = p
			v = matmul(A,phat)
			alpha = rho / sum(rhat*v)
			s = r - alpha * v

			shat = s
			t = matmul(A,shat)
			omega = sum(t*s)/sum(t*t)
			delta = alpha * p + omega * s
			relative = maxval(abs(delta) / x)
			x = x + delta		
			r = s - omega * t
			i = i + 1
		enddo	
		
		b = x
		  
	end subroutine bicgstab	
	
! ---------------------------------------------------------
! This subroutine solves a linear system of equations using the SLAP routines
! It uses an Incomplete LU BiConjugate Gradient Sparse Ax=b solver
! ---------------------------------------------------------		  
	subroutine slapsolver(a,b,initial)
	real(kind=8), INTENT(INOUT) :: a(:,:), b(:), initial(:)
	integer :: i, j, k
	real(kind=8), allocatable :: a_sparse(:), x(:), rwork(:)
	integer, allocatable :: ia(:), ja(:), iwork(:)
	integer :: n, nelt, isym, itol, itmax, iter, ierr, iunit, lenw, leniw, nsave
	real(kind=8) :: tol
	integer, allocatable :: indx(:)

		iter = 0
		ierr = 1
		n = size(b)
		allocate(x(n))
		k = 0
		do i = 1, n
			do j = 1, n
				if (a(j,i) /= 0.d0) k = k + 1
			enddo
		enddo
		nelt = k
		allocate(a_sparse(nelt))
		allocate(ia(nelt))
		allocate(ja(nelt))		
		
		k = 1
		do i = 1, n
			do j = 1, n
				if (a(j,i) /= 0.d0) then
					a_sparse(k) = a(j,i)
					ia(k) = j
					ja(k) = i
					k = k + 1
				endif
			enddo
		enddo
		
! First try with BiCGStab
		x = initial
		isym = 0
		itol = 2
		tol = 1.d-10
		itmax = 10000
		iunit = 0
		lenw = n*n
		leniw = n*n
		nsave = 10
		allocate(iwork(leniw))
		allocate(rwork(lenw))
		
!		print *, 'Solving the system...'
! 		call ds2y( n, nelt, ia, ja, a_sparse, isym )
! 		call dslucs(n, b, x, nelt, ia, ja, a_sparse, isym, itol, tol, &
! 			itmax, iter, err, ierr, iunit, rwork, lenw, iwork, leniw )
!		print *, 'Number of iterations : ', iter, 'of ', itmax
!		print *, 'Error : ', err, ierr
		if (iter >= itmax+1 .or. ierr /= 0) then
			print *, 'Couldn''t reach convergence. Trying other initialization...'
!			write(17,*) 'Couldn''t reach convergence. Trying other initialization...'
			x = 1.d0
			isym = 0
			itol = 2
			tol = 1.d-10
			itmax = 10000
			iunit = 0
			lenw = n*n
			leniw = n*n
			nsave = 10
! 			call dslucs(n, b, x, nelt, ia, ja, a_sparse, isym, itol, tol, &
! 			itmax, iter, err, ierr, iunit, rwork, lenw, iwork, leniw )
!			print *, 'Number of iterations : ', iter, 'of ', itmax
!			print *, 'Error : ', err, ierr
		endif
		
! Let's try GMRES		
		if (iter >= itmax+1 .or. ierr /= 0) then
			print *, 'Couldn''t reach convergence. Trying GMRES...'
!			write(17,*) 'Couldn''t reach convergence. Trying GMRES...'
			x = 0.d0
			isym = 0
			itol = 2
			tol = 1.d-10
			itmax = 10000
			iunit = 0
			lenw = n*n
			leniw = n*n
			nsave = 10
! 			call dsdgmr(n, b, x, nelt, ia, ja, a_sparse, isym, nsave, itol, tol, &
! 			itmax, iter, err, ierr, iunit, rwork, lenw, iwork, leniw )
!			print *, 'Number of iterations : ', iter, 'of ', itmax
!			print *, 'Error : ', err, ierr
		endif
		
		if (iter >= itmax+1 .or. ierr /= 0) then
			print *, 'Couldn''t reach convergence. Trying other initialization...'
!			write(17,*) 'Couldn''t reach convergence. Trying other initialization...'
			x = 1.d0
			isym = 0
			itol = 0
			tol = 1.d-10
			itmax = 10000
			iunit = 0
			lenw = n*n
			leniw = n*n
			nsave = 10
! 			call dsdgmr(n, b, x, nelt, ia, ja, a_sparse, isym, nsave, itol, tol, &
! 			itmax, iter, err, ierr, iunit, rwork, lenw, iwork, leniw )
!			print *, 'Number of iterations : ', iter, 'of ', itmax
!			print *, 'Error : ', err, ierr
		endif
		
		b = x
		
		deallocate(x)
		deallocate(a_sparse)
		deallocate(ia)
		deallocate(ja)
		deallocate(iwork)
		deallocate(rwork)
		
		if (iter >= itmax+1 .or. ierr /= 0) then
			print *, 'Couldn''t reach convergence. Trying LU decomposition...'
			write(17,*) 'Couldn''t reach convergence. Trying LU decomposition...'
			allocate(indx(n))
			call ludcmp(a,indx,tol)
			call lubksb(a,indx,b)
			deallocate(indx)
		endif
		
	end subroutine slapsolver		
	
	
! ---------------------------------------------------------
! Interface to the LU solver
!  INPUT:
!		- a is the system matrix (if alu=0) or the LU decomposition (if alu /= 0)
!		- b is the right hand side of the system
!		- indx is the vector of permutations
!		- alu is a flag to indicate if LU decomposition has to be done or not
! ---------------------------------------------------------	
	subroutine linsolve(a, b, indx, alu)
	real(kind=8), INTENT(INOUT) :: a(:,:), b(:)
	integer, INTENT(INOUT) :: indx(:)
	real(kind=8) :: d, v(size(b),size(b)), w(size(b)), solution(size(b))
	integer :: alu, n
		n = size(b)
		
		if (linear_solver_algorithm == 0) then
			if (alu == 0) then
				call ludcmp(a,indx,d)
				call lubksb(a,indx,b)
			else 
				call lubksb(a,indx,b)
			endif
		else if (linear_solver_algorithm == 1) then
			call svdcmp(a,n,n,n,n,w,v)
			call svbksb(a,w,v,n,n,n,n,b,solution)
			b = solution
		else if (linear_solver_algorithm == 2) then
			print *, 'Not implemented yet'
! 			call slapsolver(a,b,initial)
		else
! 			call bicgstab(a,b)
		endif
		
	end subroutine linsolve

! ---------------------------------------------------------
! Returns the SVD decomposition of a matrix
! ---------------------------------------------------------		
   subroutine svdcmp(a,m,n,mp,np,w,v)
   integer :: m,mp,n,np
   integer, parameter :: nmax=500
   real(kind=8) :: a(mp,np),v(np,np),w(np)
   integer :: i,its,j,jj,k,l,nm
   real(kind=8) :: anorm,c,f,g,h,s,scale,x,y,z,rv1(NMAX)
   	nm = 0
   	g=0.d0
   	scale=0.d0
   	anorm=0.d0
   	do 25 i=1,n
   	  l=i+1
   	  rv1(i)=scale*g
   	  g=0.d0
   	  s=0.d0
   	  scale=0.d0
   	  if(i.le.m)then
      	 do 11 k=i,m
	      	scale=scale+abs(a(k,i))
11           continue
	   	 if(scale.ne.0.d0)then
	      	do 12 k=i,m
	      	  a(k,i)=a(k,i)/scale
	      	  s=s+a(k,i)*a(k,i)
12             continue
	      	f=a(i,i)
	      	g=-sign(sqrt(s),f)
	      	h=f*g-s
	      	a(i,i)=f-g
	      	do 15 j=l,n
	      	  s=0.d0
	      	  do 13 k=i,m
	         	 s=s+a(k,i)*a(k,j)
13               continue
	      	  f=s/h
	      	  do 14 k=i,m
	         	 a(k,j)=a(k,j)+f*a(k,i)
14               continue
15             continue
	      	do 16 k=i,m
	      	  a(k,i)=scale*a(k,i)
16             continue
	   	 endif
	     endif
	     w(i)=scale *g
	     g=0.d0
	     s=0.d0
	     scale=0.d0
	     if((i.le.m).and.(i.ne.n))then
	   	 do 17 k=l,n
	      	scale=scale+abs(a(i,k))
17           continue
	   	 if(scale.ne.0.d0)then
	      	do 18 k=l,n
	      	  a(i,k)=a(i,k)/scale
	      	  s=s+a(i,k)*a(i,k)
18             continue
	      	f=a(i,l)
	      	g=-sign(sqrt(s),f)
	      	h=f*g-s
	      	a(i,l)=f-g
	      	do 19 k=l,n
	      	  rv1(k)=a(i,k)/h
19             continue
	      	do 23 j=l,m
	      	  s=0.d0
	      	  do 21 k=l,n
	         	 s=s+a(j,k)*a(i,k)
21               continue
	      	  do 22 k=l,n
	         	 a(j,k)=a(j,k)+s*rv1(k)
22               continue
23             continue
	      	do 24 k=l,n
	      	  a(i,k)=scale*a(i,k)
24             continue
	   	 endif
	     endif
	     anorm=max(anorm,(abs(w(i))+abs(rv1(i))))
25       continue
	   do 32 i=n,1,-1
	     if(i.lt.n)then
	   	 if(g.ne.0.d0)then
	      	do 26 j=l,n
	      	  v(j,i)=(a(i,j)/a(i,l))/g
26             continue
	      	do 29 j=l,n
	      	  s=0.d0
	      	  do 27 k=l,n
	         	 s=s+a(i,k)*v(k,j)
27               continue
	      	  do 28 k=l,n
	         	 v(k,j)=v(k,j)+s*v(k,i)
28               continue
29             continue
	   	 endif
	   	 do 31 j=l,n
	      	v(i,j)=0.d0
	      	v(j,i)=0.d0
31           continue
	     endif
	     v(i,i)=1.d0
	     g=rv1(i)
	     l=i
32       continue
	   do 39 i=min(m,n),1,-1
	     l=i+1
	     g=w(i)
	     do 33 j=l,n
	   	 a(i,j)=0.d0
33         continue
	     if(g.ne.0.d0)then
	   	 g=1.d0/g
	   	 do 36 j=l,n
	      	s=0.d0
	      	do 34 k=l,m
	      	  s=s+a(k,i)*a(k,j)
34             continue
	      	f=(s/a(i,i))*g
	      	do 35 k=i,m
	      	  a(k,j)=a(k,j)+f*a(k,i)
35             continue
36           continue
	   	 do 37 j=i,m
	      	a(j,i)=a(j,i)*g
37           continue
	     else
	   	 do 38 j= i,m
	      	a(j,i)=0.d0
38           continue
	     endif
	     a(i,i)=a(i,i)+1.d0
39       continue
	   do 49 k=n,1,-1
	     do 48 its=1,30
	   	 do 41 l=k,1,-1
	      	nm=l-1
	      	if((abs(rv1(l))+anorm).eq.anorm)  goto 2
	      	if((abs(w(nm))+anorm).eq.anorm)  goto 1
41           continue
1            c=0.d0
	   	 s=1.d0
	   	 do 43 i=l,k
	      	f=s*rv1(i)
	      	rv1(i)=c*rv1(i)
	      	if((abs(f)+anorm).eq.anorm) goto 2
	      	g=w(i)
	      	h=pythag(f,g)
	      	w(i)=h
	      	h=1.d0/h
	      	c= (g*h)
	      	s=-(f*h)
	      	do 42 j=1,m
	      	  y=a(j,nm)
	      	  z=a(j,i)
	      	  a(j,nm)=(y*c)+(z*s)
	      	  a(j,i)=-(y*s)+(z*c)
42             continue
43           continue
2            z=w(k)
	   	 if(l.eq.k)then
	      	if(z.lt.0.d0)then
	      	  w(k)=-z
	      	  do 44 j=1,n
	         	 v(j,k)=-v(j,k)
44               continue
	      	endif
	      	goto 3
	   	 endif
	   	 if(its.eq.30) then
				print *, 'no convergence in svdcmp'
				stop
			endif
	   	 x=w(l)
	   	 nm=k-1
	   	 y=w(nm)
	   	 g=rv1(nm)
	   	 h=rv1(k)
	   	 f=((y-z)*(y+z)+(g-h)*(g+h))/(2.d0*h*y)
	   	 g=pythag(f,1.d0)
	   	 f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
	   	 c=1.d0
	   	 s=1.d0
	   	 do 47 j=l,nm
	      	i=j+1
	      	g=rv1(i)
	      	y=w(i)
	      	h=s*g
	      	g=c*g
	      	z=pythag(f,h)
	      	rv1(j)=z
	      	c=f/z
	      	s=h/z
	      	f= (x*c)+(g*s)
	      	g=-(x*s)+(g*c)
	      	h=y*s
	      	y=y*c
	      	do 45 jj=1,n
	      	  x=v(jj,j)
	      	  z=v(jj,i)
	      	  v(jj,j)= (x*c)+(z*s)
	      	  v(jj,i)=-(x*s)+(z*c)
45             continue
	      	z=pythag(f,h)
	      	w(j)=z
	      	if(z.ne.0.d0)then
	      	  z=1.d0/z
	      	  c=f*z
	      	  s=h*z
	      	endif
	      	f= (c*g)+(s*y)
	      	x=-(s*g)+(c*y)
	      	do 46 jj=1,m
	      	  y=a(jj,j)
	      	  z=a(jj,i)
	      	  a(jj,j)= (y*c)+(z*s)
	      	  a(jj,i)=-(y*s)+(z*c)
46             continue
47           continue
	   	 rv1(l)=0.d0
	   	 rv1(k)=f
	   	 w(k)=x
48         continue
3          continue
49       continue
	   return
   end subroutine svdcmp

! ---------------------------------------------------------
! Solves a linear system of equations using the SVD decomposition
! ---------------------------------------------------------	
   subroutine svbksb(u,w,v,m,n,mp,np,b,x)
   integer :: m,mp,n,np
	integer, parameter :: NMAX=500
   real(kind=8) :: b(mp),u(mp,np),v(np,np),w(np),x(np)
   integer :: i,j,jj
   real(kind=8) :: s,tmp(NMAX)
      do 12 j=1,n
        s=0.d0
        if(w(j).ne.0.d0)then
          do 11 i=1,m
            s=s+u(i,j)*b(i)
11        continue
          s=s/w(j)
        endif
        tmp(j)=s
12    continue
      do 14 j=1,n
        s=0.d0
        do 13 jj=1,n
          s=s+v(j,jj)*tmp(jj)
13      continue
        x(j)=s
14    continue
      return
   end subroutine svbksb
		
! ---------------------------------------------------------
! Returns (a^2+b^2)^(1/2) without destructive underflow or overflow
! ---------------------------------------------------------	
	function pythag(a,b)
   real(kind=8) :: a,b,pythag
   real(kind=8) :: absa,absb
      absa=dabs(a)
      absb=dabs(b)
      if(absa.gt.absb)then
        pythag=absa*dsqrt(1.+(absb/absa)**2)
      else
        if(absb.eq.0.d0)then
          pythag=0.d0
        else
          pythag=absb*sqrt(1.d0+(absa/absb)**2)
        endif
      endif
      return
	end function pythag
	
! ---------------------------------------------------------
! Iterative improvement of a solution to a linear set of equations
! ---------------------------------------------------------
	subroutine mprove(a, alud, indx, b, x)
	real(kind=8), INTENT(IN) :: a(:,:), alud(:,:), b(:)
	real(kind=8), INTENT(INOUT) :: x(:)
	integer, INTENT(IN) :: indx(:)
	integer :: i, j, n
	real(kind=8) :: r(size(a,1)), sdp
		n = size(a,1)
		do i = 1, n
			sdp = -b(i)
			do j = 1, n
				sdp = sdp + a(i,j) * x(j)
			enddo
			r(i) = sdp
		enddo
		call lubksb(alud,indx,r)
		
		x = x - r

	end subroutine mprove

! ---------------------------------------------------------
! Solves a symmetric linear system of equations
! ---------------------------------------------------------	
	subroutine llslv(S,X,N,NR)
	integer, INTENT(IN) :: n, nr
	real(kind=8), INTENT(INOUT) :: S(nr,nr), x(nr)
	real(kind=8) :: sum
	integer :: i, j, k
	
! Decompose symmetric matrix L*LT=S
	do i = 1, n
		do j = 1, i
			sum = s(i,j)
			do k = 1, j-1
				sum = sum - S(i,k) * S(j,k)
			enddo !k
			if (j < i) then
				s(i,j) = sum / S(j,j)
			else
				S(i,j) = dsqrt(dabs(sum))
			endif
		enddo !j
	enddo !i
	
	entry llreslv(S,X,N,NR)
! Solve the system
	
	do i = 1, n
		sum = x(i)
		do j = 1, i-1
			sum = sum - S(i,j)*x(j)
		enddo !j
		x(i) = sum / S(i,i)
	enddo !i
	do i = n, 1, -1
		sum = x(i)
		do j = n, i+1, -1
			sum = sum - S(j,i)*x(j)
		enddo !j
		x(i) = sum / S(i,i)
	enddo !i
	
	end subroutine llslv
	
! ---------------------------------------------------------
! Returns a Voigt profile for the line
! ---------------------------------------------------------
      function voigt(a,vv,j)
	  real(kind=8) :: voigt, a, vv
	  integer :: j, k, i, kk, kkk
	  real(kind=8) :: d1, d2, d3, d12, d13, d23, d, y
      real(kind=8) :: a0, a1, a2, a3, a4, a5, a6, b0, b1, b2, b3, b4, b5, b6
	  real(kind=8) :: v
      complex :: z
      real(kind=8) :: xdws(28),ydws(28)

      data a0,a1,a2,a3,a4,a5,a6,b0,b1,b2,b3,b4,b5,b6/&
      122.607931777104326d0,214.382388694706425d0,181.928533092181549d0,&
      93.155580458138441d0,30.180142196210589d0,5.912626209773153d0,&
      .564189583562615d0,122.60793177387535d0,352.730625110963558d0,&
      457.334478783897737d0,348.703917719495792d0,170.354001821091472d0,&
      53.992906912940207d0,10.479857114260399d0/

      data xdws/.1d0,.2d0,.3d0,.4d0,.5d0,.6d0,.7d0,.8d0,.9d0,1.d0,1.2d0,1.4d0,1.6d0,1.8d0,2.d0,&
      3.d0,4.d0,5.d0,6.d0,7.d0,8.d0,9.d0,10.d0,12.d0,14.d0,16.d0,18.d0,20.d0/,ydws/&
      9.9335991d-02,1.9475104d-01,2.8263167d-01,3.5994348d-01,&
      4.2443639d-01,4.7476321d-01,5.1050407d-01,5.3210169d-01,&
      5.4072434d-01,5.3807950d-01,5.0727350d-01,4.5650724d-01,&
      3.9993989d-01,3.4677279d-01,3.0134040d-01,1.7827103d-01,&
      1.2934799d-01,1.0213407d-01,8.4542692d-02,7.2180972d-02,&
      6.3000202d-02,5.5905048d-02,5.0253846d-02,4.1812878d-02,&
      3.5806101d-02,3.1311397d-02,2.7820844d-02,2.5031367d-02/

      v=dabs(vv)
      IF(A.NE.0) GOTO 1 
      IF(J.NE.0) GOTO 3 
      VOIGT=DEXP(-V*V)   
      RETURN
   3  IF(V.GT.XDWS(1)) GOTO 4   
      D=V*(1.-.66666667d0*V*V)
      GOTO 8
   4  IF(V.GT.XDWS(28)) GOTO 5  
      K=27  
      DO 7 I=2,27   
      IF(XDWS(I).LT.V) GOTO 7   
      K=I   
      GOTO 6
   7  CONTINUE  
   6  KK=K-1
      KKK=K+1   
      D1=V-XDWS(KK) 
      D2=V-XDWS(K)  
      D3=V-XDWS(KKK)
      D12=XDWS(KK)-XDWS(K)  
      D13=XDWS(KK)-XDWS(KKK)
      D23=XDWS(K)-XDWS(KKK) 
      D=YDWS(KK)*D2*D3/(D12*D13)-YDWS(K)*D1*D3/(D12*D23)+YDWS(KKK)*&
        D1*D2/(D13*D23)
      GOTO 8
   5  Y=.5/V
      D=Y*(1.+Y/V)  
   8  VOIGT=5.641895836d-1*D
   9  IF(VV.LT.0.) VOIGT=-VOIGT 
      RETURN
   1  Z=CMPLX(A,-V) 
      Z=((((((A6*Z+A5)*Z+A4)*Z+A3)*Z+A2)*Z+A1)*Z+A0)/&
      (((((((Z+B6)*Z+B5)*Z+B4)*Z+B3)*Z+B2)*Z+B1)*Z+B0)
      IF(J.NE.0) GOTO 2 
      VOIGT=REAL(Z) 
      RETURN
   2  VOIGT=.5d0*AIMAG(Z) 
      GOTO 9
      END function voigt
		
! ---------------------------------------------------------
! Given x(:) and y(:) which tabulate a function and the derivative at the boundary points
! this function returns the second derivative of the spline at each point
! ---------------------------------------------------------
		subroutine splin1(x,y,yp1,ypn,y2)
		real(kind=8), INTENT(IN) :: x(:), y(:), yp1, ypn
		real(kind=8), INTENT(INOUT) :: y2(size(x))
		integer :: n, i, k
		real(kind=8) :: p, qn, sig, un, u(size(x))

			n = size(x)
			
			if (yp1 > .99d30) then
				y2(1) = 0.d0
				u(1) = 0.d0
			else
				y2(1) = -0.5d0
				u(1) = (3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
			endif

			do i = 2, n-1
				sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))				
				p = sig * y2(i-1)+2.d0
				y2(i) = (sig-1.d0)/p
				u(i) = (6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/&
					(x(i+1)-x(i-1))-sig*u(i-1))/p
			enddo
			if (ypn > .99d30) then
				qn = 0.d0
				un = 0.d0
			else
				qn = 0.5d0
				un = (3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
			endif
			
			y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.d0)

			do k = n-1, 1, -1
				y2(k) = y2(k)*y2(k+1)+u(k)
			enddo

		end subroutine splin1

! ---------------------------------------------------------
! Given xa(:) and ya(:) which tabulate a function, returns the interpolation using
! splines of vector x(:) in y(:)
! ---------------------------------------------------------
		subroutine spline(xa,ya,x,y,extrapolation)
		real(kind=8), INTENT(INOUT) :: y(:)
		real(kind=8), INTENT(IN) :: xa(:), ya(:), x(:)
		real(kind=8) :: y2a(size(xa))
		character*4 :: extrapolation
		integer :: n_x, n, i, k, khi, klo
		real(kind=8) :: a, b, h, extrap
			
			n = size(xa)
			n_x = size(x)
			call splin1(xa,ya,1.d30,1.d30,y2a)

			do i = 1, n_x					

! Downward extrapolation 
				if (x(i) < xa(1)) then
					select case(extrapolation)
						case('NO  ')
							y(i) = ya(1)
						case('SQRT')
							extrap = mean(ya/dsqrt(xa))
							y(i) = extrap * x(i)
					end select
				else 

! Upward extrapolation
				if (x(i) > xa(n)) then
					select case(extrapolation)
						case('NO  ')
							y(i) = ya(n)
						case('SQRT')
							extrap = mean(ya/dsqrt(xa))
							y(i) = extrap * x(i)
					end select
				else
! In range
						klo = 1
						khi = n
1						if(khi-klo > 1) then
							k = (khi+klo)/2
							if (xa(k) > x(i)) then
								khi = k
							else
								klo = k
							endif					
							go to 1
						endif

						h = xa(khi)-xa(klo)

						if (h == 0.d0) then
							print *, 'bad xa input in spline'
							stop
						endif
						a = (xa(khi)-x(i))/h
						b = (x(i)-xa(klo))/h

						y(i) = a*ya(klo)+b*ya(khi)+((a**3.d0-a)*y2a(klo)+(b**3.d0-b)*y2a(khi))*(h**2.d0)/6.d0		
					endif
				endif
			enddo

		end subroutine spline

! ---------------------------------------------------------
! Returns the value at x0 with an interpolation in the data (x,y)
! ---------------------------------------------------------
		function interpolate(x,y,x0)
		real(kind=8) :: interpolate
		real(kind=8) :: x0
		real(kind=8), INTENT(IN) :: x(:), y(:)
		real(kind=8) :: x1(1), y1(1)
			x1(1) = x0
			call spline(x,y,x1,y1,'NO  ')
			interpolate = y1(1)				
		end function interpolate

! ---------------------------------------------------------
! Returns the vector of linearly interpolated data for a given data x0, y0 at given values x
! ---------------------------------------------------------
		subroutine linear_interpol(x0,y0,x,y)	
		real(kind=8), INTENT(IN) :: x0(:), y0(:)
		real(kind=8), INTENT(INOUT) :: x(:), y(:)
		real(kind=8) :: a, b, ya, yb, m, n
		integer :: n1, n2, i, j, which

			which = 0
			
			n1 = size(x0)
			n2 = size(x)

			do i = 1, n2
				if (x(i) < x0(1)) then
					which = 1
				else
					do j = 1, n1-1
						if (x(i) >= x0(j) .and. x(i) <= x0(j+1)) then
							which = j
						endif
					enddo
				endif
				a = x0(which)
				b = x0(which+1)
				ya = y0(which)
				yb = y0(which+1)
				
				m = (ya-yb)/(a-b)
				n = (yb*a-ya*b)/(a-b)
				y(i) = m * x(i) + n
			enddo

						
		end subroutine linear_interpol

! ---------------------------------------------------------
! Returns the imaginary part of a complex number
! ---------------------------------------------------------
		function imaginary(x)		
		real(kind=8) :: imaginary
		complex(kind=8), INTENT(IN) :: x
		complex(kind=8) :: i
			i = cmplx(0.d0,1.d0)
			imaginary = -real(i*x)
		end function imaginary
		
		
end module maths
