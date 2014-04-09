! *********************************************************
! *********************************************************
! ZEEMAN VARIABLES
! *********************************************************
! *********************************************************

module variables_zeeman
implicit none
! wave: central wavelength of the transition
! nw: number of wavelengths
! wstep: step in wavelength in mA
! ndir: number of directions for angular integration
! zmax: maximum z
! zmin: minimum z
! dw: frequency displacement axis
! tric: value to turn on series in exp(-t)
! ditype: type of the transition (electric/magnetic dipole)
! magop: magnetooptical or birrefringence effects
! nz: number of depth points
! z: depth points
! etaz: line/continuum opacity
! aval: a of Voigt
! dwid: Doppler width in mA
! sl: source function at each point
! planck: Planck's function
! bmag: magnetic field in G
! gamma: inclination of B in degrees
! ch: azimuth of B in degrees
! vel: macroscopic velocity in km/s
! eps: epsilon
! ytr: frequency weights for the red wing
! ytp: frequency weights for the pi wing
! ytb: frequency weights for the blue wing
! kentry: tells if Zeeman shifts have to be computed at vgtgen
! d: Zeeman shifts
! s: Zeeman strengths
! a: Voigt parameter a
! vlos: line of sight velocity
! magf: magnetic field
! dopp: Doppler shift
! sline: local line source function
! bnu: local planck's function
! coz: cosine of the z angular quadrature
! cox: cosine of the x angular quadrature
! wtdir: weights
! phij: value of phij

	real(kind=8) :: wave, wstep, zmax, zmin, tric, cont, boll
	real(kind=8), allocatable :: dw(:), ytr(:), ytp(:), ytb(:)
	real(kind=8), allocatable :: wtr(:,:), wtp(:,:), wtb(:,:)
	real(kind=8), allocatable :: z(:), etaz(:), aval(:), dwid(:), Sli(:), planc(:), bmag(:)
	real(kind=8), allocatable :: gamma(:), ch(:), vel(:), eps(:)
	real(kind=8), allocatable :: coz(:), cox(:), wtdir(:), bnu(:)
   real(kind=8), allocatable :: phij(:), planus(:), plunus(:)
	real(kind=8), allocatable :: emerg(:,:,:), phi1(:), phi2(:), phi3(:), phi4(:)
	integer :: nw, ndir, nz
	real(kind=8) :: ru, lu, ju, rl, ll, jl
   character*1 :: ditype, magop
   character*132 :: nombre
	integer :: kentry
	real(kind=8) :: d(3,20), s(3,20)
	integer :: ntr(3)
   real(kind=8) :: a, vlos, magf, dopp, eta0, sline
	real(kind=8), parameter :: ONE_SQRTPI = 0.564189583546d0, PI_180 = 1.74532925199d-2
	
	integer :: magneto_optics
	real(kind=8) :: damping, mu_direction, inclination, azimuth, field_strength, mu_read
	character*8 :: type_transition
	character*6 :: type_zeeman
	character*4 :: boundary_condition
	character*40 :: output_file
	
end module variables_zeeman

! *********************************************************
! *********************************************************
! ZEEMAN SPECIFIC FUNCTIONS
! *********************************************************
! *********************************************************

module functions_delopar
use variables_zeeman
use variables
use general

implicit none
contains

!--------------------------------------------------------------
! Returns the signum of a quantity
!--------------------------------------------------------------
function sign(a)
real(kind=8) :: sign, a
	sign = a / dabs(a)
end function sign

!--------------------------------------------------------------
! Returns 0 if it is positive and 1 if it is negative
!--------------------------------------------------------------
function negative(a)
real(kind=8) :: negative, a
	if (a < 0.d0) then
		negative = 1.d0
	else
		negative = 0.d0
	endif
end function negative

!--------------------------------------------------------------
! Generate the correct opacity and distances for a plane-parallel atmosphere
!--------------------------------------------------------------
subroutine generate_zeeman_atmos_pp
integer :: temp, up, low, i
real(kind=8) :: einsteinA

! Upper and lower level of the transition	
	up = itran(1,zeeman_line)
	low = itran(2,zeeman_line)

! First calculate the opacity and source function of the transition	
	call opacity(zeeman_line,temp,1,n_radios)
	nz = n_radios
			
	allocate(z(nz))
! This rotation is due to the fact that the vector of heights have to begin at the deeper
! zones, and the r vector is flipped	
	z = rotate(r)
	
	allocate(etaz(nz))
	etaz = rotate(chil)
	
	allocate(aval(nz))
	einsteinA = dtran(1,zeeman_line) * (dlevel(2,low) / dlevel(2,up)) * &
		(PI8C2 * dtran(2,zeeman_line)**2.d0) 
	
! Natural broadening of the line
	aval = (wave*1.d-8) * einsteinA / (4.d0 * PI * thermal_v)
	aval = rotate(aval)
	
! Doppler width in mA
	allocate(dwid(nz))	
	dwid = ((wave*1.d-8)**2.d0 * doppler_width(:,zeeman_line) / PC ) * 1.d11
	dwid = rotate(dwid)
	
	
	allocate(planc(nz))
	do i = 1, nz
		planc(i) = planck(PC/(wave*1.d-8),datph(1,i))
	enddo
	planc = rotate(planc)

	allocate(Sli(nz))
	Sli = rotate(Sl)
		
	allocate(bmag(nz))
! Use the field in the atmosphere model file or use a constant field
	if (field_strength < 0.d0) then		
		bmag = rotate(datph(9,:))
	else
		bmag = field_strength
	endif
	
	allocate(gamma(nz))
! Use the field in the atmosphere model file or use a constant field inclination
	if (dabs(inclination) > 180.d0) then
		gamma = rotate(datph(10,:))
	else
		gamma = inclination
	endif
	
! Use the field in the atmosphere model file or use a constant field azimuth
	allocate(ch(nz))	
	if (dabs(azimuth) > 180.d0) then
		ch = rotate(datph(11,:))
	else
		ch = azimuth
	endif
	
! Velocity field
	allocate(vel(nz))
	vel = rotate(datph(7,:))
	
	allocate(eps(nz))
	allocate(phij(nz))
   allocate(bnu(nz))
   allocate(phi1(nz))
   allocate(phi2(nz))
   allocate(phi3(nz))
   allocate(phi4(nz))

	if (background_flag /= 0) then
! REVISAR ESTO
		cont = kappa(1)
	else
		cont = 0.d0
	endif
	cont = 0.d0
	
end subroutine generate_zeeman_atmos_pp

!--------------------------------------------------------------
! Generate the correct opacity and distances for a spherical atmosphere
!--------------------------------------------------------------
subroutine generate_zeeman_atmos_spher	
integer :: temp, up, low, i, p_corte_origen, cortes_caract
real(kind=8) :: einsteinA, p_salida

	mu_direction = mu_read
	
	p_salida = r(n_radios) * dsqrt(1.d0-mu_direction**2.d0)
	p_corte_origen = 0.d0
	
	do i = 1, n_radios
		if (r(i) <= p_salida) then
			p_corte_origen = p_corte_origen + 1
		endif
	enddo
	
	
	if (p_corte_origen /= 0) then 
		cortes_caract = 2*(n_radios - p_corte_origen)	
	else 
		cortes_caract = (n_radios - p_corte_origen)	
	endif

	
! Upper and lower level of the transition	
	up = itran(1,zeeman_line)
	low = itran(2,zeeman_line)

! First calculate the opacity and source function of the transition	
	call opacity(zeeman_line,temp,1,n_radios)
	nz = cortes_caract
	print *, 'Number of grid points for the selected radius : ', nz	
			
! Height vector
	allocate(z(nz))
	if (p_corte_origen /= 0) then
		boundary_condition = 'ZERO'
		do i = 1, nz/2
			z(i) = -dsqrt( r(n_radios-i+1)**2.d0 - p_salida**2.d0)
		enddo
		do i = 1,nz/2
			z(i+nz/2) = dsqrt( r(p_corte_origen+i)**2.d0 - p_salida**2.d0)
		enddo
		z = rotate(z)
	else
		boundary_condition = 'DIFF'	
		do i = 1, nz
			z(i) = dsqrt( r(i)**2.d0 - p_salida**2.d0)
		enddo
		z = rotate(r)
	endif

			
! Opacity vector
	allocate(etaz(nz))
	if (p_corte_origen /= 0) then
		do i = 1, nz/2
			etaz(i) = chil(n_radios-i+1)
		enddo
		do i = 1,nz/2
			etaz(i+nz/2) = chil(p_corte_origen+i)
		enddo
		etaz = rotate(etaz)
	else
		etaz = rotate(chil)
	endif

	
	allocate(aval(nz))
	einsteinA = dtran(1,zeeman_line) * (dlevel(2,low) / dlevel(2,up)) * &
		(PI8C2 * dtran(2,zeeman_line)**2.d0) 
	
! Natural broadening of the line
	if (p_corte_origen /= 0) then
		do i = 1, nz/2
			aval(i) = (wave*1.d-8) * einsteinA / (4.d0 * PI * thermal_v(n_radios-i+1))
		enddo
		do i = 1,nz/2
			aval(i+nz/2) = (wave*1.d-8) * einsteinA / (4.d0 * PI * thermal_v(p_corte_origen+i))
		enddo
		aval = rotate(aval)
	else
		aval = (wave*1.d-8) * einsteinA / (4.d0 * PI * thermal_v)
		aval = rotate(aval)
	endif
	
! Doppler width in mA
	allocate(dwid(nz))	
	if (p_corte_origen /= 0) then
		do i = 1, nz/2
			dwid(i) = ((wave*1.d-8)**2.d0 * doppler_width(n_radios-i+1,zeeman_line) / PC ) * 1.d11
		enddo
		do i = 1,nz/2
			dwid(i+nz/2) = ((wave*1.d-8)**2.d0 * doppler_width(p_corte_origen+i,zeeman_line) / PC ) * 1.d11
		enddo
		dwid = rotate(dwid)
	else
		dwid = ((wave*1.d-8)**2.d0 * doppler_width(:,zeeman_line) / PC ) * 1.d11
		dwid = rotate(dwid)
	endif
	
! Planck's function	
	allocate(planc(nz))
	if (p_corte_origen /= 0) then
		do i = 1, nz/2
			planc(i) = planck(PC/(wave*1.d-8),datph(1,n_radios-i+1))
		enddo
		do i = 1,nz/2
			planc(i+nz/2) = planck(PC/(wave*1.d-8),datph(1,p_corte_origen+i))
		enddo
		planc = rotate(planc)
	else
		do i = 1, nz
			planc(i) = planck(PC/(wave*1.d-8),datph(1,i))
		enddo
		planc = rotate(planc)
	endif	

! Source function
	allocate(Sli(nz))
	if (p_corte_origen /= 0) then
		do i = 1, nz/2
			Sli(i) = Sl(n_radios-i+1)
		enddo
		do i = 1,nz/2
			Sli(i+nz/2) = Sl(p_corte_origen+i)
		enddo
		Sli = rotate(Sli)
	else
		do i = 1, nz
			Sli(i) = Sl(i)
		enddo
		Sli = rotate(Sli)
	endif
			
	allocate(bmag(nz))
! Use the field in the atmosphere model file or use a constant field
	if (field_strength < 0.d0) then		
		if (p_corte_origen /= 0) then
			do i = 1, nz/2
				bmag(i) = datph(9,n_radios-i+1)
			enddo
			do i = 1,nz/2
				bmag(i+nz/2) = datph(9,p_corte_origen+i)
			enddo
			bmag = rotate(bmag)
		else
			do i = 1, nz
				bmag(i) = datph(9,i)
			enddo
			bmag = rotate(bmag)
		endif
	else
		bmag = field_strength
	endif
	
	allocate(gamma(nz))
! Use the field in the atmosphere model file or use a constant field inclination
	
	if (dabs(inclination) > 180.d0) then

		if (p_corte_origen /= 0) then
			do i = 1, nz/2
				gamma(i) = PI * negative(z(i)) - (datph(10,n_radios-i+1) + &
					acos(dsqrt( r(n_radios-i+1)**2.d0 - p_salida**2.d0) / r(n_radios-i+1))/PI_180)
			enddo
			do i = 1,nz/2
				gamma(i+nz/2) = PI * negative(z(i)) - (datph(10,p_corte_origen+i) + &
				  acos(dsqrt( r(p_corte_origen+i)**2.d0 - p_salida**2.d0) / r(p_corte_origen+i))/PI_180)
			enddo
			gamma = rotate(gamma)
		else
			do i = 1, nz
				gamma(i) = datph(10,i) + &
					acos(dsqrt( r(i)**2.d0 - p_salida**2.d0) / r(i))/PI_180
			enddo
			gamma = rotate(gamma)
		endif		
	else
		if (p_corte_origen /= 0) then
			do i = 1, nz/2
				gamma(i) = inclination + &
					acos(dsqrt( r(n_radios-i+1)**2.d0 - p_salida**2.d0) / r(n_radios-i+1))/PI_180
			enddo
			do i = 1,nz/2
				gamma(i+nz/2) = inclination + &
					acos(dsqrt( r(p_corte_origen+i)**2.d0 - p_salida**2.d0) / r(p_corte_origen+i))/PI_180
			enddo
			gamma = rotate(gamma)
		else
			do i = 1, nz
				gamma(i) = inclination + &
					acos(dsqrt( r(i)**2.d0 - p_salida**2.d0) / r(i))/PI_180
			enddo
			gamma = rotate(gamma)
		endif		
	endif
	
	
! Use the field in the atmosphere model file or use a constant field azimuth
	allocate(ch(nz))	
	if (dabs(azimuth) > 180.d0) then
		if (p_corte_origen /= 0) then
			do i = 1, nz/2
				ch(i) = datph(11,n_radios-i+1) 
			enddo
			do i = 1,nz/2
				ch(i+nz/2) = datph(11,p_corte_origen+i) 
			enddo
			ch = rotate(ch)
		else
			do i = 1, nz
				ch(i) = datph(11,i)
			enddo
			ch = rotate(ch)
		endif			
	else
		ch = azimuth
	endif
	
! Velocity field
	allocate(vel(nz))
	if (p_corte_origen /= 0) then
		do i = 1, nz/2
			vel(i) = datph(7,n_radios-i+1) 
		enddo
		do i = 1,nz/2
			vel(i+nz/2) = datph(7,p_corte_origen+i) 
		enddo
		vel = rotate(vel)
	else
		do i = 1, nz
			vel(i) = datph(7,i)
		enddo
		vel = rotate(vel)
	endif			
	
	allocate(eps(nz))
	allocate(phij(nz))
   allocate(bnu(nz))
   allocate(phi1(nz))
   allocate(phi2(nz))
   allocate(phi3(nz))
   allocate(phi4(nz))

	if (background_flag /= 0) then
! REVISAR ESTO
		cont = kappa(1)
	else
		cont = 0.d0
	endif
	cont = 0.d0
	
	mu_direction = 1.d0

end subroutine generate_zeeman_atmos_spher	


!--------------------------------------------------------------
! Read the atmosphere and the atomic model
!--------------------------------------------------------------
subroutine initialize_vars
	integer :: i, k, up, low
	character*1 :: lup, llo
	character*1 :: term(6)
	real(kind=8) :: divr, divb, divp, delw
	real(kind=8) :: hp, fvp, hb, fvb, hr, fvr
	data term/'S','P','D','F','G','H'/

	nombre = 'ATMOS'

! Put the frequency in Angstroms
	wave = dtran(2,zeeman_line)
	wave = PC * 1.d8 / wave
	select case(type_zeeman)
		case('ATOMIC')
			print *, 'Transition : ', wave, ' �'
		case('MOLROT')
			print *, 'Transition : ', dtran(2,zeeman_line)/1.e9,' GHz'
	end select 
	
! Upper and lower level of the transition	
	up = itran(1,zeeman_line)
	low = itran(2,zeeman_line)
		
	select case(type_zeeman) 

! Atomic transitions
		case('ATOMIC')

		! Get the multiplicity from the level label
			read(label(up)(1:1),*) ru
			read(label(low)(1:1),*) rl

		! Get the angular moment L	
			read(label(up)(2:2),*) lup
			read(label(low)(2:2),*) llo

		! Get the angular moment J
			read(label(up)(3:3),*) ju
			read(label(low)(3:3),*) jl
			ju = (ju - 1.d0) / 2.d0
			jl = (jl - 1.d0) / 2.d0

			print *, 'Upper level : ', up, '  term   ', label(up)(1:3)
			print *, '     Orbital angular moment : ', lup
			print *, '     Spin angular moment : ', (ru-1.d0)/2.d0
			print *, '     Total angular moment : ', ju
			print *, 'Lower level : ', low, '  term   ', label(low)(1:3)	
			print *, '     Orbital angular moment : ', llo
			print *, '     Spin angular moment : ', (rl-1.d0)/2.d0
			print *, '     Total angular moment : ', jl
			
! Molecular rotation transitions
		case('MOLROT')

		! Get the angular moment J
			read(label(up)(3:3),*) ju
			read(label(low)(3:3),*) jl

			print *, 'Upper level : ', up
			print *, '     Total angular moment : ', ju
			print *, 'Lower level : ', low
			print *, '     Total angular moment : ', jl				
		
	end select	
		
! Type of the transition
	if (type_transition == 'ELECTRIC') then
		ditype = 'e'
	else
		ditype = 'm'
	endif
	
! Magneto-optical effects
	if (magneto_optics == 1) then
		magop = 'y'
	else
		magop = 'n'
	endif
	
	ndir = 1
	
	allocate(dw(nw))
	allocate(cox(ndir))
	allocate(coz(ndir))
	allocate(wtdir(ndir))
	
	zmax = maxval(r)
	zmin = minval(r)
	
	if (zmax < zmin) stop

! Continuum
	dw(1) = 20000.

! Generate the frequency axis
	do i = 2, nw
		dw(i) = dble(i-2.d0)*wstep - dble((nw/2.d0)-1.d0)*wstep
	enddo

! Now put the atmosphere into the correct arrays
	
	if (geometry_type == 'PLANEP') then
		call generate_zeeman_atmos_pp
	else if (geometry_type == 'SPHERICAL') then 
		call generate_zeeman_atmos_spher	
	endif
	
		
! Test the angular momentum of the upper and lower level
	i = 0
	do while (term(i+1) /= lup)
		i = i + 1
	enddo
	lu = i
	
	i = 0
	do while (term(i+1) /= llo)
		i = i + 1
	enddo
	ll = i	
	
! Renormalization in frequency
	allocate(ytr(nw))
	allocate(ytp(nw))
	allocate(ytb(nw))
	allocate(wtr(nw,nz))
	allocate(wtp(nw,nz))
	allocate(wtb(nw,nz))
	
	ytr = wstep
	ytp = wstep
	ytb = wstep
	ytr(2) = wstep * 0.5d0
	ytr(nw) = wstep * 0.5d0
	ytp(2) = wstep * 0.5d0
	ytp(nw) = wstep * 0.5d0	
	ytb(2) = wstep * 0.5d0
	ytb(nw) = wstep * 0.5d0	
	
	do k = 1, nz
		divr = 0.d0
		divp = 0.d0
		divb = 0.d0
		a = aval(k)
		dopp = dwid(k)
		magf = bmag(k)
		vlos = vel(k)
		do i = 2, nw
			delw = dw(i)
			call vgtgen(a,delw,hp,fvp,hb,fvb,hr,fvr)

			divr = divr + ytr(i)*hr
			divp = divp + ytp(i)*hp
			divb = divb + ytb(i)*hb
		enddo

! Renormalize the weights
		wtr(:,k) = ytr / divr
		wtp(:,k) = ytp / divp
		wtb(:,k) = ytb / divb				
		
	enddo

! Get the angular quadrature	
	call angle(ndir,wtdir,coz,cox)	
		
end subroutine initialize_vars

!--------------------------------------------------------------
! Generate the angle quadrature points and weights
!--------------------------------------------------------------
subroutine angle(na, wa, cz, cx)
integer :: na
real(kind=8) :: wa(na), cz(na), cx(na), sam
	select case(na) 
! 1-point quadrature	
		case(1)
			wa(1) = 0.5d0
			cz(1) = dsqrt(1.d0/3.d0)
			cx(1) = dsqrt(1.d0/3.d0)

! 6-points quadrature			
		case(6)
			cz(1) = 0.80847425d0
			cz(2) = 0.80847425d0			
			cz(3) = 0.57735027d0
			cz(4) = 0.57735027d0
			cz(5) = 0.11417547d0
			cz(6) = 0.11417547d0
			
			cx(1) = 0.11417547d0
			cx(2) = 0.57735027d0
			cx(3) = 0.11417547d0
			cx(4) = 0.80847425d0
			cx(5) = 0.57735027d0
			cx(6) = 0.80847425d0
			
			wa = 1.d0 / 6.d0
			
! Set S(b) of Carlsson for 12 directions per octant
		case(12)
			cz(1:3) = 0.11104445d0
			cz(4:7) = 0.50307327d0
			cz(8:10) = 0.70273364d0
			cz(11:12) = 0.85708017d0
			
			cx(1) = 0.50307327d0
			cx(2) = 0.70273364d0
			cx(3) = 0.85708017d0
			cx(4) = 0.11104445d0
			cx(5) = 0.50307327d0
			cx(6) = 0.70273364d0
			cx(7) = 0.85708017d0
			cx(8) = 0.11104445d0
			cx(9) = 0.50307327d0
			cx(10) = 0.70273364d0
			cx(11) = 0.11104445d0
			cx(12) = 0.50307327d0
			
			wa(1) = 0.11434821d0
			wa(2) = 0.07937696d0
			wa(3:4) = 0.11434821d0
			wa(4:6) = 0.02525996d0
			wa(7) = 0.11434821d0
			wa(8) = 0.07937696d0
			wa(9) = 0.02525996d0
			wa(10) = 0.07937696d0
			wa(11:12) = 0.11434821d0
	end select 
	
	if (na /= 1 .and. na /= 6 .and. na /= 12) then
		print *, 'Wrong angular quadrature'
		stop
	endif
! Renormalization of the angular quadrature
	
	sam = 2.d0 * sum(wa)
	wa = wa / sam
	
end subroutine angle

!--------------------------------------------------------------
! Generates the Voigt and anomalous dispersion profiles from dipole transition
! (ru,lu,ju) -> (rl,ll,jl)  Dipole type is specified by ditype
! r(u/l) -> multiplicity
!--------------------------------------------------------------
subroutine vgtgen(a, delw, hp, fvp, hb, fvb, hr, fvr)
real(kind=8) :: a, delw, hp, fvp, hb, fvb, hr, fvr, cte
real(kind=8) :: vshift, v, hgen(3), fgen(3), h, f, dmag
integer :: i, k

	if (kentry == 0) then
		kentry = 1

		select case(type_zeeman)
			case('ATOMIC')
				call zeem_atomic(d,s,ntr)
			case('MOLROT')
				call zeem_molrot(d,s,ntr)
		end select
		

! Convert Zeeman displacements to mA/G by multiplying by e*lambda**2/(4*pi*m*c^2)	
      if (nombre /= 'MEDD') then
! Constant to give wavelength in A		
			cte = PEC / (4.d0 * PI * PME * PC**2) * 1.d-5
			d = d * cte * wave**2
      endif
	endif
	

! Velocity shift in mA
	vshift = wave*vlos/299.792458d0
   
! Doppler shift
	v = (delw + vshift) / dopp
	
! Magnetic shift in Doppler units
	dmag = magf / dopp
	hgen = 0.d0
	fgen = 0.d0
	do i = 1, 3
		do k = 1, ntr(i)
			call fvoigt(a,v-dmag*d(i,k),h,f)
			if (magop /= 'y') f = 0.d0
			hgen(i) = hgen(i) + s(i,k)*h
			fgen(i) = fgen(i) + s(i,k)*f
		enddo
	enddo

! Using Landi-type notation. i = 1,2,3 means (mu-ml)=+1,0,-1 (indices b,p,r)
! Multiply by 1/sqrt(pi) to get normalized profiles
	hp = hgen(2)*ONE_SQRTPI
	fvp = fgen(2)*ONE_SQRTPI	
	hb = hgen(1)*ONE_SQRTPI
	fvb = fgen(1)*ONE_SQRTPI	
	hr = hgen(3)*ONE_SQRTPI
	fvr = fgen(3)*ONE_SQRTPI			
	
end subroutine vgtgen

!--------------------------------------------------------------
! Generates Voigt and anomalous dispersion profiles
! See Humlicek (1982) JQSRT 27, 437
!--------------------------------------------------------------
subroutine fvoigt(da, dv, dh, dfv)
real(kind=8) :: da, dv, dh, dfv
	complex :: w4, z, t, u, v4
	real(kind=8) :: s
	z = cmplx(dv, da)
	t = cmplx(da, -dv)
	s = dabs(dv) + da
	u = t*t
	
	if (s >= 15.d0) then
		w4 = t * 0.5641896d0 / (0.5d0+u)
	elseif (s >= 5.5) then
			w4 = t*(1.410474d0+u*0.5641896d0)/(0.75d0+u*(3.d0+u))
		elseif (da >= 0.195d0*dabs(dv)-0.176d0) then
				w4 = (16.4955d0+t*(20.20933d0+t*(11.96482d0+t*(3.778987d0+t*0.5642236d0)))) / &
					(16.4955d0+t*(38.82363d0+t*(39.27121d0+t*(21.69274d0+t*(6.699398d0+t)))))
			else 
				w4 = t*(36183.31d0-u*(3321.9905d0-u*(1540.787d0-u*(219.0313d0-u*(35.76683d0-&
					u*(1.320522d0-u*0.56419d0))))))
				v4 = (32066.6d0-u*(24322.84d0-u*(9022.228d0-u*(2186.181d0-u*(364.2191d0-&
					u*(61.57037d0-u*(1.841439d0-u)))))))
				w4 = cexp(u) - w4/v4
			endif

			
	dh = dble(w4)
	dfv = aimag(w4)

end subroutine fvoigt

!--------------------------------------------------------------
! Calculates the Zeeman shifts and strengths for an atomic transition
!--------------------------------------------------------------
subroutine zeem_atomic(d, s, ntr)
real(kind=8) :: d(3,20), s(3,20)
integer :: ntr(3)
real(kind=8) :: mlo, mup, ml(20), mu(20), su, sl, gu, gl, delm, sums
integer :: nu, nl, k, n, i

! Upper level
	su = (ru-1.d0)/2.d0   ! Spin taken from the multiplicity
	gu = 0.d0
	if (ju /= 0.d0) then
		gu = 1.5d0+(su*(su+1.d0)-lu*(lu+1.d0)) / (2.d0*ju*(ju+1.d0))   ! Lande factor
	endif
	nu = int(2.d0*ju+1.d0)
	do k = 1, nu
		mu(k) = ju+1.d0-k
	enddo
	
! Lower level
	sl = (rl-1.d0)/2.d0   ! Spin taken from the multiplicity
	gl = 0.d0
	if (jl /= 0.d0) then
		gl = 1.5d0+(sl*(sl+1.d0)-ll*(ll+1.d0)) / (2.d0*jl*(jl+1.d0))   ! Lande factor
	endif
	nl = int(2.d0*jl+1.d0)
	do k = 1, nl
		ml(k) = jl+1.d0-k
	enddo
	n = min(nu,nl)
	ntr = n
	do i = 1, 3
		delm = dble(2-i)   ! delta(M) = -1,0,+1
		if ((ju == jl) .and. (delm /= 0.d0)) ntr(i) = n-1
		sums = 0.d0
		do k = 1, ntr(i)
			if ((ju > jl) .or. ((ju == jl) .and. (delm <= 0.d0))) then
				mup = ml(k) + delm
				mlo = ml(k)
			else
				mup = mu(k)
				mlo = mu(k) - delm
			endif
			d(i,k) = gl * mlo - gu * mup
			s(i,k) = str(ju,mup,jl,mlo)
		enddo	
		sums = sum(s(i,:))
		s(i,:) = s(i,:) / sums
	enddo
	
! Test for a triplet
	if ((ju == 0.d0) .or. (jl == 0.d0) .or. (gu == gl)) then	
		print *, 'Triplet detected'
		ntr = 1
		s(:,1) = 1.d0
	endif
			
end subroutine zeem_atomic

!--------------------------------------------------------------
! Calculates the Zeeman shifts and strengths for an atomic transition
!--------------------------------------------------------------
subroutine zeem_molrot(d, s, ntr)
real(kind=8) :: d(3,20), s(3,20)
integer :: ntr(3)
real(kind=8) :: mlo, mup, ml(20), mu(20), gu, gl, delm, sums
integer :: nu, nl, k, n, i

! Upper level
	gu = 0.d0
	if (ju /= 0.d0) then
! Only for CO
		gu = -0.268d0  ! Lande factor
	endif
	nu = int(2.d0*ju+1.d0)
	do k = 1, nu
		mu(k) = ju+1.d0-k
	enddo
	
! Lower level
	gl = 0.d0
	if (jl /= 0.d0) then
		gl = -0.268d0  ! Lande factor
	endif
	nl = int(2.d0*jl+1.d0)
	do k = 1, nl
		ml(k) = jl+1.d0-k
	enddo
	n = min(nu,nl)
	ntr = n
	do i = 1, 3
		delm = dble(2-i)   ! delta(M) = -1,0,+1
		if ((ju == jl) .and. (delm /= 0.d0)) ntr(i) = n-1
		sums = 0.d0
		do k = 1, ntr(i)
			if ((ju > jl) .or. ((ju == jl) .and. (delm <= 0.d0))) then
				mup = ml(k) + delm
				mlo = ml(k)
			else
				mup = mu(k)
				mlo = mu(k) - delm
			endif
			d(i,k) = gl * mlo - gu * mup
			s(i,k) = str(ju,mup,jl,mlo)
		enddo	
		sums = sum(s(i,:))
		s(i,:) = s(i,:) / sums
	enddo
	
! Test for a triplet
	if ((ju == 0.d0) .or. (jl == 0.d0) .or. (gu == gl)) then	
		print *, 'Triplet detected'
		ntr = 1
		s(:,1) = 1.d0
	endif
			
end subroutine zeem_molrot

!--------------------------------------------------------------
! Un-normalized Zeeman strengths
! ju,mu -> J and M of the upper level
! jl,ml -> J and M of the lower level
!--------------------------------------------------------------
function str(ju,mu,jl,ml)
real(kind=8) :: str, ju, mu, jl, ml
	if (ju == jl) then
		if (mu == ml) then
			str = mu**2
		else if (mu > ml) then
			str = (ju+mu)*(ju-mu+1.d0)
		else
			str = (ju-mu)*(ju+mu+1.d0)
		endif
	else if (ju < jl) then
		if (mu == ml) then
			str = (ju+1.d0)**2 - mu**2
		else if (mu > ml) then
			str = (ju-mu+1.d0)*(ju-mu+2.d0)
		else
			str = (ju+mu+1.d0)*(ju+mu+2.d0)
		endif
	else
		if (mu == ml) then
			str = ju**2-mu**2
		else if (mu > ml) then
			str = (ju+mu)*(ju+mu-1.d0)
		else
			str = (ju-mu)*(ju-mu-1.d0)
		endif	
	endif
end function str

!--------------------------------------------------------------
! Calculates the absorption matrix, the emission vector
! and the fis 
!--------------------------------------------------------------
subroutine abem(bol,k,in,idir,ila,delw,ab,em,pem,bnu,fis)
real(kind=8) :: bol, delw, ab(4,4), em(4), pem(4), bnu, fis(4)
real(kind=8) :: gam, chi, cosg, cos2g, sin2g, cos2ch, sin2ch
real(kind=8) :: hp, fvp, hb, fvb, hr, fvr, temp
real(kind=8) :: etap, rhop, etab, rhob, etar, rhor
real(kind=8) :: etai, etaq, etau, etav, rhoq, rhou, rhov, mu
integer :: k, in, idir, ila

! Local opacity
	eta0 = etaz(k)
! Set opacity to 0 to compute the continuum
	if (delw > 1.d4) eta0 = 0.d0
	
! Transform to radians
	gam = gamma(k) * PI_180
	chi = ch(k) * PI_180
	
	bnu = planc(k)
	sline = Sli(k)
	vlos = vel(k)
	magf = bmag(k)
	dopp = dwid(k)
	a = aval(k)
	
! We are considering the magnetic field vertical, so cos(gamma)=mu
!	cosg = in * coz(idir)

! Try to include the effect of B direction
! cosg is the angle between the magnetic field and the line of sight (given by mu)
	mu = in*coz(idir)
	cosg = in*(mu * dcos(gam) + dsqrt(1.d0-mu**2)*dsin(gam))

	cos2g = cosg**2
	sin2g = 1.d0 - cos2g
	
! Azimuth of the magnetic field. It is vertical, so there is invariance
	cos2ch = dcos(2.d0*chi)!1.d0
	sin2ch = dsin(2.d0*chi)!0.d0
	
! Call generalized Voigt routine to get the profiles
	call vgtgen(a,delw,hp,fvp,hb,fvb,hr,fvr)

! These are the phi_i of Rees et al ApJ 339, 1093
	etap = eta0 * hp
	rhop = eta0 * fvp
	etab = eta0 * hb
	rhob = eta0 * fvb
	etar = eta0 * hr
	rhor = eta0 * fvr
	
! Now evaluate the matrix elements of the absorption matrix
! First the eta's
	etai = 0.5d0*etap*sin2g + 0.25d0*(etab+etar)*(1.d0+cos2g)
	etaq = (0.5d0*etap - 0.25*(etab+etar))*sin2g
		
	if (ditype == 'm') etaq = -etaq    ! Check for magnetic dipole transition
	etau = etaq*sin2ch
	etaq = etaq*cos2ch
	etav = 0.5d0*(etar-etab)*cosg
	
! Then the rho's
	rhoq = (0.5d0*rhop - 0.25*(rhob+rhor))*sin2g
		
	if (ditype == 'm') rhoq = -rhoq    ! Check for magnetic dipole transition
	rhou = rhoq*sin2ch
	rhoq = rhoq*cos2ch
	rhov = 0.5d0*(rhor-rhob)*cosg	
	
! Evaluate vector fis
! It seems to be the integral of the profiles
	fis(1) = 0.5d0*hp*wtp(ila,k)*sin2g + 0.25d0*(1.d0+cos2g)*(hb*wtb(ila,k)+hr*wtr(ila,k))
	temp = (0.5d0*hp*wtp(ila,k)-0.25d0*(hb*wtb(ila,k)+hr*wtr(ila,k)))*sin2g
	fis(2) = temp*cos2ch
	fis(3) = temp*sin2ch
	fis(4) = 0.5d0*(hr*wtr(ila,k)-hb*wtb(ila,k))*cosg
	
! Evaluate absorption matrix
! Perhaps this is a good place to renormalize frequency weigths if there are velocity fields
! and the magnetic field varies with height, imposing that int(fis d(nu))=1 (not sure)

	temp = etai + cont	
	ab(1,1) = temp
	ab(2,2) = temp
	ab(3,3) = temp
	ab(4,4) = temp
	ab(1,2) = etaq
	ab(2,1) = etaq
	ab(1,3) = etau
	ab(3,1) = etau
	ab(1,4) = etav
	ab(4,1) = etav
	ab(2,3) = rhov
	ab(3,2) = -rhov
	ab(2,4) = -rhou
	ab(4,2) = rhou
	ab(3,4) = rhoq
	ab(4,3) = -rhoq
	
! Evaluate emission vector
	em(1) = cont*bnu*bol + etai*sline
	em(2) = etaq*sline
	em(3) = etau*sline
	em(4) = etav*sline
	
! Evaluate pseudoemission vector (corresponding to sline=1 and bnu=0)
	pem(1) = etai
	pem(2) = etaq
	pem(3) = etau
	pem(4) = etav
	
end subroutine abem

! --------------------------------------------------------- &
! Formal solution of the polarized RT equation for the Zeeman case
! --------------------------------------------------------- &
subroutine iquvformal(sem, direction, boundary)
real(kind=8) :: sem(4,nw-1,ndir)
real(kind=8) :: abk(4,4), emk(4), pemk(4), bnu(nz), ex(nz), step(nz)
real(kind=8) :: rhs(4,nz), ab(4,4,nz), dz(nz), s(4)
real(kind=8) :: sto(4,nz), apro(4,nz), fis(4)
real(kind=8) :: u(4,nz), clft(4,4)
real(kind=8) :: ext(4,4), extn(4,4,nz)
real(kind=8) :: dtaus
integer :: iwave, idir, in, k, k0, k1, kdel, i, j, l
real(kind=8) :: dlw, dtu, dtd, w0, w1, a, s0, su, sd
real(kind=8) :: sform, extmp, dbdt, direction
 character*4 :: boundary

	coz(1) = direction
	boll = 1.d0
	in = 1


! WAVELENGTH LOOP
  do iwave = 2, nw

    dlw = dw(iwave)

    ! ANGLE LOOP
    do idir = 1, ndir

          k0 = nz
          k1 = 1
          kdel = -1

! First calculate, for this direction and wavelength, the absorption matrix, the emission vector
! and the fis for all spatial points. Keep this information
        do k = 1, nz
          call abem(boll,k,in,idir,iwave,dlw,abk,emk,pemk,bnu(k),fis)
			 
! Some kind of normalization of the vectors and matrix to etai
          do i = 1, 4
            rhs(i,k) = emk(i) / abk(i,i)
            apro(i,k) = pemk(i) / abk(i,i)
            sto(i,k) = fis(i)
            do j = 1, 4
              if (j == i) then
                ab(i,i,k) = abk(i,i)
              else
                ab(i,j,k) = abk(i,j) / abk(i,i)
              endif
            enddo
          enddo

          ! Calculate optical and physical grids
          if (k /= 1) then
				dtaus = 0.5d0*(z(k)-z(k-1))*(ab(1,1,k)+ab(1,1,k-1)) / direction
            dz(k-1) = dabs(dtaus)
            ex(k-1) = dexp(-dz(k-1))
            step(k-1) = z(k)-z(k-1)
          endif

        enddo !k

        call matinx(abk)

        dbdt = (bnu(nz)-bnu(nz-1)) / (z(nz)-z(nz-1))
        dbdt = boll * dbdt * direction
        if (boundary == 'DIFF') then
			  u(1,nz) = boll * bnu(nz) + abk(1,1) * dbdt
   	     u(2:4,nz) = abk(2:4,1) * dbdt
		  endif
		  if (boundary == 'ZERO') then
		  	  u(1,nz) = 0.d0
			  u(2:4,nz) = 0.d0
		  endif	

        ! Identity matrix
        extn(:,:,1) = 0.d0
        do i = 1, 4
          extn(i,i,1) = 1.d0
        enddo

        ! SPATIAL LOOP
        do k = k0+kdel, k1, kdel

          dtu = dz(k)
          if (k == k1) then
            dtd = dtu
          else
            dtd = dz(k+kdel)
          endif

          a = 1.d0 - ex(k)
          w0 = (dz(k) - a) / dz(k)
          w1 = (a - dz(k) * ex(k)) / dz(k)


          do i = 1, 4
            ! Source functions at points O and M
            s0 = rhs(i,k)
            su = rhs(i,k-kdel)
	
            if (k /= k1) then
              sd = rhs(i,k+kdel)
              sform = par_sc(trick,dtu,dtd,s0,su,sd)
              s(i) = sform 
            else
              sd = s0
              dtd = dtu
              sform = lin_sc(trick,dtu,s0,su)
              s(i) = sform 
            endif

            ext(i,1:4) = -w1 * ab(i,1:4,k-kdel)
            ext(i,i) = ex(k)

          enddo   !i

          clft = w0 * ab(:,:,k)
          do j = 1, 4
            clft(j,j) = 1.d0
          enddo

          call matinx(clft)

          do i = 1, 4
            do j = 1, 4
              extmp = 0.d0
              do l = 1, 4
                extmp = extmp + clft(i,l)*ext(l,j)
              enddo
              extn(i,j,k-kdel) = extmp
            enddo
          enddo


        do i = 1, 4
          u(i,k) = 0.d0
          do j = 1, 4
            u(i,k) = u(i,k) + clft(i,j)*s(j) + extn(i,j,k-kdel)*u(j,k-kdel)
          enddo
        enddo

       
      enddo !k
      sem(1:4,iwave-1,idir) = u(1:4,1)
   
  enddo !idir

enddo !iwave

end subroutine iquvformal

! ---------------------------------------------------------
! Solves RTE using short-characteristics. Source function is supposed to vary parabolically
! between points M, O and P and RTE is integrated analitically, resulting in
!            I(O)=I(M)*exp(-tau(MO))+psi(M)*S(M)+psi(O)*S(O)+psi(P)*S(P)
! where psi are functions of the optical distance tau(MO) and tau(OP)
! It returns the value of the three psi coefficients
! ---------------------------------------------------------
   function par_sc(tric,dtm,dtp,s0,sm,sp)
   real(kind=8) :: par_sc
   real(kind=8), INTENT(IN) :: dtm, dtp, s0, sm, sp, tric
   real(kind=8) :: exu, u0, u1 ,u2, d2, d3, d4, cm, c0, cp
         if (dtm.ge.tric) then
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
			  cm=(u2-u1*(dtp+2.d0*dtm))/(dtm*(dtm+dtp))+u0
      	  c0=(u1*(dtm+dtp)-u2)/(dtm*dtp)
      	  cp=(u2-dtm*u1)/(dtp*(dtm+dtp))
	 	  else
		  	  cm = 0.d0
			  c0 = 0.d0
			  cp = 0.d0
		  endif
		  
		  par_sc = cm*sm+c0*s0+cp*sp
		  
   end function par_sc
	
! ---------------------------------------------------------
! Solves RTE using short-characteristics. Source function is supposed to vary linearly
! between points M and O and RTE is integrated analitically, resulting in
!            I(O)=I(M)*exp(-tau(MO))+psi(M)*S(M)+psi(O)*S(O)
! where psi are functions of the optical distance tau(MO)
! It returns the value of the two psi coefficients
! ---------------------------------------------------------
   function lin_sc(tric, dtm, s0, sm)
   real(kind=8) :: lin_sc
   real(kind=8), INTENT(IN) :: dtm, tric, s0, sm
   real(kind=8) :: exu, u0, u1, c0, cm, d2

      if (dtm.ge.tric) then
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

		lin_sc = cm*sm+c0*s0

   end function lin_sc	

!--------------------------------------------------------------
! Inversion of a 4x4 matrix
!--------------------------------------------------------------
subroutine matinx(a)
real(kind=8) :: a(4,4)
real(kind=8) :: b(4,4), det, maxim, fabsmax
! First some tests of singularity
	b = dabs(a)
	maxim = maxval(b)
	fabsmax = 1.d0 / maxim
	if (maxim == 0.d0) then
		print *, 'Singularity in the inversion'
		stop
	endif
	
	a = a * fabsmax
	
      b(1,1) = a(2,2) * a(3,3) * a(4,4) + a(2,3) * a(3,4) * a(4,2)&
      	+ a(2,4) * a(3,2) * a(4,3) - a(2,2) * a(3,4) * a(4,3)&
	      - a(2,3) * a(3,2) * a(4,4) - a(2,4) * a(3,3) * a(4,2)
      b(2,1) = a(2,3) * a(3,1) * a(4,4) + a(2,4) * a(3,3) * a(4,1)&
	      + a(2,1) * a(3,4) * a(4,3) - a(2,3) * a(3,4) * a(4,1)&
   	   - a(2,4) * a(3,1) * a(4,3) - a(2,1) * a(3,3) * a(4,4)
      b(3,1) = a(2,4) * a(3,1) * a(4,2) + a(2,1) * a(3,2) * a(4,4)&
	      + a(2,2) * a(3,4) * a(4,1) - a(2,4) * a(3,2) * a(4,1)&
   	   - a(2,1) * a(3,4) * a(4,2) - a(2,2) * a(3,1) * a(4,4)
      b(4,1) = a(2,1) * a(3,3) * a(4,2) + a(2,2) * a(3,1) * a(4,3)&
	      + a(2,3) * a(3,2) * a(4,1) - a(2,1) * a(3,2) * a(4,3)&
   	   - a(2,2) * a(3,3) * a(4,1) - a(2,3) * a(3,1) * a(4,2)
      b(1,2) = a(3,2) * a(4,4) * a(1,3) + a(3,3) * a(4,2) * a(1,4)&
	      + a(3,4) * a(4,3) * a(1,2) - a(3,2) * a(4,3) * a(1,4)&
   	   - a(3,3) * a(4,4) * a(1,2) - a(3,4) * a(4,2) * a(1,3)
      b(2,2) = a(3,3) * a(4,4) * a(1,1) + a(3,4) * a(4,1) * a(1,3)&
	      + a(3,1) * a(4,3) * a(1,4) - a(3,3) * a(4,1) * a(1,4)&
   	   - a(3,4) * a(4,3) * a(1,1) - a(3,1) * a(4,4) * a(1,3)
      b(3,2) = a(3,4) * a(4,2) * a(1,1) + a(3,1) * a(4,4) * a(1,2)&
	      + a(3,2) * a(4,1) * a(1,4) - a(3,4) * a(4,1) * a(1,2)&
   	   - a(3,1) * a(4,2) * a(1,4) - a(3,2) * a(4,4) * a(1,1)
      b(4,2) = a(3,1) * a(4,2) * a(1,3) + a(3,2) * a(4,3) * a(1,1)&
	      + a(3,3) * a(4,1) * a(1,2) - a(3,1) * a(4,3) * a(1,2)&
   	   - a(3,2) * a(4,1) * a(1,3) - a(3,3) * a(4,2) * a(1,1)
      b(1,3) = a(4,2) * a(1,3) * a(2,4) + a(4,3) * a(1,4) * a(2,2)&
   	   + a(4,4) * a(1,2) * a(2,3) - a(4,2) * a(1,4) * a(2,3)&
	      - a(4,3) * a(1,2) * a(2,4) - a(4,4) * a(1,3) * a(2,2)
      b(2,3) = a(4,3) * a(1,1) * a(2,4) + a(4,4) * a(1,3) * a(2,1)&
	      + a(4,1) * a(1,4) * a(2,3) - a(4,3) * a(1,4) * a(2,1)&
   	   - a(4,4) * a(1,1) * a(2,3) - a(4,1) * a(1,3) * a(2,4)
      b(3,3) = a(4,4) * a(1,1) * a(2,2) + a(4,1) * a(1,2) * a(2,4)&
	      + a(4,2) * a(1,4) * a(2,1) - a(4,4) * a(1,2) * a(2,1)&
   	   - a(4,1) * a(1,4) * a(2,2) - a(4,2) * a(1,1) * a(2,4)
      b(4,3) = a(4,1) * a(1,3) * a(2,2) + a(4,2) * a(1,1) * a(2,3)&
	      + a(4,3) * a(1,2) * a(2,1) - a(4,1) * a(1,2) * a(2,3)&
   	   - a(4,2) * a(1,3) * a(2,1) - a(4,3) * a(1,1) * a(2,2)
      b(1,4) = a(1,2) * a(2,4) * a(3,3) + a(1,3) * a(2,2) * a(3,4)&
	      + a(1,4) * a(2,3) * a(3,2) - a(1,2) * a(2,3) * a(3,4)&
   	   - a(1,3) * a(2,4) * a(3,2) - a(1,4) * a(2,2) * a(3,3)
      b(2,4) = a(1,3) * a(2,4) * a(3,1) + a(1,4) * a(2,1) * a(3,3)&
	      + a(1,1) * a(2,3) * a(3,4) - a(1,3) * a(2,1) * a(3,4)&
   	   - a(1,4) * a(2,3) * a(3,1) - a(1,1) * a(2,4) * a(3,3)
      b(3,4) = a(1,4) * a(2,2) * a(3,1) + a(1,1) * a(2,4) * a(3,2)&
	      + a(1,2) * a(2,1) * a(3,4) - a(1,4) * a(2,1) * a(3,2)&
   	   - a(1,1) * a(2,2) * a(3,4) - a(1,2) * a(2,4) * a(3,1)
      b(4,4) = a(1,1) * a(2,2) * a(3,3) + a(1,2) * a(2,3) * a(3,1)&
	      + a(1,3) * a(2,1) * a(3,2) - a(1,1) * a(2,3) * a(3,2)&
   	   - a(1,2) * a(2,1) * a(3,3) - a(1,3) * a(2,2) * a(3,1)
			
		det = a(1,1) * b(1,1) + a(1,2) * b(2,1) + a(1,3) * b(3,1) + a(1,4) * b(4,1)			
		
		a = b * (fabsmax / det)

end subroutine matinx

!--------------------------------------------------------------
! Deallocate allocated memory
!--------------------------------------------------------------
subroutine cleanall_zeeman
	if (allocated(dw)) deallocate(dw)
	if (allocated(z)) deallocate(z)
	if (allocated(etaz)) deallocate(etaz)
	if (allocated(aval)) deallocate(aval)
	if (allocated(dwid)) deallocate(dwid)
	if (allocated(Sli)) deallocate(Sli)
	if (allocated(planc)) deallocate(planc)
	if (allocated(bmag)) deallocate(bmag)
	if (allocated(gamma)) deallocate(gamma)
	if (allocated(ch)) deallocate(ch)
	if (allocated(vel)) deallocate(vel)
	if (allocated(eps)) deallocate(eps)	
	if (allocated(ytr)) deallocate(ytr)
	if (allocated(ytp)) deallocate(ytp)
	if (allocated(ytb)) deallocate(ytb)
	if (allocated(wtr)) deallocate(wtr)
	if (allocated(wtp)) deallocate(wtp)
	if (allocated(wtb)) deallocate(wtb)	
	if (allocated(cox)) deallocate(cox)
	if (allocated(coz)) deallocate(coz)
	if (allocated(wtdir)) deallocate(wtdir)
	if (allocated(bnu)) deallocate(bnu)
	if (allocated(phij)) deallocate(phij)
	if (allocated(emerg)) deallocate(emerg)
	if (allocated(bnu)) deallocate(bnu)
	if (allocated(phi1)) deallocate(phi1)
	if (allocated(phi2)) deallocate(phi2)
	if (allocated(phi3)) deallocate(phi3)
	if (allocated(phi4)) deallocate(phi4)
end subroutine cleanall_zeeman

end module functions_delopar



module zeeman
use variables
use general
use variables_zeeman
use functions_delopar
implicit none

contains

! *********************************************************
! *********************************************************
! ZEEMAN SYNTHESIS ROUTINES
! *********************************************************
! *********************************************************

! ---------------------------------------------------------
! Read the Zeeman config file
! ---------------------------------------------------------	
	subroutine zeeman_synth
	integer :: i, j, k
	
		open(unit=12,file='zeeman.dat',action='read',status='old')
		call lb(12,2)
		
! Read the line in which the Zeeman synthesis is going to be performed
		call lb(12,2)
		read(12,*) zeeman_line				
		
! Read the number of wavelength points
		call lb(12,2)
		read(12,*) nw
		
! Read the step in wavelength
		call lb(12,2)
		read(12,*) wstep
		
! Read the direction in which one is observing
		call lb(12,2)
		read(12,*) mu_read		
		mu_read = dcos(mu_read * PI / 180.d0)
		mu_direction = mu_read
		
! Read the field strength
		call lb(12,2)
		read(12,*) field_strength
		
! Read the field inclination
		call lb(12,2)
		read(12,*) inclination

! Read the field azimuth
		call lb(12,2)
		read(12,*) azimuth
								
! Read the type of transition
		call lb(12,2)
		read(12,*) type_transition
		
! Read the type of Zeeman effect
		call lb(12,2)
		read(12,*) type_zeeman		
		
! Read the type of transition
		call lb(12,2)
		read(12,*) magneto_optics
		
! Read the name of the output file
		call lb(12,2)
		read(12,*) output_file
								
		close(12)
		
! Only synthesize one line
		if (zeeman_line /= -1) then
			call initialize_vars

			allocate(emerg(4,nw-1,ndir))
	   	boll = 1.d0

			print *, 'Solving the RT equation for the polarized case with DELOPAR'
	   	call iquvformal(emerg,mu_direction,'DIFF')

			open(unit=2, file=output_file,status='replace',action='write')
	   	write(2,*) 1
			write(2,*) wave
			write(2,*) nw-1
			do i = 1, nw-1
				write(2,*) (dw(i+1),emerg(:,i,j),j=1,ndir)
			enddo
			close(2)
			call cleanall_zeeman
			
! Synthesize every line in the model
		else
			open(unit=2, file=output_file,status='replace',action='write')
			write(2,*) nr
			print *, 'Solving the RT equation for the polarized case with DELOPAR'
			do k = 1, nr
				zeeman_line = k
				print *, 'Line number ', k
				call initialize_vars
				print *

				allocate(emerg(4,nw-1,ndir))
	   		boll = 1.d0
				
	   		call iquvformal(emerg,mu_direction,'DIFF')

	   		write(2,*) wave
				write(2,*) nw-1
				do i = 1, nw-1
					write(2,*) (dw(i+1),emerg(:,i,j),j=1,ndir)
				enddo	
				
				call cleanall_zeeman
				
			enddo
			
			close(2)
		endif

			
	end subroutine zeeman_synth
		

end module zeeman

