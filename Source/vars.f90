module constantes
implicit none

! *********************************************************
! *********************************************************
! GLOBAL VARIABLES
! *********************************************************
! *********************************************************

! ---------------------------------------------------------
! PHYSICAL CONSTANTS AND RELATIONS AMONG THEM (ALL IN CGS)
! Taken from NIST (data from 1998)
!  http://physics.nist.gov/cuu/Constants/index.html
! ---------------------------------------------------------
! PI : pi
! PME : electron mass
! PK : Boltmann's constant
! PH : Planck's constant
! PC : speed of light
! UMA : atomic mass unit
! A0 : Bohr radius
! PION : hydrogen ionization potential
! ONETHIRD : 1/3
! FOURSQRTTWO : 1/(4*sqrt(2))
! CI : 1/2 * (h^2/(2*pi*me*k))^(3/2)  Constant for the Saha equation
! PI8C2 : 8 * pi / c^2
! PI4H : 4 * pi / h
! PHK : h / k
! PKUMA : k / amu
! PHC : 2 * h / c^2
! FOURSQRTTWO : 1/(4*sqrt(2))
! SQRTPI : sqrt(PI)
! FWHM : sqrt(8*ln(2))
! ---------------------------------------------------------
	real(kind=8), parameter :: PI = 3.14159265359d0, PME = 9.10938188d-28, PK = 1.3806503d-16
	real(kind=8), parameter :: PH = 6.62606876d-27, PC = 2.99792458d10, UMA = 1.66053873d-24, CI = 2.0706d-16
	real(kind=8), parameter :: A0 = 5.291772083d-9, PION = 2.17987190d-11, PEC = 4.803d-10
	real(kind=8), parameter :: PI8C2 = 8.d0 * PI / (PC**2), PI4H = 4.d0 * PI / PH, PHK = PH / PK
	real(kind=8), parameter :: PKUMA = PK / UMA, ONETHIRD = 1.d0 / 3.d0
	real(kind=8), parameter :: PHC = 2.d0 * PH / PC**2, FOURSQRTTWO = 0.176776695297d0
	real(kind=8), parameter :: SQRTPI = 1.77245385091d0, PARSEC = 3.085678d18
	real(kind=8), parameter :: FWHM = 2.35482004503d0
	real(kind=8), parameter :: C_graf = -25.16d0, C_sil = -25.11d0

! ---------------------------------------------------------
! CONSTANTS FOR LIMITS
! ---------------------------------------------------------
! nordm : maximum order of Ng acceleration
! n_ptos : maximum number of shells
! ---------------------------------------------------------

	integer, parameter :: nordm = 4
   integer, parameter :: n_ptos = 801	

end module constantes

module variables
use constantes

! ---------------------------------------------------------
! GENERAL GLOBAL VARIABLES
! ---------------------------------------------------------
!       REAL*8
! tauip1 : it is used to indicate the minimum optical depth to be considered as thin and so not
!        make formal solution
! distancia : distance to the source (pc)
! ang_size : angular size of the source (arcsec)
! tel_fwhm : full width at half maximum of the telescope (FWHM)
! omega_sor : omega parameter in SOR iterative scheme
! star_radius: radius of the central star if present 
! gaussq_n : number of points in the gaussian quadrature used in LVG 
!
!       INTEGER
! n_radios : number of shells
! nfrq : number of frequency points in each profile
! n_freq_ud_doppler : number of frequency points in each Doppler width of the profile
! n_freq_perfil : number of Doppler width which cover the profile
! nord : Ng acceleration order
! itracc : number of iterations to accelerate
! nintracc : number of iterations after accelerating to re-accelerate
! iacel : type of acceleration ( 1 -> Ng, otherwise -> no acceleration)
! iter : current iterative step
! cortes : number of intersections between a shell and the characteristics
! n_caract : total number of characteristics used
! caract_core : number of characteristics crossing the core
! iter_max : maximum number of iterations
! output_spectrum : flag to write the emerging spectrum or not
! impact_star_radius : number of impact parameters which cross the central star if present
!
!       CHARACTER*30
! nombre_salida : output file
! nombre_iteration : iteration file
! nombre_modelo : atomic/molecular model file
! nombre_atmosf : atmosphere file
! nombre_emerging : emerging spectrum output file
! nombre_atmosf_out : atmosphere output file
! nombre_flux : emerging flux output file
! nombre_background_out : background opacity output file
! nombre_radrates : radiative rates output file
! tipo_iteracion : type of iterative scheme to use
! cortes : n. cortes entre una caracteristica y una capa! 
! hor_axis : horizontal axis of the ellipse (vertical axis is equivalent to radius)
! geometry_type : type of geometry we are using
! ---------------------------------------------------------

	real(kind=8) :: tauip1, distancia, ang_size, omega_sor, maj_axis, tel_fwhm, tel_fwhm_cut, star_radius
	integer :: n_radios, nfrq, n_freq_ud_doppler, n_freq_perfil, nord, itracc, nintracc, iacel, iter
	integer :: cortes, n_caract, caract_core, iter_max, output_spectrum, impact_star_radius
	integer :: angleset_pp, gaussq_n
	character*60 :: nombre_salida, nombre_modelo, nombre_atmosf, nombre_iteration
	integer :: tipo_iteracion
	character*60 :: nombre_emerging, nombre_atmosf_out, nombre_flux, geometry_type
	character*60 :: nombre_background_out, nombre_radrates, nombre_collisrates
   character*10 :: observing_instrument, central_source_type

! ---------------------------------------------------------
! DEPTH VARYING VARIABLES
! ---------------------------------------------------------
! tau(transition,shell) : optical depth of each transition at each shell
! B(shell) : continuum source function (Planck's function at kinetic temperature)
! chil(shell) : line opacity
! kappa(shell) : continuum opacity
! chie(shell) : electron scattering opacity
! Sl(shell) : line source function
! Lstar(shell) : Lambda operator's diagonal
! Jbar(shell) : frequency and angle averaged intensity
! Jbar20(transition, shell) : frequency and angle averaged with (3*mu^2-1) intensity (for polarization)
! p(shell) : position of each characteristic
! r(shell) : shell grid
! B_dust(shell) : dust source function (Planck's function at dust temperature)
! kappa_dust(shell) : dust opacity
! ---------------------------------------------------------

	real(kind=8), allocatable :: tau(:,:), B(:), chil(:), chic(:,:), Sc(:,:), chie(:), Sl(:), Lstar(:)
	real(kind=8), allocatable :: Jbar(:), p(:), r(:), B_dust(:), kappa_dust(:), Jbar20(:,:), chic_bf(:,:,:), Sc_bf(:,:,:)
	real(kind=8), allocatable :: kappa(:), wtmu(:), mu_pp(:)
	real(kind=8), allocatable :: Jbarb(:), Lstarb(:), Lstarc(:), Ibf(:)

! ---------------------------------------------------------
! FREQUENCY VARYING VARIABLES
! ---------------------------------------------------------
! perfil(frq) : line profile
! profile(frq,caract,shell,in/out) : line profile in each point
! frqwt(frq,caract,shell,in/out) : frequency integration weights
! overl_lim(shell,caract,inout,transition_i,transition_j,begin/end) : indicates the begining
!		and the end of the overlapping between transition i and j
! In(charact,frq) : specific intensity for each characteristic and frequency
! wtnu(frq) : frequency integration weights
! eps_graphite_par(2,lambda) : graphite complex parallel dielectric function
! eps_graphite_per(2,lambda) : graphite complex perpendicular dielectric function
! eps_silicate(2,lambda) : silicate complex dielectric function
! graphite_lambda(lambda) : wavelengths where graphite's dielectric function is tabulated
! silicate_lambda(lambda) : wavelengths where silicate's dielectric function is tabulated
! ---------------------------------------------------------

	real(kind=8), allocatable :: perfil(:), In(:,:), wtnu(:), profile(:,:,:,:), frqwt(:,:,:,:)
	real(kind=8), allocatable :: eps_graphite_par(:,:), eps_graphite_per(:,:), eps_silicate(:,:)
	real(kind=8), allocatable :: graphite_lambda(:), silicate_lambda(:)
	integer, allocatable :: overl(:,:,:)


! ---------------------------------------------------------
! ANGULAR VARYING VARIABLES
! ---------------------------------------------------------
! mu(shell,charact) : angle between the shell and the characteristic for each intersection
! ---------------------------------------------------------

	real(kind=8), allocatable :: mu(:,:)
	
! ---------------------------------------------------------
! GLOBAL JBAR AND LSTAR MATRICES
! ---------------------------------------------------------
! Jbar_total(shell*level) : total Jbar matrix for all levels and shells
! Lstar_total(shell*level) : total Lambda operator diagonal for all levels and shells
! ---------------------------------------------------------	
		
   real(kind=8), allocatable :: Jbar_total(:), lstar_total(:)
	real(kind=8), allocatable :: Jbar_total_bf_a(:), Jbar_total_bf_b(:), I_total_bf(:)
	real(kind=8), allocatable :: lstar_total_bf_a(:), lstar_total_bf_b(:), lstar_total_bf_c(:)

! ---------------------------------------------------------
! OTHER VARIABLES
! ---------------------------------------------------------
! scratch(level*shells,nord) : Ng temporal storage
! precision : final maximum relative change to reach
! relat_error : current iteration maximum relative change
! relat_error_p : last iteration maximum relative change
! trick : threshold value to turn on exponential series expansion
! repeat : variable to say if we have an empty core or an absorbing/source one
! ---------------------------------------------------------	

   real(kind=8), allocatable ::  scratch(:,:)
	real(kind=8) :: precision, relat_error, trick, relat_error_p, repeat, true_error

	
! ---------------------------------------------------------
! ATOMIC MODEL
! ---------------------------------------------------------
!       INTEGER
! nl : number of levels
! ni : number of ion
! nt : number of transitions
! nr : number of radiative transitions
! nact : number of active transitions
! nli(nlm) : ion to which the level belongs to
! itran(2,ntm) : upper and lower levels of each transition
! n_partition : number of coefficients in the partition function fit
!
!       REAL*8   
! dlevel(1,nlm) : level energy in frequency units
! dlevel(2,nlm) : statistical weight of each level
! dion(1,ni) : ionization potential in frequency units
! dtran(1,nt) : effective cross-section of the transition = fij * pi * e^2 / (me * c)
!					It is converted to h*nu/(4*pi)*Blu
! dtran(2,nt) : transition frequency (Hz)
! dtran(3,nt) : collisional cross-section (only used for reading from file, see collision(:,:))
! dtran(4,nt) : Doppler width or radiation temperature
! doppler_width(shell,transition) : doppler width for each transition and point in the atmosphere
! collision(shell,transition) : collisional rates for each transition and shell
! partition(n_partition) : coefficients of the partition function
! mass : mass of the atom or molecule
!
!       CHARACTER*20
! tipo_colision : type of collisions (FIXED, FIXED_H, INTERPOL_H, etc)
! tipo_modelo : indicates if the model is an atom or a molecule
! atom_mol_nombre : indicates the model name
! ---------------------------------------------------------

	integer :: nl, ni, nact, nt, nr, n_partition, nbf, nfrq_bf, include_bound_free
	integer, allocatable :: nli(:), itran(:,:)
	real(kind=8), allocatable :: dlevel(:,:), dion(:,:), dtran(:,:), doppler_width(:,:)
	real(kind=8), allocatable :: partition(:), collision(:,:), threshold_bf(:), bf_cross_section(:,:,:)
	character*40, allocatable :: label(:)
	real(kind=8) :: mass
	character*20 :: tipo_colision, tipo_modelo, atom_mol_nombre
	
! ---------------------------------------------------------
! ATMOSPHERE MODEL
! ---------------------------------------------------------	
!       REAL*8   
! datph(1,shell) : kinetic temperature in each shell
! datph(2,shell) : electron density
! datph(3,shell) : atom/molecule abundance
! datph(4,shell) : microturbulence velocity
! datph(5,shell) : visual absorption A_v
! datph(6,shell) : dust temperature
! datph(7,shell) : macroscopic velocity
! datph(8,shell) : gas pressure
! dentro(nfrq) : inner shell boundary condition
! fuera(nfrq) : outer shell boundary condition
! tback : microwave background temperature 
! tcentral : central source temperature
! lambda_ref : dust reference wavelength
! expav : dust opacity exponent chi_dust = A_v * (lambda/lambda_ref)**(-expav)
! xhel : helium abundance referring to H2
! xhid : atomic hydrogen abundance referring to H2
!       INTEGER 
! background_flag : include background radiation or not
! central_flag : include central source or not
! dust_flag : include dust or not
! back_opac_flag : include or not background opacity
! radrates_flag : write or not radiative rates
! iterat_improv : do iterative improvement of the solution to the linearized system
! V_max : maximum velocity in the atmosphere
! lvg_flag : maximum relative change of the LVG solution (if 0, then use LTE)
! SNTB_flag : maximum relative change to turn on Socas Navarro-Trujillo Bueno acceleration
!       If 0, then don't use it
! ---------------------------------------------------------
	real(kind=8), allocatable :: datph(:,:), dentro(:), fuera(:), thermal_v(:), dentro_bf(:), fuera_bf(:)
	real(kind=8) :: tback, tcentral, lambda_ref, expav, xhel, xhid, V_max, lvg_flag, SNTB_flag
	real(kind=8) :: ref_wavelength, ref_opacity, spectral_index
	integer :: background_flag, dust_flag, central_flag, iterat_improv, back_opac_flag
	integer :: radrates_flag, converge_flag, overlap_flag, error_flag
	integer :: zeeman_line, include_collis, linear_solver_algorithm, verbose_flag
	integer :: level_highest_change, interpolate_atmosphere, n_radios_old, n_caract_old
	
! ---------------------------------------------------------
! POPULATION
! ---------------------------------------------------------
! nh(shell) : atomic hydrogen (atom model) or molecular hydrogen (molecular model) density 
! popl(shell*nl) : LTE population in each shell for all levels
! pop(shell*nl) : population of all levels in each shell
! pop(shell*nl) : initial population (LVG or LTE)
! ---------------------------------------------------------

	real(kind=8), allocatable :: nh(:), popl(:), pop(:), popi(:), pop_final(:)
	
! ---------------------------------------------------------
! CEP method
! ---------------------------------------------------------
! ---------------------------------------------------------	
	integer :: n_quadr_beta	
	real(kind=8), allocatable :: dSldn(:,:), dJbar_totaldn(:,:), dJbardn(:,:), x_e3(:), w_e3(:)
	
end module variables
