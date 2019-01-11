Program Aerosol_box

  IMPLICIT NONE

  !Definition of variables
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15,300)

  REAL(dp), PARAMETER :: pi = 2D0*ASIN(1D0)
  REAL(dp), PARAMETER :: g = 9.81D0                     ! gravitation const. m s^-2
  INTEGER, PARAMETER ::  nr_bins = 100                  ! Number of particle size bins
  INTEGER, PARAMETER ::  nr_cond = 2                    ! Number of condensable vapours
  REAL(DP), PARAMETER :: Rg=8.3145D0                    ! Universal gas constant J mol^-1 K^-1
  REAL(DP), PARAMETER :: Na=6.022D23                    ! Avogadro's number 
  REAL(DP), PARAMETER :: Mair=28.96D-3                  ! Mean molecular weight of air
  REAL(DP), PARAMETER :: kb = 1.381d-23                 ! Boltzmann constant [(m2*kg)/(s2*K)]
    
  REAL(DP), DIMENSION(nr_bins) :: diameter  ,&          ! Diameter of each size bin
       particle_mass                        ,&          ! mass of one particle in each size bin
       particle_conc                        ,&          ! number concentration in each size bin
       particle_volume                      ,&          ! volume concentration in each size bin 
       coag_loss                                        ! coagulation loss rate of particles in each size bin
       
  REAL(DP), DIMENSION(nr_cond) :: molecular_mass   ,&   ! molecular mass of the condensing vapours [kg/#]
       molecular_volume                            ,&   ! Molecule volume of condensable vapours [m^3]
       molecular_dia                               ,&   ! Molecule diameter of condensable vapours [m]
       molar_mass                                  ,&   ! Molar mass of condensable vapours [kg/m^3]
       cond_vapour                                      ! Concentration of condensable vapours [molec/m^3] 
                                        
  REAL(DP), DIMENSION(nr_cond) :: Cond_sink = 0.001     ! Assumed initial condensation sink of vapours [s^-1]

  REAL(DP) ::         PN, PM, PV                 ! Total particle number and mass concentration [cm^-3] 

       
  REAL(DP)    particle_density              ,&    ! [Kg]
              nucleation_coef               ,&    ! Nucleation coefficient 
              mass_accomm                         ! mass accomodation coefficient 

  REAL(dp) :: nucleation_rate                     ! #/(m3*s)

  REAL(DP) :: time                          ,&    ! time in simulation [s]
              simu_hours                    ,&    ! total simulation time in hours
              timestep                      ,&    ! Model time step [s]
              time_start, time_end          ,&    ! Start time and end time of model [s]
              temp                          ,&    ! Temperature [K]
              pres                                ! Pressure [Pa]
  
  integer :: i
 
  !! ======================= Programe starts ===========================================================================

  ! Assign values to parameters and initialize the simulation
  !IsNucleation = .true.
  !Iscondenstion = .false.
  !Iscoagulation = .false.


  CALL Aerosol_init(diameter, particle_mass, particle_volume, particle_conc, &
           particle_density, nucleation_coef, molecular_mass, molar_mass, &
           molecular_volume, molecular_dia, mass_accomm)
       
  PN=SUM(particle_conc)*1D-6                 ! Total particle number concentration (cm^-3)     #/cm^3
  PM=SUM(particle_conc*particle_mass)    ! Total particle mass concentration (ug/m^3)
  PV=SUM(particle_conc * particle_volume) *1D12

  simu_hours = 24D0 ! hours
  timestep = 1D1 ! s
  time_start = 0D0
  time_end = time_start + simu_hours*3600.
  time = time_start

CALL Open_files ! Open output files 


DO WHILE (time .lt. time_end) ! Main program time step loop
! Meteorological parameters:
  temp = 300D0  ! K
  pres = 1D5    ! Pa
  
  cond_vapour(1) = 1D13*sin(pi/86400 * time) ! H2SO4 molec / m^3   
  cond_vapour(2) = 1D13                      ! ELVOC molec / m^3
  
  !!! Calculate new particle formation (nucleation) here !!!
  call calculateNucleo(cond_vapour(1),1.0_dp,temp, nucleation_rate)
  particle_conc(1) = particle_conc(1) + (nucleation_rate * timestep)
  PN = SUM(particle_conc)*1D-6
  PV = SUM(particle_conc * particle_volume) *1D12
  
  
  !!! Calculate coagulation losses here !!!

  call Coagulation(timestep, particle_conc, diameter,temp,pres,particle_mass)
  PN = SUM(particle_conc)*1D-6
  PV = SUM(particle_conc * particle_volume) *1D12
  

  
  !!! Calculate condensation particle growth here !!!

  call Condensation(timestep, temp, pres, mass_accomm, molecular_mass, &
    molecular_volume, molar_mass, molecular_dia, particle_mass, particle_volume, &
    particle_conc, diameter, cond_vapour)
  PN = SUM(particle_conc)*1D-6
  PV = SUM(particle_conc * particle_volume) *1D12
    
  
  time = time + timestep

! Save output variables every 10 minutes: 
IF (MOD(INT(time-timestep), 600) .eq. 0) THEN 
    CALL Write_files
    write(*,*) 'time', time/3600.
  ENDIF
     
END DO

CONTAINS

SUBROUTINE Aerosol_init(diameter, particle_mass, particle_volume, particle_conc, &
  particle_density, nucleation_coef, molecular_mass, molar_mass, &
  molecular_volume, molecular_dia, mass_accomm)

    !! ====================== Definition of variables =====================================================================

    REAL(DP), DIMENSION(nr_bins), INTENT(OUT) :: diameter       , &    ! diamter of each size bin
                                                 particle_mass  , &    ! mass of one particle
                                                 particle_volume  , &  ! volume of one particle
                                                 particle_conc         ! number concentration

    REAL(DP), DIMENSION(nr_cond), INTENT(OUT) :: molecular_mass ,& ! molecular mass of the condensing vapours [kg/#]
                                           molecular_volume ,&     ! [m3]
                                           molecular_dia, &        ! [m]
                                           molar_mass              ! molar mass of the condensing vapours [kg/mol]

    REAL(DP), INTENT(OUT) :: nucleation_coef, mass_accomm

    REAL(DP), DIMENSION(nr_cond) :: density             ! Bulk density of condensing vapours [kg/m^3]
 
    REAL(DP)  ::  particle_density
    
    INTEGER :: i
    mass_accomm = 1D0   ! Mass accommodation coefficient 
    
    nucleation_coef = 1D-20
    
    ! Particle diameters between 2D-9 and 2.5D-6 m:
    diameter(1)=2D-9 
    DO i=2,nr_bins
    diameter(i)=diameter(i-1)*(2.5D-6/diameter(1))**(1D0/(nr_bins-1))
    END DO
      
    particle_conc = 1D0 ! Assume an initial particle number concentration of 1 m^-3
    where((abs(diameter-2D-7)-MINVAL(abs(diameter-2D-7)))<1D-20)  particle_conc=2D8 ! add 200 cm^-3 200 nm sized accumulation mode particles
    
    particle_density = 1.4D3                                        ! Assumed fixed particle density [kg/m^3]
    particle_volume = 1D0/6D0 * pi * diameter**3                      ! Single particle volume (m^3)
    particle_mass=  1D0/6D0 * pi * diameter**3 * particle_density     ! [kg]

    density = (/1.84D3, 1.4D3/)                                     ! density of sulphuric acid and SOA
    molar_mass = (/0.098D0, 0.3D0/)                                 ! H2SO4 and ELVOC
    molecular_mass = molar_mass / Na                                ! molecular mass [kg]
    molecular_volume = molecular_mass / density                     ! molecular volume [m^3]
    molecular_dia = (6D0 * molecular_volume / pi )**(1D0/3D0)       ! molecular diameter [m]

END SUBROUTINE Aerosol_init

SUBROUTINE Condensation(timestep, temp, pres, mass_accomm, molecular_mass, &
  molecular_volume, molar_mass, molecular_dia, particle_mass, particle_volume, &
  particle_conc, diameter, cond_vapour) ! Add more variables if you need it

    !REAL(DP), DIMENSION(nr_cond) :: density
  
    REAL(DP), DIMENSION(nr_bins), INTENT(IN) :: diameter, particle_mass
    REAL(DP), DIMENSION(nr_cond), INTENT(IN) :: molecular_mass, molecular_dia, &
    molecular_volume, molar_mass
    REAL(DP), INTENT(IN) :: timestep, temp, pres, mass_accomm

    REAL(DP), DIMENSION(nr_bins), INTENT(INOUT) :: particle_conc

    REAL(DP), DIMENSION(2), INTENT(IN) :: cond_vapour  ! condensing vapour concentrations, which is H2SO4 and organics (ELVOC) [#/m^3]
    
    REAL(DP), DIMENSION(nr_bins), INTENT(IN)   :: particle_volume
    REAL(DP), DIMENSION(nr_bins) :: particle_volume_new
    
    REAL(DP), DIMENSION(nr_bins)   ::  slip_correction, diffusivity, speed_p, &
    particle_conc_new, x1, x2
    
    REAL(DP), DIMENSION(nr_cond)   ::  diffusivity_gas, speed_gas
    
    REAL(DP) :: dyn_visc, l_gas, dens_air
    
    INTEGER :: j,i
    REAL(DP), PARAMETER :: K_Boltzman = 1.38E-23

    REAL(DP), DIMENSION(nr_cond) ::  &
                              nn_i, &
                              CR, &
                              Eq, &
                              C_eq 
    REAL(DP), PARAMETER :: alpha = 1.0_dp
    Real(dp), dimension(nr_bins) :: nn_p,lambda_1, lambda_2, Kn_1, Kn_2, FuchsSutugin_1, FuchsSutugin_2, CR_1, CR_2, Eq_1, Eq_2
    
    ! Add more variabels if you need it...

    dyn_visc = 1.8D-5*(temp/298D0)**0.85D0  ! dynamic viscosity of air
    dens_air=Mair*pres/(Rg*temp)        ! Air density
    l_gas=2D0*dyn_visc/(pres*SQRT(8D0*Mair/(pi*Rg*temp))) ! Gas mean free path in air (m)

    slip_correction = 1D0+(2D0*l_gas/(diameter))*&
    (1.257D0+0.4D0*exp(-1.1D0/(2D0*l_gas/diameter))) ! Cunninghams slip correction factor (Seinfeld and Pandis eq 9.34) 
    
    diffusivity = slip_correction*kb*temp/(3D0*pi*dyn_visc*diameter)   ! Diffusivity for the different particle sizes m^2/s
    speed_p = SQRT(8D0*kb*temp/(pi*particle_mass))                     ! speed of particles (m/s)
  
    diffusivity_gas=5D0/(16D0*Na*molecular_dia**2D0*dens_air)*&
    SQRT(Rg*temp*Mair/(2D0*pi)*(molar_mass+Mair)/molar_mass)            ! Diffusivity of condensable vapours (m^2 s^-1)

    ! Thermal velocity of vapour molecule
    speed_gas=SQRT(8D0*kb*temp/(pi*molecular_mass)) ! speed of H2SO4 and ELVOC molecule

    !Calculate nn_i and nn_p
    nn_i(1) = speed_gas(1)
    nn_i(2) = speed_gas(2)
    nn_p = speed_p
    !nn_p(2) = ( 8.0 * K_Boltzman * temp / (PI * particle_mass(2)) )

    !Calculate lambda
    lambda_1 = 3.0 * (diffusivity_gas(1) + diffusivity) / sqrt(nn_i(1)**2 + nn_p**2)
    lambda_2 = 3.0 * (diffusivity_gas(2) + diffusivity) / sqrt(nn_i(2)**2 + nn_p**2)
    !Calculate Kn(1) and Kn(2)
    Kn_1 = 2.0 * lambda_1 / (molecular_dia(1) + diameter)
    Kn_2 = 2.0 * lambda_2 / (molecular_dia(2) + diameter)

    !Fuchs-Sutugin correction factor for H2SO4
    FuchsSutugin_1 = ( 0.75_dp * alpha*(1.0+Kn_1) / (Kn_1**2+Kn_1 + 0.283*Kn_1*alpha + 0.75_dp*alpha) )
    FuchsSutugin_2 = ( 0.75_dp * alpha*(1.0+Kn_2) / (Kn_2**2+Kn_2 + 0.283*Kn_2*alpha + 0.75_dp*alpha) )
    
    ! Calculate the Collision rate (CR [m^3/2]) between gas molecules (H2SO4 and ELVOC) and the particles:
    CR_1 = 2.0 * PI * (molecular_dia(1) + diameter) * (diffusivity_gas(1) + diffusivity) * FuchsSutugin_1
    CR_2 = 2.0 * PI * (molecular_dia(2) + diameter) * (diffusivity_gas(2) + diffusivity) * FuchsSutugin_2

    !Calculate Kelvin effect  sigma
    !C_eq(1) = pres * exp(4.0 * 1.0 * molar_mass(1) / (Rg*temp*density(1) * 1.0))
    !C_eq(2) = pres * exp(4.0 * 1.0 * molar_mass(2) / (Rg*temp*density(2) * 1.0))
    C_eq(1) = 0.0
    C_eq(2) = 0.0

    !Calculate equilibrium
    Eq_1 = CR_1 * molecular_volume(1) * (cond_vapour(1) - C_eq(1))
    Eq_2 = CR_2 * molecular_volume(2) * (cond_vapour(2) - C_eq(2))
    
    ! Calculate the new single particle volume after condensation (particle_volume_new):
    particle_volume_new = 0.0_dp

    particle_volume_new = particle_volume + CR_1 * cond_vapour(1)*molecular_volume(1) * (timestep) !	change	of	vp	due	to	condensation	of	H2SO4
    particle_volume_new = particle_volume_new + CR_2 * cond_vapour(2)*molecular_volume(2) * (timestep) !	change	of	vp	due	to	condensation	of	ELVOC
    
    ! Use the full-stationary method to divide the particles between the existing size bins (fixed diameter grid):
    particle_conc_new=0D0 ! Initialise a new vector with the new particle concentrations
    
    particle_conc_new(nr_bins)=particle_conc(nr_bins)
    !particle_conc_new = particle_conc

    DO j = 1,nr_bins-1
      x1(j) = (particle_volume(j+1) - particle_volume_new(j)) / (particle_volume(j+1) - particle_volume(j))
      x2(j) = 1.0 - x1(j)
      particle_conc_new(j) = particle_conc_new(j) + (x1(j)*particle_conc(j))
      particle_conc_new(j+1) = particle_conc_new(j+1) + (x2(j)*particle_conc(j))
    END DO

    particle_conc=particle_conc_new
    
END SUBROUTINE Condensation

SUBROUTINE Coagulation(timestep, particle_conc, diameter, &
  temp,pres,particle_mass) ! Add more variables if you need it
  
    REAL(DP), DIMENSION(nr_bins), INTENT(IN) :: diameter            !diameter of single particla
    REAL(DP), DIMENSION(nr_bins), INTENT(INOUT) :: particle_conc    !concentration of particles for bins
    REAL(DP), INTENT(IN) :: timestep
    REAL(DP), DIMENSION(nr_bins), INTENT(IN) :: particle_mass       ! mass of one particle for bins                          
    REAL(DP), INTENT(IN) :: temp, pres
    
    REAL(DP), DIMENSION(nr_bins,nr_bins) :: coagulation_coef        ! coagulation coefficients [m^3/s]

    REAL(DP), DIMENSION(nr_bins) :: slip_correction, diffusivity, dist, speed_p, &
    Beta_Fuchs, free_path_p

    REAL(DP) ::       dyn_visc, &                                   ! dynamic viscosity, kg/(m*s)
                      l_gas                                         ! Gas mean free path in air
    INTEGER  :: i,j
    REAL(DP), DIMENSION(nr_bins) :: self_coagulation, LPcoagulation
    real(dp), dimension(nr_bins) :: particle_conc_new_coag

  ! The Coagulation coefficient is calculated according to formula 13.56 in Seinfield and Pandis (2006), Page 603

  dyn_visc = 1.8D-5*(temp/298.)**0.85                                              ! Dynamic viscosity of air

  l_gas=2D0*dyn_visc/(pres*SQRT(8D0*Mair/(pi*Rg*temp)))                        ! Gas mean free path in air (m)

  slip_correction = 1D0+(2D0*l_gas/(diameter))*&
  (1.257D0+0.4D0*exp(-1.1D0/(2D0*l_gas/diameter)))                                        ! Cunninghams slip correction factor (Seinfeld and Pandis eq 9.34)

  diffusivity = slip_correction*kb*temp/(3D0*pi*dyn_visc*diameter)                 ! Diffusivity for the different particle sizes m^2/s

  speed_p = SQRT(8D0*kb*temp/(pi*particle_mass))                                   ! Speed of particles (m/s)

  free_path_p = 8D0*diffusivity/(pi*speed_p)                                              ! Particle mean free path (m)

  dist = (1D0/(3D0*diameter*free_path_p))*((diameter+free_path_p)**3D0 &
  -(diameter**2D0+free_path_p**2D0)**(3D0/2D0))-diameter                    ! mean distance from the center of a sphere reached by particles leaving the sphere's surface (m)

  DO i = 1,nr_bins
    Beta_Fuchs = 1D0/((diameter+diameter(i))/(diameter+diameter(i)+&
    2D0*(dist**2D0+dist(i)**2D0)**0.5D0)+8D0*(diffusivity+diffusivity(i))/&
    (((speed_p**2D0+speed_p(i)**2D0)**0.5D0)*(diameter+diameter(i))))                    ! Fuchs correction factor from Seinfeld and Pandis, 2006, p. 600

    coagulation_coef(i,:) = 2D0*pi*Beta_Fuchs*(diameter*diffusivity(i)+&
    diameter*diffusivity+diameter(i)*diffusivity+diameter(i)*diffusivity(i))             ! coagulation rates between two particles of all size combinations  (m^3/s)    
  END DO

  !self_coagulation = 0.5 * SUM (( (/ (coagulation_coef(i,i), i = 1,nr_bins) /)) * particle_conc**2)

  self_coagulation = 0.0_dp
  LPcoagulation = 0.0_dp

  do i = 1, nr_bins-1
    do j = i+1, nr_bins
      LPcoagulation(i) = LPcoagulation(i) + particle_conc(i) * ( coagulation_coef(i,j) * particle_conc(j) ) 
    end do
  end do

  do i = 1, nr_bins-1
    self_coagulation(i) = 0.5_dp * ( coagulation_coef(i,i) * particle_conc(i)**2 )
  end do

  particle_conc_new_coag=0D0 ! Initialise a new vector with the new particle concentrations
  particle_conc_new_coag(nr_bins)=particle_conc(nr_bins)
  !LPcoagulation = particle_conc * SUM (( (/ (coagulation_coef(:,i), i = 2,nr_bins) /)) * (/ (particle_conc(i), i=2,nr_bins) /))
  particle_conc_new_coag = particle_conc - (LPcoagulation + self_coagulation) * timestep
  particle_conc = particle_conc_new_coag


  ! Write equations that considers how the particle number concentration in each size bin 
  !(particle_conc) is influenced by the coagulation sink (loss of smaller particles when
  ! they collide with larger ones)

  ! You can first calculate the loss (loss1) do to self-coagulation between particles in the same size bin
  ! and then calculate the loss (loss2) due to coagulation with larger particles
  ! Then add the two loss terms together loss = loss1 + loss2 
 

END SUBROUTINE Coagulation

SUBROUTINE Open_files
    OPEN(unit=101, file = 'Particles.dat')          ! Number concentrations of particles, #/cm3
    OPEN(unit=102, file = 'Nucleation_rate.dat')    ! Simulated nucleation rate, #/cm3/s
    OPEN(unit=103, file = 'Coagulation_rate.dat')   ! Simulated coagulation loss rates, #/(cm3*s)
    OPEN(unit=104, file = 'Cond_Vapours.dat')       ! Gas concentration of h2so4, oh, o3, and alpha-pinene #/cm3
    OPEN(unit=105, file = 'Diameter.dat')           ! Particle diameters 
    
    OPEN(unit=106, file = 'PN.dat')                 ! Particle diameters 
    OPEN(unit=107, file = 'PV.dat')                 ! Particle diameters
     
END SUBROUTINE Open_files   
    
SUBROUTINE Write_files
    WRITE(101, "(101E20.8)") time/86400., particle_conc / 1.d6 
    WRITE(102, "(2E20.8)") time/86400., nucleation_rate / 1.d6 
    WRITE(103, "(101E20.8)") time/86400., coag_loss / 1.d6  
    WRITE(104, "(3E20.8)") time/86400., cond_vapour / 1.d6
    WRITE(105, "(101E20.8)") time/86400., diameter*1.d9
    WRITE(106, "(101E20.8)") time/86400., PN
    WRITE(107, "(101E20.8)") time/86400., PV
END SUBROUTINE Write_files

subroutine Close_files
  close(101)
  close(102)
  close(103)
  close(104)
  close(105)
  close(106)
  close(107)
end subroutine Close_files

subroutine calculateNucleo(H2SO4Conc, RH, Temp, J2)
  implicit NONE
  real(dp), intent(in) :: H2SO4Conc
  real(dp), PARAMETER :: K = 1D-20     ![m^3 / molec / s]
  REAL(DP), intent(out) :: J2
  REAL(DP), intent(in) :: RH, Temp

  J2 = K * H2SO4Conc**2

end subroutine calculateNucleo
  
END PROGRAM Aerosol_box
