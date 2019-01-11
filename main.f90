program main

use chemf

implicit none

logical :: use_emission   = .true.
logical :: use_chemistry  = .true.
logical :: use_deposition = .false.
logical :: use_aerosol    = .true.
character(len=255), parameter :: input_dir  = './input'
character(len=255), parameter :: output_dir = './output'

!-----------------------------------------------------------------------------------------
! Constants
!-----------------------------------------------------------------------------------------
! Double precision
! http://en.wikipedia.org/wiki/Double_precision_floating-point_format
integer, parameter :: dp = selected_real_kind(15, 307)

! Physics constants
real(dp), parameter :: PI     = 2*asin(1.0_dp)                     !                  the constant pi
real(dp), parameter :: grav   = 9.81_dp                            ! [m s-2],         gravitation
real(dp), parameter :: Rgas   = 8.3144598_dp                       ! [J mol-1 K-1],   universal gas constant
real(dp), parameter :: NA     = 6.022140857e23_dp                  ! [molec mol-1],   Avogadro's number 
real(dp), parameter :: mm_air = 28.96e-3_dp                        ! [kg mol-1],      mean molar mass of air
real(dp), parameter :: kb     = 1.38064852e-23_dp                  ! [m2 kg s-2 K-1], Boltzmann constant
real(dp), parameter :: Cp     = 1012.0_dp                          ! [J kg-1 K-1],    air specific heat at constant pressure,
real(dp), parameter :: p00    = 1.01325e5_dp                       ! [Pa],            reference pressure at surface
real(dp), parameter :: nu_air = 1.59e-5_dp                         ! [m2 s-1],        kinematic viscosity of air
real(dp), parameter :: Omega  = 2*PI/(24.0_dp*60.0_dp*60.0_dp)     ! [rad s-1],       Earth angular speed
real(dp), parameter :: lambda = 300.0_dp                           ! [m],             maximum mixing length
real(dp), parameter :: vonk   = 0.4_dp                             ! [],              von Karman constant

! Latitude and longitude of Hyytiala
real(dp), parameter :: latitude_deg  = 61.8455d0                   ! [degN]
real(dp), parameter :: longitude_deg = 24.2833d0                   ! [degE]
real(dp), parameter :: latitude      = latitude_deg  * PI/180.0d0  ! [rad]
real(dp), parameter :: longitude     = longitude_deg * PI/180.0d0  ! [rad]

! Coriolis parameter at Hyytiala
real(dp), parameter :: fcor = 2*Omega*sin(latitude)  

!-----------------------------------------------------------------------------------------
! Grid parameters and variables
!-----------------------------------------------------------------------------------------
integer, parameter :: nz = 50  ! [-], number of height levels

! Model height levels, [m]
real(dp), parameter, dimension(nz) :: &
  hh = (/    0,   10,   20,   30,   40,   50,   60,   70,   80,   90, &
           100,  120,  140,  160,  180,  200,  230,  260,  300,  350, &
           400,  450,  500,  550,  600,  650,  700,  800,  900, 1000, &
          1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, &
          2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000 /)

real(dp), dimension(nz  ) :: uwind,   &             ! [m s-1], u component of wind
                             vwind,   &             ! [m s-1], v component of wind
                             theta,   &             ! [K],     potential temperature
                             temp,    &             ! [K],     air temperature
                             pres,    &             ! [Pa],    air pressure
                             O2,      &
                             N2,      &
                             Mair

  REAL(DP) :: H2O

real(dp), parameter :: hc = 10.0_dp  ! [m], canopy height
INTEGER, PARAMETER ::  nr_bins = 100                  ! Number of particle size bins
INTEGER, PARAMETER ::  nr_cond = 2                    ! Number of condensable vapours
REAL(DP), DIMENSION(nr_bins, nz) :: diameter  ,&          ! Diameter of each size bin
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
!-----------------------------------------------------------------------------------------
! Time variables
!-----------------------------------------------------------------------------------------
integer, parameter :: one_hour = 60*60   ! [s], one hour in seconds

real(dp) :: time                         ! [s], current time
real(dp) :: time_start, time_end         ! [s], start and end time

real(dp) :: dt                           ! [s], time step for main loop, usually is equal to meteorology time step
real(dp) :: dt_emis                      ! [s], time step for emission calculation
real(dp) :: dt_chem                      ! [s], time step for chemistry calculation
real(dp) :: dt_depo                      ! [s], time step for deposition calculation
real(dp) :: dt_aero                      ! [s], time step for aerosol calculation
real(dp) :: dt_output                    ! [s], time step for output

real(dp) :: time_start_emission          ! [s], time to start calculating emission
real(dp) :: time_start_chemistry         ! [s], time to start calculating chemistry
real(dp) :: time_start_deposition        ! [s], time to start calculating deposition
real(dp) :: time_start_aerosol           ! [s], time to start calculating aerosol

integer :: daynumber_start               ! [day], start day of year
integer :: daynumber                     ! [day], current day of year

integer :: counter                       ! [-], counter of time steps

integer :: i,j  ! used for loops (j)

real(dp), DIMENSION(nz-1) :: Km, Kh, Ri


real(dp) :: exp_coszen
REAL(dp) :: Conc(neq, nz), Conc_ex(neq, nz)
REAL(dp), Parameter :: ppb = 1e-9_dp
REAL(DP) :: Emi_alp,Emi_iso
REAL(DP) :: Emis_alp, Emis_iso
real(dp), PARAMETER :: alpha = 1.0_dp

!-----------------------------------------------------------------------------------------
! Initialization
!-----------------------------------------------------------------------------------------

call time_init()                 ! initialize time
call meteorology_init()          ! initialize meteorology
call chem_init(temp, pres)

call open_files()        ! open output files
call write_files(time)   ! write initial values

Emis_iso = 0.0
Emis_alp = 0.0 


!-----------------------------------------------------------------------------------------
! Start main loop
!-----------------------------------------------------------------------------------------
do while (time <= time_end)
  !---------------------------------------------------------------------------------------
  ! Meteorology
  !---------------------------------------------------------------------------------------

  ! Calculate exp_coszen for current time
  exp_coszen = get_exp_coszen(time,daynumber,latitude)
  

  ! Set lower boundary condition
  call surface_values(theta(1), time+dt)  ! theta = temperature at the surface
  
  ! Update meteorology
  call meteor_update(uwind,vwind, theta, temp, pres)

  
  !---------------------------------------------------------------------------------------
  ! Emission
  !---------------------------------------------------------------------------------------
  ! Start to calculate emission after time_start_emission
  ! Compute emission part every dt_emis, multiplying 1000 to convert s to ms to make mod easier
  if ( use_emission .and. time >= time_start_emission ) then
    if ( mod( nint((time - time_start_emission)*1000.0d0), nint(dt_emis*1000.0d0) ) == 0 ) then
      ! Calculate emission rates

      !exp_coszen = get_exp_coszen(time,daynumber,latitude)

      call calculateEmission( temp,       &   !Vector of T. We need only second layer
                              exp_coszen, &   !Mr. Sun
                              Emi_iso,   &   !Update initial values for Isoprene
                              Emi_alp )      !Update initial values for a-penene
    end if
  end if

  if ( use_emission .and. (.not. use_chemistry) ) then
    ! Add emission to the number concentrations of compoundsmak
  end if


  !---------------------------------------------------------------------------------------
  ! Deposition
  !---------------------------------------------------------------------------------------
  ! Start to calculate gas dry deposition velocity after time_start_deposition
  ! Compute deposition part every dt_depo, multiplying 1000 to convert s to ms to make mod easier
  if ( use_deposition .and. time >= time_start_deposition ) then
    if ( mod( nint((time - time_start_deposition)*1000.0d0), nint(dt_depo*1000.0d0) ) == 0 ) then
      ! Calculate deposition velocity

      ! Remove deposited concentration at level 2 which includes canopy and soil
    end if
  end if


  !---------------------------------------------------------------------------------------
  ! Chemistry
  !---------------------------------------------------------------------------------------
  ! Start to calculate chemical reactions only after some time to save the computation time
  ! Compute chemistry part every dt_chem, multiplying 1000 to convert s to ms to make mod easier
  if ( use_chemistry .and. time >= time_start_chemistry ) then
    if ( mod( nint((time - time_start_chemistry)*1000.0d0), nint(dt_chem*1000.0d0) ) == 0 ) then
      ! Solve chemical equations for each layer except boundaries
      call chem_update(temp, pres)
      !write(*,*) Conc(3,:)
      do i = 2, nz-1      
        ! for each layer call chemistry step in order to calculate conc vector for each spec.
        !Conc - concentration's vector for update
        !O2(i) - oxygen's vector
        if (i /= 2) then 
          Emis_iso = 0.0
          Emis_alp = 0.0 
        else
          Emis_iso = Emi_iso
          Emis_alp = Emi_alp  
        end if
          CALL chemistry_step(Conc(1:neq, i),time,time+dt_chem, O2(i), N2(i), Mair(i), H2O, temp(i), &
          exp_coszen, Emis_iso, Emis_alp, Cond_sink)
 

        
      end do

    end if  ! every dt_chem
  end if

  ! Update concentrations of gas phase compounds if any of these processes are considered
  ! Deposition should not be used alone because it calculates nothing in that case
  if (use_emission .or. use_chemistry) then
    ! Trick to make bottom flux zero

    Conc(1:neq,1) = Conc(1:neq,2)
    Conc(1:neq, nz) = 0.0
    
    ! Concentrations can not be lower than 0
    do i = 1, 25
      do j = 1, nz
        if (Conc(i,j) < 0.0_dp) then
              Conc(i,j) = 0.0_dp
        end if
      end do  
    end do

    ! Mixing of chemical species
    do i = 2, nz-1
      Conc_ex(1:neq,i) = Conc(1:neq,i) + (dt * (((Kh(i)*((Conc(1:neq,i+1)-Conc(1:neq,i))/(hh(i+1)-hh(i)))) &
        -(Kh(i-1)*((Conc(1:neq,i)-Conc(1:neq,i-1))/(hh(i)-hh(i-1))))) / ((hh(i+1)-hh(i-1))/2)))
    end do

    ! Set the constraints above again for output
    Conc(1:neq,2:nz-1) = Conc_ex(1:neq,2:nz-1)

  end if

  !---------------------------------------------------------------------------------------
  ! Aerosol
  !---------------------------------------------------------------------------------------
  ! Start to calculate aerosol processes only after some time to save the computation time
  ! Compute aerosol part every dt_aero, multiplying 1000 to convert s to ms to make mod easier
  if ( use_aerosol .and. time >= time_start_aerosol ) then
    if ( mod( nint((time - time_start_aerosol)*1000.0d0), nint(dt_aero*1000.0d0) ) == 0 ) then
      do i = 1, nz
        cond_vapour(1) = Conc(21,i) * 1D6
        cond_vapour(2) = Conc(25, i) * 1D6
        
        call calculateNucleo(dt_aero, cond_vapour(1), particle_conc(:,i))                           !NUCLEO

        call Coagulation(dt_aero, particle_conc(:,i), diameter, &                                   !COAG
        temp(i),pres(i),particle_mass,Mair(i))

        call Condensation(dt_aero, temp(i), pres(i),alpha, molecular_mass, &
        molecular_volume, molar_mass, molecular_dia, particle_mass, particle_volume, &
        particle_conc, cond_sink, diameter, cond_vapour)


      end do
    end if

    ! Trick to make bottom flux zero

    ! Concentrations can not be lower than 0 [molec m-3]

    ! Mixing of aerosol particles

    ! Set the constraints above again for output

    ! Update related values, e.g., total number concentration, total mass concentration

  end if

  !---------------------------------------------------------------------------------------
  ! Ending loop actions
  !---------------------------------------------------------------------------------------
  ! Advance to next time step
  time = time + dt

  ! Write data every dt_output [s]
  if ( mod( nint((time - time_start)*1000.0d0), nint(dt_output*1000.0d0) ) == 0 ) then
    write(*, '(a8,f8.3,a8)') 'time = ', time/one_hour, '   hours'
    call write_files(time)
  end if

  ! Count loop number
  counter = counter + 1

end do

!-----------------------------------------------------------------------------------------
! Finalization
!-----------------------------------------------------------------------------------------

! Close all the opened files
call close_files()

write(*,'(a,i6,a)') "Complete with ",counter, " time steps" 


contains

!-----------------------------------------------------------------------------------------
! subroutine open_files()
!
! Open needed files
!-----------------------------------------------------------------------------------------
subroutine open_files()
  logical :: dir_exist

  ! Create a new directory if it does not exist
  inquire(file=trim(adjustl(output_dir)), exist=dir_exist)
  if (.not. dir_exist) then
    ! This line may change for different operating systems
    call system('mkdir ' // trim(adjustl(output_dir)))
  end if

  ! Open files to write output results
  open( 4, file=trim(adjustl(output_dir))//'/canopy_T.dat' , status='replace', action='write')
  open( 5, file=trim(adjustl(output_dir))//'/exp_coszen.dat' , status='replace', action='write')
  open( 7, file=trim(adjustl(output_dir))//'/ground_T.dat'   , status='replace', action='write')

  open( 8, file=trim(adjustl(output_dir))//'/time.dat' , status='replace', action='write')
  open( 9, file=trim(adjustl(output_dir))//'/hh.dat'   , status='replace', action='write')
  open(10, file=trim(adjustl(output_dir))//'/uwind.dat', status='replace', action='write')
  open(11, file=trim(adjustl(output_dir))//'/vwind.dat', status='replace', action='write')
  open(12, file=trim(adjustl(output_dir))//'/theta.dat', status='replace', action='write')
  open(13, file=trim(adjustl(output_dir))//'/K.dat', status='replace', action='write')

  open(14, file=trim(adjustl(output_dir))//'/Conc_OH.dat', status='replace', action='write')
  open(15, file=trim(adjustl(output_dir))//'/Conc_H20.dat', status='replace', action='write')
  open(16, file=trim(adjustl(output_dir))//'/Conc_H2SO4.dat', status='replace', action='write')
  open(17, file=trim(adjustl(output_dir))//'/Conc_ELVOC.dat', status='replace', action='write')

  open(18, file=trim(adjustl(output_dir))//'/emi_iso.dat', status='replace', action='write')
  open(19, file=trim(adjustl(output_dir))//'/emi_alp.dat', status='replace', action='write')

  open(20, file=trim(adjustl(output_dir))//'/emi_iso_Conc.dat', status='replace', action='write')
  open(21, file=trim(adjustl(output_dir))//'/emi_alp_Conc.dat', status='replace', action='write')
  open(22, file=trim(adjustl(output_dir))//'/O3_Conc.dat', status='replace', action='write')

end subroutine open_files


!-----------------------------------------------------------------------------------------
! subroutine write_files(time)
!
! Write data to files at time
!-----------------------------------------------------------------------------------------
subroutine write_files(time)
  real(dp) :: time  ! current time
  character(255) :: outfmt_one_scalar, outfmt_two_scalar, outfmt_level, outfmt_mid_level

  ! Output real data with scientific notation with 16 decimal digits
  outfmt_one_scalar = '(es25.16e3)'                               ! for scalar
  write(outfmt_level     , '(a, i3, a)') '(', nz  , 'es25.16e3)'  ! for original levels
  write(outfmt_mid_level , '(a, i3, a)') '(', nz-1, 'es25.16e3)'  ! for middle levels
  write(outfmt_two_scalar, '(a, i3, a)') '(', 2   , 'es25.16e3)'  ! for two scalars

  ! Only output hh once
  if (time == time_start) then
    write(9, outfmt_level) hh
  end if

  ! Output every output time step
  write( 5, outfmt_level     ) exp_coszen
  write( 7, outfmt_level     ) temp(1) !
  write( 4, outfmt_level     ) temp(2) ! for canopy T
  write( 8, outfmt_one_scalar) time/(24*one_hour)  ! [day]
  write(10, outfmt_level     ) uwind
  write(11, outfmt_level     ) vwind
  write(12, outfmt_level     ) theta
  write(13, outfmt_level     ) Km
  write(14, outfmt_level     ) Conc( 3,:)               ! For OH    Concentration
  write(15, outfmt_level     ) Conc( 8,:)               ! For HO2   Conccentration
  write(16, outfmt_level     ) Conc(21,:)               ! For H2SO4 Concentration
  write(17, outfmt_level     ) Conc(25,:)               ! For ELVOC Concentration
  write(18, outfmt_level     ) Emi_iso
  write(19, outfmt_level     ) Emi_alp

  write(20, outfmt_level     ) Conc(13,:)
  write(21, outfmt_level     ) Conc(23,:)

  write(22, outfmt_level     ) Conc(1,:)               ! For ELVOC Concentration
end subroutine write_files


!-----------------------------------------------------------------------------------------
! subroutine Close_Files()
!
! Close files
!-----------------------------------------------------------------------------------------
subroutine close_files()
  close( 4)
  close( 5)
  close( 7)
  close( 8)
  close( 9)
  close(10)
  close(11)
  close(12)
  close(13)
  close(14)
  close(15)
  close(16)
  close(17)
  close(18)
  close(19)
  close(20)
  close(21)
  close(22)
end subroutine close_files


!-----------------------------------------------------------------------------------------
! subroutine time_init()
!
! Time initiation
!-----------------------------------------------------------------------------------------
subroutine time_init()
  ! Basic time variables
  time_start = 0.0d0
  time_end   = 5.0d0 * 24.0d0 * one_hour
  time       = time_start

  ! Time steps
  dt        = 0.5d0
  dt_emis   = 0.5d0
  dt_chem   = 10.0d0
  dt_depo   = 10.0d0
  dt_aero   = 10.0d0
  dt_output = 3600.0d0

  ! Day number
  daynumber_start = 31+28+31+30+31+30+31+10  ! day is Aug. 10
  daynumber       = daynumber_start

  ! Start time for each process
  time_start_emission   = 3*24*one_hour
  time_start_chemistry  = 3*24*one_hour
  time_start_deposition = 3*24*one_hour
  time_start_aerosol    = 3*24*one_hour

  ! Loop number
  counter = 0
end subroutine time_init


!-----------------------------------------------------------------------------------------
! subroutine meteorology_init()
!
! Meteorology initiation
!-----------------------------------------------------------------------------------------
subroutine meteorology_init()
  ! Wind velocity
  uwind         = 0.0d0
  uwind(nz)     = 10.0d0
  uwind(2:nz-1) = uwind(nz) * hh(2:nz-1)/hh(nz)

  vwind = 0.0d0

  ! Potential temperature
  theta     = 273.15d0 + 25.0d0
  theta(nz) = 273.15d0 + 30.0d0

  ! Air temperature and pressure
  temp = theta - (grav/Cp)*hh
  pres = barometric_law(p00, temp, hh)
end subroutine meteorology_init


subroutine chem_init(temp, pres)
  implicit none
  REAL(DP), DIMENSION(nz), intent(in) :: pres, temp
  ! Atmospheric oxygen, N2 and H2O are kept constant:
  Mair = pres*NA / (Rgas*temp) * 1e-6    ! Air molecules concentration [molecules/cm3]
  O2   = 0.21*Mair                       ! Oxygen concentration [molecules/cm3]
  N2   = 0.78*Mair                       ! Nitrogen concentration [molecules/cm3]
  H2O  = 1.0D16                          ! Water molecules [molecules/cm3]


  ! initial concentration state:
  Conc     = 0.0
  Conc(1,1:nz)  = 24.   * Mair * ppb     ! O3 concentration
  Conc(9,1:nz)  = 100.  * Mair * ppb     ! CO
  Conc(5,1:nz)  = 0.2    * Mair * ppb     ! NO2
  Conc(6,1:nz)  = 0.07   * Mair * ppb     ! NO
  Conc(20,1:nz) = 0.5   * Mair * ppb     ! SO2
  Conc(11,1:nz) = 1759. * Mair * ppb     ! CH4
  Conc(13, 1:nz) = 0.0
  Conc(23, 1:nz) = 0.0
end subroutine chem_init

subroutine chem_update(temp, pres)
  implicit none 

  real(dp), DIMENSION(nz), intent(in) :: pres 
  real(dp), DIMENSION(nz), intent(in) :: temp 

  ! Atmospheric oxygen, N2 and H2O are kept constant:
  Mair = pres*NA / (Rgas*temp) * 1e-6    ! Air molecules concentration [molecules/cm3]
  O2   = 0.21*Mair                       ! Oxygen concentration [molecules/cm3]
  N2   = 0.78*Mair                       ! Nitrogen concentration [molecules/cm3]
  H2O  = 1.0D16                          ! Water molecules [molecules/cm3]

  ! initial concentration state:
  Conc(1,1:nz)  = 24.   * Mair * ppb     ! O3 concentration
  Conc(9,1:nz)  = 100.  * Mair * ppb     ! CO
  Conc(5,1:nz)  = 0.2    * Mair * ppb     ! NO2
  Conc(6,1:nz)  = 0.07   * Mair * ppb     ! NO
  Conc(20,1:nz) = 0.5   * Mair * ppb     ! SO2
  Conc(11,1:nz) = 1759. * Mair * ppb     ! CH4

end subroutine chem_update



!-----------------------------------------------------------------------------------------
! Get the surface values from the input data file
! Now only the temperature is used.
!-----------------------------------------------------------------------------------------
subroutine surface_values(temperature, time)

  ! (Note: can also get water concentrantion, in ppt, if modify this
  ! subroutine to also use column 8)
  !
  ! Data is taken from:
  ! http://avaa.tdata.fi/web/smart

  real(dp), intent(in)            :: time ! input, in seconds
  real(dp), intent(out)           :: temperature ! output, in Kelvin
  logical, save                   :: first_time = .true.
  real(dp), dimension(8,50), save :: surface_data
  real(dp), dimension(50), save   :: temperature_data
  real(dp), parameter             :: seconds_in_day = 24*60*60
  real(dp), parameter             :: seconds_in_30min = 30*60
  integer                         :: index
  real(dp) :: time24h, time30min, time24plus15, temp1, temp2, x

  ! Only when called for the first time, read in data from file
  ! With this trick, we don't need to open the file in the main program
  IF (first_time) THEN
     open(30, file=trim(adjustl(input_dir))//'/hyytiala_2011_8_10_t_h2o.dat', status='old')
     read(30, *) surface_data
     temperature_data(1:50) = surface_data(7,1:50) ! in Celcius
     first_time = .false.
  end IF

  time24h = modulo(time, seconds_in_day) ! time modulo 24 hours
  time24plus15 = time24h + 15*60 ! time since 23:45 previous day
  time30min = modulo(time24plus15, seconds_in_30min)
  index = 1 + floor(time24plus15/seconds_in_30min)

  temp1 = temperature_data(index)
  temp2 = temperature_data(index + 1)
  x = time30min/seconds_in_30min

  ! linear interpolation between previous and next temperature data value
  temperature = temp1 + x*(temp2 - temp1) + 273.15_dp  ! now in Kelvin

end subroutine surface_values


!-----------------------------------------------------------------------------------------
! Calculate the radiation related quantities
!-----------------------------------------------------------------------------------------
REAL(dp) function get_exp_coszen(time,daynumber,latitude)
  REAL(dp), INTENT(in) :: time,latitude
  INTEGER, INTENT(in) :: daynumber
  REAL(dp) :: hourangle,zenith,coszen

  hourangle = get_hourangle(time)
  zenith = solar_zenith_angle(hourangle,daynumber,latitude)
  coszen = COS(zenith)
  IF (coszen > 0) THEN  ! sun is above horizon
    get_exp_coszen = EXP(-0.575_dp/coszen)
  ELSE
    get_exp_coszen = 0.0_dp
  ENDIF
END function get_exp_coszen


real(dp) function get_hourangle(time)
  real(dp), intent(in) :: time
  real(dp), parameter :: one_day = 24*one_hour
  get_hourangle = modulo(time,one_day)/one_day * 2 * pi - pi
end function get_hourangle

REAL(dp) FUNCTION solar_zenith_angle(hourangle,daynumber,latitude)
  ! http://en.wikipedia.org/wiki/Solar_elevation_angle
  ! http://en.wikipedia.org/wiki/Position_of_the_Sun
  INTEGER, INTENT(in) :: daynumber
  REAL(dp), INTENT(in) :: hourangle,latitude
  REAL(dp) :: declination,elevation
  REAL(dp), PARAMETER :: to_rad = pi/180.0_dp

  declination = -23.44_dp * to_rad * COS(2 * pi * (daynumber + 10)/365.0_dp)
  elevation = COS(hourangle)*COS(declination)*COS(latitude) &
              + SIN(declination)*SIN(latitude)
  solar_zenith_angle = pi/2.0_dp - elevation
  ! Notes:
  ! - Not tested near equador or on the southern hemisphere.
  ! - solar_zenith_angle can be larger than pi/2, it just means that the sun is below horizon.
  ! - solar_zenith_angle assumes time is in local solar time, which is usually not exactly true
END FUNCTION solar_zenith_angle


!-----------------------------------------------------------------------------------------
! Other functions
!-----------------------------------------------------------------------------------------
function barometric_law(p00, tempK, h) result(p)
  real(dp), intent(in) :: p00, tempK(nz), h(nz)
  real(dp) :: p(nz)
  real(dp) :: dh(nz)

  dh(2:nz) = h(2:nz) - h(1:nz-1)

  p(1) = p00
  do i=2, nz
    p(i) = p(i-1)*exp(-mm_air*grav/(Rgas*(tempK(i-1)+tempK(i))/2.0d0)*dh(i))
  end do
end function barometric_law

elemental function fm(Ri)

  real(dp), intent(in) :: Ri 
    real(dp) :: fm

    if (Ri < 0.0_dp) then 
      fm = (1.0_dp - 16.0*Ri) ** 0.5_dp
    elseif (Ri >= 0.0_dp .and. Ri < 0.2_dp) then
      fm = MAX( (1.0_dp - 5.0_dp*Ri)**2, 0.1_dp)
    else
      fm = 0.1_dp
    endif
end function fm

elemental function fh(Ri)

  real(dp), intent(in) :: Ri 
    real(dp) :: fh

    if (Ri < 0.0_dp) then 
      fh = (1.0_dp - 16.0*Ri) ** 0.75_dp
    elseif (Ri >= 0.0_dp .and. Ri < 0.2_dp) then
      fh = MAX( (1.0_dp - 5.0_dp*Ri)**2, 0.1_dp)
    else
      fh = 0.1_dp
    endif
end function fh

subroutine calc_ri(uwind, vwind, theta,Ri)

  implicit none 
  real(dp), dimension(nz), intent(in) :: uwind, vwind, theta
  real(dp), DIMENSION(nz-1), intent(out) :: Ri
  REAL(dp), DIMENSION(nz-1) :: term3, term4, term1, term2, term5

  term1 = (theta(2:nz)+theta(1:nz-1))/2.0_dp 
  term2 = ( uwind(2:nz) - uwind(1:nz-1) )**2 + ( vwind(2:nz) - vwind(1:nz-1))**2
  term3 = grav/term1
  term4 = hh(2:nz) - hh(1:nz-1)
  term5 = (theta(2:nz)-theta(1:nz-1))

  !Ri = (term1 * term4 * term5 ) / term2
  Ri = term3 * term4 * term5 / term2

end subroutine calc_ri

subroutine get_k(km,kh,Ri)

  implicit none
  REAL(dp), DIMENSION(nz-1) :: LL, term1_k, term2_k
  real(dp), DIMENSION(nz-1) :: Kz
  real(dp), dimension(nz-1), intent(out) :: Km, Kh
  real(dp), DIMENSION(nz-1), intent(in) :: Ri

  !call calc_ri(uwind, vwind, theta, Ri)
 
  Kz = (vonk * (hh(2:nz)-hh(1:nz-1)/2.0))
  LL =  Kz / ( 1.0_dp + Kz / lambda)
  term1_k = ( ( uwind(2:nz) - uwind(1:nz-1) )/ ( hh(2:nz) - hh(1:nz-1) ) )
  term2_k = ( ( vwind(2:nz) - vwind(1:nz-1) )/ ( hh(2:nz) - hh(1:nz-1) ) )


  Km = LL**2 * sqrt( (term1_k)**2 + (term2_k)**2 ) * fm(Ri)
  Kh = LL**2 * sqrt( (term1_k)**2 + (term2_k)**2 ) * fh(Ri)

end subroutine get_k

subroutine meteor_update(uwind, vwind, theta, temp, pres)

  implicit none
  real(dp), dimension(nz-1) :: flux_u , flux_v, flux_theta
  real(dp), dimension(nz) :: diffusion_u, diffusion_v, diffusion_theta 
  !real(dp), dimension(nz-1), intent(out) :: Km, Kh
  real(dp), dimension(nz), intent(inout) :: uwind, vwind, theta, temp, pres
  real(dp), dimension(nz):: du_dt, dv_dt, dtheta_dt
  
  !real(dp), DIMENSION(nz-1) :: Kz
  !real(dp), DIMENSION(nz-1) :: Ri

  call calc_ri(uwind, vwind, theta, Ri)
  !print*, ri
  call get_k(km,kh,Ri)
  
 
  flux_u = Km *( ( uwind(2:nz) - uwind(1:nz-1) ) / ( hh(2:nz) - hh(1:nz-1) ) )
  diffusion_u(2:nz-1) = (flux_u(2:nz-1) -flux_u(1:nz-2))/((hh(3:nz)-hh(1:nz-2))/2.0)
  
  diffusion_u(1) = 0.0
  diffusion_u(nz) = 0.0
 
  du_dt = fcor * (vwind-0.0) + diffusion_u
  

  
  flux_v = Km* ( ( vwind(2:nz) - vwind(1:nz-1) ) / ( hh(2:nz) - hh(1:nz-1) ) )
  diffusion_v(2:nz-1) = (flux_v(2:nz-1) -flux_v(1:nz-2))/((hh(3:nz)-hh(1:nz-2))/2.0)

  diffusion_v(1) = 0.0
  diffusion_v(nz) = 0.0
  dv_dt = -fcor * (uwind-10.0) + diffusion_v
 

  flux_theta = Kh* (theta(2:nz)-theta(1:nz-1))/(hh(2:nz)-hh(1:nz-1))
  diffusion_theta(2:nz-1) = (flux_theta(2:nz-1) - flux_theta(1:nz-2))/((hh(3:nz)-hh(1:nz-2))/2.0)

  diffusion_theta(1) = 0.0
  diffusion_theta(nz) = 0.0
  dtheta_dt =  diffusion_theta


  uwind = uwind + dt * du_dt
  vwind = vwind + dt * dv_dt
  theta = theta + dt * dtheta_dt

  uwind(1)       = 0.0d0
  uwind(nz)     = 10.0d0

  vwind(1) = 0.0d0
  theta(1)     = 273.15d0 + 25.0d0
  theta(nz) = 273.15d0 + 30.0d0

  ! Update air temperature and pressure
  temp = theta - (grav/Cp)*hh
  pres = barometric_law(p00, temp, hh)
  

end subroutine meteor_update

subroutine calculateEmission(CanopyTemp, sun, isoprene_emis, apenen_emis)
  implicit none
  REAL(DP), PARAMETER :: Dm = 0.0538

  real(dp), PARAMETER :: alpha = 0.0027
  real(dp), PARAMETER :: betta = 0.09
  real(dp), PARAMETER :: delta = 1.0
  real(dp), PARAMETER :: epss = 100.0

  REAL(DP), PARAMETER :: Ts = 303.15
  REAL(DP), PARAMETER :: Tm = 314.0

  REAL (DP), PARAMETER :: Ct1 = 95.0 * 1e3
  REAL (DP), PARAMETER :: Ct2 = 230.0*1e3
  REAL (DP), PARAMETER :: CL1 = 1.006

  REAL (DP), DIMENSION(nz), intent(in) :: CanopyTemp
  real(dp), intent(in) :: sun

  real(dp), intent(out) :: isoprene_emis, apenen_emis

  real(DP) :: PAR

  real(dp) :: gamma_isoprene, gamma_apenene 

  PAR = 1000.0_dp * sun

  gamma_isoprene = alpha * CL1 * PAR / sqrt(1.0_dp + alpha**2 * PAR**2) &
                      * EXP(Ct1 * (CanopyTemp(2) - Ts) / (Rgas * CanopyTemp(2) * Ts)) &
                      / ( 1.0_dp + EXP(Ct2 * (CanopyTemp(2) - Tm) / (Rgas * CanopyTemp(2) * Ts)) )
  isoprene_emis = Dm * epss * gamma_isoprene * delta
  isoprene_emis = isoprene_emis * NA* (1.e-14 / 36.0) / 68.0
  
  gamma_apenene = exp(betta * (CanopyTemp(2) - Ts)) ! gamma of a-pinene
  apenen_emis = Dm * epss * gamma_apenene * delta
  apenen_emis = apenen_emis * NA * (1.e-14 / 36.0) / 136.0

end subroutine calculateEmission

SUBROUTINE Aerosol_init(diameter, particle_mass, particle_volume, particle_conc, &
  particle_density, nucleation_coef, molecular_mass, molar_mass, &
  molecular_volume, molecular_dia, mass_accomm)

    !! ====================== Definition of variables =====================================================================

    REAL(DP), DIMENSION(nr_bins, nz), INTENT(OUT) :: diameter       , &    ! diamter of each size bin
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
    diameter(1,:)=2D-9 
    DO i=2,nr_bins
    diameter(i,:)=diameter(i-1,1)*(2.5D-6/diameter(1,:))**(1D0/(nr_bins-1))
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


SUBROUTINE Condensation(timestep, temp, pres,alpha, molecular_mass, &
  molecular_volume, molar_mass, molecular_dia, particle_mass, particle_volume, &
  particle_conc, cond_sink, diameter, cond_vapour) ! Add more variables if you need it

    !REAL(DP), DIMENSION(nr_cond) :: density
  
    REAL(DP), DIMENSION(nr_bins), INTENT(IN) :: diameter, particle_mass
    REAL(DP), DIMENSION(nr_cond), INTENT(IN) :: molecular_mass, molecular_dia, &
    molecular_volume, molar_mass
    REAL(DP), INTENT(IN) :: timestep, temp, pres

    REAL(DP), DIMENSION(nr_bins), INTENT(INOUT) :: particle_conc

    REAL(DP), DIMENSION(2), INTENT(IN) :: cond_vapour  ! condensing vapour concentrations, which is H2SO4 and organics (ELVOC) [#/m^3]
    
    REAL(DP), DIMENSION(nr_bins), INTENT(IN)   :: particle_volume
    REAL(DP), DIMENSION(nr_bins) :: particle_volume_new
    
    REAL(DP), DIMENSION(nr_bins)   ::  slip_correction, diffusivity, speed_p, &
    particle_conc_new, x1, x2
    
    REAL(DP), DIMENSION(nr_cond)   ::  diffusivity_gas, speed_gas
    
    REAL(DP) :: dyn_visc, l_gas, dens_air
    
    INTEGER :: j
    REAL(DP), PARAMETER :: K_Boltzman = 1.38E-23

    REAL(DP), DIMENSION(nr_cond) ::  &
                              nn_i, &
                              C_eq 
    REAL(DP), intent(IN) :: alpha
    REAL(DP), DIMENSION(nr_bins) :: cond_sink
    Real(dp), dimension(nr_bins) :: nn_p,lambda_1, lambda_2, Kn_1, Kn_2, FuchsSutugin_1, FuchsSutugin_2, CR_1, CR_2, Eq_1, Eq_2
    REAL(DP) :: Mair
    
    ! Add more variabels if you need it...

    dyn_visc = 1.8D-5*(temp/298D0)**0.85D0  ! dynamic viscosity of air
    dens_air=Mair*pres/(Rgas*temp)        ! Air density
    l_gas=2D0*dyn_visc/(pres*SQRT(8D0*Mair/(pi*Rgas*temp))) ! Gas mean free path in air (m)

    slip_correction = 1D0+(2D0*l_gas/(diameter))*&
    (1.257D0+0.4D0*exp(-1.1D0/(2D0*l_gas/diameter))) ! Cunninghams slip correction factor (Seinfeld and Pandis eq 9.34) 
    
    diffusivity = slip_correction*kb*temp/(3D0*pi*dyn_visc*diameter)   ! Diffusivity for the different particle sizes m^2/s
    speed_p = SQRT(8D0*kb*temp/(pi*particle_mass))                     ! speed of particles (m/s)
  
    diffusivity_gas=5D0/(16D0*Na*molecular_dia**2D0*dens_air)*&
    SQRT(Rgas*temp*Mair/(2D0*pi)*(molar_mass+Mair)/molar_mass)            ! Diffusivity of condensable vapours (m^2 s^-1)

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
  temp,pres,particle_mass,Mair) ! Add more variables if you need it
  
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
    REAL(DP), intent(IN) :: Mair

  ! The Coagulation coefficient is calculated according to formula 13.56 in Seinfield and Pandis (2006), Page 603

  dyn_visc = 1.8D-5*(temp/298.)**0.85                                              ! Dynamic viscosity of air

  l_gas=2D0*dyn_visc/(pres*SQRT(8D0*Mair/(pi*Rgas*temp)))                        ! Gas mean free path in air (m)

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

subroutine calculateNucleo(dt_aero, H2SO4Conc, particle_conc)
  implicit NONE
  real(dp), intent(in) :: H2SO4Conc
  real(dp), PARAMETER :: K = 1D-20     ![m^3 / molec / s]
  real(dp), dimension(nr_bins),intent(out) :: particle_conc
  REAL(DP) :: J2, dt_aero

  J2 = K * H2SO4Conc**2
  particle_conc = particle_conc + (J2 * dt_aero)

  

end subroutine calculateNucleo




end program main
