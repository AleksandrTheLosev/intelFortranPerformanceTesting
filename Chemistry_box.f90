PROGRAM chemistry

  USE chemf
  IMPLICIT NONE
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15,307)

  REAL(dp) :: Conc(neq)
  REAL(dp) :: O2, N2, Mair, H2O, TEMP
  REAL(dp), Parameter :: ppb = 1e-9_dp

  ! time management
  REAL(dp) :: time, time_start, time_end
  REAL(dp) :: dt                          ! time step, in seconds
  INTEGER  :: time_in_ms                  ! time in milliseconds
  INTEGER, PARAMETER :: one_hour = 60*60  ! in seconds
  INTEGER :: daynumber_start, daynumber
  REAL(dp) :: exp_coszen, pres

  ! latitude and longitude of Hyyt:
  REAL(dp), PARAMETER :: pi = 2*ASIN(1.)
  REAL(dp), PARAMETER :: latitude_deg = 61.8455
  REAL(dp), PARAMETER :: latitude = latitude_deg * pi/180.
  REAL(dp), PARAMETER :: longitude_deg = 24.2873
  REAL(dp), PARAMETER :: longitude = longitude_deg * pi/180.
  REAL(dp), PARAMETER :: Rgas = 8.3144598d0                   ! universal gas constant [J mol-1 K-1]
  REAL(dp), PARAMETER :: NA   = 6.022140857d23                ! Avogadro's number [molec mol-1]

  ! Start program

  ! Atmospheric oxygen, N2 and H2O are kept constant:
  pres = 1.01325e5_dp                    ! [Pa], reference pressure at surface
  TEMP = 300.                            ! Temperature in Kelvin
  Mair = pres*NA / (Rgas*temp) * 1e-6    ! Air molecules concentration [molecules/cm3]
  O2   = 0.21*Mair                       ! Oxygen concentration [molecules/cm3]
  N2   = 0.78*Mair                       ! Nitrogen concentration [molecules/cm3]
  H2O  = 1.0D16                          ! Water molecules [molecules/cm3]

  ! initial state:
  Conc     = 0.0
  Conc(1)  = 24.   * Mair * ppb     ! O3 concentration
  Conc(9)  = 100.  * Mair * ppb     ! CO
  Conc(5)  = 1.    * Mair * ppb     ! NO2
  Conc(6)  = 0.2   * Mair * ppb     ! NO
  Conc(20) = 0.5   * Mair * ppb     ! SO2
  Conc(11) = 1759. * Mair * ppb     ! CH4
  Conc(13) = 2.2   * Mair * ppb     ! C5H8
  Conc(23) = 2.2   * Mair * ppb     ! alpha-p

  daynumber_start = (31+28+31+30+31+30+31+10) ! day is 10.8.
  daynumber = daynumber_start

  time_start = 0.0
  time_end = 5 * 24 * one_hour
  dt = 10.0

  time = time_start

  exp_coszen = get_exp_coszen(time, daynumber, latitude, longitude)

  OPEN(11,file="Output/concentrations.dat",status='replace',action='write')
  OPEN(12,file="Output/radiation.dat",status='replace',action='write')
  WRITE(11,*) Conc
  WRITE(12,*) exp_coszen

  DO WHILE (time < time_end)

     exp_coszen = get_exp_coszen(time, daynumber, latitude, longitude)


     CALL chemistry_step(Conc,time,time+dt,O2,N2,Mair,H2O,TEMP,exp_coszen)
     time = time + dt
     daynumber = daynumber_start + FLOOR(time/(24 * one_hour))

     ! Write output every hour
     time_in_ms = FLOOR(1000*time)
     IF (MODULO(time_in_ms, 1000*one_hour) == 0) THEN
        WRITE(*,'(A8,F6.3,A6)') 'time = ', time/(24*one_hour), '  days'
        WRITE(11,*) Conc
        WRITE(12,*) exp_coszen
     END IF

  END DO

  CLOSE(11)
  CLOSE(12)

CONTAINS

  REAL(dp) FUNCTION get_hourangle(time)
     REAL(dp), INTENT(in) :: time
     REAL(dp), PARAMETER :: one_day = 24*one_hour
     get_hourangle = MODULO(time,one_day)/one_day * 2 * pi - pi
  END FUNCTION get_hourangle

  REAL(dp) FUNCTION solar_zenith_angle(hourangle,daynumber,latitude,longitude)
     ! http://en.wikipedia.org/wiki/Solar_elevation_angle
     ! http://en.wikipedia.org/wiki/Position_of_the_Sun
     INTEGER, INTENT(in) :: daynumber
     REAL(dp), INTENT(in) :: hourangle,latitude,longitude
     REAL(dp) :: declination,elevation
     REAL(dp), PARAMETER :: to_rad = pi/180.

     declination = -23.44 * to_rad * COS(2 * pi * (daynumber + 10)/365.)
     elevation = COS(hourangle)*COS(declination)*COS(latitude) &
                 + SIN(declination)*SIN(latitude)
     solar_zenith_angle = pi/2. - elevation
     ! Notes:
     ! - Not tested near equador or on the southern hemisphere.
     ! - solar_zenith_angle can be larger than pi/2, it just means that the sun is below horizon.
     ! - solar_zenith_angle assumes time is in local solar time, which is usually not exactly true

  END FUNCTION solar_zenith_angle

  REAL(dp) FUNCTION get_exp_coszen(time,daynumber,latitude,longitude)
     REAL(dp), INTENT(in) :: time,latitude,longitude
     INTEGER, INTENT(in) :: daynumber
     REAL(dp) :: hourangle,zenith,coszen
     hourangle = get_hourangle(time)
     zenith = solar_zenith_angle(hourangle,daynumber,latitude,longitude)
     coszen = COS(zenith)
     IF (coszen > 0) THEN  ! sun is above horizon
        get_exp_coszen = EXP(-0.575/coszen)
     ELSE
        get_exp_coszen = 0
     ENDIF
  END FUNCTION get_exp_coszen

END PROGRAM chemistry
