  ! Copyright 2020 University of Warwick

  ! Licensed under the Apache License, Version 2.0 (the "License");
  ! you may not use this file except in compliance with the License.
  ! You may obtain a copy of the License at

  !    http://www.apache.org/licenses/LICENSE-2.0

  ! Unless required by applicable law or agreed to in writing, software
  ! distributed under the License is distributed on an "AS IS" BASIS,
  ! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  ! See the License for the specific language governing permissions and
  ! limitations under the License.
  
MODULE radiative

  USE shared_data
  USE boundary

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: trac_temperature, rad_losses, user_defined_heating

  INTEGER, PARAMETER :: n = 7   !set the number of temperature boundaries used in Q(T)
  INTEGER, PARAMETER  :: kmax = n - 1
  REAL(num), DIMENSION(:), ALLOCATABLE :: t_boundary, pow, psi
  INTEGER, DIMENSION(:), ALLOCATABLE :: alpha
  REAL(num), DIMENSION(:), ALLOCATABLE :: yk, qk, ratios
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: flare
  REAL(num) :: cool, ratio

CONTAINS

  
  
  !****************************************************************************
  ! Subroutine to calculate critical TRAC temperature 
  !****************************************************************************
    
  SUBROUTINE trac_temperature(energy)
        
    REAL(num), INTENT(IN), DIMENSION(-1:,-1:) :: energy
    REAL(num), DIMENSION(:), ALLOCATABLE :: temperature
    
    REAL(num) :: grad, l_t, t_c
    REAL(num) :: t_trac, tb

    INTEGER :: iy
    
    ALLOCATE(temperature(-1: ny+2))
    

    IF (trac_method) THEN
        
        IF (eos_number /= EOS_ION) THEN
            temperature(:) = energy(1,:) * (gamma - 1.0_num) / 2.0_num
        ELSE
            temperature(:) = (energy(1,:) - (1.0_num - xi_n(1,:)) * ionise_pot) &
                * (gamma - 1.0_num) / (2.0_num - xi_n(1,:))
        END IF

        ! Set critical temperature to minimum of 10^4K
        t_trac = 2.e4_num / temp_norm
        
        ! Calculate critical temperature as in Johnston 2020
        t_c = 0.0_num
        DO iy = 1, ny
            grad = (temperature(iy+1) - temperature(iy-1)) / (dyc(iy) + dyc(iy-1))
            l_t = temperature(iy) / MAX(ABS(grad), none_zero)
            IF ((l_t < 2.0_num * dyb(iy)) .AND. (temperature(iy) > t_c)) THEN
                t_c = temperature(iy)
            END IF
        END DO
        t_trac = MIN(MAX(t_trac, t_c), 0.2_num * MAXVAL(temperature))
        
        ! Set multiplier at boundaries and centres
        DO iy = 0, ny 
            tb = 0.5_num * (temperature(iy) + temperature(iy+1))
            IF ((tb <= t_trac) .AND. (tb >= 2.e4_num / temp_norm)) THEN
                tr_factor_b(iy) = (t_trac / tb)**2.5_num
            ELSE 
                tr_factor_b(iy) = 1.0_num
            END IF
        
            tb = temperature(iy)
            IF ((tb <= t_trac) .AND. (tb >= 2.e4_num / temp_norm)) THEN
                tr_factor_c(iy) = (t_trac / tb)**2.5_num
            ELSE
                tr_factor_c(iy) = 1.0_num
            END IF
        END DO
    ELSE
        tr_factor_b(:) = 1.0_num
        tr_factor_c(:) = 1.0_num
    END IF
    
    DEALLOCATE(temperature)
 

END SUBROUTINE trac_temperature


  SUBROUTINE setup_loss_function
    ! In this subroutine specify the radiative loss Q(T) = psi_k * T^alpha_k
    ! This must be piecewise polynomial with kmax=n-1 regions
    ! bounded by n temperatures. 
    ! Q(T) must be in S.I.
    ! In the energy equation this would appear as a pressure cooling through
    ! dp/dt = -(gamma-1) n_e n_H Q(T)
    ! Only in this form for T~>10^4 so fully ionised and n_H=n_e

    REAL(num) :: frac
    
    t_boundary(:) = (/0.02_num, 0.0398_num, 0.0794_num, 0.251_num, 0.562_num, 1.995_num, 10.0_num/) * 1e6_num

    !Define power for polynomial fit.
    !Alpha is defined as integer but RTV has one none
    !integer power so have integer alpha array and real pow array defined below.
    !The integer array is needed as alpha=1 is a special case.
    frac = -2.0_num / 3.0_num
    alpha(:) = (/0, 2, 0, -2, 0, 0/)
    pow(:) = REAL(alpha, num) + (/0.0_num, 0.0_num, 0.0_num, 0.0_num, 0.0_num, frac/)

    !Usually specify RTV in cgs
    psi(:) = (/-21.85_num, -31.0_num, -21.2_num, -10.4_num, -21.94_num, -17.73_num/)
    psi = 10**psi
    !Convert to SI
    psi = 1.e-13_num * psi
     

  END SUBROUTINE setup_loss_function


  SUBROUTINE user_defined_heating
    ! Use this to define any heating needed on each timestep
    ! This is a dumb example so change to the heating needed
    ! This example specifies the heating rate (heat_in) in S.I. then
    ! converts to Lare units

    REAL(num), DIMENSION(:,:), ALLOCATABLE :: heat_in
    REAL(num) :: flare_power
    flare_power = 0.0_num
    
    ALLOCATE (heat_in(-1:nx+2, -1:ny+2))

    ! Specify heating in S.I. units W/m^3
    DO iy = 1, ny
      
      heat_in(:,iy) = 6.0e-5_num

      IF (flare_event) THEN
        CALL flare_energy_release(iy, flare_power)
      END IF
      
      heat_in(:,iy) = heat_in(:,iy) + flare_power

      ! TRAC normalisation
      heat_in(:,iy) = heat_in(:,iy) / tr_factor_c(iy)
    END DO
    
    ! Convert to internal Lare units
    heat_in = heat_in * time_norm**3 / (rho_norm * L_norm**2)
    
    energy(:,:) = energy(:,:) + heat_in * dt / rho(:,:)
    CALL energy_bcs

  END SUBROUTINE user_defined_heating



  SUBROUTINE rad_losses

    LOGICAL :: first_call = .TRUE.

    IF (first_call) THEN
      first_call = .FALSE.
      ALLOCATE (t_boundary(1:n), alpha(1:kmax), pow(1:kmax), psi(1:kmax))
      ALLOCATE (yk(1:n), qk(1:n), ratios(1:n))
      CALL setup_loss_function
      CALL set_exact_integration_arrays
    END IF

    CALL exact_integration_method
    CALL energy_bcs

  END SUBROUTINE rad_losses



  SUBROUTINE exact_integration_method

    REAL(num) :: temp_si, inverse_t_cool, yt, fac
    INTEGER :: i, k

    DO ix = 1, nx
      DO iy = 1, ny

        temp_si = 0.5_num * temp_norm * energy(ix,iy) * (gamma - 1.0_num) 

        k = -1
        DO i = 1, kmax
          IF (temp_si > t_boundary(i) .AND. temp_si <= t_boundary(i+1)) THEN
            k = i
            EXIT
          END IF
        END DO
        IF (k .LT. 0) CYCLE
        

        IF (alpha(k) .NE. 1) THEN
          fac = 1.0_num / (1.0 - pow(k))
          yt = yk(k) + fac * ratios(k) * (1.0_num - (t_boundary(k)/temp_si)**(pow(k)-1.0_num))
        ELSE
          yt = yk(k) + ratios(k) * LOG((t_boundary(k)/temp_si))
        END IF

        inverse_t_cool = cool * rho(ix,iy)
        yt = yt +  dt * inverse_t_cool *  time_norm

        k = -1
        DO i = 1, kmax
          IF ((yt > yk(i) .AND. yt <= yk(i+1)) .OR. (yt < yk(i) .AND. yt >= yk(i+1))) THEN
            k = i
            EXIT
          END IF
        END DO
        IF (k .LT. 0) THEN
          temp_si = t_boundary(1)
          energy(ix,iy) = temp_si * 2.0_num / (gamma - 1.0_num) / temp_norm
          CYCLE
        END IF

        IF (alpha(k) .NE. 1) THEN
          fac = 1.0_num / (1.0 - pow(k))
          temp_si = t_boundary(k) * (1.0_num - (1.0_num - pow(k)) / ratios(k) * (yt - yk(k)))**fac
        ELSE
          temp_si = t_boundary(k) * EXP((yt - yk(k)) / ratios(k))
        END IF
         
        energy(ix,iy) = temp_si * 2.0_num / (gamma - 1.0_num) / temp_norm 
        
      END DO
    END DO 

  END SUBROUTINE exact_integration_method



  SUBROUTINE set_exact_integration_arrays
    ! Define arrays used in Townsend exact integration method
    ! Here qk = Lambda_k from Townsend

    REAL(num) :: fac
    INTEGER :: k

    DO k = 1, kmax
      qk(k) = psi(k) * t_boundary(k)**pow(k)
    END DO
    qk(n) = qk(n-1) * (t_boundary(n)/t_boundary(n-1))**pow(n-1)

    ratio = qk(n) / t_boundary(n)
    cool = 0.5_num * ratio * (gamma - 1.0_num) * rho_norm / (kb_si * mf * mh_si) 
    DO k = 1, n
      ratios(k) = ratio * (t_boundary(k) / qk(k))
    END DO

    yk(n) = 0.0_num
    DO k = kmax, 1, -1
      IF (alpha(k) .NE. 1) THEN
        fac = 1.0_num / (1.0 - pow(k))
        yk(k) = yk(k+1) - fac * ratios(k) * (1.0_num - (t_boundary(k)/t_boundary(k+1))**(pow(k)-1.0_num))
      ELSE
        yk(k) = yk(k+1) - ratios(k) * LOG((t_boundary(k)/t_boundary(k+1)))
      END IF
    END DO
   
  END SUBROUTINE set_exact_integration_arrays


  SUBROUTINE flare_energy_release(iy, power)
    INTEGER, INTENT(IN) :: iy
    REAL(num), INTENT(OUT) :: power
    REAL(num) :: t_start, duration, time_part, spatial_part, Q_max, flare_size, flare_energy
    
    t_start = 0.50_num
    
    ! Total flare energy in J
    flare_energy = 1.e19_num 
    
    ! Flare spatial size in m
    flare_size = 3.e6_num 
    
    ! Duration of the flare in s
    duration = 60.0_num
    
    ! Give flare power in W/m^3
    Q_max = flare_energy / (0.5_num * duration) / (flare_size * sqrt(pi)) / (pi * 1.e12_num) 

    flare_size = flare_size / L_norm
    duration = duration / time_norm

    ! Triangular time profile
    IF ((time > t_start) .AND. (time < t_start + 0.5_num * duration)) THEN
      time_part = 2.0_num / duration * (time - t_start)
    
    ELSE IF ((time > t_start + 0.5_num * duration) .AND. (time < t_start + duration)) THEN
        time_part = 1 - 2.0_num / duration * (time - (t_start + 0.5_num * duration)) 

    ELSE
        time_part = 0.0_num
    END IF

    ! Gaussian spatial profile
    spatial_part = exp(-((yc(iy) - 0.5_num) / flare_size)**2)
    
    power = Q_max * time_part * spatial_part

  END SUBROUTINE flare_energy_release 


END MODULE radiative 
