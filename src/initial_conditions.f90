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
  
MODULE initial_conditions

  USE shared_data
  USE neutral
  USE diagnostics
  USE boundary

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: set_initial_conditions

CONTAINS

  !****************************************************************************
  ! This function sets up the initial condition for the code
  ! The variables which must be set are:
  !   rho - density
  !   v{x,y,z} - Velocities in x, y, z
  !   b{x,y,z} - Magnetic fields in x, y, z
  !   energy - Specific internal energy
  !   grav - Gravity
  !
  ! You may also need the neutral fraction. This can be calculated by a
  ! function call to get_neutral(temperature, rho). This routine is in
  ! core/neutral.f90 and requires the local temperature and mass density.
  ! For example to set xi_n to the neutral fraction use:
  !   xi_n = get_neutral(temperature, rho)
  ! 
  ! If using Hall_MHD then you must specific lambda_i in this routine
  !****************************************************************************

  SUBROUTINE set_initial_conditions

    REAL(num), DIMENSION(:,:), ALLOCATABLE :: temperature
    REAL(num) :: xi_v, amp, centre, width

    INTEGER :: loop
    INTEGER :: ix, iy, iy1
    REAL(num) :: a1, a2, dg, a, b, c
    REAL(num) :: legs_size, wtr, ycor, gravity_value, Tph, Tcor 
    
    REAL(num), DIMENSION(:), ALLOCATABLE :: yc_global, dyb_global, dyc_global
    REAL(num), DIMENSION(:), ALLOCATABLE :: grav_ref, temp_ref, rho_ref
    REAL(num), DIMENSION(:), ALLOCATABLE :: beta_ref, mag_ref, mu_m

    ALLOCATE(yc_global(-1:ny_global+1))
    ALLOCATE(dyb_global(-1:ny_global+1))
    ALLOCATE(dyc_global(-1:ny_global+1))
    ALLOCATE(grav_ref(-1:ny_global+2))
    ALLOCATE(temp_ref(-1:ny_global+2))
    ALLOCATE(rho_ref(-1:ny_global+2))
    ALLOCATE(mag_ref(-1:ny_global+2))
    ALLOCATE(beta_ref(-1:ny_global+2))
    ALLOCATE( mu_m(-1:ny_global+2))


    legs_size = 10.0e6_num / L_norm
    wtr = 1.e6_num / L_norm
    gravity_value = 274.0_num / (L_norm / time_norm**2)
    Tph = 2.e4_num / (mf * mh_si * L_norm**2 / time_norm**2 / kb_si)
    Tcor = 2.2e6_num / (mf * mh_si * L_norm**2 / time_norm**2 / kb_si)

    ! Below are all the variables which must be defined and their sizes

    vx = 0.0_num
    vy = 0.0_num
    vz = 0.0_num
    bx = 0.0_num
    by = 0.0_num
    bz = 0.0_num

    ! Fill in yc_global with the positions central to the yb_global points
    DO iy = -1, ny_global + 1
      yc_global(iy) = 0.5_num * (yb_global(iy-1) + yb_global(iy))
    END DO

    ! Fill in dyb_global and dyc_global
    DO iy = -1, ny_global
      dyb_global(iy) = yb_global(iy) - yb_global(iy-1)
      dyc_global(iy) = yc_global(iy+1) - yc_global(iy)
    END DO

    ! Fill in the reference gravity array 
    grav_ref = gravity_value 
    a1 = legs_size
    
    grav_ref(-1) = grav_ref(0)
    grav_ref(ny_global+1:ny_global+2) = grav_ref(ny_global)

   ! Calculate the density profile, starting from the refence density at the
   ! photosphere/chromosphere and calculating up
    rho_ref = 1.0_num
    mu_m = 1.0_num
    IF (eos_number == EOS_IDEAL .AND. (.NOT. neutral_gas)) mu_m = 0.5_num

   ! Go from photosphere/chromosphere up (along the loop)
    a = -(Tcor - Tph) / (a1 - 0.5_num * y_max)**2
    b = -2.0_num * a * 0.5_num * y_max
    c = Tcor + a * (0.5_num * y_max)**2
    DO iy = -1, ny_global + 1
!      IF (yc_global(iy) < y_max / 2.0_num) THEN
!        temp_ref(iy) = Tph + 0.5_num * (Tcor - Tph) * (TANH((yc_global(iy) - ycor) / wtr) + 1.0_num)
!      ELSE
!        temp_ref(iy) = Tph + 0.5_num * (Tcor - Tph) * (TANH((-yc_global(iy) - ycor + y_max) / wtr) + 1.0_num)
!      END IF
       IF (yc_global(iy) <= a1 .OR. yc_global(iy) >= a2) THEN
            temp_ref(iy) = Tph
       ELSE
            temp_ref(iy) = a * yc_global(iy)**2 + b * yc_global(iy) + c
       END IF
    END DO
    temp_ref(ny_global+1:ny_global+2) = temp_ref(ny_global)
    
    ! Now move from the photosphere/chromosphere up (along the loop)
    DO iy = 2, ny_global
       IF (yc_global(iy) >= 0.0_num) THEN
         iym = iy - 1
         dg = 1.0_num / (dyb_global(iy) + dyb_global(iym))

         rho_ref(iy)  = rho_ref(iym) * (temp_ref(iym) &
             * 1.0_num / dyc_global(iym) / mu_m(iym) &
             - grav_ref(iym) * dyb_global(iym) * dg)

         rho_ref(iy)  = rho_ref(iy) / (temp_ref(iy) &
             * 1.0_num  / dyc_global(iym) / mu_m(iy) &
             + grav_ref(iym) * dyb_global(iy) * dg)
       END IF
     END DO

   rho_ref(ny_global+1:ny_global+2) = rho_ref(ny_global)

  ! Fill in all the final arrays from the ref arrays
  iy1 = n_global_min(2) - 1

  DO iy = -1, ny + 2
    grav(iy) = grav_ref(iy1)
    DO ix = -1, nx + 2
      rho(ix,iy) = rho_ref(iy1)
      energy(ix,iy) = temp_ref(iy1)

      IF (eos_number /= EOS_IDEAL) THEN
        xi_v = get_neutral(energy(ix,iy), rho(ix,iy))
      ELSE
        IF (neutral_gas) THEN
          xi_v = 1.0_num
        ELSE
          xi_v = 0.0_num
        END IF
      END IF

      energy(ix,iy) = (energy(ix,iy) * (2.0_num - xi_v) &
          + (1.0_num - xi_v) * ionise_pot * (gamma - 1.0_num)) &
          / (gamma - 1.0_num)
    END DO
    iy1 = iy1 + 1
  END DO

  DO ix = -1, nx + 2
        energy(ix, ny+2) = energy(ix, ny+1) 
  END DO
  
  
!  DO iy = -1, ny + 2
!    DO ix = -1, nx + 2
!        energy(ix, iy) = energy(ix, iy) &
!             * (1 + 10 * EXP(-((yc_global(iy) - 0.5_num) / 0.005_num) **2))
!    END DO
!  END DO



  DEALLOCATE(yc_global, dyb_global, dyc_global, mu_m)
  DEALLOCATE(grav_ref, temp_ref, rho_ref, beta_ref, mag_ref)

  ! CALL add_probe(0.0_num, 0.0_num)

  END SUBROUTINE set_initial_conditions



   !  Alfven wave excitation
  ! do n1 = -1, nx+2
  !    do n2 = -1, ny+2 
  !    energy(n1,n2 = energy(n1,n2) & 
  !    + 0.01_num * energy(n1,n2) *&
  !    EXP(-((n1-nx/2.0_num)**2 + (n2-ny/2.0_num)**2)/4.0_num)
  !    end do
  ! end do



    !  Magnetoacoustic wave excitation       
!   do n1 = -2, nx+2
!      do n2 = -2, ny+2 
!     !vx(n1,n2) = EXP(-((n1-nx/2.0_num)**2 + (n2-ny/2.0_num)**2)/4.0_num)
!     
!     vy(n1,n2) = (n2-ny/2.0_num)*EXP(-((n1-nx/2.0_num)**2 + (n2-ny/2.0_num)**2)/4.0_num)
!     
!     !vz(n1,n2) = EXP(-((n1-nx/2.0_num)**2 + (n2-ny/2.0_num)**2)/4.0_num)
!      end do
!    end do
!
!
  SUBROUTINE potential_field_analytic()

        REAL(num), DIMENSION(:,:), ALLOCATABLE :: phi
        INTEGER :: i, j

        ALLOCATE(phi(-1:nx+2, -1:ny+2))


        phi(:,:) = 0.0_num
        DO i = -1, nx+2
          DO j = ny+2, -1, -1
              IF (yb_global(j) >= -10.0_num) THEN
                  phi(i,j) = sin(pi / 2 * xb_global(i) / xb_global(nx)) &
                      * exp(- pi / 2 * (yb_global(j) + 10) / xb_global(nx))
                  !phi(i,j) = - tanh(pi / 2 * xb_global(i) / xb_global(nx)) * yb_global(j)  
                  !phi(i,j) = atan((yb_global(j) - 10.0-num) / xb_global(i)) 
                
              ELSE 
                  phi(i,j) = phi(i,j+1)
              END IF
          END DO  
        END DO    
        
        DO iy = 0, ny
            DO ix = 0, nx
                bx(ix,iy) = -(phi(ix+1,iy)-phi(ix,iy))/dxc(ix)
            END DO
            bx(0, iy) = bx(1, iy) 
            bx(nx+1, iy) = bx(nx, iy)
        END DO
    
        DO ix = 0, nx
            DO iy = 0, ny
                by(ix,iy) = -(phi(ix,iy+1)-phi(ix,iy))/dyc(iy)
            END DO
            by(ix, -1) = by(ix, 0)
            by(ix, ny+1) = by(ix, ny)
        END DO
  
        CALL bfield_bcs

        DEALLOCATE(phi)

    END SUBROUTINE potential_field_analytic



  SUBROUTINE potential_field()

      REAL(num), DIMENSION(:,:), ALLOCATABLE :: phi
      REAL(num) :: w, errmax, error, residual, fractional_error
      REAL(num) :: by_min, by_min_local
      REAL(num) :: by_max, by_max_local
      INTEGER :: loop, x1, y1, redblack, i, j, n
      LOGICAL :: converged
      
      ALLOCATE(phi(-1:nx+2,-1:ny+2))
      phi(:,:) = 0.0_num
      CALL phi_mpi
      
      converged = .FALSE.
      w = 2.0_num / (1.0_num + SIN(pi / REAL(nx_global,num)))
      fractional_error = 1.e-10_num

      !Iterate to get phi^{n+1} by SOR Gauss-Seidel
      iterate: DO loop = 1, 10000000
      errmax = 0.0_num
      error = 0.0_num
      y1 = 1
      DO redblack = 1, 2
       x1 = y1
       DO iy = 1, ny 
         iym = iy - 1
         iyp = iy + 1
         DO ix = x1, nx, 2
           ixm = ix - 1
           ixp = ix + 1
           residual = &
               ((phi(ixp,iy) - phi(ix,iy))/dxc(ix) - (phi(ix,iy) - phi(ixm,iy))/dxc(ixm)) / dxb(ix) &
             + ((phi(ix,iyp) - phi(ix,iy))/dyc(iy) - (phi(ix,iy) - phi(ix,iym))/dyc(iym)) / dyb(iy)
           residual = residual / ((1.0_num/dxc(ix) +1.0_num/dxc(ixm))/dxb(ix) &
                              +  (1.0_num/dyc(iy) +1.0_num/dyc(iym))/dyb(iy))
           phi(ix,iy) = phi(ix,iy) + w * residual 
           error = ABS(residual) 
           errmax = MAX(errmax, error)
         END DO
         CALL phi_mpi
         x1 = 3 - x1
       END DO
       CALL phi_mpi
       y1 = 3 - y1
     END DO
     CALL MPI_ALLREDUCE(errmax, error, 1, mpireal, MPI_MAX, comm, errcode)
     IF (rank == 0 .AND. (MOD(loop,1000).EQ.0)) print *, 'loop, residual = ', loop, error
     IF (error < fractional_error) THEN
       converged = .TRUE.
       EXIT iterate
     END IF
    END DO iterate

    IF (rank == 0 .AND. .NOT. converged) PRINT*, 'potential_field failed'

      DO iy = 0, ny
        DO ix = 0, nx
          bx(ix,iy) = -(phi(ix+1,iy)-phi(ix,iy))/dxc(ix)
        END DO
      END DO
    
      DO ix = 0, nx
        DO iy = 0, ny
          by(ix,iy) = -(phi(ix,iy+1)-phi(ix,iy))/dyc(iy)
        END DO
      END DO
      
      
      CALL bfield_bcs

      DEALLOCATE(phi)

    CONTAINS

      SUBROUTINE phi_mpi

        REAL(num) :: total_flux, local_flux

        CALL MPI_SENDRECV( &
            phi(   1,-1), 1, cell_xface, proc_x_min, tag, &
            phi(nx+1,-1), 1, cell_xface, proc_x_max, tag, &
            comm, status, errcode)
        CALL MPI_SENDRECV( &
            phi(nx-1,-1), 1, cell_xface, proc_x_max, tag, &
            phi(  -1,-1), 1, cell_xface, proc_x_min, tag, &
            comm, status, errcode)

        CALL MPI_SENDRECV( &
            phi(-1,   1), 1, cell_yface, proc_y_min, tag, &
            phi(-1,ny+1), 1, cell_yface, proc_y_max, tag, &
            comm, status, errcode)
        CALL MPI_SENDRECV( &
            phi(-1,ny-1), 1, cell_yface, proc_y_max, tag, &
            phi(-1,  -1), 1, cell_yface, proc_y_min, tag, &
            comm, status, errcode)

        !Dipolar flux
        local_flux = 0.0_num
        IF (proc_y_min == MPI_PROC_NULL) THEN
          phi(1:nx,0) = phi(1:nx,1) + dyc(1) * EXP(-((xc(1:nx) - 30.0_num) / 5.0_num)**2) 
          phi(1:nx,0) = phi(1:nx,0) - dyc(1) * EXP(-((xc(1:nx) + 30.0_num) / 5.0_num)**2)
          local_flux = SUM(dxb(1:nx) * (EXP(-((xc(1:nx) - 30.0_num) / 5.0_num)**2) &
              - EXP(-((xc(1:nx) + 30.0_num) / 5.0_num)**2)))
          phi(1:nx,-1) = phi(1:nx,0)
         
        END IF
        CALL MPI_ALLREDUCE(local_flux, total_flux, 1, mpireal, MPI_SUM, comm, errcode)
        IF (proc_y_min == MPI_PROC_NULL) THEN
          phi(1:nx,0) = phi(1:nx,0) - dyc(1) * total_flux / length_x
          phi(1:nx,-1) = phi(1:nx,0)
        END IF        

        IF (proc_y_max == MPI_PROC_NULL) THEN
          phi(:,ny+1) = 1.0_num
          phi(:,ny+2) = 1.0_num
        END IF        
        IF (proc_x_min == MPI_PROC_NULL) THEN
          phi(0,:) = phi(1,:) 
          phi(-1,:) = phi(1,:)
        END IF 
        IF (proc_x_max == MPI_PROC_NULL) THEN
          phi(nx+1,:) = phi(nx,:) 
          phi(nx+2,:) = phi(nx,:)
        END IF 

      END SUBROUTINE phi_mpi

  END SUBROUTINE potential_field

END MODULE initial_conditions
