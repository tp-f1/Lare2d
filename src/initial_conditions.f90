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

    REAL(num) :: xi_v 
    REAL(num) :: w_tr = 0.4_num, T_ph = 10.0_num, T_cor = 1000.0_num
    REAL(num) :: y_ph = 0.0_num, y_cor = 15.0_num

    integer :: i
    integer :: j
    
    REAL(num), DIMENSION(:), ALLOCATABLE :: yc_global, dy 
    REAL(num), DIMENSION(:), ALLOCATABLE :: rho_y, temp_y, energy_y
    REAL(num), DIMENSION(:), ALLOCATABLE :: mu 
    
    ALLOCATE(yc_global(-1:ny+1))
    ALLOCATE(dy(-1:ny+1))
    ALLOCATE(rho_y(-1:ny+2))
    ALLOCATE(temp_y(-1:ny+2))
    ALLOCATE(energy_y(-1:ny+2))
    ALLOCATE(mu(-1:ny+2))
    

    DO i = -1, ny + 1
        yc_global(i) = 0.5_num * (yb_global(i-1) + yb_global(i))
        dy(i) = yb_global(i) - yb_global(i-1)
    END DO 

    ! Below are all the variables which must be defined and their sizes
    
    vx(:, :) = 0.0_num
    vy(:, :) = 0.0_num
    vz(:, :) = 0.0_num

    bx(:, :) = 0.0_num
    by(:, :) = 0.0_num
    bz(:, :) = 0.0_num

    rho(:, :) = 1.0_num
    energy(:, :) = 1.0_num

    grav(:) = 15.0_num
    DO i = -1, ny + 2
        IF (yb_global(i) > y_ph) THEN
           grav(i) = grav(i) * (y_cor - yb_global(i)) / (y_cor - y_ph)  
        END IF 

        IF (yb_global(i) > y_cor) THEN
           grav(i) = 0.0_num
        END IF
    END DO
 
 
    ! If probe points needed add them here
    CALL add_probe(0.0_num, 0.0_num)

    DEALLOCATE(yc_global)
    DEALLOCATE(dy)
    DEALLOCATE(rho_y)
    DEALLOCATE(temp_y)
    DEALLOCATE(energy_y)
    DEALLOCATE(mu)


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


  SUBROUTINE potential_field()

      REAL(num), DIMENSION(:,:), ALLOCATABLE :: phi
      REAL(num) :: w, errmax, error, residual, fractional_error
      REAL(num) :: by_min, by_min_local
      REAL(num) :: by_max, by_max_local
      INTEGER :: loop, x1, y1, redblack
      LOGICAL :: converged

      ALLOCATE(phi(-1:nx+2,-1:ny+2))
      phi = 0.0_num
      CALL phi_mpi

      converged = .FALSE.
      w = 2.0_num / (1.0_num + SIN(pi / REAL(nx_global,num)))
      fractional_error = 1.e-8_num

      ! Iterate to get phi^{n+1} by SOR Gauss-Seidel
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

      DO ix = 0, nx
        DO iy = 1, ny
          bx(ix,iy) = -(phi(ix+1,iy)-phi(ix,iy))/dxc(ix)
        END DO
      END DO

      DO ix = 1, nx
        DO iy = 0, ny
          by(ix,iy) = -(phi(ix,iy+1)-phi(ix,iy))/dyc(iy)
        END DO
      END DO

      CALL bfield_bcs

      !Only incoming flux on lower boundary
      by_min_local = MAXVAL(by)
      IF (proc_y_min == MPI_PROC_NULL) THEN
        by_min_local = MINVAL(by(1:nx,0))
      END IF
      CALL MPI_ALLREDUCE(by_min_local, by_min, 1, mpireal, MPI_MIN, comm, errcode)
      by = by - MIN(by_min, 0.0_num)

      !Find maximum By on lower boundary
      by_max_local = MINVAL(by)
      IF (proc_y_min == MPI_PROC_NULL) THEN
        by_max_local = MAXVAL(by(1:nx,0))
      END IF
      CALL MPI_ALLREDUCE(by_max_local, by_max, 1, mpireal, MPI_MAX, comm, errcode) 

      !Scale the field to maximum of 1 in normalised units  
      by = by / by_max 
      bx = bx / by_max 

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

        !Unipolar flux
        local_flux = 0.0_num
        IF (proc_y_min == MPI_PROC_NULL) THEN
          phi(1:nx,0) = phi(1:nx,1) + dyc(1) * EXP(-xc(1:nx)**2) 
          local_flux = SUM(dxb(1:nx) * EXP(-xc(1:nx)**2))
          phi(1:nx,-1) = phi(1:nx,0)
        END IF
        CALL MPI_ALLREDUCE(local_flux, total_flux, 1, mpireal, MPI_SUM, comm, errcode)
        IF (proc_y_min == MPI_PROC_NULL) THEN
          phi(1:nx,0) = phi(1:nx,0) - dyc(1) * total_flux / length_x
          phi(1:nx,-1) = phi(1:nx,0)
        END IF        

        IF (proc_y_max == MPI_PROC_NULL) THEN
          phi(:,ny+1) = 0.0_num
          phi(:,ny+2) = 0.0_num
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
