!*******************************************************************
! All the real ghost cell values are controlled by these routines.
!*******************************************************************

MODULE boundary

  USE shared_data
  USE mpiboundary

  IMPLICIT NONE

CONTAINS

  SUBROUTINE set_boundary_conditions

    LOGICAL :: first_call = .TRUE.
    LOGICAL :: second_call = .FALSE.

    IF (second_call) THEN
      IF (xbc_max == BC_OPEN) THEN
        bx(nx+2,:) = bx(nx+1,:)
        by(nx+2,:) = by(nx+1,:)
        bz(nx+2,:) = bz(nx+1,:)
      END IF
      IF (xbc_min == BC_OPEN) THEN
        bx(-2,:) = bx(-1,:)
        by(-1,:) = by(0,:)
        bz(-1,:) = bz(0,:)
      END IF
      IF (ybc_max == BC_OPEN) THEN
        bx(:,ny+2) = bx(:,ny+1)
        by(:,ny+2) = by(:,ny+1)
        bz(:,ny+2) = bz(:,ny+1)
      END IF
      IF (ybc_min == BC_OPEN) THEN
        bx(:,-1) = bx(:,0)
        by(:,-2) = by(:,-1)
        bz(:,-1) = bz(:,0)
      END IF

      second_call = .FALSE.
    END IF

    IF (first_call) THEN
      any_open = .FALSE.
      IF ((xbc_max == BC_OPEN) .OR. (xbc_min == BC_OPEN) &
          .OR. (ybc_max == BC_OPEN) .OR. (ybc_min == BC_OPEN)) any_open = .TRUE.
      first_call = .FALSE.
      second_call = .TRUE.
    END IF

  END SUBROUTINE set_boundary_conditions



  SUBROUTINE boundary_conditions

    CALL bfield_bcs
    CALL energy_bcs
    CALL density_bcs
    CALL velocity_bcs

    CALL damp_boundaries

  END SUBROUTINE boundary_conditions



  SUBROUTINE damp_boundaries
    ! These are not generic, always work damping options.
    ! Users should change the damping scheme for each problem
    REAL(num) :: a, d

    IF (damping) THEN
      ! proc_x_max boundary
      IF (proc_x_max == MPI_PROC_NULL) THEN
        d = 0.7_num * x_end
        DO iy = -1, ny + 1
          DO ix = -1, nx + 1
            IF (xb(ix) > d) THEN
              a = dt * (xb(ix) - d) / (x_end - d)
              vx(ix, iy) = vx(ix, iy) / (1.0_num + a)
              vy(ix, iy) = vy(ix, iy) / (1.0_num + a)
              vz(ix, iy) = vz(ix, iy) / (1.0_num + a)
            END IF
          END DO
        END DO
      END IF

      ! proc_x_min boundary
      IF (proc_x_min == MPI_PROC_NULL) THEN
        d = 0.7_num * x_start
        DO iy = -1, ny + 1
          DO ix = -1, nx + 1
            IF (xb(ix) < d) THEN
              a = dt * (xb(ix) - d) / (x_start - d)
              vx(ix, iy) = vx(ix, iy) / (1.0_num + a)
              vy(ix, iy) = vy(ix, iy) / (1.0_num + a)
              vz(ix, iy) = vz(ix, iy) / (1.0_num + a)
            END IF
          END DO
        END DO
      END IF

      ! top boundary
      IF (proc_y_max == MPI_PROC_NULL) THEN
        d = 0.7_num * y_end
        DO iy = -1, ny + 1
          DO ix = -1, nx + 1
            IF (yb(iy) > d) THEN
              a = dt * (yb(iy) - d) / (y_end - d)
              vx(ix, iy) = vx(ix, iy) / (1.0_num + a)
              vy(ix, iy) = vy(ix, iy) / (1.0_num + a)
              vz(ix, iy) = vz(ix, iy) / (1.0_num + a)
            END IF
          END DO
        END DO
      END IF

      ! bottom boundary
      IF(proc_y_min == MPI_PROC_NULL) THEN
      d = 0.7_num * y_start
      DO iy = -1, ny + 1
        DO ix = -1, nx + 1
          IF (yb(iy) < d) THEN
            a = dt * (yb(iy) - d) / (y_start - d)
            vx(ix, iy) = vx(ix, iy) / (1.0_num + a)
            vy(ix, iy) = vy(ix, iy) / (1.0_num + a)
            vz(ix, iy) = vz(ix, iy) / (1.0_num + a)
          END IF
        END DO
      END DO
    END IF

    END IF

  END SUBROUTINE damp_boundaries



  SUBROUTINE bfield_bcs

    CALL bfield_MPI

    IF (proc_x_max == MPI_PROC_NULL .AND. xbc_max == BC_OTHER) THEN
      bx(nx+1, :) = bx(nx-1, :)
      bx(nx+2, :) = bx(nx-2, :)
      by(nx+1, :) = by(nx  , :)
      by(nx+2, :) = by(nx-1, :)
      bz(nx+1, :) = bz(nx  , :)
      bz(nx+2, :) = bz(nx-1, :)
    END IF
    IF (proc_x_min == MPI_PROC_NULL .AND. xbc_min == BC_OTHER) THEN
      bx(-1, :) = bx(1, :)
      bx(-2, :) = bx(2, :)
      by( 0, :) = by(1, :)
      by(-1, :) = by(2, :)
      bz( 0, :) = bz(1, :)
      bz(-1, :) = bz(2, :)
    END IF

    IF (proc_y_max == MPI_PROC_NULL .AND. ybc_max == BC_OTHER) THEN
      bx(:, ny+1) = bx(:, ny  )
      bx(:, ny+2) = bx(:, ny-1)
      by(:, ny+1) = by(:, ny-1)
      by(:, ny+2) = by(:, ny-2)
      bz(:, ny+1) = bz(:, ny  )
      bz(:, ny+2) = bz(:, ny-1)
    END IF
    IF (proc_y_min == MPI_PROC_NULL .AND. ybc_min == BC_OTHER) THEN
      bx(:,  0) = bx(:, 1)
      bx(:, -1) = bx(:, 2)
      by(:, -1) = by(:, 1)
      by(:, -2) = by(:, 2)
      bz(:,  0) = bz(:, 1)
      bz(:, -1) = bz(:, 2)
    END IF

  END SUBROUTINE bfield_bcs



  SUBROUTINE bz_bcs

    CALL bz_MPI

    IF (proc_x_max == MPI_PROC_NULL .AND. xbc_max == BC_OTHER) THEN
      bz(nx+1, :) = bz(nx  , :)
      bz(nx+2, :) = bz(nx-1, :)
    END IF
    IF (proc_x_min == MPI_PROC_NULL .AND. xbc_min == BC_OTHER) THEN
      bz( 0, :) = bz(1, :)
      bz(-1, :) = bz(2, :)
    END IF
    IF (proc_y_max == MPI_PROC_NULL .AND. ybc_max == BC_OTHER) THEN
      bz(:, ny+1) = bz(:, ny  )
      bz(:, ny+2) = bz(:, ny-1)
    END IF
    IF (proc_y_min == MPI_PROC_NULL .AND. ybc_min == BC_OTHER) THEN
      bz(:,  0) = bz(:, 1)
      bz(:, -1) = bz(:, 2)
    END IF

  END SUBROUTINE bz_bcs



  SUBROUTINE energy_bcs

    CALL energy_MPI

    IF (proc_x_max == MPI_PROC_NULL .AND. xbc_max == BC_OTHER) THEN
      energy(nx+1, :) = energy(nx  , :)
      energy(nx+2, :) = energy(nx-1, :)
    END IF
    IF (proc_x_min == MPI_PROC_NULL .AND. xbc_min == BC_OTHER) THEN
      energy( 0, :) = energy(1, :)
      energy(-1, :) = energy(2, :)
    END IF

    IF (proc_y_max == MPI_PROC_NULL .AND. ybc_max == BC_OTHER) THEN
      energy(:, ny+1) = energy(:, ny  )
      energy(:, ny+2) = energy(:, ny-1)
    END IF
    IF (proc_y_min == MPI_PROC_NULL .AND. ybc_min == BC_OTHER) THEN
      energy(:,  0) = energy(:, 1)
      energy(:, -1) = energy(:, 2)
    END IF

  END SUBROUTINE energy_bcs



  SUBROUTINE velocity_bcs

    CALL velocity_MPI

    IF (proc_x_max == MPI_PROC_NULL) THEN
      IF (xbc_max == BC_OTHER) THEN
        vx(nx:nx+2, :) = 0.0_num
        vy(nx:nx+2, :) = 0.0_num
        vz(nx:nx+2, :) = 0.0_num
      END IF
    END IF
    IF (proc_x_min == MPI_PROC_NULL) THEN
      IF (xbc_min == BC_OTHER) THEN
        vx(-2:0, :) = 0.0_num
        vy(-2:0, :) = 0.0_num
        vz(-2:0, :) = 0.0_num
      END IF
    END IF

    IF (proc_y_max == MPI_PROC_NULL) THEN
      IF (ybc_max == BC_OTHER) THEN
        vx(:, ny:ny+2) = 0.0_num
        vy(:, ny:ny+2) = 0.0_num
        vz(:, ny:ny+2) = 0.0_num
      END IF
    END IF
    IF (proc_y_min == MPI_PROC_NULL) THEN
      IF (ybc_min == BC_OTHER) THEN
        vx(:, -2:0) = 0.0_num
        vy(:, -2:0) = 0.0_num
        vz(:, -2:0) = 0.0_num
      END IF
    END IF

  END SUBROUTINE velocity_bcs



  SUBROUTINE remap_v_bcs

    CALL remap_v_MPI

    IF (proc_x_max == MPI_PROC_NULL) THEN
      IF (xbc_max == BC_OTHER) THEN
        vx1(nx:nx+2, :) = 0.0_num
        vy1(nx:nx+2, :) = 0.0_num
        vz1(nx:nx+2, :) = 0.0_num
      END IF
    END IF
    IF (proc_x_min == MPI_PROC_NULL) THEN
      IF (xbc_min == BC_OTHER) THEN
        vx1(-2:0, :) = 0.0_num
        vy1(-2:0, :) = 0.0_num
        vz1(-2:0, :) = 0.0_num
      END IF
    END IF

    IF (proc_y_max == MPI_PROC_NULL) THEN
      IF (ybc_max == BC_OTHER) THEN
        vx1(:, ny:ny+2) = 0.0_num
        vy1(:, ny:ny+2) = 0.0_num
        vz1(:, ny:ny+2) = 0.0_num
      END IF
    END IF
    IF (proc_y_min == MPI_PROC_NULL) THEN
      IF (ybc_min == BC_OTHER) THEN
        vx1(:, -2:0) = 0.0_num
        vy1(:, -2:0) = 0.0_num
        vz1(:, -2:0) = 0.0_num
      END IF
    END IF

  END SUBROUTINE remap_v_bcs



  SUBROUTINE density_bcs

    REAL(num) :: a, b

    CALL density_MPI

    IF (proc_x_max == MPI_PROC_NULL .AND. xbc_max == BC_OTHER) THEN
      rho(nx+1, :) = rho(nx  , :)
      rho(nx+2, :) = rho(nx-1, :)
    END IF
    IF (proc_x_min == MPI_PROC_NULL .AND. xbc_min == BC_OTHER) THEN
      rho( 0, :) = rho(1, :)
      rho(-1, :) = rho(2, :)
    END IF

    IF (proc_y_max == MPI_PROC_NULL .AND. ybc_max == BC_OTHER) THEN
      rho(:, ny+1) = rho(:, ny  )
      rho(:, ny+2) = rho(:, ny-1)
    END IF
    IF (proc_y_min == MPI_PROC_NULL .AND. ybc_min == BC_OTHER) THEN
      rho(:, 0) = rho(:, 1)
      rho(:, -1) = rho(:, 2)
    END IF


  END SUBROUTINE density_bcs

END MODULE boundary
