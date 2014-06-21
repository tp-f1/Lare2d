!******************************************************************************
! Routines to set up the MPI routines and allocate dynamic arrays
!******************************************************************************

MODULE mpi_routines

  USE shared_data
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: mpi_initialise, mpi_close

  REAL(dbl) :: start_time, end_time

CONTAINS

  !****************************************************************************
  ! Start up the MPI layer, allocate arrays and set up MPI types
  !****************************************************************************

  SUBROUTINE mpi_initialise

    LOGICAL :: periods(c_ndims), reorder
    INTEGER :: starts(c_ndims), sizes(c_ndims), subsizes(c_ndims)
    INTEGER :: dims(c_ndims)
    INTEGER :: nx0, ny0
    INTEGER :: nxp, nyp
    INTEGER :: cx, cy

    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, errcode)

    dims = (/ nprocy, nprocx /)

    IF (PRODUCT(MAX(dims, 1)) > nproc) THEN
      dims = 0
      IF (rank == 0) THEN
        PRINT*, 'Too many processors requested in override.'
        PRINT*, 'Reverting to automatic decomposition.'
        PRINT*, '******************************************'
        PRINT*
      END IF
    END IF

    CALL MPI_DIMS_CREATE(nproc, c_ndims, dims, errcode)

    nprocx = dims(c_ndims  )
    nprocy = dims(c_ndims-1)

    periods = .TRUE.
    reorder = .TRUE.

    IF (xbc_min == BC_OTHER) periods(c_ndims  ) = .FALSE.
    IF (xbc_max == BC_OTHER) periods(c_ndims  ) = .FALSE.
    IF (ybc_min == BC_OTHER) periods(c_ndims-1) = .FALSE.
    IF (ybc_max == BC_OTHER) periods(c_ndims-1) = .FALSE.

    IF (xbc_min == BC_OPEN) periods(c_ndims  ) = .FALSE.
    IF (xbc_max == BC_OPEN) periods(c_ndims  ) = .FALSE.
    IF (ybc_min == BC_OPEN) periods(c_ndims-1) = .FALSE.
    IF (ybc_max == BC_OPEN) periods(c_ndims-1) = .FALSE.

    CALL MPI_CART_CREATE(MPI_COMM_WORLD, c_ndims, dims, periods, reorder, &
        comm, errcode)

    CALL MPI_COMM_RANK(comm, rank, errcode)
    CALL MPI_CART_COORDS(comm, rank, c_ndims, coordinates, errcode)
    CALL MPI_CART_SHIFT(comm, c_ndims-1, 1, proc_x_min, proc_x_max, errcode)
    CALL MPI_CART_SHIFT(comm, c_ndims-2, 1, proc_y_min, proc_y_max, errcode)

    cx = coordinates(c_ndims  )
    cy = coordinates(c_ndims-1)

    ! Create the subarray for this problem: subtype decribes where this
    ! process's data fits into the global picture.

    nx0 = nx_global / nprocx
    ny0 = ny_global / nprocy
    nx = nx0
    ny = ny0

    ! If the number of gridpoints cannot be exactly subdivided then fix
    ! the first nxp processors have nx0 grid points
    ! The remaining processors have nx0+1 grid points
    IF (nx0 * nprocx /= nx_global) THEN
      nxp = (nx0 + 1) * nprocx - nx_global
      IF (cx >= nxp) nx = nx0 + 1
    ELSE
      nxp = nprocx
    END IF
    IF (ny0 * nprocy /= ny_global) THEN
      nyp = (ny0 + 1) * nprocy - ny_global
      IF (cy >= nyp) ny = ny0 + 1
    ELSE
      nyp = nprocy
    END IF

    ! Set up the starting point for my subgrid (assumes arrays start at 0)

    IF (cx < nxp) THEN
      starts(1) = cx * nx0
    ELSE
      starts(1) = nxp * nx0 + (cx - nxp) * (nx0 + 1)
    END IF
    IF (cy < nyp) THEN
      starts(2) = cy * ny0
    ELSE
      starts(2) = nyp * ny0 + (cy - nyp) * (ny0 + 1)
    END IF

    ! The grid sizes
    subsizes = (/ nx+1, ny+1 /)
    sizes = (/ nx_global+1, ny_global+1 /)

    CALL MPI_TYPE_CREATE_SUBARRAY(c_ndims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, mpireal, subtype, errcode)

    CALL MPI_TYPE_COMMIT(subtype, errcode)

    ALLOCATE(energy(-1:nx+2, -1:ny+2))
    ALLOCATE(p_visc(-1:nx+2, -1:ny+2))
    ALLOCATE(rho(-1:nx+2, -1:ny+2))
    ALLOCATE(vx (-2:nx+2, -2:ny+2))
    ALLOCATE(vy (-2:nx+2, -2:ny+2))
    ALLOCATE(vz (-2:nx+2, -2:ny+2))
    ALLOCATE(vx1(-2:nx+2, -2:ny+2))
    ALLOCATE(vy1(-2:nx+2, -2:ny+2))
    ALLOCATE(vz1(-2:nx+2, -2:ny+2))
    ALLOCATE(bx (-2:nx+2, -1:ny+2))
    ALLOCATE(by (-1:nx+2, -2:ny+2))
    ALLOCATE(bz (-1:nx+2, -1:ny+2))
    ALLOCATE(eta(-1:nx+2, -1:ny+2))
    IF (rke) ALLOCATE(delta_ke(-1:nx+2, -1:ny+2))
    ALLOCATE(lambda_i(0:nx, 0:ny))

    ! Shocked and resistive need to be larger to allow offset = 4 in shock_test
    ALLOCATE(cv(-1:nx+2, -1:ny+2), cv1(-1:nx+2, -1:ny+2))
    ALLOCATE(xc(-1:nx+2), xb(-2:nx+2), dxc(-1:nx+1), dxb(-1:nx+2))
    ALLOCATE(yc(-1:ny+2), yb(-2:ny+2), dyc(-1:ny+1), dyb(-1:ny+2))
    ALLOCATE(grav(-1:ny+2))
    ALLOCATE(jx_r(0:nx+1, 0:ny+1))
    ALLOCATE(jy_r(0:nx+1, 0:ny+1))
    ALLOCATE(jz_r(0:nx+1, 0:ny+1))

    ALLOCATE(xb_global(-2:nx_global+2))
    ALLOCATE(yb_global(-2:ny_global+2))

    IF (rank == 0) start_time = MPI_WTIME()

    p_visc = 0.0_num
    eta = 0.0_num

  END SUBROUTINE mpi_initialise



  !****************************************************************************
  ! Shutdown the MPI layer, deallocate arrays and set up timing info
  !****************************************************************************

  SUBROUTINE mpi_close

    INTEGER :: seconds, minutes, hours, total

    IF (rank == 0) THEN
      end_time = MPI_WTIME()
      total = INT(end_time - start_time)
      seconds = MOD(total, 60)
      minutes = MOD(total / 60, 60)
      hours = total / 3600
      WRITE(stat_unit,*)
      WRITE(stat_unit,'(''runtime = '', i4, ''h '', i2, ''m '', i2, ''s on '', &
          & i4, '' process elements.'')') hours, minutes, seconds, nproc
    END IF

    CALL MPI_BARRIER(comm, errcode)

    DEALLOCATE(rho, energy)
    DEALLOCATE(vx, vy, vz)
    DEALLOCATE(vx1, vy1, vz1)
    DEALLOCATE(bx, by, bz)
    DEALLOCATE(p_visc)
    DEALLOCATE(eta)
    DEALLOCATE(cv, cv1)
    DEALLOCATE(xc, xb, dxb, dxc)
    DEALLOCATE(yc, yb, dyb, dyc)
    DEALLOCATE(grav)
    DEALLOCATE(jx_r, jy_r, jz_r)
    DEALLOCATE(xb_global, yb_global)
    DEALLOCATE(lambda_i)

    IF (ALLOCATED(xi_n)) DEALLOCATE(xi_n)
    IF (ALLOCATED(delta_ke)) DEALLOCATE(delta_ke)
    IF (ALLOCATED(eta_perp)) DEALLOCATE(eta_perp)
    IF (ALLOCATED(parallel_current)) DEALLOCATE(parallel_current)
    IF (ALLOCATED(perp_current)) DEALLOCATE(perp_current)

  END SUBROUTINE mpi_close

END MODULE mpi_routines
