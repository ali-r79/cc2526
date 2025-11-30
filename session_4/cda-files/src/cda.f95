PROGRAM cda
  USE kinds, ONLY: wp => dp
  USE cubes
  IMPLICIT NONE

  TYPE(cube) :: cube_ab, cube_a, cube_b, cube_ref, cube_drho
  REAL(KIND=wp), ALLOCATABLE :: cdz(:)
  INTEGER :: i

  CALL cube_get(cube_ab, '../test/CuCO+/ab.cube')
  CALL cube_get(cube_a,  '../test/CuCO+/a.cube')
  CALL cube_get(cube_b,  '../test/CuCO+/b.cube')

  cube_ref = cube_a + cube_b
  cube_drho = cube_ab - cube_ref

  cdz = cube_cdz(cube_drho)

  WRITE (*,*) 'Integral rho: ', cube_int(cube_drho)
  WRITE (*,*) 'Charge change along z:'
  DO i = 1, SIZE(cdz)
    WRITE (*,'(I6,1X,F20.10)') i, cdz(i)
  END DO

  CALL cube_write('charge_density.dat', '# index  delta_rho', cube_drho)
  CALL cube_write('electron_density.dat', '# index  electron_density', cube_ab)
  CALL write_vector('charge_z.dat', '# z_index  cumulative_delta_rho', cdz)

  CALL cube_del(cube_ab)
  CALL cube_del(cube_a)
  CALL cube_del(cube_b)
  CALL cube_del(cube_ref)
  CALL cube_del(cube_drho)

CONTAINS

  SUBROUTINE write_vector(filename, header, data)
    CHARACTER(LEN=*), INTENT(IN) :: filename, header
    REAL(KIND=wp), INTENT(IN) :: data(:)
    INTEGER :: i, unit

    OPEN(NEWUNIT=unit, FILE=filename, STATUS='replace', ACTION='write')
    WRITE(unit,'(A)') header
    DO i = 1, SIZE(data)
      WRITE(unit,'(I10,1X,ES20.10)') i, data(i)
    END DO
    CLOSE(unit)
  END SUBROUTINE write_vector

END PROGRAM cda
