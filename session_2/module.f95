MODULE lennard_jones
  USE kinds, ONLY: wp => dp
  IMPLICIT NONE
CONTAINS


  SUBROUTINE distance_sub(position_3d, r)
    IMPLICIT NONE
    REAL(wp), INTENT(IN)  :: position_3d(:,:)
    REAL(wp), INTENT(OUT) :: r(:)
    INTEGER :: n, i, j, idx

    n = SIZE(position_3d, 1)
    idx = 1

    DO i = 1, n-1
      DO j = i+1, n
        r(idx) = SQRT( (position_3d(i,1) - position_3d(j,1)) ** 2 &
                      + (position_3d(i,2) - position_3d(j,2)) ** 2 &
                      + (position_3d(i,3) - position_3d(j,3)) ** 2 )
        idx = idx + 1
      END DO
    END DO
  END SUBROUTINE distance_sub


  SUBROUTINE distance_1d_sub(position_3d, direction, r_1d)
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: position_3d(:,:)
    INTEGER, INTENT(IN) :: direction

    REAL(wp), INTENT(OUT) :: r_1d(:)

    INTEGER :: n, i, j, idx

    n = SIZE(position_3d, 1)
    idx = 1

    DO i = 1, n-1
      DO j = i+1, n
        r_1d(idx) = position_3d(i, direction) - position_3d(j, direction)
        idx = idx + 1
      END DO
    END DO
  END SUBROUTINE distance_1d_sub


  SUBROUTINE lennard_sub(sigma, epsilon, r, result)
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: sigma, epsilon
    REAL(wp), INTENT(IN) :: r(:)
    REAL(wp), INTENT(OUT) :: result(:)
    INTEGER :: i
    REAL(wp) :: sr6, sr12

    DO i = 1, SIZE(r)
      sr6  = (sigma / r(i))**6
      sr12 = sr6 * sr6

      result(i) = 4.0_wp * epsilon * (sr12 - sr6)
    END DO
  END SUBROUTINE lennard_sub


  SUBROUTINE lennard_prime_sub(sigma, epsilon, r, result)
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: sigma, epsilon
    REAL(wp), INTENT(IN) :: r(:)
    REAL(wp), INTENT(OUT) :: result(:)
    INTEGER :: i
    REAL(wp) :: s6, s12, r7, r13

    DO i = 1, SIZE(r)
      s6  = sigma ** 6
      s12 = sigma ** 12
      r7  = r(i) ** 7
      r13 = r(i) ** 13

      result(i) = 4.0_wp * epsilon * ((-12.0_wp)*(s12/r13) + 6.0_wp*(s6/r7))
    END DO

  END SUBROUTINE lennard_prime_sub


  SUBROUTINE force_1d_sub(r_1d, r, lennard_prime_result, result)

    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: r_1d(:), r(:), lennard_prime_result(:)
    REAL(wp), INTENT(OUT) :: result(:)

    result = -(r_1d / r) * lennard_prime_result

  END SUBROUTINE force_1d_sub


  SUBROUTINE position_1d_sub(position_init_1d, taw, velocity_1d_result, &
                             m, force_1d_result_k, position_1d_result)
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: position_init_1d(:), velocity_1d_result(:), &
                            m(:), force_1d_result_k(:)
    REAL(wp), INTENT(IN) :: taw

    REAL(wp), INTENT(OUT) :: position_1d_result(:)

    position_1d_result = position_init_1d + taw * velocity_1d_result + &
                         0.5_wp * taw ** 2 * (force_1d_result_k / m)

  END SUBROUTINE position_1d_sub


  SUBROUTINE velocity_1d_sub(velocity_init_1d, taw, m, force_1d_result_k, &
                             force_1d_result_kp1, velocity_1d_result)
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: velocity_init_1d(:), m(:), force_1d_result_k(:), &
                            force_1d_result_kp1(:)
    REAL(wp), INTENT(IN) :: taw

    REAL(wp), INTENT(OUT) :: velocity_1d_result(:)

    velocity_1d_result = velocity_init_1d + 0.5_wp * taw * (force_1d_result_k + &
                         force_1d_result_kp1) / m

  END SUBROUTINE velocity_1d_sub


  SUBROUTINE kin_energy_sub(m, velocity_x, velocity_y, velocity_z, kin_energy_result)

    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: m(:), velocity_x(:), velocity_y(:), velocity_z(:)
    INTEGER :: i, n
    REAL(wp), INTENT(OUT) :: kin_energy_result(:)

    n = SIZE(m)
    DO i = 1, n

      kin_energy_result(i) = 0.5_wp * m(i) * (velocity_x(i)**2 + velocity_y(i)**2 + &
                                              velocity_z(i)**2)

    END DO


  END SUBROUTINE kin_energy_sub


  SUBROUTINE diff_pot_energy_sub(force_k_x, force_k_y, force_k_z, x_kp1, &
                                 x_k, y_kp1, y_k, z_kp1, z_k, diff_pot_energy_result)

    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: force_k_x(:), force_k_y(:), force_k_z(:), x_kp1(:), x_k(:), &
                            y_kp1(:), y_k(:), z_kp1(:), z_k(:)
    INTEGER :: i, n
    REAL(wp), INTENT(OUT) :: diff_pot_energy_result(:)

    n = SIZE(x_k)
    DO i = 2, n

      diff_pot_energy_result(i) = (force_k_x(i) * (x_kp1(i) - x_k(i))) + &
                                  (force_k_y(i) * (y_kp1(i) - y_k(i))) + &
                                  (force_k_z(i)) * (z_kp1(i) - z_k(i))

    END DO

  END SUBROUTINE diff_pot_energy_sub


END MODULE lennard_jones
