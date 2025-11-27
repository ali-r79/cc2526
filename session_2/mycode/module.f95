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
    REAL(wp), INTENT(IN) :: position_3d(:,:)
    INTEGER, INTENT(IN) :: direction

    REAL(wp), INTENT(OUT) :: r_1(:)

    INTEGER :: n, i, j, idx

    n = SIZE(position_3d, 1)

    DO i = 1, n-1
      DO j = i+1, n
        r_1(idx) = ABS(position_3d(i, direction) - position_3d(i, direction))
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
      r7 = r ** 7
      r13 = r ** 13

      result(i) = 4 * epsilon * ((-12)*(s12/r13) + 6*(s6/r7))
    END DO

  END SUBROUTINE lennard_prime_sub


  SUBROUTINE force_1d_sub(r_1d, r, lennard_prime_result, result)

  IMPLICIT NONE
  REAL(wp), INTENT(IN) :: r_1d(:), r(:), lennard_prime_result(:)
  REAL(wp), INTENT(OUT) :: result(:)

  result = -SUM(r_1d / r_1 * lennard_prime_result)

  END SUBROUTINE force_1d_sub


  SUBROUTINE position_1d_sub(position_init_1d, taw, velocity_1d_result, &
                             m, force_1d_result_k, position_1d_result)
    REAL(wp), INTENT(IN) :: position_init_1d(:), velocity_1d_result_k(:), &
                            velocity_1d_result_kp1(:), force_1d_result_k(:), m(:)
    REAL(wp), INTENT(IN) :: taw

    REAL(wp), INTENT(OUT) :: position_1d_result(:)

    REAL(wp), ALLOCATABLE :: position_1d_k(:), position_1d_kp1(:)
    INTEGER :: i, n

    n = SIZE(position_init_1d)
    ALLOCATE(position_1d_k(n))
    ALLOCATE(position_1d_kp1(n))

    position_1d_k = position_init_1d
    DO i = 1, n
      position_1d_kp1 = position_1d_k + taw * velocity_1d_result_k + taw ** 2 * (force_1d_result_k/(2*m))
      position_1d_k = position_1d_kp1
    END DO

  END SUBROUTINE position_1d_sub


  SUBROUTINE velocity_1d_sub(velocity_init_1d, taw, m, force_1d_result_k, force_1d_result_kp1, velocity_1d_result)
    REAL(wp), INTENT(IN) :: velocity_init_1d(:), velocity_1d_result_k(:), &
                            velocity_1d_result_kp1(:), force_1d_result_k(:), m(:)
    REAL(wp), INTENT(IN) :: taw

    REAL(wp), INTENT(OUT) :: velocity_1d_result(:)

    REAL(wp), ALLOCATABLE :: velocity_1d_k(:), velocity_1d_kp1(:)
    INTEGER :: i, n

    n = SIZE(velocity_init_1d)
    ALLOCATE(velocity_1d_k(n))
    ALLOCATE(velocity_1d_kp1(n))

    velocity_1d_k = velocity_init_1d
    DO i = 1, n
      velocity_1d_kp1 = velocity_1d_k + (taw/2*m) * (force_1d_result_k + force_1d_result_kp1)
      velocity_1d_k = velocity_1d_kp1
    END DO
    velocity_1d_result = velocity_1d_k

  END SUBROUTINE velocity_1d_sub


END MODULE lennard_jones
