PROGRAM main
  USE kinds, ONLY: wp => dp
  USE lennard_jones, lennard_prime
  IMPLICIT NONE

  INTEGER :: n, i, nk
  REAL(wp) :: tau, epsilon, sigma
  REAL(wp), ALLOCATABLE :: m(:)
  REAL(wp), ALLOCATABLE :: position_3d(:,:), velocity_3d(:,:), force_3d(:,:)
  REAL(wp), ALLOCATABLE :: r(:), result_lennard(:)

  OPEN(11, FILE="./Ne_data.txt", STATUS="old", ACTION="read")

  READ(11,*) nk, tau
  READ(11,*) sigma, epsilon
  READ(11,*) n

  ALLOCATE(position_3d(n,3))
  ALLOCATE(velocity_3d(n,3))
  ALLOCATE(force_3d(n,3))
  ALLOCATE(m(n))
  ALLOCATE(r(n*(n-1)/2))
  ALLOCATE(result_lennard(n*(n-1)/2))

  DO i = 1, n
     READ(11,*) m(i), position_3d(i,1), position_3d(i,2), position_3d(i,3), &
                      velocity_3d(i,1), velocity_3d(i,2), velocity_3d(i,3)
  END DO

  CLOSE(11)

  CALL distance(position_3d, r)

  CALL lennard_func(sigma, epsilon, r, result_lennard)


  PRINT *, "Lennard-Jones potential values:"
  DO i = 1, SIZE(result_lennard)
     PRINT *, result_lennard(i)
  END DO

END PROGRAM main
