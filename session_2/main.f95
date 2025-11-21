PROGRAM main
  USE kinds, ONLY: wp => dp
  USE lennard_jones
  IMPLICIT NONE

  INTEGER :: n, i, j, k, nk, n_pairs, pair_idx
  REAL(wp) :: tau, epsilon, sigma, time_k
  REAL(wp), ALLOCATABLE :: m(:)
  REAL(wp), ALLOCATABLE :: position_3d(:,:), velocity_3d(:,:)
  REAL(wp), ALLOCATABLE :: r(:), r_z(:), lennard_prime(:)
  REAL(wp), ALLOCATABLE :: force_z(:), new_pos_z(:), new_vel_z(:)
  REAL(wp) :: fz_pair
  CHARACTER(len=*), PARAMETER :: atom_symbol = "Ne"

  OPEN(11, FILE="./Ne_data.txt", STATUS="old", ACTION="read")

  READ(11,*) nk, tau
  READ(11,*) sigma, epsilon
  READ(11,*) n

  ALLOCATE(position_3d(n,3))
  ALLOCATE(velocity_3d(n,3))
  ALLOCATE(m(n))
  n_pairs = n * (n - 1) / 2
  ALLOCATE(r(n_pairs))
  ALLOCATE(r_z(n_pairs))
  ALLOCATE(lennard_prime(n_pairs))
  ALLOCATE(force_z(n))
  ALLOCATE(new_pos_z(n))
  ALLOCATE(new_vel_z(n))

  DO i = 1, n
     READ(11,*) m(i), position_3d(i,1), position_3d(i,2), position_3d(i,3), &
                      velocity_3d(i,1), velocity_3d(i,2), velocity_3d(i,3)
  END DO

  CLOSE(11)

  OPEN(12, FILE="./ne_positions.xyz", STATUS="replace", ACTION="write")

  DO k = 1, nk
     time_k = (k - 1) * tau

     WRITE(12,'(I0)') n
     WRITE(12,'("time = ",F18.6)') time_k

     DO i = 1, n
        WRITE(12,'(A2,3(1X,F18.10))') atom_symbol, position_3d(i,1), &
                                       position_3d(i,2), position_3d(i,3)
     END DO

     CALL distance_sub(position_3d, r)
     CALL distance_1d_sub(position_3d, 3, r_z)
     CALL lennard_prime_sub(sigma, epsilon, r, lennard_prime)

     force_z = 0.0_wp
     pair_idx = 1
     DO i = 1, n-1
        DO j = i+1, n
           fz_pair = -lennard_prime(pair_idx) * (position_3d(i,3) - position_3d(j,3)) / r(pair_idx)
           force_z(i) = force_z(i) + fz_pair
           force_z(j) = force_z(j) - fz_pair
           pair_idx = pair_idx + 1
        END DO
     END DO

     CALL position_1d_sub(position_3d(:,3), tau, velocity_3d(:,3), m, force_z, new_pos_z)
     CALL velocity_1d_sub(velocity_3d(:,3), tau, m, force_z, force_z, new_vel_z)

     position_3d(:,3) = new_pos_z
     velocity_3d(:,3) = new_vel_z
  END DO

  CLOSE(12)

END PROGRAM main
