PROGRAM main
  USE kinds, ONLY: wp => dp
  USE lennard_jones
  IMPLICIT NONE

  INTEGER :: n, i, n_pairs, nk
  REAL(wp) :: tau, epsilon, sigma

  REAL(wp), ALLOCATABLE :: m(:)
  REAL(wp), ALLOCATABLE :: position_3d(:,:), velocity_3d(:,:), &
                           position_init(:,:), velocity_init(:,:)
  REAL(wp), ALLOCATABLE :: r(:), r_z(:), lennard_prime(:), lennard(:)
  REAL(wp), ALLOCATABLE :: force_z(:), new_pos_z(:), new_vel_z(:)
  REAL(wp), ALLOCATABLE :: kin_energy_vec(:), pot_energy_vec(:)
  REAL(wp) :: tot_kin_energy, tot_pot_energy, tot_energy
  REAL(wp) :: fz_pair
  CHARACTER(len=*), PARAMETER :: atom_symbol = "Ne"

  OPEN(11, FILE="./Ne_data.txt", STATUS="old", ACTION="read")

  READ(11,*) nk, tau
  READ(11,*) sigma, epsilon
  READ(11,*) n

  ALLOCATE(position_3d(n,3))
  ALLOCATE(velocity_3d(n,3))
  ALLOCATE(position_init(n,3))
  ALLOCATE(velocity_init(n,3))
  ALLOCATE(m(n))
  n_pairs = n * (n - 1) / 2
  ALLOCATE(r(n_pairs))
  ALLOCATE(r_z(n_pairs))
  ALLOCATE(lennard_prime(n_pairs))
  ALLOCATE(lennard(n_pairs))
  ALLOCATE(force_z(n))
  ALLOCATE(new_pos_z(n))
  ALLOCATE(new_vel_z(n))
  ALLOCATE(kin_energy_vec(n))
  ALLOCATE(pot_energy_vec(n))

  DO i = 1, n
     READ(11,*) m(i), position_3d(i,1), position_3d(i,2), position_3d(i,3), &
                      velocity_3d(i,1), velocity_3d(i,2), velocity_3d(i,3)
  END DO

  position_init = position_3d
  velocity_init = velocity_3d

  CLOSE(11)

  CALL run_simulation(tau, nk, position_init, velocity_init, "./ne_positions_tau_1.xyz", &
                      "./energy_tau_1.txt")

  CALL run_simulation(0.01_wp, nk, position_init, velocity_init, "./ne_positions_tau_0p01.xyz", &
                      "./energy_tau_0p01.txt")

CONTAINS

  SUBROUTINE run_simulation(tau_step, n_steps, position_start, velocity_start, &
                            xyz_filename, energy_filename)
    USE kinds, ONLY: wp => dp
    USE lennard_jones
    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: tau_step
    INTEGER,  INTENT(IN) :: n_steps
    REAL(wp), INTENT(IN) :: position_start(:,:), velocity_start(:,:)
    CHARACTER(len=*), INTENT(IN) :: xyz_filename, energy_filename

    INTEGER :: k, i, j, pair_idx
    INTEGER :: xyz_unit, energy_unit
    REAL(wp) :: time_k, fz_pair, pot_contribution

    position_3d = position_start
    velocity_3d = velocity_start

    OPEN(NEWUNIT=xyz_unit, FILE=xyz_filename, STATUS="replace", ACTION="write")
    OPEN(NEWUNIT=energy_unit, FILE=energy_filename, STATUS="replace", ACTION="write")

    DO k = 1, n_steps
       time_k = (k - 1) * tau_step

       WRITE(xyz_unit,'(I0)') n
       WRITE(xyz_unit,'("time = ",F18.2)') time_k

       DO i = 1, n
          WRITE(xyz_unit,'(A2,3(1X,F18.10))') atom_symbol, position_3d(i,1), &
                                             position_3d(i,2), position_3d(i,3)
       END DO

       CALL distance_sub(position_3d, r)
       CALL distance_1d_sub(position_3d, 3, r_z)
       CALL lennard_prime_sub(sigma, epsilon, r, lennard_prime)
       CALL lennard_sub(sigma, epsilon, r, lennard)

       force_z = 0.0_wp
       pot_energy_vec = 0.0_wp
       pair_idx = 1
       DO i = 1, n-1
          DO j = i+1, n
             fz_pair = -lennard_prime(pair_idx) * (position_3d(i,3) - position_3d(j,3)) &
                       / r(pair_idx)
             force_z(i) = force_z(i) + fz_pair
             force_z(j) = force_z(j) - fz_pair

             pot_contribution = lennard(pair_idx)
             pot_energy_vec(i) = pot_energy_vec(i) + 0.5_wp * pot_contribution
             pot_energy_vec(j) = pot_energy_vec(j) + 0.5_wp * pot_contribution

             pair_idx = pair_idx + 1
          END DO
       END DO

       CALL kin_energy_sub(m, velocity_3d(:,1), velocity_3d(:,2), velocity_3d(:,3), kin_energy_vec)
       tot_kin_energy = SUM(kin_energy_vec)
       tot_pot_energy = SUM(pot_energy_vec)

       IF (MOD(k, 100) == 0) THEN
          WRITE(energy_unit,'("step=",I8," time=",F12.2)') k, time_k
          WRITE(energy_unit,'("kin:",*(1X,ES16.8)," total=",ES16.8)') kin_energy_vec, tot_kin_energy
          WRITE(energy_unit,'("pot:",*(1X,ES16.8)," total=",ES16.8)') pot_energy_vec, tot_pot_energy
         !  tot_energy = tot_kin_energy + tot_pot_energy
         !  WRITE(energy_unit, '("tot:",*(1X,ES16.8)," total=",ES16.8)') tot_energy
       END IF

       CALL position_1d_sub(position_3d(:,3), tau_step, velocity_3d(:,3), m, force_z, new_pos_z)
       CALL velocity_1d_sub(velocity_3d(:,3), tau_step, m, force_z, force_z, new_vel_z)

       position_3d(:,3) = new_pos_z
       velocity_3d(:,3) = new_vel_z
    END DO

    CLOSE(xyz_unit)
    CLOSE(energy_unit)

  END SUBROUTINE run_simulation

END PROGRAM main
