PROGRAM verlet
  IMPLICIT NONE ! inline comment
  ! declare your variable here
  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND (p=13, r=300)
  INTEGER :: k
  REAL (KIND=wp) :: taw, mass, velocity_init
  REAL (KIND=wp), DIMENSION(:), ALLOCATABLE :: position_3d
  REAL (KIND=wp), DIMENSION(:), ALLOCATABLE :: force_3d
  REAL (KIND=wp), DIMENSION(:), ALLOCATABLE :: velocity_3d

  ! write your instructions here
  ALLOCATE (position_3d(3))
  ALLOCATE (force_3d(3))
  ALLOCATE (velocity_3d(3))

  force_3d = (/ 0.0_wp, 0.1_wp, 0.0_wp/)
  taw = 0.2_wp
  mass = 1.0_wp
  velocity_init = 0
  velocity_3d = velocity_init

  DO k = 1, 600
    velocity_3d =  velocity_3d + taw / mass * force_3d
    print *, "velocity 3d in ", k, "step is: ", velocity_3d
  END DO
  DEALLOCATE (position_3d, velocity_3d, force_3d)
END PROGRAM verlet