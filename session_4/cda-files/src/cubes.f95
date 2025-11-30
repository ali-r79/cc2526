MODULE cubes
  USE kinds, ONLY: wp => dp
  IMPLICIT NONE

  TYPE, PUBLIC :: cube
    PRIVATE
    CHARACTER (LEN=72) :: str1
    CHARACTER (LEN=72) :: str2
    REAL (KIND=wp) :: xmin, ymin, zmin, dx, dy, dz
    INTEGER :: nx, ny, nz, natom, npts
    INTEGER, DIMENSION(:), POINTER :: zahl
    REAL (KIND=wp), DIMENSION(:), POINTER :: chrg, x, y, z
    REAL (KIND=wp), DIMENSION(:), POINTER :: array

    ! fill
  END TYPE cube

  PUBLIC :: cube_get, &
            cube_add, &
            cube_sub, &
            cube_int, &
            cube_cdz, &
            cube_write, &
            cube_del

  INTERFACE OPERATOR(+)
    MODULE PROCEDURE cube_add
  END INTERFACE
  INTERFACE OPERATOR(-)
    MODULE PROCEDURE cube_sub
  END INTERFACE


  CONTAINS


  SUBROUTINE cube_get (mycube, infile)
    CHARACTER(LEN=*), INTENT(IN) :: infile
    TYPE (cube), INTENT(OUT) :: mycube
    INTEGER :: i
    REAL(KIND=wp) :: vx, vy, vz
    OPEN (UNIT=11, FILE=infile, STATUS="old", ACTION="read")

    READ (UNIT=11, FMT=*) mycube%str1
    READ (UNIT=11, FMT=*) mycube%str2
    READ (UNIT=11, FMT=*) mycube%natom, mycube%xmin, mycube%ymin, mycube%zmin
    READ (UNIT=11, FMT=*) mycube%nx, vx, vy, vz
    mycube%dx = vx
    READ (UNIT=11, FMT=*) mycube%ny, vx, vy, vz
    mycube%dy = vy
    READ (UNIT=11, FMT=*) mycube%nz, vx, vy, vz
    mycube%dz = vz

    mycube%npts = mycube%nx * mycube%ny * mycube%nz

    ALLOCATE(mycube%zahl(mycube%natom))
    ALLOCATE(mycube%chrg(mycube%natom))
    ALLOCATE(mycube%x(mycube%natom))
    ALLOCATE(mycube%y(mycube%natom))
    ALLOCATE(mycube%z(mycube%natom))

    ALLOCATE(mycube%array(mycube%npts))

    DO i = 1, mycube%natom
      READ (UNIT=11, FMT=*) mycube%zahl(i), mycube%chrg(i), mycube%x(i), mycube%y(i), mycube%z(i)
    END DO

    READ (UNIT=11, FMT=*) (mycube%array(i), i=1,mycube%npts)


    CLOSE (11)
    ! ...
    ! access attributes as follows:
    ! mycube%str1 = ...
  END SUBROUTINE cube_get


  FUNCTION cube_add (mycube1, mycube2)
    TYPE(cube) :: cube_add
    TYPE(cube), INTENT(IN) :: mycube1, mycube2
    INTEGER :: i

    cube_add%str1  = mycube1%str1
    cube_add%str2  = mycube1%str2
    cube_add%xmin  = mycube1%xmin
    cube_add%ymin  = mycube1%ymin
    cube_add%zmin  = mycube1%zmin
    cube_add%dx    = mycube1%dx
    cube_add%dy    = mycube1%dy
    cube_add%dz    = mycube1%dz
    cube_add%nx    = mycube1%nx
    cube_add%ny    = mycube1%ny
    cube_add%nz    = mycube1%nz
    cube_add%natom = mycube1%natom
    cube_add%npts  = mycube1%npts

    ALLOCATE(cube_add%zahl(cube_add%natom))
    ALLOCATE(cube_add%chrg(cube_add%natom))
    ALLOCATE(cube_add%x(cube_add%natom))
    ALLOCATE(cube_add%y(cube_add%natom))
    ALLOCATE(cube_add%z(cube_add%natom))
    ALLOCATE(cube_add%array(cube_add%npts))

    cube_add%zahl = mycube1%zahl
    cube_add%chrg = mycube1%chrg
    cube_add%x    = mycube1%x
    cube_add%y    = mycube1%y
    cube_add%z    = mycube1%z

    DO i = 1, cube_add%npts
      cube_add%array(i) = mycube1%array(i) + mycube2%array(i)
    END DO
  END FUNCTION cube_add


  FUNCTION cube_sub (mycube1, mycube2)
    TYPE(cube) :: cube_sub
    TYPE(cube), INTENT(IN) :: mycube1, mycube2
    INTEGER :: i

    cube_sub%str1  = mycube1%str1
    cube_sub%str2  = mycube1%str2
    cube_sub%xmin  = mycube1%xmin
    cube_sub%ymin  = mycube1%ymin
    cube_sub%zmin  = mycube1%zmin
    cube_sub%dx    = mycube1%dx
    cube_sub%dy    = mycube1%dy
    cube_sub%dz    = mycube1%dz
    cube_sub%nx    = mycube1%nx
    cube_sub%ny    = mycube1%ny
    cube_sub%nz    = mycube1%nz
    cube_sub%natom = mycube1%natom
    cube_sub%npts  = mycube1%npts

    ALLOCATE(cube_sub%zahl(cube_sub%natom))
    ALLOCATE(cube_sub%chrg(cube_sub%natom))
    ALLOCATE(cube_sub%x(cube_sub%natom))
    ALLOCATE(cube_sub%y(cube_sub%natom))
    ALLOCATE(cube_sub%z(cube_sub%natom))
    ALLOCATE(cube_sub%array(cube_sub%npts))

    cube_sub%zahl = mycube1%zahl
    cube_sub%chrg = mycube1%chrg
    cube_sub%x    = mycube1%x
    cube_sub%y    = mycube1%y
    cube_sub%z    = mycube1%z

    DO i = 1, cube_sub%npts
      cube_sub%array(i) = mycube1%array(i) - mycube2%array(i)
    END DO
  END FUNCTION cube_sub


  FUNCTION cube_int (mycube)
    REAL (KIND=wp) :: cube_int
    TYPE (cube), INTENT(IN) :: mycube
    cube_int = SUM(mycube%array) * mycube%dx * mycube%dy * mycube%dz
  END FUNCTION cube_int


  FUNCTION cube_cdz (mycube)
    REAL (KIND=wp), ALLOCATABLE :: cube_cdz(:)
    TYPE (cube), INTENT(IN) :: mycube
    INTEGER :: ix, iy, iz, idx
    REAL (KIND=wp) :: slice

    ALLOCATE(cube_cdz(mycube%nz))
    cube_cdz = 0.0_wp

    DO iz = 1, mycube%nz
      slice = 0.0_wp
      DO iy = 1, mycube%ny
        DO ix = 1, mycube%nx
          idx = ((iz - 1) * mycube%ny + (iy - 1)) * mycube%nx + ix
          slice = slice + mycube%array(idx)
        END DO
      END DO
      IF (iz == 1) THEN
        cube_cdz(iz) = slice * mycube%dx * mycube%dy * mycube%dz
      ELSE
        cube_cdz(iz) = cube_cdz(iz - 1) + slice * mycube%dx * mycube%dy * mycube%dz
      END IF
    END DO
  END FUNCTION cube_cdz


  SUBROUTINE cube_write (filename, header, mycube)
    CHARACTER(LEN=*), INTENT(IN) :: filename, header
    TYPE (cube), INTENT(IN) :: mycube
    INTEGER :: unit, i

    OPEN(NEWUNIT=unit, FILE=filename, STATUS='replace', ACTION='write')
    WRITE(unit,'(A)') header
    DO i = 1, mycube%npts
      WRITE(unit,'(I10,1X,ES20.10)') i, mycube%array(i)
    END DO
    CLOSE(unit)
  END SUBROUTINE cube_write


  SUBROUTINE cube_del (mycube)
    TYPE (cube), INTENT(INOUT) :: mycube
    IF (ASSOCIATED(mycube%zahl))  DEALLOCATE(mycube%zahl)
    IF (ASSOCIATED(mycube%chrg))  DEALLOCATE(mycube%chrg)
    IF (ASSOCIATED(mycube%x))     DEALLOCATE(mycube%x)
    IF (ASSOCIATED(mycube%y))     DEALLOCATE(mycube%y)
    IF (ASSOCIATED(mycube%z))     DEALLOCATE(mycube%z)
    IF (ASSOCIATED(mycube%array)) DEALLOCATE(mycube%array)
  END SUBROUTINE cube_del


END MODULE cubes
