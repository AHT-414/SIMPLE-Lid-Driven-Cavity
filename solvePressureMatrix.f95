subroutine solvePressureMatrix
  use variables
  implicit none

  double precision, dimension(nPoints**2)              :: rhs
  double precision, dimension(nPoints**2,nPoints**2)   :: AP
  integer                                              :: position,z
  integer :: pivot(nPoints**2),ok

  ! Set up RHS matrix
  do i = 1,nPoints
    do j = 1,nPoints
      position = i + (j-1)*nPoints
      rhs(position) = -rho*( (uStar(i+1,j)-uStar(i,j))*dy + &
        (vStar(i,j+1)-vStar(i,j))*dx )
    enddo
  enddo
  rhs(1) = 0

  ! Set up pressure correction coefficient matrix
  do i = 1,nPoints
    do j = 1,nPoints
      position = i + (j-1)*nPoints

      ! Conditionals for 4 corners
      if ((i .eq. 1) .and. (j .eq. 1)) then ! Bottom left corner
        AP(position,position) = 1 ! Special case
        cycle
      endif
      if ((i .eq. 1) .and. (j .eq. nPoints)) then ! Top left corner
        ! only set APE and APS
        ! East
          AP(position,position+1) = -rho*dy*du(i+1,j)
        ! South
          AP(position,position-nPoints) = -rho*dx*dv(i,j)
        ! Center
          AP(position,position) = -AP(position,position+1) -&
                                   AP(position,position-nPoints)
        cycle
      endif
      if ((i .eq. nPoints) .and. (j .eq. 1)) then ! bottom right corner
        ! Only need West and North
        ! North
          AP(position,position+nPoints) = -rho*dx*dv(i,j+1)
        ! West
          AP(position,position-1) = -rho*dy*du(i,j)

        ! Center
          AP(position,position) = -AP(position,position+nPoints) - &
                                   AP(position,position-1)
        cycle
      endif
      if ((i .eq. nPoints) .and. (j .eq. nPoints)) then ! top right corner
        ! Set South and West
        ! South
          AP(position,position-nPoints) = -rho*dx*dv(i,j)
        ! West
          AP(position,position-1) = -rho*dy*du(i,j)

        ! Center
          AP(position,position) = -AP(position,position-nPoints) - &
                                   AP(position,position-1)

        cycle
      endif


      ! Condionals for Walls
      if (i .eq. 1) then ! Left Wall
        ! Set north, south, and east
        ! North
          AP(position,position+nPoints) = -rho*dx*dv(i,j+1)
        ! South
          AP(position,position-nPoints) = -rho*dx*dv(i,j)
        ! East
          AP(position,position+1) = -rho*dy*du(i+1,j)

        ! Center
          AP(position,position) = -AP(position,position+nPoints) - &
                                   AP(position,position-nPoints) - &
                                   AP(position,position+1)
        cycle
      endif
      if (i .eq. nPoints) then ! Right Wall
        ! Set north, south, and west
        ! North
          AP(position,position+nPoints) = -rho*dx*dv(i,j+1)
        ! South
          AP(position,position-nPoints) = -rho*dx*dv(i,j)
        ! West
          AP(position,position-1) = -rho*dy*du(i,j)

        ! Center
          AP(position,position) = -AP(position,position+nPoints) - &
                                   AP(position,position-nPoints) - &
                                   AP(position,position-1)
        cycle
      endif
      if (j .eq. 1) then ! Bottom wall
        ! Set North, East, adn West
        ! North
          AP(position,position+nPoints) = -rho*dx*dv(i,j+1)
        ! East
          AP(position,position+1) = -rho*dy*du(i+1,j)
        ! West
          AP(position,position-1) = -rho*dy*du(i,j)

        ! Center
          AP(position,position) = -AP(position,position+nPoints) - &
                                   AP(position,position+1) - &
                                   AP(position,position-1)
        cycle
      endif
      if (j .eq. nPoints) then ! Top Wall
        ! Set East, West, South
        ! East
          AP(position,position+1) = -rho*dy*du(i+1,j)
        ! West
          AP(position,position-1) = -rho*dy*du(i,j)
        ! South
          AP(position,position-nPoints) = -rho*dx*dv(i,j)

        ! Center
          AP(position,position) = -AP(position,position-nPoints) - &
                                   AP(position,position+1) - &
                                   AP(position,position-1)
        cycle
      endif


      ! Interior points
      ! North
        AP(position,position+nPoints) = -rho*dx*dv(i,j+1)
      ! South
        AP(position,position-nPoints) = -rho*dx*dv(i,j)
      ! East
        AP(position,position+1) = -rho*dy*du(i+1,j)
      ! West
        AP(position,position-1) = -rho*dy*du(i,j)

      ! Center
      AP(position,position) = -AP(position,position+nPoints) - &
                               AP(position,position-nPoints) - &
                               AP(position,position+1) - &
                               AP(position,position-1)

    enddo
  enddo



  call DGESV(nPoints**2, 1, AP, nPoints**2,pivot,rhs,nPoints**2,ok)

  ! Update pPrime
  z = 1
  do j = 1,nPoints
    do i = 1,nPoints
      pPrime(i,j) = rhs(z)
      z = z+1
    end do
  enddo

  !do i = 1,nPoints**2
  !  print *, rhs(i)
  !enddo



end subroutine
