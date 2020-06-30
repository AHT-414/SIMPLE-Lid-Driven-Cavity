! ==============================================================================
Program cavityFlow
  use variables
  implicit none

  integer :: k
  real    :: start, finish, time

  call cpu_time(start)

  ! Define constant fluid properties
  rho = 1 ! kg/m**3
  mu = 1e-3 ! Ns / m**2

  ! The cavity will be a 1m x 1m square
  length = .5 ! meter
  uLid = 1

  Re = real(uLid)*real(Length)*real(rho)/real(mu)
  print *, "Reynolds Number = ", ceiling(Re)

  ! Define the resolution of the cavity
  nPoints = 75
  dx = length/nPoints ! meter
  dy = length/nPoints ! meter

  call allocateArrays


  ! Set initial Conditions
  uVelocity = 0
  vVelocity = 0

  ! Set Lid uVelocity BC
  uVelocity(1:nPoints+1,nPoints) = uLid
  uStar(1:nPoints+1,nPoints) = uLid

  ! Initial guess for uStar, vStar, pStar
  uStar = uVelocity
  vStar = vVelocity
  pStar = pressure


  kLoop = 1000
  ! Outer SIMPLE LOOP
  DO k = 1,kLoop
    !print *, "Calling xLinkCoeffs"
    Call xLinkCoeffs

    !print *, "Calling yLinkCoeffs"
    Call yLinkCoeffs

    !print *, "Calling calcPressureCoeffs"
    !Call pressureCoeffs

    !print *, "Calling solvePressureCorrection"
    Call solvePressureMatrix

    !print *, "Calling Update"
    Call Update

    print *, "Loop:", k
    print *, "  "



  end do

  call cpu_time(finish)

  time = finish - start
  print *, "Total Time:", time,"sec"


  ! average velocity at pressure nodes
  do i = 1,nPoints
    do j = 1,nPoints
      xVelocity(i,j) = 0.5*(uVelocity(i,j) + uVelocity(i+1,j))
      yVelocity(i,j) = 0.5*(vVelocity(i,j) + vVelocity(i,j+1))
      velMag(i,j) = sqrt(xVelocity(i,j)**2 + yVelocity(i,j)**2)
    end do
  enddo

  open(unit=12,file="xVelocity.data")
  open(unit=13,file="yVelocity.data")
  open(unit=14,file="velMag.data")
  do j = nPoints,1,-1
    write(12,*) xVelocity(:,j)
    write(13,*) yVelocity(:,j)
    write(14,*) velMag(:,j)
  end do
  close(12)
  close(13)
  close(14)


End Program
! ==============================================================================
