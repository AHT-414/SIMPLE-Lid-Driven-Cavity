Subroutine Update
  use variables
  implicit none



  do i = 2, nPoints
    do j = 2, nPoints-1
      uVelocity(i,j) = ustar(i,j) + &
        (((pPrime(i-1,j) - pPrime(i,j))*du(i,j)))
    end do
  end do

  do i = 2,nPoints-1
    do j = 2,nPoints
      vVelocity(i,j) = vstar(i,j) + &
        (((pPrime(i,j-1) - pPrime(i,j))*dv(i,j)))
    end do
  end do



  pressure = pStar + alphaP*pPrime
  pressure(1,1) = 0

  ! Apply boundary Conditions
    ! Left Wall
    uVelocity(1,1:nPoints) = -uVelocity(2,1:nPoints)
    ! Right Wall
    uVelocity(nPoints+1,1:nPoints) = -uVelocity(nPoints,1:nPoints)
    ! Bottom
    uVelocity(1:nPoints+1,1) = 0
    ! Top
    uVelocity(1:nPoints+1,nPoints) = uLid



  ! Apply Boundary Conditions
    ! top
    vVelocity(1:nPoints,nPoints+1) = -vVelocity(1:nPoints,nPoints)
    ! Bottom
    vVelocity(1:nPoints,1) = -vVelocity(1:nPoints,2)
    ! Left Wall
    vVelocity(1,1:nPoints+1) = 0
    !Right Wall
    vVelocity(nPoints,nPoints+1) = 0

  uStar = uVelocity
  vStar = vVelocity
  pStar = pressure



end Subroutine
