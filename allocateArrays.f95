Subroutine allocateArrays
  use variables
  implicit none

  ! allocate size of u, v, and p arrays
  allocate(uVelocity(nPoints+1,nPoints))
  allocate(uStar(nPoints+1,nPoints))
  allocate(uPrime(nPoints+1,nPoints))

  allocate(vVelocity(nPoints,nPoints+1))
  allocate(vStar(nPoints,nPoints+1))
  allocate(vPrime(nPoints,nPoints+1))

  allocate(pressure(nPoints,nPoints))
  allocate(pStar(nPoints,nPoints))
  allocate(pPrime(nPoints,nPoints))

  allocate(ax(nPoints-1), bx(nPoints-1), cx(nPoints-1))
  allocate(ddx(nPoints-1), solx(nPoints-1))

  allocate(az(nPoints), bz(nPoints), cz(nPoints))
  allocate(ddz(nPoints), solz(nPoints))


  allocate(ape(nPoints,nPoints),apw(nPoints,nPoints),apn(nPoints,nPoints))
  allocate(aps(nPoints,nPoints),apo(nPoints,nPoints),spo(nPoints,nPoints))


  allocate(du(nPoints+1,nPoints),dv(nPoints,nPoints+1))

  allocate(xVelocity(nPoints,nPoints))
  allocate(yVelocity(nPoints,nPoints))
  allocate(velMag(nPoints,nPoints))

  pPrime = 0

  massSweep = 2
  momentumSweep = 5

  alphaP = .2
  alphaUV = .8

  zero = 0.0

end Subroutine
