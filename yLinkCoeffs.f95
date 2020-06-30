subroutine yLinkCoeffs
  use variables
  implicit none
  real                          :: De, Dw, Dn, Ds
  real                          :: Fe, Fw, Fn, Fs
  real                          :: Ae, Aw, An, As

  De  = mu*dy / dx
  Dw  = mu*dy / dx
  Dn  = mu*dx / dy
  Ds  = mu*dx / dy

  ! Loop through interior Points
  do i = 2,nPoints-1
   do j = 2,nPoints

     Fe  = .5*rho*dy*(uVelocity(i+1,j)+uVelocity(i+1,j-1))
     Fw  = .5*rho*dy*(uVelocity(i,j)+uVelocity(i,j-1))
     Fn  = .5*rho*dx*(vVelocity(i,j)+vVelocity(i,j+1))
     Fs  = .5*rho*dx*(vVelocity(i,j-1)+vVelocity(i,j))

     Ae = max(zero, ((1-0.1 * abs(Fe/De))**5 ) )
     Aw = max(zero, ((1-0.1 * abs(Fw/Dw))**5 ) )
     An = max(zero, ((1-0.1 * abs(Fn/Dn))**5 ) )
     As = max(zero, ((1-0.1 * abs(Fs/Ds))**5 ) )

     aye = De * Ae + max(-Fe,zero)
     ayw = Dw * Aw + max(Fw,zero)
     ayn = Dn * An + max(-Fn,zero)
     ays = Ds * As + max(Fs,zero)
     ayo = aye + ayw + ayn + ays + &
       (Fe - Fw) + (Fn - Fs)

     syo = (pressure(i,j-1)-pressure(i,j)) * dx

     vStar(i,j) = alphaUV/ayo * &
       ( (aye*vVelocity(i+1,j)+ayw*vVelocity(i-1,j)+ayn*vVelocity(i,j+1)&
         +ays*vVelocity(i,j-1)) + syo ) + (1-alphaUV)*vVelocity(i,j)

     dv(i,j) = alphaUV * dy / axo

   enddo
 enddo


  ! Loop through Left and Right Walls
  ! Left Wall
  i = 1
  do j = 2,nPoints
    Fe  = .5*rho*dy*(uVelocity(i+1,j)+uVelocity(i+1,j-1))
    Fw  = 0
    Fn  = .5*rho*dx*(vVelocity(i,j)+vVelocity(i,j+1))
    Fs  = .5*rho*dx*(vVelocity(i,j-1)+vVelocity(i,j))

    Ae = max(zero, ((1-0.1 * abs(Fe/De))**5 ) )
    Aw = 0
    An = max(zero, ((1-0.1 * abs(Fn/Dn))**5 ) )
    As = max(zero, ((1-0.1 * abs(Fs/Ds))**5 ) )

    aye = De * Ae + max(-Fe,zero)
    ayw = 0
    ayn = Dn * An + max(-Fn,zero)
    ays = Ds * As + max(Fs,zero)
    ayo = aye + ayw + ayn + ays + &
     (Fe - Fw) + (Fn - Fs)

    dv(i,j) = alphaUV * dy / axo

  enddo

  ! Right Wall
  i = nPoints
  do j = 2,nPoints
    Fe  = 0
    Fw  = .5*rho*dy*(uVelocity(i,j)+uVelocity(i,j-1))
    Fn  = .5*rho*dx*(vVelocity(i,j)+vVelocity(i,j+1))
    Fs  = .5*rho*dx*(vVelocity(i,j-1)+vVelocity(i,j))

    Ae = 0
    Aw = max(zero, ((1-0.1 * abs(Fw/Dw))**5 ) )
    An = max(zero, ((1-0.1 * abs(Fn/Dn))**5 ) )
    As = max(zero, ((1-0.1 * abs(Fs/Ds))**5 ) )

    aye = 0
    ayw = De * Ae + max(-Fe,zero)
    ayn = Dn * An + max(-Fn,zero)
    ays = Ds * As + max(Fs,zero)
    ayo = aye + ayw + ayn + ays + &
     (Fe - Fw) + (Fn - Fs)

    dv(i,j) = alphaUV * dy / axo
  enddo

  ! Apply Boundary Conditions
    ! top
    vStar(1:nPoints,nPoints+1) = -vStar(1:nPoints,nPoints)
    ! Bottom
    vStar(1:nPoints,1) = -vStar(1:nPoints,2)
    ! Left Wall
    vStar(1,1:nPoints+1) = 0
    !Right Wall
    vStar(nPoints,nPoints+1) = 0


end subroutine
