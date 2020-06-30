subroutine xlinkCoeffs
  use variables
  implicit none
  real                          :: De, Dw, Dn, Ds
  real                          :: Fe, Fw, Fn, Fs
  real                          :: Ae, Aw, An, As

  De  = mu*dy / dx
  Dw  = mu*dy / dx
  Dn  = mu*dx / dy
  Ds  = mu*dx / dy

  ! Loop through all Points
    ! Except top and bottom, done in different loop
  do i = 2,nPoints
    do j = 2,nPoints-1
      Fe  = .5*rho*dy*(uVelocity(i+1,j)+uVelocity(i,j))
      Fw  = .5*rho*dy*(uVelocity(i-1,j)+uVelocity(i,j))
      Fn  = .5*rho*dx*(vVelocity(i,j+1)+vVelocity(i-1,j+1))
      Fs  = .5*rho*dx*(vVelocity(i,j)+vVelocity(i-1,j))

      Ae = max(zero, ((1-0.1 * abs(Fe/De))**5 ) )
      Aw = max(zero, ((1-0.1 * abs(Fw/Dw))**5 ) )
      An = max(zero, ((1-0.1 * abs(Fn/Dn))**5 ) )
      As = max(zero, ((1-0.1 * abs(Fs/Ds))**5 ) )

      axe = De * Ae + max(-Fe,zero)
      axw = Dw * Aw + max(Fw,zero)
      axn = Dn * An + max(-Fn,zero)
      axs = Ds * As + max(Fs,zero)
      axo = axe + axw + axn + axs + &
        (Fe - Fw) + (Fn - Fs)

      sxo = (pressure(i-1,j)-pressure(i,j)) * dy

      uStar(i,j) = alphaUV/axo * &
         ((axe*uVelocity(i+1,j)+axw*uVelocity(i-1,j)+axn*uVelocity(i,j+1) + &
           axs*uVelocity(i,j-1)) + sxo) + (1-alphaUV)*uVelocity(i,j)

      du(i,j) = alphaUV * dy / axo
    end do
  end do

  ! Top of domain
  j = nPoints
  do i = 2,nPoints
    Fe  = .5*rho*dy*(uVelocity(i+1,j)+uVelocity(i,j))
    Fw  = .5*rho*dy*(uVelocity(i-1,j)+uVelocity(i,j))
    Fn  = 0
    Fs  = .5*rho*dx*(vVelocity(i,j)+vVelocity(i-1,j))

    Ae = max(zero, ((1-0.1 * abs(Fe/De))**5 ) )
    Aw = max(zero, ((1-0.1 * abs(Fw/Dw))**5 ) )
    An = 0
    As = max(zero, ((1-0.1 * abs(Fs/Ds))**5 ) )

    axe = De * Ae + max(-Fe,zero)
    axw = Dw * Aw + max(Fw,zero)
    axn = 0
    axs = Ds * As + max(Fs,zero)
    axo = axe + axw + axn + axs + &
      (Fe - Fw) + (Fn - Fs)

    sxo = (pressure(i-1,j)-pressure(i,j)) * dy

    du(i,j) = alphaUV * dy / axo
  enddo

  ! Bottom of domain
  j = 1
  do i = 2,nPoints
    Fe  = .5*rho*dy*(uVelocity(i+1,j)+uVelocity(i,j))
    Fw  = .5*rho*dy*(uVelocity(i-1,j)+uVelocity(i,j))
    Fn  = .5*rho*dx*(vVelocity(i,j+1)+vVelocity(i-1,j+1))
    Fs  = 0

    Ae = max(zero, ((1-0.1 * abs(Fe/De))**5 ) )
    Aw = max(zero, ((1-0.1 * abs(Fw/Dw))**5 ) )
    An = max(zero, ((1-0.1 * abs(Fn/Dn))**5 ) )
    As = 0

    axe = De * Ae + max(-Fe,zero)
    axw = Dw * Aw + max(Fw,zero)
    axn = Dn * An + max(-Fn,zero)
    axs = 0
    axo = axe + axw + axn + axs + &
      (Fe - Fw) + (Fn - Fs)

    du(i,j) = alphaUV * dy / axo
  enddo

  ! Apply boundary Conditions
    ! Left Wall
    uStar(1,1:nPoints) = -uStar(2,1:nPoints)
    ! Right Wall
    uStar(nPoints+1,1:nPoints) = -uStar(nPoints,1:nPoints)
    ! Bottom
    uStar(1:nPoints+1,1) = 0
    ! Top
    uStar(1:nPoints+1,nPoints) = uLid

end subroutine
