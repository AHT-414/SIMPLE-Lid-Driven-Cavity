! ==============================================================================
Module variables
  implicit none

  ! Fluid dependent variables
  double precision, allocatable :: uVelocity(:,:), uStar(:,:), uPrime(:,:)
  double precision, allocatable :: vVelocity(:,:), vStar(:,:), vPrime(:,:)
  double precision, allocatable :: pressure(:,:), pStar(:,:), pPrime(:,:)


  ! Link coefficients
  double precision              :: axe, axw, axs, axn, axo, sxo
  double precision              :: aye, ayw, ays, ayn, ayo, syo
  double precision,allocatable  :: du(:,:), dv(:,:)
  double precision,allocatable  :: ape(:,:), apw(:,:), aps(:,:), &
                                   apn(:,:), apo(:,:), spo(:,:)


  ! Constant fluid properties
  double precision              :: rho, mu

  ! Mesh geometry and boundary Conditions
  double precision              :: length, dx, dy
  double precision              :: uLid
  real                          :: zero
  integer                       :: nPoints

  ! indexing variables
  integer                       :: i, j, l, m, n, kLoop, il

  ! Residual
  double precision              :: residual
  integer                       :: momentumSweep, massSweep

  ! Relaxation
  double precision              :: alphaUV, alphaP


  double precision              :: Re



  double precision, allocatable :: ax(:), bx(:), cx(:), ddx(:), solx(:)
  double precision, allocatable :: az(:), bz(:), cz(:), ddz(:), solz(:)


  double precision, allocatable :: xVelocity(:,:), yVelocity(:,:), velMag(:,:)

end module
! ==============================================================================
