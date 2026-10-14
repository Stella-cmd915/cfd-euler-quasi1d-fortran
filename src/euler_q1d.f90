!=================================================================
! Pseudo-1D (Quasi-1D) Euler solver for a converging-diverging duct
! Finite-Volume, RK4, Roe or FVS (Steger-Warming),
! 1st or 2nd order MUSCL with Van Leer limiter.
!
! Solves:    dU/dt + dF/dx = q
! with      U  = [rho, rho*u, rho*E]^T
!           F  = [rho*u, rho*u^2 + p, u*(rho*E + p)]^T
!           q  = -(1/S) dS/dx [rho*u, rho*u^2, u*(rho*E + p)]^T
!=================================================================
module constants
  implicit none
  double precision, parameter :: GAM   = 1.4d0
  double precision, parameter :: GAMM1 = GAM - 1.d0
  double precision, parameter :: Rgas  = 287.03d0
  double precision, parameter :: pt_in = 1.0332d5
  double precision, parameter :: Tt_in = 288.d0
  double precision :: rkconst(4)
end module constants

module variables
  implicit none
  integer          :: NX, NXC
  double precision :: XMIN, XMAX, DX, DT
end module variables

module mathfuncs
  use constants
  implicit none
contains

  ! Van Leer flux limiter Psi(r) = (r + |r|)/(1 + |r|)
  double precision function vanleer(r)
    implicit none
    double precision, intent(in) :: r
    if (r <= 0.d0) then
      vanleer = 0.d0
    else
      vanleer = (r + abs(r)) / (1.d0 + abs(r))
    endif
  end function vanleer

  ! A/A* relation as function of Mach number
  double precision function area_mach(M)
    implicit none
    double precision, intent(in) :: M
    area_mach = (1.d0/M) * &
      ( (2.d0/(GAM+1.d0)) * (1.d0 + 0.5d0*GAMM1*M*M) ) &
      **((GAM+1.d0)/(2.d0*GAMM1))
  end function area_mach

  ! Mach from A/A* (subsonic branch, Newton-Raphson)
  double precision function mach_from_area(AR)
    implicit none
    double precision, intent(in) :: AR
    double precision :: M, f, df
    integer :: it
    M = 0.2d0
    do it = 1, 80
      f  = area_mach(M) - AR
      df = (area_mach(M+1.d-6) - area_mach(M-1.d-6)) / (2.d-6)
      if (abs(df) < 1.d-14) exit
      M  = M - f/df
      if (M <= 0.d0) M = 1.d-3
      if (M >= 1.d0) M = 0.999d0
    enddo
    mach_from_area = M
  end function mach_from_area

end module mathfuncs

program Euler_Q1D
  use constants
  use variables
  use mathfuncs
  implicit none

  double precision :: CFL, TIME, xN, Smin, aa, Res, p_back
  double precision :: UCL(3), UCR(3), UPL(3), UPR(3), NORM, flux(3)
  double precision :: UPcell(3), UP(3), UC(3)
  double precision, allocatable :: q(:,:), geom(:), dSdx(:), cell_flux(:,:)
  double precision, allocatable :: U(:,:), U0(:,:), Uold(:,:)
  double precision, allocatable :: x(:), S(:)
  integer :: i, k, ntime, NTIMETOT, IORDER, K1, K2
  integer :: IRIEMANN
  character(len=80) :: tag, fres, fmach, fconv, fgeom, finp

  ! ----- Standard 4-stage low-storage RK coefficients (Jameson) -----
  rkconst(1) = 0.1084d0
  rkconst(2) = 0.2602d0
  rkconst(3) = 0.5052d0
  rkconst(4) = 1.d0

  ! ----- Read input -----
  call get_command_argument(1, finp)
  if (len_trim(finp) == 0) finp = 'eulerf.txt'

  open(1,file=trim(finp),status='old')
  read(1,*) NTIMETOT
  read(1,*) XMIN, XMAX
  read(1,*) DX, CFL
  read(1,*) IRIEMANN          ! 1 = Roe, 2 = FVS
  read(1,*) IORDER            ! 1 or 2
  read(1,*) K1, K2
  read(1,*) p_back
  read(1,*) tag               ! suffix for output files
  close(1)

  fres  = 'results_'//trim(tag)//'.dat'
  fmach = 'mach_'   //trim(tag)//'.dat'
  fconv = 'conv_'   //trim(tag)//'.dat'
  fgeom = 'geom_'   //trim(tag)//'.dat'

  ! ----- Grid (faces -> NX nodes; cells -> NXC = NX-1) -----
  NX  = int((XMAX-XMIN)/DX) + 1
  DX  = (XMAX-XMIN)/dble(NX-1)
  NXC = NX - 1

  allocate(x(1:NX), S(1:NX), dSdx(1:NX))

  do i = 1, NX
    x(i) = XMIN + dble(i-1)*DX
  enddo

  ! ----- Geometry: throat at xN, S(xN)=Smin -----
  xN   = 0.25d0 + dble(K2)/50.d0
  Smin = 0.80d0*1.d0 + dble(K1)/100.d0   ! S(0)=S(1)=1

  ! Two parabolas joined at xN, both equal to Smin there,
  ! both equal to 1 at the duct ends (x=0 and x=1).
  do i = 1, NX
    if (x(i) < 0.d0) then
      S(i)    = 1.d0
      dSdx(i) = 0.d0
    else if (x(i) > 1.d0) then
      S(i)    = 1.d0
      dSdx(i) = 0.d0
    else if (x(i) <= xN) then
      aa = (1.d0 - Smin)/(xN*xN)
      S(i)    = Smin + aa*(x(i)-xN)**2
      dSdx(i) = 2.d0*aa*(x(i)-xN)
    else
      aa = (1.d0 - Smin)/((1.d0-xN)**2)
      S(i)    = Smin + aa*(x(i)-xN)**2
      dSdx(i) = 2.d0*aa*(x(i)-xN)
    endif
  enddo

  ! Save geometry (face data + cell-center)
  open(40,file=trim(fgeom))
  write(40,'(a)') '# x  S(x)  dS/dx'
  do i = 1, NX
    write(40,'(3e22.12)') x(i), S(i), dSdx(i)
  enddo
  close(40)

  ! ----- Cell-centered geometric source factor (1/S) dS/dx -----
  allocate(geom(1:NXC))
  do i = 1, NXC
    geom(i) = 0.5d0*( dSdx(i)/S(i) + dSdx(i+1)/S(i+1) )
  enddo

  ! ----- Allocate solution and work arrays -----
  allocate(U(3,0:NXC+1), U0(3,0:NXC+1), Uold(3,1:NXC))
  allocate(cell_flux(1:3,0:NXC+1))
  allocate(q(3,0:NXC+1))

  TIME = 0.d0
  call ic(U, NXC)

  open(30,file=trim(fconv))
  write(30,'(a)') '# iter  time  residual_rho'

  ! ===========================================================
  ! Main time-marching loop
  ! ===========================================================
  do ntime = 1, NTIMETOT

    Uold(1:3,1:NXC) = U(1:3,1:NXC)

    ! Local CFL-based DT
    DT = 1.d30
    do i = 1, NXC
      call co_to_pri(U(1:3,i), UPcell)
      if (UPcell(1) <= 0.d0 .or. UPcell(3) <= 0.d0) cycle
      aa = sqrt(GAM*UPcell(3)/UPcell(1))
      DT = min(DT, CFL*DX / (abs(UPcell(2)) + aa))
    enddo
    if (DT <= 0.d0 .or. DT > 1.d2) then
      write(*,*) 'BAD DT = ', DT, ' at NTIME=', ntime
      stop
    endif

    k = 1
    do while (k <= 4)

      call bc(U, p_back)

      cell_flux = 0.d0
      do i = 1, NX
        if (IORDER == 2 .and. i >= 3 .and. i <= NX-1) then
          call muscl_recon(U(1:3,i-2:i+1), UCL, UCR)
        else
          UCL = U(1:3,i-1)
          UCR = U(1:3,i)
        endif
        call co_to_pri(UCL, UPL)
        call co_to_pri(UCR, UPR)
        NORM = 1.d0
        if (IRIEMANN == 1) then
          call calc_flux_roe(UPL, UPR, UCL, UCR, NORM, flux)
        else
          call calc_flux_fvs(UPL, UPR, UCL, UCR, flux)
        endif
        cell_flux(1:3, i-1) = cell_flux(1:3, i-1) + flux
        cell_flux(1:3, i  ) = cell_flux(1:3, i  ) - flux
      enddo

      ! Geometric source term (full, per Q1D Euler equations)
      do i = 1, NXC
        call co_to_pri(U(1:3,i), UPcell)
        q(1,i) = -U(2,i)*geom(i)
        q(2,i) = -UPcell(1)*UPcell(2)**2 * geom(i)
        q(3,i) = -(U(3,i) + UPcell(3))*UPcell(2) * geom(i)
      enddo

      if (k == 1) U0(1:3,1:NXC) = U(1:3,1:NXC)
      U(1:3,1:NXC) = U0(1:3,1:NXC) &
                   - DT/DX * rkconst(k) * cell_flux(1:3,1:NXC) &
                   + DT     * rkconst(k) * q(1:3,1:NXC)
      k = k + 1
    enddo

    TIME = TIME + DT

    call compute_residual(Uold, U, NXC, Res)
    write(30,'(i8,2(1x,e20.10))') ntime, TIME, Res
    if (mod(ntime,200) == 0) write(*,'(a,i7,a,e12.4,a,e12.4)') &
      ' it=', ntime, '   t=', TIME, '   res=', Res
    if (Res /= Res) then
      write(*,*) 'NaN at ntime=', ntime
      stop
    endif

  enddo
  close(30)

  ! Final fields
  call writeres(trim(fres), TIME)
  call writeres_mach(trim(fmach), S)

contains
  !---------------------------------------------------------
  subroutine ic(U, NXC)
    use constants
    implicit none
    integer, intent(in) :: NXC
    double precision, intent(out) :: U(3,0:NXC+1)
    double precision :: UCloc(3), UPloc(3)
    double precision :: M0, T, p, rho, vel
    integer :: i
    M0  = 0.20d0
    T   = Tt_in / (1.d0 + 0.5d0*GAMM1*M0*M0)
    p   = pt_in / (1.d0 + 0.5d0*GAMM1*M0*M0)**(GAM/GAMM1)
    rho = p / (Rgas*T)
    vel = M0 * sqrt(GAM*Rgas*T)
    UPloc(1) = rho; UPloc(2) = vel; UPloc(3) = p
    call pri_to_co(UCloc, UPloc)
    do i = 0, NXC+1
      U(1:3,i) = UCloc(1:3)
    enddo
  end subroutine ic

  !---------------------------------------------------------
  subroutine bc(U, p_back)
    use constants
    use variables
    implicit none
    double precision, intent(inout) :: U(3,0:NXC+1)
    double precision, intent(in)    :: p_back
    double precision :: UPloc(3), UCloc(3)
    double precision :: rho_i, u_i, p_i, M_in, T, p, rho, vel, E

    ! INLET — subsonic: keep total pt_in, Tt_in;
    ! infer Mach from interior cell 1 isentropically
    rho_i = U(1,1)
    u_i   = U(2,1)/rho_i
    p_i   = GAMM1*( U(3,1) - 0.5d0*rho_i*u_i*u_i )
    if (p_i <= 0.d0) p_i = 1.d0
    M_in  = u_i / sqrt(GAM*p_i/rho_i)
    if (M_in < 0.d0) M_in = 0.d0
    if (M_in > 0.99d0) M_in = 0.99d0

    T   = Tt_in / (1.d0 + 0.5d0*GAMM1*M_in*M_in)
    p   = pt_in / (1.d0 + 0.5d0*GAMM1*M_in*M_in)**(GAM/GAMM1)
    rho = p/(Rgas*T)
    vel = M_in*sqrt(GAM*Rgas*T)
    UPloc(1) = rho; UPloc(2) = vel; UPloc(3) = p
    call pri_to_co(UCloc, UPloc)
    U(1:3,0) = UCloc(1:3)

    ! OUTLET — subsonic: impose static pressure p_back,
    ! extrapolate rho, u from interior cell NXC
    rho = U(1,NXC)
    vel = U(2,NXC)/rho
    p   = p_back
    UPloc(1) = rho; UPloc(2) = vel; UPloc(3) = p
    call pri_to_co(UCloc, UPloc)
    U(1:3,NXC+1) = UCloc(1:3)
  end subroutine bc

  !---------------------------------------------------------
  subroutine writeres(fname, TIME)
    use variables
    implicit none
    character(len=*), intent(in) :: fname
    double precision, intent(in) :: TIME
    double precision :: UPloc(3), UCloc(3), aa, M
    integer :: i
    open(11,file=fname)
    write(11,'(a)') '# x   rho   u   p   Mach'
    do i = 1, NXC
      UCloc(1:3) = U(1:3,i)
      call co_to_pri(UCloc, UPloc)
      aa = sqrt(GAM*UPloc(3)/UPloc(1))
      M  = UPloc(2)/aa
      write(11,'(5e22.12)') XMIN+(dble(i)-0.5d0)*DX, UPloc(1), UPloc(2), UPloc(3), M
    enddo
    close(11)
  end subroutine writeres

  !---------------------------------------------------------
  subroutine writeres_mach(fname, Sloc)
    use constants
    use variables
    use mathfuncs
    implicit none
    character(len=*), intent(in) :: fname
    double precision, intent(in) :: Sloc(1:NX)
    double precision :: MachN, MachA, Astar, Aratio, UPloc(3), aa, Sc
    integer :: i

    ! Use the cell-center area (interpolate between face values)
    Astar = Smin   ! virtual A* = throat area (sonic-throat reference)
    open(12,file=fname)
    write(12,'(a)') '# x   S_cell   Mach_numeric   Mach_analytic   A/A*'
    do i = 1, NXC
      call co_to_pri(U(1:3,i), UPloc)
      aa    = sqrt(GAM*UPloc(3)/UPloc(1))
      MachN = UPloc(2)/aa
      Sc    = 0.5d0*(Sloc(i)+Sloc(i+1))
      Aratio = Sc/Astar
      if (Aratio < 1.d0) Aratio = 1.d0
      MachA = mach_from_area(Aratio)
      write(12,'(5e22.12)') XMIN+(dble(i)-0.5d0)*DX, Sc, MachN, MachA, Aratio
    enddo
    close(12)
  end subroutine writeres_mach

  !---------------------------------------------------------
  subroutine co_to_pri(UC, UP)
    use constants
    implicit none
    double precision, intent(in)  :: UC(3)
    double precision, intent(out) :: UP(3)
    UP(1) = UC(1)
    UP(2) = UC(2)/UC(1)
    UP(3) = GAMM1*(UC(3) - 0.5d0/UC(1) * UC(2)**2)
  end subroutine co_to_pri

  subroutine pri_to_co(UC, UP)
    use constants
    implicit none
    double precision, intent(out) :: UC(3)
    double precision, intent(in)  :: UP(3)
    UC(1) = UP(1)
    UC(2) = UP(1)*UP(2)
    UC(3) = UP(3)/GAMM1 + 0.5d0*UP(1)*UP(2)**2
  end subroutine pri_to_co

  !---------------------------------------------------------
  subroutine calc_flux_roe(UPL, UPR, UCL, UCR, NORM, F)
    implicit none
    double precision, intent(in)  :: UPL(3), UPR(3), UCL(3), UCR(3), NORM
    double precision, intent(out) :: F(3)
    double precision :: FL(3), FR(3), Aroe(3,3), DU(3)
    call flux1d(UPL, UCL, NORM, FL)
    call flux1d(UPR, UCR, NORM, FR)
    call roematrix(UPL, UCL, UPR, UCR, NORM, Aroe)
    DU = UCR - UCL
    F  = 0.5d0*(FL + FR - matmul(Aroe, DU))
  end subroutine calc_flux_roe

  subroutine flux1d(UP, UC, NORM, F)
    implicit none
    double precision, intent(in)  :: UP(3), UC(3), NORM
    double precision, intent(out) :: F(3)
    double precision :: VN
    VN   = UP(2)*NORM
    F(1) = UP(1)*VN
    F(2) = UP(1)*UP(2)*VN + UP(3)*NORM
    F(3) = (UC(3) + UP(3))*VN
  end subroutine flux1d

  subroutine roematrix(UPL, UCL, UPR, UCR, NORM, Aroe)
    use constants
    implicit none
    double precision, intent(in)  :: UPL(3), UCL(3), UPR(3), UCR(3), NORM
    double precision, intent(out) :: Aroe(3,3)
    double precision :: rL, rR, Uroe, Hl, Hr, Hroe, croe, L1, L2, L3, delta
    double precision :: left(3,3), right(3,3), eigen(3,3)
    rL = sqrt(UPL(1)); rR = sqrt(UPR(1))
    Uroe = (rL*UPL(2) + rR*UPR(2))/(rL+rR)
    Hl   = (UCL(3)+UPL(3))/UPL(1)
    Hr   = (UCR(3)+UPR(3))/UPR(1)
    Hroe = (rL*Hl + rR*Hr)/(rL+rR)
    croe = sqrt( GAMM1*(Hroe - 0.5d0*Uroe**2) )
    call leftright(rL*rR, croe, Uroe, left, right)
    L1 = abs(Uroe);  L2 = abs(Uroe+croe);  L3 = abs(Uroe-croe)
    delta = 0.05d0*croe
    if (L1 <= delta) L1 = (L1*L1+delta*delta)/(2.d0*delta)
    if (L2 <= delta) L2 = (L2*L2+delta*delta)/(2.d0*delta)
    if (L3 <= delta) L3 = (L3*L3+delta*delta)/(2.d0*delta)
    eigen = 0.d0
    eigen(1,1)=L1; eigen(2,2)=L2; eigen(3,3)=L3
    Aroe = matmul(right, eigen)
    Aroe = matmul(Aroe,  left)
  end subroutine roematrix

  subroutine leftright(D, C, U, left, right)
    use constants
    implicit none
    double precision, intent(in)  :: D, C, U
    double precision, intent(out) :: left(3,3), right(3,3)
    left(1,1) =  1.d0 - 0.5d0*GAMM1*U*U/(C*C)
    left(1,2) =  GAMM1*U/(C*C)
    left(1,3) = -GAMM1/(C*C)
    left(2,1) = (0.5d0*GAMM1*U*U - U*C)/(D*C)
    left(2,2) = (C - GAMM1*U)/(D*C)
    left(2,3) =  GAMM1/(D*C)
    left(3,1) = -(0.5d0*GAMM1*U*U + U*C)/(D*C)
    left(3,2) =  (C + GAMM1*U)/(D*C)
    left(3,3) = -GAMM1/(D*C)
    right(1,1) = 1.d0
    right(1,2) =  0.5d0*D/C
    right(1,3) = -0.5d0*D/C
    right(2,1) = U
    right(2,2) =  0.5d0*(U+C)*D/C
    right(2,3) = -0.5d0*(U-C)*D/C
    right(3,1) = 0.5d0*U*U
    right(3,2) =  (0.5d0*U*U + U*C + C*C/GAMM1)*0.5d0*D/C
    right(3,3) = -(0.5d0*U*U - U*C + C*C/GAMM1)*0.5d0*D/C
  end subroutine leftright

  !---------------------------------------------------------
  subroutine calc_flux_fvs(UPL, UPR, UCL, UCR, F)
    use constants
    implicit none
    double precision, intent(in)  :: UPL(3), UPR(3), UCL(3), UCR(3)
    double precision, intent(out) :: F(3)
    double precision :: FLp(3), FRm(3)
    call flux_fvs_pm(UPL, UCL, +1, FLp)
    call flux_fvs_pm(UPR, UCR, -1, FRm)
    F = FLp + FRm
  end subroutine calc_flux_fvs

  ! Steger-Warming FVS (subsonic + supersonic safe)
  subroutine flux_fvs_pm(UP, UC, sgn, F)
    use constants
    implicit none
    double precision, intent(in)  :: UP(3), UC(3)
    integer, intent(in)           :: sgn
    double precision, intent(out) :: F(3)
    double precision :: rho, u, p, H, a
    double precision :: l1, l2, l3, c2g
    rho = UP(1); u = UP(2); p = UP(3)
    a   = sqrt(GAM*p/rho)
    H   = (UC(3) + p)/rho
    l1 = u
    l2 = u + a
    l3 = u - a
    if (sgn == 1) then
      l1 = 0.5d0*(l1 + abs(l1))
      l2 = 0.5d0*(l2 + abs(l2))
      l3 = 0.5d0*(l3 + abs(l3))
    else
      l1 = 0.5d0*(l1 - abs(l1))
      l2 = 0.5d0*(l2 - abs(l2))
      l3 = 0.5d0*(l3 - abs(l3))
    endif
    ! Steger-Warming closed form (Anderson, CFD: principles & applications)
    c2g  = a*a/(GAM-1.d0)
    F(1) = (rho/(2.d0*GAM)) * ( 2.d0*(GAM-1.d0)*l1 + l2 + l3 )
    F(2) = (rho/(2.d0*GAM)) * ( 2.d0*(GAM-1.d0)*l1*u + l2*(u+a) + l3*(u-a) )
    F(3) = (rho/(2.d0*GAM)) * ( (GAM-1.d0)*l1*u*u &
            + 0.5d0*l2*(u+a)*(u+a) + 0.5d0*l3*(u-a)*(u-a) &
            + ((3.d0-GAM)/(2.d0*(GAM-1.d0))) * (l2 + l3)*a*a )
  end subroutine flux_fvs_pm

  !---------------------------------------------------------
  ! MUSCL reconstruction with Van Leer limiter (k=0)
  subroutine muscl_recon(Ust, UL, UR)
    use mathfuncs
    implicit none
    double precision, intent(in)  :: Ust(3,-1:2)
    double precision, intent(out) :: UL(3), UR(3)
    double precision :: dUm, dUp, r, slope
    integer :: kk
    do kk = 1, 3
      dUm = Ust(kk,0) - Ust(kk,-1)
      dUp = Ust(kk,1) - Ust(kk,0)
      if (abs(dUp) > 1.d-12) then
        r = dUm/dUp
      else
        r = 0.d0
      endif
      slope = vanleer(r)*dUp
      UL(kk) = Ust(kk,0) + 0.5d0*slope
    enddo
    do kk = 1, 3
      dUm = Ust(kk,1) - Ust(kk,0)
      dUp = Ust(kk,2) - Ust(kk,1)
      if (abs(dUp) > 1.d-12) then
        r = dUm/dUp
      else
        r = 0.d0
      endif
      slope = vanleer(r)*dUp
      UR(kk) = Ust(kk,1) - 0.5d0*slope
    enddo
  end subroutine muscl_recon

  !---------------------------------------------------------
  subroutine compute_residual(Uo, Un, NXC, Res)
    implicit none
    integer, intent(in) :: NXC
    double precision, intent(in) :: Uo(3,1:NXC), Un(3,1:NXC)
    double precision, intent(out) :: Res
    integer :: i
    Res = 0.d0
    do i = 1, NXC
      Res = Res + (Un(1,i) - Uo(1,i))**2
    enddo
    Res = sqrt(Res/dble(NXC))
  end subroutine compute_residual

end program Euler_Q1D
