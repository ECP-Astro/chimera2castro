subroutine PROBINIT (init,name,namlen,problo,probhi)

  use probdata_module
  use chimera_parser_module
  use mesa_parser_module
  use bl_constants_module
  use bl_error_module
  use bl_types
  use eos_module
  use parallel, only: parallel_IOProcessor
  use prob_params_module, only: center
  use interpolate_module, only: locate
  use model_interp_module, only: interp1d_linear
  use quadrature_module, only: quad_init

  implicit none
  integer :: init, namlen
  integer :: name(namlen)
  real (dp_t) :: problo(2), probhi(2)

  integer :: untin,i,j,k,dir
  real (dp_t) :: probhi_r, rcyl, zcyl, r, dr, dvol, volr, dvolr, gr
  real (dp_t) :: mass_chim, mass_mesa, vol_chim, vol_mesa
  real (dp_t) :: tmp1, tmp2, tmp3, domega

  namelist /fortin/ chimera_fname
  namelist /fortin/ mesa_fname
  namelist /fortin/ eos_input
  namelist /fortin/ interp_method
  namelist /fortin/ max_radius
  namelist /fortin/ radius_inner
  namelist /fortin/ rho_inner
  namelist /fortin/ temp_inner
  namelist /fortin/ i_inner
  namelist /fortin/ do_particles
  namelist /fortin/ use_quad
  namelist /fortin/ nquad

  !
  !     Build "probin" filename -- the name of file containing fortin namelist.
  !     
  integer, parameter :: maxlen = 127
  character(maxlen) :: probin
  character(maxlen) :: model
  integer :: ipp, ierr, ipp1, n

  ! assume axisymmetric
  center(1) = zero
  center(2) = half*(problo(2)+probhi(2))

  if (namlen .gt. maxlen) call bl_error("probin file name too long")

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! set namelist defaults
  chimera_fname = ""
  mesa_fname = ""
  interp_method = 1
  eos_input = eos_input_rt
  max_radius = zero
  radius_inner = zero
  rho_inner = zero
  temp_inner = zero
  i_inner = 1
  do_particles = .false.
  use_quad = .false.
  nquad = 2

  ! Read namelists
  open(newunit=untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(untin)

  if ( len_trim(chimera_fname) == 0 ) then
    call bl_error("must specify string for chimera_fname")
  end if

  if ( eos_input /= eos_input_rt .and. &
  &    eos_input /= eos_input_rp ) then
    call bl_error("invalid value for eos_input")
  end if

  if ( interp_method /= 1 .and. interp_method /= 2 ) then
    call bl_error("invalid value for interp_method")
  end if

  ! load chimera model
  call open_chimera_file(chimera_fname)
  call read_chimera_file

  ! load mesa model
  if ( len_trim(mesa_fname) > 0 ) then
    call read_mesa_file(mesa_fname)
  end if

  ! set the max radius to use chimera data
  ! if r > max_radius, use mesa data
  if ( max_radius <= zero ) then
    ! by default, use the volume-center of outer-most chimera zone (excluding ghost zones)
    ! this obviates the need for any extrapolation of chimera data
    max_radius = (three * volx_c_chim(imax_chim) )**third
!   imax_radius = imax_chim
! else
!   volr = third * max_radius**3
!   imax_radius = locate( volr, imax_chim, volx_c_chim ) - 1
  end if

  ! largest radius of castro grid in spherical coordinates
  probhi_r = sqrt(max(abs(problo(2)),abs(probhi(2)))**2+probhi(1)**2)

  ! integrate mass from chimera data on castro grid
  mass_chim = zero
  vol_chim = zero
  do j = jmin_chim, jmax_chim
    do i = imin_chim, imax_chim
      if ( x_e_chim(i) < max_radius ) then
        rcyl = x_c_chim(i) * sin( y_c_chim(j) )
        zcyl = x_c_chim(i) * cos( y_c_chim(j) )
        if ( rcyl > problo(1) .and. rcyl < probhi(1) .and. zcyl > problo(2) .and. zcyl < probhi(2) ) then
          if ( x_e_chim(i+1) <= max_radius ) then
            mass_chim = mass_chim + sum(dmass_e_chim(i,j,kmin_chim:kmax_chim))
            vol_chim = vol_chim + sum(dvol_e_chim(i,j,kmin_chim:kmax_chim))
          else
            dr = x_e_chim(i+1) - max_radius
            dvolr = dr * ( x_e_chim(i) * max_radius + dr * dr * third )
            dvol  = sum(domega_chim(j,kmin_chim:kmax_chim)) * ( dvolx_e_chim(i) - dvolr )

            vol_chim = vol_chim + dvol
            mass_chim = mass_chim + sum( rho_c_chim(i,j,kmin_chim:kmax_chim) * domega_chim(j,kmin_chim:kmax_chim) ) * &
            &                       ( dvolx_e_chim(i) - dvolr )
          end if
        end if
      end if
    end do
  end do

  ! integrate mass from mesa data on castro grid
  mass_mesa = zero
  vol_mesa = zero
  do i = 1, imax_mesa
    if ( x_e_mesa(i+1) > max_radius .and. &
    &    x_e_mesa(i)   < probhi_r ) then

      r = half * ( x_e_mesa(i) + x_e_mesa(i+1) )

      tmp1 = min( one, max( -one, probhi(2)/r ) )
      tmp2 = min( one, max( -one, problo(2)/r ) )
      tmp3 = two * sqrt( one - min( one, probhi(1)/r )**2 )

      domega = two * m_pi * ( tmp1 - tmp2 - tmp3 )
      domega = min( four * m_pi, max( zero, domega ) )

      dvolr = zero
      if ( x_e_mesa(i+1) > probhi_r ) then
        dr = x_e_mesa(i+1) - probhi_r
        dvolr = dvolr + dr * ( probhi_r * x_e_mesa(i) + dr * dr * third )
      end if

      if ( x_e_mesa(i) < max_radius ) then
        dr = max_radius - x_e_mesa(i) 
        dvolr = dvolr + dr * ( x_e_mesa(i) * max_radius + dr * dr * third )
      end if

      dvol = domega * ( dvolx_e_mesa(i) - dvolr )

      vol_mesa = vol_mesa + dvol
      mass_mesa = mass_mesa + dvol * rho_c_mesa(i)

    end if
  end do

  ! compare gravitational acceleration from edge of chimera grid with castro
  r = min( probhi_r, max_radius )
  volr = third * r**3
  if ( volr <= volx_c_chim(1) ) then
    i = 0
  else if ( volr >= volx_c_chim(imax_chim) ) then
    i = imax_chim
  else
    i = locate( volr, imax_chim, volx_c_chim ) - 1
  end if
  call interp1d_linear( i, imax_chim, volx_c_chim, gravx_c_avg_chim, volr, gr )

  if (parallel_IOProcessor()) then
    write(*,'(a,2es23.15)') 'average g(r) (chimera)        =',r,gr

    write(*,'(a,2es23.15)') 'total mass   (chimera, mesa)  =',mass_chim,mass_mesa
    write(*,'(a,es23.15)')  'total mass   (chimera + mesa) =',mass_chim+mass_mesa

    write(*,'(a,2es23.15)') 'total volume (chimera, mesa)  =',vol_chim,vol_mesa
    write(*,'(a,es23.15)')  'total volume (chimera + mesa) =',vol_chim+vol_mesa
    write(*,'(a,es23.15)')  'total volume (castro)         =',m_pi * (probhi(1)**2-problo(1)**2) * ( probhi(2)-problo(2) )
  end if

  ! set up quadrature weights and abscissae
  if ( use_quad ) call quad_init( nquad )

  return
end subroutine PROBINIT


! ::: -----------------------------------------------------------
! ::: This routine is called at problem setup time and is used
! ::: to initialize data on each grid.  
! ::: 
! ::: NOTE:  all arrays have one cell of ghost zones surrounding
! :::        the grid interior.  Values in these cells need not
! :::        be set here.
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: level     => amr level of grid
! ::: time      => time at which to init data             
! ::: lo,hi     => index limits of grid interior (cell centered)
! ::: nscal     => number of state components.
! ::: state     <= scalar array
! ::: dx        => cell size
! ::: xlo, xhi  => physical locations of lower left and upper
! :::              right hand corner of grid.  (does not include
! :::		   ghost region).
! ::: -----------------------------------------------------------

subroutine ca_initdata(level,time,lo,hi,nscal, &
                       state,state_l1,state_l2,state_h1,state_h2, &
                       delta,xlo,xhi)

  use bl_error_module
  use bl_types
  use fundamental_constants_module
  use eos_module
  use eos_type_module, only: minye, maxye
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEINT, UFS, UTEMP, UEDEN, UFX, UFA
  use network, only: nspec, naux
  use probdata_module
  use chimera_parser_module
  use mesa_parser_module
  use quadrature_module, only: xquad, wquad, quad_avg

  implicit none

  integer :: level, nscal
  integer :: lo(2), hi(2)
  integer :: state_l1,state_l2,state_h1,state_h2
  double precision :: xlo(2), xhi(2), time, delta(2)
  double precision :: state(state_l1:state_h1,state_l2:state_h2,NVAR)

  ! local variables
  real (dp_t) :: xcen(lo(1):hi(1))
  real (dp_t) :: ycen(lo(2):hi(2))
  real (dp_t) :: r(lo(1):hi(1),lo(2):hi(2))
  real (dp_t) :: theta(lo(1):hi(1),lo(2):hi(2))

  real (dp_t) :: rho_i_chim(lo(1):hi(1),lo(2):hi(2))
  real (dp_t) :: t_i_chim(lo(1):hi(1),lo(2):hi(2))
  real (dp_t) :: p_i_chim(lo(1):hi(1),lo(2):hi(2))
  real (dp_t) :: xn_i_chim(nspec,lo(1):hi(1),lo(2):hi(2))
  real (dp_t) :: u_i_chim(lo(1):hi(1),lo(2):hi(2))
  real (dp_t) :: v_i_chim(lo(1):hi(1),lo(2):hi(2))
  real (dp_t) :: w_i_chim(lo(1):hi(1),lo(2):hi(2))
  real (dp_t) :: ye_i_chim(lo(1):hi(1),lo(2):hi(2))
  real (dp_t) :: a_aux_i_chim(lo(1):hi(1),lo(2):hi(2))
  real (dp_t) :: z_aux_i_chim(lo(1):hi(1),lo(2):hi(2))

  real (dp_t) :: rho_i_mesa(lo(1):hi(1),lo(2):hi(2))
  real (dp_t) :: t_i_mesa(lo(1):hi(1),lo(2):hi(2))
  real (dp_t) :: xn_i_mesa(nspec,lo(1):hi(1),lo(2):hi(2))
  real (dp_t) :: u_i_mesa(lo(1):hi(1),lo(2):hi(2))
  real (dp_t) :: ye_i_mesa(lo(1):hi(1),lo(2):hi(2))
  real (dp_t) :: a_aux_i_mesa(lo(1):hi(1),lo(2):hi(2))
  real (dp_t) :: z_aux_i_mesa(lo(1):hi(1),lo(2):hi(2))

  integer :: i, ii, j, jj, n

  real (dp_t) :: xg, yg
  real (dp_t) :: rg(nquad,nquad)
  real (dp_t) :: tg(nquad,nquad)
  real (dp_t) :: rho_quad(nquad,nquad)
  real (dp_t) :: t_quad(nquad,nquad)
  real (dp_t) :: p_quad(nquad,nquad)
  real (dp_t) :: u_quad(nquad,nquad)
  real (dp_t) :: v_quad(nquad,nquad)
  real (dp_t) :: w_quad(nquad,nquad)
  real (dp_t) :: xn_quad(nquad,nquad)
  real (dp_t) :: ye_quad(nquad,nquad)
  real (dp_t) :: a_aux_quad(nquad,nquad)
  real (dp_t) :: z_aux_quad(nquad,nquad)

  type (eos_t) :: eos_state
 
  real (dp_t) :: drho
  real (dp_t), parameter :: dsmooth = 1.0e6_dp_t

  ! determine coordinates in r-theta and quadrature points
  do j = lo(2), hi(2)
    ycen(j) = xlo(2) + delta(2)*(dble(j-lo(2)) + half)
    do i = lo(1), hi(1)
      xcen(i) = xlo(1) + delta(1)*(dble(i-lo(1)) + half)
      r(i,j) = sqrt( xcen(i)**2 + ycen(j)**2 )
      if ( r(i,j) <= zero ) then
        theta(i,j) = zero
      else if ( ycen(j) == zero ) then
        theta(i,j) = half * m_pi
      else
        theta(i,j) = atan( xcen(i)/ycen(j) )
      end if
      if ( theta(i,j) < zero ) then
        theta(i,j) = theta(i,j) + m_pi
      end if

    end do
  end do

  if ( use_quad ) then
    do j = lo(2), hi(2)
      do i = lo(1), hi(1)
        do jj = 1, nquad
          yg = ycen(j) + half*delta(2)*xquad(jj)
          do ii = 1, nquad
            xg = xcen(i) + half*delta(1)*xquad(ii)
            rg(ii,jj) = sqrt( xg**2 + yg**2 )
            if ( rg(ii,jj) <= zero ) then
              tg(ii,jj) = zero
            else if ( yg == zero ) then
              tg(ii,jj) = half * m_pi
            else
              tg(ii,jj) = atan( xg/yg )
            end if
            if ( tg(ii,jj) < zero ) then
              tg(ii,jj) = tg(ii,jj) + m_pi
            end if
          end do
        end do

        select case (eos_input)
        case (eos_input_rt)
          call interp2d_chimera( rg, tg, rho_c_chim(:,:,1), rho_quad )
          call interp2d_chimera( rg, tg, t_c_chim(:,:,1), t_quad )
        case (eos_input_rp)
          call interp2d_chimera( rg, tg, rho_c_chim(:,:,1), rho_quad )
          call interp2d_chimera( rg, tg, p_c_chim(:,:,1), p_quad )
        case default
          call interp2d_chimera( rg, tg, rho_c_chim(:,:,1), rho_i_chim )
          call interp2d_chimera( rg, tg, t_c_chim(:,:,1), t_i_chim )
        end select
        call interp2d_chimera( rg, tg, u_c_chim(:,:,1), u_quad )
        call interp2d_chimera( rg, tg, v_c_chim(:,:,1), v_quad )
        call interp2d_chimera( rg, tg, w_c_chim(:,:,1), w_quad )
        do n = 1, nspec
          call interp2d_chimera( rg, tg, xn_c_chim(n,:,:,1), xn_quad )
          xn_i_chim(n,i,j) = quad_avg( wquad, xn_quad )
        end do
        if ( naux == 1 ) then
          call interp2d_chimera( rg, tg, ye_c_chim(:,:,1), ye_quad )
        else if ( naux == 2 ) then
          call interp2d_chimera( rg, tg, a_aux_c_chim(:,:,1), a_aux_quad )
          call interp2d_chimera( rg, tg, z_aux_c_chim(:,:,1), z_aux_quad )
        else if ( naux == 3 ) then
          call interp2d_chimera( rg, tg, ye_c_chim(:,:,1), ye_quad )
          call interp2d_chimera( rg, tg, a_aux_c_chim(:,:,1), a_aux_quad )
          call interp2d_chimera( rg, tg, z_aux_c_chim(:,:,1), z_aux_quad )
        end if
        rho_i_chim(i,j) = quad_avg( wquad, rho_quad )
        t_i_chim(i,j) = quad_avg( wquad, t_quad )
        p_i_chim(i,j) = quad_avg( wquad, p_quad )
        u_i_chim(i,j) = quad_avg( wquad, u_quad )
        v_i_chim(i,j) = quad_avg( wquad, v_quad )
        w_i_chim(i,j) = quad_avg( wquad, w_quad )
        ye_i_chim(i,j) = quad_avg( wquad, ye_quad )
        a_aux_i_chim(i,j) = quad_avg( wquad, a_aux_quad )
        z_aux_i_chim(i,j) = quad_avg( wquad, z_aux_quad )
      end do
    end do
  else
    select case (eos_input)
      case (eos_input_rt)
        call interp2d_chimera( r, theta, rho_c_chim(:,:,1), rho_i_chim )
        call interp2d_chimera( r, theta, t_c_chim(:,:,1), t_i_chim )
      case (eos_input_rp)
        call interp2d_chimera( r, theta, rho_c_chim(:,:,1), rho_i_chim )
        call interp2d_chimera( r, theta, p_c_chim(:,:,1), p_i_chim )
      case default
        call interp2d_chimera( r, theta, rho_c_chim(:,:,1), rho_i_chim )
        call interp2d_chimera( r, theta, t_c_chim(:,:,1), t_i_chim )
    end select
    do n = 1, nspec
      call interp2d_chimera( r, theta, xn_c_chim(n,:,:,1), xn_i_chim(n,:,:) )
    end do
    call interp2d_chimera( r, theta, u_c_chim(:,:,1), u_i_chim )
    call interp2d_chimera( r, theta, v_c_chim(:,:,1), v_i_chim )
    call interp2d_chimera( r, theta, w_c_chim(:,:,1), w_i_chim )
    if ( naux == 1 ) then
      call interp2d_chimera( r, theta, ye_c_chim(:,:,1), ye_i_chim )
    else if ( naux == 2 ) then
      call interp2d_chimera( r, theta, a_aux_c_chim(:,:,1), a_aux_i_chim )
      call interp2d_chimera( r, theta, z_aux_c_chim(:,:,1), z_aux_i_chim )
    else if ( naux == 3 ) then
      call interp2d_chimera( r, theta, ye_c_chim(:,:,1), ye_i_chim )
      call interp2d_chimera( r, theta, a_aux_c_chim(:,:,1), a_aux_i_chim )
      call interp2d_chimera( r, theta, z_aux_c_chim(:,:,1), z_aux_i_chim )
    end if
  end if

  if ( len_trim(mesa_fname) > 0 ) then
    call interp2d_mesa( r, rho_c_mesa, rho_i_mesa )
    call interp2d_mesa( r, t_c_mesa, t_i_mesa )
    do n = 1, nspec
      call interp2d_mesa( r, xn_c_mesa(:,n), xn_i_mesa(n,:,:) )
    end do
    call interp2d_mesa( r, u_c_mesa, u_i_mesa )
    if ( naux == 1 ) then
      call interp2d_mesa( r, ye_c_mesa, ye_i_mesa )
    else if ( naux == 2 ) then
      call interp2d_mesa( r, a_aux_c_mesa, a_aux_i_mesa )
      call interp2d_mesa( r, z_aux_c_mesa, z_aux_i_mesa )
    else if ( naux == 3 ) then
      call interp2d_mesa( r, ye_c_mesa, ye_i_mesa )
      call interp2d_mesa( r, a_aux_c_mesa, a_aux_i_mesa )
      call interp2d_mesa( r, z_aux_c_mesa, z_aux_i_mesa )
    end if
  end if

  do j = lo(2), hi(2)
    do i = lo(1), hi(1)

      xn_i_chim(:,i,j) = xn_i_chim(:,i,j) / sum( xn_i_chim(:,i,j) )
      xn_i_mesa(:,i,j) = xn_i_mesa(:,i,j) / sum( xn_i_mesa(:,i,j) )

      if ( r(i,j) <= max_radius .or. len_trim(mesa_fname) == 0 ) then

        if ( r(i,j) < radius_inner ) then
!         drho = half * ( rho_i_chim(i,j) - rho_inner ) * ( one + tanh( (radius_inner - r(i,j))/dsmooth) )
!         state(i,j,URHO) = rho_i_chim(i,j) + drho
          state(i,j,URHO) = rho_inner
          eos_state%rho = state(i,j,URHO)
          eos_state%T = t_i_chim(i,j)
          eos_state%p = p_i_chim(i,j)
        else
          select case (eos_input)
            case (eos_input_rt)
              eos_state%rho = rho_i_chim(i,j)
              eos_state%T = t_i_chim(i,j)
            case (eos_input_rp)
              eos_state%rho = rho_i_chim(i,j)
              eos_state%p = p_i_chim(i,j)
            case default
              eos_state%rho = rho_i_chim(i,j)
              eos_state%T = t_i_chim(i,j)
          end select
        end if
        eos_state%xn(:) = xn_i_chim(:,i,j)
        if ( naux == 1 ) then
          eos_state%aux = min( maxye, max( minye, ye_i_chim(i,j) ) )
        else if ( naux == 2 ) then
          eos_state%aux(1) = a_aux_i_chim(i,j)
          eos_state%aux(2) = z_aux_i_chim(i,j)
        else if ( naux == 3 ) then
          eos_state%aux(1) = min( maxye, max( minye, ye_i_chim(i,j) ) )
          eos_state%aux(2) = a_aux_i_chim(i,j)
          eos_state%aux(3) = z_aux_i_chim(i,j)
        end if

        call eos(eos_input, eos_state)

        state(i,j,UMX) = u_i_chim(i,j) * sin( theta(i,j) ) + v_i_chim(i,j) * cos( theta(i,j) )
        state(i,j,UMY) = u_i_chim(i,j) * cos( theta(i,j) ) - v_i_chim(i,j) * sin( theta(i,j) )
        state(i,j,UMZ) = w_i_chim(i,j)
        state(i,j,UFS:UFS+nspec-1) = xn_i_chim(:,i,j)
        if ( naux == 1 ) then
          state(i,j,UFX) = min( maxye, max( minye, ye_i_chim(i,j) ) )
        else if ( naux == 2 ) then
          state(i,j,UFX)   = a_aux_i_chim(i,j)
          state(i,j,UFX+1) = z_aux_i_chim(i,j)
        else if ( naux == 3 ) then
          state(i,j,UFX) = min( maxye, max( minye, ye_i_chim(i,j) ) )
          state(i,j,UFX+1) = a_aux_i_chim(i,j)
          state(i,j,UFX+2) = z_aux_i_chim(i,j)
        end if

        select case (eos_input)
          case (eos_input_rt)
            state(i,j,URHO) = rho_i_chim(i,j)
            state(i,j,UTEMP) = t_i_chim(i,j)
            state(i,j,UEINT) = eos_state%e
          case (eos_input_rp)
            state(i,j,URHO) = rho_i_chim(i,j)
            state(i,j,UTEMP) = eos_state%T
            state(i,j,UEINT) = eos_state%e
          case default
            state(i,j,URHO) = rho_i_chim(i,j)
            state(i,j,UTEMP) = t_i_chim(i,j)
            state(i,j,UEINT) = eos_state%e
        end select

      else

        eos_state%rho = rho_i_mesa(i,j)
        eos_state%T = t_i_mesa(i,j)
        eos_state%xn(:) = xn_i_mesa(:,i,j)
        if ( naux == 1 ) then
          eos_state%aux = min( maxye, max( minye, ye_i_mesa(i,j) ) )
        else if ( naux == 2 ) then
          eos_state%aux(1) = a_aux_i_mesa(i,j)
          eos_state%aux(2) = z_aux_i_mesa(i,j)
        else if ( naux == 3 ) then
          eos_state%aux(1) = min( maxye, max( minye, ye_i_mesa(i,j) ) )
          eos_state%aux(2) = a_aux_i_mesa(i,j)
          eos_state%aux(3) = z_aux_i_mesa(i,j)
        end if

        call eos(eos_input_rt, eos_state)

        state(i,j,UMX) = u_i_mesa(i,j) * sin( theta(i,j) )
        state(i,j,UMY) = u_i_mesa(i,j) * cos( theta(i,j) )
        state(i,j,UMZ) = zero
        state(i,j,UFS:UFS+nspec-1) = xn_i_mesa(:,i,j)
        if ( naux == 1 ) then
          state(i,j,UFX) = min( maxye, max( minye, ye_i_mesa(i,j) ) )
        else if ( naux == 2 ) then
          state(i,j,UFX)   = a_aux_i_mesa(i,j)
          state(i,j,UFX+1) = z_aux_i_mesa(i,j)
        else if ( naux == 3 ) then
          state(i,j,UFX) = min( maxye, max( minye, ye_i_mesa(i,j) ) )
          state(i,j,UFX+1) = a_aux_i_mesa(i,j)
          state(i,j,UFX+2) = z_aux_i_mesa(i,j)
        end if

        state(i,j,URHO) = rho_i_mesa(i,j)
        state(i,j,UTEMP) = t_i_mesa(i,j)
        state(i,j,UEINT) = eos_state%e

      end if

    end do
  end do

  do j = lo(2), hi(2)
    do i = lo(1), hi(1)
      state(i,j,UEINT)   = state(i,j,URHO) * state(i,j,UEINT)
      state(i,j,UEDEN)   = state(i,j,UEINT) + state(i,j,URHO)*sum( half*state(i,j,UMX:UMZ)**2 )
      state(i,j,UMX:UMZ) = state(i,j,URHO) * state(i,j,UMX:UMZ)
      state(i,j,UFS:UFS+nspec-1) = state(i,j,URHO) * state(i,j,UFS:UFS+nspec-1)
      if ( naux > 0 ) then
        state(i,j,UFX:UFX+naux-1) = state(i,j,URHO) * state(i,j,UFX:UFX+naux-1)
      end if
    end do
  end do

  return
end subroutine ca_initdata

