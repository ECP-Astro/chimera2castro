subroutine PROBINIT (init,name,namlen,problo,probhi)

  use probdata_module
  use chimera_parser_module
  use mesa_parser_module
  use bl_constants_module
  use bl_error_module
  use bl_types
  use eos_module
  use parallel, only: parallel_IOProcessor
  use prob_params_module, only : center
  use interpolate_module, only: locate
  use model_interp_module, only: interp1d_linear

  implicit none
  integer :: init, namlen
  integer :: name(namlen)
  real (dp_t) :: problo(2), probhi(2)

  integer :: untin,i,j,k,dir
  real (dp_t) :: probhi_r, rcyl, zcyl, r, dr, dvol, volr, dvolr, gr
  real (dp_t) :: mass_chim, mass_mesa, vol_chim, vol_mesa
  real (dp_t) :: tmp1, tmp2, tmp3, domega
  real (dp_t) :: wquad_sum

  namelist /fortin/ model_name
  namelist /fortin/ mesa_name
  namelist /fortin/ model_eos_input
  namelist /fortin/ model_interp_method
  namelist /fortin/ min_radius
  namelist /fortin/ max_radius
  namelist /fortin/ do_particles
  namelist /fortin/ use_quad

  !
  !     Build "probin" filename -- the name of file containing fortin namelist.
  !     
  integer, parameter :: maxlen = 127
  character(maxlen) :: probin
  character(maxlen) :: model
  integer :: ipp, ierr, ipp1, n

  if (namlen .gt. maxlen) call bl_error("probin file name too long")

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! set namelist defaults
  model_eos_input = eos_input_rt
  model_interp_method = 1
  min_radius = zero
  max_radius = zero
  model_name = ""
  mesa_name = ""
  do_particles = .false.
  use_quad = .false.

  ! Read namelists
  untin = 9 
  open(untin,file=probin(1:namlen),form='formatted',status='old')
  read(untin,fortin)
  close(unit=untin)

  if ( trim(model_name) == "" ) then
    call bl_error("must specify string for model_name")
  end if

  if ( model_eos_input /= eos_input_rt .and. &
  &    model_eos_input /= eos_input_rp .and. &
  &    model_eos_input /= eos_input_re .and. &
  &    model_eos_input /= eos_input_ps ) then
    call bl_error("invalid value for model_eos_input")
  end if

  if ( model_interp_method /= 1 .and. model_interp_method /= 2 ) then
    call bl_error("invalid value for model_interp_method")
  end if

  ! open initial model
  call open_chimera_file(model_name)

  ! read initial model
  call read_chimera_file

  if ( .not. trim(mesa_name) == "" ) then
    call read_mesa_file(mesa_name)
  end if

  if ( max_radius <= zero ) then
    max_radius = rad_edge_in(imax_in+1)
  end if

  ! compute largest radius in spherical coordinates
  probhi_r = sqrt(max(abs(problo(2)),abs(probhi(2)))**2+probhi(1)**2)

  ! integrate mass from CHIMERA model on CASTRO grid
  mass_chim = zero
  vol_chim = zero
  do j = jmin_in, jmax_in
    do i = imin_in, imax_in
      if ( rad_edge_in(i) < max_radius ) then
        rcyl = rad_cntr_in(i) * sin( theta_cntr_in(j) )
        zcyl = rad_cntr_in(i) * cos( theta_cntr_in(j) )
        if ( rcyl > problo(1) .and. rcyl < probhi(1) .and. zcyl > problo(2) .and. zcyl < probhi(2) ) then
          if ( rad_edge_in(i+1) <= max_radius ) then
            mass_chim = mass_chim + sum(zone_mass_in(i,j,kmin_in:kmax_in))
            vol_chim = vol_chim + sum(volume_in(i,j,kmin_in:kmax_in))
          else
            dr = rad_edge_in(i+1) - max_radius
            dvolr = dr * ( rad_edge_in(i) * max_radius + dr * dr * third )
            dvol  = sum(domega_in(j,kmin_in:kmax_in)) * ( dvol_rad_in(i) - dvolr )

            vol_chim = vol_chim + dvol
            mass_chim = mass_chim + sum( dens_in(i,j,kmin_in:kmax_in) * domega_in(j,kmin_in:kmax_in) ) * ( dvol_rad_in(i) - dvolr )
          end if
        end if
      end if
    end do
  end do

  ! integrate mass from MESA model on CASTRO grid
  mass_mesa = zero
  vol_mesa = zero
  do i = 1, nx_mesa_in
    if ( rad_edge_mesa_in(i+1) > max_radius .and. &
    &    rad_edge_mesa_in(i)   < probhi_r ) then

      r = half * ( rad_edge_mesa_in(i) + rad_edge_mesa_in(i+1) )

      tmp1 = min( one, max( -one, probhi(2)/r ) )
      tmp2 = min( one, max( -one, problo(2)/r ) )
      tmp3 = two * sqrt( one - min( one, probhi(1)/r )**2 )

      domega = two * m_pi * ( tmp1 - tmp2 - tmp3 )
      domega = min( four * m_pi, max( zero, domega ) )

      dvolr = zero
      if ( rad_edge_mesa_in(i+1) > probhi_r ) then
        dr = rad_edge_mesa_in(i+1) - probhi_r
        dvolr = dvolr + dr * ( probhi_r * rad_edge_mesa_in(i) + dr * dr * third )
      end if

      if ( rad_edge_mesa_in(i) < max_radius ) then
        dr = max_radius - rad_edge_mesa_in(i) 
        dvolr = dvolr + dr * ( rad_edge_mesa_in(i) * max_radius + dr * dr * third )
      end if

      dvol = domega * ( dvol_rad_mesa_in(i) - dvolr )

      vol_mesa = vol_mesa + dvol
      mass_mesa = mass_mesa + dvol * dens_mesa_in(i)

    end if
  end do

  ! Get the gravitational acceleration to compare with CASTRO
  r = min( probhi_r, max_radius )
  volr = third * r**3
  if ( volr <= vol_rad_cntr_in(1) ) then
    i = 0
  else if ( volr >= vol_rad_cntr_in(imax_in) ) then
    i = imax_in
  else
    i = locate( volr, imax_in, vol_rad_cntr_in ) - 1
  end if
  call interp1d_linear( i, imax_in, vol_rad_cntr_in, grad_avg_in, volr, gr )

  if (parallel_IOProcessor()) then
    write(*,'(a,2es23.15)') 'average g(r) (chimera)        =',r,gr

    write(*,'(a,2es23.15)') 'total mass   (chimera, mesa)  =',mass_chim,mass_mesa
    write(*,'(a,es23.15)')  'total mass   (chimera + mesa) =',mass_chim+mass_mesa

    write(*,'(a,2es23.15)') 'total volume (chimera, mesa)  =',vol_chim,vol_mesa
    write(*,'(a,es23.15)')  'total volume (chimera + mesa) =',vol_chim+vol_mesa
    write(*,'(a,es23.15)')  'total volume (castro)         =',m_pi * (probhi(1)**2-problo(1)**2) * ( probhi(2)-problo(2) )
  end if

  center(1) = zero
  center(2) = zero

  nquad = 2
  if ( .not. allocated( xquad ) ) allocate (xquad(nquad))
  if ( .not. allocated( wquad ) ) allocate (wquad(nquad))

  if ( nquad == 2 ) then
    xquad(1) = -one / sqrt(three)
    xquad(2) =  one / sqrt(three)
    wquad(1) = one
    wquad(2) = one
    norm_quad = fourth
  else
    call gquad( nquad, xquad, wquad, nquad )
    wquad_sum = zero
    do j = 1, nquad
      do i = 1, nquad
        wquad_sum = wquad_sum + wquad(i) * wquad(j)
      end do
    end do
    norm_quad = one / wquad_sum
  end if

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
  use chimera_parser_module
  use mesa_parser_module
  use fundamental_constants_module
  use eos_module
  use eos_type_module, only: minye, maxye
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEINT, UFS, UTEMP, UEDEN, UFX, UFA
  use network, only: nspec, naux
  use probdata_module

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
  real (dp_t) :: rho_chim(lo(1):hi(1),lo(2):hi(2))
  real (dp_t) :: temp_chim(lo(1):hi(1),lo(2):hi(2))
  real (dp_t) :: pressure_chim(lo(1):hi(1),lo(2):hi(2))
  real (dp_t) :: eint_chim(lo(1):hi(1),lo(2):hi(2))
  real (dp_t) :: entropy_chim(lo(1):hi(1),lo(2):hi(2))
  real (dp_t) :: xn_chim(lo(1):hi(1),lo(2):hi(2),nspec)
  real (dp_t) :: vrad_chim(lo(1):hi(1),lo(2):hi(2))
  real (dp_t) :: vtheta_chim(lo(1):hi(1),lo(2):hi(2))
  real (dp_t) :: vphi_chim(lo(1):hi(1),lo(2):hi(2))
  real (dp_t) :: ye_chim(lo(1):hi(1),lo(2):hi(2))
  real (dp_t) :: rho_mesa(lo(1):hi(1),lo(2):hi(2))
  real (dp_t) :: temp_mesa(lo(1):hi(1),lo(2):hi(2))
  real (dp_t) :: xn_mesa(lo(1):hi(1),lo(2):hi(2),nspec)
  real (dp_t) :: vrad_mesa(lo(1):hi(1),lo(2):hi(2))
  real (dp_t) :: ye_mesa(lo(1):hi(1),lo(2):hi(2))
  integer :: i, ii, j, jj, n

  real (dp_t) :: xg, yg
  real (dp_t) :: rg(nquad,nquad,lo(1):hi(1),lo(2):hi(2))
  real (dp_t) :: tg(nquad,nquad,lo(1):hi(1),lo(2):hi(2))

  type (eos_t) :: eos_state

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

      do jj = 1, nquad
        yg = ycen(j) + half*delta(2)*xquad(jj)
        do ii = 1, nquad
          xg = xcen(i) + half*delta(1)*xquad(ii)
          rg(ii,jj,i,j) = sqrt( xg**2 + yg**2 )
          if ( rg(ii,jj,i,j) <= zero ) then
            tg(ii,jj,i,j) = zero
          else if ( yg == zero ) then
            tg(ii,jj,i,j) = half * m_pi
          else
            tg(ii,jj,i,j) = atan( xg/yg )
          end if
          if ( tg(ii,jj,i,j) < zero ) then
            tg(ii,jj,i,j) = tg(ii,jj,i,j) + m_pi
          end if
        end do
      end do
    end do
  end do

  if ( use_quad ) then
    select case (model_eos_input)
      case (eos_input_rt)
        call interp2d_quad_chimera( rg, tg, dens_in(:,:,1), rho_chim, model_interp_method )
        call interp2d_quad_chimera( rg, tg, temp_in(:,:,1), temp_chim, model_interp_method )
      case (eos_input_rp)
        call interp2d_quad_chimera( rg, tg, dens_in(:,:,1), rho_chim, model_interp_method )
        call interp2d_quad_chimera( rg, tg, pres_in(:,:,1), pressure_chim, model_interp_method )
      case (eos_input_re)
        call interp2d_quad_chimera( rg, tg, dens_in(:,:,1), rho_chim, model_interp_method )
        call interp2d_quad_chimera( rg, tg, eint_in(:,:,1), eint_chim, model_interp_method )
      case (eos_input_ps)
        call interp2d_quad_chimera( rg, tg, enpy_in(:,:,1), entropy_chim, model_interp_method )
        call interp2d_quad_chimera( rg, tg, pres_in(:,:,1), pressure_chim, model_interp_method )
      case default
        call interp2d_quad_chimera( rg, tg, dens_in(:,:,1), rho_chim, model_interp_method )
        call interp2d_quad_chimera( rg, tg, temp_in(:,:,1), temp_chim, model_interp_method )
    end select
    do n = 1, nspec
      call interp2d_quad_chimera( rg, tg, xn_in(n,:,:,1), xn_chim(:,:,n), model_interp_method )
    end do
    call interp2d_quad_chimera( rg, tg, vrad_in(:,:,1), vrad_chim, model_interp_method )
    call interp2d_quad_chimera( rg, tg, vtheta_in(:,:,1), vtheta_chim, model_interp_method )
    call interp2d_quad_chimera( rg, tg, vphi_in(:,:,1), vphi_chim, model_interp_method )
    if ( naux > 0 ) then
      call interp2d_quad_chimera( rg, tg, ye_in(:,:,1), ye_chim, model_interp_method )
    end if
  else
    select case (model_eos_input)
      case (eos_input_rt)
        call interp2d_chimera( r, theta, dens_in(:,:,1), rho_chim, model_interp_method )
        call interp2d_chimera( r, theta, temp_in(:,:,1), temp_chim, model_interp_method )
      case (eos_input_rp)
        call interp2d_chimera( r, theta, dens_in(:,:,1), rho_chim, model_interp_method )
        call interp2d_chimera( r, theta, pres_in(:,:,1), pressure_chim, model_interp_method )
      case (eos_input_re)
        call interp2d_chimera( r, theta, dens_in(:,:,1), rho_chim, model_interp_method )
        call interp2d_chimera( r, theta, eint_in(:,:,1), eint_chim, model_interp_method )
      case (eos_input_ps)
        call interp2d_chimera( r, theta, enpy_in(:,:,1), entropy_chim, model_interp_method )
        call interp2d_chimera( r, theta, pres_in(:,:,1), pressure_chim, model_interp_method )
      case default
        call interp2d_chimera( r, theta, dens_in(:,:,1), rho_chim, model_interp_method )
        call interp2d_chimera( r, theta, temp_in(:,:,1), temp_chim, model_interp_method )
    end select
    do n = 1, nspec
      call interp2d_chimera( r, theta, xn_in(n,:,:,1), xn_chim(:,:,n), model_interp_method )
    end do
    call interp2d_chimera( r, theta, vrad_in(:,:,1), vrad_chim, model_interp_method )
    call interp2d_chimera( r, theta, vtheta_in(:,:,1), vtheta_chim, model_interp_method )
    call interp2d_chimera( r, theta, vphi_in(:,:,1), vphi_chim, model_interp_method )
    if ( naux > 0 ) then
      call interp2d_chimera( r, theta, ye_in(:,:,1), ye_chim, model_interp_method )
    end if
  end if


  if ( .not. trim(mesa_name) == "" ) then
    call interp2d_mesa( r, dens_mesa_in, rho_mesa, model_interp_method )
    call interp2d_mesa( r, temp_mesa_in, temp_mesa, model_interp_method )
    do n = 1, nspec
      call interp2d_mesa( r, xn_mesa_in(:,n), xn_mesa(:,:,n), model_interp_method )
    end do
    call interp2d_mesa( r, vrad_mesa_in, vrad_mesa, model_interp_method )
    if ( naux > 0 ) then
      call interp2d_mesa( r, ye_mesa_in, ye_mesa, model_interp_method )
    end if
  end if

  do j = lo(2), hi(2)
    do i = lo(1), hi(1)

      xn_chim(i,j,:) = xn_chim(i,j,:) / sum( xn_chim(i,j,:) )
      xn_mesa(i,j,:) = xn_mesa(i,j,:) / sum( xn_mesa(i,j,:) )

      if ( r(i,j) <= max_radius .or. trim(mesa_name) == "" ) then

        select case (model_eos_input)
          case (eos_input_rt)
            eos_state%rho = rho_chim(i,j)
            eos_state%T = temp_chim(i,j)
          case (eos_input_rp)
            eos_state%rho = rho_chim(i,j)
            eos_state%p = pressure_chim(i,j)
          case (eos_input_re)
            eos_state%rho = rho_chim(i,j)
            eos_state%e = eint_chim(i,j)
          case (eos_input_ps)
            eos_state%p = pressure_chim(i,j)
            eos_state%s = entropy_chim(i,j)
          case default
            eos_state%rho = rho_chim(i,j)
            eos_state%T = temp_chim(i,j)
        end select
        eos_state%xn(:) = xn_chim(i,j,:)
        if ( naux > 0 ) then
          eos_state%aux = min( maxye, max( minye, ye_chim(i,j) ) )
        end if

        call eos(model_eos_input, eos_state)

        state(i,j,UMX) = vrad_chim(i,j) * sin( theta(i,j) ) + vtheta_chim(i,j) * cos( theta(i,j) )
        state(i,j,UMY) = vrad_chim(i,j) * cos( theta(i,j) ) - vtheta_chim(i,j) * sin( theta(i,j) )
        state(i,j,UMZ) = vphi_chim(i,j)
        state(i,j,UFS:UFS+nspec-1) = xn_chim(i,j,:)
        if ( naux > 0 ) then
          state(i,j,UFX) = min( maxye, max( minye, ye_chim(i,j) ) )
        end if

        select case (model_eos_input)
          case (eos_input_rt)
            state(i,j,URHO) = rho_chim(i,j)
            state(i,j,UTEMP) = temp_chim(i,j)
            state(i,j,UEINT) = eos_state%e
          case (eos_input_rp)
            state(i,j,URHO) = rho_chim(i,j)
            state(i,j,UTEMP) = eos_state%T
            state(i,j,UEINT) = eos_state%e
          case (eos_input_re)
            state(i,j,URHO) = rho_chim(i,j)
            state(i,j,UTEMP) = eos_state%T
            state(i,j,UEINT) = eint_chim(i,j)
          case (eos_input_ps)
            state(i,j,URHO) = eos_state%rho
            state(i,j,UTEMP) = eos_state%T
            state(i,j,UEINT) = eos_state%e
          case default
            state(i,j,URHO) = rho_chim(i,j)
            state(i,j,UTEMP) = temp_chim(i,j)
            state(i,j,UEINT) = eos_state%e
        end select

      else

        eos_state%rho = rho_mesa(i,j)
        eos_state%T = temp_mesa(i,j)
        eos_state%xn(:) = xn_mesa(i,j,:)
        if ( naux > 0 ) then
          eos_state%aux = min( maxye, max( minye, ye_mesa(i,j) ) )
        end if
        call eos(eos_input_rt, eos_state)

        state(i,j,UMX) = vrad_mesa(i,j) * sin( theta(i,j) )
        state(i,j,UMY) = vrad_mesa(i,j) * cos( theta(i,j) )
        state(i,j,UMZ) = zero
        state(i,j,UFS:UFS+nspec-1) = xn_mesa(i,j,:)
        if ( naux > 0 ) then
          state(i,j,UFX) = min( maxye, max( minye, ye_mesa(i,j) ) )
        end if

        state(i,j,URHO) = rho_mesa(i,j)
        state(i,j,UTEMP) = temp_mesa(i,j)
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
        state(i,j,UFX)   = state(i,j,URHO) * state(i,j,UFX)
      end if
    end do
  end do

  return
end subroutine ca_initdata

