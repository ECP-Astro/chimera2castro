subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use probdata_module
  use chimera_parser_module
  use star_parser_module
  use amrex_constants_module
  use amrex_error_module
  use amrex_fort_module, only: rt => amrex_real
  use eos_module
  use eos_type_module
  use amrex_paralleldescriptor_module, only: amrex_pd_ioprocessor
  use prob_params_module, only: center
  use interpolate_module, only: locate
  use model_interp_module, only: interp1d_linear
  use quadrature_module, only: quad_init

  implicit none
  integer :: init, namlen
  integer :: name(namlen)
  real (rt) :: problo(3), probhi(3)

  integer :: untin,i,j,k,dir
  real (rt) :: probhi_r, rcyl, zcyl, r, dr, dvol, volr, dvolr, gr
  real (rt) :: mass_chim, mass_star, vol_chim, vol_star
  real (rt) :: tmp1, tmp2, tmp3, domega

  namelist /fortin/ chimera_fname
  namelist /fortin/ star_fname
  namelist /fortin/ star_type
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
  center(1) = ZERO
  center(2) = HALF*(problo(2)+probhi(2))

  if (namlen .gt. maxlen) call bl_error("probin file name too long")

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! set namelist defaults
  chimera_fname = ""
  star_fname = ""
  interp_method = 1
  eos_input = eos_input_rp
  max_radius = ZERO
  radius_inner = ZERO
  rho_inner = ZERO
  temp_inner = ZERO
  i_inner = 1
  do_particles = .false.
  use_quad = .false.
  nquad = 2
  star_type = 2

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

  ! load star model
  if ( len_trim(star_fname) > 0 ) then
    call read_star_file(star_fname)
  end if

  ! set the max radius to use chimera data
  ! if r > max_radius, use star data
  if ( max_radius <= ZERO ) then
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
  mass_chim = ZERO
  vol_chim = ZERO
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

  ! integrate mass from star data on castro grid
  mass_star = ZERO
  vol_star = ZERO
  do i = 1, imax_star
    if ( x_e_star(i+1) > max_radius .and. &
    &    x_e_star(i)   < probhi_r ) then

      r = HALF * ( x_e_star(i) + x_e_star(i+1) )

      tmp1 = min( one, max( -one, probhi(2)/r ) )
      tmp2 = min( one, max( -one, problo(2)/r ) )
      tmp3 = two * sqrt( one - min( one, probhi(1)/r )**2 )

      domega = two * m_pi * ( tmp1 - tmp2 - tmp3 )
      domega = min( four * m_pi, max( ZERO, domega ) )

      dvolr = ZERO
      if ( x_e_star(i+1) > probhi_r ) then
        dr = x_e_star(i+1) - probhi_r
        dvolr = dvolr + dr * ( probhi_r * x_e_star(i) + dr * dr * third )
      end if

      if ( x_e_star(i) < max_radius ) then
        dr = max_radius - x_e_star(i) 
        dvolr = dvolr + dr * ( x_e_star(i) * max_radius + dr * dr * third )
      end if

      dvol = domega * ( dvolx_e_star(i) - dvolr )

      vol_star = vol_star + dvol
      mass_star = mass_star + dvol * rho_c_star(i)

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

  if (amrex_pd_ioprocessor()) then
    write(*,'(a,2es23.15)') 'average g(r) (chimera)        =',r,gr

    write(*,'(a,2es23.15)') 'total mass   (chimera, star)  =',mass_chim,mass_star
    write(*,'(a,es23.15)')  'total mass   (chimera + star) =',mass_chim+mass_star

    write(*,'(a,2es23.15)') 'total volume (chimera, star)  =',vol_chim,vol_star
    write(*,'(a,es23.15)')  'total volume (chimera + star) =',vol_chim+vol_star
    write(*,'(a,es23.15)')  'total volume (castro)         =',m_pi * (probhi(1)**2-problo(1)**2) * ( probhi(2)-problo(2) )
  end if

  ! set up quadrature weights and abscissae
  if ( use_quad ) call quad_init( nquad )

  return
end subroutine amrex_probinit


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
                       state, state_lo, state_hi, &
                       delta,xlo,xhi)

  use amrex_error_module
  use amrex_fort_module, only: rt => amrex_real
  use fundamental_constants_module
  use amrex_constants_module
  use eos_module
  use eos_type_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEINT, UFS, UTEMP, UEDEN, UFX, UFA
  use network
  use probdata_module
  use chimera_parser_module
  use star_parser_module
  use quadrature_module, only: xquad, wquad, quad_avg

  implicit none

  integer :: level, nscal
  integer :: lo(3), hi(3)
  integer :: state_lo(3), state_hi(3)
  double precision :: xlo(3), xhi(3), time, delta(3)
  double precision :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),NVAR)

  ! local variables
  real (rt) :: xcen(lo(1):hi(1))
  real (rt) :: ycen(lo(2):hi(2))
  real (rt) :: r(lo(1):hi(1),lo(2):hi(2))
  real (rt) :: theta(lo(1):hi(1),lo(2):hi(2))

  real (rt) :: rho_i_chim(lo(1):hi(1),lo(2):hi(2))
  real (rt) :: t_i_chim(lo(1):hi(1),lo(2):hi(2))
  real (rt) :: p_i_chim(lo(1):hi(1),lo(2):hi(2))
  real (rt) :: xn_i_chim(nspec,lo(1):hi(1),lo(2):hi(2))
  real (rt) :: u_i_chim(lo(1):hi(1),lo(2):hi(2))
  real (rt) :: v_i_chim(lo(1):hi(1),lo(2):hi(2))
  real (rt) :: w_i_chim(lo(1):hi(1),lo(2):hi(2))
  real (rt) :: ye_i_chim(lo(1):hi(1),lo(2):hi(2))
  real (rt) :: a_aux_i_chim(lo(1):hi(1),lo(2):hi(2))
  real (rt) :: z_aux_i_chim(lo(1):hi(1),lo(2):hi(2))

  real (rt) :: rho_i_star(lo(1):hi(1),lo(2):hi(2))
  real (rt) :: t_i_star(lo(1):hi(1),lo(2):hi(2))
  real (rt) :: xn_i_star(nspec,lo(1):hi(1),lo(2):hi(2))
  real (rt) :: u_i_star(lo(1):hi(1),lo(2):hi(2))
  real (rt) :: ye_i_star(lo(1):hi(1),lo(2):hi(2))
  real (rt) :: a_aux_i_star(lo(1):hi(1),lo(2):hi(2))
  real (rt) :: z_aux_i_star(lo(1):hi(1),lo(2):hi(2))

  integer :: i, ii, j, jj, n

  real (rt) :: xg, yg
  real (rt) :: rg(nquad,nquad)
  real (rt) :: tg(nquad,nquad)
  real (rt) :: rho_quad(nquad,nquad)
  real (rt) :: t_quad(nquad,nquad)
  real (rt) :: p_quad(nquad,nquad)
  real (rt) :: u_quad(nquad,nquad)
  real (rt) :: v_quad(nquad,nquad)
  real (rt) :: w_quad(nquad,nquad)
  real (rt) :: xn_quad(nquad,nquad)
  real (rt) :: ye_quad(nquad,nquad)
  real (rt) :: a_aux_quad(nquad,nquad)
  real (rt) :: z_aux_quad(nquad,nquad)


  real (rt) :: xn_renorm(nspec)
  real (rt) :: ainv, ny, zy, zny, zzy, aa, zz, nn, alpha, beta

  type (eos_t) :: eos_state
 
  real (rt) :: drho
  real (rt), parameter :: dsmooth = 1.0e6_rt

  ! determine coordinates in r-theta and quadrature points
  do j = lo(2), hi(2)
    ycen(j) = xlo(2) + delta(2)*(dble(j-lo(2)) + HALF)
    do i = lo(1), hi(1)
      xcen(i) = xlo(1) + delta(1)*(dble(i-lo(1)) + HALF)
      r(i,j) = sqrt( xcen(i)**2 + ycen(j)**2 )
      if ( r(i,j) <= ZERO ) then
        theta(i,j) = ZERO
      else if ( ycen(j) == ZERO ) then
        theta(i,j) = HALF * m_pi
      else
        theta(i,j) = atan( xcen(i)/ycen(j) )
      end if
      if ( theta(i,j) < ZERO ) then
        theta(i,j) = theta(i,j) + m_pi
      end if

    end do
  end do

  if ( use_quad ) then
    do j = lo(2), hi(2)
      do i = lo(1), hi(1)
        do jj = 1, nquad
          yg = ycen(j) + HALF*delta(2)*xquad(jj)
          do ii = 1, nquad
            xg = xcen(i) + HALF*delta(1)*xquad(ii)
            rg(ii,jj) = sqrt( xg**2 + yg**2 )
            if ( rg(ii,jj) <= ZERO ) then
              tg(ii,jj) = ZERO
            else if ( yg == ZERO ) then
              tg(ii,jj) = HALF * m_pi
            else
              tg(ii,jj) = atan( xg/yg )
            end if
            if ( tg(ii,jj) < ZERO ) then
              tg(ii,jj) = tg(ii,jj) + m_pi
            end if
          end do
        end do

        select case (eos_input)
        case (eos_input_rt)
          call interp2dvol_chimera( rg, tg, rho_c_chim(:,:,1), rho_quad )
          call interp2drad_chimera( rg, tg, t_c_chim(:,:,1), t_quad )
        case (eos_input_rp)
          call interp2dvol_chimera( rg, tg, rho_c_chim(:,:,1), rho_quad )
          call interp2drad_chimera( rg, tg, p_c_chim(:,:,1), p_quad )
        case default
          call interp2dvol_chimera( rg, tg, rho_c_chim(:,:,1), rho_i_chim )
          call interp2drad_chimera( rg, tg, t_c_chim(:,:,1), t_i_chim )
        end select
        call interp2drad_chimera( rg, tg, u_c_chim(:,:,1), u_quad )
        call interp2drad_chimera( rg, tg, v_c_chim(:,:,1), v_quad )
        call interp2drad_chimera( rg, tg, w_c_chim(:,:,1), w_quad )
        do n = 1, nspec
          call interp2drad_chimera( rg, tg, xn_c_chim(n,:,:,1), xn_quad )
          xn_i_chim(n,i,j) = quad_avg( wquad, xn_quad )
        end do
        if ( naux == 1 ) then
          call interp2drad_chimera( rg, tg, ye_c_chim(:,:,1), ye_quad )
        else if ( naux == 2 ) then
          call interp2drad_chimera( rg, tg, a_aux_c_chim(:,:,1), a_aux_quad )
          call interp2drad_chimera( rg, tg, z_aux_c_chim(:,:,1), z_aux_quad )
        else if ( naux == 3 ) then
          call interp2drad_chimera( rg, tg, ye_c_chim(:,:,1), ye_quad )
          call interp2drad_chimera( rg, tg, a_aux_c_chim(:,:,1), a_aux_quad )
          call interp2drad_chimera( rg, tg, z_aux_c_chim(:,:,1), z_aux_quad )
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
        call interp2dvol_chimera( r, theta, rho_c_chim(:,:,1), rho_i_chim )
        call interp2drad_chimera( r, theta, t_c_chim(:,:,1), t_i_chim )
      case (eos_input_rp)
        call interp2dvol_chimera( r, theta, rho_c_chim(:,:,1), rho_i_chim )
        call interp2drad_chimera( r, theta, p_c_chim(:,:,1), p_i_chim )
      case default
        call interp2dvol_chimera( r, theta, rho_c_chim(:,:,1), rho_i_chim )
        call interp2drad_chimera( r, theta, t_c_chim(:,:,1), t_i_chim )
    end select
    do n = 1, nspec
      call interp2drad_chimera( r, theta, xn_c_chim(n,:,:,1), xn_i_chim(n,:,:) )
    end do
    call interp2drad_chimera( r, theta, u_c_chim(:,:,1), u_i_chim )
    call interp2drad_chimera( r, theta, v_c_chim(:,:,1), v_i_chim )
    call interp2drad_chimera( r, theta, w_c_chim(:,:,1), w_i_chim )
    if ( naux == 1 ) then
      call interp2drad_chimera( r, theta, ye_c_chim(:,:,1), ye_i_chim )
    else if ( naux == 2 ) then
      call interp2drad_chimera( r, theta, a_aux_c_chim(:,:,1), a_aux_i_chim )
      call interp2drad_chimera( r, theta, z_aux_c_chim(:,:,1), z_aux_i_chim )
    else if ( naux == 3 ) then
      call interp2drad_chimera( r, theta, ye_c_chim(:,:,1), ye_i_chim )
      call interp2drad_chimera( r, theta, a_aux_c_chim(:,:,1), a_aux_i_chim )
      call interp2drad_chimera( r, theta, z_aux_c_chim(:,:,1), z_aux_i_chim )
    end if
  end if

  if ( len_trim(star_fname) > 0 ) then
    call interp2dvol_star( r, rho_c_star, rho_i_star )
    call interp2drad_star( r, t_c_star, t_i_star )
    do n = 1, nspec
      call interp2drad_star( r, xn_c_star(:,n), xn_i_star(n,:,:) )
    end do
    call interp2drad_star( r, u_c_star, u_i_star )
    if ( naux == 1 ) then
      call interp2drad_star( r, ye_c_star, ye_i_star )
    else if ( naux == 2 ) then
      call interp2drad_star( r, a_aux_c_star, a_aux_i_star )
      call interp2drad_star( r, z_aux_c_star, z_aux_i_star )
    else if ( naux == 3 ) then
      call interp2drad_star( r, ye_c_star, ye_i_star )
      call interp2drad_star( r, a_aux_c_star, a_aux_i_star )
      call interp2drad_star( r, z_aux_c_star, z_aux_i_star )
    end if
  end if


  do j = lo(2), hi(2)
    do i = lo(1), hi(1)

    !  xn_i_chim(:,i,j) = xn_i_chim(:,i,j) / sum( xn_i_chim(:,i,j) )
    !  xn_i_star(:,i,j) = xn_i_star(:,i,j) / sum( xn_i_star(:,i,j) )

      if ( r(i,j) <= max_radius .or. len_trim(star_fname) == 0 ) then

        if ( r(i,j) < radius_inner ) then
!         drho = HALF * ( rho_i_chim(i,j) - rho_inner ) * ( one + tanh( (radius_inner - r(i,j))/dsmooth) )
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


        ! renormalize chimera mass fractions

        ny = 0.e0_rt
        zy = 0.e0_rt
        zny = 0.e0_rt
        zzy = 0.e0_rt

        do n = 1, nspec

          xn_renorm(n) = max(1.e-99_rt,min(1.e0_rt,xn_i_chim(n,i,j)))
          aa = aion(n)
          zz = zion(n)
          nn = nion(n)
          ainv = 1.e0_rt / aa
          ny = ny + xn_renorm(n)*nn*ainv
          zy = zy + xn_renorm(n)*zz*ainv
          zny = zny + xn_renorm(n)*zz*nn*ainv*ainv
          zzy = zzy + xn_renorm(n)*zz*zz*ainv*ainv
        enddo


        beta = (ye_i_chim(i,j)*ny-zny) / (ny*zzy - zy*zny)
        alpha = (1.e0_rt - beta*zy) / ny

       do n = 1,nspec

        aa = aion(n)
        zz = zion(n)
        nn = nion(n)

        xn_renorm(n) = xn_renorm(n)*(alpha*nn+beta*zz) / aa
        xn_renorm(n) = max(1.e-99_rt,min(1.e0_rt,xn_renorm(n)))

       enddo








        eos_state%xn(:) = xn_renorm(:)
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


        ! renormalize star mass fractions


        ny = 0.e0_rt
        zy = 0.e0_rt
        zny = 0.e0_rt
        zzy = 0.e0_rt

        do n = 1, nspec

          xn_renorm(n) = max(1.e-99_rt,min(1.e0_rt,xn_i_star(n,i,j)))
          aa = aion(n)
          zz = zion(n)
          nn = nion(n)
          ainv = 1.e0_rt / aa
          ny = ny + xn_renorm(n)*nn*ainv
          zy = zy + xn_renorm(n)*zz*ainv
          zny = zny + xn_renorm(n)*zz*nn*ainv*ainv
          zzy = zzy + xn_renorm(n)*zz*zz*ainv*ainv
        enddo


        beta = (ye_i_star(i,j)*ny-zny) / (ny*zzy - zy*zny)
        alpha = (1.e0_rt - beta*zy) / ny

       do n = 1,nspec

        aa = aion(n)
        zz = zion(n)
        nn = nion(n)

        xn_renorm(n) = xn_renorm(n)*(alpha*nn+beta*zz) / aa
        xn_renorm(n) = max(1.e-99_rt,min(1.e0_rt,xn_renorm(n)))

       enddo




        eos_state%rho = rho_i_star(i,j)
        eos_state%T = t_i_star(i,j)
        eos_state%xn(:) = xn_renorm(:)
        if ( naux == 1 ) then
          eos_state%aux = min( maxye, max( minye, ye_i_star(i,j) ) )
        else if ( naux == 2 ) then
          eos_state%aux(1) = a_aux_i_star(i,j)
          eos_state%aux(2) = z_aux_i_star(i,j)
        else if ( naux == 3 ) then
          eos_state%aux(1) = min( maxye, max( minye, ye_i_star(i,j) ) )
          eos_state%aux(2) = a_aux_i_star(i,j)
          eos_state%aux(3) = z_aux_i_star(i,j)
        end if

        call eos(eos_input_rt, eos_state)

        state(i,j,UMX) = u_i_star(i,j) * sin( theta(i,j) )
        state(i,j,UMY) = u_i_star(i,j) * cos( theta(i,j) )
        state(i,j,UMZ) = ZERO
        state(i,j,UFS:UFS+nspec-1) = xn_i_star(:,i,j)
        if ( naux == 1 ) then
          state(i,j,UFX) = min( maxye, max( minye, ye_i_star(i,j) ) )
        else if ( naux == 2 ) then
          state(i,j,UFX)   = a_aux_i_star(i,j)
          state(i,j,UFX+1) = z_aux_i_star(i,j)
        else if ( naux == 3 ) then
          state(i,j,UFX) = min( maxye, max( minye, ye_i_star(i,j) ) )
          state(i,j,UFX+1) = a_aux_i_star(i,j)
          state(i,j,UFX+2) = z_aux_i_star(i,j)
        end if

        state(i,j,URHO) = rho_i_star(i,j)
        state(i,j,UTEMP) = t_i_star(i,j)
        state(i,j,UEINT) = eos_state%e

      end if

    end do
  end do

  do j = lo(2), hi(2)
    do i = lo(1), hi(1)
      state(i,j,UEINT)   = state(i,j,URHO) * state(i,j,UEINT)
      state(i,j,UEDEN)   = state(i,j,UEINT) + state(i,j,URHO)*sum( HALF*state(i,j,UMX:UMZ)**2 )
      state(i,j,UMX:UMZ) = state(i,j,URHO) * state(i,j,UMX:UMZ)
      state(i,j,UFS:UFS+nspec-1) = state(i,j,URHO) * state(i,j,UFS:UFS+nspec-1)
      if ( naux > 0 ) then
        state(i,j,UFX:UFX+naux-1) = state(i,j,URHO) * state(i,j,UFX:UFX+naux-1)
      end if
    end do
  end do

  return
end subroutine ca_initdata





