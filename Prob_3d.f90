subroutine PROBINIT (init,name,namlen,problo,probhi)

  use probdata_module
  use chimera_parser_module
  use mesa_parser_module
  use bl_constants_module
  use bl_error_module
  use amrex_fort_module, only: rt => amrex_real
  use eos_module
  use parallel, only: parallel_IOProcessor
  use prob_params_module, only: center
  use interpolate_module, only: locate
  use model_interp_module, only: interp1d_linear
  use quadrature_module, only: quad_init

  implicit none
  integer :: init, namlen
  integer :: name(namlen)
  real (rt) :: problo(3), probhi(3)

  integer :: untin,i,j,k,dir
  real (rt) :: probhi_r, r, dr, rlo, rhi, x, y, z, theta, phi
  real (rt) :: dvol, volr, dvolr, gr
  real (rt) :: mass_chim, mass_mesa, vol_chim, vol_mesa
  real (rt) :: domega, domega_exclude, point_mass, mass_inner

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

  center(:) = zero

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
  &    eos_input /= eos_input_rp .and. &
  &    eos_input /= eos_input_re .and. &
  &    eos_input /= eos_input_ps ) then
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

  ! set the inner radius to be the max of the user-defined value or the inner grid boundary
  radius_inner = max( problo(1), radius_inner )

  ! set the interior density/temperature if not user-defined
  if ( radius_inner > zero ) then
!   volr = third * radius_inner**3
!   if ( volr <= volx_c_chim(1) ) then
!     i = 0
!   else if ( volr >= volx_c_chim(imax_chim) ) then
!     i = imax_chim
!   else
!     i = locate( volr, imax_chim, volx_c_chim ) - 1
!   end if
!   i_inner = i
!
!   if ( rho_inner <= zero ) then
!     call interp1d_linear( i, imax_chim, volx_c_chim, rhobar_c_chim, volr, rho_inner )
!   end if
!   if ( temp_inner <= zero ) then
!     call interp1d_linear( i, imax_chim, volx_c_chim, tbar_c_chim, volr, temp_inner )
!   end if
!   mass_inner = four * m_pi * rho_inner * third * radius_inner**3
    mass_inner = zero
  else
    mass_inner = zero
  end if

  ! largest radius of castro grid in spherical coordinates
  probhi_r = sqrt(sum(max(abs(problo),abs(probhi))**2))

  ! integrate mass from chimera data on castro grid
  mass_chim = zero
  vol_chim = zero
  do k = kmin_chim, kmax_chim
    phi = z_c_chim(k)
    do j = jmin_chim, jmax_chim
      theta = y_c_chim(j)
      do i = imin_chim, imax_chim
        r = x_c_chim(i)
        rlo = x_e_chim(i)
        rhi = x_e_chim(i+1)

        ! only use chimera data if zone is inside max_radius
        if ( rlo < max_radius ) then
          x = r * sin( theta ) * cos( phi )
          y = r * sin( theta ) * sin( phi )
          z = r * cos( theta )

          ! only include in sum if the zone is on castro grid
          if ( x > problo(1) .and. x < probhi(1) .and. &
               y > problo(2) .and. y < probhi(2) .and. &
               z > problo(3) .and. z < probhi(3) ) then

            if ( rhi > max_radius ) then
              dr = rhi - max_radius
              dvolr = dr * ( rlo * max_radius + dr * dr * third )
              dvol  = domega_chim(j,k) * ( dvolx_e_chim(i) - dvolr )

              vol_chim = vol_chim + dvol
              mass_chim = mass_chim + rho_c_chim(i,j,k) * dvol
            else
              mass_chim = mass_chim + dmass_e_chim(i,j,k)
              vol_chim = vol_chim + dvol_e_chim(i,j,k)
            end if
          end if
        end if

      end do
    end do
  end do

  ! integrate mass from mesa data on castro grid
  mass_mesa = zero
  vol_mesa = zero
  do i = 1, imax_mesa
    rlo = x_e_mesa(i)
    rhi = x_e_mesa(i+1)

    ! this includes mesa data on castro grid but not covered by chimera data
    if ( rhi > max_radius .and. rlo < probhi_r ) then

      r = half * ( rlo + rhi )

      dvolr = zero
      if ( rhi > probhi_r ) then
        dr = rhi - probhi_r
        dvolr = dvolr + dr * ( probhi_r * rhi + dr * dr * third )
      end if

      if ( rlo < max_radius ) then
        dr = max_radius - rlo
        dvolr = dvolr + dr * ( rlo * max_radius + dr * dr * third )
      end if

      ! TODO: calculate the solid angle of shell exterior to castro grid
      domega_exclude = zero
      domega = four * m_pi - domega_exclude

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
    write(*,'(a,es23.15)')  'total volume (castro)         =',product(probhi-problo)
  end if

  ! set up quadrature weights and abscissae
  if ( use_quad ) call quad_init( nquad )

  point_mass = zero
  do i = 1, imax_chim+1
    if ( x_e_chim(i) < radius_inner ) then
      if ( x_e_chim(i+1) <= radius_inner ) then
        point_mass = point_mass + sum( dmass_e_chim(i,jmin_chim:jmax_chim,kmin_chim:kmax_chim) )
      else
        dr = x_e_chim(i+1) - radius_inner
        dvolr = dr * ( radius_inner * x_e_chim(i+1) + dr * dr * third )
        point_mass = point_mass &
          &         + sum( rho_c_chim(i,jmin_chim:jmax_chim,kmin_chim:kmax_chim) &
          &                * domega_chim(jmin_chim:jmax_chim,kmin_chim:kmax_chim) ) &
          &         * ( dvolx_e_chim(i) - dvolr )
      end if
    end if
  end do
  point_mass = point_mass - mass_inner

  if (parallel_IOProcessor()) then
    write(*,*) 'point_mass=',point_mass
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

subroutine ca_initdata(level,time,lo,hi,nvar, &
                       state,state_l1,state_l2,state_l3,state_h1,state_h2,state_h3, &
                       dx,xlo,xhi)

  use bl_constants_module
  use bl_error_module
  use amrex_fort_module, only: rt => amrex_real
  use fundamental_constants_module
  use eos_module
  use eos_type_module, only: minye, maxye
  use meth_params_module, only: URHO, UMX, UMY, UMZ, UEINT, UFS, UTEMP, UEDEN, UFX, UFA
  use prob_params_module, only: center
  use network, only: nspec, naux
  use parallel, only: parallel_IOProcessor
  use probdata_module
  use chimera_parser_module
  use mesa_parser_module
  use quadrature_module, only: xquad, wquad, quad_avg

  implicit none

  integer :: level, nvar
  integer :: lo(3), hi(3)
  integer :: state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
  real (rt) :: xlo(3), xhi(3), time, dx(3)
  real (rt) :: state(state_l1:state_h1, &
                            state_l2:state_h2, &
                            state_l3:state_h3,nvar)

  ! local variables
  real (rt) :: xcen(lo(1):hi(1))
  real (rt) :: ycen(lo(2):hi(2))
  real (rt) :: zcen(lo(3):hi(3))

  real (rt) :: r(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
  real (rt) :: theta(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
  real (rt) :: phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))

  real (rt) :: rho_i_chim(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
  real (rt) :: t_i_chim(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
  real (rt) :: p_i_chim(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
  real (rt) :: e_i_chim(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
  real (rt) :: s_i_chim(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
  real (rt) :: xn_i_chim(nspec,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
  real (rt) :: u_i_chim(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
  real (rt) :: v_i_chim(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
  real (rt) :: w_i_chim(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
  real (rt) :: ye_i_chim(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
  real (rt) :: a_aux_i_chim(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
  real (rt) :: z_aux_i_chim(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))

  real (rt) :: rho_i_mesa(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
  real (rt) :: t_i_mesa(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
  real (rt) :: xn_i_mesa(nspec,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
  real (rt) :: u_i_mesa(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
  real (rt) :: ye_i_mesa(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
  real (rt) :: a_aux_i_mesa(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
  real (rt) :: z_aux_i_mesa(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))

  real (rt) :: xg, yg, zg
  real (rt) :: rg(nquad,nquad,nquad)
  real (rt) :: tg(nquad,nquad,nquad)
  real (rt) :: pg(nquad,nquad,nquad)
  real (rt) :: rho_quad(nquad,nquad,nquad)
  real (rt) :: t_quad(nquad,nquad,nquad)
  real (rt) :: p_quad(nquad,nquad,nquad)
  real (rt) :: e_quad(nquad,nquad,nquad)
  real (rt) :: s_quad(nquad,nquad,nquad)
  real (rt) :: u_quad(nquad,nquad,nquad)
  real (rt) :: v_quad(nquad,nquad,nquad)
  real (rt) :: w_quad(nquad,nquad,nquad)
  real (rt) :: xn_quad(nquad,nquad,nquad)
  real (rt) :: ye_quad(nquad,nquad,nquad)
  real (rt) :: a_aux_quad(nquad,nquad,nquad)
  real (rt) :: z_aux_quad(nquad,nquad,nquad)

  integer :: i, ii, j, jj, k, kk, n

  type (eos_t) :: eos_state

! real (rt), parameter :: delta = 2.9296875d6

  ! determine coordinates in r-theta-phi
  do k = lo(3), hi(3)
    zcen(k) = xlo(3) + dx(3)*(dble(k-lo(3)) + half) - center(3)
    do j = lo(2), hi(2)
      ycen(j) = xlo(2) + dx(2)*(dble(j-lo(2)) + half) - center(2)
      do i = lo(1), hi(1)
        xcen(i) = xlo(1) + dx(1)*(dble(i-lo(1)) + half) - center(1)
        r(i,j,k) = sqrt( xcen(i)**2 + ycen(j)**2 + zcen(k)**2 )
        if ( r(i,j,k) <= zero ) then
          theta(i,j,k) = zero
          phi(i,j,k) = zero
        else
          theta(i,j,k) = acos( zcen(k) / r(i,j,k) )
          phi(i,j,k) = mod( atan2( ycen(j), xcen(i) ), two*m_pi )
        end if
      end do
    end do
  end do


  ! interpolate state variables to r-theta-phi coordinates
  if ( use_quad ) then

    do k = lo(3), hi(3)
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)

          ! calculate quadrature points
          do kk = 1, nquad
            zg = zcen(k) + half*dx(3)*xquad(kk)
            do jj = 1, nquad
              yg = ycen(j) + half*dx(2)*xquad(jj)
              do ii = 1, nquad
                xg = xcen(i) + half*dx(1)*xquad(ii)
                rg(ii,jj,kk) = sqrt( xg**2 + yg**2 + zg**2 )
                if ( rg(ii,jj,kk) <= zero ) then
                  tg(ii,jj,kk) = zero
                  pg(ii,jj,kk) = zero
                else
                  tg(ii,jj,kk) = acos( zg / rg(ii,jj,kk) )
                  pg(ii,jj,kk) = mod( atan2( yg, xg ), two*m_pi )
                end if
              end do
            end do
          end do

          call interp3d_chimera( rg, tg, pg, u_c_chim(:,:,:), u_quad )
          call interp3d_chimera( rg, tg, pg, v_c_chim(:,:,:), v_quad )
          call interp3d_chimera( rg, tg, pg, w_c_chim(:,:,:), w_quad )
          call interp3d_chimera( rg, tg, pg, rho_c_chim(:,:,:), rho_quad )
          call interp3d_chimera( rg, tg, pg, t_c_chim(:,:,:), t_quad )
          call interp3d_chimera( rg, tg, pg, p_c_chim(:,:,:), p_quad )
          if ( trim(eos_name) == "stellarcollapse" ) then
            call interp3d_chimera( rg, tg, pg, ei_c_chim(:,:,:), e_quad )
          else
            call interp3d_chimera( rg, tg, pg, et_c_chim(:,:,:), e_quad )
          end if
          call interp3d_chimera( rg, tg, pg, s_c_chim(:,:,:), s_quad )
          do n = 1, nspec
            call interp3d_chimera( rg, tg, pg, xn_c_chim(n,:,:,:), xn_quad )
            xn_i_chim(n,i,j,k) = quad_avg( wquad, xn_quad )
          end do
          call interp3d_chimera( rg, tg, pg, ye_c_chim(:,:,:), ye_quad )
          call interp3d_chimera( rg, tg, pg, a_aux_c_chim(:,:,:), a_aux_quad )
          call interp3d_chimera( rg, tg, pg, z_aux_c_chim(:,:,:), z_aux_quad )

          u_i_chim(i,j,k) = quad_avg( wquad, u_quad )
          v_i_chim(i,j,k) = quad_avg( wquad, v_quad )
          w_i_chim(i,j,k) = quad_avg( wquad, w_quad )
          rho_i_chim(i,j,k) = quad_avg( wquad, rho_quad )
          t_i_chim(i,j,k) = quad_avg( wquad, t_quad )
          p_i_chim(i,j,k) = quad_avg( wquad, p_quad )
          e_i_chim(i,j,k) = quad_avg( wquad, e_quad )
          s_i_chim(i,j,k) = quad_avg( wquad, s_quad )
          ye_i_chim(i,j,k) = quad_avg( wquad, ye_quad )
          a_aux_i_chim(i,j,k) = quad_avg( wquad, a_aux_quad )
          z_aux_i_chim(i,j,k) = quad_avg( wquad, z_aux_quad )

        end do
      end do
    end do

  else

    call interp3d_chimera( r, theta, phi, u_c_chim(:,:,:), u_i_chim )
    call interp3d_chimera( r, theta, phi, v_c_chim(:,:,:), v_i_chim )
    call interp3d_chimera( r, theta, phi, w_c_chim(:,:,:), w_i_chim )
    call interp3d_chimera( r, theta, phi, rho_c_chim(:,:,:), rho_i_chim )
    call interp3d_chimera( r, theta, phi, t_c_chim(:,:,:), t_i_chim )
    call interp3d_chimera( r, theta, phi, p_c_chim(:,:,:), p_i_chim )
    if ( trim(eos_name) == "stellarcollapse" ) then
      call interp3d_chimera( r, theta, phi, ei_c_chim(:,:,:), e_i_chim )
    else
      call interp3d_chimera( r, theta, phi, et_c_chim(:,:,:), e_i_chim )
    end if
    call interp3d_chimera( r, theta, phi, s_c_chim(:,:,:), s_i_chim )
    do n = 1, nspec
      call interp3d_chimera( r, theta, phi, xn_c_chim(n,:,:,:), xn_i_chim(n,:,:,:) )
    end do
    call interp3d_chimera( r, theta, phi, ye_c_chim(:,:,:), ye_i_chim )
    call interp3d_chimera( r, theta, phi, a_aux_c_chim(:,:,:), a_aux_i_chim )
    call interp3d_chimera( r, theta, phi, z_aux_c_chim(:,:,:), z_aux_i_chim )

  end if

  if ( len_trim(mesa_fname) > 0 ) then
    call interp3d_mesa( r, u_c_mesa, u_i_mesa )
    call interp3d_mesa( r, rho_c_mesa, rho_i_mesa )
    call interp3d_mesa( r, t_c_mesa, t_i_mesa )
    do n = 1, nspec
      call interp3d_mesa( r, xn_c_mesa(:,n), xn_i_mesa(n,:,:,:) )
    end do
    call interp3d_mesa( r, ye_c_mesa, ye_i_mesa )
    call interp3d_mesa( r, a_aux_c_mesa, a_aux_i_mesa )
    call interp3d_mesa( r, z_aux_c_mesa, z_aux_i_mesa )
  end if

  do k = lo(3), hi(3)
    do j = lo(2), hi(2)
      do i = lo(1), hi(1)

        xn_i_chim(:,i,j,k) = xn_i_chim(:,i,j,k) / sum( xn_i_chim(:,i,j,k) )
        xn_i_mesa(:,i,j,k) = xn_i_mesa(:,i,j,k) / sum( xn_i_mesa(:,i,j,k) )

        if ( r(i,j,k) <= max_radius .or. len_trim(mesa_fname) == 0 ) then

          state(i,j,k,UMX) = + u_i_chim(i,j,k) * sin( theta(i,j,k) ) * cos( phi(i,j,k) ) &
                             + v_i_chim(i,j,k) * cos( theta(i,j,k) ) * cos( phi(i,j,k) ) &
                             - w_i_chim(i,j,k) * sin( theta(i,j,k) ) * sin( phi(i,j,k) )

          state(i,j,k,UMY) = + u_i_chim(i,j,k) * sin( theta(i,j,k) ) * sin( phi(i,j,k) ) &
                             + v_i_chim(i,j,k) * cos( theta(i,j,k) ) * sin( phi(i,j,k) ) &
                             + w_i_chim(i,j,k) * sin( theta(i,j,k) ) * cos( phi(i,j,k) )

          state(i,j,k,UMZ) = + u_i_chim(i,j,k) * cos( theta(i,j,k) ) &
                             - v_i_chim(i,j,k) * sin( theta(i,j,k) )

          eos_state%rho = rho_i_chim(i,j,k)
          eos_state%T = t_i_chim(i,j,k)
          eos_state%p = p_i_chim(i,j,k)
          eos_state%e = e_i_chim(i,j,k)
          eos_state%s = s_i_chim(i,j,k)
          eos_state%xn(:) = xn_i_chim(:,i,j,k)
          if ( naux == 1 ) then
            eos_state%aux(1) = min( maxye, max( minye, ye_i_chim(i,j,k) ) )
          else if ( naux == 2 ) then
            eos_state%aux(1) = a_aux_i_chim(i,j,k)
            eos_state%aux(2) = z_aux_i_chim(i,j,k)
          else if ( naux == 3 ) then
            eos_state%aux(1) = min( maxye, max( minye, ye_i_chim(i,j,k) ) )
            eos_state%aux(2) = a_aux_i_chim(i,j,k)
            eos_state%aux(3) = z_aux_i_chim(i,j,k)
          end if

          call eos(eos_input, eos_state)

        else

          state(i,j,k,UMX) = u_i_mesa(i,j,k) * sin( theta(i,j,k) ) * cos( phi(i,j,k) )
          state(i,j,k,UMY) = u_i_mesa(i,j,k) * sin( theta(i,j,k) ) * sin( phi(i,j,k) )
          state(i,j,k,UMZ) = u_i_mesa(i,j,k) * cos( theta(i,j,k) )

          eos_state%rho = rho_i_mesa(i,j,k)
          eos_state%T = t_i_mesa(i,j,k)
          eos_state%xn(:) = xn_i_mesa(:,i,j,k)
          if ( naux == 1 ) then
            eos_state%aux(1) = min( maxye, max( minye, ye_i_mesa(i,j,k) ) )
          else if ( naux == 2 ) then
            eos_state%aux(1) = a_aux_i_mesa(i,j,k)
            eos_state%aux(2) = z_aux_i_mesa(i,j,k)
          else if ( naux == 3 ) then
            eos_state%aux(1) = min( maxye, max( minye, ye_i_mesa(i,j,k) ) )
            eos_state%aux(2) = a_aux_i_mesa(i,j,k)
            eos_state%aux(3) = z_aux_i_mesa(i,j,k)
          end if

          call eos(eos_input_rt, eos_state)

        end if

        state(i,j,k,URHO)    = eos_state%rho
        state(i,j,k,UTEMP)   = eos_state%T
        state(i,j,k,UEINT)   = eos_state%rho * eos_state%e
        state(i,j,k,UEDEN)   = eos_state%rho * ( half * sum( state(i,j,k,UMX:UMZ)**2 ) + eos_state%e )
        state(i,j,k,UMX:UMZ) = eos_state%rho * state(i,j,k,UMX:UMZ)
        state(i,j,k,UFS:UFS+nspec-1) = eos_state%rho * eos_state%xn(:)
        if ( naux > 0 ) then
          state(i,j,k,UFX:UFX+naux-1) = eos_state%rho * eos_state%aux(:)
        end if

      end do
    end do
  end do

  return
end subroutine ca_initdata

