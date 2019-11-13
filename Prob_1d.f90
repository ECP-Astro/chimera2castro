subroutine PROBINIT (init,name,namlen,problo,probhi)

  use probdata_module
  use chimera_parser_module
  use star_parser_module
  use amrex_constants_module
  use amrex_error_module
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
  real (rt) :: problo(1), probhi(1)

  integer :: untin,i,dir
  real (rt) :: r, dr, dvol, volr, dvolr, gr
  real (rt) :: mass_chim, mass_star, vol_chim, vol_star
  real (rt) :: domega, point_mass, mass_inner

  namelist /fortin/ chimera_fname
  namelist /fortin/ star_fname
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

  if (namlen .gt. maxlen) call bl_error("probin file name too long")

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! set namelist defaults
  chimera_fname = ""
  star_fname = ""
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

  ! load star model
  if ( len_trim(star_fname) > 0 ) then
    call read_star_file(star_fname)
  end if

  ! set the max radius to use chimera data
  ! if r > max_radius, use star data
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

  ! integrate mass from chimera data on castro grid
  mass_chim = zero
  vol_chim = zero
  do i = imin_chim, imax_chim
    if ( x_e_chim(i) < min(probhi(1),max_radius) ) then
      if ( x_e_chim(i+1) <= min(probhi(1),max_radius) ) then
        mass_chim = mass_chim + sum(dmass_e_chim(i,jmin_chim:jmax_chim,kmin_chim:kmax_chim))
        vol_chim = vol_chim + sum(dvol_e_chim(i,jmin_chim:jmax_chim,kmin_chim:kmax_chim))
      else
        dr = x_e_chim(i+1) - min(probhi(1),max_radius)
        dvolr = dr * ( x_e_chim(i) * min(probhi(1),max_radius) + dr * dr * third )
        dvol  = four * m_pi * ( dvolx_e_chim(i) - dvolr )

        vol_chim = vol_chim + dvol
        mass_chim = mass_chim &
          &         + sum( rho_c_chim(i,jmin_chim:jmax_chim,kmin_chim:kmax_chim) &
          &                * domega_chim(jmin_chim:jmax_chim,kmin_chim:kmax_chim) ) &
          &         * ( dvolx_e_chim(i) - dvolr )
      end if
    end if
  end do

  ! integrate mass from star data on castro grid
  mass_star = zero
  vol_star = zero
  do i = 1, imax_star
    if ( x_e_star(i+1) > max_radius .and. &
    &    x_e_star(i)   < probhi(1) ) then

      domega = four * m_pi

      dvolr = zero
      if ( x_e_star(i+1) > probhi(1) ) then
        dr = x_e_star(i+1) - probhi(1)
        dvolr = dvolr + dr * ( probhi(1) * x_e_star(i) + dr * dr * third )
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
  r = min( probhi(1), max_radius )
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

    write(*,'(a,2es23.15)') 'total mass   (chimera, star)  =',mass_chim,mass_star
    write(*,'(a,es23.15)')  'total mass   (chimera + star) =',mass_chim+mass_star

    write(*,'(a,2es23.15)') 'total volume (chimera, star)  =',vol_chim,vol_star
    write(*,'(a,es23.15)')  'total volume (chimera + star) =',vol_chim+vol_star
    dr = probhi(1) - problo(1)
    write(*,'(a,es23.15)')  'total volume (castro)         =',four * m_pi * dr * ( problo(1) * probhi(1) + dr * dr * third )
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
                       state,state_l1,state_h1, &
                       dx,xlo,xhi)

  use amrex_constants_module
  use amrex_error_module
  use amrex_fort_module, only: rt => amrex_real
  use fundamental_constants_module
  use eos_module
  use eos_type_module, only: minye, maxye
  use meth_params_module, only : URHO, UMX, UMY, UMZ, UEINT, UFS, UTEMP, UEDEN, UFX, UFA
  use network, only: nspec, naux
  use parallel, only: parallel_IOProcessor
  use probdata_module
  use chimera_parser_module
  use star_parser_module
  use quadrature_module, only: xquad, wquad, quad_avg

  implicit none

  integer :: level, nvar
  integer :: lo(1), hi(1)
  integer :: state_l1,state_h1
  double precision :: xlo(1), xhi(1), time, dx(1)
  double precision :: state(state_l1:state_h1,nvar)

  ! local variables
  real (rt) :: xcen(lo(1):hi(1))

  real (rt) :: rho_i_chim(lo(1):hi(1))
  real (rt) :: t_i_chim(lo(1):hi(1))
  real (rt) :: p_i_chim(lo(1):hi(1))
  real (rt) :: e_i_chim(lo(1):hi(1))
  real (rt) :: s_i_chim(lo(1):hi(1))
  real (rt) :: xn_i_chim(nspec,lo(1):hi(1))
  real (rt) :: u_i_chim(lo(1):hi(1))
  real (rt) :: v_i_chim(lo(1):hi(1))
  real (rt) :: w_i_chim(lo(1):hi(1))
  real (rt) :: ye_i_chim(lo(1):hi(1))
  real (rt) :: a_aux_i_chim(lo(1):hi(1))
  real (rt) :: z_aux_i_chim(lo(1):hi(1))

  real (rt) :: rho_i_star(lo(1):hi(1))
  real (rt) :: t_i_star(lo(1):hi(1))
  real (rt) :: xn_i_star(nspec,lo(1):hi(1))
  real (rt) :: u_i_star(lo(1):hi(1))
  real (rt) :: ye_i_star(lo(1):hi(1))
  real (rt) :: a_aux_i_star(lo(1):hi(1))
  real (rt) :: z_aux_i_star(lo(1):hi(1))

  real (rt) :: xg(nquad)
  real (rt) :: rho_quad(nquad)
  real (rt) :: t_quad(nquad)
  real (rt) :: p_quad(nquad)
  real (rt) :: e_quad(nquad)
  real (rt) :: s_quad(nquad)
  real (rt) :: u_quad(nquad)
  real (rt) :: v_quad(nquad)
  real (rt) :: w_quad(nquad)
  real (rt) :: xn_quad(nquad)
  real (rt) :: ye_quad(nquad)
  real (rt) :: a_aux_quad(nquad)
  real (rt) :: z_aux_quad(nquad)

  integer :: i, ii, n

  type (eos_t) :: eos_state

! real (rt), parameter :: delta = 2.9296875d6

  do i = lo(1), hi(1)
    xcen(i) = xlo(1) + dx(1)*(dble(i-lo(1)) + half)
  end do

  if ( use_quad ) then
    do i = lo(1), hi(1)
      do ii = 1, nquad
        xg(ii) = xcen(i) + half*dx(1)*xquad(ii)
      end do

      call interp1drad_chimera( xg, u_c_chim(:,1,1), u_quad )
      call interp1drad_chimera( xg, v_c_chim(:,1,1), v_quad )
      call interp1drad_chimera( xg, w_c_chim(:,1,1), w_quad )
      call interp1dvol_chimera( xg, rho_c_chim(:,1,1), rho_quad )
      call interp1drad_chimera( xg, t_c_chim(:,1,1), t_quad )
      call interp1drad_chimera( xg, p_c_chim(:,1,1), p_quad )
      if ( trim(eos_name) == "stellarcollapse" ) then
        call interp1drad_chimera( xg, ei_c_chim(:,1,1), e_quad )
      else
        call interp1drad_chimera( xg, et_c_chim(:,1,1), e_quad )
      end if
      call interp1drad_chimera( xg, s_c_chim(:,1,1), s_quad )
      do n = 1, nspec
        call interp1drad_chimera( xg, xn_c_chim(n,:,1,1), xn_quad )
        xn_i_chim(n,i) = quad_avg( wquad, xn_quad )
      end do
      call interp1drad_chimera( xg, ye_c_chim(:,1,1), ye_quad )
      call interp1drad_chimera( xg, a_aux_c_chim(:,1,1), a_aux_quad )
      call interp1drad_chimera( xg, z_aux_c_chim(:,1,1), z_aux_quad )

      u_i_chim(i) = quad_avg( wquad, u_quad )
      v_i_chim(i) = quad_avg( wquad, v_quad )
      w_i_chim(i) = quad_avg( wquad, w_quad )
      rho_i_chim(i) = quad_avg( wquad, rho_quad )
      t_i_chim(i) = quad_avg( wquad, t_quad )
      p_i_chim(i) = quad_avg( wquad, p_quad )
      e_i_chim(i) = quad_avg( wquad, e_quad )
      s_i_chim(i) = quad_avg( wquad, s_quad )
      ye_i_chim(i) = quad_avg( wquad, ye_quad )
      a_aux_i_chim(i) = quad_avg( wquad, a_aux_quad )
      z_aux_i_chim(i) = quad_avg( wquad, z_aux_quad )
    end do
  else
    call interp1drad_chimera( xcen, u_c_chim(:,1,1), u_i_chim )
    call interp1drad_chimera( xcen, v_c_chim(:,1,1), v_i_chim )
    call interp1drad_chimera( xcen, w_c_chim(:,1,1), w_i_chim )
    call interp1dvol_chimera( xcen, rho_c_chim(:,1,1), rho_i_chim )
    call interp1drad_chimera( xcen, t_c_chim(:,1,1), t_i_chim )
    call interp1drad_chimera( xcen, p_c_chim(:,1,1), p_i_chim )
    if ( trim(eos_name) == "stellarcollapse" ) then
      call interp1drad_chimera( xcen, ei_c_chim(:,1,1), e_i_chim )
    else
      call interp1drad_chimera( xcen, et_c_chim(:,1,1), e_i_chim )
    end if
    call interp1drad_chimera( xcen, s_c_chim(:,1,1), s_i_chim )
    do n = 1, nspec
      call interp1drad_chimera( xcen, xn_c_chim(n,:,1,1), xn_i_chim(n,:) )
    end do
    call interp1drad_chimera( xcen, ye_c_chim(:,1,1), ye_i_chim )
    call interp1drad_chimera( xcen, a_aux_c_chim(:,1,1), a_aux_i_chim )
    call interp1drad_chimera( xcen, z_aux_c_chim(:,1,1), z_aux_i_chim )
  end if

  if ( len_trim(star_fname) > 0 ) then
    call interp1drad_star( xcen, u_c_star, u_i_star )
    call interp1dvol_star( xcen, rho_c_star, rho_i_star )
    call interp1drad_star( xcen, t_c_star, t_i_star )
    do n = 1, nspec
      call interp1drad_star( xcen, xn_c_star(:,n), xn_i_star(n,:) )
    end do
    call interp1drad_star( xcen, ye_c_star, ye_i_star )
    call interp1drad_star( xcen, a_aux_c_star, a_aux_i_star )
    call interp1drad_star( xcen, z_aux_c_star, z_aux_i_star )
  end if

  do i = lo(1), hi(1)

    xn_i_chim(:,i) = xn_i_chim(:,i) / sum( xn_i_chim(:,i) )
    xn_i_star(:,i) = xn_i_star(:,i) / sum( xn_i_star(:,i) )

    if ( xcen(i) <= max_radius .or. len_trim(star_fname) == 0 ) then

      state(i,UMX) = u_i_chim(i)
      state(i,UMY) = v_i_chim(i)
      state(i,UMZ) = w_i_chim(i)

      eos_state%rho = rho_i_chim(i)
      eos_state%T = t_i_chim(i)
      eos_state%p = p_i_chim(i)
      eos_state%e = e_i_chim(i)
      eos_state%s = s_i_chim(i)
      eos_state%xn(:) = xn_i_chim(:,i)
      if ( naux == 1 ) then
        eos_state%aux(1) = min( maxye, max( minye, ye_i_chim(i) ) )
      else if ( naux == 2 ) then
        eos_state%aux(1) = a_aux_i_chim(i)
        eos_state%aux(2) = z_aux_i_chim(i)
      else if ( naux == 3 ) then
        eos_state%aux(1) = min( maxye, max( minye, ye_i_chim(i) ) )
        eos_state%aux(2) = a_aux_i_chim(i)
        eos_state%aux(3) = z_aux_i_chim(i)
      end if

      call eos(eos_input, eos_state)

    else

      state(i,UMX) = u_i_star(i)
      state(i,UMY) = zero
      state(i,UMZ) = zero

      eos_state%rho = rho_i_star(i)
      eos_state%T = t_i_star(i)
      eos_state%xn(:) = xn_i_star(:,i)
      if ( naux == 1 ) then
        eos_state%aux(1) = min( maxye, max( minye, ye_i_star(i) ) )
      else if ( naux == 2 ) then
        eos_state%aux(1) = a_aux_i_star(i)
        eos_state%aux(2) = z_aux_i_star(i)
      else if ( naux == 3 ) then
        eos_state%aux(1) = min( maxye, max( minye, ye_i_star(i) ) )
        eos_state%aux(2) = a_aux_i_star(i)
        eos_state%aux(3) = z_aux_i_star(i)
      end if

      call eos(eos_input_rt, eos_state)

    end if

    state(i,URHO)    = eos_state%rho
    state(i,UTEMP)   = eos_state%T
    state(i,UEINT)   = eos_state%rho * eos_state%e
    state(i,UEDEN)   = eos_state%rho * ( half * sum( state(i,UMX:UMZ)**2 ) + eos_state%e )
    state(i,UMX:UMZ) = eos_state%rho * state(i,UMX:UMZ)
    state(i,UFS:UFS+nspec-1) = eos_state%rho * eos_state%xn(:)
    if ( naux > 0 ) then
      state(i,UFX:UFX+naux-1) = eos_state%rho * eos_state%aux(:)
    end if

  end do

  return
end subroutine ca_initdata

