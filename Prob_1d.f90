subroutine PROBINIT (init,name,namlen,problo,probhi)

  use probdata_module
  use chimera_parser_module
  use mesa_parser_module
  use bl_constants_module
  use bl_error_module
  use bl_types
  use eos_module
  use prob_params_module, only : center
  use interpolate_module, only: locate
  use model_interp_module, only: interp1d_linear

  implicit none
  integer :: init, namlen
  integer :: name(namlen)
  real (dp_t) :: problo(1), probhi(1)

  integer :: untin,i,j,k,dir
  real (dp_t) :: r, dr, dvol, volr, dvolr, gr
  real (dp_t) :: mass_chim, mass_mesa, vol_chim, vol_mesa

  namelist /fortin/ model_name
  namelist /fortin/ mesa_name
  namelist /fortin/ model_eos_input
  namelist /fortin/ model_interp_method
  namelist /fortin/ r_inner
  namelist /fortin/ max_radius
  namelist /fortin/ do_particles

  !
  !     Build "probin" filename -- the name of file containing fortin namelist.
  !     
  integer, parameter :: maxlen = 127
  character(maxlen) :: probin
  character(maxlen) :: model
  integer :: ipp, ierr, ipp1

  if (namlen .gt. maxlen) call bl_error("probin file name too long")

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! set namelist defaults
  model_eos_input = eos_input_rt
  model_interp_method = 1
  r_inner = zero
  max_radius = zero
  model_name = ""
  mesa_name = ""
  do_particles = .false.

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

  ! integrate mass from CHIMERA model on CASTRO grid
  mass_chim = zero
  vol_chim = zero
  do i = imin_in, imax_in
    if ( rad_edge_in(i) < min(probhi(1),max_radius) ) then
      if ( rad_edge_in(i+1) <= min(probhi(1),max_radius) ) then
        mass_chim = mass_chim + sum(zone_mass_in(i,jmin_in:jmax_in,kmin_in:kmax_in))
        vol_chim = vol_chim + sum(volume_in(i,jmin_in:jmax_in,kmin_in:kmax_in))
      else
        dr = rad_edge_in(i+1) - min(probhi(1),max_radius)
        dvolr = dr * ( rad_edge_in(i) * min(probhi(1),max_radius) + dr * dr * third )
        dvol  = four * m_pi * ( dvol_rad_in(i) - dvolr )

        vol_chim = vol_chim + dvol
        mass_chim = mass_chim + &
        &           sum( dens_in(i,jmin_in:jmax_in,kmin_in:kmax_in) * domega_in(jmin_in:jmax_in,kmin_in:kmax_in) ) &
        &           * ( dvol_rad_in(i) - dvolr )
      end if
    end if
  end do

  ! integrate mass from MESA model on CASTRO grid
  mass_mesa = zero
  vol_mesa = zero
  do i = 1, nx_mesa_in
    if ( rad_edge_mesa_in(i+1) > max_radius .and. &
    &    rad_edge_mesa_in(i)   < probhi(1) ) then

      dvolr = zero
      if ( rad_edge_mesa_in(i+1) > probhi(1) ) then
        dr = rad_edge_mesa_in(i+1) - probhi(1)
        dvolr = dvolr + dr * ( probhi(1) * rad_edge_mesa_in(i) + dr * dr * third )
      end if

      if ( rad_edge_mesa_in(i) < max_radius ) then
        dr = max_radius - rad_edge_mesa_in(i) 
        dvolr = dvolr + dr * ( rad_edge_mesa_in(i) * max_radius + dr * dr * third )
      end if

      dvol = four * m_pi * ( dvol_rad_mesa_in(i) - dvolr )

      vol_mesa = vol_mesa + dvol
      mass_mesa = mass_mesa + dvol * dens_mesa_in(i)

    end if
  end do

  ! Get the gravitational acceleration to compare with CASTRO
  r = min( probhi(1), max_radius )
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

  point_mass = zero
  do i = 1, imax_in
    if ( rad_cntr_in(i) <= problo(1) ) then
      point_mass = point_mass + sum( zone_mass_in(i,:,:) )
    else if ( vol_rad_cntr_in(i) <= third*r_inner**3 ) then
      point_mass = point_mass + sum( zone_mass_in(i,:,:) )
    end if
  end do
  if (parallel_IOProcessor()) then
    write(*,*) 'point_mass=',point_mass
  end if

  center(1) = zero

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
! ::: nvar      => number of state components.
! ::: state     <= scalar array
! ::: dx        => cell size
! ::: xlo, xhi  => physical locations of lower left and upper
! :::              right hand corner of grid.  (does not include
! :::		   ghost region).
! ::: -----------------------------------------------------------

subroutine ca_initdata(level,time,lo,hi,nvar, &
                       state,state_l1,state_h1, &
                       dx,xlo,xhi)

  use bl_error_module
  use bl_types
  use chimera_parser_module
  use mesa_parser_module
  use fundamental_constants_module
  use eos_module
  use meth_params_module, only : URHO, UMX, UMY, UMZ, UEINT, UFS, UTEMP, UEDEN, UFX, UFA
  use network, only: nspec
  use probdata_module

  implicit none

  integer :: level, nvar
  integer :: lo(1), hi(1)
  integer :: state_l1,state_h1
  double precision :: xlo(1), xhi(1), time, dx(1)
  double precision :: state(state_l1:state_h1,nvar)

  ! local variables
  real (dp_t) :: x(lo(1):hi(1))
  real (dp_t) :: rho_chim(lo(1):hi(1))
  real (dp_t) :: temp_chim(lo(1):hi(1))
  real (dp_t) :: pressure_chim(lo(1):hi(1))
  real (dp_t) :: eint_chim(lo(1):hi(1))
  real (dp_t) :: entropy_chim(lo(1):hi(1))
  real (dp_t) :: xn_chim(lo(1):hi(1),nspec)
  real (dp_t) :: vrad_chim(lo(1):hi(1))
  real (dp_t) :: vtheta_chim(lo(1):hi(1))
  real (dp_t) :: vphi_chim(lo(1):hi(1))
  real (dp_t) :: rho_mesa(lo(1):hi(1))
  real (dp_t) :: temp_mesa(lo(1):hi(1))
  real (dp_t) :: xn_mesa(lo(1):hi(1),nspec)
  real (dp_t) :: vrad_mesa(lo(1):hi(1))
  integer :: i, n

  type (eos_t) :: eos_state

  x = zero
  do i = lo(1), hi(1)
    x(i) = xlo(1) + dx(1)*(dble(i-lo(1)) + half)
  end do

  select case (model_eos_input)
    case (eos_input_rt)
      call interp1d_chimera( x, dens_in(:,1,1), rho_chim, model_interp_method )
      call interp1d_chimera( x, temp_in(:,1,1), temp_chim, model_interp_method )
    case (eos_input_rp)
      call interp1d_chimera( x, dens_in(:,1,1), rho_chim, model_interp_method )
      call interp1d_chimera( x, pres_in(:,1,1), pressure_chim, model_interp_method )
    case (eos_input_re)
      call interp1d_chimera( x, dens_in(:,1,1), rho_chim, model_interp_method )
      call interp1d_chimera( x, eint_in(:,1,1), eint_chim, model_interp_method )
    case (eos_input_ps)
      call interp1d_chimera( x, enpy_in(:,1,1), entropy_chim, model_interp_method )
      call interp1d_chimera( x, pres_in(:,1,1), pressure_chim, model_interp_method )
    case default
      call interp1d_chimera( x, dens_in(:,1,1), rho_chim, model_interp_method )
      call interp1d_chimera( x, temp_in(:,1,1), temp_chim, model_interp_method )
  end select
  do n = 1, nspec
    call interp1d_chimera( x, xn_in(n,:,1,1), xn_chim(:,n), model_interp_method )
  end do
  call interp1d_chimera( x, vrad_in(:,1,1), vrad_chim, model_interp_method )
  call interp1d_chimera( x, vtheta_in(:,1,1), vtheta_chim, model_interp_method )
  call interp1d_chimera( x, vphi_in(:,1,1), vphi_chim, model_interp_method )

  if ( .not. trim(mesa_name) == "" ) then
    call interp1d_mesa( x, dens_mesa_in, rho_mesa, model_interp_method )
    call interp1d_mesa( x, temp_mesa_in, temp_mesa, model_interp_method )
    do n = 1, nspec
      call interp1d_mesa( x, xn_mesa_in(:,n), xn_mesa(:,n), model_interp_method )
    end do
    call interp1d_mesa( x, vrad_mesa_in, vrad_mesa, model_interp_method )
  end if

  do i = lo(1), hi(1)

    xn_chim(i,:) = xn_chim(i,:) / sum( xn_chim(i,:) )
    xn_mesa(i,:) = xn_mesa(i,:) / sum( xn_mesa(i,:) )

    if ( x(i) <= max_radius .or. trim(mesa_name) == "" ) then
      select case (model_eos_input)
        case (eos_input_rt)
          eos_state%rho = rho_chim(i)
          eos_state%T = temp_chim(i)
        case (eos_input_rp)
          eos_state%rho = rho_chim(i)
          eos_state%p = pressure_chim(i)
        case (eos_input_re)
          eos_state%rho = rho_chim(i)
          eos_state%e = eint_chim(i)
        case (eos_input_ps)
          eos_state%p = pressure_chim(i)
          eos_state%s = entropy_chim(i)
        case default
          eos_state%rho = rho_chim(i)
          eos_state%T = temp_chim(i)
      end select
      eos_state%xn(:) = xn_chim(i,:)

      call eos(model_eos_input, eos_state)

      state(i,UMX) = vrad_chim(i)
      state(i,UMY) = vtheta_chim(i)
      state(i,UMZ) = vphi_chim(i)
      state(i,UFS:UFS+nspec-1) = xn_chim(i,:)

      select case (model_eos_input)
        case (eos_input_rt)
          state(i,URHO) = rho_chim(i)
          state(i,UTEMP) = temp_chim(i)
          state(i,UEINT) = eos_state%e
        case (eos_input_rp)
          state(i,URHO) = rho_chim(i)
          state(i,UTEMP) = eos_state%T
          state(i,UEINT) = eos_state%e
        case (eos_input_re)
          state(i,URHO) = rho_chim(i)
          state(i,UTEMP) = eos_state%T
          state(i,UEINT) = eint_chim(i)
        case (eos_input_ps)
          state(i,URHO) = eos_state%rho
          state(i,UTEMP) = eos_state%T
          state(i,UEINT) = eos_state%e
        case default
          state(i,URHO) = rho_chim(i)
          state(i,UTEMP) = temp_chim(i)
          state(i,UEINT) = eos_state%e
      end select

    else

      eos_state%rho = rho_mesa(i)
      eos_state%T = temp_mesa(i)
      eos_state%xn(:) = xn_mesa(i,:)
      call eos(eos_input_rt, eos_state)

      state(i,UMX) = vrad_mesa(i)
      state(i,UMY) = zero
      state(i,UMZ) = zero
      state(i,UFS:UFS+nspec-1) = xn_mesa(i,:)

      state(i,URHO) = rho_mesa(i)
      state(i,UTEMP) = temp_mesa(i)
      state(i,UEINT) = eos_state%e

!     write(*,*) 'using mesa values'

    end if

  end do

  do i = lo(1), hi(1)

    state(i,UEINT)   = state(i,URHO) * state(i,UEINT)
    state(i,UEDEN)   = state(i,UEINT) + state(i,URHO)*sum( half*state(i,UMX:UMZ)**2 )
    state(i,UMX:UMZ) = state(i,URHO) * state(i,UMX:UMZ)
    state(i,UFS:UFS+nspec-1) = state(i,URHO) * state(i,UFS:UFS+nspec-1)

  end do

  return
end subroutine ca_initdata

