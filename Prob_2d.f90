subroutine PROBINIT (init,name,namlen,problo,probhi)

  use probdata_module
  use chimera_parser_module
  use mesa_parser_module
  use bl_constants_module
  use bl_error_module
  use eos_module
  use parallel, only: parallel_IOProcessor
  use meth_params_module, only : point_mass
  use prob_params_module, only : center

  implicit none
  integer init, namlen
  integer name(namlen)
  double precision problo(2), probhi(2)

  integer untin,i,j,k,dir
  double precision tmp_r, tmp_z

  namelist /fortin/ model_name
  namelist /fortin/ mesa_name
  namelist /fortin/ model_eos_input
  namelist /fortin/ model_interp_method
  namelist /fortin/ min_radius
  namelist /fortin/ max_radius
  namelist /fortin/ do_particles

  !
  !     Build "probin" filename -- the name of file containing fortin namelist.
  !     
  integer, parameter :: maxlen = 127
  character probin*(maxlen)
  character model*(maxlen)
  integer ipp, ierr, ipp1, n

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
    max_radius = rad_cntr_in(imax_in)
  end if

  point_mass = zero
  do i = imin_in, imax_in
    do j = jmin_in, jmax_in
      tmp_r = rad_cntr_in(i) * sin( theta_cntr_in(j) )
      tmp_z = rad_cntr_in(i) * cos( theta_cntr_in(j) )
      if ( tmp_r <= problo(1) .and. tmp_z <= problo(2) ) then
        point_mass = point_mass + sum( zone_mass_in(i,j,:) )
      else if ( vol_rad_cntr_in(i) <= third*min_radius**3 ) then
        point_mass = point_mass + sum( zone_mass_in(i,j,:) )
      end if
    end do
  end do
  if (parallel_IOProcessor()) then
    write(*,*) 'point_mass=',point_mass
  end if

  center(1) = zero
  center(2) = zero

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
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEINT, UFS, UTEMP, UEDEN, UFX, UFA
  use network, only: nspec
  use probdata_module

  implicit none

  integer :: level, nscal
  integer :: lo(2), hi(2)
  integer :: state_l1,state_l2,state_h1,state_h2
  double precision :: xlo(2), xhi(2), time, delta(2)
  double precision :: state(state_l1:state_h1,state_l2:state_h2,NVAR)

  ! local variables
  real (dp_t) :: x(lo(1):hi(1))
  real (dp_t) :: y(lo(2):hi(2))
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
  real (dp_t) :: rho_mesa(lo(1):hi(1),lo(2):hi(2))
  real (dp_t) :: temp_mesa(lo(1):hi(1),lo(2):hi(2))
  real (dp_t) :: xn_mesa(lo(1):hi(1),lo(2):hi(2),nspec)
  real (dp_t) :: vrad_mesa(lo(1):hi(1),lo(2):hi(2))
  integer :: i, j, n

  type (eos_t) :: eos_state

  do i = lo(1), hi(1)
    x(i) = xlo(1) + delta(1)*(dble(i-lo(1)) + half)
  end do
  do j = lo(2), hi(2)
    y(j) = xlo(2) + delta(2)*(dble(j-lo(2)) + half)
  end do

  do j = lo(2), hi(2)
    do i = lo(1), hi(1)
      r(i,j) = sqrt( x(i)**2 + y(j)**2 )
      if ( r(i,j) <= zero ) then
        theta(i,j) = zero
      else if ( y(j) == zero ) then
        theta(i,j) = half * m_pi
      else
        theta(i,j) = atan( x(i)/y(j) )
      end if
      if ( theta(i,j) < zero ) then
        theta(i,j) = theta(i,j) + m_pi
      end if
    end do
  end do

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


  if ( .not. trim(mesa_name) == "" ) then
    call interp2d_mesa( r, dens_mesa_in, rho_mesa, model_interp_method )
    call interp2d_mesa( r, temp_mesa_in, temp_mesa, model_interp_method )
    do n = 1, nspec
      call interp2d_mesa( r, xn_mesa_in(:,n), xn_mesa(:,:,n), model_interp_method )
    end do
    call interp2d_mesa( r, vrad_mesa_in, vrad_mesa, model_interp_method )
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

        call eos(model_eos_input, eos_state)

        state(i,j,UMX) = vrad_chim(i,j) * sin( theta(i,j) ) + vtheta_chim(i,j) * cos( theta(i,j) )
        state(i,j,UMY) = vrad_chim(i,j) * cos( theta(i,j) ) - vtheta_chim(i,j) * sin( theta(i,j) )
        state(i,j,UMZ) = vphi_chim(i,j)
        state(i,j,UFS:UFS+nspec-1) = xn_chim(i,j,:)

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
        call eos(eos_input_rt, eos_state)

        state(i,j,UMX) = vrad_mesa(i,j) * sin( theta(i,j) )
        state(i,j,UMY) = vrad_mesa(i,j) * cos( theta(i,j) )
        state(i,j,UMZ) = zero
        state(i,j,UFS:UFS+nspec-1) = xn_mesa(i,j,:)

        state(i,j,URHO) = rho_mesa(i,j)
        state(i,j,UTEMP) = temp_mesa(i,j)
        state(i,j,UEINT) = eos_state%e

!       write(*,*) 'using mesa values'
        
      end if

    end do
  end do

  do j = lo(2), hi(2)
    do i = lo(1), hi(1)
        
      state(i,j,UEINT)   = state(i,j,URHO) * state(i,j,UEINT)
      state(i,j,UEDEN)   = state(i,j,UEINT) + state(i,j,URHO)*sum( half*state(i,j,UMX:UMZ)**2 )
      state(i,j,UMX:UMZ) = state(i,j,URHO) * state(i,j,UMX:UMZ)
      state(i,j,UFS:UFS+nspec-1) = state(i,j,URHO) * state(i,j,UFS:UFS+nspec-1)

    end do
  end do

  return
end subroutine ca_initdata

