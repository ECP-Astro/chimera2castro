module chimera_parser_module

  ! we read in the number of variables and their order and use this to
  ! map them into the model_state array.  We ignore anything other than
  ! density, temperature, pressure and composition.
  !
  ! composition is assumed to be in terms of mass fractions     

  use parallel, only: parallel_IOProcessor
  use network
  use amrex_fort_module, only: rt => amrex_real
  use hdf5_read_write
  use hdf5

  implicit none

  ! array dimensions
  integer, save :: nx_chim, ny_chim, nz_chim, nnc_chim, npart_chim, npart_shell_chim

  ! radial index bounds
  integer, save :: imin_chim, imax_chim

  ! theta index bounds
  integer, save :: jmin_chim, jmax_chim

  ! phi index bounds
  integer, save :: kmin_chim, kmax_chim

  ! grid coordinates
  real (rt), allocatable, save :: x_e_chim(:), dx_e_chim(:)
  real (rt), allocatable, save :: y_e_chim(:), dy_e_chim(:)
  real (rt), allocatable, save :: z_e_chim(:), dz_e_chim(:)
  real (rt), allocatable, save :: x_c_chim(:), dx_c_chim(:)
  real (rt), allocatable, save :: y_c_chim(:), dy_c_chim(:)
  real (rt), allocatable, save :: z_c_chim(:), dz_c_chim(:)

  real (rt), allocatable, save :: volx_e_chim(:), dvolx_e_chim(:)
  real (rt), allocatable, save :: voly_e_chim(:), dvoly_e_chim(:)
  real (rt), allocatable, save :: volz_e_chim(:), dvolz_e_chim(:)
  real (rt), allocatable, save :: volx_c_chim(:), dvolx_c_chim(:)
  real (rt), allocatable, save :: voly_c_chim(:), dvoly_c_chim(:)
  real (rt), allocatable, save :: volz_c_chim(:), dvolz_c_chim(:)

  real (rt), allocatable, save :: domega_chim(:,:)
  real (rt), allocatable, save :: dvol_e_chim(:,:,:)
  real (rt), allocatable, save :: dmass_e_chim(:,:,:)
  real (rt), save :: omega_chim

  ! time data
  real (rt), save :: time_chim, t_bounce

  ! fluid velocities
  real (rt), allocatable, save :: u_c_chim(:,:,:)
  real (rt), allocatable, save :: v_c_chim(:,:,:)
  real (rt), allocatable, save :: w_c_chim(:,:,:)

  ! thermo data
  real (rt), allocatable, save :: rho_c_chim(:,:,:)
  real (rt), allocatable, save :: t_c_chim(:,:,:)
  real (rt), allocatable, save :: ye_c_chim(:,:,:)
  real (rt), allocatable, save :: p_c_chim(:,:,:)
  real (rt), allocatable, save :: ei_c_chim(:,:,:)
  real (rt), allocatable, save :: et_c_chim(:,:,:)
  real (rt), allocatable, save :: eb_c_chim(:,:,:)
  real (rt), allocatable, save :: em_c_chim(:,:,:)
  real (rt), allocatable, save :: s_c_chim(:,:,:)

  real (rt), allocatable, save :: rhobar_c_chim(:)
  real (rt), allocatable, save :: tbar_c_chim(:)

  ! gravity data
  real (rt), allocatable, save :: gravx_c_chim(:,:,:)
  real (rt), allocatable, save :: gravy_c_chim(:,:,:)
  real (rt), allocatable, save :: gravz_c_chim(:,:,:)
  real (rt), allocatable, save :: gravx_c_avg_chim(:)

  ! nuclei data
  real (rt), allocatable, save :: a_nuc_chim(:)
  real (rt), allocatable, save :: z_nuc_chim(:)
  real (rt), allocatable, save :: mex_nuc_chim(:)
  real (rt), allocatable, save :: be_nuc_chim(:)
  character (5), allocatable, save :: name_nuc_chim(:)

  ! composition data
  real (rt), allocatable, save :: xn_read(:,:,:,:)
  real (rt), allocatable, save :: xn_c_chim(:,:,:,:)
  integer,     allocatable, save :: nse_c_chim(:,:,:)

  ! representative heavy/auxiliary data
  real (rt), allocatable, save :: a_aux_c_chim(:,:,:)
  real (rt), allocatable, save :: z_aux_c_chim(:,:,:)
  real (rt), allocatable, save :: be_aux_c_chim(:,:,:)

  ! particle data
  real (rt), allocatable, save :: px_chim(:)
  real (rt), allocatable, save :: py_chim(:)
  real (rt), allocatable, save :: pz_chim(:)

  ! chimera_initialized will be .true. once the model is read in and the
  ! model data arrays are initialized and filled
  logical, save :: chimera_initialized = .false.

  ! HDF5 variables
  integer(hid_t) :: file_id
  integer(hid_t) :: group_id
  integer(hid_t) :: dataset_id
  integer(hid_t) :: dataspace_id
  integer(hid_t) :: plist_id

  integer(hsize_t) :: datasize1d(1)
  integer(hsize_t) :: datasize2d(2)
  integer(hsize_t) :: datasize3d(3)
  integer(hsize_t) :: datasize4d(4)
  integer(hsize_t) :: datasize5d(5)

  integer(hsize_t) :: slab_offset2d(2)
  integer(hsize_t) :: slab_offset3d(3)
  integer(hsize_t) :: slab_offset4d(4)
  integer(hsize_t) :: slab_offset5d(5)

  public :: open_chimera_file, read_chimera_file, close_chimera_file

contains

  subroutine open_chimera_file(fname)

    use bl_constants_module
    use bl_error_module

    character (*), intent(in) :: fname

    ! local variables
    integer :: i_read(1)

    integer :: ierr

    call h5open_f(ierr)
    call h5fopen_f(fname, H5F_ACC_RDONLY_F, file_id, ierr)

    if (ierr /= 0) then
      if (parallel_IOProcessor()) print *,'Could not open file: ',fname
      call bl_error('Aborting now -- please supply filename')
    end if

    call h5gopen_f(file_id, '/mesh', group_id, ierr)
    if (ierr /= 0) call bl_error('Could not open /mesh group')

    datasize1d(1) = 1
    call read_1d_slab('nz_hyperslabs', i_read, group_id, datasize1d)

    call h5gclose_f(group_id, ierr)

    return

  end subroutine open_chimera_file

  subroutine read_chimera_file

    use actual_network, only: nspec, aion, zion
    use bl_constants_module
    use bl_error_module
    use fundamental_constants_module
    use probdata_module

    ! local variables
    real (rt) :: d_read(1)

    integer :: array_dimensions(3)
    integer :: index_bounds(2)
    integer :: eos_dimensions(2)
    integer :: part_dimensions(2)

    integer :: i, j, k, ierr
    integer :: itmp(nspec)
    integer :: net_to_castro(nspec)
    integer, allocatable :: net_in_castro(:)

    real (rt) :: dx, dvol

    ! open mesh group
    call h5gopen_f(file_id, '/mesh', group_id, ierr)
    if (ierr /= 0) call bl_error('Could not open /mesh group')

    ! read array dimensions
    datasize1d(1) = 3
    call read_1d_slab('array_dimensions', array_dimensions, group_id, datasize1d)
    nx_chim = array_dimensions(1)
    ny_chim = array_dimensions(2)
    nz_chim = array_dimensions(3)

    ! read radial index bounds
    datasize1d(1) = 2
    call read_1d_slab('radial_index_bound', index_bounds, group_id, datasize1d)
    imin_chim = index_bounds(1)
    imax_chim = index_bounds(2)

    ! read theta index bounds
    call read_1d_slab('theta_index_bound', index_bounds, group_id, datasize1d)
    jmin_chim = index_bounds(1)
    jmax_chim = index_bounds(2)
    if ( jmin_chim /= 1 .or. jmax_chim /= ny_chim ) then
       if (parallel_IOProcessor()) print *,'Error: inner theta bounds /= (1,ny_chim), (',jmin_chim,',',jmax_chim,')'
       call bl_error('Aborting now -- check theta bounds in model file')
     end if

    ! read phi index bounds
    call read_1d_slab('phi_index_bound', index_bounds, group_id, datasize1d)
    kmin_chim = index_bounds(1)
    kmax_chim = index_bounds(2)
    if ( kmin_chim /= 1 .or. kmax_chim /= nz_chim ) then
       if (parallel_IOProcessor()) print *,'Error: inner phi bounds /= (1,nz_chim), (',kmin_chim,',',kmax_chim,')'
       call bl_error('Aborting now -- check phi bounds in model file')
     end if

    ! allocate grid variables
    allocate (x_e_chim(nx_chim+1))
    allocate (y_e_chim(-5:ny_chim+7))
    allocate (z_e_chim(-5:nz_chim+7))
    allocate (x_c_chim(nx_chim))
    allocate (y_c_chim(-5:ny_chim+6))
    allocate (z_c_chim(-5:nz_chim+6))
    allocate (dx_e_chim(nx_chim))
    allocate (dy_e_chim(-5:ny_chim+6))
    allocate (dz_e_chim(-5:nz_chim+6))
    allocate (dx_c_chim(nx_chim))
    allocate (dy_c_chim(-5:ny_chim+6))
    allocate (dz_c_chim(-5:nz_chim+6))

    ! read zone edge coordinaets
    datasize1d(1) = nx_chim+1
    call read_1d_slab('x_ef', x_e_chim, group_id, datasize1d)
    datasize1d(1) = ny_chim+1
    call read_1d_slab('y_ef', y_e_chim(1:ny_chim+1), group_id, datasize1d)
    datasize1d(1) = nz_chim+1
    call read_1d_slab('z_ef', z_e_chim(1:nz_chim+1), group_id, datasize1d)

    ! read zone center coordinates
    datasize1d(1) = nx_chim
    call read_1d_slab('x_cf', x_c_chim, group_id, datasize1d)
    datasize1d(1) = ny_chim
    call read_1d_slab('y_cf', y_c_chim(1:ny_chim), group_id, datasize1d)
    datasize1d(1) = nz_chim
    call read_1d_slab('z_cf', z_c_chim(1:nz_chim), group_id, datasize1d)

    ! read zone widths
    datasize1d(1) = nx_chim
    call read_1d_slab('dx_cf', dx_c_chim, group_id, datasize1d)
    datasize1d(1) = ny_chim
    call read_1d_slab('dy_cf', dy_c_chim(1:ny_chim), group_id, datasize1d)
    datasize1d(1) = nz_chim
    call read_1d_slab('dz_cf', dz_c_chim(1:nz_chim), group_id, datasize1d)

    ! fill grid boundaries
    call grid_bc(y_e_chim, y_c_chim, dy_c_chim, 1, ny_chim, 0, 0)
    call grid_bc(z_e_chim, z_c_chim, dz_c_chim, 1, nz_chim, 3, 3)
    dx_e_chim(1:nx_chim) = x_e_chim(2:nx_chim+1) - x_e_chim(1:nx_chim)
    dy_e_chim(-5:ny_chim+6) = y_e_chim(-4:ny_chim+7) - y_e_chim(-5:ny_chim+6)
    dz_e_chim(-5:nz_chim+6) = z_e_chim(-4:nz_chim+7) - z_e_chim(-5:nz_chim+6)

    ! read times
    datasize1d(1) = 1
    call read_1d_slab('time', d_read, group_id, datasize1d)
    time_chim = d_read(1)
    call read_1d_slab('t_bounce', d_read, group_id, datasize1d)
    t_bounce = d_read(1)

    ! solid angles and volumes
    allocate (volx_e_chim(nx_chim+1))
    allocate (voly_e_chim(ny_chim+1))
    allocate (volz_e_chim(nz_chim+1))
    allocate (volx_c_chim(nx_chim))
    allocate (voly_c_chim(ny_chim))
    allocate (volz_c_chim(nz_chim))
    allocate (dvolx_e_chim(nx_chim))
    allocate (dvoly_e_chim(ny_chim))
    allocate (dvolz_e_chim(nz_chim))
    allocate (dvolx_c_chim(nx_chim))
    allocate (dvoly_c_chim(ny_chim))
    allocate (dvolz_c_chim(nz_chim))
    allocate (domega_chim(ny_chim,nz_chim))
    allocate (dvol_e_chim(nx_chim,ny_chim,nz_chim))
    allocate (dmass_e_chim(nx_chim,ny_chim,nz_chim))

    datasize2d = (/ ny_chim, nz_chim /)
    call read_2d_slab('d_omega', domega_chim, group_id, datasize2d)

    volx_e_chim(1) = third * x_e_chim(1)**3
    do i = 1, nx_chim
      dx               = dx_e_chim(i)
      dvolx_e_chim(i)  = dx * ( x_e_chim(i) * x_e_chim(i+1) + dx * dx * third )
      volx_e_chim(i+1) = volx_e_chim(i) + dvolx_e_chim(i)
    end do
    volx_c_chim(1:nx_chim) = volx_e_chim(1:nx_chim) + half*dvolx_e_chim(1:nx_chim)
    dvolx_c_chim(1:nx_chim-1) = volx_c_chim(2:nx_chim) - volx_c_chim(1:nx_chim-1)
    dvolx_c_chim(nx_chim) = dvolx_c_chim(nx_chim-1)

    voly_e_chim(1) = one - cos( y_e_chim(1) )
    do j = 1, ny_chim
      dvoly_e_chim(j)  = cos( y_e_chim(j) ) - cos( y_e_chim(j+1) )
      voly_e_chim(j+1) = voly_e_chim(j) + dvoly_e_chim(j)
    end do
    voly_c_chim(1:ny_chim) = voly_e_chim(1:ny_chim) + half*dvoly_e_chim(1:ny_chim)
    dvoly_c_chim(1:ny_chim-1) = voly_c_chim(2:ny_chim) - voly_c_chim(1:ny_chim-1)
    dvoly_c_chim(ny_chim) = dvoly_c_chim(ny_chim-1)

    volz_e_chim(1) = z_e_chim(1)
    do k = 1, nz_chim
      dvolz_e_chim(k)  = dz_e_chim(k)
      volz_e_chim(k+1) = volz_e_chim(k) + dvolz_e_chim(k)
    end do
    volz_c_chim(1:nz_chim) = volz_e_chim(1:nz_chim) + half*dvolz_e_chim(1:nz_chim)
    dvolz_c_chim(1:nz_chim-1) = volz_c_chim(2:nz_chim) - volz_c_chim(1:nz_chim-1)
    dvolz_c_chim(nz_chim) = dvolz_c_chim(nz_chim-1)

    omega_chim = sum( domega_chim(:,:) )
    do i = 1, nx_chim
      dvol_e_chim(i,:,:) = dvolx_e_chim(i) * domega_chim(:,:)
    end do

    ! close mesh group
    call h5gclose_f(group_id, ierr)

    ! open fluid group
    call h5gopen_f(file_id, '/fluid', group_id, ierr)
    if (ierr /= 0) call bl_error('Could not open /fluid group')

    ! allocate fluid variables
    allocate (u_c_chim(nx_chim,ny_chim,nz_chim))
    allocate (v_c_chim(nx_chim,ny_chim,nz_chim))
    allocate (w_c_chim(nx_chim,ny_chim,nz_chim))
    allocate (rho_c_chim(nx_chim,ny_chim,nz_chim))
    allocate (t_c_chim(nx_chim,ny_chim,nz_chim))
    allocate (ye_c_chim(nx_chim,ny_chim,nz_chim))
    allocate (p_c_chim(nx_chim,ny_chim,nz_chim))
    allocate (ei_c_chim(nx_chim,ny_chim,nz_chim))
    allocate (et_c_chim(nx_chim,ny_chim,nz_chim))
    allocate (eb_c_chim(nx_chim,ny_chim,nz_chim))
    allocate (em_c_chim(nx_chim,ny_chim,nz_chim))
    allocate (s_c_chim(nx_chim,ny_chim,nz_chim))
    allocate (gravx_c_chim(nx_chim,ny_chim,nz_chim))
    allocate (gravy_c_chim(nx_chim,ny_chim,nz_chim))
    allocate (gravz_c_chim(nx_chim,ny_chim,nz_chim))
    allocate (gravx_c_avg_chim(nx_chim))
    allocate (rhobar_c_chim(nx_chim))
    allocate (tbar_c_chim(nx_chim))

    ! read velocity variables
    datasize3d = (/ nx_chim, ny_chim, nz_chim /)
    slab_offset3d = (/ 0, 0, 0 /)
    call read_ray_hyperslab('u_c',     u_c_chim, group_id, datasize3d, slab_offset3d)
    call read_ray_hyperslab('v_c',     v_c_chim, group_id, datasize3d, slab_offset3d)
    call read_ray_hyperslab('w_c',     w_c_chim, group_id, datasize3d, slab_offset3d)

    ! read thermo variables
    call read_ray_hyperslab('rho_c',   rho_c_chim, group_id, datasize3d, slab_offset3d)
    call read_ray_hyperslab('t_c',     t_c_chim,   group_id, datasize3d, slab_offset3d)
    call read_ray_hyperslab('ye_c',    ye_c_chim,  group_id, datasize3d, slab_offset3d)
    call read_ray_hyperslab('press',   p_c_chim,   group_id, datasize3d, slab_offset3d)
    call read_ray_hyperslab('e_int',   ei_c_chim,  group_id, datasize3d, slab_offset3d)
    call read_ray_hyperslab('entropy', s_c_chim,   group_id, datasize3d, slab_offset3d)

    ! convert units to erg/K
    s_c_chim(:,:,:) = s_c_chim(:,:,:) * k_B * n_A

    dmass_e_chim(:,:,:) = rho_c_chim(:,:,:) * dvol_e_chim(:,:,:)

    ! read gravity variables
    call read_ray_hyperslab('grav_x_c', gravx_c_chim, group_id, datasize3d, slab_offset3d)
    call read_ray_hyperslab('grav_y_c', gravy_c_chim, group_id, datasize3d, slab_offset3d)
    call read_ray_hyperslab('grav_z_c', gravz_c_chim, group_id, datasize3d, slab_offset3d)
    do i = imax_chim+1, nx_chim
      gravx_c_chim(i,:,:) = gravx_c_chim(imax_chim,:,:) * (volx_c_chim(imax_chim)/volx_c_chim(i))**two3rd
    end do

    do i = 1, nx_chim
      rhobar_c_chim(i) = sum( rho_c_chim(i,:,:) * domega_chim(:,:) ) / omega_chim
      tbar_c_chim(i) = sum( t_c_chim(i,:,:) * domega_chim(:,:) ) / omega_chim
      gravx_c_avg_chim(i) = sum( gravx_c_chim(i,:,:) * domega_chim(:,:) ) / omega_chim
    end do

    ! read number of nuclei
    datasize1d(1) = 2
    call read_1d_slab('dimeos', eos_dimensions, group_id, datasize1d)
    nnc_chim = eos_dimensions(1)

    ! close fluid group
    call h5gclose_f(group_id, ierr)

    ! open abundance group
    call h5gopen_f(file_id, '/abundance', group_id, ierr)
    if (ierr /= 0) call bl_error('Could not open /abundance group')

    ! allocate composition variables
    allocate (a_nuc_chim(nnc_chim-1))
    allocate (z_nuc_chim(nnc_chim-1))
    allocate (mex_nuc_chim(nnc_chim-1))
    allocate (be_nuc_chim(nnc_chim-1))
    allocate (name_nuc_chim(nnc_chim))
    allocate (xn_read(nnc_chim,nx_chim,ny_chim,nz_chim))
    allocate (xn_c_chim(nspec,nx_chim,ny_chim,nz_chim))
    allocate (nse_c_chim(nx_chim+1,ny_chim,nz_chim))
    allocate (a_aux_c_chim(nx_chim,ny_chim,nz_chim))
    allocate (z_aux_c_chim(nx_chim,ny_chim,nz_chim))
    allocate (be_aux_c_chim(nx_chim,ny_chim,nz_chim))

    ! read nuclei values
    datasize1d(1) = nnc_chim-1
    call read_1d_slab('a_nuc',    a_nuc_chim,   group_id, datasize1d)
    call read_1d_slab('z_nuc',    z_nuc_chim,   group_id, datasize1d)
    call read_1d_slab('m_ex_nuc', mex_nuc_chim, group_id, datasize1d)

    ! binding energy
    be_nuc_chim(1:nnc_chim-1) = z_nuc_chim(1:nnc_chim-1)*7.2889705_rt + &
    &                           (a_nuc_chim(1:nnc_chim-1)-z_nuc_chim(1:nnc_chim-1))*8.0713171_rt - &
    &                           mex_nuc_chim(1:nnc_chim-1)
    
    ! read nuclei names
    datasize1d(1) = nnc_chim
    call read_1d_slab('a_name', name_nuc_chim, 5, group_id, datasize1d)

    ! read mass fractions
    datasize4d = (/ nnc_chim, nx_chim, ny_chim, nz_chim /)
    slab_offset4d = (/ 0, 0, 0, 0 /)
    call read_ray_hyperslab('xn_c', xn_read, group_id, datasize4d, slab_offset4d)

    ! create lookup tables for isotopes in castro net
    k = 0
    net_to_castro(:) = nnc_chim
    do i = 1, nspec
      do j = 1, nnc_chim-1
        if ( nint(aion(i)) == nint(a_nuc_chim(j)) .and. nint(zion(i)) == nint(z_nuc_chim(j)) ) then
          k = k + 1
          net_to_castro(i) = j
          itmp(k) = i
          exit
        end if
      end do
      if ( j > nnc_chim-1 ) then
        if (parallel_IOProcessor()) then
          write(*,'(2(a,i3),a)') ' could not find isotope (',nint(zion(i)),',',nint(aion(i)),') in CHIMERA net'
        end if
      end if
    end do

    if ( k > 0 ) then
      allocate (net_in_castro(k))
      net_in_castro(:) = itmp(1:k)
    else
      call bl_error("no species in chimera net in castro net")
    end if

    xn_c_chim(net_in_castro,:,:,:) = xn_read(net_to_castro(net_in_castro),:,:,:)
    if ( nspec > nspec_evolve ) then
      xn_c_chim(nspec,:,:,:) = xn_read(nnc_chim,:,:,:)
    end if

    ! read nse state
    datasize3d = (/ nx_chim+1, ny_chim, nz_chim /)
    slab_offset3d = (/ 0, 0, 0 /)
    call read_ray_hyperslab('nse_c', nse_c_chim, group_id, datasize3d, slab_offset3d)

    ! read representative heavy nucleus/auxiliary nucleus properties
    datasize3d = (/ nx_chim, ny_chim, nz_chim /)
    slab_offset3d = (/ 0, 0, 0 /)
    call read_ray_hyperslab('a_nuc_rep_c',  a_aux_c_chim,  group_id, datasize3d, slab_offset3d)
    call read_ray_hyperslab('z_nuc_rep_c',  z_aux_c_chim,  group_id, datasize3d, slab_offset3d)
    call read_ray_hyperslab('be_nuc_rep_c', be_aux_c_chim, group_id, datasize3d, slab_offset3d)

    ! binding energy
    eb_c_chim(:,:,:) = zero
    do i = 1, nnc_chim-1
      eb_c_chim(:,:,:) = eb_c_chim(:,:,:) - xn_read(i,:,:,:)*be_nuc_chim(i)/a_nuc_chim(i)
    end do
    eb_c_chim(:,:,:) = eb_c_chim(:,:,:) - xn_read(nnc_chim,:,:,:)*be_aux_c_chim(:,:,:)/a_aux_c_chim(:,:,:)

    ! convert from MeV/c^2 to erg/g
    eb_c_chim(:,:,:) = eb_c_chim(:,:,:) * MeV2eV * ev2erg * n_A

    ! rest mass energy (including electrons) relative to free neutrons
    em_c_chim(:,:,:) = - ( m_n - m_p - m_e ) * ye_c_chim(:,:,:) * c_light * c_light * n_A + eb_c_chim(:,:,:)

    ! remove rest mass energy and constant offset from internal energy
    et_c_chim(:,:,:) = ei_c_chim(:,:,:) - em_c_chim(:,:,:) - 8.9_rt * MeV2eV * ev2erg * n_A

    ! close abundance group
    call h5gclose_f(group_id, ierr)

    if ( do_particles ) then
      ! open particle group
      call h5gopen_f(file_id, '/particle', group_id, ierr)
      if (ierr /= 0) then
        call bl_error('Could not open /particle group')
      else
        ! read particle dimensions
        datasize1d(1) = 2
        call read_1d_slab('partdim', part_dimensions, group_id, datasize1d)
        npart_chim = part_dimensions(1)
        npart_shell_chim = part_dimensions(2)

        ! allocate particle variables
        allocate (px_chim(npart_chim))
        allocate (py_chim(npart_chim))
        allocate (pz_chim(npart_chim))

        ! read particle variables
        datasize1d(1) = npart_chim
        call read_1d_slab('px', px_chim, group_id, datasize1d)
        call read_1d_slab('py', py_chim, group_id, datasize1d)
        call read_1d_slab('pz', pz_chim, group_id, datasize1d)

        ! close particle group
        call h5gclose_f(group_id, ierr)
      end if
    end if

    ! close hdf5 file
    call h5fclose_f(file_id, ierr)
    call h5close_f(ierr)

    chimera_initialized = .true.

    return

  end subroutine read_chimera_file


  subroutine close_chimera_file

    use probdata_module
    
    if (chimera_initialized) then
       deallocate(x_e_chim)
       deallocate(y_e_chim)
       deallocate(z_e_chim)
       deallocate(x_c_chim)
       deallocate(y_c_chim)
       deallocate(z_c_chim)
       deallocate(dx_e_chim)
       deallocate(dy_e_chim)
       deallocate(dz_e_chim)
       deallocate(dx_c_chim)
       deallocate(dy_c_chim)
       deallocate(dz_c_chim)
       deallocate(volx_e_chim)
       deallocate(voly_e_chim)
       deallocate(volz_e_chim)
       deallocate(volx_c_chim)
       deallocate(voly_c_chim)
       deallocate(volz_c_chim)
       deallocate(dvolx_e_chim)
       deallocate(dvoly_e_chim)
       deallocate(dvolz_e_chim)
       deallocate(dvolx_c_chim)
       deallocate(dvoly_c_chim)
       deallocate(dvolz_c_chim)
       deallocate(domega_chim)
       deallocate(dvol_e_chim)
       deallocate(dmass_e_chim)
       deallocate(u_c_chim)
       deallocate(v_c_chim)
       deallocate(w_c_chim)
       deallocate(rho_c_chim)
       deallocate(t_c_chim)
       deallocate(ye_c_chim)
       deallocate(p_c_chim)
       deallocate(ei_c_chim)
       deallocate(et_c_chim)
       deallocate(eb_c_chim)
       deallocate(em_c_chim)
       deallocate(s_c_chim)
       deallocate(gravx_c_chim)
       deallocate(gravy_c_chim)
       deallocate(gravz_c_chim)
       deallocate(gravx_c_avg_chim)
       deallocate(rhobar_c_chim)
       deallocate(tbar_c_chim)
       deallocate(a_nuc_chim)
       deallocate(z_nuc_chim)
       deallocate(mex_nuc_chim)
       deallocate(be_nuc_chim)
       deallocate(name_nuc_chim)
       deallocate(xn_read)
       deallocate(xn_c_chim)
       deallocate(nse_c_chim)
       deallocate(a_aux_c_chim)
       deallocate(z_aux_c_chim)
       deallocate(be_aux_c_chim)
       if ( do_particles ) then
         deallocate(px_chim)
         deallocate(py_chim)
         deallocate(pz_chim)
       end if
       chimera_initialized = .false.
    endif

    return

  end subroutine close_chimera_file

  subroutine grid_bc( xag, xcg, dxg, nmin, nmax, nleft, nright )

    use bl_constants_module
    use bl_error_module

    !  boundary condition flags : nleft, nright
    !    = 0 : reflecting
    !    = 3 : periodic (eg, u(nmin-1) = u(nmax))

    ! input variables
    integer, intent(in) :: nmin         ! minimum coordinate index
    integer, intent(in) :: nmax         ! maximum coordinate index
    integer, intent(in) :: nleft        ! left boundary type
    integer, intent(in) :: nright       ! right boundary type

    ! output variables
    real (rt), intent(inout) :: xag(nmin-6:nmax+7)    ! grid edges
    real (rt), intent(inout) :: xcg(nmin-6:nmax+6)    ! grid centers
    real (rt), intent(inout) :: dxg(nmin-6:nmax+6)    ! grid widths

    ! local variables
    integer :: n        ! zone index

    ! left (inner) boundary
    select case (nleft)

    case (0,162)        ! grid symmetric with left (inner) edge

      do n = 1, 6
        dxg(nmin - n)  = dxg(nmin + n - 1)
        xag(nmin - n)  = xag(nmin - n + 1) - dxg(nmin-n)
        xcg(nmin - n)  = xag(nmin - n) + half * dxg(nmin - n)
      end do

    case (3)            ! periodic

      do n = 1, 6
        dxg(nmin - n)  = dxg(nmax + 1 - n)
        xag(nmin - n)  = xag(nmin - n + 1) - dxg(nmin - n)
        xcg(nmin - n)  = xag(nmin - n) + half * dxg(nmin - n)
      end do

    case default

      call bl_error("invalid value for boundary condition, nleft")

    end select

    ! right (outer) boundary
    select case (nright)

    case (0,162)        ! grid symmetric with right (outer) edge

      do n = 1, 6
        dxg(nmax + n)     = dxg(nmax + 1 - n)
        xag(nmax + n + 1) = xag(nmax + n) + dxg(nmax + n)
        xcg(nmax + n)     = xag(nmax + n) + half * dxg(nmax + n)
      end do ! n = 1, 6

    case (3)            ! grid symetric with left (inner) edge

      do n = 1, 6
        dxg(nmax + n)     = dxg(nmin + n - 1)
        xag(nmax + n + 1) = xag(nmax + n) + dxg(nmax + n)
        xcg(nmax + n)     = xag(nmax + n) + half * dxg(nmax + n)
      end do ! n = 1, 6

    case default

      call bl_error("invalid value for boundary condition, nright")

    end select

    return
  end subroutine grid_bc

  subroutine interp1d_chimera( x_out, state_chim, state_out )

    use bl_constants_module
    use bl_error_module
    use interpolate_module, only: locate
    use model_interp_module, only: interp1d_linear, interp1d_spline
    use probdata_module, only: interp_method

    ! input variables
    real (rt), intent(in) :: x_out(:)
    real (rt), intent(in) :: state_chim(:)

    ! output variables
    real (rt), intent(out) :: state_out(size(x_out))

    ! local variables
    integer :: ix(size(x_out))
    integer :: ix_max

    real (rt) :: dvol1, dvol2
    real (rt) :: volx_out(size(x_out))

    integer :: i, n

    volx_out = third * x_out**3
    ix_max = imax_chim

    if ( interp_method == 1 ) then
      ix(:) = 0
      do i = 1, size(x_out)
        if ( volx_out(i) <= volx_c_chim(1) ) then
          ix(i) = 0
        else if ( volx_out(i) >= volx_c_chim(ix_max) ) then
          ix(i) = ix_max
        else
          ix(i) = locate( volx_out(i), ix_max, volx_c_chim ) - 1
        end if
      end do
      call interp1d_linear( ix, ix_max, volx_c_chim, state_chim, volx_out, state_out )
    else if ( interp_method == 2 ) then
      call interp1d_spline( volx_c_chim, state_chim, volx_out, state_out )
    else
      call bl_error("invalid value for interp_method")
    end if

    return
  end subroutine interp1d_chimera

  subroutine interp2d_chimera( x_out, y_out, state_chim, state_out )

    use bl_constants_module
    use bl_error_module
    use interpolate_module, only: locate
    use model_interp_module, only: interp2d_linear, interp2d_spline
    use probdata_module, only: interp_method

    ! input variables
    real (rt), intent(in) :: x_out(:,:)
    real (rt), intent(in) :: y_out(:,:)
    real (rt), intent(in) :: state_chim(:,:)

    ! output variables
    real (rt), intent(out) :: state_out(size(x_out,1),size(x_out,2))

    ! local variables
    integer :: ix(size(x_out,1),size(x_out,2))
    integer :: ix_max

    integer :: iy(size(x_out,1),size(x_out,2))
    integer :: iy_max

    real (rt) :: volx_out(size(x_out,1),size(x_out,2))
    real (rt) :: voly_out(size(x_out,1),size(x_out,2))

    real (rt) :: dvol1, dvol2

    integer :: i, j, n

    volx_out   = third * x_out**3
    voly_out = one - cos( y_out )
    ix_max = imax_chim
    iy_max = jmax_chim

    if ( interp_method == 1 ) then

      ix(:,:) = 0
      iy(:,:) = 0
      do j = 1, size(x_out,2)
        do i = 1, size(x_out,1)
          if ( volx_out(i,j) <= volx_c_chim(1) ) then
            ix(i,j) = 0
          else if ( volx_out(i,j) >= volx_c_chim(ix_max) ) then
            ix(i,j) = ix_max
          else
            ix(i,j) = locate( volx_out(i,j), ix_max, volx_c_chim ) - 1
          end if
          if ( voly_out(i,j) <= voly_c_chim(1) ) then
            iy(i,j) = 0
          else if ( voly_out(i,j) >= voly_c_chim(iy_max) ) then
            iy(i,j) = iy_max
          else
            iy(i,j) = locate( voly_out(i,j), iy_max, voly_c_chim ) - 1
          end if
        end do
      end do

      call interp2d_linear( ix, ix_max, iy, iy_max, volx_c_chim, voly_c_chim, state_chim, &
      &                     volx_out, voly_out, state_out )

    else if ( interp_method == 2 ) then
      call interp2d_spline( volx_c_chim(1:ix_max), voly_c_chim(1:iy_max), state_chim(1:ix_max,1:iy_max), &
      &                     volx_out, voly_out, state_out )
    else
      call bl_error("invalid value for interp_method")
    end if

    return
  end subroutine interp2d_chimera

  subroutine interp2d_quad_chimera( x_out, y_out, state_chim, state_out )

    use bl_constants_module
    use bl_error_module
    use model_interp_module, only: interp2d_linear, interp2d_spline
    use interpolate_module, only: locate
    use probdata_module, only: interp_method, nquad
    use quadrature_module, only: wquad, quad_avg

    ! input variables
    real (rt), intent(in) :: x_out(:,:,:,:)
    real (rt), intent(in) :: y_out(:,:,:,:)
    real (rt), intent(in) :: state_chim(:,:)

    ! output variables
    real (rt), intent(out) :: state_out(size(x_out,3),size(x_out,4))

    ! local variables
    integer :: ix(size(x_out,1),size(x_out,2))
    integer :: ix_max

    integer :: iy(size(x_out,1),size(x_out,2))
    integer :: iy_max

    real (rt) :: volx_quad(size(x_out,1),size(x_out,2))
    real (rt) :: voly_quad(size(x_out,1),size(x_out,2))


    real (rt) :: dvol1, dvol2
    real (rt) :: state_quad(size(x_out,1),size(x_out,1))

    integer :: i, ii, j, jj, n

    ix_max = imax_chim
    iy_max = jmax_chim

    if ( interp_method == 1 ) then

      do j = 1, size(x_out,4)
        do i = 1, size(x_out,3)
          volx_quad   = third * x_out(:,:,i,j)**3
          voly_quad = one - cos( y_out(:,:,i,j) )

          ix(:,:) = 0
          iy(:,:) = 0
          do jj = 1, nquad
            do ii = 1, nquad
              if ( volx_quad(ii,jj) <= volx_c_chim(1) ) then
                ix(ii,jj) = 0
              else if ( volx_quad(ii,jj) >= volx_c_chim(ix_max) ) then
                ix(ii,jj) = ix_max
              else
                ix(ii,jj) = locate( volx_quad(ii,jj), ix_max, volx_c_chim ) - 1
              end if
              if ( voly_quad(ii,jj) <= voly_c_chim(1) ) then
                iy(ii,jj) = 0
              else if ( voly_quad(ii,jj) >= voly_c_chim(iy_max) ) then
                iy(ii,jj) = iy_max
              else
                iy(ii,jj) = locate( voly_quad(ii,jj), iy_max, voly_c_chim ) - 1
              end if
            end do
          end do

          call interp2d_linear( ix, ix_max, iy, iy_max, volx_c_chim, voly_c_chim, state_chim, &
          &                     volx_quad, voly_quad, state_quad )

          state_out(i,j) = quad_avg( wquad, state_quad )

        end do
      end do

    else if ( interp_method == 2 ) then
      do j = 1, size(x_out,4)
        do i = 1, size(x_out,3)
          volx_quad   = third * x_out(:,:,i,j)**3
          voly_quad = one - cos( y_out(:,:,i,j) )

          call interp2d_spline( volx_c_chim(1:ix_max), voly_c_chim(1:iy_max), state_chim(1:ix_max,1:iy_max), &
          &                     volx_quad, voly_quad, state_quad )

          state_out(i,j) = quad_avg( wquad, state_quad )

        end do
      end do
    else
      call bl_error("invalid value for interp_method")
    end if

    return
  end subroutine interp2d_quad_chimera

  subroutine interp3d_chimera( x_out, y_out, z_out, state_chim, state_out )

    use bl_constants_module
    use bl_error_module
    use interpolate_module, only: locate
    use model_interp_module, only: interp3d_linear, interp3d_spline
    use probdata_module, only: interp_method

    ! input variables
    real (rt), intent(in) :: x_out(:,:,:)
    real (rt), intent(in) :: y_out(:,:,:)
    real (rt), intent(in) :: z_out(:,:,:)
    real (rt), intent(in) :: state_chim(:,:,:)

    ! output variables
    real (rt), intent(out) :: state_out(size(x_out,1),size(x_out,2),size(x_out,3))

    ! local variables
    integer :: ix(size(x_out,1),size(x_out,2),size(x_out,3))
    integer :: ix_max

    integer :: iy(size(x_out,1),size(x_out,2),size(x_out,3))
    integer :: iy_max

    integer :: iz(size(x_out,1),size(x_out,2),size(x_out,3))
    integer :: iz_max

    real (rt) :: volx_out(size(x_out,1),size(x_out,2),size(x_out,3))
    real (rt) :: voly_out(size(x_out,1),size(x_out,2),size(x_out,3))
    real (rt) :: volz_out(size(x_out,1),size(x_out,2),size(x_out,3))

    real (rt) :: dvol1, dvol2

    integer :: i, j, k, n

    volx_out   = third * x_out**3
    voly_out = one - cos( y_out )
    volz_out = z_out
    ix_max = imax_chim
    iy_max = jmax_chim
    iz_max = kmax_chim

    if ( interp_method == 1 ) then

      ix = 0
      iy = 0
      iz = 0
      do k = 1, size(x_out,3)
        do j = 1, size(x_out,2)
          do i = 1, size(x_out,1)
            if ( volx_out(i,j,k) <= volx_c_chim(1) ) then
              ix(i,j,k) = 0
            else if ( volx_out(i,j,k) >= volx_c_chim(ix_max) ) then
              ix(i,j,k) = ix_max
            else
              ix(i,j,k) = locate( volx_out(i,j,k), ix_max, volx_c_chim ) - 1
            end if
            
            if ( voly_out(i,j,k) <= voly_c_chim(1) ) then
              iy(i,j,k) = 0
            else if ( voly_out(i,j,k) >= voly_c_chim(iy_max) ) then
              iy(i,j,k) = iy_max
            else
              iy(i,j,k) = locate( voly_out(i,j,k), iy_max, voly_c_chim ) - 1
            end if

            if ( volz_out(i,j,k) <= volz_c_chim(1) ) then
              iz(i,j,k) = 0
            else if ( volz_out(i,j,k) >= volz_c_chim(iz_max) ) then
              iz(i,j,k) = iz_max
            else
              iz(i,j,k) = locate( volz_out(i,j,k), iz_max, volz_c_chim ) - 1
            end if
          end do
        end do
      end do

      call interp3d_linear( ix, ix_max, iy, iy_max, iz, iz_max, &
      &                     volx_c_chim, voly_c_chim, volz_c_chim, state_chim, &
      &                     volx_out, voly_out, volz_out, state_out )

!   else if ( interp_method == 2 ) then
!     call interp3d_spline( volx_c_chim(1:ix_max), voly_c_chim(1:iy_max), volz_c_chim(1:iz_max), &
!     &                     state_chim(1:ix_max,1:iy_max,1:iz_max), &
!     &                     volx_out, voly_out, volz_out, state_out )
    else
      call bl_error("invalid value for interp_method")
    end if

    return
  end subroutine interp3d_chimera

end module chimera_parser_module
