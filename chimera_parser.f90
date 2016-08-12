module chimera_parser_module

  ! we read in the number of variables and their order and use this to
  ! map them into the model_state array.  We ignore anything other than
  ! density, temperature, pressure and composition.
  !
  ! composition is assumed to be in terms of mass fractions     

  use parallel, only: parallel_IOProcessor
  use network
  use bl_types
  use hdf5_read_write
  use hdf5

  implicit none

  ! array dimensions
  integer, save :: nx_in, ny_in, nz_in, nnc_in, npart_in, npart_shell_in

  ! radial index bounds
  integer, save :: imin_in, imax_in

  ! theta index bounds
  integer, save :: jmin_in, jmax_in

  ! phi index bounds
  integer, save :: kmin_in, kmax_in

  ! grid coordinates
  real (dp_t), allocatable, save :: rad_edge_in(:), drad_edge_in(:), rad_cntr_in(:), drad_cntr_in(:)
  real (dp_t), allocatable, save :: theta_edge_in(:), dtheta_edge_in(:), theta_cntr_in(:), dtheta_cntr_in(:)
  real (dp_t), allocatable, save :: phi_edge_in(:), dphi_edge_in(:), phi_cntr_in(:), dphi_cntr_in(:)

  real (dp_t), allocatable, save :: vol_rad_edge_in(:), vol_rad_cntr_in(:)
  real (dp_t), allocatable, save :: vol_theta_edge_in(:), vol_theta_cntr_in(:)
  real (dp_t), allocatable, save :: vol_phi_edge_in(:), vol_phi_cntr_in(:)

  real (dp_t), allocatable, save :: dvol_rad_in(:), dvol_theta_in(:), dvol_phi_in(:)
  real (dp_t), allocatable, save :: volume_in(:,:,:), zone_mass_in(:,:,:), domega_in(:,:)

  ! time data
  real (dp_t), save :: time_in, t_bounce

  ! fluid velocities
  real (dp_t), allocatable, save :: vrad_in(:,:,:)
  real (dp_t), allocatable, save :: vtheta_in(:,:,:)
  real (dp_t), allocatable, save :: vphi_in(:,:,:)

  ! thermo data
  real (dp_t), allocatable, save :: dens_in(:,:,:)
  real (dp_t), allocatable, save :: temp_in(:,:,:)
  real (dp_t), allocatable, save :: ye_in(:,:,:)
  real (dp_t), allocatable, save :: pres_in(:,:,:)
  real (dp_t), allocatable, save :: eint_in(:,:,:)
  real (dp_t), allocatable, save :: ebind_in(:,:,:)
  real (dp_t), allocatable, save :: enpy_in(:,:,:)

  ! nuclei data
  real (dp_t), allocatable, save :: a_nuc_in(:)
  real (dp_t), allocatable, save :: z_nuc_in(:)
  real (dp_t), allocatable, save :: m_ex_nuc_in(:)
  real (dp_t), allocatable, save :: be_nuc_in(:)
  character(len=5), allocatable, save :: nuc_name_in(:)

  ! composition data
  real (dp_t), allocatable, save :: xn_read(:,:,:,:)
  real (dp_t), allocatable, save :: xn_in(:,:,:,:)
  integer,     allocatable, save :: nse_in(:,:,:)

  ! representative heavy/auxiliary data
  real (dp_t), allocatable, save :: a_nuc_rep_in(:,:,:)
  real (dp_t), allocatable, save :: z_nuc_rep_in(:,:,:)
  real (dp_t), allocatable, save :: be_nuc_rep_in(:,:,:)

  ! particle data
  real (dp_t), allocatable, save :: prad_in(:)
  real (dp_t), allocatable, save :: ptheta_in(:)
  real (dp_t), allocatable, save :: pphi_in(:)

  ! model_initialized will be .true. once the model is read in and the
  ! model data arrays are initialized and filled
  logical, save :: model_initialized = .false.

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

  subroutine open_chimera_file(model_file)

    use bl_constants_module
    use bl_error_module

    character(len=*), intent(in   ) :: model_file

    ! local variables
    integer :: i_read(1)

    integer :: ierr

    call h5open_f(ierr)
    call h5fopen_f(model_file, H5F_ACC_RDONLY_F, file_id, ierr)

    if (ierr /= 0) then
      if (parallel_IOProcessor()) print *,'Couldnt open model_file: ',model_file
      call bl_error('Aborting now -- please supply model_file')
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
    real (dp_t) :: d_read(1)

    integer :: array_dimensions(3)
    integer :: index_bounds(2)
    integer :: eos_dimensions(2)
    integer :: part_dimensions(2)

    integer :: i, j, k, ierr
    integer :: itmp(nspec)
    integer :: net_to_castro(nspec)
    integer, allocatable :: net_in_castro(:)

    real (dp_t) :: dx, dvol

    ! open mesh group
    call h5gopen_f(file_id, '/mesh', group_id, ierr)
    if (ierr /= 0) call bl_error('Could not open /mesh group')

    ! read array dimensions
    datasize1d(1) = 3
    call read_1d_slab('array_dimensions', array_dimensions, group_id, datasize1d)
    nx_in = array_dimensions(1)
    ny_in = array_dimensions(2)
    nz_in = array_dimensions(3)

    ! read radial index bounds
    datasize1d(1) = 2
    call read_1d_slab('radial_index_bound', index_bounds, group_id, datasize1d)
    imin_in = index_bounds(1)
    imax_in = index_bounds(2)

    ! read theta index bounds
    call read_1d_slab('theta_index_bound', index_bounds, group_id, datasize1d)
    jmin_in = index_bounds(1)
    jmax_in = index_bounds(2)
    if ( jmin_in /= 1 .or. jmax_in /= ny_in ) then
       if (parallel_IOProcessor()) print *,'Error: inner theta bounds /= (1,ny_in), (',jmin_in,',',jmax_in,')'
       call bl_error('Aborting now -- check theta bounds in model file')
     end if

    ! read phi index bounds
    call read_1d_slab('phi_index_bound', index_bounds, group_id, datasize1d)
    kmin_in = index_bounds(1)
    kmax_in = index_bounds(2)
    if ( kmin_in /= 1 .or. kmax_in /= nz_in ) then
       if (parallel_IOProcessor()) print *,'Error: inner phi bounds /= (1,nz_in), (',kmin_in,',',kmax_in,')'
       call bl_error('Aborting now -- check phi bounds in model file')
     end if

    ! allocate grid variables
    allocate (rad_edge_in(nx_in+1), drad_edge_in(nx_in+1), rad_cntr_in(nx_in), drad_cntr_in(nx_in))
    allocate (theta_edge_in(-5:ny_in+7), dtheta_edge_in(-5:ny_in+7), theta_cntr_in(-5:ny_in+6), dtheta_cntr_in(-5:ny_in+6))
    allocate (phi_edge_in(-5:nz_in+7), dphi_edge_in(-5:nz_in+7), phi_cntr_in(-5:nz_in+6), dphi_cntr_in(-5:nz_in+6))

    ! read zone edge coordinaets
    datasize1d(1) = nx_in+1
    call read_1d_slab('x_ef', rad_edge_in, group_id, datasize1d)
    datasize1d(1) = ny_in+1
    call read_1d_slab('y_ef', theta_edge_in(1:ny_in+1), group_id, datasize1d)
    datasize1d(1) = nz_in+1
    call read_1d_slab('z_ef', phi_edge_in(1:nz_in+1), group_id, datasize1d)

    ! read zone center coordinates
    datasize1d(1) = nx_in
    call read_1d_slab('x_cf', rad_cntr_in, group_id, datasize1d)
    datasize1d(1) = ny_in
    call read_1d_slab('y_cf', theta_cntr_in(1:ny_in), group_id, datasize1d)
    datasize1d(1) = nz_in
    call read_1d_slab('z_cf', phi_cntr_in(1:nz_in), group_id, datasize1d)

    ! read zone widths
    datasize1d(1) = nx_in
    call read_1d_slab('dx_cf', drad_cntr_in, group_id, datasize1d)
    datasize1d(1) = ny_in
    call read_1d_slab('dy_cf', dtheta_cntr_in(1:ny_in), group_id, datasize1d)
    datasize1d(1) = nz_in
    call read_1d_slab('dz_cf', dphi_cntr_in(1:nz_in), group_id, datasize1d)

    drad_edge_in(1:nx_in) = rad_edge_in(2:nx_in+1) - rad_edge_in(1:nx_in)
    drad_edge_in(nx_in+1) = drad_edge_in(nx_in)

    ! fill grid boundaries
    do i = 1, 6
      dtheta_cntr_in(jmin_in - i) = dtheta_cntr_in(jmin_in + i - 1)
      theta_edge_in(jmin_in - i) = theta_edge_in(jmin_in - i + 1) - dtheta_cntr_in(jmin_in - i)
      theta_cntr_in(jmin_in - i) = theta_edge_in(jmin_in - i) + half * dtheta_cntr_in(jmin_in - i)

      dphi_cntr_in(kmin_in - i) = dphi_cntr_in(kmax_in + 1 - i)
      phi_edge_in(kmin_in - i) = phi_edge_in(kmin_in - i + 1) - dphi_cntr_in(kmin_in - i)
      phi_cntr_in(kmin_in - i) = phi_edge_in(kmin_in - i) + half * dphi_cntr_in(kmin_in - i)
    end do

    do i = 1, 6
      dtheta_cntr_in(jmax_in + i) = dtheta_cntr_in(jmax_in + 1 - i)
      theta_edge_in(jmax_in + i + 1) = theta_edge_in(jmax_in + i) + dtheta_cntr_in(jmax_in + i)
      theta_cntr_in(jmax_in + i) = theta_edge_in(jmax_in + i) + half * dtheta_cntr_in(jmax_in + i)

      dphi_cntr_in(kmax_in + i) = dphi_cntr_in(kmin_in + i - 1)
      phi_edge_in(kmax_in + i + 1) = phi_edge_in(kmax_in + i) + dphi_cntr_in(kmax_in + i)
      phi_cntr_in(kmax_in + i) = phi_edge_in(kmax_in + i) + half * dphi_cntr_in(kmax_in + i)
    end do
    dtheta_edge_in(-5:ny_in+6) = theta_edge_in(-4:ny_in+7) - theta_edge_in(-5:ny_in+6)
    dtheta_edge_in(ny_in+7) = dtheta_edge_in(ny_in+6)
    dphi_edge_in(-5:nz_in+6) = phi_edge_in(-4:nz_in+7) - phi_edge_in(-5:nz_in+6)
    dphi_edge_in(nz_in+7) = dphi_edge_in(nz_in+6)

    ! read times
    datasize1d(1) = 1
    call read_1d_slab('time', d_read, group_id, datasize1d)
    time_in = d_read(1)
    call read_1d_slab('t_bounce', d_read, group_id, datasize1d)
    t_bounce = d_read(1)

    ! solid angles and volumes
    allocate (domega_in(ny_in,nz_in))
    allocate (volume_in(nx_in,ny_in,nz_in))
    allocate (vol_rad_edge_in(nx_in+1))
    allocate (vol_theta_edge_in(ny_in+1))
    allocate (vol_phi_edge_in(nz_in+1))
    allocate (vol_rad_cntr_in(nx_in))
    allocate (vol_theta_cntr_in(ny_in))
    allocate (vol_phi_cntr_in(nz_in))
    allocate (dvol_rad_in(nx_in+1))
    allocate (dvol_theta_in(ny_in+1))
    allocate (dvol_phi_in(nz_in+1))

    datasize2d = (/ ny_in, nz_in /)
    call read_2d_slab('d_omega',domega_in, group_id, datasize2d)

    dvol_rad_in(1)       = third * rad_edge_in(1)**3
    vol_rad_edge_in(1)   = dvol_rad_in(1)
    do i = 2, nx_in+1
      dx                 = rad_edge_in(i) - rad_edge_in(i-1)
      dvol               = dx * ( rad_edge_in(i-1) * rad_edge_in(i) + dx * dx * third )
      dvol_rad_in(i)     = dvol
      vol_rad_edge_in(i) = vol_rad_edge_in(i-1) + dvol
    end do
!   vol_rad_cntr_in(1)       = vol_rad_edge_in(1) - half * dvol_rad_in(1)
    vol_rad_cntr_in(1:nx_in) = vol_rad_edge_in(1:nx_in) + half*dvol_rad_in(2:nx_in+1)

    dvol_theta_in(1) = one - cos( theta_edge_in(1) )
    vol_theta_edge_in(1) = dvol_theta_in(1)
    do j = 2, ny_in+1
      dvol                 = cos( theta_edge_in(j-1) ) - cos( theta_edge_in(j) )
      dvol_theta_in(j)     = dvol
      vol_theta_edge_in(j) = vol_theta_edge_in(j-1) + dvol
    end do
!   vol_theta_cntr_in(1)       = vol_theta_edge_in(1) - half * dvol_theta_in(1)
    vol_theta_cntr_in(1:ny_in) = vol_theta_edge_in(1:ny_in) + half*dvol_theta_in(2:ny_in+1)

    dvol_phi_in(1) = phi_edge_in(1)
    vol_phi_edge_in(1) = dvol_phi_in(1)
    do k = 2, nz_in+1
      dvol               = dphi_edge_in(k)
      dvol_phi_in(k)     = dvol
      vol_phi_edge_in(k) = vol_phi_edge_in(k-1) + dvol
    end do
!   vol_phi_cntr_in(1) = vol_phi_edge_in(1) - half * dvol_phi_in(1)
    vol_phi_cntr_in(1:nz_in) = vol_phi_edge_in(1:nz_in) + half*dvol_phi_in(2:nz_in+1)

    do i = 1, nx_in
      volume_in(i,:,:) = dvol_rad_in(i) * domega_in(:,:)
    end do

    ! close mesh group
    call h5gclose_f(group_id, ierr)

    ! open fluid group
    call h5gopen_f(file_id, '/fluid', group_id, ierr)
    if (ierr /= 0) call bl_error('Could not open /fluid group')

    ! allocate fluid variables
    allocate (vrad_in(nx_in,ny_in,nz_in))
    allocate (vtheta_in(nx_in,ny_in,nz_in))
    allocate (vphi_in(nx_in,ny_in,nz_in))
    allocate (dens_in(nx_in,ny_in,nz_in))
    allocate (temp_in(nx_in,ny_in,nz_in))
    allocate (ye_in(nx_in,ny_in,nz_in))
    allocate (pres_in(nx_in,ny_in,nz_in))
    allocate (eint_in(nx_in,ny_in,nz_in))
    allocate (ebind_in(nx_in,ny_in,nz_in))
    allocate (enpy_in(nx_in,ny_in,nz_in))

    ! read velocity variables
    datasize3d = (/ nx_in, ny_in, nz_in /)
    slab_offset3d = (/ 0, 0, 0 /)
    call read_ray_hyperslab('u_c',     vrad_in,   group_id, datasize3d, slab_offset3d)
    call read_ray_hyperslab('v_c',     vtheta_in, group_id, datasize3d, slab_offset3d)
    call read_ray_hyperslab('w_c',     vphi_in,   group_id, datasize3d, slab_offset3d)

    ! read thermo variables
    call read_ray_hyperslab('rho_c',   dens_in,   group_id, datasize3d, slab_offset3d)
    call read_ray_hyperslab('t_c',     temp_in,   group_id, datasize3d, slab_offset3d)
    call read_ray_hyperslab('ye_c',    ye_in,     group_id, datasize3d, slab_offset3d)
    call read_ray_hyperslab('press',   pres_in,   group_id, datasize3d, slab_offset3d)
    call read_ray_hyperslab('e_int',   eint_in,   group_id, datasize3d, slab_offset3d)
    call read_ray_hyperslab('entropy', enpy_in,   group_id, datasize3d, slab_offset3d)

    allocate (zone_mass_in(nx_in,ny_in,nz_in))
    zone_mass_in(:,:,:) = dens_in(:,:,:) * volume_in(:,:,:)

    ! read number of nuclei
    datasize1d(1) = 2
    call read_1d_slab('dimeos', eos_dimensions, group_id, datasize1d)
    nnc_in = eos_dimensions(1)

    ! close fluid group
    call h5gclose_f(group_id, ierr)

    ! open abundance group
    call h5gopen_f(file_id, '/abundance', group_id, ierr)
    if (ierr /= 0) call bl_error('Could not open /abundance group')

    ! allocate composition variables
    allocate (a_nuc_in(nnc_in-1))
    allocate (z_nuc_in(nnc_in-1))
    allocate (m_ex_nuc_in(nnc_in-1))
    allocate (be_nuc_in(nnc_in-1))
    allocate (nuc_name_in(nnc_in))
    allocate (xn_read(nnc_in,nx_in,ny_in,nz_in))
    allocate (xn_in(nspec,nx_in,ny_in,nz_in))
    allocate (nse_in(nx_in+1,ny_in,nz_in))
    allocate (a_nuc_rep_in(nx_in,ny_in,nz_in))
    allocate (z_nuc_rep_in(nx_in,ny_in,nz_in))
    allocate (be_nuc_rep_in(nx_in,ny_in,nz_in))

    ! read nuclei values
    datasize1d(1) = nnc_in-1
    call read_1d_slab('a_nuc',    a_nuc_in,    group_id, datasize1d)
    call read_1d_slab('z_nuc',    z_nuc_in,    group_id, datasize1d)
    call read_1d_slab('m_ex_nuc', m_ex_nuc_in, group_id, datasize1d)

    ! binding energy
    be_nuc_in(1:nnc_in-1) = z_nuc_in(1:nnc_in-1)*7.2889705d+00 + &
    &                       (a_nuc_in(1:nnc_in-1)-z_nuc_in(1:nnc_in-1))*8.0713171d+00 - &
    &                       m_ex_nuc_in(1:nnc_in-1)
    
    ! read nuclei names
    datasize1d(1) = nnc_in
    call read_1d_slab('a_name',   nuc_name_in, 5, group_id, datasize1d)

    ! read mass fractions
    datasize4d = (/ nnc_in, nx_in, ny_in, nz_in /)
    slab_offset4d = (/ 0, 0, 0, 0 /)
    call read_ray_hyperslab('xn_c', xn_read, group_id, datasize4d, slab_offset4d)

    ! create lookup tables for isotopes in castro net
    k = 0
    net_to_castro(:) = nnc_in
    do i = 1, nspec
      do j = 1, nnc_in-1
        if ( nint(aion(i)) == nint(a_nuc_in(j)) .and. nint(zion(i)) == nint(z_nuc_in(j)) ) then
          k = k + 1
          net_to_castro(i) = j
          itmp(k) = i
          exit
        end if
      end do
      if ( j > nnc_in-1 ) then
        if (parallel_IOProcessor()) then
          write(*,'(2(a,i3),a)') ' could not find isotope (',nint(zion(i)),',',nint(aion(i)),') in CHIMERA net'
        end if
      end if
    end do

    if ( k > 0 ) then
      allocate (net_in_castro(k))
      net_in_castro(:) = itmp(1:k)
    else
      call bl_error("no species in CHIMERA net in CASTRO net")
    end if

    xn_in(net_in_castro,:,:,:) = xn_read(net_to_castro(net_in_castro),:,:,:)

    ! read nse state
    datasize3d = (/ nx_in+1, ny_in, nz_in /)
    slab_offset3d = (/ 0, 0, 0 /)
    call read_ray_hyperslab('nse_c', nse_in, group_id, datasize3d, slab_offset3d)

    ! read representative heavy nucleus/auxiliary nucleus properties
    datasize3d = (/ nx_in, ny_in, nz_in /)
    slab_offset3d = (/ 0, 0, 0 /)
    call read_ray_hyperslab('a_nuc_rep_c',  a_nuc_rep_in,  group_id, datasize3d, slab_offset3d)
    call read_ray_hyperslab('z_nuc_rep_c',  z_nuc_rep_in,  group_id, datasize3d, slab_offset3d)
    call read_ray_hyperslab('be_nuc_rep_c', be_nuc_rep_in, group_id, datasize3d, slab_offset3d)

    ! binding energy
    ebind_in(:,:,:) = 0.0d0
    do i = 1, nnc_in-1
      ebind_in(:,:,:) = ebind_in(:,:,:) - xn_read(i,:,:,:)*be_nuc_in(i)/a_nuc_in(i)
    end do
    ebind_in(:,:,:) = ebind_in(:,:,:) - xn_read(nnc_in,:,:,:)*be_nuc_rep_in(:,:,:)/a_nuc_rep_in(:,:,:)
    ! convert from MeV/c^2 to erg/g
    ebind_in(:,:,:) = ebind_in(:,:,:) * MeV2eV * ev2erg * n_A

    ! remove b.e. and constant offset from internal energy
    eint_in(:,:,:) = eint_in(:,:,:) - ebind_in(:,:,:) - 8.9d0 * MeV2eV * ev2erg * n_A

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
        npart_in = part_dimensions(1)
        npart_shell_in = part_dimensions(2)

        ! allocate particle variables
        allocate (prad_in(npart_in))
        allocate (ptheta_in(npart_in))
        allocate (pphi_in(npart_in))

        ! read particle variables
        datasize1d(1) = npart_in
        call read_1d_slab('px', prad_in, group_id, datasize1d)
        call read_1d_slab('py', ptheta_in, group_id, datasize1d)
        call read_1d_slab('pz', pphi_in, group_id, datasize1d)

        ! close particle group
        call h5gclose_f(group_id, ierr)
      end if
    end if

    ! close hdf5 file
    call h5fclose_f(file_id, ierr)
    call h5close_f(ierr)

    model_initialized = .true.

    return

  end subroutine read_chimera_file


  subroutine close_chimera_file
    
    if (model_initialized) then
       deallocate(rad_edge_in)
       deallocate(rad_cntr_in)
       deallocate(drad_cntr_in)
       deallocate(theta_edge_in)
       deallocate(theta_cntr_in)
       deallocate(dtheta_cntr_in)
       deallocate(phi_edge_in)
       deallocate(phi_cntr_in)
       deallocate(dphi_cntr_in)
       deallocate(vol_rad_cntr_in)
       deallocate(vol_theta_cntr_in)
       deallocate(vol_phi_cntr_in)
       deallocate(vol_rad_edge_in)
       deallocate(vol_theta_edge_in)
       deallocate(vol_phi_edge_in)
       deallocate(dvol_rad_in)
       deallocate(dvol_theta_in)
       deallocate(dvol_phi_in)
       deallocate(domega_in)
       deallocate(volume_in)
       deallocate(zone_mass_in)
       deallocate(vrad_in)
       deallocate(vtheta_in)
       deallocate(vphi_in)
       deallocate(dens_in)
       deallocate(temp_in)
       deallocate(ye_in)
       deallocate(pres_in)
       deallocate(eint_in)
       deallocate(ebind_in)
       deallocate(enpy_in)
       deallocate(a_nuc_in)
       deallocate(z_nuc_in)
       deallocate(m_ex_nuc_in)
       deallocate(be_nuc_in)
       deallocate(nuc_name_in)
       deallocate(xn_read)
       deallocate(xn_in)
       deallocate(nse_in)
       deallocate(a_nuc_rep_in)
       deallocate(z_nuc_rep_in)
       deallocate(be_nuc_rep_in)
       deallocate(prad_in)
       deallocate(ptheta_in)
       deallocate(pphi_in)
       model_initialized = .false.
    endif

    return

  end subroutine close_chimera_file

  subroutine interp1d_chimera( rad_out, state_in, state_out, interp_method )

    use bl_constants_module
    use bl_error_module
    use model_interp_module, only: interp1d_linear, interp1d_spline
    use interpolate_module, only: locate

    ! input variables
    real (dp_t), intent(in) :: rad_out(:)
    real (dp_t), intent(in) :: state_in(:)
    integer, intent(in) :: interp_method

    ! output variables
    real (dp_t), intent(out) :: state_out(size(rad_out))

    ! local variables
    integer :: irad(size(rad_out))
    integer :: irad_max

    real (dp_t) :: dvol1, dvol2
    real (dp_t) :: vol_rad_out(size(rad_out))

    integer :: i, n

    vol_rad_out = third * rad_out**3

    if ( interp_method == 1 ) then
      irad_max = imax_in
      irad(:) = 0
      do i = 1, size(rad_out)
        if ( vol_rad_out(i) <= vol_rad_cntr_in(1) ) then
          irad(i) = 0
        else if ( vol_rad_out(i) >= vol_rad_cntr_in(irad_max) ) then
          irad(i) = irad_max
        else
          irad(i) = locate( vol_rad_out(i), irad_max, vol_rad_cntr_in ) - 1
        end if
      end do
      call interp1d_linear( irad, irad_max, vol_rad_cntr_in, state_in, vol_rad_out, state_out )
    else if ( interp_method == 2 ) then
      call interp1d_spline( vol_rad_cntr_in, state_in, vol_rad_out, state_out )
    else
      call bl_error("invalid value for interp_method")
    end if

    return
  end subroutine interp1d_chimera

  subroutine interp2d_chimera( rad_out, theta_out, state_in, state_out, interp_method )

    use bl_constants_module
    use bl_error_module
    use model_interp_module, only: interp2d_linear, interp2d_spline
    use interpolate_module, only: locate

    ! input variables
    real (dp_t), intent(in) :: rad_out(:,:)
    real (dp_t), intent(in) :: theta_out(:,:)
    real (dp_t), intent(in) :: state_in(:,:)
    integer, intent(in) :: interp_method

    ! output variables
    real (dp_t), intent(out) :: state_out(size(rad_out,1),size(rad_out,2))

    ! local variables
    integer :: irad(size(rad_out,1),size(rad_out,2))
    integer :: irad_max

    integer :: itheta(size(rad_out,1),size(rad_out,2))
    integer :: itheta_max

    real (dp_t) :: vol_rad_out(size(rad_out,1),size(rad_out,2))
    real (dp_t) :: vol_theta_out(size(rad_out,1),size(rad_out,2))

    real (dp_t) :: dvol1, dvol2

    integer :: i, j, n

    vol_rad_out   = third * rad_out**3
    vol_theta_out = one - cos( theta_out )
    irad_max = imax_in
    itheta_max = jmax_in

    if ( interp_method == 1 ) then

      irad(:,:) = 0
      itheta(:,:) = 0
      do j = 1, size(rad_out,2)
        do i = 1, size(rad_out,1)
          if ( vol_rad_out(i,j) <= vol_rad_cntr_in(1) ) then
            irad(i,j) = 0
          else if ( vol_rad_out(i,j) >= vol_rad_cntr_in(irad_max) ) then
            irad(i,j) = irad_max
          else
            irad(i,j) = locate( vol_rad_out(i,j), irad_max, vol_rad_cntr_in ) - 1
          end if
          if ( vol_theta_out(i,j) <= vol_theta_cntr_in(1) ) then
            itheta(i,j) = 0
          else if ( vol_theta_out(i,j) >= vol_theta_cntr_in(itheta_max) ) then
            itheta(i,j) = itheta_max
          else
            itheta(i,j) = locate( vol_theta_out(i,j), itheta_max, vol_theta_cntr_in ) - 1
          end if
        end do
      end do

      call interp2d_linear( irad, irad_max, itheta, itheta_max, vol_rad_cntr_in, vol_theta_cntr_in, state_in, &
      &                     vol_rad_out, vol_theta_out, state_out )

    else if ( interp_method == 2 ) then
      call interp2d_spline( vol_rad_cntr_in(1:irad_max), vol_theta_cntr_in(1:itheta_max), state_in(1:irad_max,1:itheta_max), &
      &                     vol_rad_out, vol_theta_out, state_out )
    else
      call bl_error("invalid value for interp_method")
    end if

    return
  end subroutine interp2d_chimera

end module chimera_parser_module
