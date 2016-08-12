module hdf5_read_write

  use bl_types
  use hdf5
  
  implicit none
  
  logical, public :: hyperslab_group_master
  
  public      :: write_ray_hyperslab
  private     :: write_ray_hyperslab_dbl_2d
  private     :: write_ray_hyperslab_dbl_3d
  private     :: write_ray_hyperslab_dbl_4d
  private     :: write_ray_hyperslab_dbl_5d
  private     :: write_ray_hyperslab_int_3d
  
  public      :: write_1d_slab
  private     :: write_1d_slab_int
  private     :: write_1d_slab_double
  private     :: write_1d_slab_string
  public      :: write_2d_slab
  private     :: write_2d_slab_double
  private     :: write_2d_slab_int
  
  public      :: read_ray_hyperslab
  private     :: read_ray_hyperslab_dbl_2d
  private     :: read_ray_hyperslab_dbl_3d
  private     :: read_ray_hyperslab_dbl_4d
  private     :: read_ray_hyperslab_dbl_5d
  private     :: read_ray_hyperslab_int_3d
  
  public      :: read_1d_slab
  private     :: read_1d_slab_int
  private     :: read_1d_slab_double
  private     :: read_1d_slab_string
  
  public      :: read_2d_slab
  private     :: read_2d_slab_int
  private     :: read_2d_slab_double
  
  interface write_ray_hyperslab
    module procedure write_ray_hyperslab_int_3d
    module procedure write_ray_hyperslab_dbl_2d
    module procedure write_ray_hyperslab_dbl_3d
    module procedure write_ray_hyperslab_dbl_4d
    module procedure write_ray_hyperslab_dbl_5d
  end interface write_ray_hyperslab
  
  interface write_1d_slab
    module procedure write_1d_slab_int
    module procedure write_1d_slab_double
    module procedure write_1d_slab_string
  end interface write_1d_slab

  interface write_2d_slab
    module procedure write_2d_slab_double
    module procedure write_2d_slab_int
  end interface write_2d_slab

  interface read_ray_hyperslab
    module procedure read_ray_hyperslab_int_3d
    module procedure read_ray_hyperslab_dbl_2d
    module procedure read_ray_hyperslab_dbl_3d
    module procedure read_ray_hyperslab_dbl_4d
    module procedure read_ray_hyperslab_dbl_5d
  end interface read_ray_hyperslab
  
  interface read_1d_slab
    module procedure read_1d_slab_int
    module procedure read_1d_slab_double
    module procedure read_1d_slab_string
  end interface read_1d_slab

  interface read_2d_slab
    module procedure read_2d_slab_int
    module procedure read_2d_slab_double
  end interface read_2d_slab


  contains

  subroutine write_1d_slab_int(name, value, group_id, datasize, &
               desc_option, unit_option)
    character(*), intent(in)                    :: name
    character(*), intent(in), optional          :: unit_option
    character(*), intent(in), optional          :: desc_option
    integer(hid_t)                              :: group_id
    integer(hsize_t), dimension(1), intent(in)  :: datasize
    integer, dimension(:), intent(in)           :: value
    
    integer(hid_t)                              :: dataset_id
    integer(hid_t)                              :: dataspace_id
    integer(hid_t)                              :: atype_id
    integer(hid_t)                              :: attr_id
    integer(size_t)                             :: attr_len
    integer(hsize_t), dimension(1)              :: adims = (/1/)
    integer                                     :: error
    
    call h5screate_simple_f(1, datasize, dataspace_id, error)
    if ( error /= 0 ) print *,'Could not create simple: '//name
    call h5dcreate_f(group_id, name, h5t_native_integer, &
                     dataspace_id, dataset_id, error)
    if(hyperslab_group_master) &
      call h5dwrite_f(dataset_id, h5t_native_integer, &
                      value, datasize, error)
    call h5sclose_f(dataspace_id, error)
    if(present(desc_option))then
      attr_len = len(desc_option)
      call h5screate_simple_f(1, adims, dataspace_id, error)
      if ( error /= 0 ) print *,'Could not create simple: '//name
      call h5tcopy_f(h5t_native_character, atype_id, error)
      call h5tset_size_f(atype_id, attr_len, error)
      if ( error /= 0 ) print *,'Could not set size for: '//name
      call h5acreate_f(dataset_id, 'desc', atype_id, dataspace_id, &
             attr_id, error)
      if(hyperslab_group_master) &
        call h5awrite_f(attr_id, atype_id, desc_option, adims, error)
      call h5aclose_f(attr_id, error)
      call h5sclose_f(dataspace_id, error)
    end if
    if(present(unit_option))then
      attr_len = len(unit_option)
      call h5screate_simple_f(1, adims, dataspace_id, error)
      if ( error /= 0 ) print *,'Could not create simple: '//name
      call h5tcopy_f(h5t_native_character, atype_id, error)
      call h5tset_size_f(atype_id, attr_len, error)
      if ( error /= 0 ) print *,'Could not set size for: '//name
      call h5acreate_f(dataset_id, 'unit', atype_id, dataspace_id, &
             attr_id, error)
      if(hyperslab_group_master) &
        call h5awrite_f(attr_id, atype_id, unit_option, adims, error)
      call h5aclose_f(attr_id, error)
      call h5sclose_f(dataspace_id, error)
    end if
    call h5dclose_f(dataset_id, error)

  end subroutine write_1d_slab_int

  subroutine write_1d_slab_double(name, value, group_id, datasize, &
               desc_option, unit_option)
    character(*), intent(in)                    :: name
    character(*), intent(in), optional          :: unit_option
    character(*), intent(in), optional          :: desc_option
    integer(hid_t)                              :: group_id
    integer(hsize_t), dimension(1), intent(in)  :: datasize
    real(dp_t), dimension(:), intent(in)        :: value
    
    integer(hid_t)                              :: dataset_id
    integer(hid_t)                              :: dataspace_id
    integer(hid_t)                              :: atype_id
    integer(hid_t)                              :: attr_id
    integer(size_t)                             :: attr_len
    integer(hsize_t), dimension(1)              :: adims = (/1/)
    integer                                     :: error
    
    call h5screate_simple_f(1, datasize, dataspace_id, error)
    if ( error /= 0 ) print *,'Could not create simple: '//name
    call h5dcreate_f(group_id, name, h5t_native_double, &
                     dataspace_id, dataset_id, error)
    if(hyperslab_group_master) &
      call h5dwrite_f(dataset_id, h5t_native_double, &
                      value, datasize, error)
    call h5sclose_f(dataspace_id, error)
    if(present(desc_option))then
      attr_len = len(desc_option)
      call h5screate_simple_f(1, adims, dataspace_id, error)
      if ( error /= 0 ) print *,'Could not create simple: '//name
      call h5tcopy_f(h5t_native_character, atype_id, error)
      call h5tset_size_f(atype_id, attr_len, error)
      if ( error /= 0 ) print *,'Could not set size for: '//name
      call h5acreate_f(dataset_id, 'desc', atype_id, dataspace_id, &
             attr_id, error)
      if(hyperslab_group_master) &
        call h5awrite_f(attr_id, atype_id, desc_option, adims, error)
      call h5aclose_f(attr_id, error)
      call h5sclose_f(dataspace_id, error)
    end if
    if(present(unit_option))then
      attr_len = len(unit_option)
      call h5screate_simple_f(1, adims, dataspace_id, error)
      if ( error /= 0 ) print *,'Could not create simple: '//name
      call h5tcopy_f(h5t_native_character, atype_id, error)
      call h5tset_size_f(atype_id, attr_len, error)
      if ( error /= 0 ) print *,'Could not set size for: '//name
      call h5acreate_f(dataset_id, 'unit', atype_id, dataspace_id, &
             attr_id, error)
      if(hyperslab_group_master) &
        call h5awrite_f(attr_id, atype_id, unit_option, adims, error)
      call h5aclose_f(attr_id, error)
      call h5sclose_f(dataspace_id, error)
    end if
    call h5dclose_f(dataset_id, error)

  end subroutine write_1d_slab_double

  subroutine write_1d_slab_string(name, value, strlen, group_id, datasize, &
               desc_option, unit_option)
    character(*), intent(in)                    :: name
    character(*), intent(in), optional          :: unit_option
    character(*), intent(in), optional          :: desc_option
    integer                                     :: strlen
    integer(hsize_t)                            :: sizechar
    integer(hid_t)                              :: group_id
    integer(hsize_t), dimension(1), intent(in)  :: datasize
    character(len=strlen), dimension(:), intent(in) :: value

    integer(hid_t)                              :: dataset_id
    integer(hid_t)                              :: dataspace_id
    integer(hid_t)                              :: atype_id
    integer(hid_t)                              :: attr_id
    integer(size_t)                             :: attr_len
    integer(hsize_t), dimension(1)              :: adims = (/1/)
    integer                                     :: error

    call h5screate_simple_f(1, datasize, dataspace_id, error)
    if ( error /= 0 ) print *,'Could not create simple: '//name
    sizechar = strlen
    call h5tset_size_f(h5t_native_character, sizechar, error)
    if ( error /= 0 ) print *,'Could not set size for: '//name
    call h5dcreate_f(group_id, name, h5t_native_character, &
                     dataspace_id, dataset_id, error)
    if(hyperslab_group_master) &
      call h5dwrite_f(dataset_id, h5t_native_character, &
                      value, datasize, error)
    call h5sclose_f(dataspace_id, error)
    if(present(desc_option))then
      attr_len = len(desc_option)
      call h5screate_simple_f(1, adims, dataspace_id, error)
      if ( error /= 0 ) print *,'Could not create simple: '//name
      call h5tcopy_f(h5t_native_character, atype_id, error)
      call h5tset_size_f(atype_id, attr_len, error)
      if ( error /= 0 ) print *,'Could not set size for: '//name
      call h5acreate_f(dataset_id, 'desc', atype_id, dataspace_id, &
             attr_id, error)
      if(hyperslab_group_master) &
        call h5awrite_f(attr_id, atype_id, desc_option, adims, error)
      call h5aclose_f(attr_id, error)
      call h5sclose_f(dataspace_id, error)
    end if
    if(present(unit_option))then
      attr_len = len(unit_option)
      call h5screate_simple_f(1, adims, dataspace_id, error)
      if ( error /= 0 ) print *,'Could not create simple: '//name
      call h5tcopy_f(h5t_native_character, atype_id, error)
      call h5tset_size_f(atype_id, attr_len, error)
      if ( error /= 0 ) print *,'Could not set size for: '//name
      call h5acreate_f(dataset_id, 'unit', atype_id, dataspace_id, &
             attr_id, error)
      if(hyperslab_group_master) &
        call h5awrite_f(attr_id, atype_id, unit_option, adims, error)
      call h5aclose_f(attr_id, error)
      call h5sclose_f(dataspace_id, error)
    end if
    call h5dclose_f(dataset_id, error)

  end subroutine write_1d_slab_string

  subroutine write_2d_slab_double(name, value, group_id, datasize, &
               desc_option, unit_option)
    character(*), intent(in)                    :: name
    character(*), intent(in), optional          :: unit_option
    character(*), intent(in), optional          :: desc_option
    integer(hid_t)                              :: group_id
    integer(hsize_t), dimension(2), intent(in)  :: datasize
    real(dp_t), dimension(:,:), intent(in)      :: value
    
    integer(hid_t)                              :: dataset_id
    integer(hid_t)                              :: dataspace_id
    integer(hid_t)                              :: atype_id
    integer(hid_t)                              :: attr_id
    integer(size_t)                             :: attr_len
    integer(hsize_t), dimension(1)              :: adims = (/1/)
    integer                                     :: error
    
    call h5screate_simple_f(2, datasize, dataspace_id, error)
    if ( error /= 0 ) print *,'Could not create simple: '//name
    call h5dcreate_f(group_id, name, h5t_native_double, &
                     dataspace_id, dataset_id, error)
    if(hyperslab_group_master) &
      call h5dwrite_f(dataset_id, h5t_native_double, &
                      value, datasize, error)
    call h5sclose_f(dataspace_id, error)
    if(present(desc_option))then
      attr_len = len(desc_option)
      call h5screate_simple_f(1, adims, dataspace_id, error)
      if ( error /= 0 ) print *,'Could not create simple: '//name
      call h5tcopy_f(h5t_native_character, atype_id, error)
      call h5tset_size_f(atype_id, attr_len, error)
      if ( error /= 0 ) print *,'Could not set size for: '//name
      call h5acreate_f(dataset_id, 'desc', atype_id, dataspace_id, &
             attr_id, error)
      if(hyperslab_group_master) &
        call h5awrite_f(attr_id, atype_id, desc_option, adims, error)
      call h5aclose_f(attr_id, error)
      call h5sclose_f(dataspace_id, error)
    end if
    if(present(unit_option))then
      attr_len = len(unit_option)
      call h5screate_simple_f(1, adims, dataspace_id, error)
      if ( error /= 0 ) print *,'Could not create simple: '//name
      call h5tcopy_f(h5t_native_character, atype_id, error)
      call h5tset_size_f(atype_id, attr_len, error)
      if ( error /= 0 ) print *,'Could not set size for: '//name
      call h5acreate_f(dataset_id, 'unit', atype_id, dataspace_id, &
             attr_id, error)
      if(hyperslab_group_master) &
        call h5awrite_f(attr_id, atype_id, unit_option, adims, error)
      call h5aclose_f(attr_id, error)
      call h5sclose_f(dataspace_id, error)
    end if
    call h5dclose_f(dataset_id, error)

  end subroutine write_2d_slab_double

  subroutine write_2d_slab_int(name, value, group_id, datasize, &
               desc_option, unit_option)
    character(*), intent(in)                    :: name
    character(*), intent(in), optional          :: unit_option
    character(*), intent(in), optional          :: desc_option
    integer(hid_t)                              :: group_id
    integer(hsize_t), dimension(2), intent(in)  :: datasize
    integer, dimension(:,:), intent(in)         :: value
    
    integer(hid_t)                              :: dataset_id
    integer(hid_t)                              :: dataspace_id
    integer(hid_t)                              :: atype_id
    integer(hid_t)                              :: attr_id
    integer(size_t)                             :: attr_len
    integer(hsize_t), dimension(1)              :: adims = (/1/)
    integer                                     :: error
    
    call h5screate_simple_f(2, datasize, dataspace_id, error)
    if ( error /= 0 ) print *,'Could not create simple: '//name
    call h5dcreate_f(group_id, name, h5t_native_integer, &
                     dataspace_id, dataset_id, error)
    if(hyperslab_group_master) &
      call h5dwrite_f(dataset_id, h5t_native_integer, &
                      value, datasize, error)
    call h5sclose_f(dataspace_id, error)
    if(present(desc_option))then
      attr_len = len(desc_option)
      call h5screate_simple_f(1, adims, dataspace_id, error)
      if ( error /= 0 ) print *,'Could not create simple: '//name
      call h5tcopy_f(h5t_native_character, atype_id, error)
      call h5tset_size_f(atype_id, attr_len, error)
      if ( error /= 0 ) print *,'Could not set size for: '//name
      call h5acreate_f(dataset_id, 'desc', atype_id, dataspace_id, &
             attr_id, error)
      if(hyperslab_group_master) &
        call h5awrite_f(attr_id, atype_id, desc_option, adims, error)
      call h5aclose_f(attr_id, error)
      call h5sclose_f(dataspace_id, error)
    end if
    if(present(unit_option))then
      attr_len = len(unit_option)
      call h5screate_simple_f(1, adims, dataspace_id, error)
      if ( error /= 0 ) print *,'Could not create simple: '//name
      call h5tcopy_f(h5t_native_character, atype_id, error)
      call h5tset_size_f(atype_id, attr_len, error)
      if ( error /= 0 ) print *,'Could not set size for: '//name
      call h5acreate_f(dataset_id, 'unit', atype_id, dataspace_id, &
             attr_id, error)
      if(hyperslab_group_master) &
        call h5awrite_f(attr_id, atype_id, unit_option, adims, error)
      call h5aclose_f(attr_id, error)
      call h5sclose_f(dataspace_id, error)
    end if
    call h5dclose_f(dataset_id, error)

  end subroutine write_2d_slab_int



  subroutine write_ray_hyperslab_dbl_2d(name, value, group_id, &
               global_datasize, local_datasize, slab_offset, &
               desc_option, unit_option)
    
    character(*), intent(in)                    :: name
    character(*), intent(in), optional          :: unit_option
    character(*), intent(in), optional          :: desc_option
    integer(hid_t)                              :: group_id
    integer(hsize_t), dimension(2), intent(in)  :: global_datasize
    integer(hsize_t), dimension(2), intent(in)  :: local_datasize
    integer(hsize_t), dimension(2), intent(in)  :: slab_offset
    real(dp_t), dimension(:,:), intent(in)      :: value
    
    integer(hid_t)                              :: filespace    
    integer(hid_t)                              :: memspace    
    integer(hid_t)                              :: plist_id
    integer(hid_t)                              :: dataset_id
    integer(hid_t)                              :: dataspace_id
    integer(hid_t)                              :: atype_id
    integer(hid_t)                              :: attr_id
    integer(size_t)                             :: attr_len
    integer(hsize_t), dimension(1)              :: adims = (/1/)
    integer                                     :: error
    
    call h5screate_simple_f(2, global_datasize, filespace, error)
    if ( error /= 0 ) print *,'Could not create simple: '//name
    call h5dcreate_f(group_id, name, h5t_native_double, filespace, &
                     dataset_id, error)
    call h5sclose_f(filespace, error)
    
    call h5screate_simple_f(2, local_datasize, memspace, error)
    if ( error /= 0 ) print *,'Could not create simple: '//name
    call h5dget_space_f(dataset_id, filespace, error)
    if ( error /= 0 ) print *,'Could not get space: '//name
    call h5sselect_hyperslab_f(filespace, h5s_select_set_f, slab_offset, &
           local_datasize, error)
    if ( error /= 0 ) print *,'Could not select hyperslab: '//name
    call h5pcreate_f(h5p_dataset_xfer_f, plist_id, error) 

!    select case(hdf_mpio_mode)
!      case (1)
!        call h5pset_dxpl_mpio_f(plist_id, h5fd_mpio_independent_f, error)
!      case (2)
        call h5pset_dxpl_mpio_f(plist_id, h5fd_mpio_collective_f, error)
!      case default
!        call h5pset_dxpl_mpio_f(plist_id, h5fd_mpio_independent_f, error)
!    end select

    call h5dwrite_f(dataset_id, h5t_native_double, value, local_datasize, &
           error, file_space_id=filespace, mem_space_id=memspace, &
           xfer_prp=plist_id)
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)
    if(present(desc_option))then
      attr_len = len(desc_option)
      call h5screate_simple_f(1, adims, dataspace_id, error)
      if ( error /= 0 ) print *,'Could not create simple: '//name
      call h5tcopy_f(h5t_native_character, atype_id, error)
      call h5tset_size_f(atype_id, attr_len, error)
      if ( error /= 0 ) print *,'Could not set size for: '//name
      call h5acreate_f(dataset_id, 'desc', atype_id, dataspace_id, &
             attr_id, error)
      if(hyperslab_group_master) &
        call h5awrite_f(attr_id, atype_id, desc_option, adims, error)
      call h5aclose_f(attr_id, error)
      call h5sclose_f(dataspace_id, error)
    end if
    if(present(unit_option))then
      attr_len = len(unit_option)
      call h5screate_simple_f(1, adims, dataspace_id, error)
      if ( error /= 0 ) print *,'Could not create simple: '//name
      call h5tcopy_f(h5t_native_character, atype_id, error)
      call h5tset_size_f(atype_id, attr_len, error)
      if ( error /= 0 ) print *,'Could not set size for: '//name
      call h5acreate_f(dataset_id, 'unit', atype_id, dataspace_id, &
             attr_id, error)
      if(hyperslab_group_master) &
        call h5awrite_f(attr_id, atype_id, unit_option, adims, error)
      call h5aclose_f(attr_id, error)
      call h5sclose_f(dataspace_id, error)
    end if
    call h5dclose_f(dataset_id, error)
    call h5pclose_f(plist_id, error)                     
    
  end subroutine write_ray_hyperslab_dbl_2d


  subroutine write_ray_hyperslab_dbl_3d(name, value, group_id, &
               global_datasize, local_datasize, slab_offset, &
               desc_option, unit_option)
    
    character(*), intent(in)                    :: name
    character(*), intent(in), optional          :: unit_option
    character(*), intent(in), optional          :: desc_option
    integer(hid_t)                              :: group_id
    integer(hsize_t), dimension(3), intent(in)  :: global_datasize
    integer(hsize_t), dimension(3), intent(in)  :: local_datasize
    integer(hsize_t), dimension(3), intent(in)  :: slab_offset
    real(dp_t), dimension(:,:,:), intent(in)    :: value
    
    integer(hid_t)                              :: filespace    
    integer(hid_t)                              :: memspace    
    integer(hid_t)                              :: plist_id
    integer(hid_t)                              :: dataset_id
    integer(hid_t)                              :: dataspace_id
    integer(hid_t)                              :: atype_id
    integer(hid_t)                              :: attr_id
    integer(size_t)                             :: attr_len
    integer(hsize_t), dimension(1)              :: adims = (/1/)
    integer                                     :: error
    
    call h5screate_simple_f(3, global_datasize, filespace, error)
    if ( error /= 0 ) print *,'Could not create simple: '//name
    call h5dcreate_f(group_id, name, h5t_native_double, filespace, &
                     dataset_id, error)
    call h5sclose_f(filespace, error)
    
    call h5screate_simple_f(3, local_datasize, memspace, error)
    if ( error /= 0 ) print *,'Could not create simple: '//name
    call h5dget_space_f(dataset_id, filespace, error)
    if ( error /= 0 ) print *,'Could not get space: '//name
    call h5sselect_hyperslab_f(filespace, h5s_select_set_f, slab_offset, &
           local_datasize, error)
    if ( error /= 0 ) print *,'Could not select hyperslab: '//name
    call h5pcreate_f(h5p_dataset_xfer_f, plist_id, error) 

!    select case(hdf_mpio_mode)
!      case (1)
!        call h5pset_dxpl_mpio_f(plist_id, h5fd_mpio_independent_f, error)
!      case (2)
        call h5pset_dxpl_mpio_f(plist_id, h5fd_mpio_collective_f, error)
!      case default
!        call h5pset_dxpl_mpio_f(plist_id, h5fd_mpio_independent_f, error)
!    end select

    call h5dwrite_f(dataset_id, h5t_native_double, value, local_datasize, &
           error, file_space_id=filespace, mem_space_id=memspace, &
           xfer_prp=plist_id)
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)
    if(present(desc_option))then
      attr_len = len(desc_option)
      call h5screate_simple_f(1, adims, dataspace_id, error)
      if ( error /= 0 ) print *,'Could not create simple: '//name
      call h5tcopy_f(h5t_native_character, atype_id, error)
      call h5tset_size_f(atype_id, attr_len, error)
      if ( error /= 0 ) print *,'Could not set size for: '//name
      call h5acreate_f(dataset_id, 'desc', atype_id, dataspace_id, &
             attr_id, error)
      if(hyperslab_group_master) &
        call h5awrite_f(attr_id, atype_id, desc_option, adims, error)
      call h5aclose_f(attr_id, error)
      call h5sclose_f(dataspace_id, error)
    end if
    if(present(unit_option))then
      attr_len = len(unit_option)
      call h5screate_simple_f(1, adims, dataspace_id, error)
      if ( error /= 0 ) print *,'Could not create simple: '//name
      call h5tcopy_f(h5t_native_character, atype_id, error)
      call h5tset_size_f(atype_id, attr_len, error)
      if ( error /= 0 ) print *,'Could not set size for: '//name
      call h5acreate_f(dataset_id, 'unit', atype_id, dataspace_id, &
             attr_id, error)
      if(hyperslab_group_master) &
        call h5awrite_f(attr_id, atype_id, unit_option, adims, error)
      call h5aclose_f(attr_id, error)
      call h5sclose_f(dataspace_id, error)
    end if
    call h5dclose_f(dataset_id, error)
    call h5pclose_f(plist_id, error)                     
    
  end subroutine write_ray_hyperslab_dbl_3d


  subroutine write_ray_hyperslab_dbl_4d(name, value, group_id, &
               global_datasize, local_datasize, slab_offset, &
               desc_option, unit_option)
    
    character(*), intent(in)                    :: name
    character(*), intent(in), optional          :: unit_option
    character(*), intent(in), optional          :: desc_option
    integer(hid_t)                              :: group_id
    integer(hsize_t), dimension(4), intent(in)  :: global_datasize
    integer(hsize_t), dimension(4), intent(in)  :: local_datasize
    integer(hsize_t), dimension(4), intent(in)  :: slab_offset
    real(dp_t), dimension(:,:,:,:), intent(in)  :: value
    
    integer(hid_t)                              :: filespace    
    integer(hid_t)                              :: memspace    
    integer(hid_t)                              :: plist_id
    integer(hid_t)                              :: dataset_id
    integer(hid_t)                              :: dataspace_id
    integer(hid_t)                              :: atype_id
    integer(hid_t)                              :: attr_id
    integer(size_t)                             :: attr_len
    integer(hsize_t), dimension(1)              :: adims = (/1/)
    integer                                     :: error
    
    call h5screate_simple_f(4, global_datasize, filespace, error)
    if ( error /= 0 ) print *,'Could not create simple: '//name
    call h5dcreate_f(group_id, name, h5t_native_double, filespace, &
                     dataset_id, error)
    call h5sclose_f(filespace, error)
    
    call h5screate_simple_f(4, local_datasize, memspace, error)
    if ( error /= 0 ) print *,'Could not create simple: '//name
    call h5dget_space_f(dataset_id, filespace, error)
    if ( error /= 0 ) print *,'Could not get space: '//name
    call h5sselect_hyperslab_f(filespace, h5s_select_set_f, slab_offset, &
           local_datasize, error)
    if ( error /= 0 ) print *,'Could not select hyperslab: '//name
    call h5pcreate_f(h5p_dataset_xfer_f, plist_id, error) 

!    select case(hdf_mpio_mode)
!      case (1)
!        call h5pset_dxpl_mpio_f(plist_id, h5fd_mpio_independent_f, error)
!      case (2)
        call h5pset_dxpl_mpio_f(plist_id, h5fd_mpio_collective_f, error)
!      case default
!        call h5pset_dxpl_mpio_f(plist_id, h5fd_mpio_independent_f, error)
!    end select

    call h5dwrite_f(dataset_id, h5t_native_double, value, local_datasize, &
           error, file_space_id=filespace, mem_space_id=memspace, &
           xfer_prp=plist_id)
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)
    if(present(desc_option))then
      attr_len = len(desc_option)
      call h5screate_simple_f(1, adims, dataspace_id, error)
      if ( error /= 0 ) print *,'Could not create simple: '//name
      call h5tcopy_f(h5t_native_character, atype_id, error)
      call h5tset_size_f(atype_id, attr_len, error)
      if ( error /= 0 ) print *,'Could not set size for: '//name
      call h5acreate_f(dataset_id, 'desc', atype_id, dataspace_id, &
             attr_id, error)
      if(hyperslab_group_master) &
        call h5awrite_f(attr_id, atype_id, desc_option, adims, error)
      call h5aclose_f(attr_id, error)
      call h5sclose_f(dataspace_id, error)
    end if
    if(present(unit_option))then
      attr_len = len(unit_option)
      call h5screate_simple_f(1, adims, dataspace_id, error)
      if ( error /= 0 ) print *,'Could not create simple: '//name
      call h5tcopy_f(h5t_native_character, atype_id, error)
      call h5tset_size_f(atype_id, attr_len, error)
      if ( error /= 0 ) print *,'Could not set size for: '//name
      call h5acreate_f(dataset_id, 'unit', atype_id, dataspace_id, &
             attr_id, error)
      if(hyperslab_group_master) &
        call h5awrite_f(attr_id, atype_id, unit_option, adims, error)
      call h5aclose_f(attr_id, error)
      call h5sclose_f(dataspace_id, error)
    end if
    call h5dclose_f(dataset_id, error)
    call h5pclose_f(plist_id, error)                     
    
  end subroutine write_ray_hyperslab_dbl_4d


  subroutine write_ray_hyperslab_dbl_5d(name, value, group_id, &
               global_datasize, local_datasize, slab_offset, &
               desc_option, unit_option)
    
    character(*), intent(in)                    :: name
    character(*), intent(in), optional          :: unit_option
    character(*), intent(in), optional          :: desc_option
    integer(hid_t)                              :: group_id
    integer(hsize_t), dimension(5), intent(in)  :: global_datasize
    integer(hsize_t), dimension(5), intent(in)  :: local_datasize
    integer(hsize_t), dimension(5), intent(in)  :: slab_offset
    real(dp_t), dimension(:,:,:,:,:), intent(in) :: value
    
    integer(hid_t)                              :: filespace    
    integer(hid_t)                              :: memspace    
    integer(hid_t)                              :: plist_id
    integer(hid_t)                              :: dataset_id
    integer(hid_t)                              :: dataspace_id
    integer(hid_t)                              :: atype_id
    integer(hid_t)                              :: attr_id
    integer(size_t)                             :: attr_len
    integer(hsize_t), dimension(1)              :: adims = (/1/)
    integer                                     :: error
    
    call h5screate_simple_f(5, global_datasize, filespace, error)
    if ( error /= 0 ) print *,'Could not create simple: '//name
    call h5dcreate_f(group_id, name, h5t_native_double, filespace, &
                     dataset_id, error)
    call h5sclose_f(filespace, error)
    
    call h5screate_simple_f(5, local_datasize, memspace, error)
    if ( error /= 0 ) print *,'Could not create simple: '//name
    call h5dget_space_f(dataset_id, filespace, error)
    if ( error /= 0 ) print *,'Could not get space: '//name
    call h5sselect_hyperslab_f(filespace, h5s_select_set_f, slab_offset, &
           local_datasize, error)
    if ( error /= 0 ) print *,'Could not select hyperslab: '//name
    call h5pcreate_f(h5p_dataset_xfer_f, plist_id, error) 

!    select case(hdf_mpio_mode)
!      case (1)
!        call h5pset_dxpl_mpio_f(plist_id, h5fd_mpio_independent_f, error)
!      case (2)
        call h5pset_dxpl_mpio_f(plist_id, h5fd_mpio_collective_f, error)
!      case default
!        call h5pset_dxpl_mpio_f(plist_id, h5fd_mpio_independent_f, error)
!    end select

    call h5dwrite_f(dataset_id, h5t_native_double, value, local_datasize, &
           error, file_space_id=filespace, mem_space_id=memspace, &
           xfer_prp=plist_id)
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)
    if(present(desc_option))then
      attr_len = len(desc_option)
      call h5screate_simple_f(1, adims, dataspace_id, error)
      if ( error /= 0 ) print *,'Could not create simple: '//name
      call h5tcopy_f(h5t_native_character, atype_id, error)
      call h5tset_size_f(atype_id, attr_len, error)
      if ( error /= 0 ) print *,'Could not set size for: '//name
      call h5acreate_f(dataset_id, 'desc', atype_id, dataspace_id, &
             attr_id, error)
      if(hyperslab_group_master) &
        call h5awrite_f(attr_id, atype_id, desc_option, adims, error)
      call h5aclose_f(attr_id, error)
      call h5sclose_f(dataspace_id, error)
    end if
    if(present(unit_option))then
      attr_len = len(unit_option)
      call h5screate_simple_f(1, adims, dataspace_id, error)
      if ( error /= 0 ) print *,'Could not create simple: '//name
      call h5tcopy_f(h5t_native_character, atype_id, error)
      call h5tset_size_f(atype_id, attr_len, error)
      if ( error /= 0 ) print *,'Could not set size for: '//name
      call h5acreate_f(dataset_id, 'unit', atype_id, dataspace_id, &
             attr_id, error)
      if(hyperslab_group_master) &
        call h5awrite_f(attr_id, atype_id, unit_option, adims, error)
      call h5aclose_f(attr_id, error)
      call h5sclose_f(dataspace_id, error)
    end if
    call h5dclose_f(dataset_id, error)
    call h5pclose_f(plist_id, error)                     
    
  end subroutine write_ray_hyperslab_dbl_5d


  subroutine write_ray_hyperslab_int_3d(name, value, group_id, &
               global_datasize, local_datasize, slab_offset, &
               desc_option, unit_option)
    
    character(*), intent(in)                    :: name
    character(*), intent(in), optional          :: unit_option
    character(*), intent(in), optional          :: desc_option
    integer(hid_t)                              :: group_id
    integer(hsize_t), dimension(3), intent(in)  :: global_datasize
    integer(hsize_t), dimension(3), intent(in)  :: local_datasize
    integer(hsize_t), dimension(3), intent(in)  :: slab_offset
    integer, dimension(:,:,:), intent(in)       :: value
    
    integer(hid_t)                              :: filespace    
    integer(hid_t)                              :: memspace    
    integer(hid_t)                              :: plist_id
    integer(hid_t)                              :: dataset_id
    integer(hid_t)                              :: dataspace_id
    integer(hid_t)                              :: atype_id
    integer(hid_t)                              :: attr_id
    integer(size_t)                             :: attr_len
    integer(hsize_t), dimension(1)              :: adims = (/1/)
    integer                                     :: error
    
    call h5screate_simple_f(3, global_datasize, filespace, error)
    if ( error /= 0 ) print *,'Could not create simple: '//name
    call h5dcreate_f(group_id, name, h5t_native_integer, filespace, &
                     dataset_id, error)
    call h5sclose_f(filespace, error)
    
    call h5screate_simple_f(3, local_datasize, memspace, error)
    if ( error /= 0 ) print *,'Could not create simple: '//name
    call h5dget_space_f(dataset_id, filespace, error)
    if ( error /= 0 ) print *,'Could not get space: '//name
    call h5sselect_hyperslab_f(filespace, h5s_select_set_f, slab_offset, &
           local_datasize, error)
    if ( error /= 0 ) print *,'Could not select hyperslab: '//name
    call h5pcreate_f(h5p_dataset_xfer_f, plist_id, error) 

    !!select case(hdf_mpio_mode)
      !case (1)
        !call h5pset_dxpl_mpio_f(plist_id, h5fd_mpio_independent_f, error)
      !case (2)
        call h5pset_dxpl_mpio_f(plist_id, h5fd_mpio_collective_f, error)
      !case default
        !call h5pset_dxpl_mpio_f(plist_id, h5fd_mpio_independent_f, error)
    !end select

    call h5dwrite_f(dataset_id, h5t_native_integer, value, local_datasize, &
           error, file_space_id=filespace, mem_space_id=memspace, &
           xfer_prp=plist_id)
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)
    if(present(desc_option))then
      attr_len = len(desc_option)
      call h5screate_simple_f(1, adims, dataspace_id, error)
      if ( error /= 0 ) print *,'Could not create simple: '//name
      call h5tcopy_f(h5t_native_character, atype_id, error)
      call h5tset_size_f(atype_id, attr_len, error)
      if ( error /= 0 ) print *,'Could not set size for: '//name
      call h5acreate_f(dataset_id, 'desc', atype_id, dataspace_id, &
             attr_id, error)
      if(hyperslab_group_master) &
        call h5awrite_f(attr_id, atype_id, desc_option, adims, error)
      call h5aclose_f(attr_id, error)
      call h5sclose_f(dataspace_id, error)
    end if
    if(present(unit_option))then
      attr_len = len(unit_option)
      call h5screate_simple_f(1, adims, dataspace_id, error)
      if ( error /= 0 ) print *,'Could not create simple: '//name
      call h5tcopy_f(h5t_native_character, atype_id, error)
      call h5tset_size_f(atype_id, attr_len, error)
      if ( error /= 0 ) print *,'Could not set size for: '//name
      call h5acreate_f(dataset_id, 'unit', atype_id, dataspace_id, &
             attr_id, error)
      if(hyperslab_group_master) &
        call h5awrite_f(attr_id, atype_id, unit_option, adims, error)
      call h5aclose_f(attr_id, error)
      call h5sclose_f(dataspace_id, error)
    end if
    call h5dclose_f(dataset_id, error)
    call h5pclose_f(plist_id, error)                     
    
  end subroutine write_ray_hyperslab_int_3d
  
  subroutine read_1d_slab_int(name, value, group_id, datasize)
    character(*), intent(in)                    :: name
    integer(hid_t)                              :: group_id
    integer(hsize_t), dimension(1), intent(in)  :: datasize
    integer, dimension(:), intent(out)          :: value
    
    integer(hid_t)                              :: dataset_id
    integer                                     :: error
    
    call h5dopen_f(group_id, name, dataset_id, error)
    if ( error /= 0 ) print *,'Could not open: '//name
    call h5dread_f(dataset_id, h5t_native_integer, &
                    value, datasize, error)
    if ( error /= 0 ) print *,'Could not read: '//name
    call h5dclose_f(dataset_id, error)

  end subroutine read_1d_slab_int

  
  subroutine read_1d_slab_double(name, value, group_id, datasize)
    character(*), intent(in)                    :: name
    integer(hid_t)                              :: group_id
    integer(hsize_t), dimension(1), intent(in)  :: datasize
    real(dp_t), dimension(:), intent(out)       :: value
    
    integer(hid_t)                              :: dataset_id
    integer                                     :: error
    
    call h5dopen_f(group_id, name, dataset_id, error)
    if ( error /= 0 ) print *,'Could not open: '//name
    call h5dread_f(dataset_id, h5t_native_double, &
                    value, datasize, error)
    if ( error /= 0 ) print *,'Could not read: '//name
    call h5dclose_f(dataset_id, error)

  end subroutine read_1d_slab_double
  
  subroutine read_1d_slab_string(name, value, strlen, group_id, datasize)
    character(*), intent(in)                    :: name
    integer, intent(in)                         :: strlen
    integer(hid_t)                              :: group_id
    integer(hsize_t), dimension(1), intent(in)  :: datasize
    character(len=strlen), dimension(:), intent(out) :: value
    
    integer(hid_t)                              :: dataset_id
    integer(hid_t)                              :: dataspace_id
    integer(hsize_t)                            :: sizechar
    integer                                     :: error
    
    sizechar = strlen
    call h5tset_size_f(h5t_native_character, sizechar, error)
    if ( error /= 0 ) print *,'Could not set size for: '//name
    call h5dopen_f(group_id, name, dataset_id, error)
    if ( error /= 0 ) print *,'Could not open: '//name
    call h5dread_f(dataset_id, h5t_native_character, &
                    value, datasize, error)
    if ( error /= 0 ) print *,'Could not read: '//name
    call h5dclose_f(dataset_id, error)

  end subroutine read_1d_slab_string

  subroutine read_2d_slab_int(name, value, group_id, datasize)
    character(*), intent(in)                    :: name
    integer(hid_t)                              :: group_id
    integer(hsize_t), dimension(2), intent(in)  :: datasize
    integer, dimension(:,:), intent(out)        :: value
    
    integer(hid_t)                              :: dataset_id
    integer                                     :: error
    
    call h5dopen_f(group_id, name, dataset_id, error)
    if ( error /= 0 ) print *,'Could not open: '//name
    call h5dread_f(dataset_id, h5t_native_integer, &
                    value, datasize, error)
    if ( error /= 0 ) print *,'Could not read: '//name
    call h5dclose_f(dataset_id, error)

  end subroutine read_2d_slab_int
  
  subroutine read_2d_slab_double(name, value, group_id, datasize)
    character(*), intent(in)                    :: name
    integer(hid_t)                              :: group_id
    integer(hsize_t), dimension(2), intent(in)  :: datasize
    real(dp_t), dimension(:,:), intent(out)     :: value
    
    integer(hid_t)                              :: dataset_id
    integer                                     :: error
    
    call h5dopen_f(group_id, name, dataset_id, error)
    if ( error /= 0 ) print *,'Could not open: '//name
    call h5dread_f(dataset_id, h5t_native_double, &
                    value, datasize, error)
    if ( error /= 0 ) print *,'Could not read: '//name
    call h5dclose_f(dataset_id, error)

  end subroutine read_2d_slab_double


  subroutine read_ray_hyperslab_dbl_2d(name, value, group_id, &
               datasize, slab_offset)
    
    character(*), intent(in)                    :: name
    integer(hid_t)                              :: group_id
    integer(hsize_t), dimension(2), intent(in)  :: datasize
    integer(hsize_t), dimension(2), intent(in)  :: slab_offset
    real(dp_t), dimension(:,:), intent(out)     :: value
    
    integer(hid_t)                              :: filespace    
    integer(hid_t)                              :: memspace    
    integer(hid_t)                              :: dataset_id
    integer(hsize_t), dimension(2)              :: null_offset
    integer                                     :: error
    
    null_offset = 0
    call h5dopen_f(group_id, name, dataset_id, error)
    if ( error /= 0 ) print *,'Could not open: '//name
    call h5dget_space_f(dataset_id, filespace, error)
    if ( error /= 0 ) print *,'Could not get space: '//name
    call h5sselect_hyperslab_f(filespace, h5s_select_set_f, slab_offset, &
           datasize, error)
    if ( error /= 0 ) print *,'Could not select hyperslab: '//name
    call h5screate_simple_f(2, datasize, memspace, error)
    if ( error /= 0 ) print *,'Could not create simple: '//name
    call h5sselect_hyperslab_f(memspace, h5s_select_set_f, null_offset, &
           datasize, error)
    if ( error /= 0 ) print *,'Could not select hyperslab: '//name
    call h5dread_f(dataset_id, h5t_native_double, value, datasize, &
           error, file_space_id=filespace, mem_space_id=memspace)
    if ( error /= 0 ) print *,'Could not read: '//name
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)
    call h5dclose_f(dataset_id, error)
    
  end subroutine read_ray_hyperslab_dbl_2d


  subroutine read_ray_hyperslab_dbl_3d(name, value, group_id, &
               datasize, slab_offset)
    
    character(*), intent(in)                    :: name
    integer(hid_t)                              :: group_id
    integer(hsize_t), dimension(3), intent(in)  :: datasize
    integer(hsize_t), dimension(3), intent(in)  :: slab_offset
    real(dp_t), dimension(:,:,:), intent(out)   :: value
    
    integer(hid_t)                              :: filespace    
    integer(hid_t)                              :: memspace    
    integer(hid_t)                              :: dataset_id
    integer(hsize_t), dimension(3)              :: null_offset
    integer                                     :: error
    
    null_offset = 0
    call h5dopen_f(group_id, name, dataset_id, error)
    if ( error /= 0 ) print *,'Could not open: '//name
    call h5dget_space_f(dataset_id, filespace, error)
    if ( error /= 0 ) print *,'Could not get space: '//name
    call h5sselect_hyperslab_f(filespace, h5s_select_set_f, slab_offset, &
           datasize, error)
    if ( error /= 0 ) print *,'Could not select hyperslab: '//name
    call h5screate_simple_f(3, datasize, memspace, error)
    if ( error /= 0 ) print *,'Could not create simple: '//name
    call h5sselect_hyperslab_f(memspace, h5s_select_set_f, null_offset, &
           datasize, error)
    if ( error /= 0 ) print *,'Could not select hyperslab: '//name
    call h5dread_f(dataset_id, h5t_native_double, value, datasize, &
           error, file_space_id=filespace, mem_space_id=memspace)
    if ( error /= 0 ) print *,'Could not read: '//name
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)
    call h5dclose_f(dataset_id, error)
    
  end subroutine read_ray_hyperslab_dbl_3d


  subroutine read_ray_hyperslab_dbl_4d(name, value, group_id, &
               datasize, slab_offset)
    
    character(*), intent(in)                    :: name
    integer(hid_t)                              :: group_id
    integer(hsize_t), dimension(4), intent(in)  :: datasize
    integer(hsize_t), dimension(4), intent(in)  :: slab_offset
    real(dp_t), dimension(:,:,:,:), intent(out) :: value
    
    integer(hid_t)                              :: filespace    
    integer(hid_t)                              :: memspace    
    integer(hid_t)                              :: dataset_id
    integer(hsize_t), dimension(4)              :: null_offset
    integer                                     :: error
    
    null_offset = 0
    call h5dopen_f(group_id, name, dataset_id, error)
    if ( error /= 0 ) print *,'Could not open: '//name
    call h5dget_space_f(dataset_id, filespace, error)
    if ( error /= 0 ) print *,'Could not get space: '//name
    call h5sselect_hyperslab_f(filespace, h5s_select_set_f, slab_offset, &
           datasize, error)
    if ( error /= 0 ) print *,'Could not select hyperslab: '//name
    call h5screate_simple_f(4, datasize, memspace, error)
    if ( error /= 0 ) print *,'Could not create simple: '//name
    call h5sselect_hyperslab_f(memspace, h5s_select_set_f, null_offset, &
           datasize, error)
    if ( error /= 0 ) print *,'Could not select hyperslab: '//name
    call h5dread_f(dataset_id, h5t_native_double, value, datasize, &
           error, file_space_id=filespace, mem_space_id=memspace)
    if ( error /= 0 ) print *,'Could not read: '//name
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)
    call h5dclose_f(dataset_id, error)
    
  end subroutine read_ray_hyperslab_dbl_4d


  subroutine read_ray_hyperslab_dbl_5d(name, value, group_id, &
               datasize, slab_offset)
    
    character(*), intent(in)                    :: name
    integer(hid_t)                              :: group_id
    integer(hsize_t), dimension(5), intent(in)  :: datasize
    integer(hsize_t), dimension(5), intent(in)  :: slab_offset
    real(dp_t), dimension(:,:,:,:,:), intent(out) :: value
    
    integer(hid_t)                              :: filespace    
    integer(hid_t)                              :: memspace    
    integer(hid_t)                              :: dataset_id
    integer(hsize_t), dimension(5)              :: null_offset
    integer                                     :: error
    
    null_offset = 0
    call h5dopen_f(group_id, name, dataset_id, error)
    if ( error /= 0 ) print *,'Could not open: '//name
    call h5dget_space_f(dataset_id, filespace, error)
    if ( error /= 0 ) print *,'Could not get space: '//name
    call h5sselect_hyperslab_f(filespace, h5s_select_set_f, slab_offset, &
           datasize, error)
    if ( error /= 0 ) print *,'Could not select hyperslab: '//name
    call h5screate_simple_f(5, datasize, memspace, error)
    if ( error /= 0 ) print *,'Could not create simple: '//name
    call h5sselect_hyperslab_f(memspace, h5s_select_set_f, null_offset, &
           datasize, error)
    if ( error /= 0 ) print *,'Could not select hyperslab: '//name
    call h5dread_f(dataset_id, h5t_native_double, value, datasize, &
           error, file_space_id=filespace, mem_space_id=memspace)
    if ( error /= 0 ) print *,'Could not read: '//name
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)
    call h5dclose_f(dataset_id, error)
    
  end subroutine read_ray_hyperslab_dbl_5d


  subroutine read_ray_hyperslab_int_3d(name, value, group_id, &
               datasize, slab_offset)
    
    character(*), intent(in)                    :: name
    integer(hid_t)                              :: group_id
    integer(hsize_t), dimension(3), intent(in)  :: datasize
    integer(hsize_t), dimension(3), intent(in)  :: slab_offset
    integer, dimension(:,:,:), intent(out)       :: value
    
    integer(hid_t)                              :: filespace    
    integer(hid_t)                              :: memspace    
    integer(hid_t)                              :: dataset_id
    integer(hsize_t), dimension(3)              :: null_offset
    integer                                     :: error
    
    null_offset = 0
    call h5dopen_f(group_id, name, dataset_id, error)
    if ( error /= 0 ) print *,'Could not open: '//name
    call h5dget_space_f(dataset_id, filespace, error)
    if ( error /= 0 ) print *,'Could not get space: '//name
    call h5sselect_hyperslab_f(filespace, h5s_select_set_f, slab_offset, &
           datasize, error)
    if ( error /= 0 ) print *,'Could not select hyperslab: '//name
    call h5screate_simple_f(3, datasize, memspace, error)
    if ( error /= 0 ) print *,'Could not create simple: '//name
    call h5sselect_hyperslab_f(memspace, h5s_select_set_f, null_offset, &
           datasize, error)
    if ( error /= 0 ) print *,'Could not select hyperslab: '//name
    call h5dread_f(dataset_id, h5t_native_integer, value, datasize, &
           error, file_space_id=filespace, mem_space_id=memspace)
    if ( error /= 0 ) print *,'Could not read: '//name
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)
    call h5dclose_f(dataset_id, error)
    
    
  end subroutine read_ray_hyperslab_int_3d
  

end module hdf5_read_write
