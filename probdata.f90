module probdata_module

  use bl_types

  ! chimera hdf5 filename
  character (256), save :: chimera_fname
  ! mesa model filename
  character (256), save :: mesa_fname

  ! flag for interpolation method (1=linear, 2=cubic spline)
  integer, save :: interp_method

  ! flag for eos_input (1=rt,2=rh,3=tp,4=rp,5=re,6=ps,7=ph,8=th)
  integer, save :: eos_input

  ! max radius to use chimera data
  ! if r > max_radius, use mesa data
  real (dp_t), save :: max_radius

  ! inner boundary for chimera data
  ! if r < radius_inner, use rho_inner, temp_inner
  real (dp_t), save :: radius_inner
  real (dp_t), save :: rho_inner
  real (dp_t), save :: temp_inner
  integer, save :: i_inner

  logical, save :: do_particles

  ! whether to use quadrature for interpolation
  logical, save :: use_quad
  ! number of quadrature points if using quadrature
  integer, save :: nquad

end module probdata_module
