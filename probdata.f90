module probdata_module

  use amrex_fort_module, only: rt => amrex_real

  ! chimera hdf5 filename
  character (256), save :: chimera_fname
  ! star model filename
  character (256), save :: star_fname

  ! flag for interpolation method (1=linear, 2=cubic spline)
  integer, save :: interp_method

  ! flag for eos_input (1=rt,2=rh,3=tp,4=rp,5=re,6=ps,7=ph,8=th)
  integer, save :: eos_input

  ! max radius to use chimera data
  ! if r > max_radius, use star data
  real (rt), save :: max_radius

  ! inner boundary for chimera data
  ! if r < radius_inner, use rho_inner, temp_inner
  real (rt), save :: radius_inner
  real (rt), save :: rho_inner
  real (rt), save :: temp_inner
  integer, save :: i_inner

  logical, save :: do_particles

  ! whether to use quadrature for interpolation
  logical, save :: use_quad
  ! number of quadrature points if using quadrature
  integer, save :: nquad

end module probdata_module
