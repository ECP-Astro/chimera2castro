module probdata_module

  use bl_types

  character (len=80), save :: model_name
  character (len=80), save :: mesa_name

  ! flag for interpolation method (1=linear, 2=cubic spline)
  integer, save :: model_interp_method

  ! flag for eos_input (1=rt,2=rh,3=tp,4=rp,5=re,6=ps,7=ph,8=th)
  integer, save :: model_eos_input

  real (dp_t), save :: min_radius
  real (dp_t), save :: max_radius

  logical, save :: do_particles

end module probdata_module
