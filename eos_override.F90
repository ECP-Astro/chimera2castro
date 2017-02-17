module eos_override_module

  implicit none

  public eos_override

contains

  ! This is a user hook to override the details of the EOS state.

  subroutine eos_override(state)

    !$acc routine seq

    use network, only: naux, nspec, nspec_evolve, aion, zion, aion_inv
    use bl_constants_module, only: one
    use eos_type_module, only: eos_t, composition, minye, maxye

    implicit none

    type (eos_t) :: state

    if ( naux == 1 ) then
      state % y_e = state % aux(1)
      state % mu_e = one / state % y_e
      state % zbar = state % abar / state % mu_e
    end if

    if ( nspec > nspec_evolve ) then
      if ( naux == 2 ) then
        aion(nspec) = state % aux(1)
        zion(nspec) = state % aux(2)
      else if ( naux == 3 ) then
        aion(nspec) = state % aux(2)
        zion(nspec) = state % aux(3)
      end if
      aion_inv(nspec) = one / aion(nspec)
      call composition(state)
    end if

    if ( state % y_e > maxye ) then
      state % y_e = maxye
      state % mu_e = one / state % y_e
      state % zbar = state % abar / state % mu_e
    else if ( state % y_e < minye ) then
      state % y_e = minye
      state % mu_e = one / state % y_e
      state % zbar = state % abar / state % mu_e
    end if

  end subroutine eos_override

end module eos_override_module
