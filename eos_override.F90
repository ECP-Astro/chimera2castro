module eos_override_module

  implicit none

  public eos_override

contains

  ! This is a user hook to override the details of the EOS state.

  subroutine eos_override(state)

    !$acc routine seq

    use network, only: naux, nspec, nspec_evolve, aion, zion
    use bl_constants_module, only: one
    use eos_type_module, only: eos_t, composition

    implicit none

    type (eos_t) :: state

    if ( nspec > nspec_evolve ) then
        if ( naux == 2 ) then
            aion(nspec) = state % aux(1)
            zion(nspec) = state % aux(2)
        else if ( naux == 3 ) then
            aion(nspec) = state % aux(2)
            zion(nspec) = state % aux(3)
        end if
        call composition(state)
    end if

  end subroutine eos_override

end module eos_override_module
