module problem_tagging_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  public

contains

  ! This routine prevents zones bordering center from being refined
  
  subroutine set_problem_tags(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                              state,state_l1,state_l2,state_l3,state_h1,state_h2,state_h3,&
                              set,clear,&
                              lo,hi,&
                              dx,problo,time,level) &
                              bind(C, name="set_problem_tags")

    use meth_params_module, only: NVAR

    use amrex_fort_module, only : rt => amrex_real
    use bl_constants_module, only: HALF
    use parallel, only: parallel_IOProcessor
    use prob_params_module, only: center, dim, n_error_buf, ref_ratio, blocking_factor, dx_level
    use probdata_module, only: radius_inner
    implicit none

    integer         ,intent(in   ) :: lo(3),hi(3)
    integer         ,intent(in   ) :: state_l1,state_l2,state_l3, &
                                      state_h1,state_h2,state_h3
    integer         ,intent(in   ) :: tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
    real(rt)        ,intent(in   ) :: state(state_l1:state_h1, &
                                      state_l2:state_h2, &
                                      state_l3:state_h3,NVAR)
    integer         ,intent(inout) :: tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
    real(rt)        ,intent(in   ) :: problo(3),dx(3),time
    integer         ,intent(in   ) :: level,set,clear

    integer                        :: i, j, k, n
    real(rt)                       :: loc(3), r

    logical                        :: center_test(3)
    integer                        :: boundary_buf(3), idx(3), icenter(3)

    ! We must ensure that the center zones are untagged.
    ! To do this properly we need to be aware of BoxLib's strategy
    ! for tagging, which is not cell-based, but rather chunk-based.
    ! The size of the chunk on the coarse grid is given by
    ! blocking_factor / ref_ratio -- the idea here being that
    ! blocking_factor is the smallest possible group of cells on a
    ! given level, so the smallest chunk of cells possible on the
    ! coarse grid is given by that number divided by the refinement ratio.
    ! So we cannot tag anything within that distance from the boundary.
    ! Additionally we need to stay a further amount n_error_buf away,
    ! since n_error_buf zones are always added as padding around
    ! tagged zones.

!   icenter(1:dim) = NINT((center(1:dim)-problo(1:dim)) / dx(1:dim))
!   if ( level == 0 ) then
!     boundary_buf(1:dim) = n_error_buf(level) + 1
!   else
!     boundary_buf(1:dim) = n_error_buf(level) + product(ref_ratio(1:dim,0:level),dim=2) - 1
!   end if

!   boundary_buf(1:dim) = n_error_buf(level) + max( blocking_factor(level+1) / ref_ratio(1:dim, level), 1 )
!   do n = level, 1, -1
!     boundary_buf(1:dim) = boundary_buf(1:dim) * ref_ratio(1:dim,n) + n_error_buf(n-1) + blocking_factor(n) / ref_ratio(1:dim,n-1)
!   end do
!   boundary_buf(1:dim) = n_error_buf(0) + blocking_factor(1) / ref_ratio(1:dim,0)
!   do n = 1, level
!     boundary_buf(1:dim) = boundary_buf(1:dim) * ref_ratio(1:dim,n-1) + n_error_buf(n) + blocking_factor(n+1) / ref_ratio(1:dim,n)
!     boundary_buf(1:dim) = boundary_buf(1:dim) * ref_ratio(1:dim,n-1) + 1
!   end do

!   do k = lo(3), hi(3)
!     loc(3) = problo(3) + (dble(k) + HALF) * dx(3) - center(3)
!     do j = lo(2), hi(2)
!       loc(2) = problo(2) + (dble(j) + HALF) * dx(2) - center(2)
!       do i = lo(1), hi(1)
!         loc(1) = problo(1) + (dble(i) + HALF) * dx(1) - center(1)

!         r = ( sum(loc**2) )**HALF
!         idx = [i, j, k]
!         center_test = .false.
!         center_test(1:dim) = abs(loc(1:dim)) <= dble(boundary_buf(1:dim)+1)*dx(1:dim)
!         center_test(1:dim) = abs(loc(1:dim)) <= dble(boundary_buf(1:dim)+1)*dx_level(1:dim,0)
!         center_test(1:dim) = loc(1:dim) >= center(1:dim) - (dble(boundary_buf(1:dim))+HALF)*dx(1:dim) .and. &
!                              loc(1:dim) <= center(1:dim) + (dble(boundary_buf(1:dim))+HALF)*dx(1:dim)
!         center_test(1:dim) = idx(1:dim) >= icenter(1:dim) - 2 - boundary_buf(1:dim) .and. &
!                              idx(1:dim) <= icenter(1:dim)     + boundary_buf(1:dim)

!         if ( r <= maxval( (dble(boundary_buf(1:dim))+1) * dx(1:dim) ) ) then
!         if ( all(center_test(1:dim)) ) then
!           tag(i,j,k) = clear
!           if ( all(idx==icenter) .or. all(idx==icenter-1) ) then
!             !$omp critical
!             write(*,'(a,1i23)')     'level          =', level
!             write(*,'(a,3i23)')     '  ref_ratio    =', ref_ratio(1:dim,level)
!             write(*,'(a,3i23)')     '  icenter      =', icenter(1:dim)
!             write(*,'(a,3i23)')     '  boundary_buf =', boundary_buf(1:dim)
!             write(*,'(a,3es23.15)') '  dx           =', dx(1:dim)
!             write(*,'(a,3es23.15)') '  loc          =', loc(1:dim)
!             write(*,'(a,3i23)')     '  idx          =', idx(1:dim)
!             !$omp end critical
!           end if
!         end if

!       end do
!     end do
!   end do

  end subroutine set_problem_tags

end module problem_tagging_module
