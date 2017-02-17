module model_interp_module

  use bl_types
  use bl_constants_module

  implicit none

  integer, parameter :: kx = 3, ky = 3, iknot = 0
  integer, parameter :: idx = 0, idy = 0

  interface interp1d_linear
    module procedure interp1d_linear_vec
    module procedure interp1d_linear_scalar
  end interface interp1d_linear

contains

  subroutine interp1d_linear_scalar( ix, ix_max, x_in, f_in, x_out, f_out )

    ! input variables
    integer, intent(in) :: ix, ix_max
    real (dp_t), intent(in) :: x_in(:)
    real (dp_t), intent(in) :: f_in(:)
    real (dp_t), intent(in) :: x_out

    ! output variables
    real (dp_t), intent(out) :: f_out

    ! local variables
    real (dp_t) :: dx0, dx1, dx
    integer :: j

    if ( ix < ix_max ) then
      dx  = x_in(ix+1) - x_in(ix)
      dx1 = (x_out - x_in(ix))/dx
      dx0 = one - dx1
      f_out = dx0*f_in(ix) + dx1*f_in(ix+1)
    else if ( ix == 0 ) then
      f_out = f_in(1)
    else if ( ix == ix_max ) then
      f_out = f_in(ix_max)
    else
      f_out = zero
    end if

    return
  end subroutine interp1d_linear_scalar

  subroutine interp1d_linear_vec( ix, ix_max, x_in, f_in, x_out, f_out )

    ! input variables
    integer, intent(in) :: ix(:), ix_max
    real (dp_t), intent(in) :: x_in(:)
    real (dp_t), intent(in) :: f_in(:)
    real (dp_t), intent(in) :: x_out(:)

    ! output variables
    real (dp_t), intent(out) :: f_out(size(x_out))

    ! local variables
    real (dp_t) :: dx0, dx1, dx
    integer :: j

    do j = 1, size(x_out)
      if ( ix(j) > 0 .and. ix(j) < ix_max ) then
        dx  = x_in(ix(j)+1) - x_in(ix(j))
        dx1 = (x_out(j) - x_in(ix(j)))/dx
        dx0 = one - dx1
        f_out(j) = dx0*f_in(ix(j)) + dx1*f_in(ix(j)+1)
      else if ( ix(j) == 0 ) then
        f_out(j) = f_in(1)
      else if ( ix(j) == ix_max ) then
        f_out(j) = f_in(ix_max)
      else
        f_out(j) = zero
      end if
    end do

    return
  end subroutine interp1d_linear_vec

  subroutine interp1d_spline( x_in, f_in, x_out, f_out )

    ! input variables
    real (dp_t), intent(in) :: x_in(:)
    real (dp_t), intent(in) :: f_in(:)
    real (dp_t), intent(in) :: x_out(:)

    ! output variables
    real (dp_t), intent(out) :: f_out(size(x_out))

    ! local variables
    real (dp_t) :: fpp(size(x_in)), fval, fpval, fppval
    integer :: j, nx

    nx = size(x_in)
    fpp = zero

    call spline_cubic_set( nx, x_in, f_in, 3, zero, 3, zero, fpp )
    DO j = 1, size(x_out)
      call spline_cubic_val( nx, x_in, f_in, fpp, x_out(j), fval, fpval, fppval )
      f_out(j) = fval
    END DO

    return
  end subroutine interp1d_spline

  subroutine interp2d_linear( ix, ix_max, iy, iy_max, x_in, y_in, f_in, x_out, y_out, f_out )

    ! input variables
    integer, intent(in) :: ix(:,:), ix_max
    integer, intent(in) :: iy(:,:), iy_max
    real (dp_t), intent(in) :: x_in(:)
    real (dp_t), intent(in) :: y_in(:)
    real (dp_t), intent(in) :: f_in(:,:)
    real (dp_t), intent(in) :: x_out(:,:)
    real (dp_t), intent(in) :: y_out(:,:)

    ! output variables
    real (dp_t), intent(out) :: f_out(size(x_out,1),size(x_out,2))

    ! local variables
    real (dp_t) :: dx0, dx1, dx
    real (dp_t) :: dy0, dy1, dy
    real (dp_t) :: fxy0, fxy1
    integer :: i, j

    do j = 1, size(x_out,2)
      do i = 1, size(x_out,1)
        
        if ( ix(i,j) > 0 .and. ix(i,j) < ix_max ) then

          dx  = x_in(ix(i,j)+1) - x_in(ix(i,j))
          dx1 = (x_out(i,j) - x_in(ix(i,j)))/dx
          dx0 = one - dx1

          if ( iy(i,j) > 0 .and. iy(i,j) < iy_max ) then

            dy  = y_in(iy(i,j)+1) - y_in(iy(i,j))
            dy1 = (y_out(i,j) - y_in(iy(i,j)))/dy
            dy0 = one - dy1

            fxy0 = dx0*f_in(ix(i,j),iy(i,j)) + dx1*f_in(ix(i,j)+1,iy(i,j))
            fxy1 = dx0*f_in(ix(i,j),iy(i,j)+1) + dx1*f_in(ix(i,j)+1,iy(i,j)+1)

            f_out(i,j) = dy0*fxy0 + dy1*fxy1

          else if ( iy(i,j) == 0 ) then
            f_out(i,j) = dx0*f_in(ix(i,j),1) + dx1*f_in(ix(i,j)+1,1)
          else if ( iy(i,j) == iy_max ) then
            f_out(i,j) = dx0*f_in(ix(i,j),iy_max) + dx1*f_in(ix(i,j)+1,iy_max)
          else
            f_out(i,j) = zero
          end if

        else if ( ix(i,j) == 0 ) then

          if ( iy(i,j) > 0 .and. iy(i,j) < iy_max ) then

            dy  = y_in(iy(i,j)+1) - y_in(iy(i,j))
            dy1 = (y_out(i,j) - y_in(iy(i,j)))/dy
            dy0 = one - dy1

            f_out(i,j) = dy0*f_in(1,iy(i,j)) + dy1*f_in(1,iy(i,j)+1)

          else if ( iy(i,j) == 0 ) then
            f_out(i,j) = f_in(1,1)
          else if ( iy(i,j) == iy_max ) then
            f_out(i,j) = f_in(1,iy_max)
          else
            f_out(i,j) = zero
          end if

        else if ( ix(i,j) == ix_max ) then

          if ( iy(i,j) > 0 .and. iy(i,j) < iy_max ) then

            dy  = y_in(iy(i,j)+1) - y_in(iy(i,j))
            dy1 = (y_out(i,j) - y_in(iy(i,j)))/dy
            dy0 = one - dy1

            f_out(i,j) = dy0*f_in(ix_max,iy(i,j)) + dy1*f_in(ix_max,iy(i,j)+1)

          else if ( iy(i,j) == 0 ) then
            f_out(i,j) = f_in(ix_max,1)
          else if ( iy(i,j) == iy_max ) then
            f_out(i,j) = f_in(ix_max,iy_max)
          else
            f_out(i,j) = zero
          end if

        else
          f_out(i,j) = zero
        end if
      end do
    end do

    return
  end subroutine interp2d_linear

  subroutine interp2d_spline( x_in, y_in, f_in, x_out, y_out, f_out )

    use bl_error_module
    use bspline_module, only: db2ink, db2val

    ! input variables
    real (dp_t), intent(in) :: x_in(:)
    real (dp_t), intent(in) :: y_in(:)
    real (dp_t), intent(in) :: f_in(:,:)
    real (dp_t), intent(in) :: x_out(:,:)
    real (dp_t), intent(in) :: y_out(:,:)

    ! output variables
    real (dp_t), intent(out) :: f_out(size(x_out,1),size(x_out,2))

    ! local variables
    real (dp_t) :: tx(size(x_in)+kx)
    real (dp_t) :: ty(size(y_in)+ky)
    real (dp_t) :: bcoef(size(x_in),size(y_in))

    integer :: inbvx, inbvy, iloy
    integer :: nx, ny
    real (dp_t) :: x, y, fval
    integer :: i, j, iflag

    nx = size( x_in )
    ny = size( y_in )

    inbvx = 1
    inbvy = 1
    iloy = 1

    call db2ink( x_in, nx, y_in, ny, f_in, kx, ky, iknot, tx, ty, bcoef, iflag )
    if ( iflag /= 0 ) then
      write(*,*) "error in db2ink: ", iflag
      call bl_error("error in db2ink")
    end if

    do j = 1, size(x_out,2)
      do i = 1, size(x_out,1)
        x = min( max( x_out(i,j), x_in(1) ), x_in(nx) )
        y = min( max( y_out(i,j), y_in(1) ), y_in(ny) )
        call db2val( x, y, idx, idy, tx, ty, nx, ny, kx, ky, bcoef, fval, iflag, inbvx, inbvy, iloy )
        if ( iflag /= 0 ) then
          write(*,*) "error in db2val: ", x, y, iflag
          call bl_error("error in db2val")
        end if
        f_out(i,j) = fval
      end do
    end do

    return
  end subroutine interp2d_spline

end module model_interp_module
