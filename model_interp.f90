module model_interp_module

  use bl_fort_module, only: rt => c_real
  use bl_constants_module

  implicit none

  integer, parameter :: kx = 3, ky = 3, kz = 3, iknot = 0
  integer, parameter :: idx = 0, idy = 0, idz = 0

  interface interp1d_linear
    module procedure interp1d_linear_vec
    module procedure interp1d_linear_scalar
  end interface interp1d_linear

contains

  subroutine interp1d_linear_scalar( ix, ix_max, x_in, f_in, x_out, f_out )

    ! input variables
    integer, intent(in) :: ix, ix_max
    real (rt), intent(in) :: x_in(:)
    real (rt), intent(in) :: f_in(:)
    real (rt), intent(in) :: x_out

    ! output variables
    real (rt), intent(out) :: f_out

    ! local variables
    real (rt) :: dx0, dx1, dx
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
    real (rt), intent(in) :: x_in(:)
    real (rt), intent(in) :: f_in(:)
    real (rt), intent(in) :: x_out(:)

    ! output variables
    real (rt), intent(out) :: f_out(size(x_out))

    ! local variables
    real (rt) :: dx0, dx1, dx
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
    real (rt), intent(in) :: x_in(:)
    real (rt), intent(in) :: f_in(:)
    real (rt), intent(in) :: x_out(:)

    ! output variables
    real (rt), intent(out) :: f_out(size(x_out))

    ! local variables
    real (rt) :: fpp(size(x_in)), fval, fpval, fppval
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
    real (rt), intent(in) :: x_in(:)
    real (rt), intent(in) :: y_in(:)
    real (rt), intent(in) :: f_in(:,:)
    real (rt), intent(in) :: x_out(:,:)
    real (rt), intent(in) :: y_out(:,:)

    ! output variables
    real (rt), intent(out) :: f_out(size(x_out,1),size(x_out,2))

    ! local variables
    real (rt) :: dx0, dx1, dx
    real (rt) :: dy0, dy1, dy
    real (rt) :: fxy0, fxy1
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
    real (rt), intent(in) :: x_in(:)
    real (rt), intent(in) :: y_in(:)
    real (rt), intent(in) :: f_in(:,:)
    real (rt), intent(in) :: x_out(:,:)
    real (rt), intent(in) :: y_out(:,:)

    ! output variables
    real (rt), intent(out) :: f_out(size(x_out,1),size(x_out,2))

    ! local variables
    real (rt) :: tx(size(x_in)+kx)
    real (rt) :: ty(size(y_in)+ky)
    real (rt) :: bcoef(size(x_in),size(y_in))

    integer :: inbvx, inbvy, iloy
    integer :: nx, ny
    real (rt) :: x, y, fval
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

  subroutine interp3d_linear( ix, ix_max, iy, iy_max, iz, iz_max, &
  &                           x_in, y_in, z_in, f_in, x_out, y_out, z_out, f_out )

    ! input variables
    integer, intent(in) :: ix(:,:,:), ix_max
    integer, intent(in) :: iy(:,:,:), iy_max
    integer, intent(in) :: iz(:,:,:), iz_max
    real (rt), intent(in) :: x_in(:)
    real (rt), intent(in) :: y_in(:)
    real (rt), intent(in) :: z_in(:)
    real (rt), intent(in) :: f_in(:,:,:)
    real (rt), intent(in) :: x_out(:,:,:)
    real (rt), intent(in) :: y_out(:,:,:)
    real (rt), intent(in) :: z_out(:,:,:)

    ! output variables
    real (rt), intent(out) :: f_out(size(x_out,1),size(x_out,2),size(x_out,3))

    ! local variables
    integer :: ix0, iy0, iz0
    integer :: ix1, iy1, iz1
    real (rt) :: dx, dy, dz
    real (rt) :: c(8), d(8)
    integer :: i, j, k

    do k = 1, size(x_out,3)
      do j = 1, size(x_out,2)
        do i = 1, size(x_out,1)
          ix0 = max( ix(i,j,k)  , 1 )
          ix1 = min( ix(i,j,k)+1, ix_max )
          if ( ix1 > ix0 ) then
            dx = ( x_out(i,j,k) - x_in(ix0) ) / ( x_in(ix1) - x_in(ix0) )
          else
            dx = zero
          end if

          iy0 = max( iy(i,j,k)  , 1 )
          iy1 = min( iy(i,j,k)+1, iy_max )
          if ( iy1 > iy0 ) then
            dy = ( y_out(i,j,k) - y_in(iy0) ) / ( y_in(iy1) - y_in(iy0) )
          else
            dy = zero
          end if

          iz0 = max( iz(i,j,k)  , 1 )
          iz1 = min( iz(i,j,k)+1, iz_max )
          if ( iz1 > iz0 ) then
            dz = ( z_out(i,j,k) - z_in(iz0) ) / ( z_in(iz1) - z_in(iz0) )
          else
            dz = zero
          end if

          d(1) = one
          d(2) = dx
          d(3) = dy
          d(4) = dz
          d(5) = dx*dy
          d(6) = dx*dz
          d(7) = dy*dz
          d(8) = dx*dy*dz

          c(1) = f_in(ix0,iy0,iz0)
          c(2) = f_in(ix1,iy0,iz0) - f_in(ix0,iy0,iz0)
          c(3) = f_in(ix0,iy1,iz0) - f_in(ix0,iy0,iz0)
          c(4) = f_in(ix0,iy0,iz1) - f_in(ix0,iy0,iz0)
          c(5) = f_in(ix1,iy1,iz0) - f_in(ix0,iy1,iz0) + &
                 f_in(ix0,iy0,iz0) - f_in(ix1,iy0,iz0)
          c(6) = f_in(ix1,iy0,iz1) - f_in(ix0,iy0,iz1) + &
                 f_in(ix0,iy0,iz0) - f_in(ix1,iy0,iz0)
          c(7) = f_in(ix0,iy1,iz1) - f_in(ix0,iy0,iz1) + &
                 f_in(ix0,iy0,iz0) - f_in(ix0,iy1,iz0)
          c(8) = f_in(ix1,iy1,iz1) - f_in(ix0,iy1,iz1) + &
                 f_in(ix0,iy0,iz1) - f_in(ix1,iy0,iz1) + &
                 f_in(ix0,iy1,iz0) - f_in(ix1,iy1,iz0) + &
                 f_in(ix1,iy0,iz0) - f_in(ix0,iy0,iz0)

          f_out(i,j,k) = sum( d*c )

        end do
      end do
    end do

    return
  end subroutine interp3d_linear

  subroutine interp3d_spline( x_in, y_in, z_in, f_in, x_out, y_out, z_out, f_out )

    use bl_error_module
    use bspline_module, only: db3ink, db3val

    ! input variables
    real (rt), intent(in) :: x_in(:)
    real (rt), intent(in) :: y_in(:)
    real (rt), intent(in) :: z_in(:)
    real (rt), intent(in) :: f_in(:,:,:)
    real (rt), intent(in) :: x_out(:,:,:)
    real (rt), intent(in) :: y_out(:,:,:)
    real (rt), intent(in) :: z_out(:,:,:)

    ! output variables
    real (rt), intent(out) :: f_out(size(x_out,1),size(x_out,2),size(x_out,3))

    ! local variables
    real (rt) :: tx(size(x_in)+kx)
    real (rt) :: ty(size(y_in)+ky)
    real (rt) :: tz(size(z_in)+kz)
    real (rt) :: bcoef(size(x_in),size(y_in),size(z_in))

    integer :: inbvx, inbvy, inbvz, iloy, iloz
    integer :: nx, ny, nz
    real (rt) :: x, y, z, fval
    integer :: i, j, k, iflag

    nx = size( x_in )
    ny = size( y_in )
    nz = size( z_in )

    inbvx = 1
    inbvy = 1
    inbvz = 1
    iloy = 1
    iloz = 1

    call db3ink( x_in, nx, y_in, ny, z_in, nz, f_in, kx, ky, kz, iknot, tx, ty, tz, bcoef, iflag )
    if ( iflag /= 0 ) then
      write(*,*) "error in db3ink: ", iflag
      call bl_error("error in db3ink")
    end if

    do k = 1, size(x_out,3)
      do j = 1, size(x_out,2)
        do i = 1, size(x_out,1)
          x = min( max( x_out(i,j,k), x_in(1) ), x_in(nx) )
          y = min( max( y_out(i,j,k), y_in(1) ), y_in(ny) )
          z = min( max( z_out(i,j,k), z_in(1) ), z_in(nz) )
          call db3val( x, y, z, idx, idy, idz, tx, ty, tz, nx, ny, nz, kx, ky, kz, bcoef, &
          &            fval, iflag, inbvx, inbvy, inbvz, iloy, iloz )
          if ( iflag /= 0 ) then
            write(*,*) "error in db3val: ", x, y, z, iflag
            call bl_error("error in db3val")
          end if
          f_out(i,j,k) = fval
        end do
      end do
    end do

    return
  end subroutine interp3d_spline

end module model_interp_module
