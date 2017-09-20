module star_parser_module

  use parallel, only: parallel_IOProcessor
  use amrex_fort_module, only: rt => amrex_real
  use probdata_module, only: interp_method, star_type

  implicit none

  integer, save :: nx_star, nnc_star
  integer, save :: imin_star, imax_star

  real (rt), allocatable, save :: x_e_star(:)
  real (rt), allocatable, save :: dx_e_star(:)
  real (rt), allocatable, save :: volx_e_star(:)
  real (rt), allocatable, save :: volx_c_star(:)
  real (rt), allocatable, save :: dvolx_e_star(:)
  real (rt), allocatable, save :: dvolx_c_star(:)
  real (rt), allocatable, save :: dvol_e_star(:)
  real (rt), allocatable, save :: dmass_e_star(:)

  real (rt), allocatable, save :: u_c_star(:)
  real (rt), allocatable, save :: rho_c_star(:)
  real (rt), allocatable, save :: t_c_star(:)
  real (rt), allocatable, save :: ye_c_star(:)

  real (rt), allocatable, save :: xn_c_star(:,:)
  real (rt), allocatable, save :: a_nuc_star(:)
  real (rt), allocatable, save :: z_nuc_star(:)
  real (rt), allocatable, save :: be_nuc_star(:)
  character(len=5), allocatable, save :: name_nuc_star(:)

  real (rt), allocatable, save :: a_aux_c_star(:)
  real (rt), allocatable, save :: z_aux_c_star(:)
  
  contains

  subroutine read_star_file( filename )

    use actual_network, only: nspec, nspec_evolve, aion, zion
    use bl_constants_module

    ! input variables
    character(len=*), intent(in) :: filename

    ! local variables
    character(len=1600)      :: line           ! read line for structure data

    integer                  :: zone_read      ! read in zone number
    real (rt)              :: rho_read     ! read in ln(density) [g cm^{-3}]
    real (rt)              :: t_read       ! read in ln(temperature) [k]
    real (rt)              :: r_read       ! read in ln(radius) [cm]
    real (rt)              :: u_read         ! read in velocity [cm s^{-1}]
    real (rt), allocatable :: xn_read(:)     ! read in progenitor composition array

    real (rt) :: dum1    ! dummy variable
    real (rt) :: dum2    ! dummy variable
    real (rt) :: dum3    ! dummy variable
    real (rt) :: dum4    ! dummy variable
    real (rt) :: dum5    ! dummy variable
    real (rt) :: dum6    ! dummy variable
    real (rt) :: dum7    ! dummy variable
    real (rt) :: dum8    ! dummy variable

    real (rt) :: dr, dvol

    integer :: itmp(nspec)
    integer :: net_to_castro(nspec)
    integer, allocatable :: net_in_castro(:)

    integer :: i, j, k   ! loop indices
    integer :: nread        ! file unit number
    integer :: istate       ! iostat variable

    character(len=32) :: header_name
    character(len=4)  :: nnc_string
    character(len=22) :: header_format
    character(len=44) :: nonnse_format

    ! formats

    101 format(59x,i5)
    102 format(46x,i4)


    ! open input data file
    open( newunit=nread, file=trim(adjustl(filename)), status='old', iostat=istate )
    if ( istate /= 0 ) then
      call bl_error('Aborting now -- please supply star file')
    end if

    ! read progenitor file
    do
      read(nread,'(a)',iostat=istate) line
      if ( istate /= 0 ) exit

      header_name = line(1:32)
	      if ( trim(adjustl(header_name)) == 'n_shells' ) then
	        read(line,101) nx_star
	        cycle
	      end if

	      if ( trim(adjustl(header_name)) == 'species' ) then
	        read(line,101) nnc_star
	        exit
	      end if


    end do
    read(nread,*)

    write(nnc_string,'(i4)') nnc_star

    if(star_type==1) then
    	nonnse_format = '(i5,1x,7(3x,es23.16,1x),'//trim(adjustl(nnc_string))//'(3x,es23.16,1x))'
    	header_format = '(181x,'//trim(adjustl(nnc_string))//'(18x,a5,4x))'
    else
   		nonnse_format = '(2x,i4,12(1x,es24.17),'//trim(nnc_string)//'(1x,es24.17))'
   		header_format = '(306x,'//trim(adjustl(nnc_string))//'(20x,a5))'
   	endif


    ! allocate and initialize star variables
    allocate (x_e_star(nx_star+1))
    allocate (dx_e_star(nx_star))
    allocate (volx_e_star(nx_star+1))
    allocate (volx_c_star(nx_star))
    allocate (dvolx_e_star(nx_star))
    allocate (dvol_e_star(nx_star))
    allocate (dmass_e_star(nx_star))
    allocate (u_c_star(nx_star))
    allocate (rho_c_star(nx_star))
    allocate (t_c_star(nx_star))
    allocate (ye_c_star(nx_star))
    allocate (xn_read(nnc_star))
    allocate (xn_c_star(nx_star,nspec))
    allocate (a_nuc_star(nnc_star))
    allocate (z_nuc_star(nnc_star))
    allocate (be_nuc_star(nnc_star))
    allocate (name_nuc_star(nnc_star))
    allocate (a_aux_c_star(nx_star))
    allocate (z_aux_c_star(nx_star))

    x_e_star = zero
    dx_e_star = zero
    volx_e_star = zero
    volx_c_star = zero
    dvolx_e_star = zero
    dmass_e_star = zero
    u_c_star = zero
    rho_c_star = zero
    t_c_star = zero
    ye_c_star = zero
    xn_read = zero
    xn_c_star = zero

    imin_star = 1
    imax_star = nx_star

    read(nread,header_format) (name_nuc_star(i),i=1,nnc_star)


    ! convert species names to lower case
    do i = 1, nnc_star
      call nuc_rename(name_nuc_star(i))
    end do


    ! get A and Z from names
    call nucaz_from_name( name_nuc_star, a_nuc_star, z_nuc_star, be_nuc_star, nnc_star )

    ! create lookup tables for isotopes in castro net
    k = 0
    net_to_castro(:) = nnc_star+1
    do i = 1, nspec
      do j = 1, nnc_star
        if ( nint(aion(i)) == nint(a_nuc_star(j)) .and. nint(zion(i)) == nint(z_nuc_star(j)) ) then
          k = k + 1
          net_to_castro(i) = j
          itmp(k) = i
          exit
        end if
      end do
      if ( j > nnc_star ) then
        if (parallel_IOProcessor()) then
          write(*,'(2(a,i3),a)') ' could not find isotope (',nint(zion(i)),',',nint(aion(i)),') in MESA net'
        end if
      end if
    end do

    if ( k > 0 ) then
      allocate (net_in_castro(k))
      net_in_castro(:) = itmp(1:k)
    else
      call bl_error("no species in star net in CASTRO net")
    end if

    ! read zones, MESA outputs these backwards (outer-most zones first)
    do i = 1, nx_star



      if (star_type == 1) then
      	read(nread,nonnse_format) zone_read, rho_read, t_read, r_read, dum1, dum2, u_read, dum3, (xn_read(k),k=1,nnc_star)
      else
      	read(nread,nonnse_format) zone_read, dum1, dum2, r_read, u_read, rho_read, t_read, dum3, dum4, dum5, dum6, dum7, dum8, (xn_read(k),k=1,nnc_star)
      endif

      if (star_type == 1) then
         j = nx_star-zone_read+2
         ! MESA gives these quantities in log
     	 rho_c_star(j-1) = exp( rho_read )
     	 t_c_star(j-1) = exp( t_read )
     	 x_e_star(j) = exp( r_read )
     	 u_c_star(j-1) = u_read

  	  else
  	  	 j = zone_read+1
  	  	 rho_c_star(j-1) = rho_read
  	  	 t_c_star(j-1) = t_read
  	  	 x_e_star(j) = r_read
  	  	 u_c_star(j-1) = u_read
  	  endif


       
      xn_c_star(j-1,net_in_castro) = xn_read(net_to_castro(net_in_castro))

      ye_c_star(j-1) = sum( z_nuc_star * xn_read / a_nuc_star ) / sum( xn_read )

      if ( nspec > nspec_evolve ) then
        dum1 = sum( xn_read,                       mask=(net_to_castro==nnc_star+1) )
        dum2 = sum( xn_read/a_nuc_star,            mask=(net_to_castro==nnc_star+1) )
        dum3 = sum( xn_read*z_nuc_star/a_nuc_star, mask=(net_to_castro==nnc_star+1) )
        xn_c_star(j-1,nspec)  = dum1
        a_aux_c_star(j-1) = dum1 / dum2
        z_aux_c_star(j-1) = dum1 / dum3
      else
        a_aux_c_star(j-1) = zero
        z_aux_c_star(j-1) = zero
      end if
    end do

!   if (parallel_IOProcessor()) then
!     write(*,'(a5,2f10.4)') (name_nuc_star(i), a_nuc_star(i), z_nuc_star(i), i=1,nnc_star)
!     write(*,'(f10.4)') (ye_c_star(i),i=1,nx_star)
!   end if

!   xn_c_star(:,nspec) = max( zero, min( one, one - sum( xn_c_star(:,1:nspec-1), dim=1 ) ) )

    close(nread)

    volx_e_star(1)     = third * x_e_star(1)**3
    do j = 1, nx_star
      dr               = x_e_star(j+1) - x_e_star(j)
      dvolx_e_star(j)  = dr * ( x_e_star(j) * x_e_star(j+1) + dr * dr * third )
      volx_e_star(j+1) = volx_e_star(j) + dvolx_e_star(j)
    end do
    volx_c_star(1:nx_star) = volx_e_star(1:nx_star) + half*dvolx_e_star(1:nx_star)

    dvol_e_star(:) = four * m_pi * dvolx_e_star(:)
    dmass_e_star(:) = dvol_e_star(:) * rho_c_star(:)
    
    ! overwrite the 'zeroth' zone
!   rho_c_star(1) = rho_c_star(2)
!   t_c_star(1) = t_c_star(2)
!   u_c_star(1) = u_c_star(2)
!   xn_c_star(1,:) = xn_c_star(2,:)

    return
  end subroutine read_star_file

  subroutine nuc_rename( nname )

    ! input variables
    character(len=5), intent(inout) :: nname

    ! local variables
    integer, parameter :: lc_a_ascii=iachar('a')
    integer, parameter :: uc_a_ascii=iachar('A')
    integer, parameter :: lc_z_ascii=iachar('z')
    integer, parameter :: uc_z_ascii=iachar('Z')

    character(len=5) :: amass, letters, cname
    integer :: name_len, in_ascii, j, k
    logical :: digit

    amass    = '     '
    letters  = '     '
    cname    = '     '

    ! copy to local
    cname    = trim(adjustl(nname))

    ! determine length of name
    name_len = len_trim(cname)

    ! take care of special cases
    if ( cname(1:name_len) == 'h1' .or. cname(1:name_len) == 'H1' .or. &
    &    cname(1:name_len) == 'h'  .or. cname(1:name_len) == 'H' ) then
      nname = 'p'
      name_len = 1
    else if ( cname(1:name_len) == 'h2' .or. cname(1:name_len) == 'H2' ) then
      nname = 'd'
      name_len = 1
    else if ( cname(1:name_len) == 'h3' .or. cname(1:name_len) == 'H3' ) then
      nname = 't'
      name_len = 1
    else if ( cname(1:name_len) == 'n01' .or. cname(1:name_len) == 'n' .or. &
    &         cname(1:name_len) == 'nt1' .or. cname(1:name_len) == 'neut' ) then
      nname = 'n'
      name_len = 1
    else if ( cname(1:name_len) == 'al-6' ) then
      nname = 'al-6'
      name_len = 4
    else if ( cname(1:name_len) == 'al*6' ) then
      nname = 'al*6'
      name_len = 4
    else if ( name_len > 1 ) then

      ! scan name until you find the digit
      digit = .false.
      do j = 1, name_len
        select case( cname(j:j) )
        case( '0':'9' )
          digit = .true.
        case default
          digit = .false.
        end select

        if ( digit ) then
          letters(1:(j-1)) = cname(1:(j-1)) ! pick out letters
          letters           = trim(adjustl(letters))

          amass = cname(j:name_len) ! pick out digits
          amass = trim(adjustl(amass))

          ! convert name to lower-case
          do k = 1, len(letters)
            in_ascii = iachar( letters(k:k) )
            if ( in_ascii >= uc_a_ascii .and. in_ascii <= uc_z_ascii ) then
              letters(k:k) = achar( in_ascii + (lc_a_ascii - uc_a_ascii) )
            end if
          end do


          ! copy back to species name array
          write(nname,'(a5)') trim(adjustl(letters))//trim(adjustl(amass))
          exit

        end if ! digit

      end do ! j = 1, name_len
    end if ! name_len > 1

    nname = adjustr( nname )

    return
  end subroutine nuc_rename

  subroutine nucaz_from_name(nname,a_nuc,z_nuc,be_nuc,nucnum)

    integer, intent(in)  :: nucnum

    character (len=5)    :: nname(nucnum)
    character (len=5)    :: cname
    character (len=124)  :: line
    character (len=5)    :: letters, amass       ! break name into its letters and its numbers
    logical              :: digit     ! determine if string is a number
    logical              :: begin     ! determine start of mass table

    real (rt)          :: a_nuc(nucnum), z_nuc(nucnum), be_nuc(nucnum)
    real (rt)          :: a_nuc_tmp, be_nuc_tmp
    integer              :: ia_nuc, name_len
    integer              :: i,j,k   ! loop variables
    integer              :: inuc
    integer              :: istat

    integer              :: in_ascii, lc_a_ascii, uc_a_ascii, lc_z_ascii, uc_z_ascii ! uppercase to lower case

    integer, parameter   :: n_masstable = 3218

    ! declare variables from mass table
    character(len=1)  :: cc 
    integer           :: nz, ni, zi, ai
    character(len=3)  :: el
    character(len=4)  :: oi
    real (rt)       :: binding
    character(len=11) :: binding_r

    integer           :: lun_table

    digit = .false. 

    ! declaration for turning capital into lowercase
    lc_a_ascii=iachar('a')
    uc_a_ascii=iachar('A')
    lc_z_ascii=iachar('z')
    uc_z_ascii=iachar('Z')

    ! open mass table to read
    open( newunit=lun_table, file="mass.mas12", status='old' )

    ! cycle over all species
    do i=1, nucnum

      amass    = '     '
      letters  = '     '
      cname    = '     '
      name_len = 0
      ia_nuc   = 0

      be_nuc_tmp = 0.0d0

      nname(i) = trim(adjustl(nname(i)))
      cname    = nname(i)
      name_len = len_trim(cname)

      rewind(lun_table)

      ! cycle over length of name to determine where numbers start
      if( cname(1:name_len) .eq. 'n' ) then
        a_nuc(i)  = 1.0d0
        z_nuc(i)  = 0.0d0
        be_nuc(i) = 0.0d0
        letters    = ' n' 
        amass      = '1'
        cname      = ' n1'
        name_len   = 3
      else if( cname(1:name_len) .eq. 'p' ) then
        a_nuc(i)  = 1.0d0
        z_nuc(i)  = 1.0d0
        be_nuc(i) = 0.0d0
        letters    = 'h' 
        amass      = '1'
        cname      = 'h1'
        name_len   = 2
      else if( cname(1:name_len) .eq. 'd' ) then
        a_nuc(i)  = 2.0d0
        z_nuc(i)  = 1.0d0
        be_nuc(i) = 1112.283 * 1.0d-3 * a_nuc(i)   ! to convert kev to mev
        letters    = 'h' 
        amass      = '2'
        cname      = 'h2'
        name_len   = 2
      else if( cname(1:name_len) .eq. 't' ) then
        a_nuc(i)  = 3.0d0
        z_nuc(i)  = 1.0d0
        be_nuc(i) = 2827.266d0 * 1.0d-3 * a_nuc(i)
        letters    = 'h' 
        amass      = '3'
        cname      = 'h3'
        name_len   = 2
      else if( cname(1:name_len) .eq. 'alg6' .or. cname(1:name_len) .eq. 'al-6' ) then
        a_nuc(i)  = 26.0d0
        z_nuc(i)  = 13.0d0
        be_nuc(i) = 8149.771d0 * 1.0d-3 * a_nuc(i)
        letters    = 'al'
        amass      = '26'
        cname      = 'al26'
        name_len   = 4
      else if( cname(1:name_len) .eq. 'alm6' .or. cname(1:name_len) .eq. 'al*6' ) then
        a_nuc(i)  = 26.0d0
        z_nuc(i)  = 13.0d0
        be_nuc(i) = 8149.771d0 * 1.0d-3 * a_nuc(i)
        letters    = 'al'
        amass      = '26'
        cname      = 'al26'
        name_len   = 4
      else if ( name_len > 1 ) then
        ! scan name until you find the digit
        digit = .false.
        do j=1,name_len
          select case( cname(j:j) )
          case( '0':'9' )
            digit = .true.
          case default
            digit = .false.
          end select

          if( digit ) then
            letters(1:(j-1))  = cname(1:(j-1))    ! pick out letters
            letters            = trim(adjustl(letters))

            amass = cname(j:name_len) ! pick out digits
            amass = trim(adjustl(amass))

            read( amass, '(i5)' ) ia_nuc      ! assign character mass to real variable            
            exit
          end if
        end do

        ! scan through table until we reach the desired species
        begin = .false.
        do k = 1,n_masstable
          ! read in line-by-line
          read( lun_table, '(a124)', iostat=istat ) line
          if( istat < 0 ) exit  ! end-of-file

          ! cycle through the header lines
          if( .not. begin ) then
            ! neutron is first species in mass table
            if( line(1:19) == "0  1    1    0    1" ) then
              begin = .true.
            else
              cycle
            end if
          end if
          if( istat > 0 ) cycle ! read error

          666 format( a1,i3,i5,i5,i5,1x,a3,a4,1x,13x,11x,a11,9x,1x,2x,11x,9x,1x,16x,11x,1x )
          read( line, 666 ) cc, nz, ni, zi, ai, el, oi, binding_r

          ! convert "#" character to decimal
          j = index( binding_r, "#" )
          if( j /= 0 ) binding_r(j:j) = "."
          read( binding_r, "(f11.3)" ) binding

          ! convert name to lower-case
          el = trim(adjustl(el))
          do j = 1,len(el)
            in_ascii = iachar( el(j:j) )
            if( in_ascii >= uc_a_ascii .and. in_ascii <= uc_z_ascii ) then
              el(j:j) = achar( in_ascii + (lc_a_ascii - uc_a_ascii) )
            end if
          end do

          ! match names in mass table to network species names
          if( el .eq. letters ) then
            z_nuc(i)  = dble(zi)
            a_nuc_tmp  = dble(ai)
            be_nuc_tmp = binding*a_nuc_tmp*1.0d-3 ! b.e. from table is [kev / a]
            if( ai .eq. ia_nuc ) then
              a_nuc(i)  = a_nuc_tmp
              be_nuc(i) = be_nuc_tmp
              exit
            end if
          end if

          if( k .eq. n_masstable .and. a_nuc(i) .ne. dble(ia_nuc) ) then
            a_nuc(i)  = dble(ia_nuc)
            be_nuc(i) = be_nuc_tmp
          end if
        end do ! do k = 1,n_masstable
      end if
    end do  ! do i = 1,nucnum

    close(lun_table)  ! close mass table file

    return
  end subroutine nucaz_from_name

  subroutine interp1dvol_star( x_out, state_star, state_out )

    use bl_constants_module
    use bl_error_module
    use interpolate_module, only: locate
    use model_interp_module, only: interp1d_linear, interp1d_spline
    use probdata_module, only: interp_method

    ! input variables
    real (rt), intent(in) :: x_out(:)
    real (rt), intent(in) :: state_star(:)

    ! output variables
    real (rt), intent(out) :: state_out(size(x_out))

    ! local variables
    integer :: ix(size(x_out))
    integer :: ix_max

    real (rt) :: dvol1, dvol2
    real (rt) :: volx_out(size(x_out))

    integer :: i, j, n

    volx_out = third * x_out**3

    if ( interp_method == 1 ) then

      ix_max = imax_star
      ix(:) = 1
      do i = 1, size(x_out)
        if ( volx_out(i) <= volx_c_star(1) ) then
          ix(i) = 0
        else if ( volx_out(i) >= volx_c_star(ix_max) ) then
          ix(i) = ix_max
        else
          ix(i) = locate( volx_out(i), ix_max, volx_c_star ) - 1
        end if
      end do
      call interp1d_linear( ix, ix_max, volx_c_star, state_star, volx_out, state_out )
    else if ( interp_method == 2 ) then
      call interp1d_spline( volx_c_star, state_star, volx_out, state_out )
    else
      call bl_error("invalid value for interp_method")
    end if

    return
  end subroutine interp1dvol_star


  subroutine interp1drad_star( x_out, state_star, state_out )

    use bl_constants_module
    use bl_error_module
    use interpolate_module, only: locate
    use model_interp_module, only: interp1d_linear, interp1d_spline
    use probdata_module, only: interp_method

    ! input variables
    real (rt), intent(in) :: x_out(:)
    real (rt), intent(in) :: state_star(:)

    ! output variables
    real (rt), intent(out) :: state_out(size(x_out))

    ! local variables
    integer :: ix(size(x_out))
    integer :: ix_max

    integer :: i, j, n


    if ( interp_method == 1 ) then

      ix_max = imax_star
      ix(:) = 1
      do i = 1, size(x_out)
        if ( x_out(i) <= x_e_star(1) ) then
          ix(i) = 0
        else if ( x_out(i) >= x_e_star(ix_max) ) then
          ix(i) = ix_max
        else
          ix(i) = locate( x_out(i), ix_max, x_e_star ) - 1
        end if
      end do
      call interp1d_linear( ix, ix_max, x_e_star, state_star, x_out, state_out )
    else if ( interp_method == 2 ) then
      call interp1d_spline( x_e_star, state_star, x_out, state_out )
    else
      call bl_error("invalid value for interp_method")
    end if

    return
  end subroutine interp1drad_star


  subroutine interp2dvol_star( x_out, state_star, state_out )

    use bl_constants_module
    use bl_error_module
    use interpolate_module, only: locate
    use model_interp_module, only: interp1d_linear, interp1d_spline
    use probdata_module, only: interp_method

    ! input variables
    real (rt), intent(in) :: x_out(:,:)
    real (rt), intent(in) :: state_star(:)

    ! output variables
    real (rt), intent(out) :: state_out(size(x_out,1),size(x_out,2))

    ! local variables
    integer :: ix(size(x_out,1))
    integer :: ix_max

    real (rt) :: dvol1, dvol2
    real (rt) :: volx_out(size(x_out,1),size(x_out,2))

    integer :: i, j, n

    volx_out = third * x_out**3

    if ( interp_method == 1 ) then

      ix_max = imax_star
      do j = 1, size(x_out,2)
        ix(:) = 1
        do i = 1, size(x_out,1)
          if ( volx_out(i,j) <= volx_c_star(1) ) then
            ix(i) = 0
          else if ( volx_out(i,j) >= volx_c_star(ix_max) ) then
            ix(i) = ix_max
          else
            ix(i) = locate( volx_out(i,j), ix_max, volx_c_star ) - 1
          end if
        end do
        call interp1d_linear( ix, ix_max, volx_c_star, state_star, volx_out(:,j), state_out(:,j) )
      end do
    else if ( interp_method == 2 ) then
      do j = 1, size(x_out,2)
        call interp1d_spline( volx_c_star, state_star, volx_out(:,j), state_out(:,j) )
      end do
    else
      call bl_error("invalid value for interp_method")
    end if

    return
  end subroutine interp2dvol_star

  subroutine interp2drad_star( x_out, state_star, state_out )

    use bl_constants_module
    use bl_error_module
    use interpolate_module, only: locate
    use model_interp_module, only: interp1d_linear, interp1d_spline
    use probdata_module, only: interp_method

    ! input variables
    real (rt), intent(in) :: x_out(:,:)
    real (rt), intent(in) :: state_star(:)

    ! output variables
    real (rt), intent(out) :: state_out(size(x_out,1),size(x_out,2))

    ! local variables
    integer :: ix(size(x_out,1))
    integer :: ix_max

    integer :: i, j, n


    if ( interp_method == 1 ) then

      ix_max = imax_star
      do j = 1, size(x_out,2)
        ix(:) = 1
        do i = 1, size(x_out,1)
          if ( x_out(i,j) <= x_e_star(1) ) then
            ix(i) = 0
          else if ( x_out(i,j) >= x_e_star(ix_max) ) then
            ix(i) = ix_max
          else
            ix(i) = locate( x_out(i,j), ix_max, x_e_star ) - 1
          end if
        end do
        call interp1d_linear( ix, ix_max, x_e_star, state_star, x_out(:,j), state_out(:,j) )
      end do
    else if ( interp_method == 2 ) then
      do j = 1, size(x_out,2)
        call interp1d_spline( x_e_star, state_star, x_out(:,j), state_out(:,j) )
      end do
    else
      call bl_error("invalid value for interp_method")
    end if

    return
  end subroutine interp2drad_star


  subroutine interp3dvol_star( x_out, state_star, state_out )

    use bl_constants_module
    use bl_error_module
    use interpolate_module, only: locate
    use model_interp_module, only: interp1d_linear, interp1d_spline
    use probdata_module, only: interp_method

    ! input variables
    real (rt), intent(in) :: x_out(:,:,:)
    real (rt), intent(in) :: state_star(:)

    ! output variables
    real (rt), intent(out) :: state_out(size(x_out,1),size(x_out,2),size(x_out,3))

    ! local variables
    integer :: ix(size(x_out,1))
    integer :: ix_max

    real (rt) :: volx_out(size(x_out,1),size(x_out,2),size(x_out,3))

    integer :: i, j, k

    volx_out   = third * x_out**3

    if ( interp_method == 1 ) then

      ix_max = imax_star
      do k = 1, size(x_out,3)
        do j = 1, size(x_out,2)
          ix(:) = 1
          do i = 1, size(x_out,1)
            if ( volx_out(i,j,k) <= volx_c_star(1) ) then
              ix(i) = 0
            else if ( volx_out(i,j,k) >= volx_c_star(ix_max) ) then
              ix(i) = ix_max
            else
              ix(i) = locate( volx_out(i,j,k), ix_max, volx_c_star ) - 1
            end if
          end do
          call interp1d_linear( ix, ix_max, volx_c_star, state_star, volx_out(:,j,k), state_out(:,j,k) )
        end do
      end do
    else if ( interp_method == 2 ) then
      do k = 1, size(x_out,3)
        do j = 1, size(x_out,2)
          call interp1d_spline( volx_c_star, state_star, volx_out(:,j,k), state_out(:,j,k) )
        end do
      end do
    else
      call bl_error("invalid value for interp_method")
    end if

    return
  end subroutine interp3dvol_star

  subroutine interp3drad_star( x_out, state_star, state_out )

    use bl_constants_module
    use bl_error_module
    use interpolate_module, only: locate
    use model_interp_module, only: interp1d_linear, interp1d_spline
    use probdata_module, only: interp_method

    ! input variables
    real (rt), intent(in) :: x_out(:,:,:)
    real (rt), intent(in) :: state_star(:)

    ! output variables
    real (rt), intent(out) :: state_out(size(x_out,1),size(x_out,2),size(x_out,3))

    ! local variables
    integer :: ix(size(x_out,1))
    integer :: ix_max

    integer :: i, j, k

    if ( interp_method == 1 ) then

      ix_max = imax_star
      do k = 1, size(x_out,3)
        do j = 1, size(x_out,2)
          ix(:) = 1
          do i = 1, size(x_out,1)
            if ( x_out(i,j,k) <= x_e_star(1) ) then
              ix(i) = 0
            else if ( x_out(i,j,k) >= x_e_star(ix_max) ) then
              ix(i) = ix_max
            else
              ix(i) = locate( x_out(i,j,k), ix_max, x_e_star ) - 1
            end if
          end do
          call interp1d_linear( ix, ix_max, x_e_star, state_star, x_out(:,j,k), state_out(:,j,k) )
        end do
      end do
    else if ( interp_method == 2 ) then
      do k = 1, size(x_out,3)
        do j = 1, size(x_out,2)
          call interp1d_spline( x_e_star, state_star, x_out(:,j,k), state_out(:,j,k) )
        end do
      end do
    else
      call bl_error("invalid value for interp_method")
    end if

    return
  end subroutine interp3drad_star


end module star_parser_module
