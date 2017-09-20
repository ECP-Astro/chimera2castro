module kepler_parser_module

  use parallel, only: parallel_IOProcessor
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer, save :: nx_kep, nnc_kep ! max num. zones, num. nuclei in progenitor set
  integer, save :: imin_kep, imax_kep ! min zone, num. of zones
  integer, save :: i_nse_kep !last zone in use


  character(len=5), allocatable, save :: nuc_name_kep(:) ! nuclei names
  real (rt), allocatable, save :: xn_full_kep(:) ! mass fractions
  real (rt), allocatable, save :: m0_kep(:) ! enclosed mass [g]
  real (rt), allocatable, save :: dm0_kep(:) ! cell mass [g]
  real (rt), allocatable, save :: vol_kep(:) ! cell volume [cm^3]
  real (rt), allocatable, save :: r0_kep(:) ! radius [cm]
  real (rt), allocatable, save :: rho0_kep(:) ! density [g cm^{-3}]
  real (rt), allocatable, save :: temp0_kep(:) ! temperature [K]
  real (rt), allocatable, save :: u0_kep(:) ! velocity [cm s^{-1}]
  real (rt), allocatable, save :: ye0_kep(:) ! electron fraction
  real (rt), allocatable, save :: p0_kep(:) ! pressure [erg cm^{-3}]
  real (rt), allocatable, save :: e0_kep(:) ! internal energy [erg g^{-1}]
  real (rt), allocatable, save :: s0_kep(:) ! entropy 
  real (rt), allocatable, save :: ang_vel0_kep(:) ! angular velocity 


  real (rt), allocatable, save :: a_nuc_kep(:)
  real (rt), allocatable, save :: z_nuc_kep(:)
  real (rt), allocatable, save :: a_aux_nuc_kep(:)
  real (rt), allocatable, save :: z_aux_nuc_kep(:)
  real (rt), allocatable, save :: be_aux_nuc_kep(:)


  real (rt), allocatable, save :: volx_e_kep(:)
  real (rt), allocatable, save :: volx_c_kep(:)
  real (rt), allocatable, save :: dvolx_e_kep(:)
  real (rt), allocatable, save :: dvolx_c_kep(:)





  ! ---unused MESA vars from conversion-----

  ! real (rt), allocatable, save :: x_e_kep(:)
  ! real (rt), allocatable, save :: dx_e_kep(:)
  ! real (rt), allocatable, save :: volx_e_kep(:)
  ! real (rt), allocatable, save :: volx_c_kep(:)
  ! real (rt), allocatable, save :: dvolx_e_kep(:)
  ! real (rt), allocatable, save :: dvolx_c_kep(:)
  ! real (rt), allocatable, save :: dvol_e_kep(:)
  ! real (rt), allocatable, save :: dmass_e_kep(:)

  ! real (rt), allocatable, save :: u_c_kep(:)
  ! real (rt), allocatable, save :: rho_c_kep(:)
  ! real (rt), allocatable, save :: t_c_kep(:)
  ! real (rt), allocatable, save :: ye_c_kep(:)

  ! real (rt), allocatable, save :: xn_c_kep(:,:)
  ! real (rt), allocatable, save :: a_nuc_kep(:)
  ! real (rt), allocatable, save :: z_nuc_kep(:)
  ! real (rt), allocatable, save :: be_nuc_kep(:)
  ! character(len=5), allocatable, save :: nuc_name_kep(:)

  ! real (rt), allocatable, save :: a_aux_c_kep(:)
  ! real (rt), allocatable, save :: z_aux_c_kep(:)

  ! --------------------




  
  contains

  subroutine read_kepler_file( filename )

    use actual_network, only: nspec, nspec_evolve, aion, zion
    use bl_constants_module

    ! input variables
    character(len=*), intent(in) :: filename

    ! local variables
    character(len=1600)      :: line           ! read line for structure data

    integer                :: nread     ! file unit number
    integer                :: zone_read
    real (rt)              :: dm0 ! single zone mass
    real (rt)              :: abar0 ! mean nucleon number
    real (rt), allocatable :: xn_read(:)     ! read in progenitor composition array



    real (rt) :: dr, dvol

    integer :: itmp(nspec)
    integer :: net_to_castro(nspec)
    integer, allocatable :: net_in_castro(:)

    integer :: i, j, k   ! loop indices
    integer :: nread        ! file unit number
    integer :: istate       ! iostat variable

    character(len=32) :: header_name
    character(len=4)  :: nnc_string
    character(len=20) :: header_format
    character(len=36) :: nonnse_format

    ! formats
    105 format (6x,12(1x,es24.17)) ! get cell vars (non-species info) 
    4128 format (i5,7(1x,es24.16e3))
    1012 format(39x,i5)

    ! open input data file
    open( newunit=nread, file=trim(adjustl(filename)), status='old', iostat=istate )
    if ( istate /= 0 ) then
      call bl_error('Aborting now -- please supply KEPLER file')
    end if

    ! read progenitor file
    do
      read(nread,'(a)',iostat=istate) line
      if ( istate /= 0 ) exit

      header_name = line(1:32)

      if ( trim(adjustl(header_name)) == 'num_zones' ) then
        read(line,1012) nx_kep
        cycle
      end if

      if ( trim(adjustl(header_name)) == 'num_species' ) then
        read(line,1012) nnc_kep
        exit
      end if
    end do
    read(nread,*) line
    read(nread,*) line

    write(nnc_string,'(i4)') nnc_kep
   ! nonnse_format = '(i5,1x,7(3x,es23.16,1x),'//trim(adjustl(nnc_string))//'(3x,es23.16,1x))'
   ! header_format = '(181x,'//trim(adjustl(nnc_string))//'(18x,a5,4x))'
    nonnse_format = '(6x,12(1x,es24.17),'//trim(nnc_string)//'(1x,es24.17))'
    header_format = '(306x,'//trim(nnc_string)//'(20x,a5))'


     allocate(xn_full_kep(nx_kep+1,nnc_kep)) ! mass fractions
     allocate(m0_kep(nx_kep)) ! enclosed mass [g]
     allocate(dm0_kep(nx_kep)) ! cell mass [g]
     allocate(r0_kep(nx_kep)) ! radius [cm]
     allocate(vol_kep(nx_kep)) ! cell volume [cm^3]
     allocate(rho0_kep(nx_kep)) ! density [g cm^{-3}]
     allocate(temp0_kep(nx_kep)) ! temperature [K]
     allocate(u0_kep(nx_kep)) ! velocity [cm s^{-1}]
     allocate(ye0_kep(nx_kep)) ! electron fraction
     allocate(p0_kep(nx_kep)) ! pressure [erg cm^{-3}]
     allocate(e0_kep(nx_kep)) ! internal energy [erg g^{-1}]
     allocate(s0_kep(nx_kep)) ! entropy 
     allocate(ang_vel0_kep(nx_kep)) ! angular velocity 
     allocate(a_nuc_kep(nnc_kep)) 
     allocate(z_nuc_kep(nnc_kep))
     allocate(be_nuc_kep(nnc_kep))
     allocate(a_aux_nuc_kep(nnc_kep))
     allocate(z_aux_nuc_kep(nnc_kep))
     allocate(be_aux_nuc_kep(nnc_kep))

     xn_full_kep = zero
     m0_kep = zero
     dm0_kep = zero
     vol_kep = zero
     r0_kep = zero
     rho0_kep = zero
     temp0_kep = zero
     u0_kep = zero
     ye0_kep = zero
     p0_kep = zero
     e0_kep = zero
     s0_kep = zero
     ang_vel0_kep = zero
     a_nuc_kep = zero
     z_nuc_kep = zero
     be_nuc_kep = zero



    read(nread,header_format) nuc_name_kep

    ! convert species names to lower case
    do i = 1, nnc_kep
      call nuc_rename(nuc_name_kep(i))
    end do

    ! get A and Z from names
    call nucaz_from_name( nuc_name_kep, a_nuc_kep, z_nuc_kep, be_nuc_kep, nnc_kep )

    ! create lookup tables for isotopes in castro net
    k = 0
    net_to_castro(:) = nnc_kep+1
    do i = 1, nspec
      do j = 1, nnc_kep
        if ( nint(aion(i)) == nint(a_nuc_kep(j)) .and. nint(zion(i)) == nint(z_nuc_kep(j)) ) then
          k = k + 1
          net_to_castro(i) = j
          itmp(k) = i
          exit
        end if
      end do
      if ( j > nnc_kep ) then
        if (parallel_IOProcessor()) then
          write(*,'(2(a,i3),a)') ' could not find isotope (',nint(zion(i)),',',nint(aion(i)),') in KEPLER net'
        end if
      end if
    end do

    if ( k > 0 ) then
      allocate (net_in_castro(k))
      net_in_castro(:) = itmp(1:k)
    else
      call bl_error("no species in KEPLER net in CASTRO net")
    end if

    ! read zones
    do i = 1, nx_kep

      read(nread,nonnse_format) zone_read, dm0_kep[i], m0_kep[i], r0_kep[i], u0_kep[i], roh0_kep[i], &
                                 temp0_kep[i], p0_kep[i], e0_kep[i], s0_kep[i], ang_vel0_kep[i], &
                                  dum1, ye0_kep[i], (xn_read(k),k=1,nnc_kep)

      xn_full_kep(i,net_in_castro) = xn_read(net_to_castro(net_in_castro))





      if ( nspec > nspec_evolve ) then
        dum1 = sum( xn_read,                       mask=(net_to_castro==nnc_kep+1) )
        dum2 = sum( xn_read/a_nuc_kep,             mask=(net_to_castro==nnc_kep+1) )
        dum3 = sum( xn_read*z_nuc_kep/a_nuc_kep,   mask=(net_to_castro==nnc_kep+1) )
        xn_full_kep(i,nspec)  = dum1
        a_aux_nuc_kep(i) = dum1 / dum2
        z_aux_nuc_kep(i) = dum1 / dum3
      else
        a_aux_nuc_kep(i) = zero
        z_aux_nuc_kep(i) = zero
      end if
    end do

    close(nread)

    volx_e_kep(1)     = third*r0_kep(1)**3
    do j = 1, nx_kep
      dr = r0_kep(j+1)-r0_kep(j)
      dvolx_e_kep(j) = dr*(r0_kep(j)*r0_kep(j+1)+dr*dr*third)
      volx_e_kep(j+1) = volx_e_kep(j)+dvolx_e_kep(j)
    end do

    volx_c_kep(1:nx_kep) = volx_e_kep(1:nx_kep) + half*dvolx_e_kep(1:nx_kep)

    dvol_e_kep(:) = four*m_pi*dvolx_e_kep(:)


    return
  end subroutine read_kep_file

  subroutine nuc_rename( nname )

    ! input variables
    character(len=5), intent(inout) :: nname

    ! local variables
    integer, parameter :: lc_a_ascii=iachar('a')
    integer, parameter :: uc_a_ascii=iachar('a')
    integer, parameter :: lc_z_ascii=iachar('z')
    integer, parameter :: uc_z_ascii=iachar('z')

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

  subroutine interp1d_kep( x_out, state_kep, state_out )

    use bl_constants_module
    use bl_error_module
    use interpolate_module, only: locate
    use model_interp_module, only: interp1d_linear, interp1d_spline
    use probdata_module, only: interp_method

    ! input variables
    real (rt), intent(in) :: x_out(:)
    real (rt), intent(in) :: state_kep(:)

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

      ix_max = imax_kep
      ix(:) = 1
      do i = 1, size(x_out)
        if ( volx_out(i) <= volx_c_kep(1) ) then
          ix(i) = 0
        else if ( volx_out(i) >= volx_c_kep(ix_max) ) then
          ix(i) = ix_max
        else
          ix(i) = locate( volx_out(i), ix_max, volx_c_kep ) - 1
        end if
      end do
      call interp1d_linear( ix, ix_max, volx_c_kep, state_kep, volx_out, state_out )
    else if ( interp_method == 2 ) then
      call interp1d_spline( volx_c_kep, state_kep, volx_out, state_out )
    else
      call bl_error("invalid value for interp_method")
    end if

    return
  end subroutine interp1d_kep

  subroutine interp2d_kep( x_out, state_kep, state_out )

    use bl_constants_module
    use bl_error_module
    use interpolate_module, only: locate
    use model_interp_module, only: interp1d_linear, interp1d_spline
    use probdata_module, only: interp_method

    ! input variables
    real (rt), intent(in) :: x_out(:,:)
    real (rt), intent(in) :: state_kep(:)

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

      ix_max = imax_kep
      do j = 1, size(x_out,2)
        ix(:) = 1
        do i = 1, size(x_out,1)
          if ( volx_out(i,j) <= volx_c_kep(1) ) then
            ix(i) = 0
          else if ( volx_out(i,j) >= volx_c_kep(ix_max) ) then
            ix(i) = ix_max
          else
            ix(i) = locate( volx_out(i,j), ix_max, volx_c_kep ) - 1
          end if
        end do
        call interp1d_linear( ix, ix_max, volx_c_kep, state_kep, volx_out(:,j), state_out(:,j) )
      end do
    else if ( interp_method == 2 ) then
      do j = 1, size(x_out,2)
        call interp1d_spline( volx_c_kep, state_kep, volx_out(:,j), state_out(:,j) )
      end do
    else
      call bl_error("invalid value for interp_method")
    end if

    return
  end subroutine interp2d_kep

  subroutine interp3d_kep( x_out, state_kep, state_out )

    use bl_constants_module
    use bl_error_module
    use interpolate_module, only: locate
    use model_interp_module, only: interp1d_linear, interp1d_spline
    use probdata_module, only: interp_method

    ! input variables
    real (rt), intent(in) :: x_out(:,:,:)
    real (rt), intent(in) :: state_kep(:)

    ! output variables
    real (rt), intent(out) :: state_out(size(x_out,1),size(x_out,2),size(x_out,3))

    ! local variables
    integer :: ix(size(x_out,1))
    integer :: ix_max

    real (rt) :: volx_out(size(x_out,1),size(x_out,2),size(x_out,3))

    integer :: i, j, k

    volx_out   = third * x_out**3

    if ( interp_method == 1 ) then

      ix_max = imax_kep
      do k = 1, size(x_out,3)
        do j = 1, size(x_out,2)
          ix(:) = 1
          do i = 1, size(x_out,1)
            if ( volx_out(i,j,k) <= volx_c_kep(1) ) then
              ix(i) = 0
            else if ( volx_out(i,j,k) >= volx_c_kep(ix_max) ) then
              ix(i) = ix_max
            else
              ix(i) = locate( volx_out(i,j,k), ix_max, volx_c_kep ) - 1
            end if
          end do
          call interp1d_linear( ix, ix_max, volx_c_kep, state_kep, volx_out(:,j,k), state_out(:,j,k) )
        end do
      end do
    else if ( interp_method == 2 ) then
      do k = 1, size(x_out,3)
        do j = 1, size(x_out,2)
          call interp1d_spline( volx_c_kep, state_kep, volx_out(:,j,k), state_out(:,j,k) )
        end do
      end do
    else
      call bl_error("invalid value for interp_method")
    end if

    return
  end subroutine interp3d_kep

end module kep_parser_module
