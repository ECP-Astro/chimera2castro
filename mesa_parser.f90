module mesa_parser_module

  use bl_types

  implicit none

  integer, save :: nx_mesa_in
  integer, save :: nnc_mesa_in

  real (dp_t), allocatable, save :: dens_mesa_in(:)
  real (dp_t), allocatable, save :: temp_mesa_in(:)
  real (dp_t), allocatable, save :: rad_edge_mesa_in(:)
  real (dp_t), allocatable, save :: dvol_rad_mesa_in(:)
  real (dp_t), allocatable, save :: vol_rad_edge_mesa_in(:)
  real (dp_t), allocatable, save :: vol_rad_cntr_mesa_in(:)
  real (dp_t), allocatable, save :: vrad_mesa_in(:)
  real (dp_t), allocatable, save :: xn_mesa_in(:,:)

  real (dp_t), allocatable, save :: a_nuc_mesa_in(:)
  real (dp_t), allocatable, save :: z_nuc_mesa_in(:)
  real (dp_t), allocatable, save :: be_nuc_mesa_in(:)
  character(len=5), allocatable, save :: nuc_name_mesa(:)
  
  contains

  subroutine read_mesa_file( filename )

    use actual_network, only: nspec, aion, zion
    use bl_constants_module

    ! input variables
    character(len=*), intent(in) :: filename

    ! local variables
    character(len=1600)      :: line           ! read line for structure data

    integer                  :: zone_read      ! read in zone number
    real (dp_t)              :: lnrho_read     ! read in ln(density) [g cm^{-3}]
    real (dp_t)              :: lnt_read       ! read in ln(temperature) [k]
    real (dp_t)              :: lnr_read       ! read in ln(radius) [cm]
    real (dp_t)              :: u_read         ! read in velocity [cm s^{-1}]
    real (dp_t), allocatable :: xn_read(:)     ! read in progenitor composition array

    real (dp_t) :: dum1    ! dummy variable
    real (dp_t) :: dum2    ! dummy variable
    real (dp_t) :: dum3    ! dummy variable

    real (dp_t) :: dr, dvol

    integer :: itmp(nspec)
    integer :: net_to_castro(nspec)
    integer, allocatable :: net_in_castro(:)

    integer :: ii, jj, kk   ! loop indices
    integer :: nread        ! file unit number
    integer :: istate       ! iostat variable

    character(len=32) :: header_name
    character(len=4)  :: nnc_string
    character(len=22) :: header_format
    character(len=44) :: nonnse_format

    ! formats
    1012 format(59x,i5)

    ! open input data file
    open(newunit=nread, file=trim(adjustl(filename)), status='old', iostat=istate)
    if ( istate /= 0 ) then
      call bl_error('Aborting now -- please supply MESA file')
    end if

    ! read progenitor file
    do
      read(nread,'(a)',iostat=istate) line
      if ( istate /= 0 ) exit

      header_name = line(1:32)

      if ( trim(adjustl(header_name)) == 'n_shells' ) then
        read(line,1012) nx_mesa_in
        cycle
      end if

      if ( trim(adjustl(header_name)) == 'species' ) then
        read(line,1012) nnc_mesa_in
        exit
      end if
    end do
    read(nread,*)

    write(nnc_string,'(i4)') nnc_mesa_in
    nonnse_format = '(i5,1x,7(3x,es23.16,1x),'//trim(adjustl(nnc_string))//'(3x,es23.16,1x))'
    header_format = '(181x,'//trim(adjustl(nnc_string))//'(18x,a5,4x))'

    ! allocate and initialize mesa variables
    allocate (dens_mesa_in(nx_mesa_in))
    allocate (temp_mesa_in(nx_mesa_in))
    allocate (rad_edge_mesa_in(nx_mesa_in+1))
    allocate (dvol_rad_mesa_in(nx_mesa_in+1))
    allocate (vol_rad_edge_mesa_in(nx_mesa_in+1))
    allocate (vol_rad_cntr_mesa_in(nx_mesa_in))
    allocate (vrad_mesa_in(nx_mesa_in))
    allocate (xn_mesa_in(nx_mesa_in,nspec))
    allocate (xn_read(nnc_mesa_in))
    allocate (a_nuc_mesa_in(nnc_mesa_in))
    allocate (z_nuc_mesa_in(nnc_mesa_in))
    allocate (be_nuc_mesa_in(nnc_mesa_in))
    allocate (nuc_name_mesa(nnc_mesa_in))

    dens_mesa_in = zero
    temp_mesa_in = zero
    rad_edge_mesa_in = zero
    dvol_rad_mesa_in = zero
    vol_rad_edge_mesa_in = zero
    vol_rad_cntr_mesa_in = zero
    vrad_mesa_in = zero
    xn_mesa_in = zero
    xn_read = zero

    read(nread,header_format) (nuc_name_mesa(ii),ii=1,nnc_mesa_in)

    ! convert species names to lower case
    do ii = 1, nnc_mesa_in
      call nuc_rename(nuc_name_mesa(ii))
    end do

    ! get A and Z from names
    call nucaz_from_name( nuc_name_mesa, a_nuc_mesa_in, z_nuc_mesa_in, be_nuc_mesa_in, nnc_mesa_in )

    ! create lookup tables for isotopes in castro net
    kk = 0
    net_to_castro(:) = nnc_mesa_in+1
    do ii = 1, nspec
      do jj = 1, nnc_mesa_in
        if ( nint(aion(ii)) == nint(a_nuc_mesa_in(jj)) .and. nint(zion(ii)) == nint(z_nuc_mesa_in(jj)) ) then
          kk = kk + 1
          net_to_castro(ii) = jj
          itmp(kk) = ii
          exit
        end if
      end do
      if ( jj > nnc_mesa_in ) then
        write(*,'(2(a,i3),a)') ' could not find isotope (',nint(zion(ii)),',',nint(aion(ii)),') in MESA net'
      end if
    end do

    if ( kk > 0 ) then
      allocate (net_in_castro(kk))
      net_in_castro(:) = itmp(1:kk)
    else
      call bl_error("no species in MESA net in CASTRO net")
    end if

    ! read zones, MESA outputs these backwards (outer-most zones first)
    do ii = 1, nx_mesa_in

      read(nread,nonnse_format) zone_read, lnrho_read, lnt_read, lnr_read, dum1, dum2, u_read, dum3, (xn_read(kk),kk=1,nnc_mesa_in)
      jj = nx_mesa_in-zone_read+2
      dens_mesa_in(jj-1) = exp( lnrho_read )
      temp_mesa_in(jj-1) = exp( lnt_read )
      rad_edge_mesa_in(jj) = exp( lnr_read )
      vrad_mesa_in(jj-1) = u_read
       
      xn_mesa_in(jj-1,net_in_castro) = xn_read(net_to_castro(net_in_castro))
       
    end do

    close(nread)

    vol_rad_edge_mesa_in(1) = third * rad_edge_mesa_in(1)**3
    do jj = 2, nx_mesa_in+1
      dr                       = rad_edge_mesa_in(jj) - rad_edge_mesa_in(jj-1)
      dvol                     = dr * ( rad_edge_mesa_in(jj-1) * rad_edge_mesa_in(jj) + dr * dr * third )
      dvol_rad_mesa_in(jj)     = dvol
      vol_rad_edge_mesa_in(jj) = vol_rad_edge_mesa_in(jj-1) + dvol
    end do
    vol_rad_cntr_mesa_in(1:nx_mesa_in) = vol_rad_edge_mesa_in(1:nx_mesa_in) + half*dvol_rad_mesa_in(2:nx_mesa_in+1)

    ! overwrite the 'zeroth' zone
!   dens_mesa_in(1) = dens_mesa_in(2)
!   temp_mesa_in(1) = temp_mesa_in(2)
!   vrad_mesa_in(1) = vrad_mesa_in(2)
!   xn_mesa_in(1,:) = xn_mesa_in(2,:)

    return
  end subroutine read_mesa_file

  subroutine nuc_rename( nname )

    ! input variables
    character(len=5), intent(inout) :: nname

    ! local variables
    integer, parameter :: lc_a_ascii=iachar('a')
    integer, parameter :: uc_a_ascii=iachar('a')
    integer, parameter :: lc_z_ascii=iachar('z')
    integer, parameter :: uc_z_ascii=iachar('z')

    character(len=5) :: amass, letters, cname
    integer :: name_len, in_ascii, jj, kk
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
      do jj = 1, name_len
        select case( cname(jj:jj) )
        case( '0':'9' )
          digit = .true.
        case default
          digit = .false.
        end select

        if ( digit ) then
          letters(1:(jj-1)) = cname(1:(jj-1)) ! pick out letters
          letters           = trim(adjustl(letters))

          amass = cname(jj:name_len) ! pick out digits
          amass = trim(adjustl(amass))

          ! convert name to lower-case
          do kk = 1, len(letters)
            in_ascii = iachar( letters(kk:kk) )
            if ( in_ascii >= uc_a_ascii .and. in_ascii <= uc_z_ascii ) then
              letters(kk:kk) = achar( in_ascii + (lc_a_ascii - uc_a_ascii) )
            end if
          end do

          ! copy back to species name array
          write(nname,'(a5)') trim(adjustl(letters))//trim(adjustl(amass))

          exit

        end if ! digit

      end do ! jj = 1, name_len
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

    real (dp_t)          :: a_nuc(nucnum), z_nuc(nucnum), be_nuc(nucnum)
    real (dp_t)          :: a_nuc_tmp, be_nuc_tmp
    integer              :: ia_nuc, name_len
    integer              :: ii,jj,kk   ! loop variables
    integer              :: inuc
    integer              :: istat

    integer              :: in_ascii, lc_a_ascii, uc_a_ascii, lc_z_ascii, uc_z_ascii ! uppercase to lower case

    integer, parameter   :: n_masstable = 3218

    ! declare variables from mass table
    character(len=1)  :: cc 
    integer           :: nz, ni, zi, ai
    character(len=3)  :: el
    character(len=4)  :: oi
    real (dp_t)       :: binding
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
    do ii=1, nucnum

      amass    = '     '
      letters  = '     '
      cname    = '     '
      name_len = 0
      ia_nuc   = 0

      be_nuc_tmp = 0.0d0

      nname(ii) = trim(adjustl(nname(ii)))
      cname    = nname(ii)
      name_len = len_trim(cname)

      rewind(lun_table)

      ! cycle over length of name to determine where numbers start
      if( cname(1:name_len) .eq. 'n' ) then
        a_nuc(ii)  = 1.0d0
        z_nuc(ii)  = 0.0d0
        be_nuc(ii) = 0.0d0
        letters    = ' n' 
        amass      = '1'
        cname      = ' n1'
        name_len   = 3
      else if( cname(1:name_len) .eq. 'p' ) then
        a_nuc(ii)  = 1.0d0
        z_nuc(ii)  = 1.0d0
        be_nuc(ii) = 0.0d0
        letters    = 'h' 
        amass      = '1'
        cname      = 'h1'
        name_len   = 2
      else if( cname(1:name_len) .eq. 'd' ) then
        a_nuc(ii)  = 2.0d0
        z_nuc(ii)  = 1.0d0
        be_nuc(ii) = 1112.283 * 1.0d-3 * a_nuc(ii)   ! to convert kev to mev
        letters    = 'h' 
        amass      = '2'
        cname      = 'h2'
        name_len   = 2
      else if( cname(1:name_len) .eq. 't' ) then
        a_nuc(ii)  = 3.0d0
        z_nuc(ii)  = 1.0d0
        be_nuc(ii) = 2827.266d0 * 1.0d-3 * a_nuc(ii)
        letters    = 'h' 
        amass      = '3'
        cname      = 'h3'
        name_len   = 2
      else if( cname(1:name_len) .eq. 'alg6' .or. cname(1:name_len) .eq. 'al-6' ) then
        a_nuc(ii)  = 26.0d0
        z_nuc(ii)  = 13.0d0
        be_nuc(ii) = 8149.771d0 * 1.0d-3 * a_nuc(ii)
        letters    = 'al'
        amass      = '26'
        cname      = 'al26'
        name_len   = 4
      else if( cname(1:name_len) .eq. 'alm6' .or. cname(1:name_len) .eq. 'al*6' ) then
        a_nuc(ii)  = 26.0d0
        z_nuc(ii)  = 13.0d0
        be_nuc(ii) = 8149.771d0 * 1.0d-3 * a_nuc(ii)
        letters    = 'al'
        amass      = '26'
        cname      = 'al26'
        name_len   = 4
      else if ( name_len > 1 ) then
        ! scan name until you find the digit
        digit = .false.
        do jj=1,name_len
          select case( cname(jj:jj) )
          case( '0':'9' )
            digit = .true.
          case default
            digit = .false.
          end select

          if( digit ) then
            letters(1:(jj-1))  = cname(1:(jj-1))    ! pick out letters
            letters            = trim(adjustl(letters))

            amass = cname(jj:name_len) ! pick out digits
            amass = trim(adjustl(amass))

            read( amass, '(i5)' ) ia_nuc      ! assign character mass to real variable            
            exit
          end if
        end do

        ! scan through table until we reach the desired species
        begin = .false.
        do kk = 1,n_masstable
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
          jj = index( binding_r, "#" )
          if( jj /= 0 ) binding_r(jj:jj) = "."
          read( binding_r, "(f11.3)" ) binding

          ! convert name to lower-case
          el = trim(adjustl(el))
          do jj = 1,len(el)
            in_ascii = iachar( el(jj:jj) )
            if( in_ascii >= uc_a_ascii .and. in_ascii <= uc_z_ascii ) then
              el(jj:jj) = achar( in_ascii + (lc_a_ascii - uc_a_ascii) )
            end if
          end do

          ! match names in mass table to network species names
          if( el .eq. letters ) then
            z_nuc(ii)  = dble(zi)
            a_nuc_tmp  = dble(ai)
            be_nuc_tmp = binding*a_nuc_tmp*1.0d-3 ! b.e. from table is [kev / a]
            if( ai .eq. ia_nuc ) then
              a_nuc(ii)  = a_nuc_tmp
              be_nuc(ii) = be_nuc_tmp
              exit
            end if
          end if

          if( kk .eq. n_masstable .and. a_nuc(ii) .ne. dble(ia_nuc) ) then
            a_nuc(ii)  = dble(ia_nuc)
            be_nuc(ii) = be_nuc_tmp
          end if
        end do ! do kk = 1,n_masstable
      end if
    end do  ! do ii = 1,nucnum

    close(lun_table)  ! close mass table file

    return
  end subroutine nucaz_from_name

  subroutine interp1d_mesa( rad_out, state_in, state_out, interp_method )

    use bl_constants_module
    use bl_error_module
    use model_interp_module, only: interp1d_linear, interp1d_spline
    use interpolate_module, only: locate

    ! input variables
    real (dp_t), intent(in) :: rad_out(:)
    real (dp_t), intent(in) :: state_in(:)
    integer, intent(in) :: interp_method

    ! output variables
    real (dp_t), intent(out) :: state_out(size(rad_out))

    ! local variables
    integer :: irad(size(rad_out))
    integer :: irad_max

    real (dp_t) :: dvol1, dvol2
    real (dp_t) :: vol_rad_out(size(rad_out))

    integer :: i, j, n

    vol_rad_out = third * rad_out**3

    if ( interp_method == 1 ) then

      irad_max = nx_mesa_in
      irad(:) = 1
      do i = 1, size(rad_out)
        if ( vol_rad_out(i) <= vol_rad_cntr_mesa_in(1) ) then
          irad(i) = 0
        else if ( vol_rad_out(i) >= vol_rad_cntr_mesa_in(irad_max) ) then
          irad(i) = irad_max
        else
          irad(i) = locate( vol_rad_out(i), irad_max, vol_rad_cntr_mesa_in ) - 1
        end if
      end do
      call interp1d_linear( irad, irad_max, vol_rad_cntr_mesa_in, state_in, vol_rad_out, state_out )
    else if ( interp_method == 2 ) then
      call interp1d_spline( vol_rad_cntr_mesa_in, state_in, vol_rad_out, state_out )
    else
      call bl_error("invalid value for interp_method")
    end if

    return
  end subroutine interp1d_mesa

  subroutine interp2d_mesa( rad_out, state_in, state_out, interp_method )

    use bl_constants_module
    use bl_error_module
    use model_interp_module, only: interp1d_linear, interp1d_spline
    use interpolate_module, only: locate

    ! input variables
    real (dp_t), intent(in) :: rad_out(:,:)
    real (dp_t), intent(in) :: state_in(:)
    integer, intent(in) :: interp_method

    ! output variables
    real (dp_t), intent(out) :: state_out(size(rad_out,1),size(rad_out,2))

    ! local variables
    integer :: irad(size(rad_out,1))
    integer :: irad_max

    real (dp_t) :: dvol1, dvol2
    real (dp_t) :: vol_rad_out(size(rad_out,1),size(rad_out,2))

    integer :: i, j, n

    vol_rad_out = third * rad_out**3

    if ( interp_method == 1 ) then

      irad_max = nx_mesa_in
      do j = 1, size(rad_out,2)
        irad(:) = 1
        do i = 1, size(rad_out,1)
          if ( vol_rad_out(i,j) <= vol_rad_cntr_mesa_in(1) ) then
            irad(i) = 0
          else if ( vol_rad_out(i,j) >= vol_rad_cntr_mesa_in(irad_max) ) then
            irad(i) = irad_max
          else
            irad(i) = locate( vol_rad_out(i,j), irad_max, vol_rad_cntr_mesa_in ) - 1
          end if
        end do
        call interp1d_linear( irad, irad_max, vol_rad_cntr_mesa_in, state_in, vol_rad_out(:,j), state_out(:,j) )
      end do
    else if ( interp_method == 2 ) then
      do j = 1, size(rad_out,2)
        call interp1d_spline( vol_rad_cntr_mesa_in, state_in, vol_rad_out(:,j), state_out(:,j) )
      end do
    else
      call bl_error("invalid value for interp_method")
    end if

    return
  end subroutine interp2d_mesa

end module mesa_parser_module
