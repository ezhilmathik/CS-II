
subroutine read_quakeml

#ifdef NO_FOXY
  print *,"!!! Compiled with -DNO_FOXY - cannot parse quakeML !!!"
  print *,"!!! Stopping here !!!"
  STOP 1
  
#else
    use foxy, only: xml_file
    use parameters
   
    implicit none
    type(xml_file)                :: qml_file  
    type(xml_file)                :: qml_qml   
    type(xml_file)                :: qml_eventParameters 
    type(xml_file)                :: qml_event  
    type(xml_file)                :: qml_origin, qml_magnitude   
    character(len=:), allocatable :: name       
    character(len=32)             :: input_file 
    character(len=:), allocatable :: tag_content

    call qml_file%parse(filename=trim(QuakemlFile))
    
! We have to hop down the levels - I have no idea how to do this in a more direct way...

    call qml_qml%parse(            string=qml_file%content(           'q:quakeml'))
    call qml_eventParameters%parse(string=qml_qml%content(            'eventParameters'))
    call qml_event%parse(          string=qml_eventParameters%content('event'))

    call qml_origin%parse(         string=qml_event%content(          'origin'))
    call qml_magnitude%parse(      string=qml_event%content(          'magnitude'))

! Now, extract the values we need
    call qml_extract_value(qml_origin,'longitude', eq_epi_x)
    call qml_extract_value(qml_origin,'latitude',  eq_epi_y)
    call qml_extract_value(qml_magnitude,'mag',    eq_mag)

    if (verbosity > 0) then
       print *,'QuakeML EQ longitude:', eq_epi_x
       print *,'QuakeML EQ latitude: ', eq_epi_y
       print *,'QuakeML EQ magnitude:', eq_mag
    endif

    ! If the QuakeML file should provide the fault mechanism,
    ! we are not yet able to transfer it into an initial condition.
    !  It is "idealised" (cosine bell) for the time being!

    print *,'QuakeML: only magnitude and epicentre are parsed.'
    print *,'QuakeML: Fault mechanism still to be implemented!'
    print *,'QuakeML: Choosing idealised scenario.'

    enable_idealised = .true.

 contains

  subroutine qml_extract_value(qml, tag, val)
    
    type(xml_file), intent(in)   :: qml
    character(len=*), intent(in) :: tag
    real(kind=wp), intent(out)   :: val

    integer                       :: m, n  
    type(xml_file)                :: qml_val
    character(len=:), allocatable :: content
    integer                       :: ierr

    call qml_val%parse(string=qml%content(trim(adjustl(tag))))
    content = qml_val%content('value')
    m = scan(content,"0.123456789")
    n = scan(content,"0.123456789", back=.true.)
    if (m>0 .and. n>=m) then
       read(content(m:n),*,iostat=ierr) val
    else
       ierr = 1
    endif

    if (ierr /=0 ) then
       print *,"!!! QuakeML: Error in parsing ",trim(adjustl(tag))," !!!"
       print *,"!!! - stopping here !!!"
       STOP 1
    endif

  end subroutine  qml_extract_value
#endif

end subroutine read_quakeml
