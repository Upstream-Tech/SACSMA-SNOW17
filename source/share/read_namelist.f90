subroutine read_namelist(namelist_name)
  use nrtype
  use snow17_sac

  implicit none
!input variable
  character(len=2000),intent(in)	:: namelist_name

!local variables
  integer(I4B) :: ierr

  open(UNIT=30, file=trim(namelist_name),form="FORMATTED")

  read(UNIT=30, NML=INIT_CONTROL, iostat=ierr)
  if (ierr /= 0) then
    write(*,'(/," ***** ERROR: Problem reading namelist INIT_CONTROL",/)')
    rewind(UNIT=30)
    read(UNIT=30, NML=INIT_CONTROL)
    stop " ***** ERROR: Problem reading namelist INIT_CONTROL"
  endif

  close(UNIT=30)

  return
end subroutine