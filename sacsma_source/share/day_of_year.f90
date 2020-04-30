subroutine julian_day(year,month,day,jday)
  use nrtype

  implicit none
!
! Taken from Glen Liston's SnowModel code, updated to F90
!


!input variables
  integer(I4B),dimension(:),intent(in) :: year
  integer(I4B),dimension(:),intent(in) :: month
  integer(I4B),dimension(:),intent(in) :: day

!output variables
  integer(I4B),dimension(:),intent(out) :: jday

!local variables
  integer(I4B)      :: i


! Calculate the day of year (1,...,365,366) corresponding to the date
!   iyear-imonth-iday. 
  jday = day &
          + min(1,max(0,month-1))*31 &
          + min(1,max(0,month-2))*(28+(1-min(1,mod(year,4)))) &
          + min(1,max(0,month-3))*31 &
          + min(1,max(0,month-4))*30 &
          + min(1,max(0,month-5))*31 &
          + min(1,max(0,month-6))*30 &
          + min(1,max(0,month-7))*31 &
          + min(1,max(0,month-8))*31 &
          + min(1,max(0,month-9))*30 &
          + min(1,max(0,month-10))*31 &
          + min(1,max(0,month-11))*30 &
          + min(1,max(0,month-12))*31

  return
end subroutine

!ccccccccccccccccccccccccccccccccc

subroutine julianday_scalar(iyear,imonth,iday,jday_scalar)
  use nrtype

  implicit none

!
!  Based on above subroutine from Glen Liston's SnowModel code
!

!input variables
  integer(I4B), intent(in) :: iyear
  integer(I4B), intent(in) :: imonth
  integer(I4B), intent(in) :: iday

!output variables
  integer(I4B),intent(out) :: jday_scalar



! Calculate the day of year (1...365,366) corresponding to the date

  jday_scalar = iday &
          + min(1,max(0,imonth-1))*31 &
          + min(1,max(0,imonth-2))*(28+(1-min(1,mod(iyear,4)))) &
!	  + min(1,max(0,imonth-2))*28 &
          + min(1,max(0,imonth-3))*31 &
          + min(1,max(0,imonth-4))*30 &
          + min(1,max(0,imonth-5))*31 &
          + min(1,max(0,imonth-6))*30 &
          + min(1,max(0,imonth-7))*31 &
          + min(1,max(0,imonth-8))*31 &
          + min(1,max(0,imonth-9))*30 &
          + min(1,max(0,imonth-10))*31 &
          + min(1,max(0,imonth-11))*30 &
          + min(1,max(0,imonth-12))*31

  return
end subroutine
