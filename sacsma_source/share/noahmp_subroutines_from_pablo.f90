module noahmp_subroutines
! Created by Pablo Mendoza, NCAR (August 22, 2013)
! Based on the work made by Andy Newman for calibrating SAC-SMA using
! SCE-UA (Duan et al., 1993).
! This module contains all the subroutines used to apply the Distributed
! Evaluation of Local Sensitivity Analysis, DELSA (Rakovec et al., 2013)
! using the Noah land surface model with Multiple Parameterization
! options, Noah-MP (Niu et al., 2013)
!
!
!
! Updates:
!
! - Pablo Mendoza, NCAR (September 5, 2013)
!	This new version includes some subroutines that compute scores 
!	based on signature measures.
!	New subroutines: noahmp_scores, fdcsigmeasures, sortpick


  implicit none

!Model setting variables 
  character(len = 1024) :: filelist_name
  character(len = 1024) :: cellfrac_name
  character(len = 1024) :: origsoil_name
  character(len = 1024) :: calibsoil_name
  character(len = 1024) :: origmp_name
  character(len = 1024) :: calibmp_name
  character(len = 1024) :: outputmod_name
  character(len = 1024) :: outputobs_name
  character(len = 1024) :: executable
  integer :: dt
  integer :: sim_len
  integer :: start_cal
  integer :: end_cal
  integer :: nyears
  integer :: start_year
  integer :: opt
  integer :: Ncells
  real 	  :: BasinArea

!DELSA setting variables 
  character(len = 1024) :: baseparams_name
  character(len = 1024) :: parambounds_name
  integer :: Nsamp
  integer :: Npar

  !other variables 
  real, dimension(:), allocatable   :: a      !parameter set
  real, dimension(:), allocatable   :: bl     !lower bounds on parameter set
  real, dimension(:), allocatable   :: bu     !upper bounds on parameter set
  real, dimension(:), allocatable   :: stdev    !standard deviation used for compute DELSA sensitivities
  character(len = 1024) :: sce_fname

!verification variables
  real, dimension(:),   allocatable :: qobs
  real, dimension(:),   allocatable :: qmod

!namelists

!namelist for Noah-MP
!right now read in the following parameters including initial multiplier,
! lower bound value, upper bound value and dimension
  namelist / DELSA_SETTINGS / baseparams_name, parambounds_name, Nsamp, Npar

!control namelist for running models
  namelist / MODEL_SETTINGS /  filelist_name, cellfrac_name, origsoil_name, &
			     calibsoil_name, outputmod_name, origmp_name, &
			     calibmp_name, outputobs_name, executable, dt, &
			     sim_len, start_cal, end_cal, nyears, start_year, &
			     opt, Ncells, BasinArea

  save

  contains


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc                                                               ccc
!cc                      SUBROUTINES!!!!!!!!!                     ccc
!cc 								  ccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccc
subroutine read_namelist
  
  implicit none

!local variables
  integer :: ierr

  open(UNIT=30, file="namelist_noahmp.model",form="FORMATTED")

  read(UNIT=30, NML=MODEL_SETTINGS, iostat=ierr)
  if (ierr /= 0) then
    write(*,'(/," ***** ERROR: Problem reading namelist INIT_CONTROL",/)')
    rewind(UNIT=30)
    read(UNIT=30, NML=MODEL_SETTINGS)
    stop " ***** ERROR: Problem reading namelist INIT_CONTROL"
  endif

  read(UNIT=30, NML=DELSA_SETTINGS, iostat=ierr)
  if (ierr /= 0) then
    write(*,'(/," ***** ERROR: Problem reading namelist NOAHMP",/)')
    rewind(UNIT=30)
    read(UNIT=30, NML=DELSA_SETTINGS)
    stop " ***** ERROR: Problem reading namelist NOAHMP"
  endif

  close(UNIT=30)

end subroutine

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine noahmp_scores(a,rmse,pbiasrr,pbiasfms,pbiasfhv,pbiasflv,pbiasfmm,pbiasctr)
  implicit none

!input variables
  real, dimension(:), intent(in) 	:: a		! Vector with parameter multipliers

!output variables
  real, intent(out)	:: rmse		! Root mean square error
  real, intent(out)	:: pbiasrr	! Percent bias in runoff ratio
  real, intent(out)	:: pbiasfms	! Percent bias in mid-segment slope of FDC
  real, intent(out)	:: pbiasfhv	! Percent bias in high-segment volume of FDC
  real, intent(out)	:: pbiasflv	! Percent bias in low-segment volume of FDC
  real, intent(out)	:: pbiasfmm	! Percent bias in median of streamflow
  real, intent(out)	:: pbiasctr	! Centroid of average water year (starting Oct. 1)

!local variables
  character (len=200)	:: rdline
  character (len=10)	:: auxline
  real    		:: modrealline, auxrow
  integer   		:: i,io,j,k

  print *, 'Initialize score values'

! Initialize score values
  rmse = 0.0
  pbiasrr = 0.0
  pbiasfms = 0.0
  pbiasfhv = 0.0
  pbiasflv = 0.0
  pbiasfmm = 0.0
  pbiasctr = 0.0

! Read original and modified parameter files

! Modify parameters in SOILPARM file
  open (UNIT=54,file=origsoil_name,form='formatted',status='old')
  open (UNIT=55,file=calibsoil_name,form='formatted',status='old')

  do i=1,3
    read(54,'(A)') rdline
    write(55,'(A)') TRIM(rdline)
  enddo
  ! Modify b, maxsmc, satpsi, satdk and qtz
  do j=1,19
    ! b exponent
    read(54,'(A)') rdline
    if (j.LE.9) then
      auxline = rdline(3:11)
      read(auxline,*) auxrow
      modrealline = auxrow*a(1)
      write(auxline,'(f9.2)') modrealline
      rdline(3:11) = auxline
    elseif (j.GE.10) then
      auxline = rdline(4:11)
      read(auxline,*) auxrow
      modrealline = auxrow*a(1)
      write(auxline,'(f8.2)') modrealline
      rdline(4:11) = auxline
    endif
    ! maxsmc
    auxline = rdline(34:41)
    read(auxline,*) auxrow
    modrealline = auxrow*a(2)
    write(auxline,'(f8.3)') modrealline
    rdline(34:41) = auxline
    ! psisat
    auxline = rdline(52:59)
    read(auxline,*) auxrow
    modrealline = auxrow*a(3)
    write(auxline,'(f8.3)') modrealline
    rdline(52:59) = auxline
    ! satdk
    auxline = rdline(61:69)
    read(auxline,*) auxrow
    modrealline = auxrow*a(4)
    write(auxline,'(e9.2)') modrealline
    rdline(61:69) = auxline
    ! quartz
    auxline = rdline(91:96)
    read(auxline,*) auxrow
    modrealline = auxrow*a(5)
    write(auxline,'(f6.2)') modrealline
    rdline(91:96) = auxline
    write(55,'(A)') TRIM(ADJUSTL(rdline))
  enddo
  do i=1,22 
    read(54,'(A)') rdline
    write(55,'(A)') TRIM(rdline)
  enddo

  close(UNIT=54)
  close(UNIT=55)

! Modify parameters in MPTABLE file
  open (UNIT=50,file=origmp_name,form='formatted',status='old')
  open (UNIT=51,file=calibmp_name,form='formatted',status='old')

  do i=1,45
    read(50,'(A)') rdline
    write(51,'(A)') TRIM(rdline)
  enddo
  !z0mvt (Line #45)
  read(50,'(A)') rdline
  do j=1,27
    auxline = rdline((9+(j-1)*7):(14+(j-1)*7))
    read(auxline,*) auxrow
    modrealline = auxrow*a(6)
    write(auxline,'(f6.3)') modrealline
    rdline((9+(j-1)*7):(14+(j-1)*7)) = auxline   
  enddo
  write(51,'(A)') TRIM(rdline)
  do i=1,7
    read(50,'(A)') rdline
    write(51,'(A)') TRIM(rdline)
  enddo
  !rhol (Lines #54 & 55)
  do k = 1,2
    read(50,'(A)') rdline
    do j=1,27
      auxline = rdline((9+(j-1)*7):(14+(j-1)*7))
      read(auxline,*) auxrow
      modrealline = auxrow*a(7)
      write(auxline,'(f6.3)') modrealline
      rdline((9+(j-1)*7):(14+(j-1)*7)) = auxline
    enddo
    write(51,'(A)') TRIM(rdline)
  enddo
  do i=1,3
    read(50,'(A)') rdline
    write(51,'(A)') TRIM(rdline)
  enddo
  !rhos (Lines #59 and 60)
  do k = 1,2
    read(50,'(A)') rdline
    do j=1,27
      auxline = rdline((9+(j-1)*7):(14+(j-1)*7))
      read(auxline,*) auxrow
      modrealline = auxrow*a(8)
      write(auxline,'(f6.3)') modrealline
      rdline((9+(j-1)*7):(14+(j-1)*7)) = auxline
    enddo
    write(51,'(A)') TRIM(rdline)
  enddo
  do i=1,3
    read(50,'(A)') rdline
    write(51,'(A)') TRIM(rdline)
  enddo
  !taul (Lines #64 and 65)
  do k =1,2
    read(50,'(A)') rdline
    do j=1,27
      auxline = rdline((9+(j-1)*7):(14+(j-1)*7))
      read(auxline,*) auxrow
      modrealline = auxrow*a(9)
      write(auxline,'(f6.3)') modrealline
      rdline((9+(j-1)*7):(14+(j-1)*7)) = auxline
    enddo
    write(51,'(A)') TRIM(rdline)
  enddo
  do i=1,3
    read(50,'(A)') rdline
    write(51,'(A)') TRIM(rdline)
  enddo
  !taus (Lines #69 and 70)
  do k = 1,2
    read(50,'(A)') rdline
    do j=1,27
      auxline = rdline((9+(j-1)*7):(14+(j-1)*7))
      read(auxline,*) auxrow
      modrealline = auxrow*a(10)
      write(auxline,'(f6.3)') modrealline
      rdline((9+(j-1)*7):(14+(j-1)*7)) = auxline
    enddo
    write(51,'(A)') TRIM(rdline)
  enddo
  read(50,'(A)') rdline
  write(51,'(A)') TRIM(rdline)
  !xl (Line #72)
  read(50,'(A)') rdline
  do j=1,27
    auxline = rdline((9+(j-1)*7):(14+(j-1)*7))
    read(auxline,*) auxrow
    modrealline = auxrow*a(11)
    write(auxline,'(f6.3)') modrealline
    rdline((9+(j-1)*7):(14+(j-1)*7)) = auxline   
  enddo
  write(51,'(A)') TRIM(rdline)
  read(50,'(A)') rdline
  write(51,'(A)') TRIM(rdline)
  !cwpvt (Line #74)
  read(50,'(A)') rdline
  do j=1,27
    auxline = rdline((9+(j-1)*7):(14+(j-1)*7))
    read(auxline,*) auxrow
    modrealline = auxrow*a(12)
    write(auxline,'(f6.3)') modrealline
    rdline((9+(j-1)*7):(14+(j-1)*7)) = auxline   
  enddo
  write(51,'(A)') TRIM(rdline)
  do i=1,14
    read(50,'(A)') rdline
    write(51,'(A)') TRIM(rdline)
  enddo
  !tmin (Line #89)
  read(50,'(A)') rdline
  do j=1,27
    auxline = rdline((9+(j-1)*7):(14+(j-1)*7))
    read(auxline,*) auxrow
    modrealline = auxrow*a(13)
    write(auxline,'(f6.2)') modrealline
    rdline((9+(j-1)*7):(14+(j-1)*7)) = auxline
  enddo
  write(51,'(A)') TRIM(rdline)
  !vcmx25 (Line #90)
  read(50,'(A)') rdline
  do j=1,27
    auxline = rdline((9+(j-1)*7):(14+(j-1)*7))
    read(auxline,*) auxrow
    modrealline = auxrow*a(14)
    write(auxline,'(f6.2)') modrealline
    rdline((9+(j-1)*7):(14+(j-1)*7)) = auxline
  enddo
  write(51,'(A)') TRIM(rdline)
  read(50,'(A)') rdline
  write(51,'(A)') TRIM(rdline)
  !bp is not modified now! (Line #92)
  read(50,'(A)') rdline
  write(51,'(A)') TRIM(rdline)
  !mp (Line #93)
  read(50,'(A)') rdline
  do j=1,27
    auxline = rdline((9+(j-1)*7):(14+(j-1)*7))
    read(auxline,*) auxrow
    modrealline = auxrow*a(15)
    write(auxline,'(f6.3)') modrealline
    rdline((9+(j-1)*7):(14+(j-1)*7)) = auxline
  enddo
  write(51,'(A)') TRIM(rdline)
  do i=1,10
    read(50,'(A)') rdline
    write(51,'(A)') TRIM(rdline)
  enddo
  !saimfw and saimss (Lines #104 to #115)
  do k=1,3
    read(50,'(A)') rdline
    do j=1,27
      auxline = rdline((9+(j-1)*7):(14+(j-1)*7))
      read(auxline,*) auxrow
      modrealline = auxrow*a(16)	
      write(auxline,'(f6.3)') modrealline
      rdline((9+(j-1)*7):(14+(j-1)*7)) = auxline
    enddo
    write(51,'(A)') TRIM(rdline)
  enddo
  do k=1,6
    read(50,'(A)') rdline
    do j=1,27
      auxline = rdline((9+(j-1)*7):(14+(j-1)*7))
      read(auxline,*) auxrow
      modrealline = auxrow*a(17)	
      write(auxline,'(f6.3)') modrealline
      rdline((9+(j-1)*7):(14+(j-1)*7)) = auxline
    enddo
    write(51,'(A)') TRIM(rdline)
  enddo
  do k=1,3
    read(50,'(A)') rdline
    do j=1,27
      auxline = rdline((9+(j-1)*7):(14+(j-1)*7))
      read(auxline,*) auxrow
      modrealline = auxrow*a(16)	
      write(auxline,'(f6.3)') modrealline
      rdline((9+(j-1)*7):(14+(j-1)*7)) = auxline
    enddo
    write(51,'(A)') TRIM(rdline)
  enddo

  read(50,'(A)') rdline
  write(51,'(A)') TRIM(rdline)

  !laim (Lines #117 to #128)
  do k=1,3
    read(50,'(A)') rdline
    do j=1,27
      auxline = rdline((9+(j-1)*7):(14+(j-1)*7))
      read(auxline,*) auxrow
      modrealline = auxrow*a(18)	
      write(auxline,'(f6.3)') modrealline
      rdline((9+(j-1)*7):(14+(j-1)*7)) = auxline
    enddo
    write(51,'(A)') TRIM(rdline)
  enddo
  do k=1,6
    read(50,'(A)') rdline
    do j=1,27
      auxline = rdline((9+(j-1)*7):(14+(j-1)*7))
      read(auxline,*) auxrow
      modrealline = auxrow*a(19)	
      write(auxline,'(f6.3)') modrealline
      rdline((9+(j-1)*7):(14+(j-1)*7)) = auxline
    enddo
    write(51,'(A)') TRIM(rdline)
  enddo
  do k=1,3
    read(50,'(A)') rdline
    do j=1,27
      auxline = rdline((9+(j-1)*7):(14+(j-1)*7))
      read(auxline,*) auxrow
      modrealline = auxrow*a(18)	
      write(auxline,'(f6.3)') modrealline
      rdline((9+(j-1)*7):(14+(j-1)*7)) = auxline
    enddo
    write(51,'(A)') TRIM(rdline)
  enddo

  do i=1,11
    read(50,'(A)') rdline
    write(51,'(A)') TRIM(rdline)
  enddo

  !fff (TEST1,Line #140)
  read(50,'(A)') rdline
  do j=1,27
    auxline = rdline((9+(j-1)*7):(14+(j-1)*7))
    read(auxline,*) auxrow
    modrealline = auxrow*a(20)
    write(auxline,'(f6.4)') modrealline
    rdline((9+(j-1)*7):(14+(j-1)*7)) = auxline
  enddo
  write(51,'(A)') TRIM(rdline)
  !rsbmx (TEST2,Line #141)
  read(50,'(A)') rdline
  do j=1,27
    auxline = rdline((9+(j-1)*7):(14+(j-1)*7))
    read(auxline,*) auxrow
    modrealline = auxrow*a(21)
    write(auxline,'(f6.4)') modrealline
    rdline((9+(j-1)*7):(14+(j-1)*7)) = auxline
  enddo
  write(51,'(A)') TRIM(rdline)
  !timean (TEST3,Line #142)
  read(50,'(A)') rdline
  do j=1,27
    auxline = rdline((9+(j-1)*7):(14+(j-1)*7))
    read(auxline,*) auxrow
    modrealline = auxrow*a(22)
    write(auxline,'(f6.3)') modrealline
    rdline((9+(j-1)*7):(14+(j-1)*7)) = auxline
  enddo
  write(51,'(A)') TRIM(rdline)
  !fsatmx (TEST4,Line #143)
  read(50,'(A)') rdline
  do j=1,27
    auxline = rdline((9+(j-1)*7):(14+(j-1)*7))
    read(auxline,*) auxrow
    modrealline = auxrow*a(23)
    write(auxline,'(f6.4)') modrealline
    rdline((9+(j-1)*7):(14+(j-1)*7)) = auxline
  enddo
  write(51,'(A)') TRIM(rdline)
  !mexp (TEST5,Line #144)
  read(50,'(A)') rdline
  do j=1,27
    auxline = rdline((9+(j-1)*7):(14+(j-1)*7))
    read(auxline,*) auxrow
    modrealline = auxrow*a(24)
    write(auxline,'(f6.3)') modrealline
    rdline((9+(j-1)*7):(14+(j-1)*7)) = auxline
  enddo
  write(51,'(A)') TRIM(rdline)
  !zosno (TEST6,Line #145)
  read(50,'(A)') rdline
  do j=1,27
    auxline = rdline((9+(j-1)*7):(14+(j-1)*7))
    read(auxline,*) auxrow
    modrealline = auxrow*a(25)
    write(auxline,'(f6.3)') modrealline
    rdline((9+(j-1)*7):(14+(j-1)*7)) = auxline
  enddo
  write(51,'(A)') TRIM(rdline)
  !ssi (TEST7,Line #146)
  read(50,'(A)') rdline
  do j=1,27
    auxline = rdline((9+(j-1)*7):(14+(j-1)*7))
    read(auxline,*) auxrow
    modrealline = auxrow*a(26)
    write(auxline,'(f6.4)') modrealline
    rdline((9+(j-1)*7):(14+(j-1)*7)) = auxline
  enddo
  write(51,'(A)') TRIM(rdline)
  !swemx (TEST8,Line #147)
  read(50,'(A)') rdline
  do j=1,27
    auxline = rdline((9+(j-1)*7):(14+(j-1)*7))
    read(auxline,*) auxrow
    modrealline = auxrow*a(27)
    write(auxline,'(f6.3)') modrealline
    rdline((9+(j-1)*7):(14+(j-1)*7)) = auxline
  enddo
  write(51,'(A)') TRIM(rdline)
  !albmin (TEST9,Line #148)
  read(50,'(A)') rdline
  do j=1,27
    auxline = rdline((9+(j-1)*7):(14+(j-1)*7))
    read(auxline,*) auxrow
    modrealline = auxrow*a(28)
    write(auxline,'(f6.4)') modrealline
    rdline((9+(j-1)*7):(14+(j-1)*7)) = auxline
  enddo
  write(51,'(A)') TRIM(rdline)
  !albmax (TEST10,Line #149)
  read(50,'(A)') rdline
  do j=1,27
    auxline = rdline((9+(j-1)*7):(14+(j-1)*7))
    read(auxline,*) auxrow
    modrealline = auxrow*a(29)
    write(auxline,'(f6.4)') modrealline
    rdline((9+(j-1)*7):(14+(j-1)*7)) = auxline
  enddo
  write(51,'(A)') TRIM(rdline)
  !albdec (TEST11,Line #150)
  read(50,'(A)') rdline
  do j=1,27
    auxline = rdline((9+(j-1)*7):(14+(j-1)*7))
    read(auxline,*) auxrow
    modrealline = auxrow*a(30)
    write(auxline,'(f6.3)') modrealline
    rdline((9+(j-1)*7):(14+(j-1)*7)) = auxline
  enddo
  write(51,'(A)') TRIM(rdline)

  do i=1,160
    read(50,'(A)') rdline
    write(51,'(A)') TRIM(rdline)
  enddo

  close(UNIT=50)
  close(UNIT=51)

  print *, 'model parameters have been modified'

  !call Noah-MP model
  call system(executable) !< run Noah

  print *, 'model has been run'

  call read_modflow(outputmod_name,cellfrac_name,filelist_name,sim_len,qmod)
    !    !read model output

  !read observations
  call read_obsflow(outputobs_name,sim_len,qobs)
    !    !read streamflow observations

   ! Compute RMSE
   call calc_rmse(qmod,qobs,rmse)
   ! Compute signature measures based on the FDC
   call fdcsigmeasures(qmod,qobs,pbiasrr,pbiasfms,pbiasfhv,pbiasflv,pbiasfmm)
   ! Compute percent bias in centroid of hydrograph
   call pbiascentroid(qmod,qobs,pbiasctr)

  return

end subroutine

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine calc_rmse(model,streamflow,rmse)
  implicit none

!input variables (model: simulations, streamflow: observations)
  real, dimension(:),     intent(in)  :: model
  real, dimension(:),     intent(in)  :: streamflow

!output variables
  real, intent(out)	:: rmse	! root mean square error

!local variables
  integer :: itime
  real :: sum_sqr
  

  sum_sqr = 0.0

  do itime = start_cal,end_cal
    sum_sqr = sum_sqr + (model(itime)-streamflow(itime))**2
  enddo
  rmse = sqrt(sum_sqr/real(end_cal - start_cal + 1))
!  print *,'rmse: ', model(1),streamflow(1),rmse,dt
  return
end subroutine

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine fdcsigmeasures(qmod,qobs,pbiasrr,pbiasfms,pbiasfhv,pbiasflv,pbiasfmm)
  implicit none
! This subroutine computes model scores based on signature measures of catchment
! behavior (Yilmaz et al., 2008)

!input variables (model: simulations, streamflow: observations)
  real, dimension(:),     intent(in)  :: qmod	! Modeled streamflow
  real, dimension(:),     intent(in)  :: qobs	! Simulated streamflow

!output variables
  real, intent(out)	:: pbiasrr	! Percent bias in runoff ratio
  real, intent(out)	:: pbiasfms	! Percent bias in mid-segment slope of FDC
  real, intent(out)	:: pbiasfhv	! Percent bias in high-segment volume of FDC
  real, intent(out)	:: pbiasflv	! Percent bias in low-segment volume of FDC
  real, intent(out)	:: pbiasfmm	! Percent bias in median of streamflow

!local variables
  integer :: itime
  integer :: fhvindex   ! Maximum index that contains the high flows portion in FDC
  integer :: flvindex   ! Minimum index that contains the low flows portion in FDC
  integer :: fmsmin     ! Minimum index that contains the mid-segment portion in FDC
  integer :: fmsmax     ! Maximum index that contains the mid-segment portion in FDC
  integer :: Ntime
  real, dimension(end_cal-start_cal+1)	:: probexc	! exceedance probabilities
  real, dimension(end_cal-start_cal+1)	:: modfdc	! modeled FDC
  real, dimension(end_cal-start_cal+1)	:: obsfdc	! observed FDC
  real			:: modfdcmed	! modeled median
  real			:: obsfdcmed	! observed median
  real			:: modslopemin	! modeled flow for Pexc = 0.2
  real			:: modslopemax	! modeled flow for Pexc = 0.7
  real			:: obsslopemin	! observed flow for Pexc = 0.2
  real			:: obsslopemax	! observed flow for Pexc = 0.7
  real			:: aux1,aux2	! ancillary real values

! Compute total length time
  Ntime = end_cal - start_cal + 1

 ! Compute the probabilities of exceedance and save relevant index thresholds
  do itime = 1,Ntime
    probexc(itime) = real(itime)/real(Ntime)
    if ( probexc(itime).LT.0.02 ) then
      fhvindex = itime ! Maximum index for which Pexc < 0.02
    endif
    if ( probexc(itime).LT.0.2 ) then
      fmsmin = itime ! Maximum index for which Pexc < 0.2
    endif
    if ( probexc(itime).LT.0.7 ) then
      flvindex = itime ! Maximum index for which Pexc < 0.7
    endif
  enddo

  fmsmax = flvindex
  flvindex = flvindex + 1 ! We want the minimum index for which Pexc is more than 0.7
  
! Obtain the flow arrays associated with these probabilities
  call sortpick(qmod,modfdc,modfdcmed)
  call sortpick(qobs,obsfdc,obsfdcmed)

! Print the median of observed and modeled streamflow
  print *, 'The median of modeled streamflow is'
  print *, modfdcmed
  print *, 'The median of observed streamflow is'
  print *, obsfdcmed
  print *, 'The size of the sorted vectors is'
  print *, size(modfdc)
 
! Compute the percent bias in runoff ratio
  pbiasrr = ( sum(modfdc) - sum(obsfdc) ) / sum(obsfdc) * 100.0

! Compute the percent bias in mid-segment slope of FDC
  modslopemin = modfdc(fmsmin) - (modfdc(fmsmin) - modfdc(fmsmin+1)) / ( probexc(fmsmin) - probexc(fmsmin+1) ) * ( probexc(fmsmin) - 0.2)
  obsslopemin = obsfdc(fmsmin) - (obsfdc(fmsmin) - obsfdc(fmsmin+1)) / ( probexc(fmsmin) - probexc(fmsmin+1) ) * ( probexc(fmsmin) - 0.2)
  modslopemax = modfdc(fmsmax) - (modfdc(fmsmax) - modfdc(fmsmax+1)) / ( probexc(fmsmax) - probexc(fmsmax+1) ) * ( probexc(fmsmax) - 0.7)
  obsslopemax = obsfdc(fmsmax) - (obsfdc(fmsmax) - obsfdc(fmsmax+1)) / ( probexc(fmsmax) - probexc(fmsmax+1) ) * ( probexc(fmsmax) - 0.7)
! Print the 20 % and 70 % probability exceedance flows
  print *, 'Observed 20% daily flow'
  print *, obsslopemin
  print *, 'Observed 70% daily flow'
  print *, obsslopemax
  print *, 'Modeled 20% daily flow'
  print *, modslopemin
  print *, 'Modeled 70% daily flow'
  print *, modslopemax
  pbiasfms = ( (log10(modslopemin) - log10(modslopemax)) - (log10(obsslopemin) - log10(obsslopemax)) )/ (log10(obsslopemin) - log10(obsslopemax))*100

! Compute the percent bias in high flow volumes
  aux1 =  sum(log10(modfdc(flvindex:Ntime))) - real(Ntime-flvindex+1)*log10(modfdc(Ntime))
  aux2 =  sum(log10(obsfdc(flvindex:Ntime))) - real(Ntime-flvindex+1)*log10(obsfdc(Ntime)) 
  pbiasflv = -( aux1 - aux2 ) / aux2 * 100.0

! Compute the percent bias in low flow volumes
  pbiasfhv =  ( sum(modfdc(1:fhvindex)) - sum(obsfdc(1:fhvindex)) ) / sum(obsfdc(1:fhvindex)) * 100.0

! Compute the percent bias in median
  pbiasfmm = ( log10(modfdcmed) - log10(obsfdcmed) ) / log10(obsfdcmed) * 100.0


  return
end subroutine

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine pbiascentroid(qmod,qobs,pbiasctr)
  implicit none

!input variables
  real, dimension(:), intent(in)  :: qmod	! Modeled streamflow
  real, dimension(:), intent(in)  :: qobs	! Simulated streamflow

!output variables
  real, intent(out)	:: pbiasctr	! Percent bias in centroid of signature measures

!local variables
  real			:: modcentroid	! centroid in modeled time series
  real			:: obscentroid	! centroid in observed time series

! Compute the average hydrograph centroid for both observed and modeled time series
  call avgcentroid(qmod,modcentroid)
  call avgcentroid(qobs,obscentroid)

! Compute the percent bias
  pbiasctr = (modcentroid - obscentroid)/obscentroid * 100.0

  return
end subroutine

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine avgcentroid(tseries,centroid)
  implicit none
! This subroutine computes the average centroid of a water year. A key assumption is that
! the first value in the time series provided starts on October 1, and the last value is
! for September 30.

!input variables
  real, dimension(:), intent(in)	:: tseries	! time series with the variable of interest

!output variables
  real, intent(out)			:: centroid	! mean centroid of annual hydrograph (days)

!local variables
  integer			:: iyear	! index for each water year
  integer			:: iday	! index for each day
  integer			:: year		! current year
  integer			:: startindex	! index for the beginning of each water year
  integer			:: endindex	! index for the end of each water year
  real,dimension(nyears)  	:: ctallyears   ! centroid for all water years within the time series
  real				:: sumctr	! ancillary variable

  startindex = start_cal 	! Initialize the starting year index
  ctallyears = 0.0 	! Initialize array with centroids of all water years
  centroid = 0.0

! Start loop over water years
  do iyear = 1,nyears
      
    year = start_year + iyear - 1
    sumctr = 0.0

    if ( MOD(year+1,4).EQ.0 ) then
      endindex = startindex + 366 - 1
      do iday=startindex,endindex
	sumctr = sumctr + tseries(iday)*real(iday-startindex+1)
      enddo
      ctallyears(iyear) = sumctr/sum(tseries(startindex:endindex))
      startindex = startindex + 366
    endif
    if ( MOD(year+1,4).GT.0 ) then
      endindex = startindex + 365 - 1
      do iday=startindex,endindex
	sumctr = sumctr + tseries(iday)*real(iday-startindex+1)
      enddo
      ctallyears(iyear) = sumctr/sum(tseries(startindex:endindex))
      startindex = startindex + 365
    endif

    centroid = centroid + ctallyears(iyear)/real(nyears)

  enddo ! End loop over water years



  return
end subroutine


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine sortpick(vectorin,vectorout,median)
  implicit none
! Sorts an array vector into descending numerical order, by straight insertion. 
! vectorin is replaced on output by its sorted rearrangement (vectorout)
! Additionally, the median of the array is also computed

! input variables
  real, dimension(:), intent(in)	:: vectorin 	! input vector

! output variables
  real, dimension(:), intent(out)	:: vectorout 	! output vector
  real, intent(out)			:: median 	! output mean

! local variables
    integer 	:: i,j,n
    real 	:: a

  vectorout = vectorin(start_cal:end_cal)  ! Extract input data for the time window required

  n = size(vectorout)
  print *, 'the size of the output vector (vectorout) is'
  print *, n
  
  do j=2,n ! Pick out each element in turn.
    a=vectorout(j)
    do i=j-1,1,-1 ! Look for the place to insert it.
      if (vectorout(i) >= a) exit
      vectorout(i+1)=vectorout(i)
    enddo
    vectorout(i+1)=a !Insert it.
  enddo


  if ( mod(n,2) == 0 ) then
     median = (vectorout(n/2+1) + vectorout(n/2))/2.0
  else
     median = vectorout(n/2+1)
  endif

  return
end subroutine


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine read_obsflow(stream_name,sim_len,qobs)
  implicit none
! This subroutine reads observed streamflow

!input variables
  character(len = 1024),	intent(in)  :: stream_name
  integer,			intent(in)  :: sim_len

!output variables
  real, dimension(:), intent(out) :: qobs

!local variables
  integer :: itime

!read observed streamflow
  open (UNIT=60,file=stream_name,form='formatted',status='old')

  do itime = 1,sim_len
    read (UNIT=60,*) qobs(itime)
  enddo
  close(UNIT=60)

  return
end subroutine

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine read_modflow(outputdir,cellfile,listfile,sim_len,qmod)
  use netcdf
  implicit none
!! The structure of this subroutine should be very similar to that for VIC
!input variables
  character(len = 200),	intent(in)  :: outputdir
  character(len = 200),	intent(in)  :: cellfile
  character(len = 200),	intent(in)  :: listfile
  integer,	intent(in)  :: sim_len

!output variables
  real, dimension(:), intent(out) :: qmod

!local variables
  character(len = 1000) 	:: filename
  real 				:: cellfraction
  real, dimension(sim_len)	:: qsurf, qb
  integer :: icell, itime, ncid, varid

!! Initialize qmod
  qmod = 0.0

!read observed streamflow
  open (UNIT=61,file=listfile,form='formatted',status='old')
  open (UNIT=62,file=cellfile,form='formatted',status='old')

  do icell = 1,Ncells

    read (UNIT=61,*) filename
    read (UNIT=62,*) cellfraction
    filename = TRIM(outputdir) // TRIM(filename)
    
    ! Open netcdf file for reading
    call check( nf90_open( TRIM(ADJUSTL(filename)), NF90_NOWRITE, ncid) )
    !nf90_open(filename, nf90_NoWrite, ncid) 

    ! Get surface runoff
    call check( nf90_inq_varid(ncid, "RUNSRF", varid) )
    call check( nf90_get_var(ncid, varid, qsurf) )
    !nf90_inq_varid(ncid, "SFRUNOFF", varid)
    !nf90_get_var(ncid, varid, qsurf)

    ! Get baseflow
    call check( nf90_inq_varid(ncid, "RUNSUB", varid) )
    call check( nf90_get_var(ncid, varid, qb) )
    !nf90_inq_varid(ncid, "UDRUNOFF", varid)
    !nf90_get_var(ncid, varid, qb)

    do itime = 1,sim_len

      ! The next line is the sum of baseflow + surface runoff
      qmod(itime) = qmod(itime) + (qsurf(itime) + qb(itime))*cellfraction*BasinArea*( 10**3 )      

    enddo

    ! Close the output file for grid cell icell
    call check( nf90_close(ncid) )

  enddo

  close(UNIT=61)
  close(UNIT=62)

  return
end subroutine

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine check(status)
  use netcdf
  integer	:: status
  if(status /= nf90_noerr) then
    !call handle_err(status)
    print *, trim(nf90_strerror(status))
    stop
  end if
end subroutine check  

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine sce_param_setup(a,bl,bu,stdev)
  implicit none
!output variables
  real,dimension(Npar),intent(out) :: a		! Array with parameter multipliers
  real,dimension(Npar),intent(out) :: bl	! Lower bounds for parameter multipliers
  real,dimension(Npar),intent(out) :: bu	! Upper bounds for parameter multipliers
  real,dimension(Npar),intent(out) :: stdev	! standard deviation
!local variables
  integer		:: ipar 	! index for parameters
  character (len=200) 	:: rdline
  character (len=10) 	:: auxline
  real    		:: auxrow

! Array with parameter multipliers is initialized with 1
  a = 1.0

! Read parameter bounds
  open (UNIT=70,file=parambounds_name,form='formatted',status='old')

! Read the first line
  read(70,'(A)') rdline

  do ipar = 1,Npar ! Start loop over parameters

    read(70,'(A)') rdline
    ! Extract lower bound
    auxline = rdline(10:18)
    read(auxline,*) auxrow
    bl(ipar) = auxrow
    ! Extract upper bound
    auxline = rdline(19:26)
    read(auxline,*) auxrow
    bu(ipar) = auxrow
    ! Compute standard deviation using formula for uniform distribution
    stdev(ipar) = sqrt ( 1.0 / 12.0 * (bu(ipar) - bl(ipar)) ** 2.0 )

  enddo ! End loop over parameters


  close(UNIT=70)

  return
end subroutine


end module