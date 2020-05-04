subroutine calc_rmse(model,obs,length,valid,rmse)
  use nrtype

  implicit none

!input variables
  real(dp), dimension(36500),     intent(in)  :: model
  real(dp), dimension(36500),     intent(in)  :: obs
  
  integer(I4B),	intent(in)			:: length
  logical, dimension(36500),intent(in)		:: valid

!output variables
  real(dp),		  intent(out) :: rmse

!local variables
  integer(I4B) :: i
  integer(I4B) :: num_valid
  integer(I4B) :: end_pt
  real(dp) :: sum_sqr


!!! code

!  end_pt = size(model)

  end_pt = length

  num_valid = count(valid(1:end_pt))

!print *,'num missing',end_pt-num_valid,num_valid
 
  sum_sqr = 0.0_dp

!print *,'rmse size arrays:',end_pt

!  do i = 1,end_pt
!    sum_sqr = sum_sqr + (model(i)-obs(i))**2
!  enddo

!  rmse = sqrt(sum_sqr/(real(end_pt)))

  sum_sqr = sum((model(1:end_pt)-obs(1:end_pt))**2,MASK=valid)
  rmse = sqrt(sum_sqr/real(num_valid))

!print *,num_valid,sum_sqr,rmse

  return
end subroutine calc_rmse

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calc_mse(model,obs,length,valid,mse)
  use nrtype

  implicit none

!input variables
  real(dp), dimension(36500),     intent(in)  :: model
  real(dp), dimension(36500),     intent(in)  :: obs

  integer(I4B),	intent(in)			:: length
  logical, dimension(36500),intent(in)		:: valid


!output variables
  real(dp),		  intent(out) :: mse

!local variables
  integer(I4B) :: i
  integer(I4B) :: end_pt
  integer(I4B) :: num_valid
  real(dp) :: sum_sqr


!!! code

  sum_sqr = 0.0_dp

  end_pt = length

  num_valid = count(valid(1:end_pt))

!  end_pt = size(model)
!  end_pt = length

!  do i = 1,end_pt
!    sum_sqr = sum_sqr + (model(i)-obs(i))**2
!  enddo
!  mse = sum_sqr/(real(end_pt))


  sum_sqr = sum((model(1:end_pt)-obs(1:end_pt))**2,MASK=valid)
  mse = sum_sqr/real(num_valid)

  return
end subroutine calc_mse

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calc_nse(model,obs,length,valid,nse)
  use nrtype
  use snow17_sac, only: mean_obs

  implicit none

!input variables
  real(dp), dimension(36500),     intent(in)  :: model
  real(dp), dimension(36500),     intent(in)  :: obs

  integer(I4B),	intent(in)			:: length
  logical, dimension(36500), intent(in)	:: valid

!output variables
  real(dp),		  intent(out) :: nse

!local variables
  integer(I4B) :: i
  integer(I4B) :: end_pt
  integer(I4B) :: num_valid

  real(dp) :: sum_sqr
  real(dp) :: sum_obs


!!! code

  sum_sqr = 0.0_dp
  sum_obs = 0.0_dp

!  end_pt = size(model)
  end_pt = length

  num_valid = count(valid(1:end_pt))

!print *,'num missing',end_pt,num_valid

!  do i = 1,end_pt
!    sum_sqr = sum_sqr + (model(i)-obs(i))**2
!    sum_obs = sum_obs + (obs(i)-mean_obs)**2
!  enddo
!  nse = -1.0_dp*(1.0 - sum_sqr/sum_obs)  !inverted sign for optimization

  sum_sqr = sum((model(1:end_pt)-obs(1:end_pt))**2,valid)
  sum_obs = sum((obs(1:end_pt)-mean_obs)**2,valid)

  nse = -1.0_dp*(1.0_dp - sum_sqr/sum_obs)

  return
end subroutine calc_nse

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calc_kge(model,obs,length,valid,kge)
  use nrtype
  use gauge_calib, only: pearson
  implicit none

! This subroutine calculates the negative of the Kling-Gupta Efficiency
! in order to solve a minimization problem
! from Pablo Mendoza
! modified by Andy Newman


!input variables
  real(dp), dimension(36500),     intent(in)  :: model
  real(dp), dimension(36500),     intent(in)  :: obs

  integer(I4B),	intent(in)			:: length
  logical, dimension(36500), intent(in)		:: valid

!output variables
  real(dp),     intent(out) :: kge

!local variables
  integer(I4B) :: itime
  real(dp) :: cc,alpha,betha,mu_s,mu_o,sigma_s,sigma_o
  integer(I4B) :: end_pt
  integer(I4B)	:: num_valid

!code
!  end_pt = size(model)
  end_pt = length
  num_valid = count(valid)

  mu_s = 0.0_dp
  mu_o = 0.0_dp
  sigma_s = 0.0_dp
  sigma_o = 0.0_dp
  !We first compute the mean
  mu_s = sum(model(1:end_pt),valid)/real(num_valid,kind(dp))
  mu_o = sum(obs(1:end_pt),valid)/real(num_valid,kind(dp))
  betha = mu_s/mu_o
  !Now we compute the standard deviation
!  do itime = 1,end_pt
!    sigma_s = sigma_s + (model(itime)-mu_s)**2
!    sigma_o = sigma_o + (obs(itime)-mu_o)**2
!  enddo

  sigma_s = sum((model(1:end_pt)-mu_s)**2,valid)
  sigma_o = sum((obs(1:end_pt)-mu_o)**2,valid)

  sigma_s = sqrt(mu_s/real(num_valid,kind(dp)))
  sigma_o = sqrt(mu_o/real(num_valid,kind(dp)))
  alpha = sigma_s/sigma_o

  !Compute linear correlation coefficient
  call pearson(model,obs,length,valid,cc)

!  print *,'mu',mu_s,mu_o
!  print *,'sigma',sigma_s,sigma_o
!  print *,'ab',alpha,betha
!  print *,'rho',cc

  !inverted sign for minimization
  kge = -( 1.0 - sqrt((cc-1.0)**2 + (alpha-1.0)**2 + (betha-1.0)**2) )

  return
end subroutine calc_kge

subroutine pearson(model,obs,length,valid,corr)
  use nrtype


  implicit none

  interface
    logical function isnan(a)
      use nrtype
      real(dp),intent(in) :: a
    end function isnan
  end interface

  !input variables
  real(dp), dimension(36500),intent(in)	:: model
  real(dp), dimension(36500),intent(in)	:: obs

  integer(I4B),	intent(in)		:: length
  logical, dimension(36500), intent(in)	:: valid

  !output variables
  real(dp), intent(out)			:: corr

  !local variables
  real(dp)	:: model_mean
  real(dp)	:: obs_mean
  real(dp)	:: model_var
  real(dp)	:: obs_var
  real(dp)	:: cov

  integer(dp)	:: end_pt
  integer(i4b)	:: num_valid


  !code

!  end_pt = size(model)
  end_pt = length
  num_valid = count(valid)

!compute means
  model_mean = sum(model(1:end_pt),valid)/real(num_valid,kind(dp))
  obs_mean   = sum(obs(1:end_pt),valid)/real(num_valid,kind(dp))

!compute variance,covariance
  model_var = sum((model(1:end_pt) - model_mean)**2,valid)
  obs_var   = sum((obs(1:end_pt)   - obs_mean)**2,valid)
  cov       = sum((model(1:end_pt) - model_mean)*(obs(1:end_pt) - obs_mean),valid)


!compute correlation
  corr = cov/(sqrt(model_var)*sqrt(obs_var))

  if(isnan(corr)) then
    corr = 0.0_dp
  endif

!  print *,'p',sqrt(model_var),sqrt(obs_var),cov,corr

  return
end subroutine pearson


logical function isnan(a)
  use nrtype
  implicit none

  real(dp),intent(in) :: a

  if (a.ne.a) then
    isnan = .true.
  else
    isnan = .false.
  end if

  return
end 


subroutine fdc_fms(model,obsfdc,order,length,valid,pbiasfms)
  !modified from pablo mendoza's subroutine following yilmax et al. (2008)
  !only computes bias for mid-segment of fdc

  use nrtype
  use nr, only: sort_heap

  implicit none


!input variables
  real(DP), dimension(36500),     intent(in)  :: model	! Modeled streamflow
  real(DP), dimension(36500),     intent(in)  :: obsfdc	! sorted observed streamflow 
  integer(I4B),dimension(36500),  intent(in)  :: order	! order of original observed streamflow

  integer(I4B),	intent(in)			:: length
  logical, dimension(36500),intent(in)		:: valid

!output variables
  real(DP), intent(out) :: pbiasfms

!local variables
  integer(I4B) :: itime,num_valid,cnt,i
  integer(I4B) :: fhvindex   ! Maximum index that contains the high flows portion in FDC
  integer(I4B) :: flvindex   ! Minimum index that contains the low flows portion in FDC
  integer(I4B) :: fmsmin     ! Minimum index that contains the mid-segment portion in FDC
  integer(I4B) :: fmsmax     ! Maximum index that contains the mid-segment portion in FDC
  integer(I4B) :: end_pt
  real(DP), dimension(36500)	:: probexc	! exceedance probabilities
  real(SP), dimension(36500)	:: modfdc_a	! modeled FDC
  real(SP), dimension(36500)	:: modfdc	! modeled FDC
  real(SP), dimension(36500)	:: obsfdc_sp	! obseved fdc single precision
  real(DP)			:: modslopemin	! modeled flow for Pexc = 0.2
  real(DP)			:: modslopemax	! modeled flow for Pexc = 0.7
  real(DP)			:: obsslopemin	! observed flow for Pexc = 0.2
  real(DP)			:: obsslopemax	! observed flow for Pexc = 0.7


  end_pt = length

  num_valid = count(valid(1:end_pt))

 ! Compute the probabilities of exceedance and save relevant index thresholds
  do itime = 1,num_valid
    probexc(itime) = real(itime)/real(num_valid)
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
!  call sortpick(model,modfdc,modfdcmed)
!  call sortpick(obs,obsfdc,obsfdcmed)

!print *,num_valid

  do i = 1,num_valid
    if(obsfdc(i) .gt. -0.00000001) then
      modfdc_a(i) = real(model(order(i)),kind(SP))
!print *,model(order(i)),modfdc(i)
    endif
  enddo

!   sort modelfdc into descending order
  call sort_heap(modfdc_a(1:num_valid))
  do i = 1,num_valid
    modfdc(i) = modfdc_a(num_valid-(i-1))
  enddo

!  do i = 1,num_valid
!    print *,order(i),obsfdc(i),modfdc(i)
!  enddo

  obsfdc_sp = obsfdc



  modslopemin = modfdc(fmsmin) - (modfdc(fmsmin) - modfdc(fmsmin+1)) / ( probexc(fmsmin) - probexc(fmsmin+1) ) * ( probexc(fmsmin) - 0.2)
  obsslopemin = obsfdc_sp(fmsmin) - (obsfdc_sp(fmsmin) - obsfdc_sp(fmsmin+1)) / ( probexc(fmsmin) - probexc(fmsmin+1) ) * ( probexc(fmsmin) - 0.2)
  modslopemax = modfdc(fmsmax) - (modfdc(fmsmax) - modfdc(fmsmax+1)) / ( probexc(fmsmax) - probexc(fmsmax+1) ) * ( probexc(fmsmax) - 0.7)
  obsslopemax = obsfdc_sp(fmsmax) - (obsfdc_sp(fmsmax) - obsfdc_sp(fmsmax+1)) / ( probexc(fmsmax) - probexc(fmsmax+1) ) * ( probexc(fmsmax) - 0.7)


!need to deal with zero observed or modeled streamflow when taking the log
  if(modslopemin .lt. 1e-6) then
    modslopemin = 1e-6
  endif
  if(obsslopemin .lt. 1e-6) then
    obsslopemin = 1e-6
  endif
  if(obsslopemax .lt. 1e-6) then
    modslopemax = 1e-6
  endif
  if(obsslopemax .lt. 1e-6) then
    modslopemax = 1e-6
  endif

  pbiasfms = ( (log10(modslopemin) - log10(modslopemax)) - (log10(obsslopemin) - log10(obsslopemax)) )/ (log10(obsslopemin) - log10(obsslopemax))*100

  return
end subroutine fdc_fms



subroutine fdc_fhvbias(model,obsfdc,order,length,valid,pbiasfhv)
  !modified from pablo mendoza's subroutine following yilmax et al. (2008)
  !only computes bias for high flow, top 2% of fdc

  use nrtype
  use nr, only: sort_heap

  implicit none


!input variables
  real(DP), dimension(36500),     intent(in)  :: model	! Modeled streamflow
  real(DP), dimension(36500),     intent(in)  :: obsfdc	! sorted observed streamflow 
  integer(I4B),dimension(36500),  intent(in)  :: order	! order of original observed streamflow

  integer(I4B),	intent(in)			:: length
  logical, dimension(36500),intent(in)		:: valid

!output variables
  real(DP), intent(out) :: pbiasfhv

!local variables
  integer(I4B) :: itime,num_valid,cnt,i
  integer(I4B) :: fhvindex   ! Maximum index that contains the high flows portion in FDC
  integer(I4B) :: flvindex   ! Minimum index that contains the low flows portion in FDC
  integer(I4B) :: fmsmin     ! Minimum index that contains the mid-segment portion in FDC
  integer(I4B) :: fmsmax     ! Maximum index that contains the mid-segment portion in FDC
  integer(I4B) :: end_pt
  real(DP), dimension(36500)	:: probexc	! exceedance probabilities
  real(SP), dimension(36500)	:: modfdc_a	! modeled FDC
  real(SP), dimension(36500)	:: modfdc	! modeled FDC
  real(SP), dimension(36500)	:: obsfdc_sp	! obseved fdc single precision

  end_pt = length

  num_valid = count(valid(1:end_pt))

 ! Compute the probabilities of exceedance and save relevant index thresholds
  do itime = 1,num_valid
    probexc(itime) = real(itime)/real(num_valid)
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
!  call sortpick(model,modfdc,modfdcmed)
!  call sortpick(obs,obsfdc,obsfdcmed)

!print *,num_valid

  do i = 1,num_valid
    if(obsfdc(i) .gt. -0.00000001) then
      modfdc_a(i) = real(model(order(i)),kind(SP))
!print *,model(order(i)),modfdc(i)
    endif
  enddo

!   sort modelfdc into descending order
  call sort_heap(modfdc_a(1:num_valid))
  do i = 1,num_valid
    modfdc(i) = modfdc_a(num_valid-(i-1))
  enddo

  do i = 1,num_valid
    print *,probexc(i),i,obsfdc(i),modfdc(i)
  enddo

  obsfdc_sp = obsfdc

! Compute the percent bias in high flow volumes
  pbiasfhv =  ( sum(modfdc(1:fhvindex)) - sum(obsfdc(1:fhvindex)) ) / sum(obsfdc(1:fhvindex)) * 100.0

  return
end subroutine fdc_fhvbias



subroutine sortpick(vectorin,end_pt,vectorout,order)
  use nrtype

  implicit none
! Sorts an array vector into descending numerical order, by straight insertion. 
! vectorin is replaced on output by its sorted rearrangement (vectorout)
! Additionally, the median of the array is also computed
! from pablo mendoza, from numerical recipies?

! input variables
  real(DP), dimension(36500), intent(in)	:: vectorin 	! input vector
  integer(I4B)					:: end_pt       !end of valid data

! output variables
  real(DP), dimension(36500), intent(out)	:: vectorout 	! output vector
  integer(I4B), dimension(36500), intent(out)	:: order 	! order of sorted vector (from initial input array)

! local variables
  integer(I4B) 	:: i,j,b
  real(DP) 	:: a

    
  vectorout = vectorin

! setup order of initial vector
  do i=1,end_pt
    order(i) = i
  enddo

  do j=2,end_pt ! Pick out each element in turn.
    a=vectorout(j)
    b=order(j)
    do i=j-1,1,-1 ! Look for the place to insert it.
      if (vectorout(i) >= a) exit
      vectorout(i+1)=vectorout(i)
      order(i+1)=order(i)
    enddo
    vectorout(i+1)=a !Insert it.
    order(i+1)=b
  enddo

!  do i = 1,end_pt
!    print *,order(i),vectorout(i)
!  enddo


  return
end subroutine sortpick

