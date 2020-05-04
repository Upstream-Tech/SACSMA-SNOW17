FUNCTION FUNCTN(NOPT,A)
  use nrtype
  use gauge_calib, only: snow17_sacsma_eval

  IMPLICIT NONE

!input variables
  integer(I4B),		intent(in) 	:: nopt
  real(sp), dimension(30),	intent(inout)	:: a

!local variables 
  real(dp) :: obj_val
  real(sp) :: functn


  call snow17_sacsma_eval(a,obj_val)

! save objective function value
  FUNCTN = real(obj_val,kind(sp))
! ---------------------------------------------------------------------------------------
END FUNCTION FUNCTN