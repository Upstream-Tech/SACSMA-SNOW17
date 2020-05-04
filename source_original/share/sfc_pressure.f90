subroutine sfc_pressure(elev, sfc_pres)
  use nrtype
  use constants, only: sfc_pres_a,sfc_pres_b,sfc_pres_c,&
                       sfc_pres_d,sfc_pres_e

  implicit none

  real(DP), intent(in)	 :: elev
  real(DP), intent(out)  :: sfc_pres
  
  sfc_pres = sfc_pres_a * (sfc_pres_b - (sfc_pres_c * (elev/100.0_dp)) &
             + (sfc_pres_d*((elev/100.0_dp)**sfc_pres_e)))   !sfc pres in hPa

  return
end subroutine sfc_pressure