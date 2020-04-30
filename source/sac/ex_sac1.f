      SUBROUTINE EXSAC(DTM, PCP, TMP, ETP,                  & ! FORCINGS (IN)           
            UZTWM, UZFWM, UZK, PCTIM, ADIMP, RIVA, ZPERC,   & ! PARAMETERS (IN)
            REXP, LZTWM, LZFSM, LZFPM, LZSK, LZPK, PFREE,   & ! PARAMETERS (IN)
            SIDE, RSERV,                                    & ! PARAMETERS (IN)
            UZTWC, UZFWC, LZTWC, LZFSC, LZFPC, ADIMC,       & ! STATES (INOUT)
            QS, QG, Q, ETA)                                 & ! OUTPUTS (OUT)
            
      IMPLICIT NONE  

      ! FORCINGS
      REAL, INTENT(IN)  ::    DTM ! timestep in seconds
      REAL, INTENT(IN)  ::    PCP ! precip in [units]
      REAL, INTENT(IN)  ::    TMP ! air temp in C
      REAL, INTENT(IN)  ::    ETP ! potential evapotranspiration in [units]
  
      ! PARAMETERS
      REAL, INTENT(IN)  ::  UZTWM, UZFWM, UZK, PCTIM, ADIMP, RIVA, ZPERC
      REAL, INTENT(IN)  ::  REXP, LZTWM, LZFSM, LZFPM, LZSK, LZPK, PFREE
      REAL, INTENT(IN)  ::  SIDE, RSERV
  
      ! STATES
      REAL, INTENT(INOUT)  ::  UZTWC, UZFWC, LZTWC, LZFSC, LZFPC, ADIMC
  
      ! OUTPUTS
      REAL, INTENT(OUT)  ::  QS
      REAL, INTENT(OUT)  ::  QG
      REAL, INTENT(OUT)  ::  Q
      REAL, INTENT(OUT)  ::  ETA
  
      ! INTERNALS
      REAL  ::  LZTWM, LZFSM, LZFPM, LZSK, LZPK, LZTWC, LZFSC, LZFPC
      REAL  ::  TOTAL_S1, TOTAL_S2
      REAL  ::  DT  ! timestep in days
      REAL  ::  DS
      INTEGER : ISC = 1
  
      ! SOME SHARED VARIABLES
      REAL, DIMENSION(7) :: RSUM
      REAL, DIMENSION(6) :: FGCO
      COMMON/FSMCO1/FGCO(6),RSUM(7)
   
     ! TURN OFF FROZEN GROUND PROCESS
      IFRZE = 0
  
      ! COMPUTE TOTAL INITIAL STORAGE (for debugging)
      TOTAL_S1 = UZTWC + UZFWC + LZTWC + LZFSC + LZFPC + ADIMC
  
      ! UNIT CONVERSIONS
      DT = DTM/86400.0     ! timestep from seconds to days
      EP1 = ETP
      P1 = PCP
      
      ! RUN FLAND1
      CALL SAC1(DT, P1, EP1,                                      & ! FFORCINGS
       TCI, ROIMP, SDRO, SSUR, SIF, BFS, BFP, ETA,                & ! OUTPUTS
       IFRZE, TA, LWE, WE, ISC, AESC,                             & ! FROZEN GROUND VARIABLES
       UZTWM, UZFWM, UZK, PCTIM, ADIMP, RIVA, ZPERC,              & ! PARAMETERS
       REXP, LZTWM, LZFSM, LZFPM, LZSK, LZPK, PFREE,              & ! PARAMETERS
       SIDE, RSERV,                                               & ! PARAMETERS
       UZTWC, UZFWC, LZTWC, LZFSC, LZFPC, ADIMC)                  & ! STATES
  
      ! --- COMPUTE FINAL TOTAL STORAGE AND WATER BALANCE ---
      ! SDRO:  direct runoff
      ! ROIMP: impervious area runoff
      ! SSUR:  surface runoff
      ! SIF:   interflow
      ! BFS:   non-channel baseflow
      ! BFP:   some kind of baseflow...
      ! TCI:   Total channel inflow
      QS = ROIMP + SDRO + SSUR + SIF
      QG = BFS + BFP
      Q  = TCI
  
      RETURN
      END SUBROUTINE EXSAC
  