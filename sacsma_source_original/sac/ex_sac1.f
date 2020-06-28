      SUBROUTINE EXSAC(NSOLD,DTM,PCP,TMP,ETP,
C     SAC PARAMETERS
     &                 UZTWM,UZFWM,UZK,PCTIM,ADIMP,RIVA,ZPERC,
     &                 REXP,LZTWM,LZFSM,LZFPM,LZSK,LZPK,PFREE,
     &                 SIDE,RSERV,
C     SAC State variables  ',
     &                 UZTWC,UZFWC,LZTWC,LZFSC,LZFPC,ADIMC,
C     SAC OUTPUTS
     &                 QS,QG,Q,ETA)

C      IMPLICIT NONE

C     RCS Id string, for version control
      CHARACTER*60 RCSID
C      DATA RCSID/"$Id: ex_sac1.f,v 1.1 2006/09/01 21:59:43 vicadmin Exp $"/


C     ...THIS SUBROUTINE IS THE EXECUTION ROUTINE FOR SMFLX MODEL...
C

      INTEGER NSOLD
      REAL    DTM
      REAL    PCP
      REAL    TMP
      REAL    ETP
      REAL    QS
      REAL    QG
      REAL    Q
      REAL    ETA


      REAL    LZTWM,LZFSM,LZFPM,LZSK,LZPK,LZTWC,LZFSC,LZFPC
      REAL    TOTAL_S1, TOTAL_S2
      REAL    DT
      REAL    DS

     
      COMMON/FSMCO1/FGCO(6),RSUM(7),PPE,PSC,PTA,PWE
      COMMON/FSUMS1/SROT,SIMPVT,SRODT,SROST,SINTFT,SGWFP,SGWFS,SRECHT,
     &              SETT,SE1,SE3,SE4,SE5

C    TURN OFF FROZEN GROUND PROCESS

      IFRZE = 0

C     COMPUTE TOTAL INITIAL STORAGE

      TOTAL_S1 = UZTWC + UZFWC + LZTWC + LZFSC + LZFPC + ADIMC

C     COMPUTE SURFACE MOISTURE FLUXES
      DT = DTM/86400.0
      EP1 = ETP
      P1 = PCP

      ! write(*,*) DT,P1,EP1,TCI,ROIMP,SDRO,SSUR,SIF,BFS,BFP,ETA
      ! write(*,*) 'hi', IFRZE,TA,LWE,WE,ISC,AESC
      ! write(*,*) UZTWC,UZFWC,LZTWC,LZFSC,LZFPC,ADIMC
      ! write(*,*) '-------------------------'

      CALL SAC1(DT,P1,EP1,TCI,ROIMP,SDRO,SSUR,SIF,BFS,BFP,ETA,
C     SAC FROZEN GROUND VARIABLES
     &            IFRZE,TA,LWE,WE,ISC,AESC,
C     SAC PARAMETERS
     &            UZTWM,UZFWM,UZK,PCTIM,ADIMP,RIVA,ZPERC,
     &            REXP,LZTWM,LZFSM,LZFPM,LZSK,LZPK,PFREE,
     &            SIDE,RSERV,
C     SAC State variables  ',
     &            UZTWC,UZFWC,LZTWC,LZFSC,LZFPC,ADIMC)

C     COMPUTE FINAL TOTAL STORAGE AND WATER BALANCE
C     sdro:  direct runoff
C     ROIMP: impervious area runoff
C     SSUR:  surface runoff
C     SIF:   interflow
C     BFS:   non-channel baseflow
C     BFP:   some kind of baseflow...
C     TCI:   Total channel inflow

      QS = ROIMP + SDRO + SSUR + SIF
      QG = BFS + BFP
      Q  = TCI
C
      TOTAL_S2 = UZTWC + UZFWC + LZTWC + LZFSC + LZFPC + ADIMC
      DS = (TOTAL_S2 - TOTAL_S1)

      BAL = P1-ETA-QS-QG-DS

C      PRINT*,'exsac1 -',BAL,P1,ETA,Q,QS,QG,DS,TOTAL_S1,TOTAL_S2

      RETURN
      END
