C MEMBER ZERO19
C  (from old member FCPACK19)
C
      SUBROUTINE ZERO19
C.......................................
C     THIS SUBROUTINE SETS ALL CARRYOVER VALUES TO NO SNOW CONDITIONS
C        FOR THE 'SNOW-17 ' OPERATION.
C.......................................
C     SUBROUTINE INITIALLY WRITTEN BY...
C        ERIC ANDERSON - HRL   MAY 1980

C    MODIFIED 4/20/00 BY V. KOREN TO ADD TWO MORE STATES: SNDPT & SNTMP
C.......................................
      REAL LIQW,NEGHS
C
C     COMMON BLOCK
      COMMON/SNCO19/WE,NEGHS,LIQW,TINDEX,ACCMAX,SB,SBAESC,SBWS,STORGE,

CVK     1   AEADJ,NEXLAG,EXLAG(7)
     1   AEADJ,NEXLAG,EXLAG(7),SNDPT,SNTMP
C
C.......................................
      WE=0.0
      NEGHS=0.0
      LIQW=0.0
      TINDEX=0.0
      ACCMAX=0.0
      SB=0.0
      SBAESC=0.0
      SBWS=0.0
      STORGE=0.0
      AEADJ=0.0

CVK  ADDED TWO MORE STATES
      SNDPT=0.0
      SNTMP=0.0
            
      DO 100 N=1,NEXLAG
  100 EXLAG(N)=0.0
C.......................................
      RETURN
      END

