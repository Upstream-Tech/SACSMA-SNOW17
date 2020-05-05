      SUBROUTINE AESC19(WE,LIQW,ACCMAX,SB,SBAESC,SBWS,SI,ADC,AEADJ,AESC)
C.......................................
C     THIS SUBROUTINE COMPUTES THE AREAL EXTENT OF SNOW COVER USING THE
C        AREAL DEPLETION CURVE FOR THE 'SNOW-17 ' OPERATION.
C.......................................
C     SUBROUTINE INITIALLY WRITTEN BY...
C        ERIC ANDERSON - HRL   MAY 1980
C.......................................
      REAL LIQW
      DIMENSION ADC(11)
C
C     COMMON BLOCKS
      COMMON/SNUP19/MFC,SFALLX,WINDC,SCTOL,WETOL,SNOF
C.......................................
      TWE=WE+LIQW
C      print *,'1',twe,accmax
      IF(TWE.GT.ACCMAX) ACCMAX=TWE
      IF (TWE.GE.AEADJ) AEADJ=0.0
      AI=ACCMAX
C      print *,'2',ai,twe,accmax
      IF(ACCMAX.GT.SI)AI=SI
      IF (AEADJ.GT.0.0) AI=AEADJ
C      print *,'3',ai,aeadj,twe,accmax,sb,sbws,tiny
      IF(TWE.GE.AI) GO TO 105
C      IF(TWE.LE.SB) GO TO 110
C  AJN Need tiny fudge factor for restart output
      IF(TWE.LE.SB+tiny) GO TO 110
      IF(TWE.GE.SBWS) GO TO 115
      AESC=SBAESC+((1.0-SBAESC)*((TWE-SB)/(SBWS-SB)))
      GO TO 120
  110 R=(TWE/AI)*10.0+1.0
      N=R
      FN=N
      R=R-FN
      AESC=ADC(N)+(ADC(N+1)-ADC(N))*R
C      print *,'rrrr',twi,ai,r,n,fn
      IF(AESC.GT.1.0) AESC=1.0
      SB=TWE+SNOF
      SBWS=TWE
      SBAESC=AESC
     GO TO 120
  105 SB=TWE
      SBWS=TWE
  115 AESC=1.0
  120 IF(AESC.LT.0.05) AESC=0.05
      IF(AESC.GT.1.0) AESC=1.0
C.......................................
      RETURN
      END

