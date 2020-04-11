      SUBROUTINE SNEW ( T,P,HCN,DHN )
     
C  CALCULATING SNOW DEPTH AND DENSITITY OF THE NEW SNOWFALL
C      T   - AIR TEMPERATURE, C
C      P   - NEW SNOWFALL, MM 
C      HCN - SNOW DEPTH, CM
C      DHN - SNOW DENSITY
C   NEW VALUES OF SNOW DEPTH & DENSITY WILL BE RETURNED   

C ----------------------------------------------------------------------
C ***  CALCULATING NEW SNOWFALL DENSITY DEPENDING ON TEMPERATURE      **
C ***  EQUATION FROM GOTTLIB L. 'A GENERAL RUNOFF MODEL FOR SNOWCOVERED
C ***  AND GLACIERIZED BASIN', 6TH NORDIC HYDROLOGICAL CONFERENCE, 
C ***  VEMADOLEN, SWEDEN, 1980, 172-177pp.
C-----------------------------------------------------------------------

C  CONVERT IN CALCULATION UNITS
      PX=0.1*P
      
      IF(T.LE.-15.) THEN
        DHN=0.05 
      ELSE                                                      
        DHN=0.05+0.0017*(T+15.)**1.5
      ENDIF

C **   CALCULATE CHANGE OF SNOW DEPTH       

      HCN=PX/DHN
      
      RETURN
      END      

