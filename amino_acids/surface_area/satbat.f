      SUBROUTINE SATBAT (CSTOT,NFILE)                                   GDH0897
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               GDH0897
      INCLUDE 'KEYS.i'                                                  GDH0897
      INCLUDE 'SIZES.i'                                                 GDH0897
C***********************************************************************GDH0897
C     THIS SUBROUTINE PRINTS OUT THE ATOM BY ATOM SUMMARY OF THE SOLVATIGDH0897
C     CALCULATION.                                                      GDH0897
C ON INPUT                                                              GDH0897
C     ELEMNT = ELEMENT SYMBOL                                           GDH0897
C     SRFACT = CDS ENERGY CONTRIBUTION      (kcal)                      GDH0897
C                                                                       GDH0897
C***********************************************************************GDH0897
      COMMON /ELEMTS/ ELEMNT(107)                                       GDH0897
      COMMON /MLKSTI/ NUMAT,NAT(NUMATM)                                 GDH0897
      COMMON /SURF  / SURFCT,SRFACT(NUMATM),ATAR(NUMATM),               GDH0897
     1                HEXLGS,ATLGAR,CSAREA(100)                         GDH0897
      COMMON /TRADCM/ IRAD                                              GDH0897
      CHARACTER*2 ELEMNT                                                GDH0897
      SAVE                                                              GDH0897
C                                                                       GDH0897
         TOTPFE = 0.D0                                                  GDH0897
         TOTAR  = 0.D0                                                  GDH0897
         TOTCDS = 0.D0                                                  GDH0897
         TOTSFE = 0.D0                                                  GDH0897
         DO 150 L = 1, NUMAT                                            GDH0897
             TOTAR  = TOTAR  + ATAR(L)                                  GDH0897
             TOTCDS = TOTCDS + SRFACT(L)                                GDH0897
             TOTSFE = TOTSFE + SRFACT(L)                                GDH0897
150      CONTINUE                                                       GDH0897
         IF (IRAD.EQ.21)THEN                                            GDH0897
            CSTOT=ATLGAR*HEXLGS                                         GDH0897
            TOTCDS=TOTCDS + CSTOT                                       GDH0897
            TOTSFE=TOTSFE + CSTOT                                       GDH0897
         ENDIF                                                          GDH0897
            WRITE(NFILE,1000)                                           GDH0897
      DO 201 L = 1, NUMAT                                               GDH0897
         WRITE(NFILE,1100)L, ELEMNT(NAT(L)), ATAR(L), SRFACT(L)         GDH0897
201   CONTINUE                                                          GDH0897
         IF (IRAD.EQ.21)                                                GDH0897
     1       WRITE(NFILE,1200) ATLGAR, CSTOT                            GDH0897
         WRITE(NFILE, 1400) TOTAR, TOTSFE                               GDH0897
1000  FORMAT(/,2X,'Atom', T12, 'Chemical', T25, 'Area', T35,            GDH0897
     1       'SM5.0R Solv.',                                            GDH0897
     2       /,1X,'number', T13,'symbol',                               GDH0897
     3        T23,'(Ang**2)',T35, 'free energy',                        GDH0897
     4        /, 1X, T38,                                               GDH0897
     5       '(kcal)',/)                                                GDH0897
1100  FORMAT(1X,I3,T15,A2,T23,F9.5,T35,F9.5)                            GDH0897
1200  FORMAT(/,1X,'CS Contribution',T23,F7.2,T35,F8.2)                  GDH0897
1400  FORMAT(/,1X,'Total:',T23,F7.2,T35,F8.2)                           GDH0897
      RETURN                                                            GDH0897
      END                                                               GDH0897
