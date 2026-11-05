C     Reproduces the CR-176883 4-digit sample case (sample_naca4.txt).
      LOGICAL LF
      DIMENSION ALPHA(10)
      DATA ALPHA / -10.0,-8.0,-6.0,-4.0,-2.0,0.0,2.0,4.0,6.0,8.0/
      LF = .FALSE.
      DO 1 NN = 1,10
      RE = 0.
      ALT = 10000.
      CLD = 0.2
      CALL NACA4(0.,0.12,CLD,Z,ALPHA(NN),0.5,RE,ALT,1.0,LF,
     &CL,CD,CLMX,CLOPT,X,CDMIN,DELCD)
      WRITE(6,'(F6.1,2F10.6,E14.6,5F10.6)') ALPHA(NN),CL,CD,RE,Z,X,
     &CDMIN,DELCD
    1 CONTINUE
      END
