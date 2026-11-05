C     Reproduces the CR-176883 16-series sample case (sample_naca16.txt).
      DIMENSION AM(2),AN(3)
      DATA AM/.45,.7/,AN/0.,2.,4./
      ALTUDE=0.
      TC=.04
      CLD=0.228150
      CALL ATMCON(ALTUDE,TEMP,PRES,DEN,VIS,WA)
      ALPHAO=-2.
      DO 1 K=1,2
      DO 1 J=1,3
      A=AN(J)-ALPHAO
      CALL CLCD16(AM(K),A,CL,CD,1,1,1.0,TC,RN,CLD,.5,WA,VIS,DEN,
     &CF,CDL,CDC,CDZ)
      WRITE(6,'(F6.2,F6.1,2F10.6,E14.6)') AM(K),AN(J),CL,CD,RN
    1 CONTINUE
      END
