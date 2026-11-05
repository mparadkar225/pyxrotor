C     CR-176883 Clark Y sample cases (report Appendix II, PDF pp. 74-86):
C     the main program's arithmetic for case 1 (T/C 0.117, chord 0.75 ft,
C     10000 ft, Mach 0, 0.1, 0.2, 0.4, alpha -1..5 deg) and case 2 (T/C
C     0.30, Mach 0.6, 0.85, alpha -2..10 deg). Prints case, Mach, RN,
C     alpha, CL, CD as in sample_clarky.txt.
      DIMENSION HM1(4),A1(7),HM2(2),A2(7)
      REAL NCLD
      DATA HM1/0.0,0.1,0.2,0.4/,A1/-1.,0.,1.,2.,3.,4.,5./
      DATA HM2/0.6,0.85/,A2/-2.,0.,2.,4.,6.,8.,10./
      CHORD=0.75
      DO 30 IC=1,2
      ALTUDE=10000.
      CALL ATMCON(ALTUDE,TEMP,PRES,DEN,VIS,WA)
      NCLD=0.
      T=0.117
      IF(IC.EQ.2) T=0.30
      CALL CLTMAX(T,NCLD,ALPHLO,DNCLD)
      NM=4
      IF(IC.EQ.2) NM=2
      DO 20 K=1,NM
      HM=HM1(K)
      IF(IC.EQ.2) HM=HM2(K)
      RN=DEN*HM*WA*CHORD/VIS
      DO 10 J=1,7
      A=A1(J)
      IF(IC.EQ.2) A=A2(J)
      ALPHAO=A-ALPHLO
      CALL CLCDCY(ALPHAO,T,HM,CL,CD,1,1)
      WRITE(6,'(I2,F5.2,E14.6,F6.1,2F11.6)') IC,HM,RN,A,CL,CD
   10 CONTINUE
   20 CONTINUE
   30 CONTINUE
      END
