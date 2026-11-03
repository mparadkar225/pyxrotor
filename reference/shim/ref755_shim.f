C=======================================================================
C     Non-interactive entry points for the XROTOR 7.55 reference library.
C
C     Each routine transcribes an interactive XROTOR routine with its
C     terminal prompts replaced by arguments. Everything else is copied
C     statement for statement, so the reference computes exactly what
C     the interactive program does. Keep these in sync with the source
C     cited in each header; do not "improve" them.
C
C     Everything not covered here is called directly from Python
C     (SETX, APER, DESGEN, PUTAERO, ...), with inputs set in COMMON.
C=======================================================================

      SUBROUTINE REFSTART
C-----------------------------------------------------------------------
C     Program start-up, from PROGRAM XROTOR (xrotor.f), up to the menu.
C-----------------------------------------------------------------------
      INCLUDE 'XROTOR.INC'
C
      VERSION = 7.55
      LUREAD = 5    ! terminal read
      LUWRIT = 6    ! terminal write
      LUTEMP = 3    ! general-use disk I/O unit  (usually available)
      LUSAVE = 4    ! save file                  (usually open)
C
      CALL INIT
      RETURN
      END ! REFSTART



      SUBROUTINE REFARBI(NB, VELIN, SGNRAD, ROOTIN, NSTIN, X, C, B)
C-----------------------------------------------------------------------
C     Arbitrary rotor geometry input, from SUBROUTINE ARBI (xrotor.f).
C
C     NB      number of blades
C     VELIN   flight speed (m/s)
C     SGNRAD  tip radius (m); negative means X and C are dimensional
C     ROOTIN  hub radius (m)
C     NSTIN   number of input stations
C     X, C, B r/R, chord/R, blade angle (deg) at each station
C-----------------------------------------------------------------------
      INCLUDE 'XROTOR.INC'
      DIMENSION X(NSTIN), C(NSTIN), B(NSTIN)
C
      GREEK = .FALSE.
C
      CONV = .FALSE.
      NBLDS = NB
      VEL = VELIN
      RAD = ABS(SGNRAD)
      ROOT = ROOTIN
      IF(RAD.LE.ROOT) STOP 'REFARBI: RAD <= ROOT'
      XI0 = ROOT/RAD
      XW0 = ROOT/RAD
      CALL SETX
      CALL XWINIT
C
      NST = NSTIN
      IF(NST.LT.3) STOP 'REFARBI: NST < 3'
C
      DO 20 N=1, NST
        W1(N) = X(N)
        W2(N) = C(N)
        W3(N) = B(N)
        IF(SGNRAD.LT.0.0) THEN
         W1(N) = W1(N)/RAD
         W2(N) = W2(N)/RAD
        ENDIF
        IF(W1(N).LT.XI0 .OR. W1(N).GT.XITIP)
     &    STOP 'REFARBI: radii must be between XI0 and 1.0'
        IF(N.GT.1 .AND. W1(N).LE.W1(N-1))
     &    STOP 'REFARBI: radii must monotonically increase'
   20 CONTINUE
C
      DO N=1, NST
        W3(N) = W3(N)*PI/180.
        W1(N) = TINVRT(W1(N))
      END DO
C
      CALL SPLINE(W2,W4,W1,NST)
      CALL SPLINE(W3,W5,W1,NST)
C
      DO I=1, II
        CH(I)   = SEVAL(T(I),W2,W4,W1,NST)
        BETA(I) = SEVAL(T(I),W3,W5,W1,NST)
        BETA0(I) = BETA(I)
        CH(I) = ABS(CH(I))                   ! negative chords are a no-no
        UBODY(I) = 0.
      END DO
C
      IF(NAME.EQ.' ') NAME = 'Arbitrary blade'
C
C---- estimate reasonable advance ratio to start iterative routines
      IS = II/2 + 1
C---HHY had to set A0 to 0.0 as A0 is now section property
      A0  = 0.0
      ANG = BETA(IS) - A0
      ADV = XI(IS)*SIN(ANG)/COS(ANG)
      ADV = MAX(0.1,ADV)
      ADW = ADV
C
      CALL SETIAERO
      RETURN
      END ! REFARBI



C=======================================================================
C     OPER menu command handlers (xoper.f SUBROUTINE OPER), with the
C     prompts replaced by arguments. Plotting calls are omitted; they
C     do not change XROTOR state.
C=======================================================================

      SUBROUTINE REFATMO(ALTIN)
C---- OPER/top-level ATMO r (xoper.f label 35, xrotor.f)
      INCLUDE 'XROTOR.INC'
      ALT = ALTIN
      CALL ATMO(ALT,VSO,RHO,RMU)
      RETURN
      END


      SUBROUTINE REFVELO(VELIN)
C---- OPER VELO r (xoper.f label 38)
      INCLUDE 'XROTOR.INC'
      VELOLD = VEL
      VEL = VELIN
C--- Change CT,CQ,CP to give same thrust,torque,power
      THR  = TTOT *(RHO*VELOLD**2*RAD**2)
      TTOT = THR / (RHO*VEL**2*RAD**2)
      TRQ  = QTOT *(RHO*VELOLD**2*RAD**3)
      QTOT = TRQ / (RHO*VEL**2*RAD**3)
      PWR  = PTOT *(RHO*VELOLD**3*RAD**2)
      PTOT = PWR / (RHO*VEL**3*RAD**2)
      CONV = .FALSE.
      RETURN
      END


      SUBROUTINE REFANGL(DELB)
C---- OPER ANGL r (xoper.f label 40)
      INCLUDE 'XROTOR.INC'
      DO I=1, II
        BETA(I)  = BETA(I)  + DELB*PI/180.
        BETA0(I) = BETA0(I) + DELB*PI/180.
      ENDDO
      CONV = .FALSE.
      RETURN
      END


      SUBROUTINE REFADVA(ADVIN)
C---- OPER ADVA r (xoper.f label 42)
      INCLUDE 'XROTOR.INC'
      ADV = ADVIN
      CONV = .FALSE.
      CALL APER(4,2,LOPRINI)
      RETURN
      END


      SUBROUTINE REFRPM(RPM)
C---- OPER RPM r (xoper.f label 45)
      INCLUDE 'XROTOR.INC'
      ADV = VEL / (RAD*RPM*PI/30.)
      CONV = .FALSE.
      CALL APER(4,2,LOPRINI)
      RETURN
      END


      SUBROUTINE REFTPQ(ISPEC,SPEC,IFIX,RPMIN)
C---- OPER THRU / TORQ / POWE (xoper.f labels 50, 60, 70)
C     ISPEC 1,2,3 = thrust, torque, power; IFIX 1 = fix pitch (P),
C     2 = fix rpm (R) at RPMIN (RPMIN <= 0: keep the current rpm)
      INCLUDE 'XROTOR.INC'
      IF(ISPEC.EQ.1) TSPEC = SPEC
      IF(ISPEC.EQ.2) QSPEC = SPEC
      IF(ISPEC.EQ.3) PSPEC = SPEC
      RPM = VEL / (RAD*ADV*PI/30.0)
      CONV = .FALSE.
      IF(IFIX.EQ.1) CALL APER(ISPEC,2,LOPRINI)
      IF(IFIX.EQ.2) THEN
       IF(RPMIN.GT.0.0) RPM = RPMIN
       ADV = VEL / (RAD*RPM*PI/30.0)
       CALL APER(ISPEC,1,LOPRINI)
      ENDIF
C---- Check for valid blade angle change
      IF(IFIX.NE.1) THEN
       IF(.NOT.CONV) THEN
C----- convergence failed: restore clobbered blade angles
        DO I=1, II
         BETA(I)  = BETA(I)  - DBETA
         BETA0(I) = BETA0(I) - DBETA
        ENDDO
       ENDIF
      ENDIF
      RETURN
      END



      SUBROUTINE REFDESGEN(OKIN,LFIXIA)
C-----------------------------------------------------------------------
C     DESGEN (xdesi.f) with the Betz-limit prompt answered by OKIN.
C     Returns CONV unchanged if DESGEN returns early.
C     LFIXIA (not in 7.55): CALL SETIAERO on the new grid after XWINIT,
C     the fix for the stale aero-section indices.
C-----------------------------------------------------------------------
      INCLUDE 'XROTOR.INC'
      LOGICAL OK, OKIN, LFIXIA
C
      DO I=1, II
        N = 1
        CALL GETAERO(N,XISECT,A0,CLMAX,CLMIN,
     &              DCLDA,DCLDA_STALL,DCL_STALL,
     &              CDMIN,CLDMIN,DCDCL2,CMCON,MCRIT,REREF,REXP)
        IF(CLDES(I).GT.CLMAX .OR. CLDES(I).LT.CLMIN) RETURN
      ENDDO
C
C---- Set up operating parameters
      CONV = .FALSE.
      BLDS = FLOAT(NBLDS)
      RAD = RADDES
      VEL = VELDES
      IF(ADVDES.EQ.0.0) THEN
       ADV = VELDES/(RADDES*RPMDES*PI/30.0)
      ELSE
       ADV = ADVDES
      ENDIF
C
C---- Check for reasonable thrust or power inputs
      TMIN = -0.5*PI*RHO*VEL**2*RAD**2
      PMIN = -0.5*PI*RHO*VEL**3*RAD**2 * 16.0/27.0
      IF((DEST .AND. TDDES.LT.TMIN).OR.(DESP .AND. PDDES.LT.PMIN)) THEN
       OK = OKIN
       IF(.NOT.OK) RETURN
      ENDIF
C
C---- Non-dimensionalize design thrust or power
      TDES = TDDES/(RHO*VEL**2*RAD**2)
      PDES = PDDES/(RHO*VEL**3*RAD**2)
C
C---- Set up discrete coordinates
      XI0 = R0DES/RAD
      XW0 = RWDES/RAD
      CALL SETX
      CALL XWINIT
      IF(LFIXIA) CALL SETIAERO
C
C---- Fill CL array with design-CL array
      DO I=1, II
        CL(I) = CLDES(I)
        STALL(I) = .FALSE.
      ENDDO
C
C---- Initialize wake advance ratios
      EFFINV = 1.0
      IF(FREE.AND.DEST) EFFINV = 0.50*SQRT(ABS(1.0+TDES*2./PI)) + 0.50
      IF(FREE.AND.DESP) EFFINV = 0.25*SQRT(ABS(1.0+PDES*4./PI)) + 0.75
      UDUCT     = 0.0
      IF(DUCT) THEN
        UDUCT = URDUCT-1.0
      ENDIF
      ADW = EFFINV*ADV
C
C---- Initialize rotor for iteration
      IF(LDESINI .OR. .NOT.LROTOR) CALL DESINI
C
C---- Calculate MIL rotor (varying chord, varying twist, fixed CL case)
      CALL DESMIL(2)
      RETURN
      END ! REFDESGEN



      SUBROUTINE REFDINPU(NBIN,RTIN,RHIN,RWIN,VIN,AIN,RPIN,
     &                    TIN,PIN,CLIN,OKIN,LFIXIA)
C-----------------------------------------------------------------------
C     DESI INPU (xdesi.f) with prompts as arguments, then the EDIT
C     <Return> that follows it (LDESINI, DESGEN). TIN = PIN = 0 is
C     re-prompted in XROTOR; here it returns without designing.
C-----------------------------------------------------------------------
      INCLUDE 'XROTOR.INC'
      LOGICAL OKIN, LFIXIA
C
      NBLDS = NBIN
      RADDES = RTIN
      R0DES = RHIN
      RWDES = RWIN
      VELDES = VIN
      ADVDES = 0.
      RPMDES = 0.
      ADVDES = AIN
      IF(ADVDES.EQ.0.0) RPMDES = RPIN
      TDDES = 0.
      PDDES = 0.
      TDDES = TIN
      IF(TDDES.EQ.0.0) PDDES = PIN
      IF(TDDES.EQ.0.0 .AND. PDDES.EQ.0.0) RETURN
      DEST = TDDES .NE. 0.0
      DESP = .NOT.DEST
      CLDES0 = CLIN
      IF(CLDES0 .EQ. 0.0) THEN
        IF(DEST) CLDES1 = SIGN(0.8,TDDES)
        IF(DESP) CLDES1 = SIGN(0.8,PDDES)
        CALL SETCLD(CLDES1,CLDES1)
      ELSE
        CALL SETCLD(CLDES0,CLDES0)
      ENDIF
      IF(NAME.EQ.' ') NAME = 'Designed blade'
      LROTOR = .FALSE.
C
      LDESINI = .TRUE.
      CALL REFDESGEN(OKIN,LFIXIA)
      RETURN
      END ! REFDINPU



      SUBROUTINE REFDESN(IINEW)
C-----------------------------------------------------------------------
C     DESI N (xdesi.f): change the number of radial stations.
C-----------------------------------------------------------------------
      INCLUDE 'XROTOR.INC'
C
      IF(LROTOR) THEN
       IISAV = II
       DO I = 1, IISAV
         W1(I) = XI(I)
         W2(I) = CH(I)
         W4(I) = BETA(I)
         W6(I) = UBODY(I)
         W8(I) = CLDES(I)
       ENDDO
       CALL SPLINE(W2,W3,W1,II)
       CALL SPLINE(W4,W5,W1,II)
       CALL SPLINE(W6,W7,W1,II)
       CALL SPLINE(W8,W9,W1,II)
      ENDIF
C
      II = IINEW
      IINF = II + II/2
      CALL SETX
      IF(LROTOR) THEN
       DO I = 1, II
         CH(I)    = SEVAL(XI(I),W2,W3,W1,IISAV)
         BETA(I)  = SEVAL(XI(I),W4,W5,W1,IISAV)
         UBODY(I) = SEVAL(XI(I),W6,W7,W1,IISAV)
         CLDES(I) = SEVAL(XI(I),W8,W9,W1,IISAV)
         BETA0(I) = BETA(I)
       ENDDO
      ENDIF
      RETURN
      END ! REFDESN



      SUBROUTINE REFFMT(IFMT,X,IW,ID)
C-----------------------------------------------------------------------
C     Formatted write of X into COMMON FNAME, for testing pyxrotor's
C     Fortran output: IFMT 1 Gw.d, 2 Ew.d, 3 Fw.d, 4 Iw (of INT(X)).
C-----------------------------------------------------------------------
      INCLUDE 'XROTOR.INC'
      CHARACTER*20 F
      IF(IFMT.EQ.1) WRITE(F,'(''(G'',I2,''.'',I2,'')'')') IW, ID
      IF(IFMT.EQ.2) WRITE(F,'(''(E'',I2,''.'',I2,'')'')') IW, ID
      IF(IFMT.EQ.3) WRITE(F,'(''(F'',I2,''.'',I2,'')'')') IW, ID
      IF(IFMT.EQ.4) WRITE(F,'(''(I'',I2,'')'')') IW
      FNAME = ' '
      IF(IFMT.LE.3) WRITE(FNAME,F) X
      IF(IFMT.EQ.4) WRITE(FNAME,F) INT(X)
      RETURN
      END ! REFFMT



      SUBROUTINE REFRDL(N,A,IERR)
C-----------------------------------------------------------------------
C     List-directed READ of N reals from COMMON FNAME (IERR=1: error).
C-----------------------------------------------------------------------
      INCLUDE 'XROTOR.INC'
      DIMENSION A(*)
      IERR = 1
      READ(FNAME,*,ERR=10,END=10) (A(I),I=1,N)
      IERR = 0
 10   RETURN
      END ! REFRDL



      SUBROUTINE REFGETFLT(A,N,IERR)
C-----------------------------------------------------------------------
C     GETFLT on COMMON FNAME (N in: max count, 0 = all; out: count).
C-----------------------------------------------------------------------
      INCLUDE 'XROTOR.INC'
      DIMENSION A(*)
      LOGICAL ERROR
      CALL GETFLT(FNAME,A,N,ERROR)
      IERR = 0
      IF(ERROR) IERR = 1
      RETURN
      END ! REFGETFLT



      SUBROUTINE REFSETCAS(ITYPE,R1,R2,R3,YESIN,IANS,CSIN,RPMIN,
     &                     LPIN,PIN)
C-----------------------------------------------------------------------
C     SETCAS (xoper.f) with first, last, increment given (NINPUT = 3)
C     and its prompts answered by arguments:
C       YESIN  'Use engine rpm/power line ?' (asked only if LPWRVAR)
C       IANS   RSEQ 'Fix power P or thrust T or blade pitch A ?':
C              0 A (or anything else), 1 T, 2 Q, 3 P
C       CSIN   VSEQ: .TRUE. constant-speed (CS), .FALSE. fixed (FP)
C       RPMIN  VSEQ CS rpm; 0 keeps the prompt's default
C       LPIN,PIN  VSEQ CS power; LPIN .FALSE. keeps the default
C-----------------------------------------------------------------------
      INCLUDE 'XROTOR.INC'
      LOGICAL YES, YESIN, CSIN, LPIN
      CHARACTER*1 ANS, ANS4*4
C
      KCASE = 0
C
      IF(ITYPE.EQ.1) THEN
        KCASE = 1
        ADV1 = R1
        ADV2 = R2
        DADV = R3
        IF(ADV1.EQ.ADV2) RETURN
        DADV = SIGN(DADV,ADV2-ADV1)
        NP = 1
        IF(DADV .NE. 0.0) NP = INT((ADV2-ADV1)/DADV + 0.5) + 1
        IF(NP.LE.0) RETURN
        YES = .FALSE.
        XANS = 0.
        IF(LPWRVAR) YES = YESIN
        IF(YES) XANS = 100.0
        IF(NCASE+NP .GT. ICASX) THEN
         NP = ICASX - NCASE
        ENDIF
        DO IP=1, NP
          NCASE = NCASE + 1
          CASPAR(0,NCASE) = XANS + FLOAT(KCASE)
          CASPAR(1,NCASE) = ADV1 + DADV*FLOAT(IP-1)
          CASPAR(2,NCASE) = VEL
          CASPAR(3,NCASE) = BETA(II)
          CASPAR(4,NCASE) = ALT
          CASPAR(5,NCASE) = RHO
          CASPAR(6,NCASE) = RMU
          CASPAR(7,NCASE) = VSO
          CASPAR(8,NCASE) = 999.
          CASPAR(9,NCASE) = 999.
          CASPAR(10,NCASE)= 999.
          CASPAR(11,NCASE)= 999.
        ENDDO
C
      ELSEIF(ITYPE.EQ.2) THEN
        KCASE = 2
        RPM1 = R1
        RPM2 = R2
        DRPM = R3
        IF(RPM1.EQ.RPM2) RETURN
        DRPM = SIGN(DRPM,RPM2-RPM1)
        NP = 1
        IF(DRPM .NE. 0.0) NP = INT((RPM2-RPM1)/DRPM + 0.5) + 1
        IF(NP.LE.0) RETURN
        YES = .FALSE.
        XANS = 0.
        IF(LPWRVAR) YES = YESIN
        IF(YES) XANS = 100.0
        ANS = ' '
        IF(IANS.EQ.1) ANS = 'T'
        IF(IANS.EQ.2) ANS = 'Q'
        IF(IANS.EQ.3) ANS = 'P'
        IF(ANS.EQ.'T') XANS = XANS + 1000.0
        IF(ANS.EQ.'Q') XANS = XANS + 2000.0
        IF(ANS.EQ.'P') XANS = XANS + 3000.0
        IF(NCASE+NP .GT. ICASX) THEN
         NP = ICASX - NCASE
        ENDIF
        DO IP=1, NP
          NCASE = NCASE + 1
          RPM = RPM1 + DRPM*FLOAT(IP-1)
          CASPAR(0,NCASE) = XANS + FLOAT(KCASE)
          CASPAR(1,NCASE) = VEL/(RPM*RAD) * 30.0/PI
          CASPAR(2,NCASE) = VEL
          CASPAR(3,NCASE) = BETA(II)
          CASPAR(4,NCASE) = ALT
          CASPAR(5,NCASE) = RHO
          CASPAR(6,NCASE) = RMU
          CASPAR(7,NCASE) = VSO
          CASPAR(8,NCASE) = 999.
          CASPAR(9,NCASE) = 999.
          CASPAR(10,NCASE)= 999.
          CASPAR(11,NCASE)= 999.
        ENDDO
C
      ELSEIF(ITYPE.EQ.3) THEN
        KCASE = 3
        BET1 = R1
        BET2 = R2
        DBET = R3
        IF(BET1.EQ.BET2) RETURN
        DBET = SIGN(DBET,BET2-BET1)
        NP = 1
        IF(DBET .NE. 0.0) NP = INT((BET2-BET1)/DBET + 0.5) + 1
        IF(NP.LE.0) RETURN
        YES = .FALSE.
        XANS = 0.
        IF(LPWRVAR) YES = YESIN
        IF(YES) XANS = 100.0
        IF(NCASE+NP .GT. ICASX) THEN
         NP = ICASX - NCASE
        ENDIF
        DO IP=1, NP
          NCASE = NCASE + 1
          BET = BET1 + DBET*FLOAT(IP-1)
          CASPAR(0,NCASE) = XANS + FLOAT(KCASE)
          CASPAR(1,NCASE) = ADV
          CASPAR(2,NCASE) = VEL
          CASPAR(3,NCASE) = BET * PI/180.0
          CASPAR(4,NCASE) = ALT
          CASPAR(5,NCASE) = RHO
          CASPAR(6,NCASE) = RMU
          CASPAR(7,NCASE) = VSO
          CASPAR(8,NCASE) = 999.
          CASPAR(9,NCASE) = 999.
          CASPAR(10,NCASE)= 999.
          CASPAR(11,NCASE)= 999.
        ENDDO
C
      ELSEIF(ITYPE.EQ.4) THEN
        KCASE = 4
        VEL1 = R1
        VEL2 = R2
        DVEL = R3
        IF(VEL1.EQ.VEL2) RETURN
        DVEL = SIGN(DVEL,VEL2-VEL1)
        NP = 1
        IF(DVEL .NE. 0.0) NP = INT((VEL2-VEL1)/DVEL + 0.5) + 1
        IF(NP.LE.0) RETURN
        YES = .FALSE.
        XANS = 0.
        IF(LPWRVAR) YES = YESIN
        IF(YES) XANS = 100.0
        ANS4 = 'FP'
        IF(CSIN) ANS4 = 'CS'
        IF(ANS4.EQ.'CS') THEN
          RPM  = VEL/(RAD*ADV*PI/30.)
          IF(RPMIN.NE.0.0) RPM = RPMIN
          ADV = VEL / (RAD*RPM*PI/30.)
          KCASE = 5
          IF(XANS.NE.100.0) THEN
            IF(PSPEC.LE.0.0 .AND. PTOT.GT.0.0) 
     &         PSPEC = PTOT * (RHO*VEL**3*RAD**2)
            IF(LPIN) PSPEC = PIN
          ENDIF
        ENDIF
        IF(NCASE+NP .GT. ICASX) THEN
         NP = ICASX - NCASE
        ENDIF
        DO IP=1, NP
          NCASE = NCASE + 1
          VVEL  = VEL1 + DVEL*FLOAT(IP-1)
          CASPAR(0,NCASE) = XANS + FLOAT(KCASE)
          CASPAR(1,NCASE) = VVEL * ADV/VEL
          CASPAR(2,NCASE) = VVEL
          CASPAR(3,NCASE) = BETA(II)
          CASPAR(4,NCASE) = ALT
          CASPAR(5,NCASE) = RHO
          CASPAR(6,NCASE) = RMU
          CASPAR(7,NCASE) = VSO
          CASPAR(8,NCASE) = 999.
          CASPAR(9,NCASE) = 999.
          CASPAR(10,NCASE)= 999.
          CASPAR(11,NCASE)= 999.
        ENDDO
      ENDIF
      RETURN
      END ! REFSETCAS



      SUBROUTINE REFCASRUN
C-----------------------------------------------------------------------
C     The computing part of CASPLT (xcasepl.f), plotting removed:
C     solve every stored case, then restore and re-solve the current
C     operating point. Scratch W1-W4 (plot data) are not filled.
C-----------------------------------------------------------------------
      INCLUDE 'XROTOR.INC'
      DIMENSION BETSAV(IX), BET0SAV(IX)
C
      IF(NCASE.EQ.0) RETURN
C
      ADVSAV = ADV
      VELSAV = VEL
      ALTSAV = ALT
      RHOSAV = RHO
      VSOSAV = VSO
      RMUSAV = RMU
      DO I=1, II
        BETSAV(I)  = BETA(I)
        BET0SAV(I) = BETA0(I)
      ENDDO
C
      DO 10 ICASE=1, NCASE
        ITQP = IFIX(CASPAR(0,ICASE)/1000.0)
        IPWR = IFIX(CASPAR(0,ICASE)/100.0) - 10*ITQP
        ITYP = IFIX(CASPAR(0,ICASE)) - 100*IPWR - 1000*ITQP
        ADV = CASPAR(1,ICASE)
        VEL = CASPAR(2,ICASE)
        BET = CASPAR(3,ICASE)
        ALT = CASPAR(4,ICASE)
        RHO = CASPAR(5,ICASE)
        RMU = CASPAR(6,ICASE)
        VSO = CASPAR(7,ICASE)
C
        DELB = BET - BETSAV(II)
        DO I=1, II
          BETA(I)  = BETSAV(I)  + DELB
          BETA0(I) = BET0SAV(I) + DELB
        ENDDO
C
        CONV = .FALSE.
        IF(IPWR.EQ.0) THEN
          IF(ITYP.NE.5) THEN
            CALL APER(4,2,.TRUE.)
           ELSE
            CALL APER(3,1,.TRUE.)
          ENDIF
C
          IF(ITYP.EQ.1) THEN
            CALL APER(4,2,.TRUE.)
          ELSEIF(ITYP.EQ.2) THEN
            IF(ITQP.EQ.1) THEN
             CALL APER(1,1,.TRUE.)
            ELSEIF(ITQP.EQ.2) THEN
             CALL APER(2,1,.TRUE.)
            ELSEIF(ITQP.EQ.3) THEN
             CALL APER(3,1,.TRUE.)
            ELSE
             CALL APER(4,2,.TRUE.)
            ENDIF
          ELSEIF(ITYP.EQ.3) THEN
           CALL APER(4,2,.TRUE.)
          ELSEIF(ITYP.EQ.4) THEN
           CALL APER(4,2,.TRUE.)
          ELSEIF(ITYP.EQ.5) THEN
           CALL APER(3,1,.TRUE.)
          ENDIF
         ELSE
          IF(ITYP.EQ.1) THEN
             CALL APER(5,1,.TRUE.)
           ELSEIF(ITYP.EQ.2) THEN
             CALL APER(5,1,.TRUE.)
           ELSEIF(ITYP.EQ.3) THEN
             CALL APER(5,2,.TRUE.)
           ELSEIF(ITYP.EQ.4) THEN
             CALL APER(5,2,.TRUE.)
           ELSEIF(ITYP.EQ.5) THEN
             CALL APER(5,1,.TRUE.)
          ENDIF
        ENDIF
C
        CASPAR(1,ICASE) = ADV
        CASPAR(3,ICASE) = BETA(II)
        IF(CONV) THEN
         CASPAR(8,ICASE) = PTOT*RHO*VEL**3*RAD**2
         CASPAR(9,ICASE) = TTOT*RHO*VEL**2*RAD**2
         CASPAR(10,ICASE)= QTOT*RHO*VEL**2*RAD**3
         CASPAR(11,ICASE)= TTOT/PTOT
        ELSE
         CASPAR(8,ICASE) = 999.
         CASPAR(9,ICASE) = 999.
         CASPAR(10,ICASE)= 999.
         CASPAR(11,ICASE)= 999.
        ENDIF
   10 CONTINUE
C
      ADV = ADVSAV
      VEL = VELSAV
      ALT = ALTSAV
      VSO = VSOSAV
      RHO = RHOSAV
      RMU = RMUSAV
      DO I=1, II
        BETA(I) = BETSAV(I)
        BETA0(I) = BET0SAV(I)
      ENDDO
      CONV = .FALSE.
      CALL APER(4,2,.TRUE.)
      RETURN
      END ! REFCASRUN



      SUBROUTINE REFADDC
C---- OPER ADDC (xoper.f)
      INCLUDE 'XROTOR.INC'
      IF(NCASE.GE.ICASX) RETURN
      NCASE = NCASE+1
      CASPAR(1,NCASE) = ADV
      CASPAR(2,NCASE) = VEL
      CASPAR(3,NCASE) = BETA(II)
      CASPAR(4,NCASE) = ALT
      CASPAR(5,NCASE) = RHO
      CASPAR(6,NCASE) = RMU
      CASPAR(7,NCASE) = VSO
      CASPAR(8,NCASE) = 999.
      CASPAR(9,NCASE) = 999.
      CASPAR(10,NCASE)= 999.
      CASPAR(11,NCASE)= 999.
      RETURN
      END ! REFADDC



      SUBROUTINE REFCASE(ICASE)
C---- OPER CASE i (xoper.f): rerun stored case ICASE
      INCLUDE 'XROTOR.INC'
      IF(NCASE.LE.0) RETURN
      IF(ICASE.LE.0 .OR. ICASE.GT.NCASE) RETURN
C
      ADV = CASPAR(1,ICASE)
      VEL = CASPAR(2,ICASE)
      BET = CASPAR(3,ICASE)
      ALT = CASPAR(4,ICASE)
      RHO = CASPAR(5,ICASE)
      RMU = CASPAR(6,ICASE)
      VSO = CASPAR(7,ICASE)
      POW = CASPAR(8,ICASE)
      THR = CASPAR(9,ICASE)
      TRQ = CASPAR(10,ICASE)
      EFF = CASPAR(11,ICASE)
C
      DELB = BET - BETA(II)
      DO I=1, II
        BETA(I)  = BETA(I)  + DELB
        BETA0(I) = BETA0(I) + DELB
      ENDDO
C
      CONV = .FALSE.
      CALL APER(4,2,.TRUE.)
C
      IF(CONV) THEN
       CASPAR(8,ICASE) = PTOT*RHO*VEL**3*RAD**2
       CASPAR(9,ICASE) = TTOT*RHO*VEL**2*RAD**2
       CASPAR(10,ICASE)= QTOT*RHO*VEL**2*RAD**3
       CASPAR(11,ICASE)= TTOT/PTOT
      ENDIF
      RETURN
      END ! REFCASE



      SUBROUTINE REFCPUT
C---- OPER CPUT f: case table to file FNAME (new file; no prompt)
      INCLUDE 'XROTOR.INC'
      LU = 12
      OPEN(LU,FILE=FNAME,STATUS='UNKNOWN')
      CALL SHOCAS(LU,NPARX,NCASE,CASPAR,RAD,NAME)
      CLOSE(LU)
      RETURN
      END ! REFCPUT



      SUBROUTINE REFCGET
C---- OPER CGET f: case table from file FNAME
      INCLUDE 'XROTOR.INC'
      LU = 12
      OPEN(LU,FILE=FNAME,STATUS='OLD')
      CALL GETCAS(LU,NPARX,NCASE,CASPAR)
      CLOSE(LU)
      RETURN
      END ! REFCGET



      SUBROUTINE REFSAVE
C---- top-level SAVE to FNAME (must not exist: SAVE prompts otherwise)
      INCLUDE 'XROTOR.INC'
      CHARACTER*80 F
      F = FNAME
      CALL SAVE(F)
      RETURN
      END ! REFSAVE



      SUBROUTINE REFLOAD
C---- top-level LOAD from FNAME (prompts on stdin if the file's
C---- station count differs from II)
      INCLUDE 'XROTOR.INC'
      CHARACTER*80 F
      F = FNAME
      CALL LOAD(F)
      RETURN
      END ! REFLOAD
