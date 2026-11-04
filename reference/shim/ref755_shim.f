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
C---- DESGEN, 7.55 order (see REFDESGENX)
      LOGICAL OKIN, LFIXIA
      CALL REFDESGENX(OKIN,LFIXIA,.FALSE.,0.0D0,0.0D0)
      RETURN
      END ! REFDESGEN



      SUBROUTINE REFDESGENX(OKIN,LFIXIA,LFIXCL,CL1,CL2)
C-----------------------------------------------------------------------
C     DESGEN (xdesi.f) with the Betz-limit prompt answered by OKIN.
C     Returns CONV unchanged if DESGEN returns early.
C     Not in 7.55, both for the bug-fix options:
C     LFIXIA  CALL SETIAERO on the new grid after XWINIT (stale
C             aero-section indices)
C     LFIXCL  re-evaluate the linear CL (CL1 root, CL2 tip) with SETCLD
C             on the new grid, and do the stall check there instead
C-----------------------------------------------------------------------
      INCLUDE 'XROTOR.INC'
      LOGICAL OK, OKIN, LFIXIA, LFIXCL
C
      IF(.NOT.LFIXCL) THEN
      DO I=1, II
        N = 1
        CALL GETAERO(N,XISECT,A0,CLMAX,CLMIN,
     &              DCLDA,DCLDA_STALL,DCL_STALL,
     &              CDMIN,CLDMIN,DCDCL2,CMCON,MCRIT,REREF,REXP)
        IF(CLDES(I).GT.CLMAX .OR. CLDES(I).LT.CLMIN) RETURN
      ENDDO
      ENDIF
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
      IF(LFIXCL) THEN
       CALL SETCLD(CL1,CL2)
       DO I=1, II
        N = 1
        CALL GETAERO(N,XISECT,A0,CLMAX,CLMIN,
     &              DCLDA,DCLDA_STALL,DCL_STALL,
     &              CDMIN,CLDMIN,DCDCL2,CMCON,MCRIT,REREF,REXP)
        IF(CLDES(I).GT.CLMAX .OR. CLDES(I).LT.CLMIN) RETURN
       ENDDO
      ENDIF
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
      END ! REFDESGENX



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



      SUBROUTINE REFNACELL(UCONSIN,LFIXSRC)
C-----------------------------------------------------------------------
C     Top-level NACE (NACELL, xrotor.f) on the file named in COMMON
C     FNAME, with the 'Enter constant Ubody' prompt answered by UCONSIN
C     and the plot (NACPLT) removed.
C     LFIXSRC (not in 7.55): put each source at the midpoint of the
C     interval whose area change it represents.
C-----------------------------------------------------------------------
      INCLUDE 'XROTOR.INC'
      CHARACTER     FNTEMP*80
      LOGICAL LFIXSRC
      PARAMETER (IZX=200)
      DIMENSION ABODY(IZX), ZBODY(IZX), RBODY(IZX), SOURCE(IZX)
C
      GREEK = .FALSE.
      LU = LUTEMP
      IF(RAD.EQ.0.0) RETURN
C
      FNTEMP = FNAME
      CALL LC2UC(FNTEMP)
      IF(FNTEMP(1:5).EQ.'CONST') THEN
        UCONS = 0.0
        VCONS = 0.0
        UCONS = UCONSIN
        DO I=1, II
          UBODY(I) = UCONS
          VBODY(I) = VCONS
        END DO
        ZPROP = 0.0
        AWAKE = 0.0
        NZ = 2
        ZBODY(1) = 0.0
        ABODY(1) = 0.0
        ZBODY(2) = 1.0
        ABODY(2) = 0.0
        GO TO 20
      ELSE
      OPEN(LU,FILE=FNAME,STATUS='OLD',ERR=80)
      READ(LU,*,ERR=90) ZPROP, AWAKE
      IF(AWAKE .LT. 0.0) AWAKE = PI*AWAKE**2
      DO K=1, IZX
        READ(LU,*,END=11,ERR=90) ZBODY(K), ABODY(K)
        IF(ABODY(K) .LT. 0.0) ABODY(K) = PI*ABODY(K)**2
      END DO
   11 NZ = K-1
      CLOSE(LU)
      ENDIF
C
      ZPROP = ZPROP/RAD
      AWAKE = AWAKE/RAD**2
      DO K=1, NZ
        ZBODY(K) = ZBODY(K)/RAD    - ZPROP
        ABODY(K) = ABODY(K)/RAD**2
        RBODY(K) = SQRT( ABODY(K)/PI )
      END DO
C
      SOURCE(1) = 0.0
      DO I=2, NZ
        SOURCE(I) = ABODY(I) - ABODY(I-1)
      END DO
      IF(LFIXSRC) THEN
        DO I=1, NZ-1
          SOURCE(I) = ABODY(I+1) - ABODY(I)
        END DO
      ENDIF
C
      DO I=1, II
        CALL SRCVEL(NZ,ZBODY,SOURCE, 0.0,XI(I),UBODY(I),VBODY(I))
      ENDDO
C
 20   XW0 = SQRT(AWAKE/PI)
      CALL XWINIT
      CONV = .FALSE.
      RETURN
C
   80 DO 85 I=1, II
        UBODY(I) = 0.
   85 CONTINUE
      RETURN
C
   90 DO 95 I=1, II
        UBODY(I) = 0.
   95 CONTINUE
      RETURN
      END ! REFNACELL



      SUBROUTINE REFSAVVEL
C---- top-level VPUT (SAVVEL) to the new file FNAME
      INCLUDE 'XROTOR.INC'
      CHARACTER*80 F
      F = FNAME
      CALL SAVVEL(F)
      RETURN
      END ! REFSAVVEL



      SUBROUTINE REFGETVEL(UWTIN,VWTIN)
C-----------------------------------------------------------------------
C     Top-level VGET (GETVEL, xrotor.f) on the file FNAME, with the
C     axial and tangential weight prompts answered by UWTIN, VWTIN.
C-----------------------------------------------------------------------
      INCLUDE 'XROTOR.INC'
C
      GREEK = .FALSE.
      LU = LUTEMP
      OPEN(LU,FILE=FNAME,STATUS='OLD',ERR=95)
      DO 10 I=1, IX
        READ(LU,*,END=11,ERR=95) W1(I), W2(I), W3(I)
 10   CONTINUE
      I = IX+1
 11   IF(I.LT.3) GO TO 96
      NADD = I-1
      CLOSE(LU)
      DO 20 I=1, NADD
        RADD(I) = W1(I)
        UADD(I) = W2(I)
        VADD(I) = W3(I)
 20   CONTINUE
      UWT = UWTIN
      VWT = VWTIN
      DO 30 I=1, NADD
        UADD(I) = UWT*UADD(I)
        VADD(I) = VWT*VADD(I)
 30   CONTINUE
      CALL SPLINE(UADD,UADDR,RADD,NADD)
      CALL SPLINE(VADD,VADDR,RADD,NADD)
      RETURN
 95   CONTINUE
 96   CONTINUE
      RETURN
      END ! REFGETVEL



      SUBROUTINE REFLDW(N,A)
C---- list-directed WRITE of A(1..N) into COMMON SAVFIL//FNAME (probe)
      INCLUDE 'XROTOR.INC'
      DIMENSION A(*)
      CHARACTER*160 L
      L = ' '
      WRITE(L,*) (A(I),I=1,N)
      SAVFIL = L(1:80)
      FNAME = L(81:160)
      RETURN
      END ! REFLDW



C=======================================================================
C     CROTOR counter-rotation iteration (Philip Carter, Esotec,
C     crotor.f v0.85 in CROTOR 755-ES1.3, GPL v2 or later), transcribed
C     for the RPM command: LOADF/LOADA, STORCR, LOADCR, PUTVEL and the
C     iteration loop (labels 400-520). CROTOR's own COMMON variables
C     (ES XROTOR.INC) live in /CRSHIM/ here. Not transcribed: NAME and
C     the file names (no effect on the numbers), output and plots, and
C     the convergence menu (label 590): a solve that fails returns
C     LOK = .FALSE.
C=======================================================================
      SUBROUTINE REFCRLOAD(ICR)
C---- LFWD (ICR=1) / LAFT (ICR=2) from the file FNAME (LOADF/LOADA)
      INCLUDE 'XROTOR.INC'
      LOGICAL LFWD
      COMMON /CRSHIM/ RADCR(2),RAKECR(2),XI0CR(2),XW0CR(2),
     &  XICR(2,IX),CHCR(2,IX),BETACR(2,IX),BETA0CR(2,IX),
     &  UBODYCR(2,IX),
     &  RADDF(IX),UADDF(IX),VADDF(IX),UADDRF(IX),VADDRF(IX),
     &  RADDA(IX),UADDA(IX),VADDA(IX),UADDRA(IX),VADDRA(IX),
     &  VELCR,IICR(2),NBLDCR(2),NADDF,NADDA,LFWD
      CHARACTER*80 F
C
      F = FNAME
      CALL LOAD(F)
      VEL = VELCR
C
      IF(ICR.EQ.1) THEN
        NADDF = NADD
        DO 100 I=1,NADD
          RADDF(I)  = RADD(I)
          UADDF(I)  = UADD(I)
          VADDF(I)  = VADD(I)
          UADDRF(I) = UADDR(I)
          VADDRF(I) = VADDR(I)
 100    CONTINUE
        CALL CRSTOR(1)
        LFWD = .TRUE.
      ELSE
        NADDA = NADD
        DO 200 I=1,NADD
          RADDA(I)  = RADD(I)
          UADDA(I)  = UADD(I)
          VADDA(I)  = VADD(I)
          UADDRA(I) = UADDR(I)
          VADDRA(I) = VADDR(I)
 200    CONTINUE
        CALL CRSTOR(2)
        LFWD = .FALSE.
      ENDIF
      RETURN
      END ! REFCRLOAD



      SUBROUTINE REFCRVCLR
C---- VCLR: zero the stored and current external slipstreams
      INCLUDE 'XROTOR.INC'
      LOGICAL LFWD
      COMMON /CRSHIM/ RADCR(2),RAKECR(2),XI0CR(2),XW0CR(2),
     &  XICR(2,IX),CHCR(2,IX),BETACR(2,IX),BETA0CR(2,IX),
     &  UBODYCR(2,IX),
     &  RADDF(IX),UADDF(IX),VADDF(IX),UADDRF(IX),VADDRF(IX),
     &  RADDA(IX),UADDA(IX),VADDA(IX),UADDRA(IX),VADDRA(IX),
     &  VELCR,IICR(2),NBLDCR(2),NADDF,NADDA,LFWD
      NADD  = 0
      NADDF = 0
      NADDA = 0
      CONV = .FALSE.
      RETURN
      END ! REFCRVCLR



      SUBROUTINE CRSTOR(J)
C---- STORCR: store blade geometry for the CR iteration
      INCLUDE 'XROTOR.INC'
      LOGICAL LFWD
      COMMON /CRSHIM/ RADCR(2),RAKECR(2),XI0CR(2),XW0CR(2),
     &  XICR(2,IX),CHCR(2,IX),BETACR(2,IX),BETA0CR(2,IX),
     &  UBODYCR(2,IX),
     &  RADDF(IX),UADDF(IX),VADDF(IX),UADDRF(IX),VADDRF(IX),
     &  RADDA(IX),UADDA(IX),VADDA(IX),UADDRA(IX),VADDRA(IX),
     &  VELCR,IICR(2),NBLDCR(2),NADDF,NADDA,LFWD
      RADCR(J)  = RAD
      RAKECR(J) = RAKE
      XI0CR(J)  = XI0
      XW0CR(J)  = XW0
      IICR(J)   = II
      NBLDCR(J) = NBLDS
      DO I=1, II
        XICR(J,I)    = XI(I)
        CHCR(J,I)    = CH(I)
        BETACR(J,I)  = BETA(I)
        BETA0CR(J,I) = BETA0(I)
        UBODYCR(J,I) = UBODY(I)
      END DO
      RETURN
      END ! CRSTOR



      SUBROUTINE CRLOAD(J)
C---- LOADCR: load geometry and slipstream for the CR iteration
      INCLUDE 'XROTOR.INC'
      LOGICAL LFWD
      COMMON /CRSHIM/ RADCR(2),RAKECR(2),XI0CR(2),XW0CR(2),
     &  XICR(2,IX),CHCR(2,IX),BETACR(2,IX),BETA0CR(2,IX),
     &  UBODYCR(2,IX),
     &  RADDF(IX),UADDF(IX),VADDF(IX),UADDRF(IX),VADDRF(IX),
     &  RADDA(IX),UADDA(IX),VADDA(IX),UADDRA(IX),VADDRA(IX),
     &  VELCR,IICR(2),NBLDCR(2),NADDF,NADDA,LFWD
C
      IF(J .EQ. 1) THEN
         K = 2
      ELSE
         K = 1
      ENDIF
C--- Store blade angles of outgoing rotor
      IF(LFWD) THEN
        IF(J.NE.1) THEN
          DO I=1,II
            BETACR(K,I) = BETA(I)
            BETA0CR(K,I)= BETA0(I)
          ENDDO
        ELSE
          DO I=1,II
            BETACR(J,I) = BETA(I)
            BETA0CR(J,I)= BETA0(I)
          ENDDO
        ENDIF
      ELSE
        IF(J.NE.2) THEN
          DO I=1,II
            BETACR(K,I) = BETA(I)
            BETA0CR(K,I)= BETA0(I)
          ENDDO
        ELSE
          DO I=1,II
            BETACR(J,I) = BETA(I)
            BETA0CR(J,I)= BETA0(I)
          ENDDO
        ENDIF
      ENDIF
C--- Read in new geometry
      RAD   = RADCR(J)
      RAKE  = RAKECR(J)
      XI0   = XI0CR(J)
      XW0   = XW0CR(J)
      II    = IICR(J)
      NBLDS = NBLDCR(J)
      DO I=1, II
        XI(I)    = XICR(J,I)
        CH(I)    = CHCR(J,I)
        BETA(I)  = BETACR(J,I)
        BETA0(I) = BETA0CR(J,I)
        UBODY(I) = UBODYCR(J,I)
      END DO
C--- load imposed velocity profiles
      IF(J .EQ. 1) THEN
         LFWD = .TRUE.
         NADD = NADDF
         DO I=1,NADD
           RADD(I)  = RADDF(I)
           UADD(I)  = UADDF(I)
           VADD(I)  = VADDF(I)
           UADDR(I) = UADDRF(I)
           VADDR(I) = VADDRF(I)
         ENDDO
      ELSE
         LFWD = .FALSE.
         NADD = NADDA
         DO I=1,NADD
           RADD(I)  = RADDA(I)
           UADD(I)  = UADDA(I)
           VADD(I)  = VADDA(I)
           UADDR(I) = UADDRA(I)
           VADDR(I) = VADDRA(I)
         ENDDO
      ENDIF
      IINF = II + II/2
      CALL SETIAERO
      RETURN
      END ! CRLOAD



      SUBROUTINE CRPUTV(UWTF,VWTF,UWTA,VWTA)
C---- PUTVEL: weighted slipstream of the current rotor, stored for
C---- the other one
      INCLUDE 'XROTOR.INC'
      LOGICAL LFWD
      COMMON /CRSHIM/ RADCR(2),RAKECR(2),XI0CR(2),XW0CR(2),
     &  XICR(2,IX),CHCR(2,IX),BETACR(2,IX),BETA0CR(2,IX),
     &  UBODYCR(2,IX),
     &  RADDF(IX),UADDF(IX),VADDF(IX),UADDRF(IX),VADDRF(IX),
     &  RADDA(IX),UADDA(IX),VADDA(IX),UADDRA(IX),VADDRA(IX),
     &  VELCR,IICR(2),NBLDCR(2),NADDF,NADDA,LFWD
C
      NADD=II
      NADDA=NADD
      NADDF=NADD
C
      BLDS = FLOAT(NBLDS)
      DO 10 I=1, II
        VT = BLDS*GAM(I)/(4.0*PI*XI(I))
        VA = VT*XI(I)/ADW
        UDUCT     = 0.0
        VADUCT_VA = 1.0
        IF(DUCT) THEN
          UDUCT = URDUCT-1.0
          VADUCT_VA = 2.0*URDUCT
        ENDIF
        VA = VA * VADUCT_VA
        UTOT = 1.0 + UDUCT + UBODY(I)
        CALL UVADD(XI(I),WA,WT)
        CI = XI(I)/ADV - WT  -  VT
        SI = UTOT      + WA  +  VA
        RDIM = XI(I)*RAD
        UDIM =  2.0*(SI - (UTOT      + WA))*VEL
        VDIM = -2.0*(CI - (XI(I)/ADV - WT))*VEL
        IF(LFWD) THEN
          RADDA(I)=RDIM
          UADDA(I)=UDIM * UWTA
          VADDA(I)=VDIM * VWTA
        ELSE
          RADDF(I)=RDIM
          UADDF(I)=UDIM * UWTF
          VADDF(I)=VDIM * VWTF
        ENDIF
 10   CONTINUE
      IF(LFWD) THEN
        CALL SPLINE(UADDA,UADDRA,RADDA,NADDA)
        CALL SPLINE(VADDA,VADDRA,RADDA,NADDA)
      ELSE
        CALL SPLINE(UADDF,UADDRF,RADDF,NADDF)
        CALL SPLINE(VADDF,VADDRF,RADDF,NADDF)
      ENDIF
      RETURN
      END ! CRPUTV



      SUBROUTINE REFCRRPM(RPMF,RRAT,VELIN,ICRITL,CRCON,
     &                    UWTF,VWTF,UWTA,VWTA,
     &                    NSOLV,TDHIST,LCONO,LOK,GAMF,TTF,PTF)
C-----------------------------------------------------------------------
C     INPU (rpm and flight speed; synchronized rpm, ratio RRAT), then
C     the RPM command (crotor.f label 400, COMAND = 'RPM ', no input).
C     Outputs: NSOLV rotor solves (ICRITR), TDHIST(k) the thrust change
C     DEL_THRUST/THRUST after solve k, LCONO = LCON (system converged),
C     LOK = every solve converged, and the forward rotor's final GAM,
C     TTOT, PTOT (COMMON holds the aft rotor at the end).
C-----------------------------------------------------------------------
      INCLUDE 'XROTOR.INC'
      LOGICAL LFWD
      COMMON /CRSHIM/ RADCR(2),RAKECR(2),XI0CR(2),XW0CR(2),
     &  XICR(2,IX),CHCR(2,IX),BETACR(2,IX),BETA0CR(2,IX),
     &  UBODYCR(2,IX),
     &  RADDF(IX),UADDF(IX),VADDF(IX),UADDRF(IX),VADDRF(IX),
     &  RADDA(IX),UADDA(IX),VADDA(IX),UADDRA(IX),VADDRA(IX),
     &  VELCR,IICR(2),NBLDCR(2),NADDF,NADDA,LFWD
      DIMENSION TDHIST(*), GAMF(IX)
      LOGICAL LCONO, LOK, LCONF, LCONA, LCON, LNCON, LFORWARD
C
C---- INPU (label 100), rpm and speed part
      RPMA = RPMF * RRAT
      VELOLD = VEL
      VELCR = VELIN
      VEL = VELCR
      THR  = TTOT *(RHO*VELOLD**2*RAD**2)
      TTOT = THR / (RHO*VEL**2*RAD**2)
      TRQ  = QTOT *(RHO*VELOLD**2*RAD**3)
      QTOT = TRQ / (RHO*VEL**2*RAD**3)
      PWR  = PTOT *(RHO*VELOLD**3*RAD**2)
      PTOT = PWR / (RHO*VEL**3*RAD**2)
C
C---- label 400
      LCONF  = .FALSE.
      LCONA  = .FALSE.
      LCON   = .FALSE.
      LNCON  = .FALSE.
      TTOTF  = 0.0
      TTOTA  = 0.0
      ICRITR = 0
      RPMINF = RPMF
      RPMINA = RPMA
      LOK = .TRUE.
C
      LFORWARD=.FALSE.
 410  ICRITR = ICRITR + 1
      LFORWARD = .NOT.LFORWARD
      IF(LFORWARD) THEN
         CALL CRLOAD(1)
      ELSE
         CALL CRLOAD(2)
      ENDIF
C
C---- label 440: RPM, constant pitch and rpm
      IF(LFWD)THEN
         RPM = RPMINF
      ELSE
         RPM = RPMINA
      ENDIF
      ADV = VEL / (RAD*RPM*PI/30.)
      CONV = .FALSE.
      CALL APER(4,2,LOPRINI)
      IF(.NOT.CONV) THEN
        LOK = .FALSE.
        GO TO 600
      ENDIF
C
      CALL CRPUTV(UWTF,VWTF,UWTA,VWTA)
C
      IF(LFWD) THEN
         TDIFF = ABS((TTOT-TTOTF)/TTOT)
         TTOTF = TTOT
         DO I=1, II
           GAMF(I) = GAM(I)
         ENDDO
         TTF = TTOT
         PTF = PTOT
      ELSE
         TDIFF = ABS((TTOT-TTOTA)/TTOT)
         TTOTA = TTOT
      END IF
      TDHIST(ICRITR) = TDIFF
C
      IF(TDIFF .LE. CRCON) THEN
         IF (LFWD) THEN
            LCONF = .TRUE.
         ELSE
            LCONA = .TRUE.
         END IF
      END IF
C
      IF(LCONF .AND. LCONA) GO TO 500
      IF(ICRITR .GE. ICRITL) GO TO 510
      GO TO 410
C
 500  IF (.NOT.LFWD .AND. .NOT.LCON)  GO TO 410
      LCON  = .TRUE.
      GO TO 520
C
 510  IF (.NOT.LFWD .AND. .NOT.LNCON) GO TO 410
      LNCON = .TRUE.
C
 520  IF(LFWD) GO TO 410
C
 600  NSOLV = ICRITR
      LCONO = LCON
      RETURN
      END ! REFCRRPM



      SUBROUTINE REFCRPOWE(POWF,POWA,RPMF,RRAT,VELIN,LPITCH,LSRPM,
     &                     ICRITL,CRCON,UWTF,VWTF,UWTA,VWTA,
     &                     NSOLV,TDHIST,LCONO,LOK,GAMF,TTF,PTF)
C-----------------------------------------------------------------------
C     INPU (rpm and flight speed; the aft rpm is RPMF*RRAT) and the POWE
C     command (crotor.f label 400, COMAND = 'POWE', no input) with the
C     power of each rotor POWF, POWA (POWINF, POWINA; INPU's POWERT/2
C     each by default). LPITCH (FIX): blade angles fixed; LSRPM (SYNC):
C     rpms synchronized at RRAT. Outputs as REFCRRPM. A failed solve
C     restores its blade angles (label 590) and returns LOK = .FALSE.;
C     the convergence menu and AUTO are not transcribed.
C-----------------------------------------------------------------------
      INCLUDE 'XROTOR.INC'
      LOGICAL LFWD
      COMMON /CRSHIM/ RADCR(2),RAKECR(2),XI0CR(2),XW0CR(2),
     &  XICR(2,IX),CHCR(2,IX),BETACR(2,IX),BETA0CR(2,IX),
     &  UBODYCR(2,IX),
     &  RADDF(IX),UADDF(IX),VADDF(IX),UADDRF(IX),VADDRF(IX),
     &  RADDA(IX),UADDA(IX),VADDA(IX),UADDRA(IX),VADDRA(IX),
     &  VELCR,IICR(2),NBLDCR(2),NADDF,NADDA,LFWD
      DIMENSION TDHIST(*), GAMF(IX)
      LOGICAL LCONO, LOK, LCONF, LCONA, LCON, LNCON, LFORWARD
      LOGICAL LPITCH, LSRPM, LSWIT
C
C---- INPU (label 100), rpm and speed part
      RPMA = RPMF * RRAT
      VELOLD = VEL
      VELCR = VELIN
      VEL = VELCR
      THR  = TTOT *(RHO*VELOLD**2*RAD**2)
      TTOT = THR / (RHO*VEL**2*RAD**2)
      TRQ  = QTOT *(RHO*VELOLD**2*RAD**3)
      QTOT = TRQ / (RHO*VEL**2*RAD**3)
      PWR  = PTOT *(RHO*VELOLD**3*RAD**2)
      PTOT = PWR / (RHO*VEL**3*RAD**2)
      RPMFWD = 0.0
      RPMAFT = 0.0
C
C---- label 400
      LCONF  = .FALSE.
      LCONA  = .FALSE.
      LCON   = .FALSE.
      LNCON  = .FALSE.
      LSWIT  = .FALSE.
      TTOTF  = 0.0
      TTOTA  = 0.0
      ICRITR = 0
      POWINF = POWF
      POWINA = POWA
      RPMINF = RPMF
      RPMINA = RPMA
      LOK = .TRUE.
C
      LFORWARD=.FALSE.
 410  ICRITR = ICRITR + 1
      LFORWARD = .NOT.LFORWARD
      IF(LFORWARD) THEN
         CALL CRLOAD(1)
      ELSE
         CALL CRLOAD(2)
      ENDIF
C
C---- label 430: POWE, specified power
      IF(LFWD) THEN
         PSPEC = POWINF
         RPM   = RPMINF
      ELSE
         PSPEC = POWINA
         RPM   = RPMINA
      ENDIF
C
      CONV = .FALSE.
      IF(.NOT.LPITCH) THEN
         ADV = VEL/(RAD*RPM*PI/30.0)
         CALL APER(3,1,LOPRINI)
      ELSE
         IF(.NOT.LSRPM) THEN
           CALL APER(3,2,LOPRINI)
           IF(LFWD) THEN
             RPMFWD = VEL/(RAD*ADV*PI/30.)
           ELSE
             RPMAFT = VEL/(RAD*ADV*PI/30.)
           ENDIF
         ELSE
           IF((LCONA.OR.LCONF) .AND. .NOT.LSWIT) THEN
             RPMTF = (RPMFWD+RPMAFT)/(1.+RRAT)
             RPMTA =  RPMTF*RRAT
             LSWIT = .TRUE.
             LCONF = .FALSE.
             LCONA = .FALSE.
             LCON  = .FALSE.
           ENDIF
           IF(LSWIT) THEN
             IF(LFWD)THEN
               RPM = RPMTF
             ELSE
               RPM = RPMTA
             ENDIF
             ADV = VEL / (RAD*RPM*PI/30.)
             CONV = .FALSE.
             CALL APER(4,2,LOPRINI)
           ELSE
             CALL APER(3,2,LOPRINI)
             IF(LFWD) THEN
               RPMFWD = VEL/(RAD*ADV*PI/30.)
             ELSE
               RPMAFT = VEL/(RAD*ADV*PI/30.)
             ENDIF
           ENDIF
         ENDIF
      ENDIF
C
      IF(CONV) THEN
        IF(.NOT.LPITCH) THEN
          IF(ABS(DBETA).GT.PI/15.0 .AND. TTOT.LT.0.0) GO TO 590
        ENDIF
        GO TO 445
      ENDIF
      GO TO 590
C
 445  CALL CRPUTV(UWTF,VWTF,UWTA,VWTA)
C
      IF(LFWD) THEN
         TDIFF = ABS((TTOT-TTOTF)/TTOT)
         TTOTF = TTOT
         DO I=1, II
           GAMF(I) = GAM(I)
         ENDDO
         TTF = TTOT
         PTF = PTOT
      ELSE
         TDIFF = ABS((TTOT-TTOTA)/TTOT)
         TTOTA = TTOT
      END IF
      TDHIST(ICRITR) = TDIFF
C
      IF(TDIFF .LE. CRCON) THEN
         IF (LFWD) THEN
            LCONF = .TRUE.
         ELSE
            LCONA = .TRUE.
         END IF
      END IF
C
      IF(LCONF .AND. LCONA) THEN
        IF(LPITCH.AND.LSRPM.AND..NOT.LSWIT) THEN
          GO TO 410
        ELSE
          GO TO 500
        ENDIF
      ENDIF
      IF(ICRITR .GE. ICRITL) GO TO 510
      GO TO 410
C
 500  IF (.NOT.LFWD .AND. .NOT.LCON)  GO TO 410
      LCON  = .TRUE.
      GO TO 520
C
 510  IF (.NOT.LFWD .AND. .NOT.LNCON) GO TO 410
      LNCON = .TRUE.
C
 520  IF(LFWD) GO TO 410
      GO TO 600
C
C---- label 590: convergence failed, return blade angles
 590  DO I=1, II
        BETA(I)  = BETA(I)  - DBETA
        BETA0(I) = BETA0(I) - DBETA
      ENDDO
      LOK = .FALSE.
C
 600  NSOLV = ICRITR
      LCONO = LCON
      RETURN
      END ! REFCRPOWE
