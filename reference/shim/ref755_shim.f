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
