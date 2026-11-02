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
