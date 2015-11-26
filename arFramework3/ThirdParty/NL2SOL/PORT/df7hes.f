      SUBROUTINE DF7HES(D, G, IRT, IV, LIV, LV, P, V, X)
C
C  ***  COMPUTE FINITE-DIFFERENCE HESSIAN, STORE IT IN V STARTING
C  ***  AT V(IV(FDH)) = V(-IV(H)).
C
C  ***  IF IV(COVREQ) .GE. 0 THEN DF7HES USES GRADIENT DIFFERENCES,
C  ***  OTHERWISE FUNCTION DIFFERENCES.  STORAGE IN V IS AS IN DG7LIT.
C
C IRT VALUES...
C     1 = COMPUTE FUNCTION VALUE, I.E., V(F).
C     2 = COMPUTE G.
C     3 = DONE.
C
C
C  ***  PARAMETER DECLARATIONS  ***
C
      INTEGER IRT, LIV, LV, P
      INTEGER IV(LIV)
      DOUBLE PRECISION D(P), G(P), V(LV), X(P)
C
C  ***  LOCAL VARIABLES  ***
C
      INTEGER GSAVE1, HES, HMI, HPI, HPM, I, K, KIND, L, M, MM1, MM1O2,
     1        PP1O2, STPI, STPM, STP0
      DOUBLE PRECISION DEL, HALF, NEGPT5, ONE, TWO, ZERO
C
C  ***  EXTERNAL SUBROUTINES  ***
C
      EXTERNAL DV7CPY
C
C DV7CPY.... COPY ONE VECTOR TO ANOTHER.
C
C  ***  SUBSCRIPTS FOR IV AND V  ***
C
      INTEGER COVREQ, DELTA, DELTA0, DLTFDC, F, FDH, FX, H, KAGQT, MODE,
     1        NFGCAL, SAVEI, SWITCH, TOOBIG, W, XMSAVE
C
C/6
C     DATA HALF/0.5D+0/, NEGPT5/-0.5D+0/, ONE/1.D+0/, TWO/2.D+0/,
C    1     ZERO/0.D+0/
C/7
      PARAMETER (HALF=0.5D+0, NEGPT5=-0.5D+0, ONE=1.D+0, TWO=2.D+0,
     1     ZERO=0.D+0)
C/
C
C/6
C     DATA COVREQ/15/, DELTA/52/, DELTA0/44/, DLTFDC/42/, F/10/,
C    1     FDH/74/, FX/53/, H/56/, KAGQT/33/, MODE/35/, NFGCAL/7/,
C    2     SAVEI/63/, SWITCH/12/, TOOBIG/2/, W/65/, XMSAVE/51/
C/7
      PARAMETER (COVREQ=15, DELTA=52, DELTA0=44, DLTFDC=42, F=10,
     1           FDH=74, FX=53, H=56, KAGQT=33, MODE=35, NFGCAL=7,
     2           SAVEI=63, SWITCH=12, TOOBIG=2, W=65, XMSAVE=51)
C/
C
C+++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++++
C
      IRT = 4
      KIND = IV(COVREQ)
      M = IV(MODE)
      IF (M .GT. 0) GO TO 10
         IV(H) = -IABS(IV(H))
         IV(FDH) = 0
         IV(KAGQT) = -1
         V(FX) = V(F)
 10   IF (M .GT. P) GO TO 999
      IF (KIND .LT. 0) GO TO 110
C
C  ***  COMPUTE FINITE-DIFFERENCE HESSIAN USING BOTH FUNCTION AND
C  ***  GRADIENT VALUES.
C
      GSAVE1 = IV(W) + P
      IF (M .GT. 0) GO TO 20
C        ***  FIRST CALL ON DF7HES.  SET GSAVE = G, TAKE FIRST STEP  ***
         CALL DV7CPY(P, V(GSAVE1), G)
         IV(SWITCH) = IV(NFGCAL)
         GO TO 90
C
 20   DEL = V(DELTA)
      X(M) = V(XMSAVE)
      IF (IV(TOOBIG) .EQ. 0) GO TO 40
C
C     ***  HANDLE OVERSIZE V(DELTA)  ***
C
         IF (DEL*X(M) .GT. ZERO) GO TO 30
C             ***  WE ALREADY TRIED SHRINKING V(DELTA), SO QUIT  ***
              IV(FDH) = -2
              GO TO 220
C
C        ***  TRY SHRINKING V(DELTA)  ***
 30      DEL = NEGPT5 * DEL
         GO TO 100
C
 40   HES = -IV(H)
C
C  ***  SET  G = (G - GSAVE)/DEL  ***
C
      DO 50 I = 1, P
         G(I) = (G(I) - V(GSAVE1)) / DEL
         GSAVE1 = GSAVE1 + 1
 50      CONTINUE
C
C  ***  ADD G AS NEW COL. TO FINITE-DIFF. HESSIAN MATRIX  ***
C
      K = HES + M*(M-1)/2
      L = K + M - 2
      IF (M .EQ. 1) GO TO 70
C
C  ***  SET  H(I,M) = 0.5 * (H(I,M) + G(I))  FOR I = 1 TO M-1  ***
C
      MM1 = M - 1
      DO 60 I = 1, MM1
         V(K) = HALF * (V(K) + G(I))
         K = K + 1
 60      CONTINUE
C
C  ***  ADD  H(I,M) = G(I)  FOR I = M TO P  ***
C
 70   L = L + 1
      DO 80 I = M, P
         V(L) = G(I)
         L = L + I
 80      CONTINUE
C
 90   M = M + 1
      IV(MODE) = M
      IF (M .GT. P) GO TO 210
C
C  ***  CHOOSE NEXT FINITE-DIFFERENCE STEP, RETURN TO GET G THERE  ***
C
      DEL = V(DELTA0) * DMAX1(ONE/D(M), DABS(X(M)))
      IF (X(M) .LT. ZERO) DEL = -DEL
      V(XMSAVE) = X(M)
 100  X(M) = X(M) + DEL
      V(DELTA) = DEL
      IRT = 2
      GO TO 999
C
C  ***  COMPUTE FINITE-DIFFERENCE HESSIAN USING FUNCTION VALUES ONLY.
C
 110  STP0 = IV(W) + P - 1
      MM1 = M - 1
      MM1O2 = M*MM1/2
      IF (M .GT. 0) GO TO 120
C        ***  FIRST CALL ON DF7HES.  ***
         IV(SAVEI) = 0
         GO TO 200
C
 120  I = IV(SAVEI)
      HES = -IV(H)
      IF (I .GT. 0) GO TO 180
      IF (IV(TOOBIG) .EQ. 0) GO TO 140
C
C     ***  HANDLE OVERSIZE STEP  ***
C
         STPM = STP0 + M
         DEL = V(STPM)
         IF (DEL*X(XMSAVE) .GT. ZERO) GO TO 130
C             ***  WE ALREADY TRIED SHRINKING THE STEP, SO QUIT  ***
              IV(FDH) = -2
              GO TO 220
C
C        ***  TRY SHRINKING THE STEP  ***
 130     DEL = NEGPT5 * DEL
         X(M) = X(XMSAVE) + DEL
         V(STPM) = DEL
         IRT = 1
         GO TO 999
C
C  ***  SAVE F(X + STP(M)*E(M)) IN H(P,M)  ***
C
 140  PP1O2 = P * (P-1) / 2
      HPM = HES + PP1O2 + MM1
      V(HPM) = V(F)
C
C  ***  START COMPUTING ROW M OF THE FINITE-DIFFERENCE HESSIAN H.  ***
C
      HMI = HES + MM1O2
      IF (MM1 .EQ. 0) GO TO 160
      HPI = HES + PP1O2
      DO 150 I = 1, MM1
         V(HMI) = V(FX) - (V(F) + V(HPI))
         HMI = HMI + 1
         HPI = HPI + 1
 150     CONTINUE
 160  V(HMI) = V(F) - TWO*V(FX)
C
C  ***  COMPUTE FUNCTION VALUES NEEDED TO COMPLETE ROW M OF H.  ***
C
      I = 1
C
 170  IV(SAVEI) = I
      STPI = STP0 + I
      V(DELTA) = X(I)
      X(I) = X(I) + V(STPI)
      IF (I .EQ. M) X(I) = V(XMSAVE) - V(STPI)
      IRT = 1
      GO TO 999
C
 180  X(I) = V(DELTA)
      IF (IV(TOOBIG) .EQ. 0) GO TO 190
C        ***  PUNT IN THE EVENT OF AN OVERSIZE STEP  ***
         IV(FDH) = -2
         GO TO 220
C
C  ***  FINISH COMPUTING H(M,I)  ***
C
 190  STPI = STP0 + I
      HMI = HES + MM1O2 + I - 1
      STPM = STP0 + M
      V(HMI) = (V(HMI) + V(F)) / (V(STPI)*V(STPM))
      I = I + 1
      IF (I .LE. M) GO TO 170
      IV(SAVEI) = 0
      X(M) = V(XMSAVE)
C
 200  M = M + 1
      IV(MODE) = M
      IF (M .GT. P) GO TO 210
C
C  ***  PREPARE TO COMPUTE ROW M OF THE FINITE-DIFFERENCE HESSIAN H.
C  ***  COMPUTE M-TH STEP SIZE STP(M), THEN RETURN TO OBTAIN
C  ***  F(X + STP(M)*E(M)), WHERE E(M) = M-TH STD. UNIT VECTOR.
C
      DEL = V(DLTFDC) * DMAX1(ONE/D(M), DABS(X(M)))
      IF (X(M) .LT. ZERO) DEL = -DEL
      V(XMSAVE) = X(M)
      X(M) = X(M) + DEL
      STPM = STP0 + M
      V(STPM) = DEL
      IRT = 1
      GO TO 999
C
C  ***  RESTORE V(F), ETC.  ***
C
 210  IV(FDH) = HES
 220  V(F) = V(FX)
      IRT = 3
      IF (KIND .LT. 0) GO TO 999
         IV(NFGCAL) = IV(SWITCH)
         GSAVE1 = IV(W) + P
         CALL DV7CPY(P, G, V(GSAVE1))
         GO TO 999
C
 999  RETURN
C  ***  LAST CARD OF DF7HES FOLLOWS  ***
      END
