      SUBROUTINE DF7DHB(B, D, G, IRT, IV, LIV, LV, P, V, X)
C
C  ***  COMPUTE FINITE-DIFFERENCE HESSIAN, STORE IT IN V STARTING
C  ***  AT V(IV(FDH)) = V(-IV(H)).  HONOR SIMPLE BOUNDS IN B.
C
C  ***  IF IV(COVREQ) .GE. 0 THEN DF7DHB USES GRADIENT DIFFERENCES,
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
      DOUBLE PRECISION B(2,P), D(P), G(P), V(LV), X(P)
C
C  ***  LOCAL VARIABLES  ***
C
      LOGICAL OFFSID
      INTEGER GSAVE1, HES, HMI, HPI, HPM, I, K, KIND, L, M, MM1, MM1O2,
     1        NEWM1, PP1O2, STPI, STPM, STP0
      DOUBLE PRECISION DEL, DEL0, T, XM, XM1
      DOUBLE PRECISION HALF, HLIM, ONE, TWO, ZERO
C
C  ***  EXTERNAL SUBROUTINES  ***
C
      EXTERNAL DV7CPY, DV7SCP
C
C DV7CPY.... COPY ONE VECTOR TO ANOTHER.
C DV7SCP... COPY SCALAR TO ALL COMPONENTS OF A VECTOR.
C
C  ***  SUBSCRIPTS FOR IV AND V  ***
C
      INTEGER COVREQ, DELTA, DELTA0, DLTFDC, F, FDH, FX, H, KAGQT, MODE,
     1        NFGCAL, SAVEI, SWITCH, TOOBIG, W, XMSAVE
C
C/6
C     DATA HALF/0.5D+0/, HLIM/0.1D+0/, ONE/1.D+0/, TWO/2.D+0/,
C    1     ZERO/0.D+0/
C/7
      PARAMETER (HALF=0.5D+0, HLIM=0.1D+0, ONE=1.D+0, TWO=2.D+0,
     1           ZERO=0.D+0)
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
         HES = IABS(IV(H))
         IV(H) = -HES
         IV(FDH) = 0
         IV(KAGQT) = -1
         V(FX) = V(F)
C        *** SUPPLY ZEROS IN CASE B(1,I) = B(2,I) FOR SOME I ***
         CALL DV7SCP(P*(P+1)/2, V(HES), ZERO)
 10   IF (M .GT. P) GO TO 999
      IF (KIND .LT. 0) GO TO 120
C
C  ***  COMPUTE FINITE-DIFFERENCE HESSIAN USING BOTH FUNCTION AND
C  ***  GRADIENT VALUES.
C
      GSAVE1 = IV(W) + P
      IF (M .GT. 0) GO TO 20
C        ***  FIRST CALL ON DF7DHB.  SET GSAVE = G, TAKE FIRST STEP  ***
         CALL DV7CPY(P, V(GSAVE1), G)
         IV(SWITCH) = IV(NFGCAL)
         GO TO 80
C
 20   DEL = V(DELTA)
      X(M) = V(XMSAVE)
      IF (IV(TOOBIG) .EQ. 0) GO TO 30
C
C     ***  HANDLE OVERSIZE V(DELTA)  ***
C
         DEL0 = V(DELTA0) * DMAX1(ONE/D(M), DABS(X(M)))
         DEL = HALF * DEL
         IF (DABS(DEL/DEL0) .LE. HLIM) GO TO 140
C
 30   HES = -IV(H)
C
C  ***  SET  G = (G - GSAVE)/DEL  ***
C
      DEL = ONE / DEL
      DO 40 I = 1, P
         G(I) = DEL * (G(I) - V(GSAVE1))
         GSAVE1 = GSAVE1 + 1
 40      CONTINUE
C
C  ***  ADD G AS NEW COL. TO FINITE-DIFF. HESSIAN MATRIX  ***
C
      K = HES + M*(M-1)/2
      L = K + M - 2
      IF (M .EQ. 1) GO TO 60
C
C  ***  SET  H(I,M) = 0.5 * (H(I,M) + G(I))  FOR I = 1 TO M-1  ***
C
      MM1 = M - 1
      DO 50 I = 1, MM1
         IF (B(1,I) .LT. B(2,I)) V(K) = HALF * (V(K) + G(I))
         K = K + 1
 50      CONTINUE
C
C  ***  ADD  H(I,M) = G(I)  FOR I = M TO P  ***
C
 60   L = L + 1
      DO 70 I = M, P
         IF (B(1,I) .LT. B(2,I)) V(L) = G(I)
         L = L + I
 70      CONTINUE
C
 80   M = M + 1
      IV(MODE) = M
      IF (M .GT. P) GO TO 340
      IF (B(1,M) .GE. B(2,M)) GO TO 80
C
C  ***  CHOOSE NEXT FINITE-DIFFERENCE STEP, RETURN TO GET G THERE  ***
C
      DEL = V(DELTA0) * DMAX1(ONE/D(M), DABS(X(M)))
      XM = X(M)
      IF (XM .LT. ZERO) GO TO 90
         XM1 = XM + DEL
         IF (XM1 .LE. B(2,M)) GO TO 110
           XM1 = XM - DEL
           IF (XM1 .GE. B(1,M)) GO TO 100
           GO TO 280
 90    XM1 = XM - DEL
       IF (XM1 .GE. B(1,M)) GO TO 100
       XM1 = XM + DEL
       IF (XM1 .LE. B(2,M)) GO TO 110
       GO TO 280
C
 100  DEL = -DEL
 110  V(XMSAVE) = XM
      X(M) = XM1
      V(DELTA) = DEL
      IRT = 2
      GO TO 999
C
C  ***  COMPUTE FINITE-DIFFERENCE HESSIAN USING FUNCTION VALUES ONLY.
C
 120  STP0 = IV(W) + P - 1
      MM1 = M - 1
      MM1O2 = M*MM1/2
      HES = -IV(H)
      IF (M .GT. 0) GO TO 130
C        ***  FIRST CALL ON DF7DHB.  ***
         IV(SAVEI) = 0
         GO TO 240
C
 130  IF (IV(TOOBIG) .EQ. 0) GO TO 150
C        ***  PUNT IN THE EVENT OF AN OVERSIZE STEP  ***
 140     IV(FDH) = -2
         GO TO 350
 150  I = IV(SAVEI)
      IF (I .GT. 0) GO TO 190
C
C  ***  SAVE F(X + STP(M)*E(M)) IN H(P,M)  ***
C
      PP1O2 = P * (P-1) / 2
      HPM = HES + PP1O2 + MM1
      V(HPM) = V(F)
C
C  ***  START COMPUTING ROW M OF THE FINITE-DIFFERENCE HESSIAN H.  ***
C
      NEWM1 = 1
      GO TO 260
 160  HMI = HES + MM1O2
      IF (MM1 .EQ. 0) GO TO 180
      HPI = HES + PP1O2
      DO 170 I = 1, MM1
         T = ZERO
         IF (B(1,I) .LT. B(2,I)) T = V(FX) - (V(F) + V(HPI))
         V(HMI) = T
         HMI = HMI + 1
         HPI = HPI + 1
 170     CONTINUE
 180  V(HMI) = V(F) - TWO*V(FX)
      IF (OFFSID) V(HMI) = V(FX) - TWO*V(F)
C
C  ***  COMPUTE FUNCTION VALUES NEEDED TO COMPLETE ROW M OF H.  ***
C
      I = 0
      GO TO 200
C
 190  X(I) = V(DELTA)
C
C  ***  FINISH COMPUTING H(M,I)  ***
C
      STPI = STP0 + I
      HMI = HES + MM1O2 + I - 1
      STPM = STP0 + M
      V(HMI) = (V(HMI) + V(F)) / (V(STPI)*V(STPM))
 200  I = I + 1
      IF (I .GT. M) GO TO 230
         IF (B(1,I) .LT. B(2,I)) GO TO 210
         GO TO 200
C
 210  IV(SAVEI) = I
      STPI = STP0 + I
      V(DELTA) = X(I)
      X(I) = X(I) + V(STPI)
      IRT = 1
      IF (I .LT. M) GO TO 999
      NEWM1 = 2
      GO TO 260
 220  X(M) = V(XMSAVE) - DEL
      IF (OFFSID) X(M) = V(XMSAVE) + TWO*DEL
      GO TO 999
C
 230  IV(SAVEI) = 0
      X(M) = V(XMSAVE)
C
 240  M = M + 1
      IV(MODE) = M
      IF (M .GT. P) GO TO 330
      IF (B(1,M) .LT. B(2,M)) GO TO 250
      GO TO 240
C
C  ***  PREPARE TO COMPUTE ROW M OF THE FINITE-DIFFERENCE HESSIAN H.
C  ***  COMPUTE M-TH STEP SIZE STP(M), THEN RETURN TO OBTAIN
C  ***  F(X + STP(M)*E(M)), WHERE E(M) = M-TH STD. UNIT VECTOR.
C
 250  V(XMSAVE) = X(M)
      NEWM1 = 3
 260  XM = V(XMSAVE)
      DEL = V(DLTFDC) * DMAX1(ONE/D(M), DABS(XM))
      XM1 = XM + DEL
      OFFSID = .FALSE.
      IF (XM1 .LE. B(2,M)) GO TO 270
         OFFSID = .TRUE.
         XM1 = XM - DEL
         IF (XM - TWO*DEL .GE. B(1,M)) GO TO 300
         GO TO 280
 270   IF (XM-DEL .GE. B(1,M)) GO TO 290
       OFFSID = .TRUE.
       IF (XM + TWO*DEL .LE. B(2,M)) GO TO 310
C
 280  IV(FDH) = -2
      GO TO 350
C
 290  IF (XM .GE. ZERO) GO TO 310
      XM1 = XM - DEL
 300  DEL = -DEL
 310  GO TO (160, 220, 320), NEWM1
 320  X(M) = XM1
      STPM = STP0 + M
      V(STPM) = DEL
      IRT = 1
      GO TO 999
C
C  ***  HANDLE SPECIAL CASE OF B(1,P) = B(2,P) -- CLEAR SCRATCH VALUES
C  ***  FROM LAST ROW OF FDH...
C
 330  IF (B(1,P) .LT. B(2,P)) GO TO 340
         I = HES + P*(P-1)/2
         CALL DV7SCP(P, V(I), ZERO)
C
C  ***  RESTORE V(F), ETC.  ***
C
 340  IV(FDH) = HES
 350  V(F) = V(FX)
      IRT = 3
      IF (KIND .LT. 0) GO TO 999
         IV(NFGCAL) = IV(SWITCH)
         GSAVE1 = IV(W) + P
         CALL DV7CPY(P, G, V(GSAVE1))
         GO TO 999
C
 999  RETURN
C  ***  LAST LINE OF DF7DHB FOLLOWS  ***
      END
