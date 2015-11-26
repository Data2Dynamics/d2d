      SUBROUTINE DS7BQN(B, D, DST, IPIV, IPIV1, IPIV2, KB, L, LV, NS,
     1                  P, P1, STEP, TD, TG, V, W, X, X0)
C
C  ***  COMPUTE BOUNDED MODIFIED NEWTON STEP  ***
C
      INTEGER KB, LV, NS, P, P1
      INTEGER IPIV(P), IPIV1(P), IPIV2(P)
      DOUBLE PRECISION B(2,P), D(P), DST(P), L(1),
     1                 STEP(P), TD(P), TG(P), V(LV), W(P), X(P),
     2                 X0(P)
C     DIMENSION L(P*(P+1)/2)
C
      DOUBLE PRECISION DD7TPR, DR7MDC, DV2NRM
      EXTERNAL DD7TPR, I7SHFT, DL7ITV, DL7IVM, DQ7RSH, DR7MDC, DV2NRM,
     1        DV2AXY,DV7CPY, DV7IPR, DV7SCP, DV7SHF
C
C  ***  LOCAL VARIABLES  ***
C
      INTEGER I, J, K, P0, P1M1
      DOUBLE PRECISION ALPHA, DST0, DST1, DSTMAX, DSTMIN, DX, GTS, T,
     1                 TI, T1, XI
      DOUBLE PRECISION FUDGE, HALF, MEPS2, ONE, TWO, ZERO
C
C  ***  V SUBSCRIPTS  ***
C
      INTEGER DSTNRM, GTSTEP, PHMNFC, PHMXFC, PREDUC, RADIUS, STPPAR
C
C/6
C     DATA DSTNRM/2/, GTSTEP/4/, PHMNFC/20/, PHMXFC/21/, PREDUC/7/,
C    1     RADIUS/8/, STPPAR/5/
C/7
      PARAMETER (DSTNRM=2, GTSTEP=4, PHMNFC=20, PHMXFC=21, PREDUC=7,
     1           RADIUS=8, STPPAR=5)
      SAVE MEPS2
C/
C
      DATA FUDGE/1.0001D+0/, HALF/0.5D+0/, MEPS2/0.D+0/,
     1     ONE/1.0D+0/, TWO/2.D+0/, ZERO/0.D+0/
C
C+++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++++
C
      DSTMAX = FUDGE * (ONE + V(PHMXFC)) * V(RADIUS)
      DSTMIN = (ONE + V(PHMNFC)) * V(RADIUS)
      DST1 = ZERO
      IF (MEPS2 .LE. ZERO) MEPS2 = TWO * DR7MDC(3)
      P0 = P1
      NS = 0
      DO 10 I = 1, P
         IPIV1(I) = I
         IPIV2(I) = I
 10      CONTINUE
      DO 20 I = 1, P1
 20      W(I) = -STEP(I) * TD(I)
      ALPHA = DABS(V(STPPAR))
      V(PREDUC) = ZERO
      GTS = -V(GTSTEP)
      IF (KB .LT. 0) CALL DV7SCP(P, DST, ZERO)
      KB = 1
C
C     ***  -W = D TIMES RESTRICTED NEWTON STEP FROM X + DST/D.
C
C     ***  FIND T SUCH THAT X - T*W IS STILL FEASIBLE.
C
 30   T = ONE
      K = 0
      DO 60 I = 1, P1
         J = IPIV(I)
         DX = W(I) / D(J)
         XI = X(J) - DX
         IF (XI .LT. B(1,J)) GO TO 40
         IF (XI .LE. B(2,J)) GO TO 60
              TI = ( X(J)  -  B(2,J) ) / DX
              K = I
              GO TO 50
 40      TI = ( X(J)  -  B(1,J) ) / DX
              K = -I
 50      IF (T .LE. TI) GO TO 60
              T = TI
 60      CONTINUE
C
      IF (P .GT. P1) CALL DV7CPY(P-P1, STEP(P1+1), DST(P1+1))
      CALL DV2AXY(P1, STEP, -T, W, DST)
      DST0 = DST1
      DST1 = DV2NRM(P, STEP)
C
C  ***  CHECK FOR OVERSIZE STEP  ***
C
      IF (DST1 .LE. DSTMAX) GO TO 80
      IF (P1 .GE. P0) GO TO 70
         IF (DST0 .LT. DSTMIN) KB = 0
         GO TO 110
C
 70   K = 0
C
C  ***  UPDATE DST, TG, AND V(PREDUC)  ***
C
 80   V(DSTNRM) = DST1
      CALL DV7CPY(P1, DST, STEP)
      T1 = ONE - T
      DO 90 I = 1, P1
 90      TG(I) = T1 * TG(I)
      IF (ALPHA .GT. ZERO) CALL DV2AXY(P1, TG, T*ALPHA, W, TG)
      V(PREDUC) = V(PREDUC) + T*((ONE - HALF*T)*GTS +
     1                        HALF*ALPHA*T*DD7TPR(P1,W,W))
      IF (K .EQ. 0) GO TO 110
C
C     ***  PERMUTE L, ETC. IF NECESSARY  ***
C
      P1M1 = P1 - 1
      J = IABS(K)
      IF (J .EQ. P1) GO TO 100
         NS = NS + 1
         IPIV2(P1) = J
         CALL DQ7RSH(J, P1, .FALSE., TG, L, W)
         CALL I7SHFT(P1, J, IPIV)
         CALL I7SHFT(P1, J, IPIV1)
         CALL DV7SHF(P1, J, TG)
         CALL DV7SHF(P1, J, DST)
 100  IF (K .LT. 0) IPIV(P1) = -IPIV(P1)
      P1 = P1M1
      IF (P1 .LE. 0) GO TO 110
      CALL DL7IVM(P1, W, L, TG)
      GTS = DD7TPR(P1, W, W)
      CALL DL7ITV(P1, W, L, W)
      GO TO 30
C
C     ***  UNSCALE STEP  ***
C
 110  DO 120 I = 1, P
         J = IABS(IPIV(I))
         STEP(J) = DST(I) / D(J)
 120     CONTINUE
C
C  ***  FUDGE STEP TO ENSURE THAT IT FORCES APPROPRIATE COMPONENTS
C  ***  TO THEIR BOUNDS  ***
C
      IF (P1 .GE. P0) GO TO 150
      K = P1 + 1
      DO 140 I = K, P0
         J = IPIV(I)
         T = MEPS2
         IF (J .GT. 0) GO TO 130
            T = -T
            J = -J
            IPIV(I) = J
 130     T = T * DMAX1(DABS(X(J)), DABS(X0(J)))
         STEP(J) = STEP(J) + T
 140     CONTINUE
C
 150  CALL DV2AXY(P, X, ONE, STEP, X0)
      IF (NS .GT. 0) CALL DV7IPR(P0, IPIV1, TD)
 999  RETURN
C  ***  LAST LINE OF DS7BQN FOLLOWS  ***
      END
