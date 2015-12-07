      SUBROUTINE DG7QSB(B, D, DIHDI, G, IPIV, IPIV1, IPIV2, KA, L, LV,
     1                  P, P0, PC, STEP, TD, TG, V, W, X, X0)
C
C  ***  COMPUTE HEURISTIC BOUNDED NEWTON STEP  ***
C
      INTEGER KA, LV, P, P0, PC
      INTEGER IPIV(P), IPIV1(P), IPIV2(P)
      DOUBLE PRECISION B(2,P), D(P), DIHDI(1), G(P), L(1),
     1                 STEP(P,2), TD(P), TG(P), V(LV), W(P), X0(P), X(P)
C     DIMENSION DIHDI(P*(P+1)/2), L(P*(P+1)/2)
C
      DOUBLE PRECISION DD7TPR
      EXTERNAL DD7TPR,DG7QTS, DS7BQN, DS7IPR,DV7CPY, DV7IPR,
     1         DV7SCP, DV7VMP
C
C  ***  LOCAL VARIABLES  ***
C
      INTEGER K, KB, KINIT, NS, P1, P10
      DOUBLE PRECISION DS0, NRED, PRED, RAD
      DOUBLE PRECISION ZERO
C
C  ***  V SUBSCRIPTS  ***
C
      INTEGER DST0, DSTNRM, GTSTEP, NREDUC, PREDUC, RADIUS
C
C/6
C     DATA DST0/3/, DSTNRM/2/, GTSTEP/4/, NREDUC/6/, PREDUC/7/,
C    1     RADIUS/8/
C/7
      PARAMETER (DST0=3, DSTNRM=2, GTSTEP=4, NREDUC=6, PREDUC=7,
     1           RADIUS=8)
C/
      DATA ZERO/0.D+0/
C
C+++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++++
C
      P1 = PC
      IF (KA .LT. 0) GO TO 10
         NRED = V(NREDUC)
         DS0 = V(DST0)
         GO TO 20
 10   P0 = 0
      KA = -1
C
 20   KINIT = -1
      IF (P0 .EQ. P1) KINIT = KA
      CALL DV7CPY(P, X, X0)
      PRED = ZERO
      RAD = V(RADIUS)
      KB = -1
      V(DSTNRM) = ZERO
      IF (P1 .GT. 0) GO TO 30
         NRED = ZERO
         DS0 = ZERO
         CALL DV7SCP(P, STEP, ZERO)
         GO TO 60
C
 30   CALL DV7CPY(P, TD, D)
      CALL DV7IPR(P, IPIV, TD)
      CALL DV7VMP(P, TG, G, D, -1)
      CALL DV7IPR(P, IPIV, TG)
 40   K = KINIT
      KINIT = -1
      V(RADIUS) = RAD - V(DSTNRM)
      CALL DG7QTS(TD, TG, DIHDI, K, L, P1, STEP, V, W)
      P0 = P1
      IF (KA .GE. 0) GO TO 50
         NRED = V(NREDUC)
         DS0 = V(DST0)
C
 50   KA = K
      V(RADIUS) = RAD
      P10 = P1
      CALL DS7BQN(B, D, STEP(1,2), IPIV, IPIV1, IPIV2, KB, L, LV,
     1            NS, P, P1, STEP, TD, TG, V, W, X, X0)
      IF (NS .GT. 0) CALL DS7IPR(P10, IPIV1, DIHDI)
      PRED = PRED + V(PREDUC)
      IF (NS .NE. 0) P0 = 0
      IF (KB .LE. 0) GO TO 40
C
 60   V(DST0) = DS0
      V(NREDUC) = NRED
      V(PREDUC) = PRED
      V(GTSTEP) = DD7TPR(P, G, STEP)
C
 999  RETURN
C  ***  LAST LINE OF DG7QSB FOLLOWS  ***
      END
