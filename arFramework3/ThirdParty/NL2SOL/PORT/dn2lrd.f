      SUBROUTINE DN2LRD(DR, IV, L, LH, LIV, LV, ND, NN, P, R, RD, V)
C
C  ***  COMPUTE REGRESSION DIAGNOSTIC AND DEFAULT COVARIANCE MATRIX FOR
C        DRN2G  ***
C
C  ***  PARAMETERS  ***
C
      INTEGER LH, LIV, LV, ND, NN, P
      INTEGER IV(LIV)
      DOUBLE PRECISION DR(ND,P), L(LH), R(NN), RD(NN), V(LV)
C
C  ***  CODED BY DAVID M. GAY (WINTER 1982, FALL 1983)  ***
C
C  ***  EXTERNAL FUNCTIONS AND SUBROUTINES  ***
C
      DOUBLE PRECISION DD7TPR
      EXTERNAL DD7TPR, DL7ITV, DL7IVM,DO7PRD, DV7SCP
C
C  ***  LOCAL VARIABLES  ***
C
      INTEGER COV, I, J, M, STEP1
      DOUBLE PRECISION A, FF, S, T
C
C  ***  CONSTANTS  ***
C
      DOUBLE PRECISION NEGONE, ONE, ONEV(1), ZERO
C
C  ***  INTRINSIC FUNCTIONS  ***
C/+
      DOUBLE PRECISION DSQRT
C/
C
C  ***  IV AND V SUBSCRIPTS  ***
C
      INTEGER F, H, MODE, RDREQ, STEP
C/6
C     DATA F/10/, H/56/, MODE/35/, RDREQ/57/, STEP/40/
C/7
      PARAMETER (F=10, H=56, MODE=35, RDREQ=57, STEP=40)
C/
C/6
C     DATA NEGONE/-1.D+0/, ONE/1.D+0/, ZERO/0.D+0/
C/7
      PARAMETER (NEGONE=-1.D+0, ONE=1.D+0, ZERO=0.D+0)
C/
      DATA ONEV(1)/1.D+0/
C
C++++++++++++++++++++++++++++++++  BODY  +++++++++++++++++++++++++++++++
C
      STEP1 = IV(STEP)
      I = IV(RDREQ)
      IF (I .LE. 0) GO TO 999
      IF (MOD(I,4) .LT. 2) GO TO 30
      FF = ONE
      IF (V(F) .NE. ZERO) FF = ONE / DSQRT(DABS(V(F)))
      CALL DV7SCP(NN, RD, NEGONE)
      DO 20 I = 1, NN
         A = R(I)**2
         M = STEP1
         DO 10 J = 1, P
            V(M) = DR(I,J)
            M = M + 1
 10         CONTINUE
         CALL DL7IVM(P, V(STEP1), L, V(STEP1))
         S = DD7TPR(P, V(STEP1), V(STEP1))
         T = ONE - S
         IF (T .LE. ZERO) GO TO 20
         A = A * S / T
         RD(I) = DSQRT(A) * FF
 20      CONTINUE
C
 30   IF (IV(MODE) - P .LT. 2) GO TO 999
C
C  ***  COMPUTE DEFAULT COVARIANCE MATRIX  ***
C
      COV = IABS(IV(H))
      DO 50 I = 1, NN
         M = STEP1
         DO 40 J = 1, P
            V(M) = DR(I,J)
            M = M + 1
 40         CONTINUE
         CALL DL7IVM(P, V(STEP1), L, V(STEP1))
         CALL DL7ITV(P, V(STEP1), L, V(STEP1))
         CALL DO7PRD(1, LH, P, V(COV), ONEV, V(STEP1), V(STEP1))
 50      CONTINUE
C
 999  RETURN
C  ***  LAST LINE OF DN2LRD FOLLOWS  ***
      END
