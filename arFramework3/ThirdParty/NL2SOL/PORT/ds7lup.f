      SUBROUTINE DS7LUP(A, COSMIN, P, SIZE, STEP, U, W, WCHMTD, WSCALE,
     1                  Y)
C
C  ***  UPDATE SYMMETRIC  A  SO THAT  A * STEP = Y  ***
C  ***  (LOWER TRIANGLE OF  A  STORED ROWWISE       ***
C
C  ***  PARAMETER DECLARATIONS  ***
C
      INTEGER P
      DOUBLE PRECISION A(1), COSMIN, SIZE, STEP(P), U(P), W(P),
     1                 WCHMTD(P), WSCALE, Y(P)
C     DIMENSION A(P*(P+1)/2)
C
C  ***  LOCAL VARIABLES  ***
C
      INTEGER I, J, K
      DOUBLE PRECISION DENMIN, SDOTWM, T, UI, WI
C
C     ***  CONSTANTS  ***
      DOUBLE PRECISION HALF, ONE, ZERO
C
C  ***  EXTERNAL FUNCTIONS AND SUBROUTINES  ***
C
      DOUBLE PRECISION DD7TPR, DV2NRM
      EXTERNAL DD7TPR, DS7LVM, DV2NRM
C
C/6
C     DATA HALF/0.5D+0/, ONE/1.D+0/, ZERO/0.D+0/
C/7
      PARAMETER (HALF=0.5D+0, ONE=1.D+0, ZERO=0.D+0)
C/
C
C-----------------------------------------------------------------------
C
      SDOTWM = DD7TPR(P, STEP, WCHMTD)
      DENMIN = COSMIN * DV2NRM(P,STEP) * DV2NRM(P,WCHMTD)
      WSCALE = ONE
      IF (DENMIN .NE. ZERO) WSCALE = DMIN1(ONE, DABS(SDOTWM/DENMIN))
      T = ZERO
      IF (SDOTWM .NE. ZERO) T = WSCALE / SDOTWM
      DO 10 I = 1, P
 10      W(I) = T * WCHMTD(I)
      CALL DS7LVM(P, U, A, STEP)
      T = HALF * (SIZE * DD7TPR(P, STEP, U)  -  DD7TPR(P, STEP, Y))
      DO 20 I = 1, P
 20      U(I) = T*W(I) + Y(I) - SIZE*U(I)
C
C  ***  SET  A = A + U*(W**T) + W*(U**T)  ***
C
      K = 1
      DO 40 I = 1, P
         UI = U(I)
         WI = W(I)
         DO 30 J = 1, I
              A(K) = SIZE*A(K) + UI*W(J) + WI*U(J)
              K = K + 1
 30           CONTINUE
 40      CONTINUE
C
 999  RETURN
C  ***  LAST CARD OF DS7LUP FOLLOWS  ***
      END
