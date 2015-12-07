      SUBROUTINE DL7SQR(N, A, L)
C
C  ***  COMPUTE  A = LOWER TRIANGLE OF  L*(L**T),  WITH BOTH
C  ***  L  AND  A  STORED COMPACTLY BY ROWS.  (BOTH MAY OCCUPY THE
C  ***  SAME STORAGE.
C
C  ***  PARAMETERS  ***
C
      INTEGER N
      DOUBLE PRECISION A(1), L(1)
C     DIMENSION A(N*(N+1)/2), L(N*(N+1)/2)
C
C  ***  LOCAL VARIABLES  ***
C
      INTEGER I, II, IJ, IK, IP1, I0, J, JJ, JK, J0, K, NP1
      DOUBLE PRECISION T
C
      NP1 = N + 1
      I0 = N*(N+1)/2
      DO 30 II = 1, N
         I = NP1 - II
         IP1 = I + 1
         I0 = I0 - I
         J0 = I*(I+1)/2
         DO 20 JJ = 1, I
              J = IP1 - JJ
              J0 = J0 - J
              T = 0.0D0
              DO 10 K = 1, J
                   IK = I0 + K
                   JK = J0 + K
                   T = T + L(IK)*L(JK)
 10                CONTINUE
              IJ = I0 + J
              A(IJ) = T
 20           CONTINUE
 30      CONTINUE
 999  RETURN
      END
