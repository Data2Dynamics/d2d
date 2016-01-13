      SUBROUTINE DL7NVR(N, LIN, L)
C
C  ***  COMPUTE  LIN = L**-1,  BOTH  N X N  LOWER TRIANG. STORED   ***
C  ***  COMPACTLY BY ROWS.  LIN AND L MAY SHARE THE SAME STORAGE.  ***
C
C  ***  PARAMETERS  ***
C
      INTEGER N
      DOUBLE PRECISION L(1), LIN(1)
C     DIMENSION L(N*(N+1)/2), LIN(N*(N+1)/2)
C
C  ***  LOCAL VARIABLES  ***
C
      INTEGER I, II, IM1, JJ, J0, J1, K, K0, NP1
      DOUBLE PRECISION ONE, T, ZERO
C/6
C     DATA ONE/1.D+0/, ZERO/0.D+0/
C/7
      PARAMETER (ONE=1.D+0, ZERO=0.D+0)
C/
C
C  ***  BODY  ***
C
      NP1 = N + 1
      J0 = N*(NP1)/2
      DO 30 II = 1, N
         I = NP1 - II
         LIN(J0) = ONE/L(J0)
         IF (I .LE. 1) GO TO 999
         J1 = J0
         IM1 = I - 1
         DO 20 JJ = 1, IM1
              T = ZERO
              J0 = J1
              K0 = J1 - JJ
              DO 10 K = 1, JJ
                   T = T - L(K0)*LIN(J0)
                   J0 = J0 - 1
                   K0 = K0 + K - I
 10                CONTINUE
              LIN(J0) = T/L(K0)
 20           CONTINUE
         J0 = J0 - 1
 30      CONTINUE
 999  RETURN
C  ***  LAST CARD OF DL7NVR FOLLOWS  ***
      END
