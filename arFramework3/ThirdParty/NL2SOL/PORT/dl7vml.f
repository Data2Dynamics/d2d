      SUBROUTINE DL7VML(N, X, L, Y)
C
C  ***  COMPUTE  X = L*Y, WHERE  L  IS AN  N X N  LOWER TRIANGULAR
C  ***  MATRIX STORED COMPACTLY BY ROWS.  X AND Y MAY OCCUPY THE SAME
C  ***  STORAGE.  ***
C
      INTEGER N
      DOUBLE PRECISION X(N), L(1), Y(N)
C     DIMENSION L(N*(N+1)/2)
      INTEGER I, II, IJ, I0, J, NP1
      DOUBLE PRECISION T, ZERO
C/6
C     DATA ZERO/0.D+0/
C/7
      PARAMETER (ZERO=0.D+0)
C/
C
      NP1 = N + 1
      I0 = N*(N+1)/2
      DO 20 II = 1, N
         I = NP1 - II
         I0 = I0 - I
         T = ZERO
         DO 10 J = 1, I
              IJ = I0 + J
              T = T + L(IJ)*Y(J)
 10           CONTINUE
         X(I) = T
 20      CONTINUE
 999  RETURN
C  ***  LAST CARD OF DL7VML FOLLOWS  ***
      END
