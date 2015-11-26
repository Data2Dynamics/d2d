      SUBROUTINE DL7TVM(N, X, L, Y)
C
C  ***  COMPUTE  X = (L**T)*Y, WHERE  L  IS AN  N X N  LOWER
C  ***  TRIANGULAR MATRIX STORED COMPACTLY BY ROWS.  X AND Y MAY
C  ***  OCCUPY THE SAME STORAGE.  ***
C
      INTEGER N
      DOUBLE PRECISION X(N), L(1), Y(N)
C     DIMENSION L(N*(N+1)/2)
      INTEGER I, IJ, I0, J
      DOUBLE PRECISION YI, ZERO
C/6
C     DATA ZERO/0.D+0/
C/7
      PARAMETER (ZERO=0.D+0)
C/
C
      I0 = 0
      DO 20 I = 1, N
         YI = Y(I)
         X(I) = ZERO
         DO 10 J = 1, I
              IJ = I0 + J
              X(J) = X(J) + YI*L(IJ)
 10           CONTINUE
         I0 = I0 + I
 20      CONTINUE
 999  RETURN
C  ***  LAST CARD OF DL7TVM FOLLOWS  ***
      END
