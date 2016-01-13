      SUBROUTINE DV2AXY(P, W, A, X, Y)
C
C  ***  SET W = A*X + Y  --  W, X, Y = P-VECTORS, A = SCALAR  ***
C
      INTEGER P
      DOUBLE PRECISION A, W(P), X(P), Y(P)
C
      INTEGER I
C
      DO 10 I = 1, P
 10      W(I) = A*X(I) + Y(I)
      RETURN
      END
