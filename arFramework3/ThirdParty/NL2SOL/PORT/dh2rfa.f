      SUBROUTINE DH2RFA(N, A, B, X, Y, Z)
C
C  ***  APPLY 2X2 HOUSEHOLDER REFLECTION DETERMINED BY X, Y, Z TO
C  ***  N-VECTORS A, B  ***
C
      INTEGER N
      DOUBLE PRECISION A(N), B(N), X, Y, Z
      INTEGER I
      DOUBLE PRECISION T
      DO 10 I = 1, N
         T = A(I)*X + B(I)*Y
         A(I) = A(I) + T
         B(I) = B(I) + T*Z
 10      CONTINUE
 999  RETURN
C  ***  LAST LINE OF DH2RFA FOLLOWS  ***
      END
