      SUBROUTINE I7COPY(P, Y, X)
C
C  ***  SET Y = X, WHERE X AND Y ARE INTEGER P-VECTORS  ***
C
      INTEGER P
      INTEGER X(P), Y(P)
C
      INTEGER I
C
      DO 10 I = 1, P
 10      Y(I) = X(I)
 999  RETURN
      END
