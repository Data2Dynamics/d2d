      SUBROUTINE I7PNVR(N, X, Y)
C
C  ***  SET PERMUTATION VECTOR X TO INVERSE OF Y  ***
C
      INTEGER N
      INTEGER X(N), Y(N)
C
      INTEGER I, J
      DO 10 I = 1, N
         J = Y(I)
         X(J) = I
 10      CONTINUE
C
 999  RETURN
C  ***  LAST LINE OF I7PNVR FOLLOWS  ***
      END
