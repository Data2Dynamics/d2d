      SUBROUTINE I7SHFT(N, K, X)
C
C  ***  SHIFT X(K),...,X(N) LEFT CIRCULARLY ONE POSITION IF K .GT. 0.
C  ***  SHIFT X(-K),...,X(N) RIGHT CIRCULARLY ONE POSITION IF K .LT. 0.
C
      INTEGER N, K
      INTEGER X(N)
C
      INTEGER I, II, K1, NM1, T
C
      IF (K .LT. 0) GO TO 20
      IF (K .GE. N) GO TO 999
      NM1 = N - 1
      T = X(K)
      DO 10 I = K, NM1
 10      X(I) = X(I+1)
      X(N) = T
      GO TO 999
C
 20   K1 = -K
      IF (K1 .GE. N) GO TO 999
      T = X(N)
      NM1 = N - K1
      DO 30 II = 1, NM1
         I = N - II
         X(I+1) = X(I)
 30      CONTINUE
      X(K1) = T
 999  RETURN
C  ***  LAST LINE OF I7SHFT FOLLOWS  ***
      END
