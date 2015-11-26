      SUBROUTINE DD7MLP(N, X, Y, Z, K)
C
C ***  SET X = DIAG(Y)**K * Z
C ***  FOR X, Z = LOWER TRIANG. MATRICES STORED COMPACTLY BY ROW
C ***  K = 1 OR -1.
C
      INTEGER N, K
C/6
C     DOUBLE PRECISION X(1), Y(N), Z(1)
C/7
      DOUBLE PRECISION X(*), Y(N), Z(*)
C/
      INTEGER I, J, L
      DOUBLE PRECISION ONE, T
      DATA ONE/1.D+0/
C
      L = 1
      IF (K .GE. 0) GO TO 30
      DO 20 I = 1, N
         T = ONE / Y(I)
         DO 10 J = 1, I
            X(L) = T * Z(L)
            L = L + 1
 10         CONTINUE
 20      CONTINUE
      GO TO 999
C
 30   DO 50 I = 1, N
         T = Y(I)
         DO 40 J = 1, I
            X(L) = T * Z(L)
            L = L + 1
 40         CONTINUE
 50      CONTINUE
 999  RETURN
C  ***  LAST CARD OF DD7MLP FOLLOWS  ***
      END
