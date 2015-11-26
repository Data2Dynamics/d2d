      SUBROUTINE DS7DMP(N, X, Y, Z, K)
C
C ***  SET X = DIAG(Z)**K * Y * DIAG(Z)**K
C ***  FOR X, Y = COMPACTLY STORED LOWER TRIANG. MATRICES
C ***  K = 1 OR -1.
C
      INTEGER N, K
C/6S
C     DOUBLE PRECISION X(1), Y(1), Z(N)
C/7S
      DOUBLE PRECISION X(*), Y(*), Z(N)
C/
      INTEGER I, J, L
      DOUBLE PRECISION ONE, T
      DATA ONE/1.D+0/
C
      L = 1
      IF (K .GE. 0) GO TO 30
      DO 20 I = 1, N
         T = ONE / Z(I)
         DO 10 J = 1, I
            X(L) = T * Y(L) / Z(J)
            L = L + 1
 10         CONTINUE
 20      CONTINUE
      GO TO 999
C
 30   DO 50 I = 1, N
         T = Z(I)
         DO 40 J = 1, I
            X(L) = T * Y(L) * Z(J)
            L = L + 1
 40         CONTINUE
 50      CONTINUE
 999  RETURN
C  ***  LAST CARD OF DS7DMP FOLLOWS  ***
      END
