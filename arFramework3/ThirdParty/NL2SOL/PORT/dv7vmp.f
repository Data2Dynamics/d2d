      SUBROUTINE DV7VMP(N, X, Y, Z, K)
C
C ***  SET X(I) = Y(I) * Z(I)**K, 1 .LE. I .LE. N (FOR K = 1 OR -1)  ***
C
      INTEGER N, K
      DOUBLE PRECISION X(N), Y(N), Z(N)
      INTEGER I
C
      IF (K .GE. 0) GO TO 20
      DO 10 I = 1, N
 10      X(I) = Y(I) / Z(I)
      GO TO 999
C
 20   DO 30 I = 1, N
 30      X(I) = Y(I) * Z(I)
 999  RETURN
C  ***  LAST CARD OF DV7VMP FOLLOWS  ***
      END
