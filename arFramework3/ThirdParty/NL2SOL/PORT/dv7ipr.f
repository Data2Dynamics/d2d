      SUBROUTINE DV7IPR(N, IP, X)
C
C     PERMUTE X SO THAT X.OUTPUT(I) = X.INPUT(IP(I)).
C     IP IS UNCHANGED ON OUTPUT.
C
      INTEGER N
      INTEGER IP(N)
      DOUBLE PRECISION X(N)
C
      INTEGER I, J, K
      DOUBLE PRECISION T
      DO 30 I = 1, N
         J = IP(I)
         IF (J .EQ. I) GO TO 30
         IF (J .GT. 0) GO TO 10
            IP(I) = -J
            GO TO 30
 10      T = X(I)
         K = I
 20      X(K) = X(J)
         K = J
         J = IP(K)
         IP(K) = -J
         IF (J .GT. I) GO TO 20
         X(K) = T
 30      CONTINUE
 999  RETURN
C  ***  LAST LINE OF DV7IPR FOLLOWS  ***
      END
