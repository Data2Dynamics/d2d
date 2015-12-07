      SUBROUTINE DS7IPR(P, IP, H)
C
C  APPLY THE PERMUTATION DEFINED BY IP TO THE ROWS AND COLUMNS OF THE
C  P X P SYMMETRIC MATRIX WHOSE LOWER TRIANGLE IS STORED COMPACTLY IN H.
C  THUS H.OUTPUT(I,J) = H.INPUT(IP(I), IP(J)).
C
      INTEGER P
      INTEGER IP(P)
      DOUBLE PRECISION H(1)
C
      INTEGER I, J, J1, JM, K, K1, KK, KM, KMJ, L, M
      DOUBLE PRECISION T
C
C ***  BODY  ***
C
      DO 90 I = 1, P
         J = IP(I)
         IF (J .EQ. I) GO TO 90
         IP(I) = IABS(J)
         IF (J .LT. 0) GO TO 90
         K = I
 10         J1 = J
            K1 = K
            IF (J .LE. K) GO TO 20
               J1 = K
               K1 = J
 20         KMJ = K1-J1
            L = J1-1
            JM = J1*L/2
            KM = K1*(K1-1)/2
            IF (L .LE. 0) GO TO 40
               DO 30 M = 1, L
                  JM = JM+1
                  T = H(JM)
                  KM = KM+1
                  H(JM) = H(KM)
                  H(KM) = T
 30               CONTINUE
 40         KM = KM+1
            KK = KM+KMJ
            JM = JM+1
            T = H(JM)
            H(JM) = H(KK)
            H(KK) = T
            J1 = L
            L = KMJ-1
            IF (L .LE. 0) GO TO 60
               DO 50 M = 1, L
                  JM = JM+J1+M
                  T = H(JM)
                  KM = KM+1
                  H(JM) = H(KM)
                  H(KM) = T
 50               CONTINUE
 60         IF (K1 .GE. P) GO TO 80
               L = P-K1
               K1 = K1-1
               KM = KK
               DO 70 M = 1, L
                  KM = KM+K1+M
                  JM = KM-KMJ
                  T = H(JM)
                  H(JM) = H(KM)
                  H(KM) = T
 70               CONTINUE
 80         K = J
            J = IP(K)
            IP(K) = -J
            IF (J .GT. I) GO TO 10
 90      CONTINUE
 999  RETURN
C  ***  LAST LINE OF DS7IPR FOLLOWS  ***
      END
