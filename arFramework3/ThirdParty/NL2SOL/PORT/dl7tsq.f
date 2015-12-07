      SUBROUTINE DL7TSQ(N, A, L)
C
C  ***  SET A TO LOWER TRIANGLE OF (L**T) * L  ***
C
C  ***  L = N X N LOWER TRIANG. MATRIX STORED ROWWISE.  ***
C  ***  A IS ALSO STORED ROWWISE AND MAY SHARE STORAGE WITH L.  ***
C
      INTEGER N
      DOUBLE PRECISION A(1), L(1)
C     DIMENSION A(N*(N+1)/2), L(N*(N+1)/2)
C
      INTEGER I, II, IIM1, I1, J, K, M
      DOUBLE PRECISION LII, LJ
C
      II = 0
      DO 50 I = 1, N
         I1 = II + 1
         II = II + I
         M = 1
         IF (I .EQ. 1) GO TO 30
         IIM1 = II - 1
         DO 20 J = I1, IIM1
              LJ = L(J)
              DO 10 K = I1, J
                   A(M) = A(M) + LJ*L(K)
                   M = M + 1
 10                CONTINUE
 20           CONTINUE
 30      LII = L(II)
         DO 40 J = I1, II
 40           A(J) = LII * L(J)
 50      CONTINUE
C
 999  RETURN
C  ***  LAST CARD OF DL7TSQ FOLLOWS  ***
      END
