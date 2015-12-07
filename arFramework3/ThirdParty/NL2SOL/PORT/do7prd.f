      SUBROUTINE DO7PRD(L, LS, P, S, W, Y, Z)
C
C  ***  FOR I = 1..L, SET S = S + W(I)*Y(.,I)*(Z(.,I)**T), I.E.,
C  ***        ADD W(I) TIMES THE OUTER PRODUCT OF Y(.,I) AND Z(.,I).
C
      INTEGER L, LS, P
      DOUBLE PRECISION S(LS), W(L), Y(P,L), Z(P,L)
C     DIMENSION S(P*(P+1)/2)
C
      INTEGER I, J, K, M
      DOUBLE PRECISION WK, YI, ZERO
      DATA ZERO/0.D+0/
C
      DO 30 K = 1, L
         WK = W(K)
         IF (WK .EQ. ZERO) GO TO 30
         M = 1
         DO 20 I = 1, P
              YI = WK * Y(I,K)
              DO 10 J = 1, I
                   S(M) = S(M) + YI*Z(J,K)
                   M = M + 1
 10                CONTINUE
 20           CONTINUE
 30      CONTINUE
C
 999  RETURN
C  ***  LAST CARD OF DO7PRD FOLLOWS  ***
      END
