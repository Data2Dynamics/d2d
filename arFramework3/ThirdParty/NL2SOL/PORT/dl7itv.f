      SUBROUTINE DL7ITV(N, X, L, Y)
C
C  ***  SOLVE  (L**T)*X = Y,  WHERE  L  IS AN  N X N  LOWER TRIANGULAR
C  ***  MATRIX STORED COMPACTLY BY ROWS.  X AND Y MAY OCCUPY THE SAME
C  ***  STORAGE.  ***
C
      INTEGER N
      DOUBLE PRECISION X(N), L(1), Y(N)
      INTEGER I, II, IJ, IM1, I0, J, NP1
      DOUBLE PRECISION XI, ZERO
C/6
C     DATA ZERO/0.D+0/
C/7
      PARAMETER (ZERO=0.D+0)
C/
C
      DO 10 I = 1, N
 10      X(I) = Y(I)
      NP1 = N + 1
      I0 = N*(N+1)/2
      DO 30 II = 1, N
         I = NP1 - II
         XI = X(I)/L(I0)
         X(I) = XI
         IF (I .LE. 1) GO TO 999
         I0 = I0 - I
         IF (XI .EQ. ZERO) GO TO 30
         IM1 = I - 1
         DO 20 J = 1, IM1
              IJ = I0 + J
              X(J) = X(J) - XI*L(IJ)
 20           CONTINUE
 30      CONTINUE
 999  RETURN
C  ***  LAST CARD OF DL7ITV FOLLOWS  ***
      END
