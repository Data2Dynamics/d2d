      SUBROUTINE DV7SCP(P, Y, S)
C
C  ***  SET P-VECTOR Y TO SCALAR S  ***
C
      INTEGER P
      DOUBLE PRECISION S, Y(P)
C
      INTEGER I
C
      DO 10 I = 1, P
 10      Y(I) = S
      RETURN
      END
