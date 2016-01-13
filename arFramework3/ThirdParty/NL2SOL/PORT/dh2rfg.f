      DOUBLE PRECISION FUNCTION DH2RFG(A, B, X, Y, Z)
C
C  ***  DETERMINE X, Y, Z SO  I + (1,Z)**T * (X,Y)  IS A 2X2
C  ***  HOUSEHOLDER REFLECTION SENDING (A,B)**T INTO (C,0)**T,
C  ***  WHERE  C = -SIGN(A)*SQRT(A**2 + B**2)  IS THE VALUE DH2RFG
C  ***  RETURNS.
C
      DOUBLE PRECISION A, B, X, Y, Z
C
      DOUBLE PRECISION A1, B1, C, T
C/+
      DOUBLE PRECISION DSQRT
C/
      DOUBLE PRECISION ZERO
      DATA ZERO/0.D+0/
C
C  ***  BODY  ***
C
      IF (B .NE. ZERO) GO TO 10
         X = ZERO
         Y = ZERO
         Z = ZERO
         DH2RFG = A
         GO TO 999
 10   T = DABS(A) + DABS(B)
      A1 = A / T
      B1 = B / T
      C = DSQRT(A1**2 + B1**2)
      IF (A1 .GT. ZERO) C = -C
      A1 = A1 - C
      Z = B1 / A1
      X = A1 / C
      Y = B1 / C
      DH2RFG = T * C
 999  RETURN
C  ***  LAST LINE OF DH2RFG FOLLOWS  ***
      END
