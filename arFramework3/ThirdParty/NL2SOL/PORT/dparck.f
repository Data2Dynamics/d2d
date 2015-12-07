      SUBROUTINE DPARCK(ALG, D, IV, LIV, LV, N, V)
C
C  ***  CHECK ***SOL (VERSION 2.3) PARAMETERS, PRINT CHANGED VALUES  ***
C
C  ***  ALG = 1 FOR REGRESSION, ALG = 2 FOR GENERAL UNCONSTRAINED OPT.
C
      INTEGER ALG, LIV, LV, N
      INTEGER IV(LIV)
      DOUBLE PRECISION D(N), V(LV)
C
      DOUBLE PRECISION DR7MDC
      EXTERNAL DIVSET, DR7MDC,DV7CPY,DV7DFL
C DIVSET  -- SUPPLIES DEFAULT VALUES TO BOTH IV AND V.
C DR7MDC -- RETURNS MACHINE-DEPENDENT CONSTANTS.
C DV7CPY  -- COPIES ONE VECTOR TO ANOTHER.
C DV7DFL  -- SUPPLIES DEFAULT PARAMETER VALUES TO V ALONE.
C
C  ***  LOCAL VARIABLES  ***
C
      INTEGER ALG1, I, II, IV1, J, K, L, M, MIV1, MIV2, NDFALT, PARSV1,
     1        PU
      INTEGER IJMP, JLIM(4), MINIV(4), NDFLT(4)
C/6S
C     INTEGER VARNM(2), SH(2)
C     REAL CNGD(3), DFLT(3), VN(2,34), WHICH(3)
C/7S
      CHARACTER*1 VARNM(2), SH(2)
      CHARACTER*4 CNGD(3), DFLT(3), VN(2,34), WHICH(3)
C/
      DOUBLE PRECISION BIG, MACHEP, TINY, VK, VM(34), VX(34), ZERO
C
C  ***  IV AND V SUBSCRIPTS  ***
C
      INTEGER ALGSAV, DINIT, DTYPE, DTYPE0, EPSLON, INITS, IVNEED,
     1        LASTIV, LASTV, LMAT, NEXTIV, NEXTV, NVDFLT, OLDN,
     2        PARPRT, PARSAV, PERM, PRUNIT, VNEED
C
C
C/6
C     DATA ALGSAV/51/, DINIT/38/, DTYPE/16/, DTYPE0/54/, EPSLON/19/,
C    1     INITS/25/, IVNEED/3/, LASTIV/44/, LASTV/45/, LMAT/42/,
C    2     NEXTIV/46/, NEXTV/47/, NVDFLT/50/, OLDN/38/, PARPRT/20/,
C    3     PARSAV/49/, PERM/58/, PRUNIT/21/, VNEED/4/
C/7
      PARAMETER (ALGSAV=51, DINIT=38, DTYPE=16, DTYPE0=54, EPSLON=19,
     1           INITS=25, IVNEED=3, LASTIV=44, LASTV=45, LMAT=42,
     2           NEXTIV=46, NEXTV=47, NVDFLT=50, OLDN=38, PARPRT=20,
     3           PARSAV=49, PERM=58, PRUNIT=21, VNEED=4)
      SAVE BIG, MACHEP, TINY
C/
C
      DATA BIG/0.D+0/, MACHEP/-1.D+0/, TINY/1.D+0/, ZERO/0.D+0/
C/6S
C     DATA VN(1,1),VN(2,1)/4HEPSL,4HON../
C     DATA VN(1,2),VN(2,2)/4HPHMN,4HFC../
C     DATA VN(1,3),VN(2,3)/4HPHMX,4HFC../
C     DATA VN(1,4),VN(2,4)/4HDECF,4HAC../
C     DATA VN(1,5),VN(2,5)/4HINCF,4HAC../
C     DATA VN(1,6),VN(2,6)/4HRDFC,4HMN../
C     DATA VN(1,7),VN(2,7)/4HRDFC,4HMX../
C     DATA VN(1,8),VN(2,8)/4HTUNE,4HR1../
C     DATA VN(1,9),VN(2,9)/4HTUNE,4HR2../
C     DATA VN(1,10),VN(2,10)/4HTUNE,4HR3../
C     DATA VN(1,11),VN(2,11)/4HTUNE,4HR4../
C     DATA VN(1,12),VN(2,12)/4HTUNE,4HR5../
C     DATA VN(1,13),VN(2,13)/4HAFCT,4HOL../
C     DATA VN(1,14),VN(2,14)/4HRFCT,4HOL../
C     DATA VN(1,15),VN(2,15)/4HXCTO,4HL.../
C     DATA VN(1,16),VN(2,16)/4HXFTO,4HL.../
C     DATA VN(1,17),VN(2,17)/4HLMAX,4H0.../
C     DATA VN(1,18),VN(2,18)/4HLMAX,4HS.../
C     DATA VN(1,19),VN(2,19)/4HSCTO,4HL.../
C     DATA VN(1,20),VN(2,20)/4HDINI,4HT.../
C     DATA VN(1,21),VN(2,21)/4HDTIN,4HIT../
C     DATA VN(1,22),VN(2,22)/4HD0IN,4HIT../
C     DATA VN(1,23),VN(2,23)/4HDFAC,4H..../
C     DATA VN(1,24),VN(2,24)/4HDLTF,4HDC../
C     DATA VN(1,25),VN(2,25)/4HDLTF,4HDJ../
C     DATA VN(1,26),VN(2,26)/4HDELT,4HA0../
C     DATA VN(1,27),VN(2,27)/4HFUZZ,4H..../
C     DATA VN(1,28),VN(2,28)/4HRLIM,4HIT../
C     DATA VN(1,29),VN(2,29)/4HCOSM,4HIN../
C     DATA VN(1,30),VN(2,30)/4HHUBE,4HRC../
C     DATA VN(1,31),VN(2,31)/4HRSPT,4HOL../
C     DATA VN(1,32),VN(2,32)/4HSIGM,4HIN../
C     DATA VN(1,33),VN(2,33)/4HETA0,4H..../
C     DATA VN(1,34),VN(2,34)/4HBIAS,4H..../
C/7S
      DATA VN(1,1),VN(2,1)/'EPSL','ON..'/
      DATA VN(1,2),VN(2,2)/'PHMN','FC..'/
      DATA VN(1,3),VN(2,3)/'PHMX','FC..'/
      DATA VN(1,4),VN(2,4)/'DECF','AC..'/
      DATA VN(1,5),VN(2,5)/'INCF','AC..'/
      DATA VN(1,6),VN(2,6)/'RDFC','MN..'/
      DATA VN(1,7),VN(2,7)/'RDFC','MX..'/
      DATA VN(1,8),VN(2,8)/'TUNE','R1..'/
      DATA VN(1,9),VN(2,9)/'TUNE','R2..'/
      DATA VN(1,10),VN(2,10)/'TUNE','R3..'/
      DATA VN(1,11),VN(2,11)/'TUNE','R4..'/
      DATA VN(1,12),VN(2,12)/'TUNE','R5..'/
      DATA VN(1,13),VN(2,13)/'AFCT','OL..'/
      DATA VN(1,14),VN(2,14)/'RFCT','OL..'/
      DATA VN(1,15),VN(2,15)/'XCTO','L...'/
      DATA VN(1,16),VN(2,16)/'XFTO','L...'/
      DATA VN(1,17),VN(2,17)/'LMAX','0...'/
      DATA VN(1,18),VN(2,18)/'LMAX','S...'/
      DATA VN(1,19),VN(2,19)/'SCTO','L...'/
      DATA VN(1,20),VN(2,20)/'DINI','T...'/
      DATA VN(1,21),VN(2,21)/'DTIN','IT..'/
      DATA VN(1,22),VN(2,22)/'D0IN','IT..'/
      DATA VN(1,23),VN(2,23)/'DFAC','....'/
      DATA VN(1,24),VN(2,24)/'DLTF','DC..'/
      DATA VN(1,25),VN(2,25)/'DLTF','DJ..'/
      DATA VN(1,26),VN(2,26)/'DELT','A0..'/
      DATA VN(1,27),VN(2,27)/'FUZZ','....'/
      DATA VN(1,28),VN(2,28)/'RLIM','IT..'/
      DATA VN(1,29),VN(2,29)/'COSM','IN..'/
      DATA VN(1,30),VN(2,30)/'HUBE','RC..'/
      DATA VN(1,31),VN(2,31)/'RSPT','OL..'/
      DATA VN(1,32),VN(2,32)/'SIGM','IN..'/
      DATA VN(1,33),VN(2,33)/'ETA0','....'/
      DATA VN(1,34),VN(2,34)/'BIAS','....'/
C/
C
      DATA VM(1)/1.0D-3/, VM(2)/-0.99D+0/, VM(3)/1.0D-3/, VM(4)/1.0D-2/,
     1     VM(5)/1.2D+0/, VM(6)/1.D-2/, VM(7)/1.2D+0/, VM(8)/0.D+0/,
     2     VM(9)/0.D+0/, VM(10)/1.D-3/, VM(11)/-1.D+0/, VM(13)/0.D+0/,
     3     VM(15)/0.D+0/, VM(16)/0.D+0/, VM(19)/0.D+0/, VM(20)/-10.D+0/,
     4     VM(21)/0.D+0/, VM(22)/0.D+0/, VM(23)/0.D+0/, VM(27)/1.01D+0/,
     5     VM(28)/1.D+10/, VM(30)/0.D+0/, VM(31)/0.D+0/, VM(32)/0.D+0/,
     6     VM(34)/0.D+0/
      DATA VX(1)/0.9D+0/, VX(2)/-1.D-3/, VX(3)/1.D+1/, VX(4)/0.8D+0/,
     1     VX(5)/1.D+2/, VX(6)/0.8D+0/, VX(7)/1.D+2/, VX(8)/0.5D+0/,
     2     VX(9)/0.5D+0/, VX(10)/1.D+0/, VX(11)/1.D+0/, VX(14)/0.1D+0/,
     3     VX(15)/1.D+0/, VX(16)/1.D+0/, VX(19)/1.D+0/, VX(23)/1.D+0/,
     4     VX(24)/1.D+0/, VX(25)/1.D+0/, VX(26)/1.D+0/, VX(27)/1.D+10/,
     5     VX(29)/1.D+0/, VX(31)/1.D+0/, VX(32)/1.D+0/, VX(33)/1.D+0/,
     6     VX(34)/1.D+0/
C
C/6S
C     DATA VARNM(1)/1HP/, VARNM(2)/1HP/, SH(1)/1HS/, SH(2)/1HH/
C     DATA CNGD(1),CNGD(2),CNGD(3)/4H---C,4HHANG,4HED V/,
C    1     DFLT(1),DFLT(2),DFLT(3)/4HNOND,4HEFAU,4HLT V/
C/7S
      DATA VARNM(1)/'P'/, VARNM(2)/'P'/, SH(1)/'S'/, SH(2)/'H'/
      DATA CNGD(1),CNGD(2),CNGD(3)/'---C','HANG','ED V'/,
     1     DFLT(1),DFLT(2),DFLT(3)/'NOND','EFAU','LT V'/
C/
      DATA IJMP/33/, JLIM(1)/0/, JLIM(2)/24/, JLIM(3)/0/, JLIM(4)/24/,
     1     NDFLT(1)/32/, NDFLT(2)/25/, NDFLT(3)/32/, NDFLT(4)/25/
      DATA MINIV(1)/82/, MINIV(2)/59/, MINIV(3)/103/, MINIV(4)/103/
C
C...............................  BODY  ................................
C
      PU = 0
      IF (PRUNIT .LE. LIV) PU = IV(PRUNIT)
      IF (ALGSAV .GT. LIV) GO TO 20
      IF (ALG .EQ. IV(ALGSAV)) GO TO 20
         IF (PU .NE. 0) WRITE(PU,10) ALG, IV(ALGSAV)
 10      FORMAT(/40H THE FIRST PARAMETER TO DIVSET SHOULD BE,I3,
     1          12H RATHER THAN,I3)
         IV(1) = 67
         GO TO 999
 20   IF (ALG .LT. 1 .OR. ALG .GT. 4) GO TO 340
      MIV1 = MINIV(ALG)
      IF (IV(1) .EQ. 15) GO TO 360
      ALG1 = MOD(ALG-1,2) + 1
      IF (IV(1) .EQ. 0) CALL DIVSET(ALG, IV, LIV, LV, V)
      IV1 = IV(1)
      IF (IV1 .NE. 13 .AND. IV1 .NE. 12) GO TO 30
      IF (PERM .LE. LIV) MIV1 = MAX0(MIV1, IV(PERM) - 1)
      IF (IVNEED .LE. LIV) MIV2 = MIV1 + MAX0(IV(IVNEED), 0)
      IF (LASTIV .LE. LIV) IV(LASTIV) = MIV2
      IF (LIV .LT. MIV1) GO TO 300
      IV(IVNEED) = 0
      IV(LASTV) = MAX0(IV(VNEED), 0) + IV(LMAT) - 1
      IV(VNEED) = 0
      IF (LIV .LT. MIV2) GO TO 300
      IF (LV .LT. IV(LASTV)) GO TO 320
 30   IF (IV1 .LT. 12 .OR. IV1 .GT. 14) GO TO 60
         IF (N .GE. 1) GO TO 50
              IV(1) = 81
              IF (PU .EQ. 0) GO TO 999
              WRITE(PU,40) VARNM(ALG1), N
 40           FORMAT(/8H /// BAD,A1,2H =,I5)
              GO TO 999
 50      IF (IV1 .NE. 14) IV(NEXTIV) = IV(PERM)
         IF (IV1 .NE. 14) IV(NEXTV) = IV(LMAT)
         IF (IV1 .EQ. 13) GO TO 999
         K = IV(PARSAV) - EPSLON
         CALL DV7DFL(ALG1, LV-K, V(K+1))
         IV(DTYPE0) = 2 - ALG1
         IV(OLDN) = N
         WHICH(1) = DFLT(1)
         WHICH(2) = DFLT(2)
         WHICH(3) = DFLT(3)
         GO TO 110
 60   IF (N .EQ. IV(OLDN)) GO TO 80
         IV(1) = 17
         IF (PU .EQ. 0) GO TO 999
         WRITE(PU,70) VARNM(ALG1), IV(OLDN), N
 70      FORMAT(/5H /// ,1A1,14H CHANGED FROM ,I5,4H TO ,I5)
         GO TO 999
C
 80   IF (IV1 .LE. 11 .AND. IV1 .GE. 1) GO TO 100
         IV(1) = 80
         IF (PU .NE. 0) WRITE(PU,90) IV1
 90      FORMAT(/13H ///  IV(1) =,I5,28H SHOULD BE BETWEEN 0 AND 14.)
         GO TO 999
C
 100  WHICH(1) = CNGD(1)
      WHICH(2) = CNGD(2)
      WHICH(3) = CNGD(3)
C
 110  IF (IV1 .EQ. 14) IV1 = 12
      IF (BIG .GT. TINY) GO TO 120
         TINY = DR7MDC(1)
         MACHEP = DR7MDC(3)
         BIG = DR7MDC(6)
         VM(12) = MACHEP
         VX(12) = BIG
         VX(13) = BIG
         VM(14) = MACHEP
         VM(17) = TINY
         VX(17) = BIG
         VM(18) = TINY
         VX(18) = BIG
         VX(20) = BIG
         VX(21) = BIG
         VX(22) = BIG
         VM(24) = MACHEP
         VM(25) = MACHEP
         VM(26) = MACHEP
         VX(28) = DR7MDC(5)
         VM(29) = MACHEP
         VX(30) = BIG
         VM(33) = MACHEP
 120  M = 0
      I = 1
      J = JLIM(ALG1)
      K = EPSLON
      NDFALT = NDFLT(ALG1)
      DO 150 L = 1, NDFALT
         VK = V(K)
         IF (VK .GE. VM(I) .AND. VK .LE. VX(I)) GO TO 140
              M = K
              IF (PU .NE. 0) WRITE(PU,130) VN(1,I), VN(2,I), K, VK,
     1                                    VM(I), VX(I)
 130          FORMAT(/6H ///  ,2A4,5H.. V(,I2,3H) =,D11.3,7H SHOULD,
     1               11H BE BETWEEN,D11.3,4H AND,D11.3)
 140     K = K + 1
         I = I + 1
         IF (I .EQ. J) I = IJMP
 150     CONTINUE
C
      IF (IV(NVDFLT) .EQ. NDFALT) GO TO 170
         IV(1) = 51
         IF (PU .EQ. 0) GO TO 999
         WRITE(PU,160) IV(NVDFLT), NDFALT
 160     FORMAT(/13H IV(NVDFLT) =,I5,13H RATHER THAN ,I5)
         GO TO 999
 170  IF ((IV(DTYPE) .GT. 0 .OR. V(DINIT) .GT. ZERO) .AND. IV1 .EQ. 12)
     1                  GO TO 200
      DO 190 I = 1, N
         IF (D(I) .GT. ZERO) GO TO 190
              M = 18
              IF (PU .NE. 0) WRITE(PU,180) I, D(I)
 180     FORMAT(/8H ///  D(,I3,3H) =,D11.3,19H SHOULD BE POSITIVE)
 190     CONTINUE
 200  IF (M .EQ. 0) GO TO 210
         IV(1) = M
         GO TO 999
C
 210  IF (PU .EQ. 0 .OR. IV(PARPRT) .EQ. 0) GO TO 999
      IF (IV1 .NE. 12 .OR. IV(INITS) .EQ. ALG1-1) GO TO 230
         M = 1
         WRITE(PU,220) SH(ALG1), IV(INITS)
 220     FORMAT(/22H NONDEFAULT VALUES..../5H INIT,A1,14H..... IV(25) =,
     1          I3)
 230  IF (IV(DTYPE) .EQ. IV(DTYPE0)) GO TO 250
         IF (M .EQ. 0) WRITE(PU,260) WHICH
         M = 1
         WRITE(PU,240) IV(DTYPE)
 240     FORMAT(20H DTYPE..... IV(16) =,I3)
 250  I = 1
      J = JLIM(ALG1)
      K = EPSLON
      L = IV(PARSAV)
      NDFALT = NDFLT(ALG1)
      DO 290 II = 1, NDFALT
         IF (V(K) .EQ. V(L)) GO TO 280
              IF (M .EQ. 0) WRITE(PU,260) WHICH
 260          FORMAT(/1H ,3A4,9HALUES..../)
              M = 1
              WRITE(PU,270) VN(1,I), VN(2,I), K, V(K)
 270          FORMAT(1X,2A4,5H.. V(,I2,3H) =,D15.7)
 280     K = K + 1
         L = L + 1
         I = I + 1
         IF (I .EQ. J) I = IJMP
 290     CONTINUE
C
      IV(DTYPE0) = IV(DTYPE)
      PARSV1 = IV(PARSAV)
      CALL DV7CPY(IV(NVDFLT), V(PARSV1), V(EPSLON))
      GO TO 999
C
 300  IV(1) = 15
      IF (PU .EQ. 0) GO TO 999
      WRITE(PU,310) LIV, MIV2
 310  FORMAT(/10H /// LIV =,I5,17H MUST BE AT LEAST,I5)
      IF (LIV .LT. MIV1) GO TO 999
      IF (LV .LT. IV(LASTV)) GO TO 320
      GO TO 999
C
 320  IV(1) = 16
      IF (PU .NE. 0) WRITE(PU,330) LV, IV(LASTV)
 330  FORMAT(/9H /// LV =,I5,17H MUST BE AT LEAST,I5)
      GO TO 999
C
 340  IV(1) = 67
      IF (PU .NE. 0) WRITE(PU,350) ALG
 350  FORMAT(/10H /// ALG =,I5,21H MUST BE 1 2, 3, OR 4)
      GO TO 999
 360  IF (PU .NE. 0) WRITE(PU,370) LIV, MIV1
 370  FORMAT(/10H /// LIV =,I5,17H MUST BE AT LEAST,I5,
     1       37H TO COMPUTE TRUE MIN. LIV AND MIN. LV)
      IF (LASTIV .LE. LIV) IV(LASTIV) = MIV1
      IF (LASTV .LE. LIV) IV(LASTV) = 0
C
 999  RETURN
C  ***  LAST LINE OF DPARCK FOLLOWS  ***
      END
