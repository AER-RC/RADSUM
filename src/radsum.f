C     path:      %P%
C     revision:  $Revision$
C     created:   $Date$  
C     presently: %H%  %T%
C------------------------------------------------------------------------------
      PROGRAM RADSUM 
C
C******************************************************************************
C
C     REVSD. FROM RADDR TO COMPUTE RADIANCE SUMS OVER MODEL LAYERS BY 
C
C                                       M. J. IACONO, S. A. CLOUGH
C                                       A.E.R. INC., 26 MARCH 1990
C
C     GENERALIZED, CORRECTED - R. D. WORSHAM   30 AUGUST 1991
C     
C     CORRECTED FOR THREE ANGLE SUMMATION  - M. J. IACONO  6 DECEMBER 1991  
C
C     ***REMINDER**************************************************************
C     ** THIS CODE MUST BE COMPILED THE SAME WAY AS LBLRTM, WHETHER SINGLE OR *
C     ** DOUBLE PRECISION IS USED, IN ORDER FOR THE INPUT TO BE READ PROPERLY.*
C     *************************************************************************
C
C******************************************************************************
C
C     READS PANELS OF LBLRTM RADIANCE OUTPUT AT INTERVALS OF DV WAVENUMBERS,
C     AND AVERAGES OVER GROUPS OF FACT*DV WAVENUMBERS FOR FLUX OUTPUT. 
C     FOR EXAMPLE: THERE ARE 21 PANELS (BOXES) FOR THE WAVENUMBER RANGE
C                  FROM 299.75 TO 310.25. (WITH DV=0.5)
C                  THE FLUXES WOULD BE OUTPUT FOR 300-305 AND 305-310 CM-1.
C
C     NANG ANGLES AND FIRST MOMENT QUADRATURE ARE USED FOR THE FLUX SUMMATION. 
C
C     INPUT:
C           I) TAPE31, TAPE32, ...
C              - LBLRTM OUTPUT FILES OF DOWNWELLING RADIANCE 
C          II) TAPE61, TAPE62, ...
C              - LBLRTM OUTPUT FILES OF UPWELLING RADIANCE
C         III) RADIN.DAT
C              CONTAINS:
C              - THE WAVENUMBER VALUE AT THE MIDDLE OF THE FIRST PANEL TO READ,
C              - THE WAVENUMBER OFFSET (IN DV) FOR THE DESIRED OUTPUT TO
C                BE OFFSET FROM THE VALUE V1P (IOFF)
C              - THE FACTOR (JDEL) BY WHICH TO MULTIPLY THE INCOMING DV TO
C                DETERMINE THE OUTPUT DV
C              - THE NUMBER OF ANGLES (NANG, CURRENTLY NANG<=3)
C              - THE NUMBER OF LEVELS (NLEV, CURRENTLY NLEV<=61)
C              - THE MODEL'S PRESSURE LEVELS (NLEV PRESSURES)
C
C     OUTPUT:
C           I) flxupdn.dat
C              CONTAINS (AT FACT*DV WAVENUMBER INTERVALS): 
C              - LBLRTM PRESSURE LEVELS
C              - THE UPWARD, DOWNWARD, AND NET (UP-DOWN) FLUXES FOR ALL
C                LBLRTM LEVELS
C              - THE HEATING RATE FOR EACH LAYER, COMPUTED FROM NET FLUX,
C                INDEXED AND WRITTEN WITH THE BOTTOM LEVEL OF THE LAYER
C******************************************************************************
C 
      PARAMETER (MXFSC=200,LIMMAX=2500,MXANGL=3)
C
      IMPLICIT DOUBLE PRECISION (V)
C
      DOUBLE PRECISION XID,SEC,HMOLID,XALTZ,YID
      INTEGER KFILD(MXANGL),KFILU(MXANGL),LFILE
      REAL HTR,NETFLX
      CHARACTER*8 HVRRAD
C
      DIMENSION RNUMW(LIMMAX)
      DIMENSION SRADD(LIMMAX,MXANGL),SRADU(LIMMAX,MXANGL)
      DIMENSION FLXTTD(MXFSC,LIMMAX),FLXTTU(MXFSC,LIMMAX),PRESLV(MXFSC)
      DIMENSION HTR(MXFSC,LIMMAX),NETFLX(MXFSC,LIMMAX),PRETHK(MXFSC)
      DIMENSION FILHDR(2),PNLHDR(2),IWDF(2),IWDP(2)
C
      COMMON /MANE/ RADD(LIMMAX,MXANGL),RADU(LIMMAX,MXANGL)
      COMMON /FILHDR/ XID(10),SEC,P0,T0,HMOLID(60),XALTZ(4),WK(60),
     *                PZL,PZU,TZL,TZU,WBROAD,DVT,V1V,V2V,TBOUND,
     *                EMISIV,FSCDID(17),NMOL,LAYRS,YID1,YID(10),LSTWDF
      COMMON /PNLHDR/ V1P,V2P,DVP,NLIM,LSTWDP 
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOG,RADCN1,RADCN2
C
      EQUIVALENCE (FILHDR(1),XID(1)), (PNLHDR(1),V1P),
     *            (XID(1),IWDF(1))  , (V1P,IWDP(1)),
     *            (FSCDID(1),IHIRAC), (FSCDID(2),ILBLF4),
     *            (FSCDID(3),ICNTNM), (FSCDID(4),IAERSL),
     *            (FSCDID(5),IEMIT) , (FSCDID(6),ISCAN),
     *            (FSCDID(7),IPLOT) , (FSCDID(8),IPATHL), 
     *            (FSCDID(9),JRAD)  , (FSCDID(10),ITEST),
     *            (FSCDID(11),IMRG) , (FSCDID(12),SCNID),
     *            (FSCDID(13),HWHM) , (FSCDID(14),IDABS),
     *            (FSCDID(15),IATM) , (FSCDID(16),LAYR1),
     *            (FSCDID(17),NLAYFS)
C 
C     Assign SCCS version number to module
C
      DATA HVRRAD / '$Revision$' /
C
C     This version of the code processes radiances computed for two
C     direction cosines:  0.84494897, and 0.35505103 for K=1.
C     The appropriate Gaussian weights (first moment quadrature)
C     for these angles are:
C
      DATA GWGD1,GWGD2/0.31804138,0.18195862/
C
C     For one angle (0.66666666667)
C
      DATA GWGO1/0.50/
C
C     For three angles  (0.91141204,0.59053314,0.21234054)
C
      DATA GWGT1,GWGT2,GWGT3/0.20093191,0.22924111,0.06982698/
C 
C     For 3 angle Lacis: secants = 1.219512195,2.43902439,3.658536585
C
      DATA WTLAC1,WTLAC2,WTLAC3/0.34914738,0.04243534,0.10841767/
C
C     HEATFC is the factor one must multiply DELTA-FLUX/DELTA-PRESSURE, 
C     with flux in W/M-2 and pressure in Mb, to get the heating rate in
C     units of Degrees/day.  It is equal to 
C
C           (g)x(#sec/day)x(1e-5)/(specific heat of air at const. p)
C        =  (9.8066)(3600)(1e-5)/(1.004)
C
      DATA HEATFC /8.4391/
C
C     Initialize variables
C
      DATA PRESLV / MXFSC*0.0 /
C
C******************************************************************************
C 
      NFHDRF = NWDL(IWDF,LSTWDF)
      NPHDRF = NWDL(IWDP,LSTWDP)
      PI = 2.*ASIN(1.)
      RADCN1 = 2.*PLANCK*CLIGHT*CLIGHT*1.E-07
      RADCN2 = PLANCK*CLIGHT/BOLTZ
C
C     Input wavenumber value of first panel to be read in,
c     number of pressure levels, number of angles, and 
c     output factor for DV
C  
      OPEN(UNIT=44,FILE='RADINIT')
      READ(44,442) IOFF,JDEL,NANG,NLEV,TBND,ILACIS
      IF (ILACIS .EQ. 1) NANG = 3
C
      IOPT = 0
 5    CALL OPNFIL(KFILD,KFILU,NANG,LFILE,IOPT)
      IF (IOPT.EQ.1) THEN
         READ(44,444,END=9999) IOFF,JDEL
      ELSE
         IOPT = 1
      ENDIF
C
C     Initialize array element increment for layer boundary pressures
C
      NLYR = 1
C
  10  CONTINUE
C
      DO 20 I = 1,NANG
         CALL BUFIN(KFILD(I),KEOF,FILHDR,NFHDRF)
         CALL BUFIN(KFILU(I),KEOF,FILHDR,NFHDRF)
  20  CONTINUE
      PRESLV(NLYR) = PZL
C 
      IF(KEOF .LT. 0) GO TO 10 
      IF(KEOF .EQ. 0) GO TO 1000 
C
  30  CONTINUE
C
      DO 40 I = 1,NANG
         CALL BUFIN(KFILD(I),KEOF,PNLHDR,NPHDRF)
         CALL BUFIN(KFILU(I),KEOF,PNLHDR,NPHDRF)
 40   CONTINUE
C      
      IF(KEOF .EQ. 0) GO TO 1000 
      IF(KEOF .LT. 0) GO TO 999
C
      DO 50 I = 1,NANG
         CALL BUFIN(KFILD(I),KEOF,RADD(1,I),NLIM) 
         CALL BUFIN(KFILU(I),KEOF,RADU(1,I),NLIM) 
  50  CONTINUE
C 
      KNUMW = (NLIM-1)/JDEL
      FACDV = FLOAT(JDEL)*DVP
C     
      IF (IOFF.LT.0) THEN
         JOFF = ABS(IOFF)
         KNUMW = 1
      ELSE
         JOFF = IOFF
      ENDIF
C     
      DO 60 K = 1,KNUMW+1
         RNUMW(K) = V1P+FLOAT(JOFF)*DVP+FACDV*FLOAT(K-1)
 60   CONTINUE
C
      L = 0
      DO 100 J = JOFF+1,NLIM,JDEL
         L = L+1
         IF (L.GT.KNUMW) GOTO 100
         DO 90 I = 1,NANG
C
            SRADD(L,I) = 0.5*RADD(J,I)
            SRADU(L,I) = 0.5*RADU(J,I)
C
            DO 80 K = J+1,J+JDEL-1
               SRADD(L,I) = SRADD(L,I)+RADD(K,I)
               SRADU(L,I) = SRADU(L,I)+RADU(K,I)
 80         CONTINUE
C
            SRADD(L,I) = SRADD(L,I)+0.5*RADD(J+JDEL,I)
            SRADU(L,I) = SRADU(L,I)+0.5*RADU(J+JDEL,I)
            SRADD(L,I) = SRADD(L,I)/FLOAT(JDEL)
            SRADU(L,I) = SRADU(L,I)/FLOAT(JDEL)
  90     CONTINUE
C
         IF (NANG.EQ.1) THEN
            FLXTTD(NLEV-NLYR,L) = GWGO1*SRADD(L,1)*FACDV*1.E04*2.*PI
            FLXTTU(NLYR+1,L)    = GWGO1*SRADU(L,1)*FACDV*1.E04*2.*PI
         ELSEIF (NANG.EQ.2) THEN
            FLXTTD(NLEV-NLYR,L) = (GWGD1*SRADD(L,1)+GWGD2*SRADD(L,2))
     *                            *FACDV*1.E04*2.*PI
            FLXTTU(NLYR+1,L)    = (GWGD1*SRADU(L,1)+GWGD2*SRADU(L,2))
     *                            *FACDV*1.E04*2.*PI
         ELSEIF (NANG.EQ.3 .AND. ILACIS .EQ. 0) THEN
            FLXTTD(NLEV-NLYR,L) = (GWGT1*SRADD(L,1)+GWGT2*SRADD(L,2)
     *                            +GWGT3*SRADD(L,3))*FACDV*1.E04*2.*PI
            FLXTTU(NLYR+1,L)    = (GWGT1*SRADU(L,1)+GWGT2*SRADU(L,2)
     *                            +GWGT3*SRADU(L,3))*FACDV*1.E04*2.*PI
         ELSEIF (NANG.EQ.3 .AND. ILACIS .EQ. 1) THEN
            FLXTTD(NLEV-NLYR,L) = (WTLAC1*SRADD(L,1)+WTLAC2*SRADD(L,2)
     *                            +WTLAC3*SRADD(L,3))*FACDV*1.E04*2.*PI
            FLXTTU(NLYR+1,L)    = (WTLAC1*SRADU(L,1)+WTLAC2*SRADU(L,2)
     *                            +WTLAC3*SRADU(L,3))*FACDV*1.E04*2.*PI
         ELSE
            STOP ' ERROR IN NANG '
         ENDIF
C
 100  CONTINUE
C
      GOTO 30
C 
999   CONTINUE
C
C     Return and do next level
C
      NLYR = NLYR+1
      GO TO 10
C
 1000 CONTINUE
C
C     Apply correction to upward sfc flux (by extrapolation)
C
C      DO 150 K = 1,KNUMW
C   150 FLXTTU(1,K) = 3*(FLXTTU(2,K)-FLXTTU(3,K))+FLXTTU(4,K)
C     
C     Calculate sfc flux from surface temperature at 0.5 wavenumber
C     intervals and sum over DVP intervals (here, 5 cm-1) 
C
      XKT = TBND/RADCN2
      DO 150 K = 1,KNUMW
         FSUM = 0.
         RVFRC = (RNUMW(K+1)-RNUMW(K))/10.0
         RVBAR1 = RNUMW(K)
         DO 160 KVF = 1,10
            RVBAR2 = RVBAR1+RVFRC
            RVBAR = (RVBAR1+RVBAR2)/2.0
            FSUM = FSUM+BBFCN(RVBAR,XKT)*FACDV*1.E04*PI
            RVBAR1 = RVBAR2
 160     CONTINUE
         FLXTTU(1,K) = FSUM/10.0
 150  CONTINUE
C
C  Compute net fluxes and heating rates, then output fluxes and heating
C  rates from top of atmosphere down for each level:
C
      DO 200 K = 1,KNUMW
         WRITE(LFILE,900)
         WRITE(LFILE,901) RNUMW(K)-0.25,RNUMW(K+1)-0.25,HVRRAD
         WRITE(LFILE,902)
         DO 200 N = NLEV,1,-1
            NETFLX(N,K) = FLXTTU(N,K)-FLXTTD(N,K)
            IF (N.EQ.NLEV) THEN
               HTR(N,K) = 0.
               PRESLV(N) = PZU
            ELSE
               HTR(N,K) = NETFLX(N,K)-NETFLX(N+1,K)
               PRETHK(N) = PRESLV(N)-PRESLV(N+1)
               HTR(N,K) = HEATFC*HTR(N,K)/PRETHK(N)
            ENDIF
            WRITE(LFILE,903) N-1,PRESLV(N),FLXTTU(N,K),
     *           FLXTTD(N,K),NETFLX(N,K),HTR(N,K)
 200     CONTINUE
C     
      GOTO 5
C
C     Formats:
C     
 442  FORMAT(4I5,F8.1,I5)
 444  FORMAT(2I5)
 900  FORMAT('1')
 901  FORMAT('WAVENUMBER BAND: ',F8.2,'-',F8.2,' CM -1',15X,
     *       'RADSUM SCCS version ',A8)
 902  FORMAT(' LEV   PRESSURE       FLUX UP       FLUX DOWN',
     *       '      NET FLUX     HEATING RATE',/,
     *       '         MB           W/M2            W/M2  ',
     *       '       W/M2         DEG/DAY   ')
 903  FORMAT(1X,0P,I2,1P,4(2X,E13.6),2X,E18.11) 
C
 9999 STOP
C
      END 
C
C------------------------------------------------------------------------------
C
      SUBROUTINE OPNFIL(KFILD,KFILU,NANG,LFILE,IOPT)
C 
C******************************************************************************
C     THIS SUBROUTINE OPENS THE NEEDED FILES. 
C 
      INTEGER KFILD(NANG), KFILU(NANG), LFILE
      CHARACTER TAPE*4,KFIL*6,CFORM*11
C
      DATA CFORM/'UNFORMATTED'/
C
      DATA TAPE/'TAPE'/
C 
      DO 10 I = 1,NANG
         IF (IOPT.EQ.0) THEN
            KFILD(I) = 30+I
            KFILU(I) = 60+I
            WRITE(KFIL,'(A4,I2.2)') TAPE,KFILD(I)
            OPEN (KFILD(I),FILE=KFIL,FORM=CFORM)
            WRITE(KFIL,'(A4,I2.2)') TAPE,KFILU(I)
            OPEN (KFILU(I),FILE=KFIL,FORM=CFORM)
         ENDIF  
         REWIND KFILD(I)
         REWIND KFILU(I)
  10  CONTINUE
C     
      IF (IOPT.EQ.0) THEN
         LFILE = 21
         OPEN (LFILE,FILE='flxupdn.dat')
         REWIND LFILE
      ENDIF
C
      RETURN
      END 
C
C------------------------------------------------------------------------------
C
      SUBROUTINE RHEAD(XID,SV,DSV,NTOT,NVAR)
C 
C******************************************************************************
C     This subroutine extracts information from XID.
C 
      REAL XID(10), SV, DSV 
      INTEGER NTOT, NVAR
      CHARACTER*80 HEADER 
C 
      PRINT 10,XID
      WRITE(HEADER,20) XID
      READ(HEADER,30) SV,DSV,NTOT,NVAR
C
      RETURN
C
C     Formats:
C     
 10   FORMAT(' XID',15A8) 
 20   FORMAT(10A8)
 30   FORMAT(8X,1PE12.4,17X,1PE12.4,8X,I5,5X,I5,8X) 
C
      END 
C
C------------------------------------------------------------------------------
C
      FUNCTION NWDL(IWD,ILAST)                                           
C******************************************************************************
      DIMENSION IWD(*)                                                   
      ILAST = -654321                                                      
      DO 10 I = 1,9000
         IF(IWD(I).NE.ILAST) GO TO 10                                       
         NWDL = I-1
         GO TO 12
 10   CONTINUE
 12   RETURN
      END
C
C------------------------------------------------------------------------------
C
      SUBROUTINE ENDFIL(IFILE)
C******************************************************************************
      DIMENSION IDUM(6)
      DATA IDUM /6*-99/
      CALL BUFOUT(IFILE,IDUM,6)
      RETURN
      END
C
C------------------------------------------------------------------------------
C
      FUNCTION BBFCN(XVI,XKT)
C
C     Function BBFCN calculates the black body function for
C     wavenumber value XVI
C
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOG,RADCN1,RADCN2
C
      IF (XKT.GT.0.0) THEN
C
         XVIOKT = XVI/XKT
C
         IF (XVIOKT.LE.0.01) THEN
            BBFCN = RADCN1*(XVI**2) *XKT/(1.+0.5*XVIOKT)
         ELSEIF (XVIOKT.LE.80.0) THEN
            BBFCN = RADCN1*(XVI**3)/(EXP(XVIOKT)-1.)
         ELSE
            BBFCN = 0.
         ENDIF
      ELSE
         BBFCN = 0.
      ENDIF
C
      RETURN
      END
C     
C
C------------------------------------------------------------------------------
C
      BLOCK DATA
C
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOG,RADCN1,RADCN2
      DATA PLANCK/6.626176E-27/, BOLTZ/1.380662E-16/,
     *     CLIGHT/2.99792458E10/, AVOG/6.022045E23/
C
      END
C
C------------------------------------------------------------------------------
C




