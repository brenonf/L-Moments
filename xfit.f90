PROGRAM XFIT
!C***********************************************************************
!C*                                                                     *
!C*  FORTRAN CODE WRITTEN FOR INCLUSION IN IBM RESEARCH REPORT RC20525, *
!C*  'FORTRAN ROUTINES FOR USE WITH THE METHOD OF L-MOMENTS, VERSION 3' *
!C*                                                                     *
!C*  J. R. M. HOSKING                                                   *
!C*  IBM RESEARCH DIVISION                                              *
!C*  T. J. WATSON RESEARCH CENTER                                       *
!C*  YORKTOWN HEIGHTS                                                   *
!C*  NEW YORK 10598, U.S.A.                                             *
!C*                                                                     *
!C*  VERSION 3     AUGUST 1996                                          *
!C*                                                                     *
!C***********************************************************************
!C
!C  EXAMPLE PROGRAM FOR REGIONAL FREQUENCY ANALYSIS USING THE METHOD OF
!C  L-MOMENTS. THE PROGRAM FITS A DISTRIBUTION TO REGIONAL DATA AND USES
!C  IT TO ESTIMATE QUANTILES AT EACH SITE.
!C
!C  THIS EXAMPLE FITS A WAKEBY DISTRIBUTION, USING A VARIANT (PLOTTING
!C  POSITION ESTIMATORS INSTEAD OF UNBIASED) OF THE REGIONAL L-MOMENT
!C  ALGORITHM DESCRIBED BY HOSKING AND WALLIS ("REGIONAL FREQUENCY
!C  ANALYSIS: AN APPROACH BASED ON L-MOMENTS", CAMBRIDGE UNIV. PRESS,
!C  1997).  TO FIT A DIFFERENT DISTRIBUTION, REPLACE THE CALLS TO
!C  SUBROUTINES PELWAK AND QUAWAK BY THE APPROPRIATE PEL... AND QUA...
!C  ROUTINES, CHANGE THE 'WAKEBY' IN FORMAT STATEMENT 6030 AND CHANGE
!C  THE VALUE OF PARAMETER NPAR.
!C
!C  PARAMETERS OF PROGRAM:
!C  MAXNS  - SHOULD BE AT LEAST AS LARGE AS THE NUMBER OF SITES IN THE
!C           REGION
!C  MAXN   - SHOULD BE AT LEAST AS LARGE AS THE LARGEST RECORD LENGTH
!C           AT ANY SITE IN THE REGION
!C  NPAR   - NUMBER OF PARAMETERS IN THE DISTRIBUTION TO BE FITTED
!C           (5 FOR WAKEBY, OF COURSE)
!C  NPROB  - NUMBER OF FLOOD QUANTILES TO BE ESTIMATED AT EACH SITE
!C  INFILE - STREAM NUMBER TO WHICH INPUT FILE IS ATTACHED

!C  ARRAYS TO BE INITIALIZED IN DATA STATEMENTS:
!C  PROB(NPROB) - PROBABILITIES FOR WHICH QUANTILES ARE TO BE ESTIMATED
!C
!C  VARIABLES TO BE INITIALIZED IN DATA STATEMENTS:
!C  A      - ) PARAMETERS OF
!C  B      - ) PLOTTING POSITION

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXNS=100,MAXN=200,NPAR=5,NPROB=10,INFILE=7)
      CHARACTER*32 SITEID
      DOUBLE PRECISION PROB(NPROB),QUANT(NPROB),PARA(5),RMOM(5),RQUANT(NPROB),WEIGHT(MAXNS),X(MAXN),XMOM(5,MAXNS)

      DATA PROB/0.1D0,0.2D0,0.5D0,0.8D0,0.9D0,0.95D0,0.98D0,0.99D0,0.999D0,0.9999D0/
      DATA A,B/-0.35D0,0D0/

      IF(A.EQ.0D0.AND.B.EQ.0D0)WRITE(6,6000)
      IF(A.NE.0D0.OR.B.NE.0D0)WRITE(6,6010)A,B

!         READ THE DATA AND CALCULATE AT-SITE L-MOMENTS.
!         ASSUMED STRUCTURE OF DATA FILE IS AS FOLLOWS.
!         1. ONE RECORD CONTAINING THE NUMBER OF SITES IN THE REGION.
!         2. FOR EACH SITE:
!            A  ONE RECORD CONTAINING AN IDENTIFYING LABEL FOR THE SITE;
!            B. ONE RECORD CONTAINING THE RECORD LENGTH AT THE SITE;
!            C. THE DATA VALUES, IN FREE FORMAT.
      OPEN(FILE="MAXWIND.DAT",UNIT=15,status='unknown')
      READ(15,*)NSITE
      DO 10 ISITE=1,NSITE
      READ(15,'(A32)')SITEID
      READ(15,*)N
      READ(15,*)(X(I),I=1,N)
      WEIGHT(ISITE)=N
      CALL SORT(X,N)
      CALL SAMLMR(X,N,XMOM(1,ISITE),NPAR,A,B)
      WRITE(6,6020)ISITE,SITEID,N,(XMOM(I,ISITE),I=1,NPAR)
   10 CONTINUE

      CALL REGLMR(NSITE,NPAR,5,XMOM,WEIGHT,RMOM)
      WRITE(6,6030)(RMOM(I),I=1,NPAR)

!         FIT REGIONAL FREQUENCY DISTRIBUTION

      CALL PELWAK(RMOM,PARA,IFAIL)
      IF(IFAIL.NE.0)WRITE(6,6040)IFAIL
      WRITE(6,6050)(PARA(I),I=1,NPAR)

!         CALCULATE QUANTILES OF REGIONAL FREQUENCY DISTRIBUTION

      WRITE(6,6060)(PROB(IQ),IQ=1,NPROB)
      DO 20 IQ=1,NPROB
   20 RQUANT(IQ)=QUAWAK(PROB(IQ),PARA)
      WRITE(6,6070)(RQUANT(IQ),IQ=1,NPROB)

!         CALCULATE QUANTILE ESTIMATES FOR EACH SITE

      DO 40 ISITE=1,NSITE
      DO 30 IQ=1,NPROB
   30 QUANT(IQ)=XMOM(1,ISITE)*RQUANT(IQ)
      WRITE(6,6080)ISITE,(QUANT(IQ),IQ=1,NPROB)
   40 CONTINUE

      STOP

 6000 FORMAT(' REGIONAL ANALYSIS, UNBIASED L-MOMENTS'/)
 6010 FORMAT(' REGIONAL ANALYSIS,',   ' L-MOMENT PLOTTING POSITION PARAMETERS ',2F8.4/)
 6020 FORMAT(' SITE',I3,1X,A32,'N=',I3,'   L-MOMENT RATIOS', F9.2,4F9.4)
 6030 FORMAT(//' REGIONAL AVERAGE L-MOMENT RATIOS',5F9.4)
 6040 FORMAT(/' PARAMETER ESTIMATION: FAIL FLAG',I2)
 6050 FORMAT(/' REGIONAL WAKEBY PARAMETERS',5F12.4)
 6060 FORMAT(///'  SITE',25X,'QUANTILES'/' NUMBER',10F10.4/1X,106('-'))
 6070 FORMAT(' REGION',10F10.2/)
 6080 FORMAT(1X,I4,2X,10F10.2)
      END PROGRAM XFIT




!***************************************************************************************
!-------------------------------  Subrotinas  -----------------------------------------
!**************************************************************************************
SUBROUTINE SORT(X,N)
!C***********************************************************************
!C*                                                                     *
!C*  FORTRAN CODE WRITTEN FOR INCLUSION IN IBM RESEARCH REPORT RC20525, *
!C*  'FORTRAN ROUTINES FOR USE WITH THE METHOD OF L-MOMENTS, VERSION 3' *
!C*                                                                     *
!C*  J. R. M. HOSKING                                                   *
!C*  IBM RESEARCH DIVISION                                              *
!C*  T. J. WATSON RESEARCH CENTER                                       *
!C*  YORKTOWN HEIGHTS                                                   *
!C*  NEW YORK 10598, U.S.A.                                             *
!C*                                                                     *
!C*  VERSION 3     AUGUST 1996                                          *
!C*                                                                     *
!C***********************************************************************
!C
!C  SORTS THE ARRAY X INTO ASCENDING ORDER
!
!C  PARAMETERS OF ROUTINE:
!C  X      *IN/OUT* ARRAY OF LENGTH N. CONTAINS THE NUMBERS TO BE SORTED.
!C                  ON EXIT, CONTAINS THE SORTED NUMBERS.
!C  N      * INPUT* NUMBER OF ELEMENTS TO BE SORTED
!C
!C  METHOD USED IS SHELL SORT WITH SEQUENCE OF INCREMENTS AS IN
!C  D.F.KNUTH (1969) 'THE ART OF COMPUTER PROGRAMMING', VOL.3, P.95
!C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X(N)
      IF(N.LE.1)RETURN
      J=4
      DO 10 I=1,100
      J=3*J+1
      IF(J.GE.N)GOTO 20
   10 CONTINUE
   20 CONTINUE
      M=(J/3)
      DO 60 MM=1,100
      M=M/3
      IF(M.EQ.0)RETURN
      DO 50 I=M+1,N
      TEST=X(I)
      J=I
      DO 30 JJ=1,100
      J=J-M
      IF(J.LE.0)GOTO 40
      IF(TEST.GE.X(J))GOTO 40
   30 X(J+M)=X(J)
   40 CONTINUE
   50 X(J+M)=TEST
   60 CONTINUE
      END


!*************************************************************************************
!_------------------------------------------------------------------------------------
!*************************************************************************************

SUBROUTINE SAMLMR(X,N,XMOM,NMOM,A,B)
!C***********************************************************************
!C*                                                                     *
!C*  FORTRAN CODE WRITTEN FOR INCLUSION IN IBM RESEARCH REPORT RC20525, *
!C*  'FORTRAN ROUTINES FOR USE WITH THE METHOD OF L-MOMENTS, VERSION 3' *
!C*                                                                     *
!C*  J. R. M. HOSKING                                                   *
!C*  IBM RESEARCH DIVISION                                              *
!C*  T. J. WATSON RESEARCH CENTER                                       *
!C*  YORKTOWN HEIGHTS                                                   *
!C*  NEW YORK 10598, U.S.A.                                             *
!C*                                                                     *
!C*  VERSION 3     AUGUST 1996                                          *
!C*                                                                     *
!C***********************************************************************
!C
!C  SAMPLE L-MOMENTS OF A DATA ARRAY
!C
!C  PARAMETERS OF ROUTINE:
!C  X      * INPUT* ARRAY OF LENGTH N. CONTAINS THE DATA, IN ASCENDING
!C                  ORDER.
!C  N      * INPUT* NUMBER OF DATA VALUES
!C  XMOM   *OUTPUT* ARRAY OF LENGTH NMOM. ON EXIT, CONTAINS THE SAMPLE
!C                  L-MOMENTS L-1, L-2, T-3, T-4, ... .
!C  NMOM   * INPUT* NUMBER OF L-MOMENTS TO BE FOUND. AT MOST MAX(N,20).
!C  A      * INPUT* ) PARAMETERS OF PLOTTING
!C  B      * INPUT* ) POSITION (SEE BELOW)
!C
!C  FOR UNBIASED ESTIMATES (OF THE LAMBDA'S) SET A=B=ZERO. OTHERWISE,
!C  PLOTTING-POSITION ESTIMATORS ARE USED, BASED ON THE PLOTTING POSITION
!C  (J+A)/(N+B)  FOR THE J'TH SMALLEST OF N OBSERVATIONS. FOR EXAMPLE,
!C  A=-0.35D0 AND B=0.0D0 YIELDS THE ESTIMATORS RECOMMENDED BY
!C  HOSKING ET AL. (1985, TECHNOMETRICS) FOR THE GEV DISTRIBUTION.
!C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X(N),XMOM(NMOM),SUM(20)
      DATA ZERO/0D0/,ONE/1D0/
      IF(NMOM.GT.20.OR.NMOM.GT.N)GOTO 1000
      DO 10 J=1,NMOM
   10 SUM(J)=ZERO
      IF(A.EQ.ZERO.AND.B.EQ.ZERO)GOTO 50
      IF(A.LE.-ONE.OR.A.GE.B)GOTO 1010

!C         PLOTTING-POSITION ESTIMATES OF PWM'S

      DO 30 I=1,N
      PPOS=(I+A)/(N+B)
      TERM=X(I)
      SUM(1)=SUM(1)+TERM
      DO 20 J=2,NMOM
      TERM=TERM*PPOS
   20 SUM(J)=SUM(J)+TERM
   30 CONTINUE
      DO 40 J=1,NMOM
   40 SUM(J)=SUM(J)/N
      GOTO 100

!C         UNBIASED ESTIMATES OF PWM'S

   50 DO 70 I=1,N
      Z=I
      TERM=X(I)
      SUM(1)=SUM(1)+TERM
      DO 60 J=2,NMOM
      Z=Z-ONE
      TERM=TERM*Z
   60 SUM(J)=SUM(J)+TERM
   70 CONTINUE
      Y=N
      Z=N
      SUM(1)=SUM(1)/Z
      DO 80 J=2,NMOM
      Y=Y-ONE
      Z=Z*Y
   80 SUM(J)=SUM(J)/Z

!C         L-MOMENTS

  100 K=NMOM
      P0=ONE
      IF(NMOM-NMOM/2*2.EQ.1)P0=-ONE
      DO 120 KK=2,NMOM
      AK=K
      P0=-P0
      P=P0
      TEMP=P*SUM(1)
      DO 110 I=1,K-1
      AI=I
      P=-P*(AK+AI-ONE)*(AK-AI)/(AI*AI)
  110 TEMP=TEMP+P*SUM(I+1)
      SUM(K)=TEMP
  120 K=K-1
      XMOM(1)=SUM(1)
      IF(NMOM.EQ.1)RETURN
      XMOM(2)=SUM(2)
      IF(SUM(2).EQ.ZERO)GOTO 1020
      IF(NMOM.EQ.2)RETURN
      DO 130 K=3,NMOM
  130 XMOM(K)=SUM(K)/SUM(2)
      RETURN

 1000 WRITE(6,7000)
      RETURN
 1010 WRITE(6,7010)
      RETURN
 1020 WRITE(6,7020)
      RETURN

 7000 FORMAT(' *** ERROR *** ROUTINE SAMLMR : PARAMETER NMOM INVALID')
 7010 FORMAT(' *** ERROR *** ROUTINE SAMLMR :', ' PLOTTING-POSITION PARAMETERS INVALID')
 7020 FORMAT(' *** ERROR *** ROUTINE SAMLMR : ALL DATA VALUES EQUAL')
      END

!*************************************************************************************
!_------------------------------------------------------------------------------------
!*************************************************************************************


SUBROUTINE REGLMR(NSITE,NMOM,NXMOM,XMOM,WEIGHT,RMOM)
!C***********************************************************************
!C*                                                                     *
!C*  FORTRAN CODE WRITTEN FOR INCLUSION IN IBM RESEARCH REPORT RC20525, *
!C*  'FORTRAN ROUTINES FOR USE WITH THE METHOD OF L-MOMENTS, VERSION 3' *
!C*                                                                     *
!C*  J. R. M. HOSKING                                                   *
!C*  IBM RESEARCH DIVISION                                              *
!C*  T. J. WATSON RESEARCH CENTER                                       *
!C*  YORKTOWN HEIGHTS                                                   *
!C*  NEW YORK 10598, U.S.A.                                             *
!C*                                                                     *
!C*  VERSION 3     AUGUST 1996                                          *
!C*                                                                     *
!C***********************************************************************
!C
!C  REGIONAL WEIGHTED AVERAGE OF L-MOMENTS
!C
!C  PARAMETERS OF ROUTINE:
!C  NSITE  * INPUT* NUMBER OF SITES IN REGION
!C  NMOM   * INPUT* NUMBER OF L-MOMENTS TO BE FOUND.
!C  NXMOM  * INPUT* THE FIRST DIMENSION OF ARRAY XMOM, AS DECLARED IN THE
!C                  CALLING PROGRAM.
!C  XMOM   * INPUT* ARRAY OF DIMENSION (NXMOM,NSITE). X(I,J) CONTAINS
!C                  THE I'TH L-MOMENT RATIO FOR SITE J.
!C  WEIGHT * INPUT* ARRAY OF LENGTH NSITE. CONTAINS THE WEIGHTS TO BE
!C                  APPLIED TO EACH SITE.
!C  RMOM   *OUTPUT* ARRAY OF LENGTH NMOM. ON EXIT, CONTAINS THE REGIONAL
!C                  WEIGHTED AVERAGE L-MOMENT RATIOS.
!C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION XMOM(NXMOM,NSITE),WEIGHT(NSITE),RMOM(NMOM)
      DATA ZERO/0D0/,ONE/1D0/
      IF(NMOM.LT.2.OR.NMOM.GT.NXMOM)GOTO 1000
      DO 10 J=1,NMOM
   10 RMOM(J)=ZERO
      WSUM=ZERO
      DO 30 ISITE=1,NSITE
      SMEAN=XMOM(1,ISITE)
      IF(SMEAN.EQ.ZERO)GOTO 1010
      W=WEIGHT(ISITE)
      WSUM=WSUM+W
      RMOM(2)=RMOM(2)+W*XMOM(2,ISITE)/SMEAN
      IF(NMOM.EQ.2)GOTO 30
      DO 20 J=3,NMOM
   20 RMOM(J)=RMOM(J)+W*XMOM(J,ISITE)
   30 CONTINUE
      IF(WSUM.LE.ZERO)GOTO 1020
      RMOM(1)=ONE
      RMOM(2)=RMOM(2)/WSUM
      IF(NMOM.EQ.2)RETURN
      DO 40 J=3,NMOM
   40 RMOM(J)=RMOM(J)/WSUM
      RETURN

 1000 WRITE(6,7000)
      RETURN
 1010 WRITE(6,7010)ISITE
      RETURN
 1020 WRITE(6,7020)
      RETURN

 7000 FORMAT(' *** ERROR *** ROUTINE REGLMR : PARAMETER NMOM INVALID')
 7010 FORMAT(' *** ERROR *** ROUTINE REGLMR : ZERO MEAN AT SITE',I4)
 7020 FORMAT(' *** ERROR *** ROUTINE REGLMR :',' SUM OF WEIGHTS IS NEGATIVE OR ZERO')
      END

!**********************************************************************************
!----------------------------------------------------------------------------------
!**********************************************************************************

SUBROUTINE PELWAK(XMOM,PARA,IFAIL)
!C***********************************************************************
!C*                                                                     *
!C*  FORTRAN CODE WRITTEN FOR INCLUSION IN IBM RESEARCH REPORT RC20525, *
!C*  'FORTRAN ROUTINES FOR USE WITH THE METHOD OF L-MOMENTS, VERSION 3' *
!C*                                                                     *
!C*  J. R. M. HOSKING                                                   *
!C*  IBM RESEARCH DIVISION                                              *
!C*  T. J. WATSON RESEARCH CENTER                                       *
!C*  YORKTOWN HEIGHTS                                                   *
!C*  NEW YORK 10598, U.S.A.                                             *
!C*                                                                     *
!C*  VERSION 3     AUGUST 1996                                          *
!C*                                                                     *
!C*  VERSION 3.04  JULY 2005                                            *
!C*  * Minor bug fix in test for validity of L-moments.                 *
!C*                                                                     *
!C***********************************************************************
!C
!C  PARAMETER ESTIMATION VIA L-MOMENTS FOR THE WAKEBY DISTRIBUTION
!C
!C  PARAMETERS OF ROUTINE:
!C  XMOM   * INPUT* ARRAY OF LENGTH 5. CONTAINS THE L-MOMENTS LAMBDA-1,
!C                  LAMBDA-2, TAU-3, TAU-4, TAU-5.
!C  PARA   *OUTPUT* ARRAY OF LENGTH 5. ON EXIT, CONTAINS THE PARAMETERS
!C                  IN THE ORDER XI, ALPHA, BETA, GAMMA, DELTA.
!C  IFAIL  *OUTPUT* FAIL FLAG. ON EXIT, IT IS SET AS FOLLOWS.
!C                  0 SUCCESSFUL EXIT
!C                  1 ESTIMATES COULD ONLY BE OBTAINED BY SETTING XI=0
!C                  2 ESTIMATES COULD ONLY BE OBTAINED BY FITTING A
!C                    GENERALIZED PARETO DISTRIBUTION
!C                  3 L-MOMENTS INVALID
!C
!C  PROCEDURE:
!C  1. LOOK FOR A SOLUTION WITH XI UNCONSTRAINED;
!C  2. IF NONE FOUND, LOOK FOR A SOLUTION WITH XI=0;
!C  3. IF NONE FOUND, FIT A GENERALIZED PARETO DISTRIBUTION TO THE
!C     FIRST 3 L-MOMENTS.
!C  ESTIMATES ARE CALCULATED USING THE FORMULAS GIVEN BY GREENWOOD ET AL.
!C  (1979, WATER RESOUR. RES., TABLE 5), BUT EXPRESSED IN TERMS OF
!C  L-MOMENTS RATHER THAN PROBABILITY WEIGHTED MOMENTS.
!C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION XMOM(5),PARA(5)
      DATA ZERO/0D0/,HALF/0.5D0/,ONE/1D0/,TWO/2D0/,THREE/3D0/,FOUR/4D0/
      DATA X2/2D0/,X3/3D0/,X4/4D0/,X5/5D0/,X7/7D0/,X8/8D0/,X9/9D0/,X10/10D0/,&
      X11/11D0/,X16/16D0/,X25/25D0/,X29/29D0/,X32/32D0/,X35/35D0/,X85/85D0/,X125/125D0/,X203/203D0/
!C
      IF(XMOM(2).LE.ZERO)GOTO 1000
      IF(DABS(XMOM(3)).GE.ONE)GOTO 1000
      IF(DABS(XMOM(4)).GE.ONE)GOTO 1000
      IF(DABS(XMOM(5)).GE.ONE)GOTO 1000
      IFAIL=0
!C
!C         CALCULATE THE L-MOMENTS (LAMBDA'S)
!C
      ALAM1=XMOM(1)
      ALAM2=XMOM(2)
      ALAM3=XMOM(3)*ALAM2
      ALAM4=XMOM(4)*ALAM2
      ALAM5=XMOM(5)*ALAM2
!C
!C         ESTIMATE N1,N2,N3,C1,C2,C3 WHEN XI.NE.0
!C
      XN1= X3*ALAM2-X25*ALAM3 +X32*ALAM4
      XN2=-X3*ALAM2 +X5*ALAM3  +X8*ALAM4
      XN3= X3*ALAM2 +X5*ALAM3  +X2*ALAM4
      XC1= X7*ALAM2-X85*ALAM3+X203*ALAM4-X125*ALAM5
      XC2=-X7*ALAM2+X25*ALAM3  +X7*ALAM4 -X25*ALAM5
      XC3= X7*ALAM2 +X5*ALAM3  -X7*ALAM4  -X5*ALAM5
!C
!C         ESTIMATE B AND D
!C
      XA=XN2*XC3-XC2*XN3
      XB=XN1*XC3-XC1*XN3
      XC=XN1*XC2-XC1*XN2
      DISC=XB*XB-FOUR*XA*XC
      IF(DISC.LT.ZERO)GOTO 10
      DISC=DSQRT(DISC)
      ROOT1=HALF*(-XB+DISC)/XA
      ROOT2=HALF*(-XB-DISC)/XA
      B= DMAX1(ROOT1,ROOT2)
      D=-DMIN1(ROOT1,ROOT2)
      IF(D.GE.ONE)GOTO 10
!C
!C         ESTIMATE A, C AND XI
!C
      A=(ONE+B)*(TWO+B)*(THREE+B)/(FOUR*(B+D))*((ONE+D)*ALAM2-(THREE-D)*ALAM3)
      C=-(ONE-D)*(TWO-D)*(THREE-D)/(FOUR*(B+D))*((ONE-B)*ALAM2-(THREE+B)*ALAM3)
      XI=ALAM1-A/(ONE+B)-C/(ONE-D)
!C
!C         CHECK FOR VALID PARAMETERS
!C
      IF(C.GE.ZERO.AND.A+C.GE.ZERO)GOTO 30
!C
!C         CAN'T FIND VALID ESTIMATES FOR XI UNRESTRICTED, SO TRY XI=0
!C
!C         ESTIMATE B AND D FOR XI=0
!C
   10 IFAIL=1
      XI=ZERO
      ZN1=X4*ALAM1-X11*ALAM2+X9*ALAM3
      ZN2=-ALAM2+X3*ALAM3
      ZN3=ALAM2+ALAM3
      ZC1=X10*ALAM1-X29*ALAM2+X35*ALAM3-X16*ALAM4
      ZC2=-ALAM2+X5*ALAM3-X4*ALAM4
      ZC3=ALAM2-ALAM4
      ZA=ZN2*ZC3-ZC2*ZN3
      ZB=ZN1*ZC3-ZC1*ZN3
      ZC=ZN1*ZC2-ZC1*ZN2
      DISC=ZB*ZB-FOUR*ZA*ZC
      IF(DISC.LT.ZERO)GOTO 20
      DISC=DSQRT(DISC)
      ROOT1=HALF*(-ZB+DISC)/ZA
      ROOT2=HALF*(-ZB-DISC)/ZA
      B= DMAX1(ROOT1,ROOT2)
      D=-DMIN1(ROOT1,ROOT2)
      IF(D.GE.ONE)GOTO 20
!C
!C         ESTIMATE A AND C
!C
      A= (ONE+B)*(TWO+B)/(B+D)*(ALAM1-(TWO-D)*ALAM2)
      C=-(ONE-D)*(TWO-D)/(B+D)*(ALAM1-(TWO+B)*ALAM2)
      IF(C.GE.ZERO.AND.A+C.GE.ZERO)GOTO 30
!C
!C         CAN'T FIND VALID ESTIMATES EVEN WITH XI=0 -
!C         FIT GENERALIZED PARETO DISTRIBUTION INSTEAD
!C
   20 IFAIL=2
      D=-(ONE-THREE*XMOM(3))/(ONE+XMOM(3))
      C=(ONE-D)*(TWO-D)*XMOM(2)
      B=ZERO
      A=ZERO
      XI=XMOM(1)-C/(ONE-D)
      IF(D.GT.ZERO)GOTO 30
      A=C
      B=-D
      C=ZERO
      D=ZERO

!C         COPY RESULTS INTO ARRAY PARA

   30 PARA(1)=XI
      PARA(2)=A
      PARA(3)=B
      PARA(4)=C
      PARA(5)=D
      RETURN

 1000 IFAIL=3
      DO 1010 I=1,5
 1010 PARA(I)=ZERO
      END


!*****************************************************************************************************
!-----------------------------------------------------------------------------------------------------
!*****************************************************************************************************
DOUBLE PRECISION FUNCTION QUAWAK(F,PARA)
!C***********************************************************************
!C*                                                                     *
!C*  FORTRAN CODE WRITTEN FOR INCLUSION IN IBM RESEARCH REPORT RC20525, *
!C*  'FORTRAN ROUTINES FOR USE WITH THE METHOD OF L-MOMENTS, VERSION 3' *
!C*                                                                     *
!C*  J. R. M. HOSKING                                                   *
!C*  IBM RESEARCH DIVISION                                              *
!C*  T. J. WATSON RESEARCH CENTER                                       *
!C*  YORKTOWN HEIGHTS                                                   *
!C*  NEW YORK 10598, U.S.A.                                             *
!C*                                                                     *
!C*  VERSION 3     AUGUST 1996                                          *
!C*                                                                     *
!C***********************************************************************
!C
!C  QUANTILE FUNCTION OF THE WAKEBY DISTRIBUTION
!C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION PARA(5)
      DATA ZERO/0D0/,ONE/1D0/
!C
!C         UFL SHOULD BE CHOSEN SO THAT EXP(UFL) JUST DOES NOT CAUSE
!C         UNDERFLOW
!C
      DATA UFL/-170D0/
!C
      XI=PARA(1)
      A=PARA(2)
      B=PARA(3)
      C=PARA(4)
      D=PARA(5)
!C
!C         TEST FOR VALID PARAMETERS
!!C
      IF(B+D.LE.ZERO.AND.(B.NE.ZERO.OR.C.NE.ZERO.OR.D.NE.ZERO))GOTO 1000
      IF(A.EQ.ZERO.AND.B.NE.ZERO)GOTO 1000
      IF(C.EQ.ZERO.AND.D.NE.ZERO)GOTO 1000
      IF(C.LT.ZERO.OR.A+C.LT.ZERO)GOTO 1000
      IF(A.EQ.ZERO.AND.C.EQ.ZERO)GOTO 1000

      IF(F.LE.ZERO.OR.F.GE.ONE)GOTO 10
      Z=-DLOG(ONE-F)
      Y1=Z
      IF(B.EQ.ZERO)GOTO 5
      TEMP=-B*Z
      IF(TEMP.LT.UFL)Y1=ONE/B
      IF(TEMP.GE.UFL)Y1=(ONE-DEXP(TEMP))/B
    5 CONTINUE
      Y2=Z
      IF(D.NE.ZERO)Y2=(ONE-DEXP(D*Y2))/(-D)
      QUAWAK=XI+A*Y1+C*Y2
      RETURN

   10 IF(F.EQ.ZERO)GOTO 20
      IF(F.EQ.ONE)GOTO 30
      GOTO 1010
   20 QUAWAK=XI
      RETURN
   30 IF(D.GT.ZERO)GOTO 1010
      IF(D.LT.ZERO)QUAWAK=XI+A/B-C/D
      IF(D.EQ.ZERO.AND.C.GT.ZERO)GOTO 1010
      IF(D.EQ.ZERO.AND.C.EQ.ZERO.AND.B.EQ.ZERO)GOTO 1010
      IF(D.EQ.ZERO.AND.C.EQ.ZERO.AND.B.GT.ZERO)QUAWAK=XI+A/B
      RETURN

 1000 WRITE(6,7000)
      QUAWAK=ZERO
      RETURN
 1010 WRITE(6,7010)
      QUAWAK=ZERO
      RETURN

 7000 FORMAT(' *** ERROR *** ROUTINE QUAWAK : PARAMETERS INVALID')
 7010 FORMAT(' *** ERROR *** ROUTINE QUAWAK :', ' ARGUMENT OF FUNCTION INVALID')
      END