 PROGRAM XCLUST
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
!C*  VERSION 3.02  MARCH 1997                                           *
!C*  * Minor change to FORMAT statement 6080                            *
!C*                                                                     *
!C*  VERSION 3.04  JULY 2005                                            *
!C*  * Removed declarations of unused variables                         *
!C*                                                                     *
!C***********************************************************************
!C
!C  EXAMPLE PROGRAM FOR CLUSTER ANALYSIS.  THE PROGRAM READS IN
!C  ATTRIBUTES FOR A NUMBER OF SITES, TRANSFORMS THE ATTRIBUTES, FORMS
!C  CLUSTERS BY WARD'S METHOD, PRINTS INFORMATION ABOUT THE CLUSTERS,
!C  AND REFINES THE CLUSTERS USING THE K-MEANS ALGORITHM.
!C    THE ANALYSIS FOLLOWS HOSKING AND WALLIS ("REGIONAL FREQUENCY
!C  ANALYSIS: AN APPROACH BASED ON L-MOMENTS", CAMBRIDGE UNIV. PRESS,
!C  1997, SEC. 9.2).
!C
!C  PARAMETERS OF PROGRAM:
!C  MAXNS  - MAXIMUM NUMBER OF SITES
!C  NATMAX - MAXIMUM NUMBER OF ATTRIBUTES
!C  INFILE - STREAM NUMBER TO WHICH INPUT FILE IS ATTACHED
!C  KOUT   - STREAM NUMBER TO WHICH OUTPUT FILE IS ATTACHED
!C
!C  ROUTINES USED: CLUAGG,CLUINF,CLUKM,KMNS,OPTRA,QTRAN
!C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXNS=104,NATMAX=10,INFILE=7,KOUT=6)
!C
      PARAMETER (NWORK=MAXNS*(MAXNS-1)/2,NRWORK=NWORK)
      CHARACTER*12 ID(MAXNS)
      INTEGER MERGE(2,MAXNS),IWORK(MAXNS),IASSGN(MAXNS),LIST(MAXNS),NUM(MAXNS)
      DOUBLE PRECISION X(MAXNS,NATMAX),WGSS(MAXNS),WORK(NWORK),WEIGHT(NATMAX),Z(NATMAX),CENT(NATMAX,MAXNS)
      REAL RWORK(NRWORK)
      DATA ZERO/0D0/,ONE/1D0/,THREE/3D0/
!C

      N=104
      OPEN(FILE="APPALACH.DAT",UNIT=20,status='unknown')
      !C
!C         FOR EACH DATA POINT ...APPALACH.DAT
!C
      DO 10 I=1,N
!C
!C         ... READ THE DATA ...
!C
      READ(20,5000)ID(I),XLAT,XLONG,AREA,ELEV
!C
!C         ... AND COMPUTE THE TRANSFORMED ATTRIBUTES
!C
      X(I,1)=DLOG(AREA)
      X(I,2)=DSQRT(ELEV)
      X(I,3)=XLAT
      X(I,4)=XLONG
!C
!C         END OF LOOP OVER DATA POINTS
!C
   10 CONTINUE
!C
!C         SET WEIGHTS FOR EACH ATTRIBUTE
!C
      NATT=4
      WEIGHT(1)=THREE
      WEIGHT(2)=ONE
      WEIGHT(3)=ONE
      WEIGHT(4)=ONE
!C
!C         FOR EACH ATTRIBUTE ...
!C
      DO 40 J=1,NATT
!C
!C         ... CALCULATE ITS STANDARD DEVIATION ACROSS THE DATA POINTS ...
!C
      SUM1=ZERO
      SUM2=ZERO
      DO 20 I=1,N
      SUM1=SUM1+X(I,J)
      SUM2=SUM2+X(I,J)**2
   20 CONTINUE
      SD=DSQRT((SUM2-SUM1*SUM1/N)/(N-ONE))
!C
!C         ... DIVIDE THE WEIGHT BY THIS STANDARD DEVIATION ...
!C
      WEIGHT(J)=WEIGHT(J)/SD
!C
!C         ... AND APPLY THE WEIGHT TO EACH DATA POINT
!C
      DO 30 I=1,N
   30 X(I,J)=X(I,J)*WEIGHT(J)
!C
!C         END OF LOOP OVER ATTRIBUTES
!C
   40 CONTINUE
!C
!C         WARD'S ALGORITHM
!C
      NX=MAXNS
      NW=NWORK
      CALL CLUAGG(3,X,NX,N,NATT,MERGE,WGSS,IWORK,WORK,NW)
      WRITE(KOUT,6000)
      DO 50 I=1,N-1
      WRITE(KOUT,6010)I,N-I,MERGE(1,I),MERGE(2,I),WGSS(I)
   50 CONTINUE
!C
!C         PRINT INFORMATION ABOUT THE 7-CLUSTER GROUPING
!C
      NCLUST=7
      CALL CLUINF(NCLUST,N,MERGE,IASSGN,LIST,NUM)
      WRITE(KOUT,6020)
      WRITE(KOUT,6030)(IASSGN(I),I=1,N)
      WRITE(KOUT,6040)
      IORIG=0
      DO 60 ICL=1,NCLUST
      NN=NUM(ICL)
      WRITE(KOUT,6050)ICL,NN
      WRITE(KOUT,6060)(IABS(LIST(I)),I=IORIG+1,IORIG+NN)
      IORIG=IORIG+NN
   60 CONTINUE
!C
!C         ADJUST CLUSTERS BY K-MEANS ALGORITHM
!C
      MAXIT=10
      CALL CLUKM(X,NX,N,NATT,NCLUST,IASSGN,LIST,NUM,SS,MAXIT,IWORK,RWORK,NRWORK)
!C
!C         PRINT INFORMATION ABOUT ADJUSTED CLUSTERS
!C
      WRITE(KOUT,6070)SS
      WRITE(KOUT,6020)
      WRITE(KOUT,6030)(IASSGN(I),I=1,N)
      WRITE(KOUT,6040)
      IORIG=0
      DO 70 ICL=1,NCLUST
      NN=NUM(ICL)
      WRITE(KOUT,6050)ICL,NN
      WRITE(KOUT,6060)(IABS(LIST(I)),I=IORIG+1,IORIG+NN)
      IORIG=IORIG+NN
   70 CONTINUE
!C
!C         FIND CLUSTER CENTERS, IN SPACE OF TRANSFORMED ATTRIBUTES
!C
      WRITE(KOUT,6080)
      DO 80 ICL=1,NCLUST
      DO 80 IATT=1,NATT
   80 CENT(IATT,ICL)=ZERO
      ICL=1
      DO 100 I=1,N
      L=LIST(I)
      ISITE=IABS(L)
      DO 90 IATT=1,NATT
      CENT(IATT,ICL)=CENT(IATT,ICL)+X(ISITE,IATT)
   90 CONTINUE
      IF(L.LT.0)ICL=ICL+1
  100 CONTINUE
      DO 110 ICL=1,NCLUST
      NN=NUM(ICL)
      DO 110 IATT=1,NATT
  110 CENT(IATT,ICL)=CENT(IATT,ICL)/NN
!C
!C         TRANSFORM BACK TO ORIGINAL ATTRIBUTES
!C
      DO 120 ICL=1,NCLUST
      Z(1)=DEXP(CENT(1,ICL)/WEIGHT(1))
      Z(2)=(CENT(2,ICL)/WEIGHT(2))**2
      Z(3)=CENT(3,ICL)/WEIGHT(3)
      Z(4)=CENT(4,ICL)/WEIGHT(4)
      WRITE(KOUT,6090)ICL,(Z(J),J=1,NATT)
  120 CONTINUE
!C
      STOP
!C
 5000 FORMAT(A8,4F8.0)
 6000 FORMAT(' MERGING SEQUENCE FROM WARD''S ALGORITHM'//)!Aqui eliminei algumas coisas p tentar rodar, porque de outra forma não tava dando p corrigir.
 6010 FORMAT(1X,I3,I9,I10,I5,F12.2)
 6020 FORMAT(/' ASSIGNMENT OF SITES TO CLUSTERS')
 6030 FORMAT(1X,10I4)
 6040 FORMAT(//' CLUSTER MEMBERSHIP')
 6050 FORMAT(/' CLUSTER',I4,'  HAS',I4,' MEMBERS:')
 6060 FORMAT(1X,10I4)
 6070 FORMAT(///' ADJUSTED CLUSTERS FROM K-MEANS ALGORITHM'/ ' (SUM OF SQUARES =',F12.2,')')
 6080 FORMAT(/' CLUSTER CENTERS'/ '          AREA      ELEV       LAT      LONG')
 6090 FORMAT(1X,I3,6F10.2)
      END PROGRAM XCLUST


!***************************************************************************************
!-------------------------------  Subrotinas  -----------------------------------------
!**************************************************************************************

SUBROUTINE CLUAGG(METHOD,X,NX,N,NATT,MERGE,DISP,IWORK,WORK,NW)
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
!C*  VERSION 3.02  MARCH 1997                                           *
!C*  * Implement single-link and complete-link clustering               *
!C*                                                                     *
!C***********************************************************************
!C
!C  CLUSTER ANALYSIS BY ANY OF SEVERAL AGGLOMERATIVE HIERARCHICAL METHODS
!C
!C  PARAMETERS OF ROUTINE:
!C  METHOD * INPUT* CLUSTERING METHOD. SHOULD BE SET TO:
!C                   1 FOR SINGLE-LINK CLUSTERING
!C                   2 FOR COMPLETE-LINK CLUSTERING
!C                   3 FOR WARD'S PROCEDURE
!C  X      * INPUT* ARRAY OF DIMENSION (NX,NATT).  X(I,J) SHOULD CONTAIN
!C                  THE J'TH ATTRIBUTE FOR THE I'TH DATA POINT.
!C  NX     * INPUT* THE FIRST DIMENSION OF ARRAY X, AS DECLARED IN THE
!C                  CALLING PROGRAM.
!C  N      * INPUT* NUMBER OF DATA POINTS
!C  NATT   * INPUT* NUMBER OF ATTRIBUTES FOR EACH DATA POINT
!C  MERGE  *OUTPUT* ARRAY OF DIMENSION (2,N). MERGE(1,I) AND MERGE(2,I)
!C                  ARE THE LABELS OF THE CLUSTERS MERGED AT THE I'TH
!C                  STAGE.  MERGE(1,N) AND MERGE(2,N) ARE NOT USED.
!C  DISP   *OUTPUT* ARRAY OF LENGTH N.  DISP(I) IS A MEASURE OF THE
!C                  WITHIN-CLUSTER DISPERSION AFTER THE I'TH MERGE.
!C                  DISPERSION IS DEFINED DIFFERENTLY FOR EACH METHOD:
!C                  SEE BELOW.  DISP(N) IS NOT USED.
!C  IWORK  * LOCAL* WORK ARRAY OF LENGTH N
!C  WORK   * LOCAL* WORK ARRAY OF LENGTH NW
!C  NW     * INPUT* LENGTH OF ARRAY WORK. MUST BE AT LEAST N*(N-1)/2.
!C
!C  Agglomerative hierarchical clustering: general description.
!C  Initially there are N clusters, each containing one data point,
!C  labeled 1 through N in the same order as the data points.  At each
!C  stage of clustering, two clusters are merged.  Their labels are saved
!C  in the MERGE array.  The smaller of the two labels is used as the
!C  label of the merged cluster.  After the Mth stage of clustering
!C  there are N-M clusters.  To find which data points belong to which
!C  clusters, use routine CLUINF.
!C
!C  Single-link clustering: the distance between two clusters A and B is
!C  defined to be the minimum of the Euclidean distances between pairs of
!C  points with one point in A and one in B.  At each stage, the two
!C  clusters separated by the smallest distance are merged.  The square
!C  of this distance is saved in the corresponding element of array DISP.
!C
!C  Complete-link clustering: the distance between two clusters A and B
!C  is defined to be the maximum of the Euclidean distances between pairs
!C  of points with one point in A and one in B.  At each stage, the two
!C  clusters separated by the smallest distance are merged.  The square
!C  of this distance is saved in the corresponding element of array DISP.
!C  DISP(I) is therefore the largest squared Euclidean distance between
!C  two points that are in the same cluster after the Ith merge.
!C
!C  Ward's procedure: at each stage, the clusters that are merged are
!C  chosen to minimize the within-cluster sum of squared deviations of
!C  each attribute about the cluster mean.  This sum of squares is saved
!C  in the corresponding element of array DISP.
!C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X(NX,NATT),DISP(N),WORK(NW)
      INTEGER MERGE(2,N),IWORK(N)
      DATA ZERO/0D0/,HALF/0.5D0/
!C
!C         BIG IS A LARGE NUMBER, USED TO INITIALIZE THE SEARCH CRITERION
!C
      DATA BIG/1D72/
!C
      NWREQ=N*(N-1)/2
      IF(NW.LT.NWREQ)GOTO 1000
!C
!C         INITIALLY THERE ARE N CLUSTERS, EACH CONTAINING ONE DATA
!C         POINT.  COMPUTE THE COST (INCREASE IN DISPERSION) OF MERGING
!C         EACH PAIR OF CLUSTERS.
!C
      IW=0
      DO 20 J=2,N
      DO 20 I=1,J-1
      SUM=ZERO
      DO 10 IATT=1,NATT
   10 SUM=SUM+(X(I,IATT)-X(J,IATT))**2
      IW=IW+1
      WORK(IW)=SUM
   20 CONTINUE
      DO 30 I=1,N
   30 IWORK(I)=1
      CCOST=ZERO
!C
!C         START OF MAIN LOOP
!C
      DO 100 IMERGE=1,N-1
!C
!C         FIND THE PAIR OF CLUSTERS WITH THE LOWEST COST OF MERGING
!C
      COST=BIG
      DO 50 J=2,N
      IF(IWORK(J).EQ.0)GOTO 50
      IORIG=(J-1)*(J-2)/2
      DO 40 I=1,J-1
      IF(IWORK(I).EQ.0)GOTO 40
      LIJ=IORIG+I
      IF(WORK(LIJ).GE.COST)GOTO 40
      COST=WORK(LIJ)
      II=I
      JJ=J
   40 CONTINUE
   50 CONTINUE
!C
!C         MERGE THEM
!C
      MERGE(1,IMERGE)=II
      MERGE(2,IMERGE)=JJ
      IF(METHOD.EQ.1.OR.METHOD.EQ.2)DISP(IMERGE)=COST
      IF(METHOD.EQ.3)CCOST=CCOST+COST
      IF(METHOD.EQ.3)DISP(IMERGE)=HALF*CCOST
!C
!C         COMPUTE THE COST OF MERGING THE NEW CLUSTER WITH EACH OF THE
!C         OTHERS
!C
      NI=IWORK(II)
      NJ=IWORK(JJ)
      NIJ=NI+NJ
      DO 60 KK=1,N
      NK=IWORK(KK)
      IF(NK.EQ.0)GOTO 60
      IF(KK.EQ.II.OR.KK.EQ.JJ)GOTO 60
      MM=MAX0(II,KK)
      M =MIN0(II,KK)
      IK=(MM-1)*(MM-2)/2+M
      MM=MAX0(JJ,KK)
      M =MIN0(JJ,KK)
      JK=(MM-1)*(MM-2)/2+M
      IF(METHOD.EQ.1)WORK(IK)=DMIN1(WORK(IK),WORK(JK))
      IF(METHOD.EQ.2)WORK(IK)=DMAX1(WORK(IK),WORK(JK))
      IF(METHOD.EQ.3)WORK(IK)=((NI+NK)*WORK(IK)+(NJ+NK)*WORK(JK)-NK*COST)/(NIJ+NK)
   60 CONTINUE
      IWORK(II)=NIJ
      IWORK(JJ)=0
!C
!C         END OF MAIN LOOP
!C
  100 CONTINUE
!C
      RETURN
!C
 1000 WRITE(6,7000)NWREQ
      RETURN
!C
 7000 FORMAT(' *** ERROR *** ROUTINE CLUAGG : INSUFFICIENT WORKSPACE.', ' LENGTH OF WORK ARRAY SHOULD BE AT LEAST ',I8)
!C
      END
      
!C=====================================================

 
SUBROUTINE CLUINF(NCLUST,N,MERGE,IASSGN,LIST,NUM)
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
!C*  VERSION 3.02  MARCH 1997                                           *
!C*  * Check for N.LT.NCLUST                                            *
!C*  * Minor changes to comments                                        *
!C*                                                                     *
!C***********************************************************************
!C
!C  OBTAINS INFORMATION ABOUT CLUSTERS ARISING FROM AGGLOMERATIVE
!C  HIERARCHICAL CLUSTERING
!C
!C  AGGLOMERATIVE HIERARCHICAL CLUSTERING PROCEDURES TYPICALLY PRODUCE A
!C  LIST OF THE CLUSTERS MERGED AT EACH STAGE OF THE CLUSTERING.  THIS
!C  ROUTINE USES THIS LIST TO CONSTRUCT ARRAYS THAT EXPLICITLY SHOW
!C  WHICH CLUSTER A GIVEN DATA POINT BELONGS TO, AND WHICH DATA POINTS
!C  BELONG TO A GIVEN CLUSTER.
!C
!C  PARAMETERS OF ROUTINE:
!C  NCLUST * INPUT* NUMBER OF CLUSTERS
!C  N      * INPUT* NUMBER OF DATA POINTS
!C  MERGE  * INPUT* ARRAY OF DIMENSION (2,N). MERGE(1,I) AND MERGE(2,I)
!C                  IDENTIFY THE CLUSTERS MERGED AT THE I'TH STEP.
!C                  THIS IS THE ARRAY MERGE RETURNED BY ROUTINE CLUAGG,
!C                  AND SHOULD BE LEFT UNCHANGED AFTER EXIT FROM THAT
!C                  ROUTINE.
!C  IASSGN *OUTPUT* ARRAY OF LENGTH N. ITS I'TH ELEMENT IS THE NUMBER
!C                  OF THE CLUSTER TO WHICH THE I'TH DATA POINT BELONGS.
!C  LIST   *OUTPUT* ARRAY OF LENGTH N. CONTAINS THE DATA POINTS IN
!C                  CLUSTER 1, FOLLOWED BY THE DATA POINTS IN CLUSTER 2,
!C                  ETC.  DATA POINTS IN EACH CLUSTER ARE LISTED IN
!C                  INCREASING ORDER.  THE LAST DATA POINT IN EACH
!C                  CLUSTER IS INDICATED BY A NEGATIVE NUMBER.
!C                  SEE THE EXAMPLE BELOW.
!C  NUM    *OUTPUT* ARRAY OF LENGTH NCLUST.  NUMBER OF DATA POINTS IN
!C                  EACH CLUSTER.
!C
!C  CLUSTER NUMBERS USED IN ARRAYS IASSGN, LIST AND NUM RANGE FROM 1 TO
!C  NCLUST.  THEY ARE ARBITRARY, BUT ARE UNIQUELY DEFINED: CLUSTER 1
!C  CONTAINS DATA POINT 1, CLUSTER M (M.GE.2) CONTAINS DATA POINT J,
!C  WHERE J=MERGE(2,N-M).
!C
!C  EXAMPLE OF THE LIST ARRAY.  SUPPOSE THAT THERE ARE 8 DATA POINTS
!C  AND 3 CLUSTERS, AND THAT THE ELEMENTS OF THE LIST ARRAY ARE
!C  1, -4, 3, 6, -8, 2, 5, -7.  THEN THE CLUSTERS ARE AS FOLLOWS:
!C  CLUSTER 1 CONTAINS POINTS 1 AND 4; CLUSTER 2 CONTAINS POINTS
!C  3, 6 AND 8; CLUSTER 3 CONTAINS POINTS 2, 5 AND 7.
!C
      INTEGER MERGE(2,N),IASSGN(N),LIST(N),NUM(NCLUST)
      IF(N.LT.NCLUST)GOTO 1000
!C
!C       CONSTRUCT THE IASSGN ARRAY
!C
      IASSGN(1)=1
      DO 10 I=1,NCLUST-1
      ITEMP=MERGE(2,N-I)
      IASSGN(ITEMP)=I+1
   10 CONTINUE
      DO 20 I=NCLUST,N-1
      ICL=N-I
      II=MERGE(1,ICL)
      JJ=MERGE(2,ICL)
      IASSGN(JJ)=IASSGN(II)
   20 CONTINUE
!C
!C       CONSTRUCT THE LIST AND NUM ARRAYS
!C
      LASTI=0
      I=0
      DO 70 ICL=1,NCLUST
      DO 60 K=1,N
      IF(IASSGN(K).NE.ICL)GOTO 60
      I=I+1
      LIST(I)=K
   60 CONTINUE
      LIST(I)=-LIST(I)
      NUM(ICL)=I-LASTI
      LASTI=I
   70 CONTINUE
      RETURN
!C
 1000 WRITE(6,7000)
      RETURN
!C
 7000 FORMAT(' *** ERROR *** ROUTINE CLUINF :', ' NUMBER OF CLUSTERS EXCEEDS NUMBER OF DATA POINTS')
!C
      END
      
!C=====================================================
      
SUBROUTINE CLUKM(X,NX,N,NATT,NCLUST,IASSGN,LIST,NUM,SS,MAXIT,IWORK,RW,NW)
!C***********************************************************************
!!C*                                                                    *
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
!C  CLUSTER ANALYSIS BY THE K-MEANS ALGORITHM
!C
!C  PARAMETERS OF ROUTINE:
!C  X      * INPUT* ARRAY OF DIMENSION (NX,NATT).  X(I,J) SHOULD
!C                  CONTAIN THE J'TH ATTRIBUTE FOR THE I'TH DATA POINT.
!C  NX     * INPUT* THE FIRST DIMENSION OF ARRAY X, AS DECLARED IN THE
!C                  CALLING PROGRAM.
!C  N      * INPUT* NUMBER OF DATA POINTS
!C  NATT   * INPUT* NUMBER OF ATTRIBUTES FOR EACH DATA POINT
!C  NCLUST * INPUT* NUMBER OF CLUSTERS
!C  IASSGN *IN/OUT* ARRAY OF LENGTH N.  ON ENTRY, SHOULD CONTAIN THE
!C                  INITIAL ASSIGNMENT OF SITES TO CLUSTERS.  ON EXIT,
!C                  CONTAINS THE FINAL ASSIGNMENT.  THE I'TH ELEMENT OF
!C                  THE ARRAY CONTAINS THE LABEL OF THE CLUSTER TO WHICH
!C                  THE I'TH DATA POINT BELONGS.  LABELS MUST BE BETWEEN
!C                  1 AND NCLUST, AND EACH OF THE VALUES 1 THROUGH NCLUST
!C                  MUST OCCUR AT LEAST ONCE.
!C  LIST   *OUTPUT* ARRAY OF LENGTH N. CONTAINS THE DATA POINTS IN
!C                  CLUSTER 1, FOLLOWED BY THE DATA POINTS IN CLUSTER 2,
!C                  ETC.  DATA POINTS IN EACH CLUSTER ARE LISTED IN
!C                  INCREASING ORDER.  THE LAST DATA POINT IN EACH
!C                  CLUSTER IS INDICATED BY A NEGATIVE NUMBER.
!C  NUM    *OUTPUT* ARRAY OF LENGTH NCLUST.  NUMBER OF DATA POINTS IN
!C                  EACH CLUSTER.
!C  SS     *OUTPUT* WITHIN-GROUP SUM OF SQUARES OF THE FINAL CLUSTERS.
!C  MAXIT  * INPUT* MAXIMUM NUMBER OF ITERATIONS FOR THE K-MEANS
!C                  CLUSTERING ALGORITHM
!C  IWORK  * LOCAL* (INTEGER) WORK ARRAY OF LENGTH NCLUST*3
!C  RW     * LOCAL* REAL WORK ARRAY OF LENGTH NW.  N.B. THIS ARRAY IS OF
!C                  TYPE REAL, NOT DOUBLE PRECISION!
!C  NW     * INPUT* LENGTH OF ARRAY RW.  MUST BE AT LEAST
!C                  (N+NCLUST)*(NATT+1)+2*NCLUST
!C
!C  OTHER ROUTINES USED: APPLIED STATISTICS ALGORITHM AS136 (ROUTINES
!C                       KMNS,OPTRA,QTRAN), AVAILABLE FROM
!C                       HTTP://STAT.LIB.CMU.EDU/APSTAT/136
!C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X(NX,NATT)
      INTEGER IASSGN(N),LIST(N),NUM(NCLUST),IWORK(NCLUST,3)
      REAL RW(NW)
      DATA ZERO/0D0/
!C
!C         SET ADDRESSES FOR SUBDIVISION OF WORK ARRAY
!C
      MC=1
      MA=MC+NCLUST*NATT
      MAN1=MA+N*NATT
      MAN2=MAN1+NCLUST
      MWSS=MAN2+NCLUST
      MD=MWSS+NCLUST
      NWREQ=MD+N-1
      IF(NW.LT.NWREQ)GOTO 1000
      LA=MA-1
      LWSS=MWSS-1
!C
!C         COPY ATTRIBUTES TO WORK ARRAY
!C
      IW=LA
      DO 5 IATT=1,NATT
      DO 5 I=1,N
      IW=IW+1
    5 RW(IW)=X(I,IATT)
!C
!C         COMPUTE CLUSTER CENTERS
!C
      DO 10 ICL=1,NCLUST
   10 NUM(ICL)=0
      IWMAX=NCLUST*NATT
      DO 20 IW=1,IWMAX
   20 RW(IW)=ZERO
      DO 40 I=1,N
      ICL=IASSGN(I)
      IF(ICL.LE.0.OR.ICL.GT.NCLUST)GOTO 1010
      NUM(ICL)=NUM(ICL)+1
      IW=ICL
      DO 30 IATT=1,NATT
      RW(IW)=RW(IW)+X(I,IATT)
      IW=IW+NCLUST
   30 CONTINUE
   40 CONTINUE
      DO 60 ICL=1,NCLUST
      NSIZE=NUM(ICL)
      IF(NSIZE.EQ.0)GOTO 1020
      IW=ICL
      DO 50 IATT=1,NATT
      RW(IW)=RW(IW)/NSIZE
      IW=IW+NCLUST
   50 CONTINUE
   60 CONTINUE
!C
!C         CALL ALGORITHM AS136
!C
      CALL KMNS(RW(MA),N,NATT,RW(MC),NCLUST,IASSGN,LIST,NUM,RW(MAN1),RW(MAN2),&
      &IWORK(1,1),RW(MD),IWORK(1,2),IWORK(1,3),MAXIT,RW(MWSS),IFAULT)
      IF(IFAULT.EQ.2)WRITE(6,7030)
!C
!C         COMPUTE LIST ARRAY AND FINAL SUM OF SQUARES
!C
      I=0
      DO 80 ICL=1,NCLUST
      DO 70 K=1,N
      IF(IASSGN(K).NE.ICL)GOTO 70
      I=I+1
      LIST(I)=K
   70 CONTINUE
      LIST(I)=-LIST(I)
   80 CONTINUE
      SS=ZERO
      DO 90 ICL=1,NCLUST
   90 SS=SS+RW(LWSS+ICL)
!C
      RETURN
!C
 1000 WRITE(6,7000)NWREQ
      RETURN
 1010 WRITE(6,7010)I
      RETURN
 1020 WRITE(6,7020)ICL
      RETURN
!C
 7000 FORMAT(' *** ERROR *** ROUTINE CLUKM  : INSUFFICIENT WORKSPACE.',  ' LENGTH OF WORK ARRAY SHOULD BE AT LEAST ',I8)
 7010 FORMAT(' *** ERROR *** ROUTINE CLUKM  :',  ' INVALID INITIAL CLUSTER NUMBER FOR DATA POINT ',I5)
 7020 FORMAT(' *** ERROR *** ROUTINE CLUKM  :',  ' INITIAL CLUSTERS INVALID.  CLUSTER ',I4,' HAS NO MEMBERS.')
 7030 FORMAT(' ** WARNING ** ROUTINE CLUKM  :',  ' ITERATION HAS NOT CONVERGED. RESULTS MAY BE UNRELIABLE.')
!C
      END

 !*********************************************************************************

SUBROUTINE KMNS(A, M, N, C, K, IC1, IC2, NC, AN1, AN2, NCP, D,ITRAN, LIVE, ITER, WSS, IFAULT)

!C     ALGORITHM AS 136  APPL. STATIST. (1979) VOL.28, NO.1
!C
!C     Divide M points in N-dimensional space into K clusters so that
!C     the within cluster sum of squares is minimized.
!C
      INTEGER IC1(M), IC2(M), NC(K), NCP(K), ITRAN(K), LIVE(K)
      REAL    A(M,N), D(M), C(K,N), AN1(K), AN2(K), WSS(K), DT(2)
      REAL    ZERO, ONE
!C
!C     Define BIG to be a very large positive number
!C
      DATA BIG /1.E30/, ZERO /0.0/, ONE /1.0/
!C
      IFAULT = 3
      IF (K .LE. 1 .OR. K .GE. M) RETURN
!C
!C     For each point I, find its two closest centres, IC1(I) and
!C     IC2(I).     Assign it to IC1(I).
!C
      DO 50 I = 1, M
	IC1(I) = 1
	IC2(I) = 2
	DO 10 IL = 1, 2
	  DT(IL) = ZERO
	  DO 10 J = 1, N
	    DA = A(I,J) - C(IL,J)
	    DT(IL) = DT(IL) + DA*DA
   10   CONTINUE
	IF (DT(1) .GT. DT(2)) THEN
	  IC1(I) = 2
	  IC2(I) = 1
	  TEMP = DT(1)
	  DT(1) = DT(2)
	  DT(2) = TEMP
	END IF
	DO 50 L = 3, K
	  DB = ZERO
	  DO 30 J = 1, N
	    DC = A(I,J) - C(L,J)
	    DB = DB + DC*DC
	    IF (DB .GE. DT(2)) GO TO 50
   30     CONTINUE
	  IF (DB .LT. DT(1)) GO TO 40
	  DT(2) = DB
	  IC2(I) = L
	  GO TO 50
   40     DT(2) = DT(1)
	  IC2(I) = IC1(I)
	  DT(1) = DB
	  IC1(I) = L
   50 CONTINUE
!C
!C     Update cluster centres to be the average of points contained
!C     within them.
!C
      DO 70 L = 1, K
	NC(L) = 0
	DO 60 J = 1, N
   60   C(L,J) = ZERO
   70 CONTINUE
      DO 90 I = 1, M
	L = IC1(I)
	NC(L) = NC(L) + 1
	DO 80 J = 1, N
   80   C(L,J) = C(L,J) + A(I,J)
   90 CONTINUE
!C
!C     Check to see if there is any empty cluster at this stage
!C
      DO 120 L = 1, K
	IF (NC(L) .EQ. 0) THEN
	  IFAULT = 1
	  RETURN
	END IF
	AA = NC(L)
	DO 110 J = 1, N
  110   C(L,J) = C(L,J) / AA
!C
!C     Initialize AN1, AN2, ITRAN & NCP
!C     AN1(L) = NC(L) / (NC(L) - 1)
!C     AN2(L) = NC(L) / (NC(L) + 1)
!C     ITRAN(L) = 1 if cluster L is updated in the quick-transfer stage,
!C              = 0 otherwise
!C     In the optimal-transfer stage, NCP(L) stores the step at which
!C     cluster L is last updated.
!C     In the quick-transfer stage, NCP(L) stores the step at which
!C     cluster L is last updated plus M.
!C
	AN2(L) = AA / (AA + ONE)
	AN1(L) = BIG
	IF (AA .GT. ONE) AN1(L) = AA / (AA - ONE)
	ITRAN(L) = 1
	NCP(L) = -1
  120 CONTINUE
      INDX = 0
      DO 140 IJ = 1, ITER
!C
!C     In this stage, there is only one pass through the data.   Each
!C     point is re-allocated, if necessary, to the cluster that will
!C     induce the maximum reduction in within-cluster sum of squares.
!C
	CALL OPTRA(A, M, N, C, K, IC1, IC2, NC, AN1, AN2, NCP, D, ITRAN, LIVE, INDX)
!C
!C     Stop if no transfer took place in the last M optimal transfer
!C     steps.
!C
	IF (INDX .EQ. M) GO TO 150
!C
!C     Each point is tested in turn to see if it should be re-allocated
!C     to the cluster to which it is most likely to be transferred,
!C     IC2(I), from its present cluster, IC1(I).   Loop through the
!C     data until no further change is to take place.
!C
	CALL QTRAN(A, M, N, C, K, IC1, IC2, NC, AN1, AN2, NCP, D,ITRAN, INDX)
!C
!C     If there are only two clusters, there is no need to re-enter the
!C     optimal transfer stage.
!C
	IF (K .EQ. 2) GO TO 150
!C
!C     NCP has to be set to 0 before entering OPTRA.
!C
	DO 130 L = 1, K
  130   NCP(L) = 0
  140 CONTINUE
!C
!C     Since the specified number of iterations has been exceeded, set
!C     IFAULT = 2.   This may indicate unforeseen looping.
!C
      IFAULT = 2
!C
!C     Compute within-cluster sum of squares for each cluster.
!C
  150 DO 160 L = 1, K
	WSS(L) = ZERO
	DO 160 J = 1, N
	  C(L,J) = ZERO
  160 CONTINUE
      DO 170 I = 1, M
	II = IC1(I)
	DO 170 J = 1, N
	  C(II,J) = C(II,J) + A(I,J)
  170 CONTINUE
      DO 190 J = 1, N
	DO 180 L = 1, K
  180   C(L,J) = C(L,J) / FLOAT(NC(L))
	DO 190 I = 1, M
	  II = IC1(I)
	  DA = A(I,J) - C(II,J)
	  WSS(II) = WSS(II) + DA*DA
  190 CONTINUE
!C
      RETURN
      END
!C
!C
      SUBROUTINE OPTRA(A, M, N, C, K, IC1, IC2, NC, AN1, AN2, NCP, D,ITRAN, LIVE, INDX)
!C
!C     ALGORITHM AS 136.1  APPL. STATIST. (1979) VOL.28, NO.1
!C
!C     This is the optimal transfer stage.
!C
!C     Each point is re-allocated, if necessary, to the cluster that
!C     will induce a maximum reduction in the within-cluster sum of
!C     squares.
!C
      INTEGER IC1(M), IC2(M), NC(K), NCP(K), ITRAN(K), LIVE(K)
      REAL    A(M,N), D(M), C(K,N), AN1(K), AN2(K), ZERO, ONE
!C
!C     Define BIG to be a very large positive number.
!C
      DATA BIG /1.0E30/, ZERO /0.0/, ONE/1.0/
!C
!C     If cluster L is updated in the last quick-transfer stage, it
!C     belongs to the live set throughout this stage.   Otherwise, at
!C     each step, it is not in the live set if it has not been updated
!C     in the last M optimal transfer steps.
!C
      DO 10 L = 1, K
	IF (ITRAN(L) .EQ. 1) LIVE(L) = M + 1
   10 CONTINUE
      DO 100 I = 1, M
	INDX = INDX + 1
	L1 = IC1(I)
	L2 = IC2(I)
	LL = L2
!C
!C     If point I is the only member of cluster L1, no transfer.
!C
	IF (NC(L1) .EQ. 1) GO TO 90
!C
!C     If L1 has not yet been updated in this stage, no need to
!C     re-compute D(I).
!C
	IF (NCP(L1) .EQ. 0) GO TO 30
	DE = ZERO
	DO 20 J = 1, N
	  DF = A(I,J) - C(L1,J)
	  DE = DE + DF*DF
   20   CONTINUE
	D(I) = DE * AN1(L1)
!C
!C     Find the cluster with minimum R2.
!C
   30   DA = ZERO
	DO 40 J = 1, N
	  DB = A(I,J) - C(L2,J)
	  DA = DA + DB*DB
   40   CONTINUE
	R2 = DA * AN2(L2)
	DO 60 L = 1, K
!C
!C     If I >= LIVE(L1), then L1 is not in the live set.   If this is
!C     true, we only need to consider clusters that are in the live set
!C     for possible transfer of point I.   Otherwise, we need to consider
!C     all possible clusters.
!C
	  IF (I .GE. LIVE(L1) .AND. I .GE. LIVE(L) .OR. L .EQ. L1 .OR. L .EQ. LL) GO TO 60
	  RR = R2 / AN2(L)
	  DC = ZERO
	  DO 50 J = 1, N
	    DD = A(I,J) - C(L,J)
	    DC = DC + DD*DD
	    IF (DC .GE. RR) GO TO 60
   50     CONTINUE
	  R2 = DC * AN2(L)
	  L2 = L
   60     CONTINUE
	  IF (R2 .LT. D(I)) GO TO 70
!C
!C     If no transfer is necessary, L2 is the new IC2(I).
!C
	  IC2(I) = L2
	  GO TO 90
!C
!C     Update cluster centres, LIVE, NCP, AN1 & AN2 for clusters L1 and
!C     L2, and update IC1(I) & IC2(I).
!C
   70     INDX = 0
	  LIVE(L1) = M + I
	  LIVE(L2) = M + I
	  NCP(L1) = I
	  NCP(L2) = I
	  AL1 = NC(L1)
	  ALW = AL1 - ONE
	  AL2 = NC(L2)
	  ALT = AL2 + ONE
	  DO 80 J = 1, N
	    C(L1,J) = (C(L1,J) * AL1 - A(I,J)) / ALW
	    C(L2,J) = (C(L2,J) * AL2 + A(I,J)) / ALT
   80     CONTINUE
	  NC(L1) = NC(L1) - 1
	  NC(L2) = NC(L2) + 1
	  AN2(L1) = ALW / AL1
	  AN1(L1) = BIG
	  IF (ALW .GT. ONE) AN1(L1) = ALW / (ALW - ONE)
	  AN1(L2) = ALT / AL2
	  AN2(L2) = ALT / (ALT + ONE)
	  IC1(I) = L2
	  IC2(I) = L1
   90   CONTINUE
	IF (INDX .EQ. M) RETURN
  100 CONTINUE
      DO 110 L = 1, K
!C
!C     ITRAN(L) = 0 before entering QTRAN.   Also, LIVE(L) has to be
!C     decreased by M before re-entering OPTRA.
!C
	ITRAN(L) = 0
	LIVE(L) = LIVE(L) - M
  110 CONTINUE
!C
      RETURN
      END
!C
!C
SUBROUTINE QTRAN(A, M, N, C, K, IC1, IC2, NC, AN1, AN2, NCP, D,ITRAN, INDX)
!C
!C     ALGORITHM AS 136.2  APPL. STATIST. (1979) VOL.28, NO.1
!C
!C     This is the quick transfer stage.
!C     IC1(I) is the cluster which point I belongs to.
!C     IC2(I) is the cluster which point I is most likely to be
!C         transferred to.
!C     For each point I, IC1(I) & IC2(I) are switched, if necessary, to
!C     reduce within-cluster sum of squares.  The cluster centres are
!C     updated after each step.
!C
      INTEGER IC1(M), IC2(M), NC(K), NCP(K), ITRAN(K)
      REAL    A(M,N), D(M), C(K,N), AN1(K), AN2(K), ZERO, ONE
!C
!C     Define BIG to be a very large positive number
!C
      DATA BIG /1.0E30/, ZERO /0.0/, ONE /1.0/
!C
!C     In the optimal transfer stage, NCP(L) indicates the step at which
!C     cluster L is last updated.   In the quick transfer stage, NCP(L)
!C     is equal to the step at which cluster L is last updated plus M.
!C
      ICOUN = 0
      ISTEP = 0
   10 DO 70 I = 1, M
	ICOUN = ICOUN + 1
	ISTEP = ISTEP + 1
	L1 = IC1(I)
	L2 = IC2(I)
!C
!C     If point I is the only member of cluster L1, no transfer.
!C
	IF (NC(L1) .EQ. 1) GO TO 60
!C
!C     If ISTEP > NCP(L1), no need to re-compute distance from point I to
!C     cluster L1.   Note that if cluster L1 is last updated exactly M
!C     steps ago, we still need to compute the distance from point I to
!C     cluster L1.
!C
	IF (ISTEP .GT. NCP(L1)) GO TO 30
	DA = ZERO
	DO 20 J = 1, N
	  DB = A(I,J) - C(L1,J)
	  DA = DA + DB*DB
   20   CONTINUE
	D(I) = DA * AN1(L1)
!C
!C     If ISTEP >= both NCP(L1) & NCP(L2) there will be no transfer of
!C     point I at this step.
!C
   30   IF (ISTEP .GE. NCP(L1) .AND. ISTEP .GE. NCP(L2)) GO TO 60
	R2 = D(I) / AN2(L2)
	DD = ZERO
	DO 40 J = 1, N
	  DE = A(I,J) - C(L2,J)
	  DD = DD + DE*DE
	  IF (DD .GE. R2) GO TO 60
   40   CONTINUE
!C
!C     Update cluster centres, NCP, NC, ITRAN, AN1 & AN2 for clusters
!C     L1 & L2.   Also update IC1(I) & IC2(I).   Note that if any
!C     updating occurs in this stage, INDX is set back to 0.
!C
	ICOUN = 0
	INDX = 0
	ITRAN(L1) = 1
	ITRAN(L2) = 1
	NCP(L1) = ISTEP + M
	NCP(L2) = ISTEP + M
	AL1 = NC(L1)
	ALW = AL1 - ONE
	AL2 = NC(L2)
	ALT = AL2 + ONE
	DO 50 J = 1, N
	  C(L1,J) = (C(L1,J) * AL1 - A(I,J)) / ALW
	  C(L2,J) = (C(L2,J) * AL2 + A(I,J)) / ALT
   50   CONTINUE
	NC(L1) = NC(L1) - 1
	NC(L2) = NC(L2) + 1
	AN2(L1) = ALW / AL1
	AN1(L1) = BIG
	IF (ALW .GT. ONE) AN1(L1) = ALW / (ALW - ONE)
	AN1(L2) = ALT / AL2
	AN2(L2) = ALT / (ALT + ONE)
	IC1(I) = L2
	IC2(I) = L1
!C
!C     If no re-allocation took place in the last M steps, return.
!C
   60   IF (ICOUN .EQ. M) RETURN
   70 CONTINUE
      GO TO 10
END