 
      SUBROUTINE DB2INK(X,NX,Y,NY,FCN,LDF,KX,KY,TX,TY,BCOEF,WORK,IFLAG)
C***BEGIN PROLOGUE  DB2INK
C***DATE WRITTEN   25 MAY 1982
C***REVISION DATE  25 MAY 1982
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***CATEGORY NO.  E1A
C***KEYWORDS  INTERPOLATION, TWO-DIMENSIONS, GRIDDED DATA, SPLINES,
C             PIECEWISE POLYNOMIALS
C***AUTHOR  BOISVERT, RONALD, NBS
C             SCIENTIFIC COMPUTING DIVISION
C             NATIONAL BUREAU OF STANDARDS
C             WASHINGTON, DC 20234
C***PURPOSE  DOUBLE PRECISION VERSION OF B2INK.
C            DB2INK DETERMINES A PIECEWISE POLYNOMIAL FUNCTION THAT
C            INTERPOLATES TWO-DIMENSIONAL GRIDDED DATA. USERS SPECIFY
C            THE POLYNOMIAL ORDER (DEGREE+1) OF THE INTERPOLANT AND
C            (OPTIONALLY) THE KNOT SEQUENCE.
C***DESCRIPTION
C
C   DB2INK determines the parameters of a  function  that  interpolates
C   the two-dimensional gridded data (X(i),Y(j),FCN(i,j)) for i=1,..,NX
C   and j=1,..,NY. The interpolating function and its  derivatives  may
C   subsequently be evaluated by the function DB2VAL.
C
C   The interpolating  function  is  a  piecewise  polynomial  function
C   represented as a tensor product of one-dimensional  B-splines.  The
C   form of this function is
C
C                          NX   NY
C              S(x,y)  =  SUM  SUM  a   U (x) V (y)
C                         i=1  j=1   ij  i     j
C
C   where the functions U(i)  and  V(j)  are  one-dimensional  B-spline
C   basis functions. The coefficients a(i,j) are chosen so that
C
C         S(X(i),Y(j)) = FCN(i,j)   for i=1,..,NX and j=1,..,NY
C
C   Note that  for  each  fixed  value  of  y  S(x,y)  is  a  piecewise
C   polynomial function of x alone, and for each fixed value of x  S(x,
C   y) is a piecewise polynomial function of y alone. In one  dimension
C   a piecewise polynomial may  be  created  by  partitioning  a  given
C   interval into subintervals and defining a distinct polynomial piece
C   on each one. The points where adjacent subintervals meet are called
C   knots. Each of the functions U(i) and V(j)  above  is  a  piecewise
C   polynomial.
C
C   Users of DB2INK choose  the  order  (degree+1)  of  the  polynomial
C   pieces used to define the piecewise polynomial in each of the x and
C   y directions (KX and KY). Users also  may  define  their  own  knot
C   sequence in x and y separately (TX and TY).  If  IFLAG=0,  however,
C   DB2INK will choose sequences of knots that result  in  a  piecewise
C   polynomial interpolant with KX-2 continuous partial derivatives  in
C   x and KY-2 continuous partial derivatives in y. (KX knots are taken
C   near each endpoint in the x direction,  not-a-knot  end  conditions
C   are used, and the remaining knots are placed at data points  if  KX
C   is even or at midpoints between data points if KX  is  odd.  The  y
C   direction is treated similarly.)
C
C   After a call to DB2INK, all information  necessary  to  define  the
C   interpolating function are contained in the parameters NX, NY,  KX,
C   KY, TX, TY, and BCOEF. These quantities should not be altered until
C   after the last call of the evaluation routine DB2VAL.
C
C
C   I N P U T
C   ---------
C
C   X       Double precision 1D array (size NX)
C           Array of x abcissae. Must be strictly increasing.
C
C   NX      Integer scalar (.GE. 3)
C           Number of x abcissae.
C
C   Y       Double precision 1D array (size NY)
C           Array of y abcissae. Must be strictly increasing.
C
C   NY      Integer scalar (.GE. 3)
C           Number of y abcissae.
C
C   FCN     Double precision 2D array (size LDF by NY)
C           Array of function values to interpolate. FCN(I,J) should
C           contain the function value at the point (X(I),Y(J))
C
C   LDF     Integer scalar (.GE. NX)
C           The actual leading dimension of FCN used in the calling
C           calling program.
C
C   KX      Integer scalar (.GE. 2, .LT. NX)
C           The order of spline pieces in x.
C           (Order = polynomial degree + 1)
C
C   KY      Integer scalar (.GE. 2, .LT. NY)
C           The order of spline pieces in y.
C           (Order = polynomial degree + 1)
C
C
C   I N P U T   O R   O U T P U T
C   -----------------------------
C
C   TX      Double precision 1D array (size NX+KX)
C           The knots in the x direction for the spline interpolant.
C           If IFLAG=0 these are chosen by DB2INK.
C           If IFLAG=1 these are specified by the user.
C                      (Must be non-decreasing.)
C
C   TY      Double precision 1D array (size NY+KY)
C           The knots in the y direction for the spline interpolant.
C           If IFLAG=0 these are chosen by DB2INK.
C           If IFLAG=1 these are specified by the user.
C                      (Must be non-decreasing.)
C
C
C   O U T P U T
C   -----------
C
C   BCOEF   Double precision 2D array (size NX by NY)
C           Array of coefficients of the B-spline interpolant.
C           This may be the same array as FCN.
C
C
C   M I S C E L L A N E O U S
C   -------------------------
C
C   WORK    Double precision 1D array (size NX*NY + max( 2*KX*(NX+1),
C                                             2*KY*(NY+1) ))
C           Array of working storage.
C
C   IFLAG   Integer scalar.
C           On input:  0 == knot sequence chosen by DB2INK
C                      1 == knot sequence chosen by user.
C           On output: 1 == successful execution
C                      2 == IFLAG out of range
C                      3 == NX out of range
C                      4 == KX out of range
C                      5 == X not strictly increasing
C                      6 == TX not non-decreasing
C                      7 == NY out of range
C                      8 == KY out of range
C                      9 == Y not strictly increasing
C                     10 == TY not non-decreasing
C
C***REFERENCES  CARL DE BOOR, A PRACTICAL GUIDE TO SPLINES,
C                 SPRINGER-VERLAG, NEW YORK, 1978.
C               CARL DE BOOR, EFFICIENT COMPUTER MANIPULATION OF TENSOR
C                 PRODUCTS, ACM TRANSACTIONS ON MATHEMATICAL SOFTWARE,
C                 VOL. 5 (1979), PP. 173-182.
C***ROUTINES CALLED  DBTPCF,DBKNOT
C***END PROLOGUE  DB2INK
C
C  ------------
C  DECLARATIONS
C  ------------
C
C  PARAMETERS
C
      INTEGER
     *        NX, NY, LDF, KX, KY, IFLAG
      DOUBLE PRECISION
     *     X(NX), Y(NY), FCN(LDF,NY), TX(*), TY(*), BCOEF(NX,NY),
     *     WORK(*)
C
C  LOCAL VARIABLES
C
      INTEGER
     *        I, IW, NPK
C
C  -----------------------
C  CHECK VALIDITY OF INPUT
C  -----------------------
C
C***FIRST EXECUTABLE STATEMENT
      IF ((IFLAG .LT. 0) .OR. (IFLAG .GT. 1))  GO TO 920
      IF (NX .LT. 3)  GO TO 930
      IF (NY .LT. 3)  GO TO 970
      IF ((KX .LT. 2) .OR. (KX .GE. NX))  GO TO 940
      IF ((KY .LT. 2) .OR. (KY .GE. NY))  GO TO 980
      DO 10 I=2,NX
         IF (X(I) .LE. X(I-1))  GO TO 950
   10 CONTINUE
      DO 20 I=2,NY
         IF (Y(I) .LE. Y(I-1))  GO TO 990
   20 CONTINUE
      IF (IFLAG .EQ. 0)  GO TO 50
         NPK = NX + KX
         DO 30 I=2,NPK
            IF (TX(I) .LT. TX(I-1))  GO TO 960
   30    CONTINUE
         NPK = NY + KY
         DO 40 I=2,NPK
            IF (TY(I) .LT. TY(I-1))  GO TO 1000
   40    CONTINUE
   50 CONTINUE
C
C  ------------
C  CHOOSE KNOTS
C  ------------
C
      IF (IFLAG .NE. 0)  GO TO 100
         CALL DBKNOT(X,NX,KX,TX)
         CALL DBKNOT(Y,NY,KY,TY)
  100 CONTINUE
C
C  -------------------------------
C  CONSTRUCT B-SPLINE COEFFICIENTS
C  -------------------------------
C
      IFLAG = 1
      IW = NX*NY + 1
      CALL DBTPCF(X,NX,FCN,LDF,NY,TX,KX,WORK,WORK(IW))
      CALL DBTPCF(Y,NY,WORK,NY,NX,TY,KY,BCOEF,WORK(IW))
      GO TO 9999
C
C  -----
C  EXITS
C  -----
C
  920 CONTINUE
      CALL XERRWV('DB2INK -  IFLAG=I1 IS OUT OF RANGE.',
     *            36,2,1,1,IFLAG,I2,0,R1,R2)
      IFLAG = 2
      GO TO 9999
C
  930 CONTINUE
      IFLAG = 3
      CALL XERRWV('DB2INK -  NX=I1 IS OUT OF RANGE.',
     *            32,IFLAG,1,1,NX,I2,0,R1,R2)
      GO TO 9999
C
  940 CONTINUE
      IFLAG = 4
      CALL XERRWV('DB2INK -  KX=I1 IS OUT OF RANGE.',
     *            32,IFLAG,1,1,KX,I2,0,R1,R2)
      GO TO 9999
C
  950 CONTINUE
      IFLAG = 5
      CALL XERRWV('DB2INK -  X ARRAY MUST BE STRICTLY INCREASING.',
     *            46,IFLAG,1,0,I1,I2,0,R1,R2)
      GO TO 9999
C
  960 CONTINUE
      IFLAG = 6
      CALL XERRWV('DB2INK -  TX ARRAY MUST BE NON-DECREASING.',
     *            42,IFLAG,1,0,I1,I2,0,R1,R2)
      GO TO 9999
C
  970 CONTINUE
      IFLAG = 7
      CALL XERRWV('DB2INK -  NY=I1 IS OUT OF RANGE.',
     *            32,IFLAG,1,1,NY,I2,0,R1,R2)
      GO TO 9999
C
  980 CONTINUE
      IFLAG = 8
      CALL XERRWV('DB2INK -  KY=I1 IS OUT OF RANGE.',
     *            32,IFLAG,1,1,KY,I2,0,R1,R2)
      GO TO 9999
C
  990 CONTINUE
      IFLAG = 9
      CALL XERRWV('DB2INK -  Y ARRAY MUST BE STRICTLY INCREASING.',
     *            46,IFLAG,1,0,I1,I2,0,R1,R2)
      GO TO 9999
C
 1000 CONTINUE
      IFLAG = 10
      CALL XERRWV('DB2INK -  TY ARRAY MUST BE NON-DECREASING.',
     *            42,IFLAG,1,0,I1,I2,0,R1,R2)
      GO TO 9999
C
 9999 CONTINUE
      RETURN
      END

      DOUBLE PRECISION FUNCTION DB2VAL(XVAL,YVAL,IDX,IDY,TX,TY,NX,NY,
     *  KX,KY,BCOEF,WORK)
C***BEGIN PROLOGUE  DB2VAL
C***DATE WRITTEN   25 MAY 1982
C***REVISION DATE  25 MAY 1982
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***CATEGORY NO.  E1A
C***KEYWORDS  INTERPOLATION, TWO-DIMENSIONS, GRIDDED DATA, SPLINES,
C             PIECEWISE POLYNOMIALS
C***AUTHOR  BOISVERT, RONALD, NBS
C             SCIENTIFIC COMPUTING DIVISION
C             NATIONAL BUREAU OF STANDARDS
C             WASHINGTON, DC 20234
C***PURPOSE  DB2VAL EVALUATES THE PIECEWISE POLYNOMIAL INTERPOLATING
C            FUNCTION CONSTRUCTED BY THE ROUTINE DB2INK OR ONE OF ITS
C            PARTIAL DERIVATIVES.
C            DOUBLE PRECISION VERSION OF B2VAL.
C***DESCRIPTION
C
C   DB2VAL  evaluates   the   tensor   product   piecewise   polynomial
C   interpolant constructed  by  the  routine  DB2INK  or  one  of  its
C   derivatives at the point (XVAL,YVAL). To evaluate  the  interpolant
C   itself, set IDX=IDY=0, to evaluate the first partial  with  respect
C   to x, set IDX=1,IDY=0, and so on.
C
C   DB2VAL returns 0.0E0 if (XVAL,YVAL) is out of range. That is, if
C            XVAL.LT.TX(1) .OR. XVAL.GT.TX(NX+KX) .OR.
C            YVAL.LT.TY(1) .OR. YVAL.GT.TY(NY+KY)
C   If the knots TX  and  TY  were  chosen  by  DB2INK,  then  this  is
C   equivalent to
C            XVAL.LT.X(1) .OR. XVAL.GT.X(NX)+EPSX .OR.
C            YVAL.LT.Y(1) .OR. YVAL.GT.Y(NY)+EPSY
C   where EPSX = 0.1*(X(NX)-X(NX-1)) and EPSY = 0.1*(Y(NY)-Y(NY-1)).
C
C   The input quantities TX, TY, NX, NY, KX, KY, and  BCOEF  should  be
C   unchanged since the last call of DB2INK.
C
C
C   I N P U T
C   ---------
C
C   XVAL    Double precision scalar
C           X coordinate of evaluation point.
C
C   YVAL    Double precision scalar
C           Y coordinate of evaluation point.
C
C   IDX     Integer scalar
C           X derivative of piecewise polynomial to evaluate.
C
C   IDY     Integer scalar
C           Y derivative of piecewise polynomial to evaluate.
C
C   TX      Double precision 1D array (size NX+KX)
C           Sequence of knots defining the piecewise polynomial in
C           the x direction.  (Same as in last call to DB2INK.)
C
C   TY      Double precision 1D array (size NY+KY)
C           Sequence of knots defining the piecewise polynomial in
C           the y direction.  (Same as in last call to DB2INK.)
C
C   NX      Integer scalar
C           The number of interpolation points in x.
C           (Same as in last call to DB2INK.)
C
C   NY      Integer scalar
C           The number of interpolation points in y.
C           (Same as in last call to DB2INK.)
C
C   KX      Integer scalar
C           Order of polynomial pieces in x.
C           (Same as in last call to DB2INK.)
C
C   KY      Integer scalar
C           Order of polynomial pieces in y.
C           (Same as in last call to DB2INK.)
C
C   BCOEF   Double precision 2D array (size NX by NY)
C           The B-spline coefficients computed by DB2INK.
C
C   WORK    Double precision 1D array (size 3*max(KX,KY) + KY)
C           A working storage array.
C
C***REFERENCES  CARL DE BOOR, A PRACTICAL GUIDE TO SPLINES,
C                 SPRINGER-VERLAG, NEW YORK, 1978.
C***ROUTINES CALLED  DINTRV,DBVALU
C***END PROLOGUE  DB2VAL
C
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C
C   MODIFICATION
C   ------------
C
C   ADDED CHECK TO SEE IF X OR Y IS OUT OF RANGE, IF SO, RETURN 0.0
C
C   R.F. BOISVERT, NIST
C   22 FEB 00
C
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C  ------------
C  DECLARATIONS
C  ------------
C
C  PARAMETERS
C
      INTEGER
     *        IDX, IDY, NX, NY, KX, KY
      DOUBLE PRECISION
     *     XVAL, YVAL, TX(*), TY(*), BCOEF(NX,NY), WORK(*)
C
C  LOCAL VARIABLES
C
      INTEGER
     *        ILOY, INBVX, INBV, K, LEFTY, MFLAG, KCOL, IW
      DOUBLE PRECISION
     *     DBVALU
C
      DATA ILOY /1/,  INBVX /1/
C     SAVE ILOY    ,  INBVX
C
C
C***FIRST EXECUTABLE STATEMENT
      DB2VAL = 0.0D0
C  NEXT STATEMENT - RFB MOD
      IF (XVAL.LT.TX(1) .OR. XVAL.GT.TX(NX+KX) .OR.
     +    YVAL.LT.TY(1) .OR. YVAL.GT.TY(NY+KY))      GO TO 100
      CALL DINTRV(TY,NY+KY,YVAL,ILOY,LEFTY,MFLAG)
      IF (MFLAG .NE. 0)  GO TO 100
         IW = KY + 1
         KCOL = LEFTY - KY
         DO 50 K=1,KY
            KCOL = KCOL + 1
            WORK(K) = DBVALU(TX,BCOEF(1,KCOL),NX,KX,IDX,XVAL,INBVX,
     *                      WORK(IW))
   50    CONTINUE
         INBV = 1
         KCOL = LEFTY - KY + 1
         DB2VAL = DBVALU(TY(KCOL),WORK,KY,KY,IDY,YVAL,INBV,WORK(IW))
  100 CONTINUE
      RETURN
      END
      SUBROUTINE DB3INK(X,NX,Y,NY,Z,NZ,FCN,LDF1,LDF2,KX,KY,KZ,TX,TY,TZ,
     *  BCOEF,WORK,IFLAG)
C***BEGIN PROLOGUE  DB3INK
C***DATE WRITTEN   25 MAY 1982
C***REVISION DATE  25 MAY 1982
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***CATEGORY NO.  E1A
C***KEYWORDS  INTERPOLATION, THREE-DIMENSIONS, GRIDDED DATA, SPLINES,
C             PIECEWISE POLYNOMIALS
C***AUTHOR  BOISVERT, RONALD, NBS
C             SCIENTIFIC COMPUTING DIVISION
C             NATIONAL BUREAU OF STANDARDS
C             WASHINGTON, DC 20234
C***PURPOSE  DOUBLE PRECISION VERSION OF DB3INK
C            DB3INK DETERMINES A PIECEWISE POLYNOMIAL FUNCTION THAT
C            INTERPOLATES THREE-DIMENSIONAL GRIDDED DATA. USERS SPECIFY
C            THE POLYNOMIAL ORDER (DEGREE+1) OF THE INTERPOLANT AND
C            (OPTIONALLY) THE KNOT SEQUENCE.
C***DESCRIPTION
C
C   DB3INK determines the parameters of a  function  that  interpolates
C   the three-dimensional gridded data (X(i),Y(j),Z(k),FCN(i,j,k))  for
C   i=1,..,NX, j=1,..,NY, and k=1,..,NZ. The interpolating function and
C   its derivatives may  subsequently  be  evaluated  by  the  function
C   DB3VAL.
C
C   The interpolating  function  is  a  piecewise  polynomial  function
C   represented as a tensor product of one-dimensional  B-splines.  The
C   form of this function is
C
C                      NX   NY   NZ
C        S(x,y,z)  =  SUM  SUM  SUM  a   U (x) V (y) W (z)
C                     i=1  j=1  k=1   ij  i     j     k
C
C   where the functions U(i), V(j), and  W(k)  are  one-dimensional  B-
C   spline basis functions. The coefficients a(i,j) are chosen so that
C
C   S(X(i),Y(j),Z(k)) = FCN(i,j,k)  for i=1,..,NX, j=1,..,NY, k=1,..,NZ
C
C   Note that for fixed values of y  and  z  S(x,y,z)  is  a  piecewise
C   polynomial function of x alone, for fixed values of x and z  S(x,y,
C   z) is a piecewise polynomial function of y  alone,  and  for  fixed
C   values of x and y S(x,y,z)  is  a  function  of  z  alone.  In  one
C   dimension a piecewise polynomial may be created by  partitioning  a
C   given interval into subintervals and defining a distinct polynomial
C   piece on each one. The points where adjacent subintervals meet  are
C   called knots. Each of the functions U(i), V(j), and W(k) above is a
C   piecewise polynomial.
C
C   Users of DB3INK choose  the  order  (degree+1)  of  the  polynomial
C   pieces used to define the piecewise polynomial in each of the x, y,
C   and z directions (KX, KY, and KZ). Users also may define their  own
C   knot sequence in x, y, and z separately (TX, TY, and TZ). If IFLAG=
C   0, however, DB3INK will choose sequences of knots that result in  a
C   piecewise  polynomial  interpolant  with  KX-2  continuous  partial
C   derivatives in x, KY-2 continuous partial derivatives in y, and KZ-
C   2 continuous partial derivatives in z. (KX  knots  are  taken  near
C   each endpoint in x, not-a-knot end conditions  are  used,  and  the
C   remaining knots are placed at data points  if  KX  is  even  or  at
C   midpoints between data points if KX is odd. The y and z  directions
C   are treated similarly.)
C
C   After a call to DB3INK, all information  necessary  to  define  the
C   interpolating function are contained in the parameters NX, NY,  NZ,
C   KX, KY, KZ, TX, TY, TZ, and BCOEF. These quantities should  not  be
C   altered until after the last call of the evaluation routine DB3VAL.
C
C
C   I N P U T
C   ---------
C
C   X       Double precision 1D array (size NX)
C           Array of x abcissae. Must be strictly increasing.
C
C   NX      Integer scalar (.GE. 3)
C           Number of x abcissae.
C
C   Y       Double precision 1D array (size NY)
C           Array of y abcissae. Must be strictly increasing.
C
C   NY      Integer scalar (.GE. 3)
C           Number of y abcissae.
C
C   Z       Double precision 1D array (size NZ)
C           Array of z abcissae. Must be strictly increasing.
C
C   NZ      Integer scalar (.GE. 3)
C           Number of z abcissae.
C
C   FCN     Double precision 3D array (size LDF1 by LDF2 by NY)
C           Array of function values to interpolate. FCN(I,J,K) should
C           contain the function value at the point (X(I),Y(J),Z(K))
C
C   LDF1    Integer scalar (.GE. NX)
C           The actual first dimension of FCN used in the
C           calling program.
C
C   LDF2    Integer scalar (.GE. NY)
C           The actual second dimension of FCN used in the calling
C           program.
C
C   KX      Integer scalar (.GE. 2, .LT. NX)
C           The order of spline pieces in x.
C           (Order = polynomial degree + 1)
C
C   KY      Integer scalar (.GE. 2, .LT. NY)
C           The order of spline pieces in y.
C           (Order = polynomial degree + 1)
C
C   KZ      Integer scalar (.GE. 2, .LT. NZ)
C           The order of spline pieces in z.
C           (Order = polynomial degree + 1)
C
C
C   I N P U T   O R   O U T P U T
C   -----------------------------
C
C   TX      Double precision 1D array (size NX+KX)
C           The knots in the x direction for the spline interpolant.
C           If IFLAG=0 these are chosen by DB3INK.
C           If IFLAG=1 these are specified by the user.
C                      (Must be non-decreasing.)
C
C   TY      Double precision 1D array (size NY+KY)
C           The knots in the y direction for the spline interpolant.
C           If IFLAG=0 these are chosen by DB3INK.
C           If IFLAG=1 these are specified by the user.
C                      (Must be non-decreasing.)
C
C   TZ      Double precision 1D array (size NZ+KZ)
C           The knots in the z direction for the spline interpolant.
C           If IFLAG=0 these are chosen by DB3INK.
C           If IFLAG=1 these are specified by the user.
C                      (Must be non-decreasing.)
C
C
C   O U T P U T
C   -----------
C
C   BCOEF   Double precision 3D array (size NX by NY by NZ)
C           Array of coefficients of the B-spline interpolant.
C           This may be the same array as FCN.
C
C
C   M I S C E L L A N E O U S
C   -------------------------
C
C   WORK    Double precision 1D array (size NX*NY*NZ + max( 2*KX*(NX+1),
C                             2*KY*(NY+1), 2*KZ*(NZ+1) )
C           Array of working storage.
C
C   IFLAG   Integer scalar.
C           On input:  0 == knot sequence chosen by B2INK
C                      1 == knot sequence chosen by user.
C           On output: 1 == successful execution
C                      2 == IFLAG out of range
C                      3 == NX out of range
C                      4 == KX out of range
C                      5 == X not strictly increasing
C                      6 == TX not non-decreasing
C                      7 == NY out of range
C                      8 == KY out of range
C                      9 == Y not strictly increasing
C                     10 == TY not non-decreasing
C                     11 == NZ out of range
C                     12 == KZ out of range
C                     13 == Z not strictly increasing
C                     14 == TY not non-decreasing
C
C***REFERENCES  CARL DE BOOR, A PRACTICAL GUIDE TO SPLINES,
C                 SPRINGER-VERLAG, NEW YORK, 1978.
C               CARL DE BOOR, EFFICIENT COMPUTER MANIPULATION OF TENSOR
C                 PRODUCTS, ACM TRANSACTIONS ON MATHEMATICAL SOFTWARE,
C                 VOL. 5 (1979), PP. 173-182.
C***ROUTINES CALLED  DBTPCF,DBKNOT
C***END PROLOGUE  DB3INK
C
C  ------------
C  DECLARATIONS
C  ------------
C
C  PARAMETERS
C
      INTEGER
     *        NX, NY, NZ, LDF1, LDF2, KX, KY, KZ, IFLAG
      DOUBLE PRECISION
     *     X(NX), Y(NY), Z(NZ), FCN(LDF1,LDF2,NZ), TX(*), TY(*), TZ(*),
     *     BCOEF(NX,NY,NZ), WORK(*)
C
C  LOCAL VARIABLES
C
      INTEGER
     *        I, J, LOC, IW, NPK
C
C  -----------------------
C  CHECK VALIDITY OF INPUT
C  -----------------------
C
C***FIRST EXECUTABLE STATEMENT
      IF ((IFLAG .LT. 0) .OR. (IFLAG .GT. 1))  GO TO 920
      IF (NX .LT. 3)  GO TO 930
      IF (NY .LT. 3)  GO TO 970
      IF (NZ .LT. 3)  GO TO 1010
      IF ((KX .LT. 2) .OR. (KX .GE. NX))  GO TO 940
      IF ((KY .LT. 2) .OR. (KY .GE. NY))  GO TO 980
      IF ((KZ .LT. 2) .OR. (KZ .GE. NZ))  GO TO 1020
      DO 10 I=2,NX
         IF (X(I) .LE. X(I-1))  GO TO 950
   10 CONTINUE
      DO 20 I=2,NY
         IF (Y(I) .LE. Y(I-1))  GO TO 990
   20 CONTINUE
      DO 30 I=2,NZ
         IF (Z(I) .LE. Z(I-1))  GO TO 1030
   30 CONTINUE
      IF (IFLAG .EQ. 0)  GO TO 70
         NPK = NX + KX
         DO 40 I=2,NPK
            IF (TX(I) .LT. TX(I-1))  GO TO 960
   40    CONTINUE
         NPK = NY + KY
         DO 50 I=2,NPK
            IF (TY(I) .LT. TY(I-1))  GO TO 1000
   50    CONTINUE
         NPK = NZ + KZ
         DO 60 I=2,NPK
            IF (TZ(I) .LT. TZ(I-1))  GO TO 1040
   60    CONTINUE
   70 CONTINUE
C
C  ------------
C  CHOOSE KNOTS
C  ------------
C
      
      IF (IFLAG .NE. 0)  GO TO 100
         CALL DBKNOT(X,NX,KX,TX)
         CALL DBKNOT(Y,NY,KY,TY)
         CALL DBKNOT(Z,NZ,KZ,TZ)
  100 CONTINUE
C
C  -------------------------------
C  CONSTRUCT B-SPLINE COEFFICIENTS
C  -------------------------------
C
      IFLAG = 1
      IW = NX*NY*NZ + 1
C
C     COPY FCN TO WORK IN PACKED FOR DBTPCF
C      LOC = 0
C      DO 510 K=1,NZ
C         DO 510 J=1,NY
C            DO 510 I=1,NX
C               LOC = LOC + 1
C               WORK(LOC) = FCN(I,J,K)
C  510 CONTINUE
C
C      CALL DBTPCF(X,NX,WORK,NX,NY*NZ,TX,KX,BCOEF,WORK(IW))
C      CALL DBTPCF(Y,NY,BCOEF,NY,NX*NZ,TY,KY,WORK,WORK(IW))
C      CALL DBTPCF(Z,NZ,WORK,NZ,NX*NY,TZ,KZ,BCOEF,WORK(IW))
      
C  KBW CHANGE ===> eliminate unnecessary copy of FCN to WORK
      CALL DBTPCF(X,NX,FCN,NX,NY*NZ,TX,KX,BCOEF,WORK(IW))
      CALL DBTPCF(Y,NY,BCOEF,NY,NX*NZ,TY,KY,WORK,WORK(IW))
      CALL DBTPCF(Z,NZ,WORK,NZ,NX*NY,TZ,KZ,BCOEF,WORK(IW))
      
      
      GO TO 9999
C
C  -----
C  EXITS
C  -----
C
  920 CONTINUE
      CALL XERRWV('DB3INK -  IFLAG=I1 IS OUT OF RANGE.',
     *            35,2,1,1,IFLAG,I2,0,R1,R2)
      IFLAG = 2
      GO TO 9999
C
  930 CONTINUE
      IFLAG = 3
      CALL XERRWV('DB3INK -  NX=I1 IS OUT OF RANGE.',
     *            32,IFLAG,1,1,NX,I2,0,R1,R2)
      GO TO 9999
C
  940 CONTINUE
      IFLAG = 4
      CALL XERRWV('DB3INK -  KX=I1 IS OUT OF RANGE.',
     *            32,IFLAG,1,1,KX,I2,0,R1,R2)
      GO TO 9999
C
  950 CONTINUE
      IFLAG = 5
      CALL XERRWV('DB3INK -  X ARRAY MUST BE STRICTLY INCREASING.',
     *            46,IFLAG,1,0,I1,I2,0,R1,R2)
      GO TO 9999
C
  960 CONTINUE
      IFLAG = 6
      CALL XERRWV('DB3INK -  TX ARRAY MUST BE NON-DECREASING.',
     *            42,IFLAG,1,0,I1,I2,0,R1,R2)
      GO TO 9999
C
  970 CONTINUE
      IFLAG = 7
      CALL XERRWV('DB3INK -  NY=I1 IS OUT OF RANGE.',
     *            32,IFLAG,1,1,NY,I2,0,R1,R2)
      GO TO 9999
C
  980 CONTINUE
      IFLAG = 8
      CALL XERRWV('DB3INK -  KY=I1 IS OUT OF RANGE.',
     *            32,IFLAG,1,1,KY,I2,0,R1,R2)
      GO TO 9999
C
  990 CONTINUE
      IFLAG = 9
      CALL XERRWV('DB3INK -  Y ARRAY MUST BE STRICTLY INCREASING.',
     *            46,IFLAG,1,0,I1,I2,0,R1,R2)
      GO TO 9999
C
 1000 CONTINUE
      IFLAG = 10
      CALL XERRWV('DB3INK -  TY ARRAY MUST BE NON-DECREASING.',
     *            42,IFLAG,1,0,I1,I2,0,R1,R2)
      GO TO 9999
C
 1010 CONTINUE
      IFLAG = 11
      CALL XERRWV('DB3INK -  NZ=I1 IS OUT OF RANGE.',
     *            32,IFLAG,1,1,NZ,I2,0,R1,R2)
      GO TO 9999
C
 1020 CONTINUE
      IFLAG = 12
      write(0,*) '*******************',myid, kz, nz
      CALL XERRWV('DB3INK -  KZ=I1 IS OUT OF RANGE.',
     *            32,IFLAG,1,1,KZ,I2,0,R1,R2)
      GO TO 9999
C
 1030 CONTINUE
      IFLAG = 13
      CALL XERRWV('DB3INK -  Z ARRAY MUST BE STRICTLY INCREASING.',
     *            46,IFLAG,1,0,I1,I2,0,R1,R2)
      GO TO 9999
C
 1040 CONTINUE
      IFLAG = 14
      CALL XERRWV('DB3INK -  TZ ARRAY MUST BE NON-DECREASING.',
     *            42,IFLAG,1,0,I1,I2,0,R1,R2)
      GO TO 9999
C
 9999 CONTINUE
      RETURN
      END
      SUBROUTINE DBKNOT(X,N,K,T)
C***BEGIN PROLOGUE  DBKNOT
C***REFER TO  DB2INK,DB3INK
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***END PROLOGUE  DBKNOT
C
C  --------------------------------------------------------------------
C  DBKNOT CHOOSES A KNOT SEQUENCE FOR INTERPOLATION OF ORDER K AT THE
C  DATA POINTS X(I), I=1,..,N.  THE N+K KNOTS ARE PLACED IN THE ARRAY
C  T.  K KNOTS ARE PLACED AT EACH ENDPOINT AND NOT-A-KNOT END
C  CONDITIONS ARE USED.  THE REMAINING KNOTS ARE PLACED AT DATA POINTS
C  IF N IS EVEN AND BETWEEN DATA POINTS IF N IS ODD.  THE RIGHTMOST
C  KNOT IS SHIFTED SLIGHTLY TO THE RIGHT TO INSURE PROPER INTERPOLATION
C  AT X(N) (SEE PAGE 350 OF THE REFERENCE).
C  DOUBLE PRECISION VERSION OF BKNOT.
C  --------------------------------------------------------------------
C
C  ------------
C  DECLARATIONS
C  ------------
C
C  PARAMETERS
C
      INTEGER
     *        N, K
      DOUBLE PRECISION
     *     X(N), T(*)
C
C  LOCAL VARIABLES
C
      INTEGER
     *        I, J, IPJ, NPJ, IP1
      DOUBLE PRECISION
     *     RNOT
C
C
C  ----------------------------
C  PUT K KNOTS AT EACH ENDPOINT
C  ----------------------------
C
C     (SHIFT RIGHT ENPOINTS SLIGHTLY -- SEE PG 350 OF REFERENCE)
      RNOT = X(N) + 0.10D0*( X(N)-X(N-1) )
      DO 110 J=1,K
         T(J) = X(1)
         NPJ = N + J
         T(NPJ) = RNOT
  110 CONTINUE
C
C  --------------------------
C  DISTRIBUTE REMAINING KNOTS
C  --------------------------
C
      IF (MOD(K,2) .EQ. 1)  GO TO 150
C
C     CASE OF EVEN K --  KNOTS AT DATA POINTS
C
      I = (K/2) - K
      JSTRT = K+1
      DO 120 J=JSTRT,N
         IPJ = I + J
         T(J) = X(IPJ)
  120 CONTINUE
      GO TO 200
C
C     CASE OF ODD K --  KNOTS BETWEEN DATA POINTS
C
  150 CONTINUE
      I = (K-1)/2 - K
      IP1 = I + 1
      JSTRT = K + 1
      DO 160 J=JSTRT,N
         IPJ = I + J
         T(J) = 0.50D0*( X(IPJ) + X(IPJ+1) )
  160 CONTINUE
  200 CONTINUE
C
      RETURN
      END
      SUBROUTINE DBTPCF(X,N,FCN,LDF,NF,T,K,BCOEF,WORK)
C***BEGIN PROLOGUE  DBTPCF
C***REFER TO  DB2INK,DB3INK
C***ROUTINES CALLED  DBINTK,DBNSLV
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***END PROLOGUE  DBTPCF
C
C  -----------------------------------------------------------------
C  DBTPCF COMPUTES B-SPLINE INTERPOLATION COEFFICIENTS FOR NF SETS
C  OF DATA STORED IN THE COLUMNS OF THE ARRAY FCN. THE B-SPLINE
C  COEFFICIENTS ARE STORED IN THE ROWS OF BCOEF HOWEVER.
C  EACH INTERPOLATION IS BASED ON THE N ABCISSA STORED IN THE
C  ARRAY X, AND THE N+K KNOTS STORED IN THE ARRAY T. THE ORDER
C  OF EACH INTERPOLATION IS K. THE WORK ARRAY MUST BE OF LENGTH
C  AT LEAST 2*K*(N+1).
C  DOUBLE PRECISION VERSION OF BTPCF.
C  -----------------------------------------------------------------
C
C  ------------
C  DECLARATIONS
C  ------------
C
C  PARAMETERS
C
      INTEGER
     *        N, LDF, K
      DOUBLE PRECISION
     *     X(N), FCN(LDF,NF), T(*), BCOEF(NF,N), WORK(*)
C
C  LOCAL VARIABLES
C
      INTEGER
     *        I, J, K1, K2, IQ, IW
C
C  ---------------------------------------------
C  CHECK FOR NULL INPUT AND PARTITION WORK ARRAY
C  ---------------------------------------------
C
C***FIRST EXECUTABLE STATEMENT
      IF (NF .LE. 0)  GO TO 500
      K1 = K - 1
      K2 = K1 + K
      IQ = 1 + N
      IW = IQ + K2*N+1
C
C  -----------------------------
C  COMPUTE B-SPLINE COEFFICIENTS
C  -----------------------------
C
C
C   FIRST DATA SET
C
      CALL DBINTK(X,FCN,T,N,K,WORK,WORK(IQ),WORK(IW))
      DO 20 I=1,N
         BCOEF(1,I) = WORK(I)
   20 CONTINUE
C
C  ALL REMAINING DATA SETS BY BACK-SUBSTITUTION
C
      IF (NF .EQ. 1)  GO TO 500
      DO 100 J=2,NF
         DO 50 I=1,N
            WORK(I) = FCN(I,J)
   50    CONTINUE
         CALL DBNSLV(WORK(IQ),K2,N,K1,K1,WORK)
         DO 60 I=1,N
            BCOEF(J,I) = WORK(I)
   60    CONTINUE
  100 CONTINUE
C
C  ----
C  EXIT
C  ----
C
  500 CONTINUE
      RETURN
      END
      SUBROUTINE XERRWV(MESSG,NMESSG,NERR,LEVEL,NI,I1,I2,NR,R1,R2)
C***BEGIN PROLOGUE  XERRWV
C***DATE WRITTEN   800319   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Processes error message allowing 2 integer and two real
C            values to be included in the message.
C***DESCRIPTION
C     Abstract
C        XERRWV processes a diagnostic message, in a manner
C        determined by the value of LEVEL and the current value
C        of the library error control flag, KONTRL.
C        (See subroutine XSETF for details.)
C        In addition, up to two integer values and two real
C        values may be printed along with the message.
C
C     Description of Parameters
C      --Input--
C        MESSG - the Hollerith message to be processed.
C        NMESSG- the actual number of characters in MESSG.
C        NERR  - the error number associated with this message.
C                NERR must not be zero.
C        LEVEL - error category.
C                =2 means this is an unconditionally fatal error.
C                =1 means this is a recoverable error.  (I.e., it is
C                   non-fatal if XSETF has been appropriately called.)
C                =0 means this is a warning message only.
C                =-1 means this is a warning message which is to be
C                   printed at most once, regardless of how many
C                   times this call is executed.
C        NI    - number of integer values to be printed. (0 to 2)
C        I1    - first integer value.
C        I2    - second integer value.
C        NR    - number of real values to be printed. (0 to 2)
C        R1    - first real value.
C        R2    - second real value.
C
C     Examples
C        CALL XERRWV('SMOOTH -- NUM (=I1) WAS ZERO.',29,1,2,
C    1   1,NUM,0,0,0.,0.)
C        CALL XERRWV('QUADXY -- REQUESTED ERROR (R1) LESS THAN MINIMUM (
C    1R2).,54,77,1,0,0,0,2,ERRREQ,ERRMIN')
C
C     Latest revision ---  19 MAR 1980
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  FDUMP,I1MACH,J4SAVE,XERABT,XERCTL,XERPRT,XERSAV,
C                    XGETUA
C***END PROLOGUE  XERRWV
      CHARACTER*(*) MESSG
      CHARACTER*20 LFIRST
      CHARACTER*37 FORM
      DIMENSION LUN(5)
C     GET FLAGS
C***FIRST EXECUTABLE STATEMENT  XERRWV
      LKNTRL = J4SAVE(2,0,.FALSE.)
      MAXMES = J4SAVE(4,0,.FALSE.)
C     CHECK FOR VALID INPUT
      IF ((NMESSG.GT.0).AND.(NERR.NE.0).AND.
     1    (LEVEL.GE.(-1)).AND.(LEVEL.LE.2)) GO TO 10
         IF (LKNTRL.GT.0) CALL XERPRT('FATAL ERROR IN...',17)
         CALL XERPRT('XERROR -- INVALID INPUT',23)
         IF (LKNTRL.GT.0) CALL FDUMP
         IF (LKNTRL.GT.0) CALL XERPRT('JOB ABORT DUE TO FATAL ERROR.',
     1  29)
         IF (LKNTRL.GT.0) CALL XERSAV(' ',0,0,0,KDUMMY)
         CALL XERABT('XERROR -- INVALID INPUT',23)
         RETURN
   10 CONTINUE
C     RECORD MESSAGE
      JUNK = J4SAVE(1,NERR,.TRUE.)
      CALL XERSAV(MESSG,NMESSG,NERR,LEVEL,KOUNT)
C     LET USER OVERRIDE
      LFIRST = MESSG
      LMESSG = NMESSG
      LERR = NERR
      LLEVEL = LEVEL
      CALL XERCTL(LFIRST,LMESSG,LERR,LLEVEL,LKNTRL)
C     RESET TO ORIGINAL VALUES
      LMESSG = NMESSG
      LERR = NERR
      LLEVEL = LEVEL
      LKNTRL = MAX0(-2,MIN0(2,LKNTRL))
      MKNTRL = IABS(LKNTRL)
C     DECIDE WHETHER TO PRINT MESSAGE
      IF ((LLEVEL.LT.2).AND.(LKNTRL.EQ.0)) GO TO 100
      IF (((LLEVEL.EQ.(-1)).AND.(KOUNT.GT.MIN0(1,MAXMES)))
     1.OR.((LLEVEL.EQ.0)   .AND.(KOUNT.GT.MAXMES))
     2.OR.((LLEVEL.EQ.1)   .AND.(KOUNT.GT.MAXMES).AND.(MKNTRL.EQ.1))
     3.OR.((LLEVEL.EQ.2)   .AND.(KOUNT.GT.MAX0(1,MAXMES)))) GO TO 100
         IF (LKNTRL.LE.0) GO TO 20
            CALL XERPRT(' ',1)
C           INTRODUCTION
            IF (LLEVEL.EQ.(-1)) CALL XERPRT
     1('WARNING MESSAGE...THIS MESSAGE WILL ONLY BE PRINTED ONCE.',57)
            IF (LLEVEL.EQ.0) CALL XERPRT('WARNING IN...',13)
            IF (LLEVEL.EQ.1) CALL XERPRT
     1      ('RECOVERABLE ERROR IN...',23)
            IF (LLEVEL.EQ.2) CALL XERPRT('FATAL ERROR IN...',17)
   20    CONTINUE
C        MESSAGE
         CALL XERPRT(MESSG,LMESSG)
         CALL XGETUA(LUN,NUNIT)
         ISIZEI = LOG10(FLOAT(I1MACH(9))) + 1.0
         ISIZEF = LOG10(FLOAT(I1MACH(10))**I1MACH(11)) + 1.0
         DO 50 KUNIT=1,NUNIT
            IUNIT = LUN(KUNIT)
            IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
            DO 22 I=1,MIN(NI,2)
               WRITE (FORM,21) I,ISIZEI
   21          FORMAT ('(11X,21HIN ABOVE MESSAGE, I',I1,'=,I',I2,')   ')
               IF (I.EQ.1) WRITE (IUNIT,FORM) I1
               IF (I.EQ.2) WRITE (IUNIT,FORM) I2
   22       CONTINUE
            DO 24 I=1,MIN(NR,2)
               WRITE (FORM,23) I,ISIZEF+10,ISIZEF
   23          FORMAT ('(11X,21HIN ABOVE MESSAGE, R',I1,'=,E',
     1         I2,'.',I2,')')
               IF (I.EQ.1) WRITE (IUNIT,FORM) R1
               IF (I.EQ.2) WRITE (IUNIT,FORM) R2
   24       CONTINUE
            IF (LKNTRL.LE.0) GO TO 40
C              ERROR NUMBER
               WRITE (IUNIT,30) LERR
   30          FORMAT (15H ERROR NUMBER =,I10)
   40       CONTINUE
   50    CONTINUE
C        TRACE-BACK
         IF (LKNTRL.GT.0) CALL FDUMP
  100 CONTINUE
      IFATAL = 0
      IF ((LLEVEL.EQ.2).OR.((LLEVEL.EQ.1).AND.(MKNTRL.EQ.2)))
     1IFATAL = 1
C     QUIT HERE IF MESSAGE IS NOT FATAL
      IF (IFATAL.LE.0) RETURN
      IF ((LKNTRL.LE.0).OR.(KOUNT.GT.MAX0(1,MAXMES))) GO TO 120
C        PRINT REASON FOR ABORT
         IF (LLEVEL.EQ.1) CALL XERPRT
     1   ('JOB ABORT DUE TO UNRECOVERED ERROR.',35)
         IF (LLEVEL.EQ.2) CALL XERPRT
     1   ('JOB ABORT DUE TO FATAL ERROR.',29)
C        PRINT ERROR SUMMARY
         CALL XERSAV(' ',-1,0,0,KDUMMY)
  120 CONTINUE
C     ABORT
      IF ((LLEVEL.EQ.2).AND.(KOUNT.GT.MAX0(1,MAXMES))) LMESSG = 0
      CALL XERABT(MESSG,LMESSG)
      RETURN
      END
      SUBROUTINE XERSAV(MESSG,NMESSG,NERR,LEVEL,ICOUNT)
C***BEGIN PROLOGUE  XERSAV
C***DATE WRITTEN   800319   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  Z
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Records that an error occurred.
C***DESCRIPTION
C     Abstract
C        Record that this error occurred.
C
C     Description of Parameters
C     --Input--
C       MESSG, NMESSG, NERR, LEVEL are as in XERROR,
C       except that when NMESSG=0 the tables will be
C       dumped and cleared, and when NMESSG is less than zero the
C       tables will be dumped and not cleared.
C     --Output--
C       ICOUNT will be the number of times this message has
C       been seen, or zero if the table has overflowed and
C       does not contain this message specifically.
C       When NMESSG=0, ICOUNT will not be altered.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C     Latest revision ---  19 Mar 1980
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  I1MACH,S88FMT,XGETUA
C***END PROLOGUE  XERSAV
      INTEGER LUN(5)
      CHARACTER*(*) MESSG
      CHARACTER*20 MESTAB(10),MES
      DIMENSION NERTAB(10),LEVTAB(10),KOUNT(10)
      SAVE MESTAB,NERTAB,LEVTAB,KOUNT,KOUNTX
C     NEXT TWO DATA STATEMENTS ARE NECESSARY TO PROVIDE A BLANK
C     ERROR TABLE INITIALLY
      DATA KOUNT(1),KOUNT(2),KOUNT(3),KOUNT(4),KOUNT(5),
     1     KOUNT(6),KOUNT(7),KOUNT(8),KOUNT(9),KOUNT(10)
     2     /0,0,0,0,0,0,0,0,0,0/
      DATA KOUNTX/0/
C***FIRST EXECUTABLE STATEMENT  XERSAV
      IF (NMESSG.GT.0) GO TO 80
C     DUMP THE TABLE
         IF (KOUNT(1).EQ.0) RETURN
C        PRINT TO EACH UNIT
         CALL XGETUA(LUN,NUNIT)
         DO 60 KUNIT=1,NUNIT
            IUNIT = LUN(KUNIT)
            IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
C           PRINT TABLE HEADER
            WRITE (IUNIT,10)
   10       FORMAT (32H0          ERROR MESSAGE SUMMARY/
     1      51H MESSAGE START             NERR     LEVEL     COUNT)
C           PRINT BODY OF TABLE
            DO 20 I=1,10
               IF (KOUNT(I).EQ.0) GO TO 30
               WRITE (IUNIT,15) MESTAB(I),NERTAB(I),LEVTAB(I),KOUNT(I)
   15          FORMAT (1X,A20,3I10)
   20       CONTINUE
   30       CONTINUE
C           PRINT NUMBER OF OTHER ERRORS
            IF (KOUNTX.NE.0) WRITE (IUNIT,40) KOUNTX
   40       FORMAT (41H0OTHER ERRORS NOT INDIVIDUALLY TABULATED=,I10)
            WRITE (IUNIT,50)
   50       FORMAT (1X)
   60    CONTINUE
         IF (NMESSG.LT.0) RETURN
C        CLEAR THE ERROR TABLES
         DO 70 I=1,10
   70       KOUNT(I) = 0
         KOUNTX = 0
         RETURN
   80 CONTINUE
C     PROCESS A MESSAGE...
C     SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG,
C     OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL.
      MES = MESSG
      DO 90 I=1,10
         II = I
         IF (KOUNT(I).EQ.0) GO TO 110
         IF (MES.NE.MESTAB(I)) GO TO 90
         IF (NERR.NE.NERTAB(I)) GO TO 90
         IF (LEVEL.NE.LEVTAB(I)) GO TO 90
         GO TO 100
   90 CONTINUE
C     THREE POSSIBLE CASES...
C     TABLE IS FULL
         KOUNTX = KOUNTX+1
         ICOUNT = 1
         RETURN
C     MESSAGE FOUND IN TABLE
  100    KOUNT(II) = KOUNT(II) + 1
         ICOUNT = KOUNT(II)
         RETURN
C     EMPTY SLOT FOUND FOR NEW MESSAGE
  110    MESTAB(II) = MES
         NERTAB(II) = NERR
         LEVTAB(II) = LEVEL
         KOUNT(II)  = 1
         ICOUNT = 1
         RETURN
      END
      SUBROUTINE XGETUA(IUNITA,N)
C***BEGIN PROLOGUE  XGETUA
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Returns unit number(s) to which error messages are being
C            sent.
C***DESCRIPTION
C     Abstract
C        XGETUA may be called to determine the unit number or numbers
C        to which error messages are being sent.
C        These unit numbers may have been set by a call to XSETUN,
C        or a call to XSETUA, or may be a default value.
C
C     Description of Parameters
C      --Output--
C        IUNIT - an array of one to five unit numbers, depending
C                on the value of N.  A value of zero refers to the
C                default unit, as defined by the I1MACH machine
C                constant routine.  Only IUNIT(1),...,IUNIT(N) are
C                defined by XGETUA.  The values of IUNIT(N+1),...,
C                IUNIT(5) are not defined (for N .LT. 5) or altered
C                in any way by XGETUA.
C        N     - the number of units to which copies of the
C                error messages are being sent.  N will be in the
C                range from 1 to 5.
C
C     Latest revision ---  19 MAR 1980
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  J4SAVE
C***END PROLOGUE  XGETUA
      DIMENSION IUNITA(5)
C***FIRST EXECUTABLE STATEMENT  XGETUA
      N = J4SAVE(5,0,.FALSE.)
      DO 30 I=1,N
         INDEX = I+4
         IF (I.EQ.1) INDEX = 3
         IUNITA(I) = J4SAVE(INDEX,0,.FALSE.)
   30 CONTINUE
      RETURN
      END

      
      SUBROUTINE DBINTK(X,Y,T,N,K,BCOEF,Q,WORK)
C***BEGIN PROLOGUE  DBINTK
C***DATE WRITTEN   800901   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***CATEGORY NO.  E1A
C***KEYWORDS  B-SPLINE,DATA FITTING,DOUBLE PRECISION,INTERPOLATION,
C             SPLINE
C***AUTHOR  AMOS, D. E., (SNLA)
C***PURPOSE  Produces the B-spline coefficients, BCOEF, of the
C            B-spline of order K with knots T(I), I=1,...,N+K, which
C            takes on the value Y(I) at X(I), I=1,...,N.
C***DESCRIPTION
C
C     Written by Carl de Boor and modified by D. E. Amos
C
C     References
C
C         A Practical Guide to Splines by C. de Boor, Applied
C         Mathematics Series 27, Springer, 1979.
C
C     Abstract    **** a double precision routine ****
C
C         DBINTK is the SPLINT routine of the reference.
C
C         DBINTK produces the B-spline coefficients, BCOEF, of the
C         B-spline of order K with knots T(I), I=1,...,N+K, which
C         takes on the value Y(I) at X(I), I=1,...,N.  The spline or
C         any of its derivatives can be evaluated by calls to DBVALU.
C
C         The I-th equation of the linear system A*BCOEF = B for the
C         coefficients of the interpolant enforces interpolation at
C         X(I), I=1,...,N.  Hence, B(I) = Y(I), for all I, and A is
C         a band matrix with 2K-1 bands if A is invertible.  The matrix
C         A is generated row by row and stored, diagonal by diagonal,
C         in the rows of Q, with the main diagonal going into row K.
C         The banded system is then solved by a call to DBNFAC (which
C         constructs the triangular factorization for A and stores it
C         again in Q), followed by a call to DBNSLV (which then
C         obtains the solution BCOEF by substitution).  DBNFAC does no
C         pivoting, since the total positivity of the matrix A makes
C         this unnecessary.  The linear system to be solved is
C         (theoretically) invertible if and only if
C                 T(I) .LT. X(I) .LT. T(I+K),        for all I.
C         Equality is permitted on the left for I=1 and on the right
C         for I=N when K knots are used at X(1) or X(N).  Otherwise,
C         violation of this condition is certain to lead to an error.
C
C         DBINTK calls DBSPVN, DBNFAC, DBNSLV, XERROR
C
C     Description of Arguments
C
C         Input       X,Y,T are double precision
C           X       - vector of length N containing data point abscissa
C                     in strictly increasing order.
C           Y       - corresponding vector of length N containing data
C                     point ordinates.
C           T       - knot vector of length N+K
C                     Since T(1),..,T(K) .LE. X(1) and T(N+1),..,T(N+K)
C                     .GE. X(N), this leaves only N-K knots (not nec-
C                     essarily X(I) values) interior to (X(1),X(N))
C           N       - number of data points, N .GE. K
C           K       - order of the spline, K .GE. 1
C
C         Output      BCOEF,Q,WORK are double precision
C           BCOEF   - a vector of length N containing the B-spline
C                     coefficients
C           Q       - a work vector of length (2*K-1)*N, containing
C                     the triangular factorization of the coefficient
C                     matrix of the linear system being solved.  The
C                     coefficients for the interpolant of an
C                     additional data set (X(I),yY(I)), I=1,...,N
C                     with the same abscissa can be obtained by loading
C                     YY into BCOEF and then executing
C                         CALL DBNSLV(Q,2K-1,N,K-1,K-1,BCOEF)
C           WORK    - work vector of length 2*K
C
C     Error Conditions
C         Improper input is a fatal error
C         Singular system of equations is a fatal error
C***REFERENCES  C. DE BOOR, *A PRACTICAL GUIDE TO SPLINES*, APPLIED
C                 MATHEMATICS SERIES 27, SPRINGER, 1979.
C               D.E. AMOS, *COMPUTATION WITH SPLINES AND B-SPLINES*,
C                 SAND78-1968,SANDIA LABORATORIES,MARCH,1979.
C***ROUTINES CALLED  DBNFAC,DBNSLV,DBSPVN,XERROR
C***END PROLOGUE  DBINTK
C
C
      INTEGER IFLAG, IWORK, K, N, I, ILP1MX, J, JJ, KM1, KPKM2, LEFT,
     1 LENQ, NP1
      DOUBLE PRECISION BCOEF(N), Y(N), Q(*), T(*), X(N), XI, WORK(*)
C     DIMENSION Q(2*K-1,N), T(N+K)
C***FIRST EXECUTABLE STATEMENT  DBINTK
      IF(K.LT.1) GO TO 100
      IF(N.LT.K) GO TO 105
      JJ = N - 1
      IF(JJ.EQ.0) GO TO 6
      DO 5 I=1,JJ
      IF(X(I).GE.X(I+1)) GO TO 110
    5 CONTINUE
    6 CONTINUE
      NP1 = N + 1
      KM1 = K - 1
      KPKM2 = 2*KM1
      LEFT = K
C                ZERO OUT ALL ENTRIES OF Q
      LENQ = N*(K+KM1)
      DO 10 I=1,LENQ
        Q(I) = 0.0D0
   10 CONTINUE
C
C  ***   LOOP OVER I TO CONSTRUCT THE  N  INTERPOLATION EQUATIONS
      DO 50 I=1,N
        XI = X(I)
        ILP1MX = MIN0(I+K,NP1)
C        *** FIND  LEFT  IN THE CLOSED INTERVAL (I,I+K-1) SUCH THAT
C                T(LEFT) .LE. X(I) .LT. T(LEFT+1)
C        MATRIX IS SINGULAR IF THIS IS NOT POSSIBLE
        LEFT = MAX0(LEFT,I)
        IF (XI.LT.T(LEFT)) GO TO 80
   20   IF (XI.LT.T(LEFT+1)) GO TO 30
        LEFT = LEFT + 1
        IF (LEFT.LT.ILP1MX) GO TO 20
        LEFT = LEFT - 1
        IF (XI.GT.T(LEFT+1)) GO TO 80
C        *** THE I-TH EQUATION ENFORCES INTERPOLATION AT XI, HENCE
C        A(I,J) = B(J,K,T)(XI), ALL J. ONLY THE  K  ENTRIES WITH  J =
C        LEFT-K+1,...,LEFT ACTUALLY MIGHT BE NONZERO. THESE  K  NUMBERS
C        ARE RETURNED, IN  BCOEF (USED FOR TEMP.STORAGE HERE), BY THE
C        FOLLOWING
   30   CALL DBSPVN(T, K, K, 1, XI, LEFT, BCOEF, WORK, IWORK)
C        WE THEREFORE WANT  BCOEF(J) = B(LEFT-K+J)(XI) TO GO INTO
C        A(I,LEFT-K+J), I.E., INTO  Q(I-(LEFT+J)+2*K,(LEFT+J)-K) SINCE
C        A(I+J,J)  IS TO GO INTO  Q(I+K,J), ALL I,J,  IF WE CONSIDER  Q
C        AS A TWO-DIM. ARRAY , WITH  2*K-1  ROWS (SEE COMMENTS IN
C        DBNFAC). IN THE PRESENT PROGRAM, WE TREAT  Q  AS AN EQUIVALENT
C        ONE-DIMENSIONAL ARRAY (BECAUSE OF FORTRAN RESTRICTIONS ON
C        DIMENSION STATEMENTS) . WE THEREFORE WANT  BCOEF(J) TO GO INTO
C        ENTRY
C            I -(LEFT+J) + 2*K + ((LEFT+J) - K-1)*(2*K-1)
C                   =  I-LEFT+1 + (LEFT -K)*(2*K-1) + (2*K-2)*J
C        OF  Q .
        JJ = I - LEFT + 1 + (LEFT-K)*(K+KM1)
        DO 40 J=1,K
          JJ = JJ + KPKM2
          Q(JJ) = BCOEF(J)
   40   CONTINUE
   50 CONTINUE
C
C     ***OBTAIN FACTORIZATION OF  A  , STORED AGAIN IN  Q.
      CALL DBNFAC(Q, K+KM1, N, KM1, KM1, IFLAG)
      GO TO (60, 90), IFLAG
C     *** SOLVE  A*BCOEF = Y  BY BACKSUBSTITUTION
   60 DO 70 I=1,N
        BCOEF(I) = Y(I)
   70 CONTINUE
      CALL DBNSLV(Q, K+KM1, N, KM1, KM1, BCOEF)
      RETURN
C
C
   80 CONTINUE
      CALL XERROR( ' DBINTK,  SOME ABSCISSA WAS NOT IN THE SUPPORT OF TH
     1E CORRESPONDING BASIS FUNCTION AND THE SYSTEM IS SINGULAR.',109,2,
     21)
      RETURN
   90 CONTINUE
      CALL XERROR( ' DBINTK,  THE SYSTEM OF SOLVER DETECTS A SINGULAR SY
     1STEM ALTHOUGH THE THEORETICAL CONDITIONS FOR A SOLUTION WERE SATIS
     2FIED.',123,8,1)
      RETURN
  100 CONTINUE
      CALL XERROR( ' DBINTK,  K DOES NOT SATISFY K.GE.1', 35, 2, 1)
      RETURN
  105 CONTINUE
      CALL XERROR( ' DBINTK,  N DOES NOT SATISFY N.GE.K', 35, 2, 1)
      RETURN
  110 CONTINUE
      CALL XERROR( ' DBINTK,  X(I) DOES NOT SATISFY X(I).LT.X(I+1) FOR S
     1OME I', 57, 2, 1)
      RETURN
      END
      SUBROUTINE DBNFAC(W,NROWW,NROW,NBANDL,NBANDU,IFLAG)
C***BEGIN PROLOGUE  DBNFAC
C***REFER TO  DBINT4,DBINTK
C
C  DBNFAC is the BANFAC routine from
C        * A Practical Guide to Splines *  by C. de Boor
C
C  DBNFAC is a double precision routine
C
C  Returns in  W  the LU-factorization (without pivoting) of the banded
C  matrix  A  of order  NROW  with  (NBANDL + 1 + NBANDU) bands or diag-
C  onals in the work array  W .
C
C *****  I N P U T  ****** W is double precision
C  W.....Work array of size  (NROWW,NROW)  containing the interesting
C        part of a banded matrix  A , with the diagonals or bands of  A
C        stored in the rows of  W , while columns of  A  correspond to
C        columns of  W . This is the storage mode used in  LINPACK  and
C        results in efficient innermost loops.
C           Explicitly,  A  has  NBANDL  bands below the diagonal
C                            +     1     (main) diagonal
C                            +   NBANDU  bands above the diagonal
C        and thus, with    MIDDLE = NBANDU + 1,
C          A(I+J,J)  is in  W(I+MIDDLE,J)  for I=-NBANDU,...,NBANDL
C                                              J=1,...,NROW .
C        For example, the interesting entries of A (1,2)-banded matrix
C        of order  9  would appear in the first  1+1+2 = 4  rows of  W
C        as follows.
C                          13 24 35 46 57 68 79
C                       12 23 34 45 56 67 78 89
C                    11 22 33 44 55 66 77 88 99
C                    21 32 43 54 65 76 87 98
C
C        All other entries of  W  not identified in this way with an en-
C        try of  A  are never referenced .
C  NROWW.....Row dimension of the work array  W .
C        must be  .GE.  NBANDL + 1 + NBANDU  .
C  NBANDL.....Number of bands of  A  below the main diagonal
C  NBANDU.....Number of bands of  A  above the main diagonal .
C
C *****  O U T P U T  ****** W is double precision
C  IFLAG.....Integer indicating success( = 1) or failure ( = 2) .
C     If  IFLAG = 1, then
C  W.....contains the LU-factorization of  A  into a unit lower triangu-
C        lar matrix  L  and an upper triangular matrix  U (both banded)
C        and stored in customary fashion over the corresponding entries
C        of  A . This makes it possible to solve any particular linear
C        system  A*X = B  for  X  by a
C              CALL DBNSLV ( W, NROWW, NROW, NBANDL, NBANDU, B )
C        with the solution X  contained in  B  on return .
C     If  IFLAG = 2, then
C        one of  NROW-1, NBANDL,NBANDU failed to be nonnegative, or else
C        one of the potential pivots was found to be zero indicating
C        that  A  does not have an LU-factorization. This implies that
C        A  is singular in case it is totally positive .
C
C *****  M E T H O D  ******
C     Gauss elimination  W I T H O U T  pivoting is used. The routine is
C  intended for use with matrices  A  which do not require row inter-
C  changes during factorization, especially for the  T O T A L L Y
C  P O S I T I V E  matrices which occur in spline calculations.
C     The routine should NOT be used for an arbitrary banded matrix.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DBNFAC
C
      INTEGER IFLAG, NBANDL, NBANDU, NROW, NROWW, I, IPK, J, JMAX, K,
     1 KMAX, MIDDLE, MIDMK, NROWM1
      DOUBLE PRECISION W(NROWW,NROW), FACTOR, PIVOT
C
C***FIRST EXECUTABLE STATEMENT  DBNFAC
      IFLAG = 1
      MIDDLE = NBANDU + 1
C                         W(MIDDLE,.) CONTAINS THE MAIN DIAGONAL OF  A .
      NROWM1 = NROW - 1
      IF (NROWM1) 120, 110, 10
   10 IF (NBANDL.GT.0) GO TO 30
C                A IS UPPER TRIANGULAR. CHECK THAT DIAGONAL IS NONZERO .
      DO 20 I=1,NROWM1
        IF (W(MIDDLE,I).EQ.0.0D0) GO TO 120
   20 CONTINUE
      GO TO 110
   30 IF (NBANDU.GT.0) GO TO 60
C              A IS LOWER TRIANGULAR. CHECK THAT DIAGONAL IS NONZERO AND
C                 DIVIDE EACH COLUMN BY ITS DIAGONAL .
      DO 50 I=1,NROWM1
        PIVOT = W(MIDDLE,I)
        IF (PIVOT.EQ.0.0D0) GO TO 120
        JMAX = MIN0(NBANDL,NROW-I)
        DO 40 J=1,JMAX
          W(MIDDLE+J,I) = W(MIDDLE+J,I)/PIVOT
   40   CONTINUE
   50 CONTINUE
      RETURN
C
C        A  IS NOT JUST A TRIANGULAR MATRIX. CONSTRUCT LU FACTORIZATION
   60 DO 100 I=1,NROWM1
C                                  W(MIDDLE,I)  IS PIVOT FOR I-TH STEP .
        PIVOT = W(MIDDLE,I)
        IF (PIVOT.EQ.0.0D0) GO TO 120
C                 JMAX  IS THE NUMBER OF (NONZERO) ENTRIES IN COLUMN  I
C                     BELOW THE DIAGONAL .
        JMAX = MIN0(NBANDL,NROW-I)
C              DIVIDE EACH ENTRY IN COLUMN  I  BELOW DIAGONAL BY PIVOT .
        DO 70 J=1,JMAX
          W(MIDDLE+J,I) = W(MIDDLE+J,I)/PIVOT
   70   CONTINUE
C                 KMAX  IS THE NUMBER OF (NONZERO) ENTRIES IN ROW  I  TO
C                     THE RIGHT OF THE DIAGONAL .
        KMAX = MIN0(NBANDU,NROW-I)
C                  SUBTRACT  A(I,I+K)*(I-TH COLUMN) FROM (I+K)-TH COLUMN
C                  (BELOW ROW  I ) .
        DO 90 K=1,KMAX
          IPK = I + K
          MIDMK = MIDDLE - K
          FACTOR = W(MIDMK,IPK)
          DO 80 J=1,JMAX
            W(MIDMK+J,IPK) = W(MIDMK+J,IPK) - W(MIDDLE+J,I)*FACTOR
   80     CONTINUE
   90   CONTINUE
  100 CONTINUE
C                                       CHECK THE LAST DIAGONAL ENTRY .
  110 IF (W(MIDDLE,NROW).NE.0.0D0) RETURN
  120 IFLAG = 2
      RETURN
      END
      SUBROUTINE DBNSLV(W,NROWW,NROW,NBANDL,NBANDU,B)
C***BEGIN PROLOGUE  DBNSLV
C***REFER TO  DBINT4,DBINTK
C
C  DBNSLV is the BANSLV routine from
C        * A Practical Guide to Splines *  by C. de Boor
C
C  DBNSLV is a double precision routine
C
C  Companion routine to  DBNFAC . It returns the solution  X  of the
C  linear system  A*X = B  in place of  B , given the LU-factorization
C  for  A  in the work array  W from DBNFAC.
C
C *****  I N P U T  ****** W,B are DOUBLE PRECISION
C  W, NROWW,NROW,NBANDL,NBANDU.....Describe the LU-factorization of a
C        banded matrix  A  of order  NROW  as constructed in  DBNFAC .
C        For details, see  DBNFAC .
C  B.....Right side of the system to be solved .
C
C *****  O U T P U T  ****** B is DOUBLE PRECISION
C  B.....Contains the solution  X , of order  NROW .
C
C *****  M E T H O D  ******
C     (With  A = L*U, as stored in  W,) the unit lower triangular system
C  L(U*X) = B  is solved for  Y = U*X, and  Y  stored in  B . Then the
C  upper triangular system  U*X = Y  is solved for  X  . The calcul-
C  ations are so arranged that the innermost loops stay within columns.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DBNSLV
C
      INTEGER NBANDL, NBANDU, NROW, NROWW, I, J, JMAX, MIDDLE, NROWM1
      DOUBLE PRECISION W(NROWW,NROW), B(NROW)
C***FIRST EXECUTABLE STATEMENT  DBNSLV
      MIDDLE = NBANDU + 1
      IF (NROW.EQ.1) GO TO 80
      NROWM1 = NROW - 1
      IF (NBANDL.EQ.0) GO TO 30
C                                 FORWARD PASS
C            FOR I=1,2,...,NROW-1, SUBTRACT  RIGHT SIDE(I)*(I-TH COLUMN
C            OF  L )  FROM RIGHT SIDE  (BELOW I-TH ROW) .
      DO 20 I=1,NROWM1
        JMAX = MIN0(NBANDL,NROW-I)
        DO 10 J=1,JMAX
          B(I+J) = B(I+J) - B(I)*W(MIDDLE+J,I)
   10   CONTINUE
   20 CONTINUE
C                                 BACKWARD PASS
C            FOR I=NROW,NROW-1,...,1, DIVIDE RIGHT SIDE(I) BY I-TH DIAG-
C            ONAL ENTRY OF  U, THEN SUBTRACT  RIGHT SIDE(I)*(I-TH COLUMN
C            OF  U)  FROM RIGHT SIDE  (ABOVE I-TH ROW).
   30 IF (NBANDU.GT.0) GO TO 50
C                                A  IS LOWER TRIANGULAR .
      DO 40 I=1,NROW
        B(I) = B(I)/W(1,I)
   40 CONTINUE
      RETURN
   50 I = NROW
   60 B(I) = B(I)/W(MIDDLE,I)
      JMAX = MIN0(NBANDU,I-1)
      DO 70 J=1,JMAX
        B(I-J) = B(I-J) - B(I)*W(MIDDLE-J,I)
   70 CONTINUE
      I = I - 1
      IF (I.GT.1) GO TO 60
   80 B(1) = B(1)/W(MIDDLE,1)
      RETURN
      END
      SUBROUTINE DBSPVN(T,JHIGH,K,INDEX,X,ILEFT,VNIKX,WORK,IWORK)
C***BEGIN PROLOGUE  DBSPVN
C***DATE WRITTEN   800901   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***CATEGORY NO.  E3,K6
C***KEYWORDS  B-SPLINE,DATA FITTING,DOUBLE PRECISION,INTERPOLATION,
C             SPLINE
C***AUTHOR  AMOS, D. E., (SNLA)
C***PURPOSE  Calculates the value of all (possibly) nonzero basis
C            functions at X.
C***DESCRIPTION
C
C     Written by Carl de Boor and modified by D. E. Amos
C
C     Reference
C         SIAM J. Numerical Analysis, 14, No. 3, June, 1977, pp.441-472.
C
C     Abstract    **** a double precision routine ****
C         DBSPVN is the BSPLVN routine of the reference.
C
C         DBSPVN calculates the value of all (possibly) nonzero basis
C         functions at X of order MAX(JHIGH,(J+1)*(INDEX-1)), where T(K)
C         .LE. X .LE. T(N+1) and J=IWORK is set inside the routine on
C         the first call when INDEX=1.  ILEFT is such that T(ILEFT) .LE.
C         X .LT. T(ILEFT+1).  A call to DINTRV(T,N+1,X,ILO,ILEFT,MFLAG)
C         produces the proper ILEFT.  DBSPVN calculates using the basic
C         algorithm needed in DBSPVD.  If only basis functions are
C         desired, setting JHIGH=K and INDEX=1 can be faster than
C         calling DBSPVD, but extra coding is required for derivatives
C         (INDEX=2) and DBSPVD is set up for this purpose.
C
C         Left limiting values are set up as described in DBSPVD.
C
C     Description of Arguments
C
C         Input      T,X are double precision
C          T       - knot vector of length N+K, where
C                    N = number of B-spline basis functions
C                    N = sum of knot multiplicities-K
C          JHIGH   - order of B-spline, 1 .LE. JHIGH .LE. K
C          K       - highest possible order
C          INDEX   - INDEX = 1 gives basis functions of order JHIGH
C                          = 2 denotes previous entry with work, IWORK
C                              values saved for subsequent calls to
C                              DBSPVN.
C          X       - argument of basis functions,
C                    T(K) .LE. X .LE. T(N+1)
C          ILEFT   - largest integer such that
C                    T(ILEFT) .LE. X .LT.  T(ILEFT+1)
C
C         Output     VNIKX, WORK are double precision
C          VNIKX   - vector of length K for spline values.
C          WORK    - a work vector of length 2*K
C          IWORK   - a work parameter.  Both WORK and IWORK contain
C                    information necessary to continue for INDEX = 2.
C                    When INDEX = 1 exclusively, these are scratch
C                    variables and can be used for other purposes.
C
C     Error Conditions
C         Improper input is a fatal error.
C***REFERENCES  C. DE BOOR, *PACKAGE FOR CALCULATING WITH B-SPLINES*,
C                 SIAM JOURNAL ON NUMERICAL ANALYSIS, VOLUME 14, NO. 3,
C                 JUNE 1977, PP. 441-472.
C***ROUTINES CALLED  XERROR
C***END PROLOGUE  DBSPVN
C
C
      INTEGER ILEFT, IMJP1, INDEX, IPJ, IWORK, JHIGH, JP1, JP1ML, K, L
      DOUBLE PRECISION T, VM, VMPREV, VNIKX, WORK, X
C     DIMENSION T(ILEFT+JHIGH)
      DIMENSION T(*), VNIKX(K), WORK(*)
C     CONTENT OF J, DELTAM, DELTAP IS EXPECTED UNCHANGED BETWEEN CALLS.
C     WORK(I) = DELTAP(I), WORK(K+I) = DELTAM(I), I = 1,K
C***FIRST EXECUTABLE STATEMENT  DBSPVN
      IF(K.LT.1) GO TO 90
      IF(JHIGH.GT.K .OR. JHIGH.LT.1) GO TO 100
      IF(INDEX.LT.1 .OR. INDEX.GT.2) GO TO 105
      IF(X.LT.T(ILEFT) .OR. X.GT.T(ILEFT+1)) GO TO 110
      GO TO (10, 20), INDEX
   10 IWORK = 1
      VNIKX(1) = 1.0D0
      IF (IWORK.GE.JHIGH) GO TO 40
C
   20 IPJ = ILEFT + IWORK
      WORK(IWORK) = T(IPJ) - X
      IMJP1 = ILEFT - IWORK + 1
      WORK(K+IWORK) = X - T(IMJP1)
      VMPREV = 0.0D0
      JP1 = IWORK + 1
      DO 30 L=1,IWORK
        JP1ML = JP1 - L
        VM = VNIKX(L)/(WORK(L)+WORK(K+JP1ML))
        VNIKX(L) = VM*WORK(L) + VMPREV
        VMPREV = VM*WORK(K+JP1ML)
   30 CONTINUE
      VNIKX(JP1) = VMPREV
      IWORK = JP1
      IF (IWORK.LT.JHIGH) GO TO 20
C
   40 RETURN
C
C
   90 CONTINUE
      CALL XERROR( ' DBSPVN,  K DOES NOT SATISFY K.GE.1', 35, 2, 1)
      RETURN
  100 CONTINUE
      CALL XERROR( ' DBSPVN,  JHIGH DOES NOT SATISFY 1.LE.JHIGH.LE.K',
     1 48, 2, 1)
      RETURN
  105 CONTINUE
      CALL XERROR( ' DBSPVN,  INDEX IS NOT 1 OR 2',29,2,1)
      RETURN
  110 CONTINUE
      CALL XERROR( ' DBSPVN,  X DOES NOT SATISFY T(ILEFT).LE.X.LE.T(ILEF
     1T+1)', 56, 2, 1)
      RETURN
      END
      SUBROUTINE FDUMP
C***BEGIN PROLOGUE  FDUMP
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  Z
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Symbolic dump (should be locally written).
C***DESCRIPTION
C        ***Note*** Machine Dependent Routine
C        FDUMP is intended to be replaced by a locally written
C        version which produces a symbolic dump.  Failing this,
C        it should be replaced by a version which prints the
C        subprogram nesting list.  Note that this dump must be
C        printed on each of up to five files, as indicated by the
C        XGETUA routine.  See XSETUA and XGETUA for details.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C     Latest revision ---  23 May 1979
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  FDUMP
C***FIRST EXECUTABLE STATEMENT  FDUMP
      RETURN
      END
      FUNCTION J4SAVE(IWHICH,IVALUE,ISET)
C***BEGIN PROLOGUE  J4SAVE
C***REFER TO  XERROR
C     Abstract
C        J4SAVE saves and recalls several global variables needed
C        by the library error handling routines.
C
C     Description of Parameters
C      --Input--
C        IWHICH - Index of item desired.
C                = 1 Refers to current error number.
C                = 2 Refers to current error control flag.
C                 = 3 Refers to current unit number to which error
C                    messages are to be sent.  (0 means use standard.)
C                 = 4 Refers to the maximum number of times any
C                     message is to be printed (as set by XERMAX).
C                 = 5 Refers to the total number of units to which
C                     each error message is to be written.
C                 = 6 Refers to the 2nd unit for error messages
C                 = 7 Refers to the 3rd unit for error messages
C                 = 8 Refers to the 4th unit for error messages
C                 = 9 Refers to the 5th unit for error messages
C        IVALUE - The value to be set for the IWHICH-th parameter,
C                 if ISET is .TRUE. .
C        ISET   - If ISET=.TRUE., the IWHICH-th parameter will BE
C                 given the value, IVALUE.  If ISET=.FALSE., the
C                 IWHICH-th parameter will be unchanged, and IVALUE
C                 is a dummy parameter.
C      --Output--
C        The (old) value of the IWHICH-th parameter will be returned
C        in the function value, J4SAVE.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C    Adapted from Bell Laboratories PORT Library Error Handler
C     Latest revision ---  23 MAY 1979
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  J4SAVE
      LOGICAL ISET
      INTEGER IPARAM(9)
      SAVE IPARAM
      DATA IPARAM(1),IPARAM(2),IPARAM(3),IPARAM(4)/0,2,0,10/
      DATA IPARAM(5)/1/
      DATA IPARAM(6),IPARAM(7),IPARAM(8),IPARAM(9)/0,0,0,0/
C***FIRST EXECUTABLE STATEMENT  J4SAVE
      J4SAVE = IPARAM(IWHICH)
      IF (ISET) IPARAM(IWHICH) = IVALUE
      RETURN
      END
      SUBROUTINE XERABT(MESSG,NMESSG)
C***BEGIN PROLOGUE  XERABT
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Aborts program execution and prints error message.
C***DESCRIPTION
C     Abstract
C        ***Note*** machine dependent routine
C        XERABT aborts the execution of the program.
C        The error message causing the abort is given in the calling
C        sequence, in case one needs it for printing on a dayfile,
C        for example.
C
C     Description of Parameters
C        MESSG and NMESSG are as in XERROR, except that NMESSG may
C        be zero, in which case no message is being supplied.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C     Latest revision ---  19 MAR 1980
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  XERABT
      CHARACTER*(*) MESSG
C***FIRST EXECUTABLE STATEMENT  XERABT
      STOP
      END
      SUBROUTINE XERCTL(MESSG1,NMESSG,NERR,LEVEL,KONTRL)
C***BEGIN PROLOGUE  XERCTL
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Allows user control over handling of individual errors.
C***DESCRIPTION
C     Abstract
C        Allows user control over handling of individual errors.
C        Just after each message is recorded, but before it is
C        processed any further (i.e., before it is printed or
C        a decision to abort is made), a call is made to XERCTL.
C        If the user has provided his own version of XERCTL, he
C        can then override the value of KONTROL used in processing
C        this message by redefining its value.
C        KONTRL may be set to any value from -2 to 2.
C        The meanings for KONTRL are the same as in XSETF, except
C        that the value of KONTRL changes only for this message.
C        If KONTRL is set to a value outside the range from -2 to 2,
C        it will be moved back into that range.
C
C     Description of Parameters
C
C      --Input--
C        MESSG1 - the first word (only) of the error message.
C        NMESSG - same as in the call to XERROR or XERRWV.
C        NERR   - same as in the call to XERROR or XERRWV.
C        LEVEL  - same as in the call to XERROR or XERRWV.
C        KONTRL - the current value of the control flag as set
C                 by a call to XSETF.
C
C      --Output--
C        KONTRL - the new value of KONTRL.  If KONTRL is not
C                 defined, it will remain at its original value.
C                 This changed value of control affects only
C                 the current occurrence of the current message.
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  XERCTL
      CHARACTER*20 MESSG1
C***FIRST EXECUTABLE STATEMENT  XERCTL
      RETURN
      END
      SUBROUTINE XERPRT(MESSG,NMESSG)
C***BEGIN PROLOGUE  XERPRT
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  Z
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Prints error messages.
C***DESCRIPTION
C     Abstract
C        Print the Hollerith message in MESSG, of length NMESSG,
C        on each file indicated by XGETUA.
C     Latest revision ---  19 MAR 1980
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  I1MACH,S88FMT,XGETUA
C***END PROLOGUE  XERPRT
      INTEGER LUN(5)
      CHARACTER*(*) MESSG
C     OBTAIN UNIT NUMBERS AND WRITE LINE TO EACH UNIT
C***FIRST EXECUTABLE STATEMENT  XERPRT
      CALL XGETUA(LUN,NUNIT)
      LENMES = LEN(MESSG)
      DO 20 KUNIT=1,NUNIT
         IUNIT = LUN(KUNIT)
         IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
         DO 10 ICHAR=1,LENMES,72
            LAST = MIN0(ICHAR+71 , LENMES)
            WRITE (IUNIT,'(1X,A)') MESSG(ICHAR:LAST)
   10    CONTINUE
   20 CONTINUE
      RETURN
      END
      SUBROUTINE XERROR(MESSG,NMESSG,NERR,LEVEL)
C***BEGIN PROLOGUE  XERROR
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Processes an error (diagnostic) message.
C***DESCRIPTION
C     Abstract
C        XERROR processes a diagnostic message, in a manner
C        determined by the value of LEVEL and the current value
C        of the library error control flag, KONTRL.
C        (See subroutine XSETF for details.)
C
C     Description of Parameters
C      --Input--
C        MESSG - the Hollerith message to be processed, containing
C                no more than 72 characters.
C        NMESSG- the actual number of characters in MESSG.
C        NERR  - the error number associated with this message.
C                NERR must not be zero.
C        LEVEL - error category.
C                =2 means this is an unconditionally fatal error.
C                =1 means this is a recoverable error.  (I.e., it is
C                   non-fatal if XSETF has been appropriately called.)
C                =0 means this is a warning message only.
C                =-1 means this is a warning message which is to be
C                   printed at most once, regardless of how many
C                   times this call is executed.
C
C     Examples
C        CALL XERROR('SMOOTH -- NUM WAS ZERO.',23,1,2)
C        CALL XERROR('INTEG  -- LESS THAN FULL ACCURACY ACHIEVED.',
C                    43,2,1)
C        CALL XERROR('ROOTER -- ACTUAL ZERO OF F FOUND BEFORE INTERVAL F
C    1ULLY COLLAPSED.',65,3,0)
C        CALL XERROR('EXP    -- UNDERFLOWS BEING SET TO ZERO.',39,1,-1)
C
C     Latest revision ---  19 MAR 1980
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  XERRWV
C***END PROLOGUE  XERROR
      CHARACTER*(*) MESSG
C***FIRST EXECUTABLE STATEMENT  XERROR
      CALL XERRWV(MESSG,NMESSG,NERR,LEVEL,0,0,0,0,0.,0.)
      RETURN
      END


      DOUBLE PRECISION FUNCTION DB3VAL(XVAL,YVAL,ZVAL,IDX,IDY,IDZ,
     *  TX,TY,TZ,NX,NY,NZ,KX,KY,KZ,BCOEF,WORK)
C***BEGIN PROLOGUE  DB3VAL
C***DATE WRITTEN   25 MAY 1982
C***REVISION DATE  25 MAY 1982
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***CATEGORY NO.  E1A
C***KEYWORDS  INTERPOLATION, THREE-DIMENSIONS, GRIDDED DATA, SPLINES,
C             PIECEWISE POLYNOMIALS
C***AUTHOR  BOISVERT, RONALD, NBS
C             SCIENTIFIC COMPUTING DIVISION
C             NATIONAL BUREAU OF STANDARDS
C             WASHINGTON, DC 20234
C***PURPOSE  DB3VAL EVALUATES THE PIECEWISE POLYNOMIAL INTERPOLATING
C            FUNCTION CONSTRUCTED BY THE ROUTINE B3INK OR ONE OF ITS
C            PARTIAL DERIVATIVES.
C            DOUBLE PRECISION VERSION OF B3VAL.
C***DESCRIPTION
C
C   DB3VAL  evaluates   the   tensor   product   piecewise   polynomial
C   interpolant constructed  by  the  routine  DB3INK  or  one  of  its
C   derivatives  at  the  point  (XVAL,YVAL,ZVAL).  To   evaluate   the
C   interpolant  itself,  set  IDX=IDY=IDZ=0,  to  evaluate  the  first
C   partial with respect to x, set IDX=1,IDY=IDZ=0, and so on.
C
C   DB3VAL returns 0.0D0 if (XVAL,YVAL,ZVAL) is out of range. That is,
C            XVAL.LT.TX(1) .OR. XVAL.GT.TX(NX+KX) .OR.
C            YVAL.LT.TY(1) .OR. YVAL.GT.TY(NY+KY) .OR.
C            ZVAL.LT.TZ(1) .OR. ZVAL.GT.TZ(NZ+KZ)
C   If the knots TX, TY, and TZ were chosen by  DB3INK,  then  this  is
C   equivalent to
C            XVAL.LT.X(1) .OR. XVAL.GT.X(NX)+EPSX .OR.
C            YVAL.LT.Y(1) .OR. YVAL.GT.Y(NY)+EPSY .OR.
C            ZVAL.LT.Z(1) .OR. ZVAL.GT.Z(NZ)+EPSZ
C   where EPSX = 0.1*(X(NX)-X(NX-1)), EPSY =  0.1*(Y(NY)-Y(NY-1)),  and
C   EPSZ = 0.1*(Z(NZ)-Z(NZ-1)).
C
C   The input quantities TX, TY, TZ, NX, NY, NZ, KX, KY, KZ, and  BCOEF
C   should remain unchanged since the last call of DB3INK.
C
C
C   I N P U T
C   ---------
C
C   XVAL    Double precision scalar
C           X coordinate of evaluation point.
C
C   YVAL    Double precision scalar
C           Y coordinate of evaluation point.
C
C   ZVAL    Double precision scalar
C           Z coordinate of evaluation point.
C
C   IDX     Integer scalar
C           X derivative of piecewise polynomial to evaluate.
C
C   IDY     Integer scalar
C           Y derivative of piecewise polynomial to evaluate.
C
C   IDZ     Integer scalar
C           Z derivative of piecewise polynomial to evaluate.
C
C   TX      Double precision 1D array (size NX+KX)
C           Sequence of knots defining the piecewise polynomial in
C           the x direction.  (Same as in last call to DB3INK.)
C
C   TY      Double precision 1D array (size NY+KY)
C           Sequence of knots defining the piecewise polynomial in
C           the y direction.  (Same as in last call to DB3INK.)
C
C   TZ      Double precision 1D array (size NZ+KZ)
C           Sequence of knots defining the piecewise polynomial in
C           the z direction.  (Same as in last call to DB3INK.)
C
C   NX      Integer scalar
C           The number of interpolation points in x.
C           (Same as in last call to DB3INK.)
C
C   NY      Integer scalar
C           The number of interpolation points in y.
C           (Same as in last call to DB3INK.)
C
C   NZ      Integer scalar
C           The number of interpolation points in z.
C           (Same as in last call to DB3INK.)
C
C   KX      Integer scalar
C           Order of polynomial pieces in x.
C           (Same as in last call to DB3INK.)
C
C   KY      Integer scalar
C           Order of polynomial pieces in y.
C           (Same as in last call to DB3INK.)
C
C   KZ      Integer scalar
C           Order of polynomial pieces in z.
C           (Same as in last call to DB3INK.)
C
C   BCOEF   Double precision 2D array (size NX by NY by NZ)
C           The B-spline coefficients computed by DB3INK.
C
C   WORK    Double precision 1D array (size KY*KZ+3*max(KX,KY,KZ)+KZ)
C           A working storage array.
C
C***REFERENCES  CARL DE BOOR, A PRACTICAL GUIDE TO SPLINES,
C                 SPRINGER-VERLAG, NEW YORK, 1978.
C***ROUTINES CALLED  DINTRV,DBVALU
C***END PROLOGUE  DB3VAL
C
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C
C   MODIFICATION
C   ------------
C
C   ADDED CHECK TO SEE IF X OR Y IS OUT OF RANGE, IF SO, RETURN 0.0
C
C   R.F. BOISVERT, NIST
C   22 FEB 00
C
C<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
C  ------------
C  DECLARATIONS
C  ------------
C
C  PARAMETERS
C
      INTEGER
     *        IDX, IDY, IDZ, NX, NY, NZ, KX, KY, KZ
      DOUBLE PRECISION
     *     XVAL, YVAL, ZVAL, TX(*), TY(*), TZ(*), BCOEF(NX,NY,NZ),
     *     WORK(*)
C
C  LOCAL VARIABLES
C
      INTEGER
     *        ILOY, ILOZ, INBVX, INBV1, INBV2, LEFTY, LEFTZ, MFLAG,
     *        KCOLY, KCOLZ, IZ, IZM1, IW, I, J, K
      DOUBLE PRECISION
     *     DBVALU
C
      DATA ILOY /1/,  ILOZ /1/,  INBVX /1/
C     SAVE ILOY    ,  ILOZ    ,  INBVX
C
C
C***FIRST EXECUTABLE STATEMENT
      DB3VAL = 0.0D0
C  NEXT STATEMENT - RFB MOD
      IF (XVAL.LT.TX(1) .OR. XVAL.GT.TX(NX+KX) .OR.
     +    YVAL.LT.TY(1) .OR. YVAL.GT.TY(NY+KY) .OR.
     +    ZVAL.LT.TZ(1) .OR. ZVAL.GT.TZ(NZ+KZ)) then
      write(0,*) 'db3.f LINE 1785, NX,NY,NZ: ',NX,NY,NZ
      write(0,*) 'db3.f LINE 1785, KX,KY,KZ: ',KX,KY,KZ
      write(0,*) 'db3.f LINE 1785, X,Y,Z: ',XVAL,YVAL,ZVAL
      write(0,*) '1786, TX(1),TY(1),TZ(1): ',TX(1),TY(1),TZ(1)
      write(0,*) '1787, TX(end),TY(end),TZ(end): ',
     +                  TX(NX+KX),TY(NY+KY),TZ(NZ+KZ)
      stop 'DB3VAL value out of range -- stopping'
      ENDIF
      
      CALL DINTRV(TY,NY+KY,YVAL,ILOY,LEFTY,MFLAG)
      IF (MFLAG .NE. 0)  stop 'problem calling DINTRV TY '
      CALL DINTRV(TZ,NZ+KZ,ZVAL,ILOZ,LEFTZ,MFLAG)
      IF (MFLAG .NE. 0)  stop 'problem calling DINTRV TZ '
         IZ = 1 + KY*KZ
         IW = IZ + KZ
         KCOLZ = LEFTZ - KZ
         I = 0
         DO 50 K=1,KZ
            KCOLZ = KCOLZ + 1
            KCOLY = LEFTY - KY
            DO 50 J=1,KY
               I = I + 1
               KCOLY = KCOLY + 1
               WORK(I) = DBVALU(TX,BCOEF(1,KCOLY,KCOLZ),NX,KX,IDX,XVAL,
     *                           INBVX,WORK(IW))
   50    CONTINUE
         INBV1 = 1
         IZM1 = IZ - 1
         KCOLY = LEFTY - KY + 1
         DO 60 K=1,KZ
            I = (K-1)*KY + 1
            J = IZM1 + K
            WORK(J) = DBVALU(TY(KCOLY),WORK(I),KY,KY,IDY,YVAL,
     *                           INBV1,WORK(IW))
  60     CONTINUE
         INBV2 = 1
         KCOLZ = LEFTZ - KZ + 1
         DB3VAL = DBVALU(TZ(KCOLZ),WORK(IZ),KZ,KZ,IDZ,ZVAL,INBV2,
     *                  WORK(IW))
  100 CONTINUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION DBVALU(T,A,N,K,IDERIV,X,INBV,WORK)
C***BEGIN PROLOGUE  DBVALU
C***DATE WRITTEN   800901   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***CATEGORY NO.  E3,K6
C***KEYWORDS  B-SPLINE,DATA FITTING,DOUBLE PRECISION,INTERPOLATION,
C             SPLINE
C***AUTHOR  AMOS, D. E., (SNLA)
C***PURPOSE  Evaluates the B-representation of a B-spline at X for the
C            function value or any of its derivatives.
C***DESCRIPTION
C
C     Written by Carl de Boor and modified by D. E. Amos
C
C     Reference
C         SIAM J. Numerical Analysis, 14, No. 3, June, 1977, pp.441-472.
C
C     Abstract   **** a double precision routine ****
C         DBVALU is the BVALUE function of the reference.
C
C         DBVALU evaluates the B-representation (T,A,N,K) of a B-spline
C         at X for the function value on IDERIV=0 or any of its
C         derivatives on IDERIV=1,2,...,K-1.  Right limiting values
C         (right derivatives) are returned except at the right end
C         point X=T(N+1) where left limiting values are computed.  The
C         spline is defined on T(K) .LE. X .LE. T(N+1).  DBVALU returns
C         a fatal error message when X is outside of this interval.
C
C         To compute left derivatives or left limiting values at a
C         knot T(I), replace N by I-1 and set X=T(I), I=K+1,N+1.
C
C         DBVALU calls DINTRV
C
C     Description of Arguments
C
C         Input      T,A,X are double precision
C          T       - knot vector of length N+K
C          A       - B-spline coefficient vector of length N
C          N       - number of B-spline coefficients
C                    N = sum of knot multiplicities-K
C          K       - order of the B-spline, K .GE. 1
C          IDERIV  - order of the derivative, 0 .LE. IDERIV .LE. K-1
C                    IDERIV = 0 returns the B-spline value
C          X       - argument, T(K) .LE. X .LE. T(N+1)
C          INBV    - an initialization parameter which must be set
C                    to 1 the first time DBVALU is called.
C
C         Output     WORK,DBVALU are double precision
C          INBV    - INBV contains information for efficient process-
C                    ing after the initial call and INBV must not
C                    be changed by the user.  Distinct splines require
C                    distinct INBV parameters.
C          WORK    - work vector of length 3*K.
C          DBVALU  - value of the IDERIV-th derivative at X
C
C     Error Conditions
C         An improper input is a fatal error
C***REFERENCES  C. DE BOOR, *PACKAGE FOR CALCULATING WITH B-SPLINES*,
C                 SIAM JOURNAL ON NUMERICAL ANALYSIS, VOLUME 14, NO. 3,
C                 JUNE 1977, PP. 441-472.
C***ROUTINES CALLED  DINTRV,XERROR
C***END PROLOGUE  DBVALU
C
C
      INTEGER I,IDERIV,IDERP1,IHI,IHMKMJ,ILO,IMK,IMKPJ, INBV, IPJ,
     1 IP1, IP1MJ, J, JJ, J1, J2, K, KMIDER, KMJ, KM1, KPK, MFLAG, N
      DOUBLE PRECISION A, FKMJ, T, WORK, X
      DIMENSION T(*), A(N), WORK(*)
C***FIRST EXECUTABLE STATEMENT  DBVALU
      DBVALU = 0.0D0
      IF(K.LT.1) GO TO 102
      IF(N.LT.K) GO TO 101
      IF(IDERIV.LT.0 .OR. IDERIV.GE.K) GO TO 110
      KMIDER = K - IDERIV
C
C *** FIND *I* IN (K,N) SUCH THAT T(I) .LE. X .LT. T(I+1)
C     (OR, .LE. T(I+1) IF T(I) .LT. T(I+1) = T(N+1)).
      KM1 = K - 1
      CALL DINTRV(T, N+1, X, INBV, I, MFLAG)
      IF (X.LT.T(K)) GO TO 120
      IF (MFLAG.EQ.0) GO TO 20
      IF (X.GT.T(I)) GO TO 130
   10 IF (I.EQ.K) GO TO 140
      I = I - 1
      IF (X.EQ.T(I)) GO TO 10
C
C *** DIFFERENCE THE COEFFICIENTS *IDERIV* TIMES
C     WORK(I) = AJ(I), WORK(K+I) = DP(I), WORK(K+K+I) = DM(I), I=1.K
C
   20 IMK = I - K
      DO 30 J=1,K
        IMKPJ = IMK + J
        WORK(J) = A(IMKPJ)
   30 CONTINUE
      IF (IDERIV.EQ.0) GO TO 60
      DO 50 J=1,IDERIV
        KMJ = K - J
        FKMJ = DBLE(FLOAT(KMJ))
        DO 40 JJ=1,KMJ
          IHI = I + JJ
          IHMKMJ = IHI - KMJ
          WORK(JJ) = (WORK(JJ+1)-WORK(JJ))/(T(IHI)-T(IHMKMJ))*FKMJ
   40   CONTINUE
   50 CONTINUE
C
C *** COMPUTE VALUE AT *X* IN (T(I),(T(I+1)) OF IDERIV-TH DERIVATIVE,
C     GIVEN ITS RELEVANT B-SPLINE COEFF. IN AJ(1),...,AJ(K-IDERIV).
   60 IF (IDERIV.EQ.KM1) GO TO 100
      IP1 = I + 1
      KPK = K + K
      J1 = K + 1
      J2 = KPK + 1
      DO 70 J=1,KMIDER
        IPJ = I + J
        WORK(J1) = T(IPJ) - X
        IP1MJ = IP1 - J
        WORK(J2) = X - T(IP1MJ)
        J1 = J1 + 1
        J2 = J2 + 1
   70 CONTINUE
      IDERP1 = IDERIV + 1
      DO 90 J=IDERP1,KM1
        KMJ = K - J
        ILO = KMJ
        DO 80 JJ=1,KMJ
          WORK(JJ) = (WORK(JJ+1)*WORK(KPK+ILO)+WORK(JJ)
     1              *WORK(K+JJ))/(WORK(KPK+ILO)+WORK(K+JJ))
          ILO = ILO - 1
   80   CONTINUE
   90 CONTINUE
  100 DBVALU = WORK(1)
      RETURN
C
C
  101 CONTINUE
      CALL XERROR( ' DBVALU,  N DOES NOT SATISFY N.GE.K',35,2,1)
      RETURN
  102 CONTINUE
      CALL XERROR( ' DBVALU,  K DOES NOT SATISFY K.GE.1',35,2,1)
      RETURN
  110 CONTINUE
      CALL XERROR( ' DBVALU,  IDERIV DOES NOT SATISFY 0.LE.IDERIV.LT.K',
     1 50, 2, 1)
      RETURN
  120 CONTINUE
      CALL XERROR( ' DBVALU,  X IS N0T GREATER THAN OR EQUAL TO T(K)',
     1 48, 2, 1)
      RETURN
  130 CONTINUE
      CALL XERROR( ' DBVALU,  X IS NOT LESS THAN OR EQUAL TO T(N+1)',
     1 47, 2, 1)
      RETURN
  140 CONTINUE
      CALL XERROR( ' DBVALU,  A LEFT LIMITING VALUE CANN0T BE OBTAINED A
     1T T(K)',    58, 2, 1)
      RETURN
      END
      SUBROUTINE DINTRV(XT,LXT,X,ILO,ILEFT,MFLAG)
C***BEGIN PROLOGUE  DINTRV
C***DATE WRITTEN   800901   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  E3,K6
C***KEYWORDS  B-SPLINE,DATA FITTING,DOUBLE PRECISION,INTERPOLATION,
C             SPLINE
C***AUTHOR  AMOS, D. E., (SNLA)
C***PURPOSE  Computes the largest integer ILEFT in 1.LE.ILEFT.LE.LXT
C            such that XT(ILEFT).LE.X where XT(*) is a subdivision of
C            the X interval.
C***DESCRIPTION
C
C     Written by Carl de Boor and modified by D. E. Amos
C
C     Reference
C         SIAM J.  Numerical Analysis, 14, No. 3, June 1977, pp.441-472.
C
C     Abstract    **** a double precision routine ****
C         DINTRV is the INTERV routine of the reference.
C
C         DINTRV computes the largest integer ILEFT in 1 .LE. ILEFT .LE.
C         LXT such that XT(ILEFT) .LE. X where XT(*) is a subdivision of
C         the X interval.  Precisely,
C
C                      X .LT. XT(1)                1         -1
C         if  XT(I) .LE. X .LT. XT(I+1)  then  ILEFT=I  , MFLAG=0
C           XT(LXT) .LE. X                         LXT        1,
C
C         That is, when multiplicities are present in the break point
C         to the left of X, the largest index is taken for ILEFT.
C
C     Description of Arguments
C
C         Input      XT,X are double precision
C          XT      - XT is a knot or break point vector of length LXT
C          LXT     - length of the XT vector
C          X       - argument
C          ILO     - an initialization parameter which must be set
C                    to 1 the first time the spline array XT is
C                    processed by DINTRV.
C
C         Output
C          ILO     - ILO contains information for efficient process-
C                    ing after the initial call and ILO must not be
C                    changed by the user.  Distinct splines require
C                    distinct ILO parameters.
C          ILEFT   - largest integer satisfying XT(ILEFT) .LE. X
C          MFLAG   - signals when X lies out of bounds
C
C     Error Conditions
C         None
C***REFERENCES  C. DE BOOR, *PACKAGE FOR CALCULATING WITH B-SPLINES*,
C                 SIAM JOURNAL ON NUMERICAL ANALYSIS, VOLUME 14, NO. 3,
C                 JUNE 1977, PP. 441-472.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DINTRV
C
C
      INTEGER IHI, ILEFT, ILO, ISTEP, LXT, MFLAG, MIDDLE
      DOUBLE PRECISION X, XT
      DIMENSION XT(LXT)
C***FIRST EXECUTABLE STATEMENT  DINTRV
      IHI = ILO + 1
      IF (IHI.LT.LXT) GO TO 10
      IF (X.GE.XT(LXT)) GO TO 110
      IF (LXT.LE.1) GO TO 90
      ILO = LXT - 1
      IHI = LXT
C
   10 IF (X.GE.XT(IHI)) GO TO 40
      IF (X.GE.XT(ILO)) GO TO 100
C
C *** NOW X .LT. XT(IHI) . FIND LOWER BOUND
      ISTEP = 1
   20 IHI = ILO
      ILO = IHI - ISTEP
      IF (ILO.LE.1) GO TO 30
      IF (X.GE.XT(ILO)) GO TO 70
      ISTEP = ISTEP*2
      GO TO 20
   30 ILO = 1
      IF (X.LT.XT(1)) GO TO 90
      GO TO 70
C *** NOW X .GE. XT(ILO) . FIND UPPER BOUND
   40 ISTEP = 1
   50 ILO = IHI
      IHI = ILO + ISTEP
      IF (IHI.GE.LXT) GO TO 60
      IF (X.LT.XT(IHI)) GO TO 70
      ISTEP = ISTEP*2
      GO TO 50
   60 IF (X.GE.XT(LXT)) GO TO 110
      IHI = LXT
C
C *** NOW XT(ILO) .LE. X .LT. XT(IHI) . NARROW THE INTERVAL
   70 MIDDLE = (ILO+IHI)/2
      IF (MIDDLE.EQ.ILO) GO TO 100
C     NOTE. IT IS ASSUMED THAT MIDDLE = ILO IN CASE IHI = ILO+1
      IF (X.LT.XT(MIDDLE)) GO TO 80
      ILO = MIDDLE
      GO TO 70
   80 IHI = MIDDLE
      GO TO 70
C *** SET OUTPUT AND RETURN
   90 MFLAG = -1
      ILEFT = 1
      RETURN
  100 MFLAG = 0
      ILEFT = ILO
      RETURN
  110 MFLAG = 1
      ILEFT = LXT
      RETURN
      END


      SUBROUTINE DBINT4(X,Y,NDATA,IBCL,IBCR,FBCL,FBCR,KNTOPT,T,BCOEF,N,
     1   K,W)
C***BEGIN PROLOGUE  DBINT4
C***DATE WRITTEN   800901   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***CATEGORY NO.  E1A
C***KEYWORDS  B-SPLINE,DATA FITTING,DOUBLE PRECISION,INTERPOLATION,
C             SPLINE
C***AUTHOR  AMOS, D. E., (SNLA)
C***PURPOSE  Computes the B-representation of a cubic spline
C            which interpolates data (X(I),Y(I)),I=1,NDATA.
C***DESCRIPTION
C
C     Written by D. E. Amos, August, 1979.
C
C     References
C         SAND-78-1968
C
C         A Practical Guide to Splines by C. de Boor, Applied
C         Mathematics Series 27, Springer, 1979.
C
C         SIAM J. Numerical Analysis, 14, No. 3, June, 1977, pp.441-472.
C
C     Abstract    **** a double precision routine ****
C
C         DBINT4 computes the B representation (T,BCOEF,N,K) of a
C         cubic spline (K=4) which interpolates data (X(I),Y(I)),
C         I=1,NDATA.  Parameters IBCL, IBCR, FBCL, FBCR allow the
C         specification of the spline first or second derivative at
C         both X(1) and X(NDATA).  When this data is not specified
C         by the problem, it is common practice to use a natural
C         spline by setting second derivatives at X(1) and X(NDATA)
C         to zero (IBCL=IBCR=2,FBCL=FBCR=0.0).  The spline is defined
C         on T(4) .LE. X .LE. T(N+1) with (ordered) interior knots at
C         X(I) values where N=NDATA+2.  The knots T(1),T(2),T(3) lie to
C         the left of T(4)=X(1) and the knots T(N+2), T(N+3), T(N+4)
C         lie to the right of T(N+1)=X(NDATA) in increasing order.  If
C         no extrapolation outside (X(1),X(NDATA)) is anticipated, the
C         knots T(1)=T(2)=T(3)=T(4)=X(1) and T(N+2)=T(N+3)=T(N+4)=
C         T(N+1)=X(NDATA) can be specified by KNTOPT=1.  KNTOPT=2
C         selects a knot placement for T(1), T(2), T(3) to make the
C         first 7 knots symmetric about T(4)=X(1) and similarly for
C         T(N+2), T(N+3), T(N+4) about T(N+1)=X(NDATA).  KNTOPT=3
C         allows the user to make his own selection, in increasing
C         order, for T(1), T(2), T(3) to the left of X(1) and T(N+2),
C         T(N+3), T(N+4) to the right of X(NDATA) in the work array
C         W(1) through W(6).  In any case, the interpolation on
C         T(4) .LE. X .LE. T(N+1) by using function DBVALU is unique
C         for given boundary conditions.
C
C         DBINT4 calls DBSPVD, DBNFAC, DBNSLV, D1MACH, XERROR
C
C     Description of Arguments
C
C         Input      X,Y,FBCL,FBCR,W are double precision
C           X      - X vector of abscissae of length NDATA, distinct
C                    and in increasing order
C           Y      - Y vector of ordinates of length NDATA
C           NDATA  - number of data points, NDATA .GE. 2
C           IBCL   - selection parameter for left boundary condition
C                    IBCL = 1 constrain the first derivative at
C                             X(1) to FBCL
C                         = 2 constrain the second derivative at
C                             X(1) to FBCL
C           IBCR   - selection parameter for right boundary condition
C                    IBCR = 1 constrain first derivative at
C                             X(NDATA) to FBCR
C                    IBCR = 2 constrain second derivative at
C                             X(NDATA) to FBCR
C           FBCL   - left boundary values governed by IBCL
C           FBCR   - right boundary values governed by IBCR
C           KNTOPT - knot selection parameter
C                    KNTOPT = 1 sets knot multiplicity at T(4) and
C                               T(N+1) to 4
C                           = 2 sets a symmetric placement of knots
C                               about T(4) and T(N+1)
C                           = 3 sets T(I)=W(I) and T(N+1+I)=W(3+I),I=1,3
C                               where W(I),I=1,6 is supplied by the user
C           W      - work array of dimension at least 5*(NDATA+2)
C                    If KNTOPT=3, then W(1),W(2),W(3) are knot values to
C                    the left of X(1) and W(4),W(5),W(6) are knot
C                    values to the right of X(NDATA) in increasing
C                    order to be supplied by the user
C
C         Output     T,BCOEF are double precision
C           T      - knot array of length N+4
C           BCOEF  - B spline coefficient array of length N
C           N      - number of coefficients, N=NDATA+2
C           K      - order of spline, K=4
C
C     Error Conditions
C         Improper  input is a fatal error
C         Singular system of equations is a fatal error
C***REFERENCES  D.E. AMOS, *COMPUTATION WITH SPLINES AND B-SPLINES*,
C                 SAND78-1968, SANDIA LABORATORIES, MARCH 1979.
C               C. DE BOOR, *PACKAGE FOR CALCULATING WITH B-SPLINES*,
C                 SIAM JOURNAL ON NUMERICAL ANALYSIS, VOLUME 14, NO. 3,
C                 JUNE 1977, PP. 441-472.
C               C. DE BOOR, *A PRACTICAL GUIDE TO SPLINES*, APPLIED
C                 MATHEMATICS SERIES 27, SPRINGER, 1979.
C***ROUTINES CALLED  D1MACH,DBNFAC,DBNSLV,DBSPVD,XERROR
C***END PROLOGUE  DBINT4
C
C
      INTEGER I, IBCL, IBCR, IFLAG, ILB, ILEFT, IT, IUB, IW, IWP, J,
     1 JW, K, KNTOPT, N, NDATA, NDM, NP, NWROW
      DOUBLE PRECISION BCOEF,FBCL,FBCR,T,TOL,TXN,TX1,VNIKX,W,WDTOL,
     1 WORK,X,XL,Y
      DOUBLE PRECISION D1MACH
      DIMENSION X(*), Y(*), T(*), BCOEF(*), W(5,*), VNIKX(4,4), WORK(15)
C***FIRST EXECUTABLE STATEMENT  DBINT4
      WDTOL = D1MACH(4)
      TOL = DSQRT(WDTOL)
      IF (NDATA.LT.2) GO TO 200
      NDM = NDATA - 1
      DO 10 I=1,NDM
        IF (X(I).GE.X(I+1)) GO TO 210
   10 CONTINUE
      IF (IBCL.LT.1 .OR. IBCL.GT.2) GO TO 220
      IF (IBCR.LT.1 .OR. IBCR.GT.2) GO TO 230
      IF (KNTOPT.LT.1 .OR. KNTOPT.GT.3) GO TO 240
      K = 4
      N = NDATA + 2
      NP = N + 1
      DO 20 I=1,NDATA
        T(I+3) = X(I)
   20 CONTINUE
      GO TO (30, 50, 90), KNTOPT
C     SET UP KNOT ARRAY WITH MULTIPLICITY 4 AT X(1) AND X(NDATA)
   30 CONTINUE
      DO 40 I=1,3
        T(4-I) = X(1)
        T(NP+I) = X(NDATA)
   40 CONTINUE
      GO TO 110
C     SET UP KNOT ARRAY WITH SYMMETRIC PLACEMENT ABOUT END POINTS
   50 CONTINUE
      IF (NDATA.GT.3) GO TO 70
      XL = (X(NDATA)-X(1))/3.0D0
      DO 60 I=1,3
        T(4-I) = T(5-I) - XL
        T(NP+I) = T(NP+I-1) + XL
   60 CONTINUE
      GO TO 110
   70 CONTINUE
      TX1 = X(1) + X(1)
      TXN = X(NDATA) + X(NDATA)
      DO 80 I=1,3
        T(4-I) = TX1 - X(I+1)
        T(NP+I) = TXN - X(NDATA-I)
   80 CONTINUE
      GO TO 110
C     SET UP KNOT ARRAY LESS THAN X(1) AND GREATER THAN X(NDATA) TO BE
C     SUPPLIED BY USER IN WORK LOCATIONS W(1) THROUGH W(6) WHEN KNTOPT=3
   90 CONTINUE
      DO 100 I=1,3
        T(4-I) = W(4-I,1)
        JW = MAX0(1,I-1)
        IW = MOD(I+2,5)+1
        T(NP+I) = W(IW,JW)
        IF (T(4-I).GT.T(5-I)) GO TO 250
        IF (T(NP+I).LT.T(NP+I-1)) GO TO 250
  100 CONTINUE
  110 CONTINUE
C
      DO 130 I=1,5
        DO 120 J=1,N
          W(I,J) = 0.0D0
  120   CONTINUE
  130 CONTINUE
C     SET UP LEFT INTERPOLATION POINT AND LEFT BOUNDARY CONDITION FOR
C     RIGHT LIMITS
      IT = IBCL + 1
      CALL DBSPVD(T, K, IT, X(1), K, 4, VNIKX, WORK)
      IW = 0
      IF (DABS(VNIKX(3,1)).LT.TOL) IW = 1
      DO 140 J=1,3
        W(J+1,4-J) = VNIKX(4-J,IT)
        W(J,4-J) = VNIKX(4-J,1)
  140 CONTINUE
      BCOEF(1) = Y(1)
      BCOEF(2) = FBCL
C     SET UP INTERPOLATION EQUATIONS FOR POINTS I=2 TO I=NDATA-1
      ILEFT = 4
      IF (NDM.LT.2) GO TO 170
      DO 160 I=2,NDM
        ILEFT = ILEFT + 1
        CALL DBSPVD(T, K, 1, X(I), ILEFT, 4, VNIKX, WORK)
        DO 150 J=1,3
          W(J+1,3+I-J) = VNIKX(4-J,1)
  150   CONTINUE
        BCOEF(I+1) = Y(I)
  160 CONTINUE
C     SET UP RIGHT INTERPOLATION POINT AND RIGHT BOUNDARY CONDITION FOR
C     LEFT LIMITS(ILEFT IS ASSOCIATED WITH T(N)=X(NDATA-1))
  170 CONTINUE
      IT = IBCR + 1
      CALL DBSPVD(T, K, IT, X(NDATA), ILEFT, 4, VNIKX, WORK)
      JW = 0
      IF (DABS(VNIKX(2,1)).LT.TOL) JW = 1
      DO 180 J=1,3
        W(J+1,3+NDATA-J) = VNIKX(5-J,IT)
        W(J+2,3+NDATA-J) = VNIKX(5-J,1)
  180 CONTINUE
      BCOEF(N-1) = FBCR
      BCOEF(N) = Y(NDATA)
C     SOLVE SYSTEM OF EQUATIONS
      ILB = 2 - JW
      IUB = 2 - IW
      NWROW = 5
      IWP = IW + 1
      CALL DBNFAC(W(IWP,1), NWROW, N, ILB, IUB, IFLAG)
      IF (IFLAG.EQ.2) GO TO 190
      CALL DBNSLV(W(IWP,1), NWROW, N, ILB, IUB, BCOEF)
      RETURN
C
C
  190 CONTINUE
      CALL XERROR( ' DBINT4,  THE SYSTEM OF EQUATIONS IS SINGULAR',
     1 45, 2, 1)
      RETURN
  200 CONTINUE
      CALL XERROR( ' DBINT4,  NDATA IS LESS THAN 2', 30, 2, 1)
      RETURN
  210 CONTINUE
      CALL XERROR( ' DBINT4,  X VALUES ARE NOT DISTINCT OR NOT ORDERED',
     1 50, 2, 1)
      RETURN
  220 CONTINUE
      CALL XERROR( ' DBINT4,  IBCL IS NOT 1 OR 2', 28, 2, 1)
      RETURN
  230 CONTINUE
      CALL XERROR( ' DBINT4,  IBCR IS NOT 1 OR 2', 28, 2, 1)
      RETURN
  240 CONTINUE
      CALL XERROR( ' DBINT4,  KNTOPT IS NOT 1, 2, OR 3', 34, 2, 1)
      RETURN
  250 CONTINUE
      CALL XERROR(  ' DBINT4,  KNOT INPUT THROUGH W ARRAY IS NOT ORDERED
     1 PROPERLY',  60, 2, 1)
      RETURN
      END



      SUBROUTINE DBSPVD(T,K,NDERIV,X,ILEFT,LDVNIK,VNIKX,WORK)
C***BEGIN PROLOGUE  DBSPVD
C***DATE WRITTEN   800901   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***CATEGORY NO.  E3,K6
C***KEYWORDS  B-SPLINE,DATA FITTING,DOUBLE PRECISION,INTERPOLATION,
C             SPLINE
C***AUTHOR  AMOS, D. E., (SNLA)
C***PURPOSE  Calculates the value and all derivatives of order less than
C            NDERIV of all basis functions which do not vanish at X.
C***DESCRIPTION
C
C     Written by Carl de Boor and modified by D. E. Amos
C
C     Reference
C         SIAM J. Numerical Analysis, 14, No. 3, June, 1977, pp.441-472.
C
C     Abstract    **** a double precision routine ****
C
C         DBSPVD is the BSPLVD routine of the reference.
C
C         DBSPVD calculates the value and all derivatives of order
C         less than NDERIV of all basis functions which do not
C         (possibly) vanish at X.  ILEFT is input such that
C         T(ILEFT) .LE. X .LT. T(ILEFT+1).  A call to INTRV(T,N+1,X,
C         ILO,ILEFT,MFLAG) will produce the proper ILEFT.  The output of
C         DBSPVD is a matrix VNIKX(I,J) of dimension at least (K,NDERIV)
C         whose columns contain the K nonzero basis functions and
C         their NDERIV-1 right derivatives at X, I=1,K, J=1,NDERIV.
C         These basis functions have indices ILEFT-K+I, I=1,K,
C         K .LE. ILEFT .LE. N.  The nonzero part of the I-th basis
C         function lies in (T(I),T(I+K)), I=1,N).
C
C         If X=T(ILEFT+1) then VNIKX contains left limiting values
C         (left derivatives) at T(ILEFT+1).  In particular, ILEFT = N
C         produces left limiting values at the right end point
C         X=T(N+1).  To obtain left limiting values at T(I), I=K+1,N+1,
C         set X= next lower distinct knot, call INTRV to get ILEFT,
C         set X=T(I), and then call DBSPVD.
C
C         DBSPVD calls DBSPVN
C
C     Description of Arguments
C         Input      T,X are double precision
C          T       - knot vector of length N+K, where
C                    N = number of B-spline basis functions
C                    N = sum of knot multiplicities-K
C          K       - order of the B-spline, K .GE. 1
C          NDERIV  - number of derivatives = NDERIV-1,
C                    1 .LE. NDERIV .LE. K
C          X       - argument of basis functions,
C                    T(K) .LE. X .LE. T(N+1)
C          ILEFT   - largest integer such that
C                    T(ILEFT) .LE. X .LT.  T(ILEFT+1)
C          LDVNIK  - leading dimension of matrix VNIKX
C
C         Output     VNIKX,WORK are double precision
C          VNIKX   - matrix of dimension at least (K,NDERIV) contain-
C                    ing the nonzero basis functions at X and their
C                    derivatives columnwise.
C          WORK    - a work vector of length (K+1)*(K+2)/2
C
C     Error Conditions
C         Improper input is a fatal error
C***REFERENCES  C. DE BOOR, *PACKAGE FOR CALCULATING WITH B-SPLINES*,
C                 SIAM JOURNAL ON NUMERICAL ANALYSIS, VOLUME 14, NO. 3,
C                 JUNE 1977, PP. 441-472.
C***ROUTINES CALLED  DBSPVN,XERROR
C***END PROLOGUE  DBSPVD
C
C
      INTEGER I,IDERIV,ILEFT,IPKMD,J,JJ,JLOW,JM,JP1MID,K,KMD, KP1, L,
     1 LDUMMY, M, MHIGH, NDERIV
      DOUBLE PRECISION FACTOR, FKMD, T, V, VNIKX, WORK, X
C     DIMENSION T(ILEFT+K), WORK((K+1)*(K+2)/2)
C     A(I,J) = WORK(I+J*(J+1)/2),  I=1,J+1  J=1,K-1
C     A(I,K) = W0RK(I+K*(K-1)/2)  I=1.K
C     WORK(1) AND WORK((K+1)*(K+2)/2) ARE NOT USED.
      DIMENSION T(*), VNIKX(LDVNIK,NDERIV), WORK(*)
C***FIRST EXECUTABLE STATEMENT  DBSPVD
      IF(K.LT.1) GO TO 200
      IF(NDERIV.LT.1 .OR. NDERIV.GT.K) GO TO 205
      IF(LDVNIK.LT.K) GO TO 210
      IDERIV = NDERIV
      KP1 = K + 1
      JJ = KP1 - IDERIV
      CALL DBSPVN(T, JJ, K, 1, X, ILEFT, VNIKX, WORK, IWORK)
      IF (IDERIV.EQ.1) GO TO 100
      MHIGH = IDERIV
      DO 20 M=2,MHIGH
        JP1MID = 1
        DO 10 J=IDERIV,K
          VNIKX(J,IDERIV) = VNIKX(JP1MID,1)
          JP1MID = JP1MID + 1
   10   CONTINUE
        IDERIV = IDERIV - 1
        JJ = KP1 - IDERIV
        CALL DBSPVN(T, JJ, K, 2, X, ILEFT, VNIKX, WORK, IWORK)
   20 CONTINUE
C
      JM = KP1*(KP1+1)/2
      DO 30 L = 1,JM
        WORK(L) = 0.0D0
   30 CONTINUE
C     A(I,I) = WORK(I*(I+3)/2) = 1.0       I = 1,K
      L = 2
      J = 0
      DO 40 I = 1,K
        J = J + L
        WORK(J) = 1.0D0
        L = L + 1
   40 CONTINUE
      KMD = K
      DO 90 M=2,MHIGH
        KMD = KMD - 1
        FKMD = FLOAT(KMD)
        I = ILEFT
        J = K
        JJ = J*(J+1)/2
        JM = JJ - J
        DO 60 LDUMMY=1,KMD
          IPKMD = I + KMD
          FACTOR = FKMD/(T(IPKMD)-T(I))
          DO 50 L=1,J
            WORK(L+JJ) = (WORK(L+JJ)-WORK(L+JM))*FACTOR
   50     CONTINUE
          I = I - 1
          J = J - 1
          JJ = JM
          JM = JM - J
   60   CONTINUE
C
        DO 80 I=1,K
          V = 0.0D0
          JLOW = MAX0(I,M)
          JJ = JLOW*(JLOW+1)/2
          DO 70 J=JLOW,K
            V = WORK(I+JJ)*VNIKX(J,M) + V
            JJ = JJ + J + 1
   70     CONTINUE
          VNIKX(I,M) = V
   80   CONTINUE
   90 CONTINUE
  100 RETURN
C
C
  200 CONTINUE
      CALL XERROR( ' DBSPVD,  K DOES NOT SATISFY K.GE.1', 35, 2, 1)
      RETURN
  205 CONTINUE
      CALL XERROR( ' DBSPVD,  NDERIV DOES NOT SATISFY 1.LE.NDERIV.LE.K',
     1 50, 2, 1)
      RETURN
  210 CONTINUE
      CALL XERROR( ' DBSPVD,  LDVNIK DOES NOT SATISFY LDVNIK.GE.K',45,
     1 2, 1)
      RETURN
      END

