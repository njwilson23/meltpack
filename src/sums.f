c****************************************************************************
C $Id: sums.f,v 1.1 2000-01-26 14:12:55-07 braup Exp $
C
C $Log: sums.f,v $
C Revision 1.1  2000-01-26 14:12:55-07  braup
C Initial revision
C
C
C*    Subroutine sums
C
C     Purpose:
C	Compute the sums of terms for the matrix elements and column
C	vector used for determining the best-fit surface to the
C	cross-correlation values in the neighborhood of the peak.
C
C     Arguments:
C	Name	Type	   I/O	Description
C	----	----	   ---	-----------
C	CPVAL	Real Array  I	5 by 5 array of cross-correlation
C				values, in units of standard deviations
C				above background, centered on
C				correlation peak.
C	MFIT	Int.	    I	Method of fitting surface 
C				 1 - Elliptical paraboloid
C				 2 - Elliptical Gaussian
C				 3 - Reciprocal Paraboloid
C	Z	Real Array  O	Function of cross-correlation value to
C				which quadratic surface is to be fit
C	WGHTS	Real Array  O	Weights to be assigned each point in 5
C				by 5 array of values.
C	B	Real Array  O	Matrix of sums of polynomial terms 
C	VECTOR	Real Array  O	Column vector obtained by summing
C				products of polynomial terms and an
C				appropriate function of cross-
C				correlation value.
C
C     Local Variables and Arrays:
C	Name	Type	    Description
C	----	----	    -----------
C	I	Integer	    Loop index:  matrix row.
C	IC	Integer	    Loop index:  column of peak-neighborhood
C			    array.
C	IR	Integer	    Loop index:  row of peak-neighborhood array.
C	IVALPT	Integer	    Pointer into arrays of cross-correlation
C			    values, functions of values, and weights in
C			    neighborhood of peak.
C	J	Integer	    Loop index:  matrix column.
C	TERM	Real Array  Terms in approximating quadratic polynomial.
C	VAL	Real	    Value of one element in CPVAL array,
C			    increased to 1.0 standard deviations above
C			    background if necessary.
C
C     Program History:
C	Version	  Date	      Author	Code/Contr.	Change Request
C	-------	  ----	      ------	-----------	--------------
C	5.0       Jan 1990    D. Steinwand   EDC 	LAS 5.0 Conversion
C         2       SEPT 83     R. White       CSC        SCN 20-66
C	001	  June 1983   R. White	     CSC	N/A
C
C*    Language:  Fortran 77.
c
C**   ALGORITHM DESCRIPTION
C
C	    Clear array of matrix elements
C	    Clear elements in column vector 
C	    DO-FOR each value in correlation peak array
C		IF value less than 1 standard deviation above background
C		 THEN
C		    Replace value by 1.0
C		END-IF
C		IF elliptical paraboloid specified THEN
C		    Use cross-correlation value for z 
C		    Set weight to one
C		ELSE-IF elliptical Gaussian specified THEN
C		    Use negative logatithm of correlation value as z
C		    Use square of correlation value as weight
C		ELSE (reciprocal paraboloid)
C		    Use reciprocal of correlation value as z
C		    Use fourth power of correlation value as weight
C		END-IF
C		DO-FOR each of terms 1, x, y, x**2, x*y, y**2
C		    Add z times term to column-vector entry
C		    DO-FOR each of terms
C			Add product of terms to matrix entry
C		    END-DO
C		END-DO
C	    END-DO
C**         Return
c
c****************************************************************************
      SUBROUTINE sums(cpval, mfit, z, wghts, b, vector)

      integer	mfit
      real	b(6,6),		cpval(25),	vector(6),
     .		wghts(25),	z(25)
 
C  Specification of Local Variables and Arrays
C---------------------------------------------
      integer	i,	ic,	ir,	ivalpt,		j
      real	term(6)
      real	val
 
C  Initialize arrays and constants
C---------------------------------
      do 100 i=1,6
	  do 200 j=1,6
	      b(i,j) = 0.0
 200      continue
	  vector(i) = 0.0
 100  continue
      term(1) = 1.0
      ivalpt = 0
 
C  Compute function of correlation coefficient and assign weight
C  for each location in neighborhood of peak
C-------------------------------------------
      do 300 ir=1,5
	  do 400 ic=1,5
	      ivalpt = ivalpt + 1
	      val = amax1(cpval(ivalpt), 1.0)
	      if (mfit.eq.1) then
		  z(ivalpt) = val
		  wghts(ivalpt) = 1.0
	      else if (mfit.eq.2) then
		  z(ivalpt) = - alog(val)
		  wghts(ivalpt) = val**2
	      else
		  z(ivalpt) = 1.0/val
		  wghts(ivalpt) = val**4
	      end if
 
C  Compute matrix and vector elements in linear equations
C  for best-fit coefficients
C---------------------------
	      term(2) = ic - 3
	      term(3) = ir - 3
	      term(4) = term(2)*term(2)
	      term(5) = term(2)*term(3)
	      term(6) = term(3)*term(3)
	      do 500 i=1,6
		  vector(i) = vector(i)
     &				+ wghts(ivalpt)*term(i)*z(ivalpt)
		  do 600 j=1,6
		      b(i,j) = b(i,j) + wghts(ivalpt)*term(i)*term(j)
 600              continue
 500          continue
 400      continue
 300  continue
      return
      end
