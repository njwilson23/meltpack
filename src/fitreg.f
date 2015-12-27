c****************************************************************************
C $Id: fitreg.f,v 1.1 2000-01-26 14:12:55-07 braup Exp $
C
C $Log: fitreg.f,v $
C Revision 1.1  2000-01-26 14:12:55-07  braup
C Initial revision
C
C
C*    Subroutine fitreg
C
C     Purpose:
C	Fit a quadratic surface to the neighborhood of the correlation
C	peak and from it determine the best-fit registration offsets and
C	their estimated errors.
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
C	PKOFFS	Real Array  O	Best-fit horizontal (pixel) and vertical
C				(line) offsets of correlation peak
C				relative to center of 5 by 5 array.
C	TLERRS	Real Array  O	Estimated horizontal error, vertical
C				error, and h-v cross term in best-fit
C				offsets.
C
C     Local Variables and Arrays:
C	Name	Type	    Description
C	----	----	    -----------
C	B	Real Array  Matrix of sums of polynomial terms.
C	COEFFS	Real Array  Coefficients of best-fit quadratic surface
C			    to values in neighborhood of correlation
C			    peak.
C	DENOM	Real	    Denominator of equation for best-fit
C			    offsets.
C	IERROR	Integer	    Status return from matrix-inversion routine.
C	VECTOR	Real Array  Column vector forming right-hand side of
C			    equation to be solved for coefficients.
C	WGHTS	Real Array  Weight assigned to each value in
C			    neighborhood of peak.
C	WORK	Real Array  Working buffer required by matrix inversion
C			    routine.
C	Z	Real Array  Values of the function of cross-correlation
C			    coefficient to which surface is to be fit.
C
C     Program History:
C	Version	  Date	      Author	Code/Contr.	Change Request
C	-------	  ----	      ------	-----------	--------------
C	5.0	  Jan 1990    D. Steinwand  EDC		LAS 5.0 Conversion
C         2       Sept 83     R. White      CSC         SCN 20-66
C	001	  June 1983   R. White	    CSC		N/A
C
C*    Language:  Fortran 77.
C
C**   ALGORITHM DESCRIPTION
C
C	    Compute elements of matrix and column vector
C	    Invert matrix
C	    Compute coefficients of best-fit quadratic surface 
C	    Compute best-fit peak offsets 
C	    Estimate errors in best-fit offsets	    
C**         Return
C
c****************************************************************************
      SUBROUTINE fitreg(cpval, mfit, pkoffs, tlerrs)
C
C  Specification of Subroutine Arguments
C---------------------------------------
      integer	mfit
      real	cpval(1),	pkoffs(2),	tlerrs(3)
C
C  Specification of Local Variables and Arrays
C---------------------------------------------
      integer	i, j
      real	b(6,6),		coeffs(6),
     .		vector(6),	wghts(25),	work(12),
     .		z(25)
      real	denom
C
C  Compute elements of matrix and column vector
C----------------------------------------------
      call sums(cpval, mfit, z, wghts, b, vector)
C
C  Invert matrix to get best-fit peak offsets
C--------------------------------------------
      call kvert(b, 6, 6, work)
      do 100 i=1,6
	  coeffs(i) = 0.0
	  do 200 J=1,6
	      coeffs(i) = coeffs(i) + b(i,j)*vector(j)
 200      continue
 100  continue
      denom = 4.0*coeffs(4)*coeffs(6) - coeffs(5)*coeffs(5)
      pkoffs(1) = (coeffs(3)*coeffs(5) - 2.0*coeffs(2)*coeffs(6))/denom
      pkoffs(2) = (coeffs(2)*coeffs(5) - 2.0*coeffs(3)*coeffs(4))/denom
C
C  Estimate errors in best-fit offsets
C-------------------------------------
      call esterr(z, wghts, b, coeffs, pkoffs, tlerrs)
C
      return
      end
