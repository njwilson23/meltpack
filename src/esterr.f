c****************************************************************************
C $Id: esterr.f,v 1.1 2000-01-26 14:12:55-07 braup Exp $
C
C $Log: esterr.f,v $
C Revision 1.1  2000-01-26 14:12:55-07  braup
C Initial revision
C
C
C*    Subroutine esterr
C
C     Purpose:
C	Estimate the probable errors in the best-fit location of the
C	correlation peak.
C
C     Arguments:
C	Name	Type	   I/O	Description
C	----	----	   ---	-----------
C	Z	Real Array  I	Function of cross-correlation value to
C				which quadratic surface is to be fit
C	WGHTS	Real Array  I	Weights to be assigned each point in 5
C				by 5 array of values.
C	BNVRS	Real Array  I	Inverse of matrix B.
C	COEFFS	Real Array  I	Coefficients of best-fit quadratic
C				surface to values of Z.
C	PKOFFS	Real Array  O	Best-fit horizontal (pixel) and vertical
C				(line) offsets of coordinate peak
C				relative to center of 5 by 5 array.
C	TLERRS	Real Array  O	Estimated horizontal error, vertical
C				error, and h-v cross term in best-fit
C				offsets.
C
C     Local Variables and Arrays:
C	Name	Type	    Description
C	----	----	    -----------
C	C	Real	    Constant for computing weights or errors.
C	DENOM	Real	    Denominator used for partial derivatives.
C	DU	Real Array  Partial derivatives of horizontal best-fit
C			    peak offset with respect to coefficients of
C			    polynomial surface.
C	DV	Real Array  Partial derivatives of vertical offset.
C	F	Real	    Value of function of correlation coefficient
C			    computed using best-fit coefficients.
C	I	Integer	    Loop index.
C	IVALPT	Integer	    Pointer into arrays of values determined
C			    from measured correlation coefficients.
C	J	Integer	    Loop index.
C	USUM	Real	    Accumulator for sum of terms to determine
C			    estimated errors in horizontal offset.
C	VAR	Real	    Residual variance of fit.
C	VSUM	Real	    Sum of terms to estimate vertical offset
C			    error.
C	X	Real	    Loop index:  horizontal offset from center
C			    of neighborhood of peak correlation values.
C	XSUM	Real	    Accumlator for sum of cross terms to
C			    estimate errors when search chip was
C			    rotated.
C	Y	Real	    Loop index:  vertical offset from peak.
C
C     Program History:
C	Version	  Date	      Author	Code/Contr.	Change Request
C	-------	  ----	      ------	-----------	--------------
C       5.0       Jan 1990    D. Steinwand  EDC         LAS 5.0 Conversion
C         2       SEpt 83     R. White      CSC         SCN 20-66
C	001	  June 1983   R. White	    CSC		N/A
C
C*    Language:  Fortran 77.
C
C**   ALGORITHM DESCRIPTION
C
C	    Compute residual variance for elements of peak array 
C	    Compute value of constant needed for error estimation 
C	    Compute partial derivatives of peak offsets 
C	    Compute estimated errors in offsets and cross term needed to 
C    		recompute errors if search subimage is rotated.
C**    Return 
C
c****************************************************************************
      SUBROUTINE esterr(z, wghts, bnvrs, coeffs, pkoffs, tlerrs)
C
C  Specification of Subroutine Arguments
C---------------------------------------
      real	pkoffs(2),	bnvrs(6,6),	coeffs(6),
     .		tlerrs(3),	wghts(25),	z(25)
C
C  Specification of Local Variables and Arrays
C---------------------------------------------
      integer	i,	ivalpt,		j
      real	du(6),	dv(6)
      real	c,	denom,	f,	usum,	var,	vsum,	x,
     .		xsum,	y
C
C  Compute residual variance for elements of peak array
C------------------------------------------------------
      ivalpt = 1
      var = 0.0
      do 100 y=-2.0,2.0,1.0
	  do 200 x=-2.0,2.0,1.0
	      f = coeffs(1) + coeffs(2)*x + coeffs(3)*y
     .		   + coeffs(4)*x*x + coeffs(5)*x*y + coeffs(6)*y*y
	      var = var + wghts(ivalpt)*(f - z(ivalpt))**2
	      ivalpt = ivalpt + 1
 200      continue
 100  continue
C
C  Compute constant to use with weights
C--------------------------------------
      c = var/19.0
C
C  Compute partial derivatives of peak offsets with respect to
C  polynomial coefficients
C-------------------------
      denom = 4.0*coeffs(4)*coeffs(6) - coeffs(5)*coeffs(5)
      du(1) = 0.0
      du(2) = -2.0*coeffs(6)/denom
      du(3) = coeffs(5)/denom
      du(4) = -4.0*coeffs(6)*pkoffs(1)/denom
      du(5) = (coeffs(3) + 2.0*coeffs(5)*pkoffs(1))/denom
      du(6) = (-2.0*coeffs(2) - 4.0*coeffs(4)*pkoffs(1))/denom
      dv(1) = 0.0
      dv(2) = du(3)
      dv(3) = -2.0*coeffs(4)/denom
      dv(4) = (-2.0*coeffs(3) - 4.0*coeffs(6)*pkoffs(2))/denom
      dv(5) = (coeffs(2) + 2.0*coeffs(5)*pkoffs(2))/denom
      dv(6) = -4.0*coeffs(4)*pkoffs(2)/denom
C
C  Compute estimated errors in offsets
C-------------------------------------
      usum = 0.0
      vsum = 0.0
      xsum = 0.0
      do 300 i=1,6
	  do 400 j=1,6
	      usum = usum + du(i)*du(j)*bnvrs(i,j)
	      vsum = vsum + dv(i)*dv(j)*bnvrs(i,j)
	      xsum = xsum + du(i)*dv(j)*bnvrs(i,j)
 400      continue
 300  continue
      tlerrs(1) = sqrt(abs(c*usum))
      tlerrs(2) = sqrt(abs(c*vsum))
      tlerrs(3) = c*xsum
C
      return
      end
