c****************************************************************************
C $Id: eval.f,v 1.1 2000-01-26 14:12:55-07 braup Exp $
C
C $Log: eval.f,v $
C Revision 1.1  2000-01-26 14:12:55-07  braup
C Initial revision
C
C
C*    Subroutine eval
C
C     Purpose:
C	Evaluate various measures of correlation validity and extract
C	a subarea of the cross correlation array centered on the peak.
C
C     Arguments:
C	Name	Type	   I/O	Description
C	----	----	   ---	-----------
C	NCOL	Integer	    I	Number of columns in coincidence-count array.
C	NROW	Integer	    I	Number of rows in coincindence-count array.
C	CCNORM	Real Array  I	Array of normalized cross-correlation
C				coefficients for each alignment of
C				reference edge image relative to search
C				edge image.
C	PKVAL	Real Array  I	Table of top 32 normalized values.
C	IPKCOL	Int. Array  I	Table of column numbers for top 32 values. 
C	IPKROW	Int. Array  I	Table of row numbers for top 32 values.
C	SUMS	Real Array  I   Sum of all normalized values and sum of squares
C	CSMIN	Real	    I	Minimum acceptable correlation strength.
C	STRENG	Real	    O	Strength of correlation 
C	IACREJ	Integer	    O	Accept/Reject code (see geompak.h)
C	CPVAL	Real Array  O	5 by 5 array of cross-correlation
C				values, in units of standard deviations
C				above background, centered on
C				correlation peak.
C
C     Local Variables and Arrays:
C	Name	Type	    Description
C	----	----	    -----------
C	BMEAN	Real	    Mean of background cross-correlation values.
C	BSIGMA	Real	    Standard deviation of background about mean.
C	I	Integer	    Loop index.
C	ICOL	Integer	    Initial column in 9 by 9 neighborhood
C			    (unless size reduced by proximity to edges
C			    of correlatin array) of correlation peak.
C	IDIST	Integer	    Maximum horizontal or vertical distance from
C			    array peak.
C	IPTR	Integer	    Pointer into 5 by 5 array of output values.
C	IPT5	Int. Array  Pointers to highest two values more than 2
C			    columns or rows from peak value.
C	IPT7	Integer	    Pointer to largest value more than 3 columns
C			    or rows from peak.
C	IROW	Integer	    Initial row in 9 by 9 neighborhood.
C	KRBASE	Integer	    Offset to beginning of current row in array.
C	LCOL	Integer	    Last column in 9 by 9 neighborhood centered
C			    on correlation peak.
C	LROW	Integer	    Last row in 9 by 9 neighborhood.
C	NPTS	Integer	    Number of array values outside 9 by 9
C			    neighborhood.
C	N5	Integer	    Number of large values outside 5 by 5
C			    neighborhood of peak.
C	N7	Integer	    Number of large values outside 7 by 7
C			    neighborhood of peak.
C
C     Program History:
C	Version	  Date	      Author	Code/Contr.	Change Request
C	-------	  ----	      ------	-----------	--------------
C  	5.0       Jan 1990    D. Steinwand   EDC        LAS 5.0 Conversion
C	001	  June 1983   R. White	     CSC		N/A
C
C*    Language:  Fortran 77.
C
C**   ALGORITHM DESCRIPTION
C
C	    Initialize accept/reject code to 1
C	    Initialize correlation strength to 0.0
C	    Initialize pointers to two largest values in 5 by 5 array
C	     centered on peak and largest value in 7 by 7
C	    IF largest value occurred less than two rows or columns from
C	     edge of correlation array THEN
C		Set accept/reject code to 2
C		Return to caller
C	    END-IF
C	    Find largest values outside 5 by 5 and 7 by 7 neighborhoods
C	     of peak
C	    IF any value among top 3 or more than 1 value among top
C	     5 is more than 2 rows or columns from peak THEN
C		Set accept/reject code to 3
C		Return to caller
C	    END-IF
C	    Find edges of 9 by 9 neighborhood of peak, reduced if
C	     necessary by proximity to edges of full array
C	    DO-FOR all alignments within "9 by 9" neighborhood
C		Subtract normalized cross-correlation value from sum
C		Subtract square of value from sum of squares
C	    END-DO
C	    Compute mean and standard deviation for all values
C	     outside the 9 by 9 subarray
C	    Compute correlation strength (FD Sect. 4.3)
C	    IF correlation strength less than minimum THEN
C		Set accept/reject code to 4
C		Return to caller
C	    END-IF
C	    Define 5 by 5 subarray centered on peak
C	    DO-FOR each value in subarray
C		Convert value to standard deviations above
C		 background
C	    END-DO
C**     Return
C
c****************************************************************************
      SUBROUTINE eval(ncol, nrow, ccnorm, pkval, ipkcol, ipkrow, sums,
     .		    csmin, streng, iacrej, cpval)
C
C  Specification of Subroutine Arguments
C---------------------------------------
      integer	iacrej,		ncol,		nrow
      integer	ipkcol(32),	ipkrow(32)
      real	ccnorm(1),	cpval(25),	pkval(32),	sums(2)
      real	csmin,		streng
C
C  Specification of Local Variables and Arrays
C---------------------------------------------
      integer	ipt5(2)
      integer	i,	icol,	idist,	iptr,	ipt7,	irow,	krbase,
     .		lcol,	lrow,	npts,	n5,	n7, j
      real	bmean,	bsigma
C
C  Initialize accept/reject code, strength, and outlier pointers
C---------------------------------------------------------------
      iacrej = 1
      streng = 0.0
      ipt5(1) = 32
      ipt5(2) = 32
      ipt7 = 1
C
C  Check for peak value within two rows or columns of edge
C---------------------------------------------------------
      if ((ipkcol(1).le.2) .or. (ipkcol(1).ge.(ncol-1)) .or.
     .	  (ipkrow(1).le.2) .or. (ipkrow(1).ge.(nrow-1))) then
	  iacrej = 2
	  return
      end if
C
C  Find largest values outside 5 by 5 and 7 by 7 neighborhoods of peak
C---------------------------------------------------------------------
      n5 = 0
      n7 = 0
      i = 2
  50  if ((n5.ge.2).or.(i.gt.32)) goto 100
	  idist = max0(iabs(ipkcol(1)-ipkcol(i)),
     .		      iabs(ipkrow(1)-ipkrow(i)))
	  if (idist.gt.2) then
	      n5 = n5 + 1
	      ipt5(n5) = i
	      if ((n7.eq.0).and.(idist.gt.3)) then
		  ipt7 = i
		  n7 = 1
	      end if
	  end if
	  i = i + 1
	  goto 50
 100  continue
      if ((ipt5(1).le.3) .or. (ipt5(2).le.5)) then
	  iacrej = 3
	  return
      end if
C
C  Find edges of 9 by 9 array centered on peak
C---------------------------------------------
      icol = max0(1, ipkcol(1) - 4)
      irow = max0(1, ipkcol(1) - 4)
      lcol = min0(ncol, ipkcol(1) + 4)
      lrow = min0(nrow, ipkrow(1) + 4)
C
C  Eliminate points within 9 by 9 array from background statistics
C-----------------------------------------------------------------
      krbase = ncol*(irow-1)
      do 200 i=irow,lrow
	  do 225 j=icol,lcol
	      sums(1) = sums(1) - ccnorm(krbase+j)
	      sums(2) = sums(2) - ccnorm(krbase+j)*ccnorm(krbase+j)
 225      continue
	  krbase = krbase + ncol
 200  continue
C
C  Compute background mean and standard deviation
C------------------------------------------------
      npts = ncol*nrow - (lcol - icol + 1) *(lrow - irow + 1)
      bmean = sums(1)/npts
      bsigma = sqrt(sums(2)/npts - bmean*bmean)
C
C  Compute correlation strength and check against minimum.  Note that the
C   LAS 4.x version of this code contains an error...it often tries to
C   access pkval array element zero when there are no signifcant 
C   subsidiary peaks, resulting in unpredictable results
C--------------------------------------------------------
      if (n7 .eq. 0) then
         streng = 2 * ((pkval(1) - bmean)/bsigma) - 0.2
      else
         streng = (pkval(1) - bmean)/bsigma
     .		+ (pkval(1) - pkval(ipt7))/bsigma
     .		+ 0.2*(n7 - 1.0)
      end if

      if (streng.lt.csmin) then
	  iacrej = 4
	  return
      end if
C
C  Convert 5 by 5 neighborhood of peak to standard deviations above mean
C-----------------------------------------------------------------------
      krbase = ncol*(ipkrow(1) - 3)
      iptr = 1
      do 300 I=1,5
	  do 325 j=ipkcol(1)-2,ipkcol(1)+2
	      cpval(iptr) = (ccnorm(krbase+j) - bmean)/bsigma
	      iptr = iptr + 1
 325      continue
	  krbase = krbase + ncol
 300  continue
C
      return
      end
