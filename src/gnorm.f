c****************************************************************************
C $Id: gnorm.f,v 1.1 2000-01-26 14:12:55-07 braup Exp $
C
C $Log: gnorm.f,v $
C Revision 1.1  2000-01-26 14:12:55-07  braup
C Initial revision
C
C
C*    Subroutine gnorm
C
C     Purpose:
C	Convert raw cross-product sums to normalized cross-correlation
C	coefficients, while tabulating statistics needed for
C	subsequent evaluation.
C
C     Arguments:
C	Name	Type	   I/O	Description
C	----	----	   ---	-----------
C	IMAGER	I*2 Array   I   Reference subimage.
C	IMAGES	I*2 Array   I	Search subimage.
C	NPLR    real Array  I	Actual size of reference subimage:
C				number of pixels per line, and number of
C				lines.
C	NPLS	real Array  I   Actual size of search subimage:  number
C				of pixels per line, and number of lines.
C	UNORMC	Real Array  I	Array of unnormalized (raw) cross-
C				product sums for each alignment of
C				reference subimage relative to search
C				subimage. 
C	CCNORM	Real Array  O	Array of normalized cross-correlation
C				coefficients for each alignment of
C				reference subimage relative to search
C				subimage.
C	PKVAL	Real Array  O	Table of top 32 normalized values.
C	IPKCOL	Int. Array  O	Table of column numbers for top 32
C				values. 
C	IPKROW	Int. Array  O	Table of row numbers for top 32 values.
C	SUMS	Real Array  O   Sum of all normalized values, and sum of
C				squares.
C
C     Local Variables and Arrays:
C	Name	Type	    Description
C	----	----	    -----------
C	COLSQR	Real Array  Sum of squares of search-subimage pixel
C			    values for one column of subarea overlain by
C			    reference subimage in a particular
C			    alignment.
C	COLSUM	Real Array  Sum of search-subimage pixel values for one
C			    column of currently overlain subarea.
C	I	Integer	    Loop index:  pointer to current row of
C			    alignments.
C	IPFREE	Integer	    Index to elements in value and location
C			    arrays which are freed by deleting 32-nd
C			    smallest.
C	IPT	Integer	    Pointer into image array.
C	IPTR	Int. Array  Pointers to value and location arrays.
C	IXCOL	Int. Array  Column locations of 32 largest values, with
C			    sequencing specified by IPTR array.
C	IXROW	Int. Array  Row locations of 32 largest values.
C	J	Integer	    Loop index:  pointer to element in current
C			    row of alignments.
C	JSTART	Integer	    Pointer to first element in current row.
C	JSTOP	Integer	    Pointer to last element in current row.
C	K	Integer	    Loop index:  pointer into table of peak
C			    values.
C	KOL	Integer	    Loop index:  column of search subimage.
C	KOLADD	Integer	    Number of search-subimage column being
C			    added, on going to new horizontal alignment,
C			    to subarea overlain by reference subimage.
C	KOLSUB	Integer	    Number of search-subimage column being
C			    subtracted from overlain subarea.
C	LINE	Integer	    Loop index:  current search-subimage line.
C	LNADD	Integer	    Number of search-subimage line being added,
C			    as result of new vertical alignment, to 
C			    subarea overlain by reference subimage.
C	LNSUB	Integer	    Number of search-subimage line being
C			    subtacted from overlain subarea.
C	NCOL	Integer	    Number of columns in cross-product-sum
C			    array.
C	NROW	Integer	    Number of rows in cross-product-sum	array.
C	NRTOT	Integer	    Total number of pixels in reference
C			    subimage.
C	REFSUM	Real	    Sum of all reference-subimage pixel values.
C	REFSQR	Real	    Sum of squares of reference-subimage pixels.
C	RHO	Real	    Temporary storage for cross-correlation
C			    coefficient.
C	RMEAN	Real	    Mean of reference-subimage pixel values.
C	RSIGMA	Real	    Standard deviation about the mean of
C			    reference-subimage pixels values.
C	RTOTAL	Real	    Number of pixels in reference subimage.
C	SIGMA1	Real	    Value obtained after first iteration of
C			    Newton's method to approximate square-root
C			    to get search-subarea standard deviation.
C	SIGMAS	Real	    Standard deviation, multiplied by RTOTAL,
C			    of search-subarea pixel values about mean;
C			    also serves as initial estimate of value
C			    for next alignment in row.
C	SRCHSM	Real	    Sum of pixel values for search subarea
C			    currently overlain by reference subimage.
C	SRCHSQ	Real	    Sum of squares of search-subarea pixels.
C	TEMP	Real	    Temporary variable to simplify computation
C			    of search-subarea standard deviation.
C	TEMPMN	Real	    Minimum allowed value of TEMP to avoid
C			    division by very small values if search
C			    subarea is nearly uniform.
C	XVAL	Real Array  Temporary storage for largest values, with
C			    sequencing specified by IPTR array.
C
C     Program History:
C	Version	  Date	      Author	Code/Contr.	Change Request
C	-------	  ---------   ------	-----------	--------------
C	5.0       Jan 1990    D. Steinwand   EDC 	LAS 5.0 Conversion
C       2.0       JULY 1984   Y. JUN         SAR        SCN#20-305
C	001	  Aug. 1983   R. White	     CSC	N/A
C
C*    Language:  Fortran 77.
C
C**   ALGORITHM DESCRIPTION
C
C	    Tabulate sum and sum of squares of pixel values in reference
C	     subimage
C	    Compute mean, standard deviation, and reciprocal of standard
C	     deviation for reference subimage
C	    Clear sum of normalized cross-correlation values and sum of
C	     squares
C	    Initialize table of 32 largest values
C	    DO-FOR each row of unnormalized cross-product values
C		DO-FOR each column of search subimage
C		    IF topmost row THEN
C			Sum pixel values and values squared for as many
C			 lines as there are lines in reference subimage
C		    ELSE
C			Adjust column sum and sum of squares for new
C			 pixel added at bottom of column and old pixel
C			 lost at top
C		    END-IF
C		END-DO
C		DO-FOR each unnormalized value in row
C		    IF leftmost alignment THEN
C			Add up columns sums of pixel values and sums of 
C			 squares for all search columns overlain by
C			 reference subimage
C			Compute mean and exact standard deviation
C		    ELSE
C			Adjust pixel sums for new column added and old
C			 column lost
C			Compute mean grey-level value
C			Approximate the standard deviation using two
C			 iterations of Newton's method
C		    END-IF
C		    Save standard deviation as initial estimate for next
C		     alignment in row
C		    Compute normalized cross-correlation coefficient
C		    Add normalized value and its square to sums
C		    IF value exceeds 32-nd largest so far THEN
C			Insert value into ordered list of top 32, using
C			 pointer array to reduce data shuffling
C			Insert horizontal and vertical offsets for this
C			 value into table of offsets
C		    END-IF
C		END-DO
C	    END-DO
C	    Copy top 32 values and coordinates into output arrays
C**         Return
C
c****************************************************************************
      SUBROUTINE gnorm(imager, images, nplr, npls, unormc, ccnorm,
     .			   pkval, ipkcol, ipkrow, sums)
      real	imager(1),	images(1)
      integer	ipkcol(1),	ipkrow(1)
      real	ccnorm(1),	pkval(1),	sums(2),
     .		unormc(1),      nplr(2),        npls(2)
C
C  Specification of Local Variables and Arrays
C---------------------------------------------
      integer	iptr(32),	ixcol(32),
     .          nsnew(2),       nrnew(2),	ixrow(32)
      integer	i,	ipfree,		ipt,		j,
     .	jstart,		jstop,		k,		kol,
     .	koladd,		kolsub,		line,		lnadd,
     .	lnsub,		ncol,		nrow,		nrtot,
     .  igl, iglnew, iglold
      real colsqr(256),	colsum(256),	xval(32)
      real refsum,	refsqr,		rho,		rmean,
     .     rsigma,	rtotal,		sigma1,		sigmas,
     .     srchsm,	srchsq,		temp,		tempmn
C
C  Compute pixel-value statistics for reference subimage
C-------------------------------------------------------
      nsnew(1) = nint(npls(1))
      nsnew(2) = nint(npls(2))
      nrnew(1) = nint(nplr(1))
      nrnew(2) = nint(nplr(2))
      refsum = 0
      refsqr = 0
      nrtot = nrnew(1) * nrnew(2)
      do 100 ipt=1,nrtot
	  igl = imager(ipt)
	  refsum = refsum + igl
	  refsqr = refsqr + igl*igl
 100  continue
      rtotal = nrtot
      tempmn = .01 * rtotal**2
      rmean = refsum/rtotal
      rsigma = sqrt(amax1((refsqr/rtotal - rmean*rmean), 0.01))
C
C  Clear sums and sums of squares of normalized correlation values
C-----------------------------------------------------------------
      sums(1) = 0.0
      sums(2) = 0.0
      do 200 k=1,32
	  xval(k) = -1.0
	  iptr(k) = k
 200  continue
      ncol = nsnew(1) - nrnew(1) + 1
      nrow = nsnew(2) - nrnew(2) + 1
C
C  Compute normalized cross-corr values for one row of alignemnts at a time
C--------------------------------------------------------------------------
      jstart = 1
      jstop = ncol
      do 300 i=1,nrow
C
C  Get column sums and sums of squares for portion of search
C  subimage overlain by reference in current row of alignments
C-------------------------------------------------------------
	  if (i .eq. 1) then
	      do 310 kol=1,nsnew(1)
		  colsum(kol) = 0.0
		  colsqr(kol) = 0.0
 310          continue
	      ipt = 1
	      do 320 line=1,nrnew(2)
		  do 330 kol=1,nsnew(1)
		      igl = images(ipt)
		      colsum(kol) = colsum(kol) + igl
		      colsqr(kol) = colsqr(kol) + igl*igl
		      ipt = ipt + 1
 330    	  continue
 320          continue
	  else
	      lnsub = (i - 2) * nsnew(1)
	      lnadd = lnsub + nsnew(1)*nrnew(2)
	      do 340 kol=1,nsnew(1)
		  iglnew = images(lnadd+kol)
		  iglold = images(lnsub+kol)
		  colsum(kol) = colsum(kol) + iglnew - iglold
		  colsqr(kol) = colsqr(kol) + iglnew*iglnew
     .				- iglold*iglold
 340          continue
	  end if
C
C  Complete comutation of search-subarea pixel statistics
C--------------------------------------------------------
	  do 350 j=jstart,jstop
	      if (j .eq. jstart) then
		  srchsm = 0.0
		  srchsq = 0.0
		  do 360 kol=1,nrnew(1)
		      srchsm = srchsm + colsum(kol)
		      srchsq = srchsq + colsqr(kol)
 360              continue
		  temp = amax1(tempmn,(rtotal*srchsq - srchsm*srchsm))
		  sigmas = sqrt(temp)
	      else
		  kolsub = j - jstart
		  koladd = kolsub + nrnew(1)
		  srchsm = srchsm + colsum(koladd) - colsum(kolsub)
		  srchsq = srchsq + colsqr(koladd) - colsqr(kolsub)
		  temp = amax1(tempmn,(rtotal*srchsq - srchsm*srchsm))
		  sigma1 = 0.5 * (sigmas + temp/sigmas)
		  sigmas = 0.5 * (sigma1 + temp/sigma1)
	      end if
C
C  Compute normalized cross-correlation value
C--------------------------------------------
	      rho = (unormc(j) - rmean*srchsm) / (rsigma*sigmas)
	      ccnorm(j) = rho
	      sums(1) = sums(1) + rho
	      sums(2) = sums(2) + rho*rho
C
C  Check whether value among top 32
C----------------------------------
	      if (rho.gt.xval(iptr(32))) then
		  k = 32
		  ipfree = iptr(32)
 700                  if((k.le.1).or.(rho.le.xval(iptr(k-1))))goto 750
    		      iptr(k) = iptr(k-1)
		      k = k - 1
		      goto 700
 750		  continue
		  iptr(k) = ipfree
		  xval(ipfree) = rho
		  ixcol(ipfree) = j - jstart + 1
		  ixrow(ipfree) = i
	      end if
 350      continue
	  jstart = jstart + ncol
	  jstop =jstop + ncol
 300  continue
C
C  Copy peak values and coordinates in correct sequence
C------------------------------------------------------
      do 400 k=1,32
	  pkval(k) = xval(iptr(k))
	  ipkcol(k) = ixcol(iptr(k))
	  ipkrow(k) = ixrow(iptr(k))
 400  continue
C
      return
      end
