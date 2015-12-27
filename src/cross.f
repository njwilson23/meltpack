c****************************************************************************
C $Id: cross.f,v 1.1 2000-01-26 14:12:55-07 braup Exp $
C
C $Log: cross.f,v $
C Revision 1.1  2000-01-26 14:12:55-07  braup
C Initial revision
C
C
C*    Subroutine cross
C
C     Purpose:
C	Compute the unnormalized (raw) sum of pixel-by-pixel cross
C	products between the reference and search images for every 
C	combination of horizontal and vertical offsets of the reference
C	relative to the search image.
C
C     Arguments:
C	Name	Type	   I/O	Description
C	----	----	   ---	-----------
C	IMAGES	I*2 Array   I	Search subimage.
C
C	IMAGER	I*2 Array   I   Reference subimage.
C
C	NOMOFF	real Array  I	Nominal horizontal and vertical offsets
C				of upper left corner of reference
C				subimage relative to search subimage.
C	NPLS	real Array  I   Actual size of search subimage:  number
C				of pixels per line, and number of lines.
C	NPLR    real Array  I	Actual size of reference subimage:
C				number of pixels per line, and number of
C				lines.
C	UNORMC	Real Array  O	Array of unnormalized (raw) counts of
C				edge coincidences for each alignment of
C				reference image relative to search
C				image. 
C
C     Local Variables and Arrays:
C	Name	  Type	      Description
C	----	  ----	      -----------
C	IMGPTR	Integer	    Pointer into subimage arrays.
C	LINE	Integer	    Loop index:  current buffer line
C	LNCONT	Integer	    Pointer to start of next segment of line.
C	LNSTRT	Integer	    Pointer to first segment of line.
C	MEMBFR	I*2 Array   Buffer in which reference subimage is
C			    properly aligned and zero padded for AP.
C	MEMBFS	I*2 Array   Buffer for aligning and padding search
C			    subimage.
C	MEMDIM	Int. Array  Power-of-2 dimensions for AP arrays. 
C	MEMPTR	Integer	    Loop index:  pointer into AP transfer
C			    buffers.
C	NCOL	Integer	    Number of columns in cross-product sums
C			    array.
C	NDXOUT	Integer	    Pointer into array for correlation output.
C	NROW	Integer	    Number of rows in cross-product sums array.
C       CSER CMPLX ARRAY    COMPLEX OF SEARCH IMAGE, PRODUCT, AND COORELATION
C       CREF CMPLX ARRAY    COMPLEX OF REFRENCE IMAGE
C
C     File References:  None.
C
C     Subroutines Called:
C	Name	  Purpose of Call
C	----	  ---------------
C       FFT2D     Two dimensional FFT of a complex array
C       
C     Program History:
C	Version	  Date	      Author	Code/Contr.	Change Request
C	-------	  ----	      ------	-----------	--------------
C       5.0       Jan 1990    D. Steinwand  EDC         LAS 5.0 Conversion
C       3.0       JULY 1985   S .KONDAL     SAR         NON AP VERSION
C       2.0       JULY 1984   Y. JUN        SAR         SCN#20-305
C	1	  Aug. 1983   R. White/G. Neal	CSC	N/A
C
C*    Language:  Fortran 77.
c
C**   ALGORITHM DESCRIPTION
C
C       IF not successful after preset time THEN
C       	Report fatal error and terminate execution
C       End-if
C       Determine lowest powers of 2 not less than dimensions of search image
C       Zero fill images to power-of-2 dimesnions.
C       Take complex of subimages
C       Take FFT of complex data
C       Take complex conjugate of one image
C       Make point by point product of FFT'S
C       Take inverse FFT of product
C       Extract the part of correlation array which is valid
C**     Return
c
c****************************************************************************
      SUBROUTINE cross(images, imager, nomoff, npls, nplr, unormc)
C
C  Specification of Subroutine Arguments
C---------------------------------------
      real	imager(1),	images(1)
      real	nomoff(2),	nplr(2),	npls(2)
      real	unormc(1)
C
C  Specification of Local Variables and Arrays
C---------------------------------------------
      real	membfr(65536),	membfs(65536)
      integer	imgptr,
     .		line,		lncont,
     .		lnstrt,	memptr,
     .		ncol,		ndxout,
     .		nrow
      integer   memdim(2),      nsnew(2),       nrnew(2)
      integer i,j,n1,n2,n3
      integer nel(2)
      real denom
      complex cser(65536),    cref(65536)
C
      nsnew(1) = nint(npls(1))
      nsnew(2) = nint(npls(2))
      nrnew(1) = nint(nplr(1))
      nrnew(2) = nint(nplr(2))
C
C  Zero extend binary search image to next higher power of 2
C-----------------------------------------------------------
      do 100 i = 1,2
         memdim(i) = 64
  50     if (nsnew(i) .le. memdim(i)) goto 100
      	 memdim(i) = memdim(i) * 2
	 goto 50
 100  continue
C
C  Fix for incorrect handling of non-square arrays
C-------------------------------------------------
      memdim(1) = max(memdim(1), memdim(2))
      memdim(2) = memdim(1)
      lnstrt = 1
      imgptr = 1
      do 200 line = 1,memdim(2)
	  if (line.le.nsnew(2)) then
	      lncont = lnstrt + nsnew(1)
	      do 225 memptr=lnstrt,lncont-1
		  membfs(memptr) = images(imgptr)
		  imgptr = imgptr + 1
 225          continue
	  else
	      lncont = lnstrt
	  end if
	  lnstrt = lnstrt + memdim(1)
	  if (lncont.lt.lnstrt) then
	      do 250 memptr=lncont,lnstrt-1
		  membfs(memptr) = 0.0
 250          continue
	  end if
 200  continue
C
C  Make comlex 2-d of search data
C--------------------------------
      memptr = 1
      do 300 i = 1, memdim(2)
         do 350 j = 1, memdim(1)
           cser((j-1)*memdim(2)+i) = cmplx(membfs(memptr))
           memptr = memptr + 1
 350     continue
 300  continue
C
C  Now zero extend binary reference image
C----------------------------------------
      lnstrt = 1
      imgptr = 1
      do 400 line = 1,memdim(2)
	  if (line.le.nrnew(2)) then
	      lncont = lnstrt + nrnew(1)
	      do 425 memptr=lnstrt,lncont-1
		  membfr(memptr) = imager(imgptr)
		  imgptr = imgptr + 1
 425          continue
	  else
	      lncont = lnstrt
	  end if
	  lnstrt = lnstrt + memdim(1)
	  if (lncont.lt.lnstrt) then
	      do 450 memptr=lncont,lnstrt-1
		  membfr(memptr) = 0.0
 450          continue
	  end if
 400  continue
C
C  Make comlex 2-d of refrence data
C----------------------------------
      memptr = 1
      do 500 i = 1, memdim(2)
         do 525 j = 1, memdim(1)
           cref((j-1)*memdim(2)+i) = cmplx(membfr(memptr))
           memptr = memptr + 1
 525     continue
 500  continue
C
      n1 = memdim(2)
      n2 = memdim(1)
      n3 = 1
      nel(1) = n1
      nel(2) = n2
C
C  Take fft of search and refrence data
C--------------------------------------
      call fft2d(cser,nel,1)
      call fft2d(cref,nel,1)
C
C  Make point by multiplication of ft of search ft with conjogate of
C  refrence subimage ft
C----------------------
      do 600 i = 1, n1*n2
           cser(i) = cser(i)*conjg(cref(i))
 600  continue
C
C  Take inverse fft of cser
C--------------------------
      call fft2d(cser,nel,-1)
C
C  Extract useful valid correlation
C----------------------------------
      ncol = nsnew(1) - nrnew(1) + 1
      nrow = nsnew(2) - nrnew(2) + 1
      denom = nel(1)*nel(2)
      ndxout = 1
      do 700 i = 1, nrow
        do 725 j = 1, ncol
	  unormc(ndxout) = real(cser((j-1)*memdim(2)+i)) / denom
          ndxout = ndxout + 1
 725    continue
 700  continue
C
      return
      end
