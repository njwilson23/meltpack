c****************************************************************************
C $Id: gcorr.f,v 1.1 2000-01-26 14:12:55-07 braup Exp $
C
C $Log: gcorr.f,v $
C Revision 1.1  2000-01-26 14:12:55-07  braup
C Initial revision
C
C
C*    Subroutine gcorr
C
C     Purpose:
C    Correlate a reference subimage with a search subimage and
C    evaluate the results.
C
C     Arguments:
C    Name    Type       I/O    Description
C    ----    ----       ---    -----------
C    IMAGES    I*2 Array   I    Search subimage.
C    IMAGER    I*2 Array   I   Reference subimage.
C    NPLS    real Array  I   Actual size of search subimage:  number
C                of pixels per line, and number of lines.
C    NPLR    real Array  I    Actual size of reference subimage:
C                number of pixels per line, and number of
C                lines.
C    CSMIN    Real        I    Minimum acceptable correlation strength.
C    MFIT    Integer        I    Method of fitting surface
C                 1 - Elliptical paraboloid
C                 2 - Elliptical Gaussian
C                 3 - Reciprocal Paraboloid
C                 4 - Round to nearest integers
C    DDMAX    Real        I    Maximum allowed diagonal displacement
C                from nominal tiepoint location to
C                location found by correlation.
C    IOFFRQ    real Array  I    Requested maximum horizontal and
C                vertical search offsets.
C    NOMOFF    real Array  I    Nominal horizontal and vertical offsets
C                of upper left corner of reference
C                subimage relative to search subimage.
C    IACREJ    Integer        O    Accept/Reject code (see geompak.h)
C    STRENG    Real        O    Strength of correlation
C    BFOFFS    Real Array  O    Best-fit horizontal (pixel) and vertical
C                (line) offsets of correlation peak
C                relative to nominal input location.
C    TLERRS    Real Array  O    Estimated horizontal error, vertical
C                error, and h-v cross term in best-fit
C                offsets.
C    DDACT    Real        O    Actual diagonal displacement from
C                nominal tiepoint location to location
C                found by correlation.
C
C     Local Variables and Arrays:
C    Name    Type        Description
C    ----    ----        -----------
C    CCNORM    Real Array  Normalized cross-correlation coefficient for
C                each search alignment.
C    CPVAL    Real Array  5 by 5 array of cross-correlation values, in
C                units of standard deviations above
C                background, centered on correlation peak.
C    IPKCOL    Int. Array  Column number for each of top 32 correlation
C                values.
C    IPKROW    Int. Array  Row number for each of top 32 correlation
C                values.
C    PKOFFS    Real Array  Fractional horizontal and vertical offsets
C                of best-fit peak location from nearest
C                integer location.
C    PKVAL    Real Array  Largest 32 values of normalized cross-
C                correlation coefficients.
C    SUMS    Real Array  Sum and sum of squares of all cross-
C                correlation values.
C    UNORMC    Real Array  Unnormalized cross-product sums for each
C                alignment.
C
C     Program History:
C    Version      Date          Author    Code/Contr.    Change Request
C    -------      ----          ------    -----------    --------------
C       5.0      Jan 1990    D. Steinwand  EDC     LAS 5.0 Conversion
C       2         JULY 1984   Y. JUN        SAR         SCN#20-305
C       1      Aug. 1983   R. White        CSC        N/A
C
C*    Language:  Fortran 77.
c
C**   ALGORITHM DESCRIPTION
C
C        Compute raw cross-product sums
C        Compute and tabulate normalized cross-correlation values
C        Evaluate strength of correlation peak
C        IF peak not rejected as too weak or too close to edge THEN
C        IF surface fitting specified THEN
C            Do surface fit
C            Add offsets returned by surface fit to the offsets from
C             edges of full array to center of 5 by 5 subarray
C        ELSE
C            Enter integer offsets at array peak as best-fit
C             values
C            Enter (0.5, 0.5) as estimated errors
C        END-IF
C        Compute diagonal displacement from nominal location
C        IF maximum offset parameter was defaulted THEN
C            Use 2.0 less than lesser of two requested search
C             offsets as maximum allowed
C        END-IF
C        IF displacement exceeds maximum allowed THEN
C            Set accept/reject code to 5
C        END-IF
C        END-IF
C**         Return
c
c****************************************************************************
      SUBROUTINE gcorr(images, imager, npls, nplr, csmin, mfit,
     .      ddmx, ioffrq, nomoff,
     .      iacrej, streng, bfoffs, tlerrs, ddact)
C
C
      real      imager(1), images(1)
      real      ioffrq(2), nomoff(2), nplr(2), npls(2)
      integer   iacrej, mfitunormc
      real      bfoffs(2), tlerrs(3)
      real      csmin, ddact, ddmx, streng
C
C  Specification of Local Variables and Arrays
C---------------------------------------------
      integer   ipkcol(32), ipkrow(32), jnpls(2), jnplr(2)
      real      ccnorm(16641), cpval(25), pkoffs(2),
     .          pkval(32), sums(2), unormc(16641)
      integer   ncol, nrow
C
C  Compute raw cross-product sums
C--------------------------------
      call cross(images, imager, nomoff, npls, nplr, unormc)
C
C  Compute normalized cross-correlation values and compile statistics
C--------------------------------------------------------------------
      call gnorm(imager, images, nplr, npls, unormc,
     .           ccnorm, pkval, ipkcol, ipkrow, sums)
C
C  Evaluate strength of correlation peak
C---------------------------------------
      ncol=nint(npls(1))-nint(nplr(1))+1
      nrow=nint(npls(2))-nint(nplr(2))+1
      call eval(ncol, nrow, ccnorm, pkval, ipkcol, ipkrow, sums, csmin,
     .          streng, iacrej, cpval)
C
C  Determine offsets of peak relative to nominal location
C--------------------------------------------------------
      if (iacrej.eq.1) then
          if (mfit.ne.4) then
              call fitreg(cpval, mfit, pkoffs, tlerrs)
              bfoffs(1) = (ipkcol(1) - 1) - nomoff(1) + pkoffs(1)
              bfoffs(2) = (ipkrow(1) - 1) - nomoff(2) + pkoffs(2)
          else
              bfoffs(1) = (ipkcol(1) - 1) - nomoff(1)
              bfoffs(2) = (ipkrow(1) - 1) - nomoff(2)
              tlerrs(1) = 0.5
              tlerrs(2) = 0.5
          end if
C
C  Determine diagonal displacement from nominal and check against maximum
C  acceptable value
C------------------
          ddact = sqrt(bfoffs(1)**2 + bfoffs(2)**2)
          if (ddmx.gt.0.0) then
              if (ddact.gt.ddmx) then
                      iacrej = 5
              end if
          else
              if ( (bfoffs(1)**2 .gt. ioffrq(1)**2) .or.
     1           (bfoffs(2)**2 .gt. ioffrq(2)**2) ) then
                  iacrej = 5
              end if
          end if
C
      end if
C
      return
      end
