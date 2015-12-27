c****************************************************************************
C $Id: fft2d.f,v 1.1 2000-01-26 14:12:55-07 braup Exp $
C
C $Log: fft2d.f,v $
C Revision 1.1  2000-01-26 14:12:55-07  braup
C Initial revision
C
c  
c  NAME:			FFT2D
c
c  PURPOSE:  Performs a two dimensional Fast Fourier Transformation
c  
c  PROGRAM HISTORY:
c  VERSION	 DATE	AUTHOR	   CODE/CONT   REASON
c  -------	 ----	------	   ---------   -----------------------------
c    5.0	1/90	D. Steinwand  CSB      Adaption from text (see below)
c  
c  COMPUTER HARDWARE AND/OR SOFTWARE LIMITATIONS:   must be run under TAE
c  
c  PROJECT:  LAS
c
c  ALGORITHM REFERENCES:
c	This subroutine has been derrived from the FOURN routine in
c
c	Press, William H., Flannery, Brian P., Teukolsky, Saul A., and
c	 Vetterling, William T.  1986, "Numerical Recipes--The Art of
c  	 Scientific Computing" (New York:  Cambridge University Press).
c
c	If isign = 1, the DATA array will be replaced by its two dimensional
c	discrete fourier transform.
c	If isign = -1, the DATA array will be replaced by its two dimensional
c	inverse transform times the product of the lengths of all dimensions.
c	DATA should be input as a multidimensional complex array.
c****************************************************************************
      subroutine fft2d(data,nel,isign)
      double precision wr,wi,wpr,wpi,wtemp,theta
      integer nel(2)
      dimension data(*)

      ntot=nel(1)*nel(2)
c
c  Loop over both dimensions
c---------------------------
      nprev=1
      do 18 idim=1,2
        n=nel(idim)
        nrem=ntot/(n*nprev)
        ip1=2*nprev
        ip2=ip1*n
        ip3=ip2*nrem
        i2rev=1
c
c  Perform bit reversal
c----------------------
        do 14 i2=1,ip2,ip1
          if(i2.lt.i2rev)then
            do 13 i1=i2,i2+ip1-2,2
              do 12 i3=i1,ip3,ip2
                i3rev=i2rev+i3-i2
                tempr=data(i3)
                tempi=data(i3+1)
                data(i3)=data(i3rev)
                data(i3+1)=data(i3rev+1)
                data(i3rev)=tempr
                data(i3rev+1)=tempi
12            continue
13          continue
          endif
          ibit=ip2/2
1         if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
            i2rev=i2rev-ibit
            ibit=ibit/2
            go to 1
          endif
          i2rev=i2rev+ibit
14      continue
c
c  Danielson-Lanczos section
c---------------------------
        ifp1=ip1
2       if(ifp1.lt.ip2)then
          ifp2=2*ifp1
c
c  Initialize for the trigonometric recurrence
c---------------------------------------------
          theta=isign*6.28318530717959d0/(ifp2/ip1)
          wpr=-2.d0*dsin(0.5d0*theta)**2
          wpi=dsin(theta)
          wr=1.D0
          wi=0.D0
          do 17 i3=1,ifp1,ip1
            do 16 i1=i3,i3+ip1-2,2
              do 15 i2=i1,ip3,ifp2
c
c  The Danielson-Lanczos formula:
c-------------------------------
                k1=i2
                k2=k1+ifp1
                tempr=sngl(wr)*data(k2)-sngl(wi)*data(k2+1)
                tempi=sngl(wr)*data(k2+1)+sngl(wi)*data(k2)
                data(k2)=data(k1)-tempr
                data(k2+1)=data(k1+1)-tempi
                data(k1)=data(k1)+tempr
                data(k1+1)=data(k1+1)+tempi
15            continue
16          continue
c
c  The trigonometric recurrence
c------------------------------
            wtemp=wr
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
17        continue
          ifp1=ifp2
          go to 2
        endif
        nprev=n*nprev
18    continue
      return
      end
