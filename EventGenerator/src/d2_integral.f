C
C  Part of the Daya Bay Cosmic code.
C
C  $Id: d2_integral.f,v 1.3 2014/03/22 21:40:44 ehrlich Exp $
C  $Author: ehrlich $
C  $Date: 2014/03/22 21:40:44 $
C
C  Original Mu2e author Yury Kolomensky
C

*********************************************************************
* do integral of the Funcation ||
*  method Divide to sub mini bins & add content ||
*********************************************************************
      Subroutine d2_integral(a1,a2,b1,b2,sum,m,n,vert)
      implicit none
      real*8 gaisser
      EXTERNAL gaisser
      integer i,j, m,n
      real*8 a1,a2,b1,b2
      real*8 sum,sumbox,binx,biny
      real*8 x,y
      Real*8 fff(m,n)
      logical vert

**********************************
      binx=(a2-a1)/m
      biny=(b2-b1)/n

*--> book 2d matrix, Get the sum
**********************************
      sumbox=0.d0
      do i=1,m
       do j=1,n
         x=a1+binx*(i-0.5)
         y=b1+biny*(j-0.5)
         IF(vert) THEN
           fff(i,j)=gaisser(x,y)*SIN(ACOS(y))
         ELSE
           fff(i,j)=gaisser(x,y)*y
         ENDIF
         sumbox=fff(i,j)+sumbox
       enddo
      enddo

**********************************

*--> do integral
**********************************
      SUM = sumbox*binx*biny
*
*      print*,'binx,y,z; sum: ',binx,biny,binz,sum
**********************************

      end


