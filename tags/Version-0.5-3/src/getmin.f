c  Copyright (C) 1995-2010 Berwin A. Turlach <berwin@maths.uwa.edu.au>
c
c  This program is free software; you can redistribute it and/or modify
c  it under the terms of the GNU General Public License as published by
c  the Free Software Foundation; either version 2 of the License, or
c  (at your option) any later version.
c
c  This program is distributed in the hope that it will be useful,
c  but WITHOUT ANY WARRANTY; without even the implied warranty of
c  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c  GNU General Public License for more details.
c
c  You should have received a copy of the GNU General Public License
c  along with this program; if not, write to the 
c     Free Software Foundation
c     59 Temple Place - Suite 330, 
c     Boston, MA 02111-1307, USA.
c 
c  Given function values on a grid, this routine tries to estimate the
c  minima of the function.
c  If the minimal value is found in the interior, the minimum is
c  determined by a quadratic fit through the point where the minimal
c  value occurred plus the points left and right to it.
c
c  Input:
c     x    nx1 vector, dp, the x-values at which the function is evaluated.
c     y    nx1 vector, dp, the value of the functions at the x-coordinates.
c     n    integer, length of x and y.
c     flag integer, if 0 the global minima is searched.  If <0 we search
c                   for the left-most local minima, otherwise for the
c                   right-most local minima.
c     nmin integer, if not 0 the number of local minima is counted.
c
c  output:
c     xmin  dp, the x-value at which the minima is approx. located.
c     ymin  dp, the approx. minima
c     flag  integer, if -1 the minima is at the left end.
c                    if  1 the minima is at the right end.
c                    if  0 the minima is in the middle of the data.
c                    if  5 the minima is in the middle, but the
c                          quadratic fit went wrong.
c     nmin  integer, either 0 or the number of local minima found.
c
      subroutine getmin(x,y,n,xmin,ymin,flag,nmin)
      implicit none
      integer n, flag, i, imin, nmin
      double precision x(*), y(*), xmin, ymin, a1, a2, a3, c0, c1, c2, 
     +     cmin
      if( flag .EQ. 0 )then
         cmin = y(1)
         imin = 1
         do 10 i=2,n
            if( y(i) .LT. cmin )then
               cmin = y(i)
               imin = i
            endif
 10      continue
      elseif( flag .LT. 0)then
         cmin = y(1)
         imin = 1
         do 20 i=2,n
            if( y(i) .LT. cmin )then
               cmin = y(i)
               imin = i
            else
               goto 21
            endif
 20      continue
 21      continue
      else
         cmin = y(n)
         imin = n
         do 30 i=n-1,1,-1
            if( y(i) .LT. cmin )then
               cmin = y(i)
               imin = i
            else
               goto 31
            endif
 30      continue
 31      continue
      endif
      if( imin .EQ. 1 )then
         xmin = x(1)
         ymin = y(1)
         flag = -1
      elseif( imin .EQ. n )then
         xmin = x(n)
         ymin = y(n)
         flag = 1
      else
         a1 = y(imin-1)/((x(imin-1)-x(imin  ))*(x(imin-1)-x(imin+1)))
         a2 = y(imin  )/((x(imin-1)-x(imin  ))*(x(imin  )-x(imin+1)))
         a3 = y(imin+1)/((x(imin-1)-x(imin+1))*(x(imin  )-x(imin+1)))
         c2 = a1-a2+a3
         c1 = -(x(imin  )+x(imin+1))*a1
     +        +(x(imin-1)+x(imin+1))*a2
     +        -(x(imin-1)+x(imin  ))*a3
         c0 = a1*x(imin)*x(imin+1) - a2*x(imin-1)*x(imin+1) +
     +        a3*x(imin-1)*x(imin)
         xmin = -c1/(2*c2)
         if( x(imin-1) .LE. xmin .AND. xmin .LE. x(imin+1) )then
            flag = 0
            ymin = (c2*xmin+c1)*xmin+c0
         else
            xmin = x(imin)
            ymin = y(imin)
            flag = 5
         endif
      endif
      if( nmin .EQ. 0 ) return
      nmin = 0
      if( y(1) .LE. y(2) ) nmin = nmin+1
      do 40 i=2,n-1
         if( y(i) .LE. y(i-1) .AND. y(i) .LE. y(i+1) ) nmin = nmin+1
 40   continue
      if( y(n) .LE. y(n-1) ) nmin = nmin+1
      return
      end
      
