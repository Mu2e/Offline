// $Id: SquareMapper.cc,v 1.4 2013/07/25 23:56:46 echenard Exp $
// $Author: echenard $
// $Date: 2013/07/25 23:56:46 $
//
// Sqaure position map generator: 
//   tesselate a plane with squares starting from the center of the plane
//
// Use basis vector, l and k, defined as  
// l = right
// k = up
//
/*

 -----------------------
 |       |       |      |
 | -1 1  | 0 1   | 1 1  |
 |       |       |      |
 -----------------------
 |       |       |      |  l,k coordinates
 | -1 0  | 0 0   | 1 0  |
 |       |       |      |
 -----------------------
 |       |       |      |
 | -1 -1 | 0 -1  | 1 -1 |
 |       |       |      |
 -----------------------

  steps :  (1,0), (0,-1), (-1,0), (0,1) (clockwise from top left corner)
  segment: top=0, right=1, bottom=2,left=3

*/
// Tesselation algorithm: tessalate in "rings" from the center
//   for each ring, start at -l,+l (top left corner), 
//   then go n time each step to create the ring 
//
// Neighbors add (+-1,0) or (0,+-1)
// next ring of neighbours, add +(-2,2) and go around the ring,...
//


// C++ includes
#include <iostream>
#include <map>
#include <cmath>

// Mu2e includes
#include "CalorimeterGeom/inc/SquareMapper.hh"

//CLHEP includes
#include "CLHEP/Vector/TwoVector.h"


namespace mu2e {


      SquareMapper::SquareMapper(void) : _step(),
                                         _apexX({-0.5,0.5,0.5,-0.5,-0.5}),
			                 _apexY({-0.5,-0.5,0.5,0.5,-0.5}) 
   
      {
         _step.push_back( SquLK( 1, 0) );  //right
	 _step.push_back( SquLK( 0,-1) ); //down
	 _step.push_back( SquLK(-1, 0) );  //left
	 _step.push_back( SquLK( 0, 1) ); //up 
      }
                       




      CLHEP::Hep2Vector SquareMapper::xyFromIndex(int thisIndex) const
      {        
         SquLK thisLK = lk(thisIndex);
         double x = thisLK._l;
         double y = thisLK._k;
	 return CLHEP::Hep2Vector(x,y);
      }
  
      int SquareMapper::indexFromXY(double x0, double y0) const
      {        
         int l = int( std::abs(x0)-0.5);
         int k = int( std::abs(y0)-0.5 );
	 if (x0<0) l *= -1;
	 if (y0<0) k *= -1;

	 SquLK lk(l,k);	 
	 return index(lk);
      }


  

      std::vector<int> SquareMapper::neighbors(int thisIndex, int level)  const
      {	 
	 std::vector<int> thisNeighbour;
	 thisNeighbour.reserve(100);

         SquLK init = lk(thisIndex);
	 SquLK lk(init._l - level, init._k + level);

         for (unsigned int i=0;i<_step.size();++i)
	 {       	     
	     for (int iseg=0;iseg<2*level;++iseg)
	     {	  
		lk += _step[i];  
		thisNeighbour.push_back( index(lk) );
	     }
	 }
         return thisNeighbour;
      }



      SquLK SquareMapper::lk(int thisIndex) const
      {         
	 if (thisIndex==0) return SquLK(0,0);

	 int nRing = int(0.5*sqrt(thisIndex) + 0.5);

         int nSeg  = (thisIndex - (2*nRing-1)*(2*nRing-1)) / (2*nRing);
         int nPos  = (thisIndex - (2*nRing-1)*(2*nRing-1)) % (2*nRing);
	 
	 if (nSeg==0) return SquLK( -nRing+nPos ,  nRing     ); 
	 if (nSeg==1) return SquLK( nRing       ,  nRing-nPos); 
	 if (nSeg==2) return SquLK( nRing-nPos  , -nRing     ); 	 
	              return SquLK(-nRing       , -nRing+nPos); 
      } 

      int SquareMapper::index(SquLK& thisLK) const
      {
         if (thisLK._l==0 && thisLK._k==0) return 0;
	 
	 int nRing = ring(thisLK);
	 int pos   = (2*nRing-1)*(2*nRing-1);

	 //add position along segment 
	 if ( thisLK._k ==  nRing && thisLK._l < nRing)   pos +=           nRing + thisLK._l;
	 if ( thisLK._l ==  nRing && thisLK._k > -nRing)  pos += 2*nRing + nRing - thisLK._k;
	 if ( thisLK._k == -nRing && thisLK._l > -nRing)  pos += 4*nRing + nRing - thisLK._l;
	 if ( thisLK._l == -nRing && thisLK._k < nRing)   pos += 6*nRing + nRing + thisLK._k;

	 return pos;      
      }

      int SquareMapper::ring(SquLK& thisLK) const
      {         
	  return std::max(std::abs(thisLK._l),std::abs(thisLK._k));
      }




}


             
