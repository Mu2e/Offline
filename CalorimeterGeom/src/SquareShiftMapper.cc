// $Id: SquareShiftMapper.cc,v 1.4 2013/07/25 23:56:46 echenard Exp $
// $Author: echenard $
// $Date: 2013/07/25 23:56:46 $
//
// Sqaure position map generator: 
//   tesselate a plane with squares, every row shifted horizontaly by 0.5 square size, starting from the center of the plane
//
// Use basis vector, l and k, defined as  
// l = up right
// k = down right
//
/*

       --------------------
       |         |        |
       |  0 -1   |   1 0  |
       |         |        |
       |         |        |
 ------------------------------
 |         |         |        |  
 |  -1 -1  |   0 0   |  1 1   |   l,k coordinates
 |         |         |        |
 |         |         |        |
 ------------------------------
       |         |        |
       |  -1 0   |  0 1   |
       |         |        |
       |         |        |
       --------------------

  steps :  (1,1) (0,1) (-1,0) (-1,-1) (0,-1) (1,0) (clockwise from top left corner)

*/
// Tesselation algorithm: tessalate in "rings" from the center
//   for each ring, start at 0,-l (top left corner), 
//   then go n time each step to create the ring 
//
// Neighbors add (0,-1) and go around the ring
// next ring of neighbours, add (0,-2) and go around the ring,...
//

#include "CalorimeterGeom/inc/SquareShiftMapper.hh"

#include "CLHEP/Vector/TwoVector.h"

#include <iostream>
#include <map>
#include <cmath>



namespace mu2e {


      SquareShiftMapper::SquareShiftMapper() : 
        _step(),
        _apexX({-0.5,0.5,0.5,-0.5,-0.5}),
        _apexY({-0.5,-0.5,0.5,0.5,-0.5})    
      {
          _step.push_back( SquShiftLK(  1,  1) ); //right
	  _step.push_back( SquShiftLK(  0 , 1) ); //down right
	  _step.push_back( SquShiftLK( -1,  0) ); //down left
	  _step.push_back( SquShiftLK( -1, -1) ); //left
	  _step.push_back( SquShiftLK(  0, -1) ); //up left 
	  _step.push_back( SquShiftLK(  1,  0) ); //up right
      }
                       




      CLHEP::Hep2Vector SquareShiftMapper::xyFromIndex(int thisIndex) const
      {        
          SquShiftLK thisLK = lk(thisIndex);
	  return CLHEP::Hep2Vector( (thisLK._l+thisLK._k)/2.0, (thisLK._l-thisLK._k) );
      }
  
      int SquareShiftMapper::indexFromXY(double x0, double y0) const
      {        
	  int l,k;
	  int ny = (y0>0) ? int(std::abs(y0)+0.5) : -int(std::abs(y0)+0.5);

	  if (ny%2==0)
	  {
	      int nx = (x0>0) ? int(std::abs(x0)+0.5) : -int(std::abs(x0)+0.5);
	      l = nx+ny/2;
	      k = nx-ny/2;
	  } else {
	      int nx = (x0>0) ? int(std::abs(x0)) : (-int(std::abs(x0))-1);
	      l = nx + (ny+1)/2;
	      k = nx - (ny-1)/2;
	  }

	  SquShiftLK lk(l,k);	 
	  return index(lk);
      }

  

      std::vector<int> SquareShiftMapper::neighbors(int thisIndex, unsigned int level)  const
      {	 
	  std::vector<int> thisNeighbour;
	  thisNeighbour.reserve(100);

	  SquShiftLK init = lk(thisIndex);
	  SquShiftLK lk(init._l, init._k - level);

	  for (unsigned int i=0;i<_step.size();++i)
	  {       	     
	      for (unsigned int iseg=0;iseg<level;++iseg)
	      {	  
		 lk.add(_step[i]);  
		 thisNeighbour.push_back( index(lk) );
	      }
	  }
	  return thisNeighbour;
      }


      SquShiftLK SquareShiftMapper::lk(int thisIndex) const
      {         
	 if (thisIndex==0) return SquShiftLK(0,0);

	 int nRing = int(0.5+sqrt(0.25+float(thisIndex-1)/3.0));
	 int nSeg  = (thisIndex - 3*nRing*(nRing-1)-1) / nRing;
	 int nPos  = (thisIndex - 3*nRing*(nRing-1)-1) % nRing;

	 int l =          nPos*_step[nSeg]._l;
	 int k = -nRing + nPos*_step[nSeg]._k;
	 	 
	 for (int i=0;i<nSeg;++i)
	 {
	    l += _step[i]._l*nRing;
	    k += _step[i]._k*nRing;
	 }

	 return SquShiftLK(l,k);
      } 


      int SquareShiftMapper::index(SquShiftLK const &thisLK) const
      {
	 if (thisLK._l==0 && thisLK._k==0) return 0;

	 int nRing = ring(thisLK);
	 int pos   = 3*nRing*(nRing-1)+1;
	  
	 //add position along segment -- a drawing helps...
	 if (thisLK._k == -nRing && thisLK._l==0) return pos;
	 if (thisLK._l ==  nRing)                 return pos + nRing   + thisLK._k;
	 if (thisLK._l == -nRing)                 return pos + 4*nRing + std::abs(thisLK._k);
	 if (thisLK._k ==  nRing)                 return pos + 3*nRing - thisLK._l;
	 if (thisLK._k == -nRing)                 return pos + 6*nRing - std::abs(thisLK._l);
	 if (thisLK._l > thisLK._k)               return pos + thisLK._l;
                                                  return pos + 3*nRing +std::abs(thisLK._l);
      }


      int SquareShiftMapper::ring(const SquShiftLK &thisLK) const
      {         
	  if (thisLK._l*thisLK._k>0) return std::max(std::abs(thisLK._l),std::abs(thisLK._k));
	  return std::abs(thisLK._l-thisLK._k);
      }

}


             
