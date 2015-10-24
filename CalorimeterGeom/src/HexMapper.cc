// $Id: HexMapper.cc,v 1.4 2013/07/25 23:56:46 echenard Exp $
// $Author: echenard $
// $Date: 2013/07/25 23:56:46 $
//
// Hexagon position map generator: 
//   tesselate a plane with hexagons starting at the center of the plane
//
// Use basis vector, l and k, defined as  
// l = up right
// k = down right
//
// so 
// l-k = up
// l+k = right 
/*
           ____
          /    \ 
     ____/ 1 -1 \____
    /    \      /    \
   / 0 -1 \____/  1 0 \   l,k coordinates
   \      /    \      /  (for whatever reason, l seems to come before k...)
    \____/  0 0 \____/   (hope you appreciate this piece of art...)
    /    \      /    \
   / -1 0 \____/ 0 1  \
   \      /    \      /
    \____/ -1 1 \____/
         \      /    
          \____/      
*/
// Tesselation algorithm: tessalate in "rings" from the center
//   for each ring, start at n+,-n (one step down right from the upmost hexagon), 
//   then go n time each step to create the ring (see example function at the end of the file)
//
// Neighbors of (l0,k0), just add +-(1,0), +-(0,1) or +-(1,-1) to (l0,k0)
//                       next ring of neighbours, add +(2,2) and go around the ring,...
//


// C++ / CLHEP includes
#include <iostream>
#include <map>
#include <cmath>
#include "CLHEP/Vector/TwoVector.h"

// Mu2e includes
#include "CalorimeterGeom/inc/HexMapper.hh"



namespace mu2e {


      HexMapper::HexMapper(void) : _step(),
                                   _apexX( {-0.2886751,+0.2886751,+0.5773502,+0.2886751,-0.2886751,-0.5773502,-0.2886751} ),
			           _apexY( {-0.5,-0.5,0,0.5,0.5,0,-0.5} ) 
      {
          _step.push_back( HexLK( 0, 1) ); //down right
	  _step.push_back( HexLK(-1, 1) ); //down
	  _step.push_back( HexLK(-1, 0) ); //down left
	  _step.push_back( HexLK( 0,-1) ); //up left
	  _step.push_back( HexLK( 1,-1) ); //up
	  _step.push_back( HexLK( 1, 0) ); //up right      
      }
                       




      CLHEP::Hep2Vector HexMapper::xyFromIndex(int thisIndex) const
      {        
          HexLK thisLK = lk(thisIndex);
          double x = (thisLK._l+thisLK._k)*sqrt(3.0)/2.0;
          double y = (thisLK._l-thisLK._k)/2.0;
	  return CLHEP::Hep2Vector(x,y);
      }
  

      int HexMapper::indexFromXY(double x, double y) const
      {        
	  int nl = int(x/sqrt(3)) + int(y);
	  int nk = int(x/sqrt(3)) - int(y);

	  double x0 = fmod(x,sqrt(3))*sqrt(3);  //sqrt(3) for convenience
	  double y0 = fmod(y,1.0);

	  if (x0>2.9999999) x0=2.99999999;

	  if (x0<0) { --nl;--nk; x0+=3.0;}
	  if (y0<0) { --nl;++nk; y0+=1.0;}

	  ++nl;
	  switch (int(2*x0))
	  {
	       case 0:        
        	   --nl; --nk;
        	   if (y0 > 0.5) ++nl; 
        	   else          ++nk;
		   break;

	       case 1:
		   if (y0 < 0.5 && y0 < (1.0-x0) ) --nl;
		   if (y0 > 0.5 && y0 > x0) --nk;
        	   break;

	       case 2: break;

	       case 3: break;

	       case 4:
		  if (y0<0.5 && y0 < (x0-2.0)) ++nk;
		  if (y0>0.5 && y0 > (3.0-x0)) ++nl;
		  break;

	       case 5:
        	   if (y0>0.5) ++nl; 
        	   else        ++nk;
		   break;   

	      default: std::cout<<"Weird behaviour in HexMapper::indexFromXY for case "<<
	                   int(2*x0)<<"  x0,y0="<<x0<<" "<<y0<<"   nl,nk="<<nl<<" "<<nk<<std::endl;
	  }    

	  HexLK lk(nl,nk);	 
	  return index(lk);
      }


      std::vector<int> HexMapper::neighbors(int thisIndex, int level)  const
      {	 
	  std::vector<int> thisNeighbour;
	  thisNeighbour.reserve(100);

          HexLK init = lk(thisIndex);
	  HexLK lk(init._l + level, init._k - level);

          for (unsigned int i=0;i<_step.size();++i)
	  {       	     
	      for (int iseg=0;iseg<level;++iseg)
	      {	  
		 lk += _step[i];
		 thisNeighbour.push_back( index(lk) );
	      }
	  }

          return thisNeighbour;
      }


      HexLK HexMapper::lk(int thisIndex) const
      {         
	  if (thisIndex==0) return HexLK(0,0);

	  int nRing = int(0.5+sqrt(0.25+(float(thisIndex)-1.0)/3.0));
          int nSeg  = (thisIndex -1 -3*nRing*(nRing-1))/nRing;
          int nPos  = (thisIndex -1 -3*nRing*(nRing-1))%nRing;

	  int l = nRing+(nPos+1)*_step[nSeg]._l;
	  int k = -nRing+(nPos+1)*_step[nSeg]._k;


	  //pre-calculate this
	  
	  for (int i=0;i<nSeg;++i) {
	     l += _step[i]._l*nRing;
	     k += _step[i]._k*nRing;
	  }

	  return HexLK(l,k);
      } 


      int HexMapper::index(HexLK& thisLK) const
      {
          if (thisLK._l==0 && thisLK._k==0) return 0;

	  int nring = ring(thisLK);
	  int pos = (nring>0) ? 1+3*nring*(nring-1): 0;

	  //find segment along the ring
	  int segment(0);	 
	  if ( std::abs(thisLK._l+thisLK._k) == nring && thisLK._k!=0) segment += 1;
	  if ( std::abs(thisLK._k) == nring && thisLK._l!=0)           segment += 2;
	  if ( (thisLK._l+thisLK._k) <=0 && thisLK._k<nring)           segment += 3;
	  pos += segment*nring;

	  //add position along segment	 
	  if (segment==0 || segment==3)  pos += nring - std::abs(thisLK._k)-1;
	  if (segment==1 || segment==4)  pos += nring - std::abs(thisLK._l)-1;
	  if (segment==2 || segment==5)  pos += std::abs(thisLK._l)-1;

	  return pos;      
      }
           

      int HexMapper::ring(HexLK& thisLK) const
      {         
	  if (thisLK._l*thisLK._k > 0)                     return std::abs(thisLK._l+thisLK._k);
	  if ( std::abs(thisLK._l) > std::abs(thisLK._k) ) return std::abs(thisLK._l);
	  return std::abs(thisLK._k);
      }




}


             
