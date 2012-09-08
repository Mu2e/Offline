// Hexagon position map generator, tesselate a disk with hexagons
// Basis vector, k and l, defined as  
// l = up right
// k = down right
//
// so 
// l-k = up
// l+k = right 
/*
           ____
          /    \ 
     ____/ 1 -1 \ ____
    /    \      /    \
   / 0 -1 \____/  1 0 \   l,k coordinates
   \      /    \      /  (for whatever reason, l seems to come before k...)
    \____/  0 0 \____/   (bonus points for drawing this btw...)
    /    \      /    \
   / -1 0 \____/ 0 1  \
   \      /    \      /
    \____/ -1 1 \____/
         \      /    
          \____/      
*/
// For each ring, start at n+,-n (one step down right from the upmost hexagon), 
// then go n time each step to create the ring.
//
// Neighbors of (l0,k0), just add +-(1,0), +-(0,1) or +-(1,-1) to (l0,k0)
//                       next ring of neighbours, add +(2,2) and go around the ring,...
//
// Hexagon dimensions: _hexsize refers to distance between flats sides of hexagon
//                     !!!!!N.B: in mu2e, the hexagon = crystal+wrapper+shell, not only the crystal!
//
// Disk dimensions:    _radiusIn and _radiusOut are the inner and outer radius of the disks
//


// C++ includes
#include <iostream>
#include <map>
#include <cassert>
#include "CLHEP/Vector/TwoVector.h"


// Mu2e includes
#include "CalorimeterGeom/inc/HexPosition.hh"
#include "CalorimeterGeom/inc/HexPositionMap.hh"


namespace mu2e {


      HexPositionMap::HexPositionMap(double hexsize, double radiusIn, double radiusOut) : 
      _step(), _positions(), _lkToIdx(), _hexsize(hexsize), _radiusIn(radiusIn), _radiusOut(radiusOut)
      {

	  _step.clear();           
	  _step.push_back( std::pair<int,int>(0,1)  ); 
	  _step.push_back( std::pair<int,int>(-1,1) ); 
	  _step.push_back( std::pair<int,int>(-1,0) ); 
	  _step.push_back( std::pair<int,int>(0,-1) ); 
	  _step.push_back( std::pair<int,int>(1,-1) ); 
	  _step.push_back( std::pair<int,int>(1,0)  ); 

	  generate();
      }

      
      
      
      void HexPositionMap::generate(void)
      {

	 int nring = int(1.5*_radiusOut/_hexsize)+1;
	 int ncrys = 1+nring*(nring-1)/2*6;

	 _positions.reserve(ncrys);
	 _positions.clear();
	 _lkToIdx.clear();

	 if (_radiusIn<0.001) {
	    HexPosition center(0,0,0);
	    _positions.push_back(center);
	    _lkToIdx[0] = 0;
	 }

	 for(int iring=1;iring<nring;iring++){  

	   std::pair<int,int> cp(iring,-iring); 

	   for(unsigned int istep =0;istep<_step.size();istep++){

	     for(int iseg=0;iseg<iring;iseg++){
        	cp.first  += _step[istep].first;
        	cp.second += _step[istep].second;

        	int thisId = getId(cp.first,cp.second);
		HexPosition thisPos(cp.first,cp.second,_hexsize);
		
		if (thisPos.getDistMin()<_radiusIn || thisPos.getDistMax()>_radiusOut) continue;
		
		_positions.push_back(thisPos);
        	_lkToIdx[thisId] = _positions.size()-1;

	     }
	   }//end step
	 }//end ring

	 return;

      } 


      
      // get neighbours for or crystal, include only crystals in the disk
      // level indicate how many rings around the crystal to consider: 
      // 1=immediate surrounding, 2 = next-to-immediate,...  (see top picture)      
      std::vector<int> HexPositionMap::getNeighbors(unsigned int id, unsigned int level) const
      {
	 std::vector<int> neighbors;
	 neighbors.reserve(100);

	 if (level < 0 || id > _positions.size()) return neighbors;

	 std::pair<int,int> lkPosition(_positions[id].getl()+level,_positions[id].getk()-level);
	 for (unsigned int i=0;i<_step.size();++i)
	 {       
	     for (unsigned int iseg=0;iseg<level;++iseg)
	     {	  
		lkPosition.first  += _step[i].first;
		lkPosition.second += _step[i].second;

		std::map<int,int>::const_iterator it = _lkToIdx.find(getId(lkPosition.first,lkPosition.second));
		if (it == _lkToIdx.end()) continue; //take only neighbours inside the disk
		neighbors.push_back(it->second);
	     }
	 }

	 return neighbors;
      }

      
      
      
      
      int HexPositionMap::findIdFromPosition(double x, double y) const
      {

	  double ld = (1.15470053837*x+2.0*y)/2.0/_hexsize;
	  double kd = (1.15470053837*x-2.0*y)/2.0/_hexsize;

	  int l = (ld > 0) ? int(ld+0.5) : int(ld-0.5);
	  int k = (kd > 0) ? int(kd+0.5) : int(kd-0.5);

	  std::map<int,int>::const_iterator it = _lkToIdx.find(getId(l,k));
	  if (it != _lkToIdx.end()) return it->second-1;
	  return -1;

      }
      
      
      int HexPositionMap::getId(int l, int k) const
      {
	 int id = 100*std::abs(l)+std::abs(k);
	 if (k<0) id +=10000;
	 if (l<0) id +=20000;
	 return id;
      }


}

