#ifndef CalorimeterGeom_HexPositionMap_hh
#define CalorimeterGeom_HexPositionMap_hh
//
// $Id: HexPositionMap.hh,v 1.2 2013/03/05 20:33:25 aluca Exp $
// $Author: aluca $
// $Date: 2013/03/05 20:33:25 $
//
// Hexagonal map generator for disks
//
// Original author B Echenard 
//

// C++ includes
#include <vector>
#include <map>

// Mu2e includes
#include "CalorimeterGeom/inc/HexPosition.hh"
#include "CLHEP/Vector/TwoVector.h"


namespace mu2e {


     class HexPositionMap {


       public:

	 HexPositionMap(double scale, double radiusIn, double radiusOut);
	 virtual ~HexPositionMap(){}


	 int nCrystals(void) const                       {return _positions.size();}
	 HexPosition const& crystalPos(int idx) const {return _positions.at(idx);}

	 int  findIdFromPosition(double x, double y) const;
	 std::vector<int> neighbors(unsigned int id, unsigned int level = 1) const;


       private:

	 std::vector< std::pair<int,int> > _step; 
	 std::vector<HexPosition>          _positions;
	 std::map<int,int>                 _lkToIdx;
	 double                            _hexsize;
	 double                            _radiusIn;
	 double                            _radiusOut;

	 void generate(void);    
	 int  id(int l, int k) const;
     };


} //namespace mu2e

#endif
