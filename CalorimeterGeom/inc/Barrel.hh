#ifndef CalorimeterGeom_Barrel_hh
#define CalorimeterGeom_Barrel_hh
// $Id: Barrel.hh,v 1.1 2013/09/05 17:11:28 gianipez Exp $
// $Author: gianipez $
// $Date: 2013/09/05 17:11:28 $
//
// Hold information about a Gianipez in the calorimter.
//
// Original author Gianipez
//

// C++ includes
#include <vector>

// Mu2e includes
#include "CalorimeterGeom/inc/CaloSection.hh"
#include "CalorimeterGeom/inc/Crystal.hh"

//CLHEP includess
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/ThreeVector.h"


namespace mu2e {


  class Barrel : public CaloSection {


  public:

    Barrel(int id
	   , double rin, double rout
	   , double cellSize
	   , int nWheels, int nCrystalWheel
	   , CLHEP::Hep3Vector crystalOffset); 

    std::vector<int> neighbors(int crystalId, int level=1) const; 
    int              idxFromPosition(CLHEP::Hep3Vector pos) const;

    double           innerRadius()   const {return _radiusIn;}
    double           outerRadius()   const {return _radiusOut;}
    
    int              nCrystalWheel() const {return _nCrystalWheel;}
    int              nWheels()       const {return _nWheels;}
	   
  private:

    void             fillCrystals();
    bool             isInsideBarrel(CLHEP::Hep3Vector const& pos) const;
    std::vector<int> findNeighbors(int crystalId, int level=1) const; 


    double             _radiusIn;
    double             _radiusOut;
    double             _cellSize;	 
    int                _nWheels;
    int                _nCrystalWheel;
  };
}

#endif /* CalorimeterGeom_Barrel_hh */
