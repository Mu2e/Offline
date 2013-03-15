#ifndef CalorimeterGeom_CalorimeterMaker_hh
#define CalorimeterGeom_CalorimeterMaker_hh
// $Id: CalorimeterMaker.hh,v 1.17 2013/03/15 15:52:03 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/15 15:52:03 $

// original authors Julie Managan and Robert Bernstein

//
// C++ includes
#include <iomanip>
#include <memory>
#include <iostream>
#include <cmath>
#include <vector>
#include <string>

//
//Mu2e includes
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/Vane.hh"
#include "ConfigTools/inc/SimpleConfig.hh"

//
// other includes
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e{

  class SimpleConfig;

    //
    //forward declarations
    class Calorimeter;

    class CalorimeterMaker{

    public:
      CalorimeterMaker(SimpleConfig const& config, double solenoidOffset);

      ~CalorimeterMaker();

      //
      // simple dumb vector

      // Accessor and unique_ptr to calorimeter needed by GeometryService.
      std::unique_ptr<Calorimeter> _calo;
      std::unique_ptr<Calorimeter> calorimeterPtr() { return std::move(_calo); }

    private:

      void BuildIt();
      void MakeVanes();
      //void MakeCalorimeter();

      //void FillNearestNeighbours();

    };

} //namespace mu2e

#endif /* CalorimeterGeom_CalorimeterMaker_hh */
