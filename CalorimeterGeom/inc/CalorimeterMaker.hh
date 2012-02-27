#ifndef CalorimeterGeom_CalorimeterMaker_hh
#define CalorimeterGeom_CalorimeterMaker_hh
// $Id: CalorimeterMaker.hh,v 1.14 2012/02/27 06:05:35 gandr Exp $
// $Author: gandr $
// $Date: 2012/02/27 06:05:35 $

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
#include "Mu2eUtilities/inc/SimpleConfig.hh"

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

      // Accessor and auto_ptr to calorimeter needed by GeometryService.
      std::auto_ptr<Calorimeter> _calo;
      std::auto_ptr<Calorimeter> getCalorimeterPtr() { return _calo; }

    private:

      void BuildIt();
      void MakeVanes();
      //void MakeCalorimeter();

      //void FillNearestNeighbours();

    };

} //namespace mu2e

#endif /* CalorimeterGeom_CalorimeterMaker_hh */
