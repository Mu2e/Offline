#ifndef CalorimeterGeom_HybridCalorimeterMaker_hh
#define CalorimeterGeom_HybridCalorimeterMaker_hh
//
// $Id: HybridCalorimeterMaker.hh,v 1.1 2013/09/05 17:11:28 gianipez Exp $
// $Author: gianipez $
// $Date: 2013/09/05 17:11:28 $
//
// original authors G. Pezzullo

//
// C++ includes
#include <iomanip>
#include <memory>
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include "CLHEP/Vector/ThreeVector.h"

//
//Mu2e includes
#include "CalorimeterGeom/inc/HybridCalorimeter.hh"
#include "CalorimeterGeom/inc/Disk.hh"
#include "ConfigTools/inc/SimpleConfig.hh"


namespace mu2e{


    class SimpleConfig;
    class HybridCalorimeter;

    class HybridCalorimeterMaker{

    public:
 
       HybridCalorimeterMaker(SimpleConfig const& config, double solenoidOffset);
      ~HybridCalorimeterMaker();

      // Accessor and unique_ptr to calorimeter needed by GeometryService.
      std::unique_ptr<HybridCalorimeter> _calo;
      std::unique_ptr<HybridCalorimeter> calorimeterPtr() { return std::move(_calo); }

    private:

      void MakeDisk(void);
      void MakeBarrel(void);
      void CheckIt(void);

      int _verbosityLevel;

    };

} //namespace mu2e

#endif /* CalorimeterGeom_HybridCalorimeterMaker_hh */
