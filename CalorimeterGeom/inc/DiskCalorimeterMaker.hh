#ifndef CalorimeterGeom_DiskCalorimeterMaker_hh
#define CalorimeterGeom_DiskCalorimeterMaker_hh
//
// $Id: DiskCalorimeterMaker.hh,v 1.2 2013/03/05 20:33:25 aluca Exp $
// $Author: aluca $
// $Date: 2013/03/05 20:33:25 $
//
// original authors B. Echenard

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
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "CalorimeterGeom/inc/Disk.hh"
#include "ConfigTools/inc/SimpleConfig.hh"


namespace mu2e{


    class SimpleConfig;
    class DiskCalorimeter;

    class DiskCalorimeterMaker{

    public:
 
       DiskCalorimeterMaker(SimpleConfig const& config, double solenoidOffset);
      ~DiskCalorimeterMaker();

      // Accessor and auto_ptr to calorimeter needed by GeometryService.
      std::auto_ptr<DiskCalorimeter> _calo;
      std::auto_ptr<DiskCalorimeter> calorimeterPtr() { return _calo; }

    private:

      void MakeDisks(void);
      void CheckIt(void);

      int verbosityLevel;

    };

} //namespace mu2e

#endif /* CalorimeterGeom_DiskCalorimeterMaker_hh */
