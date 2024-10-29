#ifndef CalorimeterGeom_DiskCalorimeterMaker_hh
#define CalorimeterGeom_DiskCalorimeterMaker_hh
//
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
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "Offline/CalorimeterGeom/inc/Disk.hh"
#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/DataProducts/inc/CaloConst.hh"


namespace mu2e{


    class SimpleConfig;
    class DiskCalorimeter;

    class DiskCalorimeterMaker{

    public:

       DiskCalorimeterMaker(SimpleConfig const& config, double solenoidOffset);

      std::unique_ptr<DiskCalorimeter> calorimeterPtr() { return std::move(calo_); }

    private:

      // unique_ptr to calorimeter needed by GeometryService.
      std::unique_ptr<DiskCalorimeter> calo_;

      void CheckIt(void);
      void MakeIt(void);

      int verbosityLevel_;
      double FPHalfZLength_;
      double diskCaseHalfZLength_;
      double BPHalfZLength_;
      double diskHalfZLength_;
      double FEBHalfZLength_;
      double motherHalfZ_;
      double crateToDiskDeltaZ_;
      CLHEP::Hep3Vector diskOriginToCrystalOrigin_;


    };

}

#endif
