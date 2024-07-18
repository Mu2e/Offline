#ifndef CalorimeterGeom_DiskCalorimeterMaker_hh
#define CalorimeterGeom_DiskCalorimeterMaker_hh
//
// original authors B. Echenard
//
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "Offline/CalorimeterGeom/inc/Disk.hh"
#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/DataProducts/inc/CaloConst.hh"
#include "CLHEP/Vector/ThreeVector.h"


namespace mu2e{

    class SimpleConfig;
    class DiskCalorimeter;

    class DiskCalorimeterMaker{

      public:
         DiskCalorimeterMaker(SimpleConfig const& config, double solenoidOffset);

        // Accessor and unique_ptr to calorimeter needed by GeometryService.
        std::unique_ptr<DiskCalorimeter> calorimeterPtr() { return std::move(calo_); }

      private:
        void   checkIt();
        void   makeIt();

        int    verbosityLevel_;
        std::unique_ptr<DiskCalorimeter> calo_;
    };

}

#endif
