//
// Interface to the disk calorimeter
// Original author B. Echenard
//

#ifndef CalorimeterGeom_Calorimeter_hh
#define CalorimeterGeom_Calorimeter_hh

#include "Offline/Mu2eInterfaces/inc/Detector.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/CalorimeterGeom/inc/CaloGeomUtil.hh"
#include "Offline/CalorimeterGeom/inc/CaloInfo.hh"
#include "Offline/CalorimeterGeom/inc/Disk.hh"
#include "Offline/CalorimeterGeom/inc/Crystal.hh"

#include "CLHEP/Vector/ThreeVector.h"
#include <vector>

namespace mu2e {

    class Calorimeter: virtual public Detector {

        public:

           //no constructor for this interface
           virtual ~Calorimeter() = default;

           virtual size_t                        nDisks()          const = 0;
           virtual const Disk&                   disk(size_t i)    const = 0;
           virtual const DiskPtrs&               diskPtrs()        const = 0;

           virtual size_t                        nCrystals()       const = 0;
           virtual const Crystal&                crystal(size_t i) const = 0;
           virtual const CrystalPtrs&            crystalPtrs()     const = 0;

           virtual const CaloInfo&               caloInfo()        const = 0;
           virtual const CaloGeomUtil&           geomUtil()        const = 0;

           virtual const std::vector<int>&       neighbors(int crystalId, bool rawMap=false)                     const = 0;
           virtual const std::vector<int>&       nextNeighbors(int crystalId, bool rawMap=false)                 const = 0;
           virtual       std::vector<int>        neighborsByLevel(int crystalId, int level, bool rawMap = false) const = 0;
           virtual int                           crystalIdxFromPosition(const CLHEP::Hep3Vector& pos)            const = 0;
           virtual int                           nearestIdxFromPosition(const CLHEP::Hep3Vector& pos)            const = 0;

           virtual void                          print(std::ostream &os = std::cout)                             const = 0;
    };
}

#endif
