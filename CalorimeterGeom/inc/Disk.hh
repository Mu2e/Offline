#ifndef CalorimeterGeom_CaloSection_hh
#define CalorimeterGeom_CaloSection_hh
//
// Hold information about a disk in the calorimter.
//
// Original author B Echenard
//

#include "Offline/CalorimeterGeom/inc/DiskGeomInfo.hh"
#include "Offline/CalorimeterGeom/inc/Crystal.hh"
#include "Offline/CalorimeterGeom/inc/CrystalMapper.hh"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/Rotation.h"
#include <vector>
#include <memory>

namespace mu2e {

   class Disk {

       public:
           Disk(int id, double rCrystalIn, double rCrystalOut, double nominalCellSize,
                double nominalCellLength, int offset, const CLHEP::Hep3Vector& diskOriginToCrystalOrigin);

           int                             id()              const {return id_;}

           size_t                          nCrystals()       const {return crystalList_.size();}
           int                             crystalOffset()   const {return globalCrystalOffset_;}
           const Crystal&                  crystal(int i)    const {return crystalList_.at(i);}
                 Crystal&                  crystal(int i)          {return crystalList_.at(i);}

           const DiskGeomInfo&             geomInfo()        const {return geomInfo_;}
                 DiskGeomInfo&             geomInfo()              {return geomInfo_;}

           int                             idxFromPosition        (double x, double y)       const;
           std::vector<int>                nearestNeighborsFromPos(double x, double y)       const;
           std::vector<int>                findLocalNeighbors     (int crystalId, int level) const;
           int                             idMinCrystalInside     (int row)                  const;
           int                             idMaxCrystalInside     (int row)                  const;
           void                            boundingBoxes          (int thisRow,std::vector<double>& params) const;

           void                            print(std::ostream& os = std::cout) const;

       private:
           void                            fillCrystalsIdeal (const CLHEP::Hep3Vector& crystalOriginInDisk,
                                                              double nominalCellLength);
           void                            fillCrystals      (const CLHEP::Hep3Vector& crystalOriginInDisk);
           void                            fixCrystalPosition();
           bool                            isInsideDisk      (double x, double y, double widthX, double widthY) const;
           bool                            isInsideCrystal   (int icry, double x, double y) const;
           const bool                      isCrystalIdxValid (int i) const;
           const bool                      isMapIdxValid     (int i) const;

           int                             id_;
           std::vector<Crystal>            crystalList_;
           DiskGeomInfo                    geomInfo_;
           double                          radiusInCrystal_;
           double                          radiusOutCrystal_;
           double                          nominalCellSize_;
           int                             globalCrystalOffset_;
           std::shared_ptr<CrystalMapper>  crystalMap_;
           std::vector<int>                mapToCrystal_;
           std::vector<int>                crystalToMap_;
           int                             rowMax_;

           static constexpr int            invalidID_ = -1;
   };
}

#endif
