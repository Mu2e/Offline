#ifndef CalorimeterGeom_CaloSection_hh
#define CalorimeterGeom_CaloSection_hh
//
// Hold information about a disk in the calorimter. Note that there are two radii. The physics dimensions of the disk, including the
// inner / outer rings (radiusIn_/radiusOut_), and the radii inside which the crystals are placed (radiusInCrystal_/radiusOutCrystal_).
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
           Disk(int id, double rin, double rout, double rCrystalIn, double rCrystalOut, double nominalCellSize,
                int offset, const CLHEP::Hep3Vector& diskOriginToCrystalOrigin);

           int                             id()              const {return id_;}
           size_t                          nCrystals()       const {return crystalList_.size();}
           int                             crystalOffset()   const {return globalCrystalOffset_;}
           const Crystal&                  crystal(int i)    const {return crystalList_.at(i);}
                 Crystal&                  crystal(int i)          {return crystalList_.at(i);}

           const DiskGeomInfo&             geomInfo()        const {return geomInfo_;}
                 DiskGeomInfo&             geomInfo()              {return geomInfo_;}

           double                          innerRadius()     const {return radiusIn_;}
           double                          outerRadius()     const {return radiusOut_;}

           int                             idxFromPosition       (double x, double y) const;
           std::vector<int>                findLocalNeighbors    (int crystalId, int level, bool raw=false) const;
           std::vector<int>                nearestIdxFromPosition(double x, double y) const;
           int                             idMinCrystalInside    (int row);
           int                             idMaxCrystalInside    (int row);
           void                            boundingBoxes         (int thisRow,std::vector<double>& params);

           void                            print(std::ostream& os = std::cout) const;


       private:
           void                            fillCrystalsIdeal (const CLHEP::Hep3Vector &crystalOriginInDisk);
           void                            fillCrystals      (const CLHEP::Hep3Vector &crystalOriginInDisk);
           bool                            isInsideDisk      (double x, double y, double widthX, double widthY) const;
           bool                            isInsideCrystal   (int icry, double x, double y) const;
           void                            fixCrystalPosition();

           int                             id_;
           std::vector<Crystal>            crystalList_;
           DiskGeomInfo                    geomInfo_;
           double                          radiusIn_;
           double                          radiusOut_;
           double                          radiusInCrystal_;
           double                          radiusOutCrystal_;
           double                          nominalCellSize_;
           int                             globalCrystalOffset_;
           std::shared_ptr<CrystalMapper>  crystalMap_;
           std::vector<int>                mapToCrystal_;
           std::vector<int>                crystalToMap_;

   };

   using DiskPtr  = std::shared_ptr<Disk>;
   using DiskPtrs = std::vector<DiskPtr>;
}

#endif
