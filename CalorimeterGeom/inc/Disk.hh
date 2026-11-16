#ifndef CalorimeterGeom_Disk_hh
#define CalorimeterGeom_Disk_hh
//
// Hold information about a disk in the calorimter.
//
// The crystal numbering scheme start at the center of the disk
// to facilitate navigation
//
// Original author B Echenard
//
#include "Offline/CalorimeterGeom/inc/DiskInfo.hh"
#include "Offline/CalorimeterGeom/inc/Crystal.hh"
#include "Offline/CalorimeterGeom/inc/CrystalMapper.hh"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/Rotation.h"
#include <vector>
#include <memory>
#include <iostream>
#include <string>

namespace mu2e {

  class Disk {
    public:
        Disk(unsigned id, double rCrystalIn, double rCrystalOut,
             double nominalCellSize, double nominalCellLength, unsigned offset,
             const CLHEP::Hep3Vector& diskOriginToCrystalOrigin,
             const std::string& crystalFileName);

        unsigned                        id()                const {return id_;}

        size_t                          nCrystals()         const {return crystals_.size();}
        const Crystal&                  crystal(unsigned i) const {return crystals_.at(i);}
              Crystal&                  crystal(unsigned i)       {return crystals_.at(i);}
        unsigned                        crystalOffset()     const {return crystalOffset_;}

        const DiskInfo&                 diskInfo()          const {return diskInfo_;}
              DiskInfo&                 diskInfo()                {return diskInfo_;}

        std::vector<unsigned>           neighbors      (unsigned localId, unsigned level=1) const;
        std::vector<unsigned>           idxFromRow     (int thisRow)                  const;
        unsigned                        idxFromPosition(const CLHEP::Hep3Vector& pos) const;
        bool                            isInsideDisk   (const CLHEP::Hep3Vector& pos) const;
        bool                            isInsideCrystal(const CLHEP::Hep3Vector& pos) const;

        CLHEP::Hep3Vector               mu2eToDisk     (const CLHEP::Hep3Vector& p) const;
        CLHEP::Hep3Vector               diskToMu2e     (const CLHEP::Hep3Vector& p) const;
        CLHEP::Hep3Vector               mu2eToDiskFF   (const CLHEP::Hep3Vector& p) const;
        CLHEP::Hep3Vector               diskFFToMu2e   (const CLHEP::Hep3Vector& p) const;
        CLHEP::Hep3Vector               mu2eToCrystal  (unsigned localId, const CLHEP::Hep3Vector& p) const;
        CLHEP::Hep3Vector               crystalToMu2e  (unsigned localId, const CLHEP::Hep3Vector& p) const;

        void                            moveCrystal    (unsigned localId, const CLHEP::Hep3Vector& disp);
        void                            moveDisk       (const CLHEP::Hep3Vector& disp,
                                                        const CLHEP::HepRotation& rotation);

        void                            print          (std::ostream& os = std::cout) const;


    private:
        void                            fillCrystals      (const CLHEP::Hep3Vector&, double nominalCellLength,
                                                           const std::string& filename);
        void                            checkPosition     (const CLHEP::Hep3Vector& pos, const CLHEP::Hep3Vector& size,
                                                           unsigned crystalID)                            const;
        bool                            isInsideCrystal   (unsigned icry, const CLHEP::Hep3Vector& pos)   const;
        bool                            isCrystalIdxValid (unsigned i)                                    const;

        unsigned                        id_;
        std::vector<Crystal>            crystals_;
        DiskInfo                        diskInfo_;
        double                          radiusInCrystal_;
        double                          radiusOutCrystal_;
        double                          nominalCellSize_;
        unsigned                        crystalOffset_;
        std::shared_ptr<CrystalMapper>  crystalMap_;
        std::vector<unsigned>           mapToCrystal_;
        std::vector<unsigned>           crystalToMap_;
  };
}
#endif
