#ifndef CalorimeterGeom_DiskCalorimeter_hh
#define CalorimeterGeom_DiskCalorimeter_hh
//
// Hold informations about the disk calorimeter
//
// Original author B. Echenard
//
#include "Offline/Mu2eInterfaces/inc/ProditionsEntity.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/CalorimeterGeom/inc/CaloG4Info.hh"
#include "Offline/CalorimeterGeom/inc/Disk.hh"
#include "Offline/CalorimeterGeom/inc/Crystal.hh"
#include "cetlib_except/exception.h"

#include "CLHEP/Vector/ThreeVector.h"
#include <vector>
#include <memory>
#include <iostream>


namespace mu2e {



  class DiskCalorimeter: public Calorimeter, public ProditionsEntity {

    friend class DiskCalorimeterMaker;

    public:
      using Hep3Vector = CLHEP::Hep3Vector;

      DiskCalorimeter();
      ~DiskCalorimeter() = default;
      DiskCalorimeter(const DiskCalorimeter& rhs);
      DiskCalorimeter(DiskCalorimeter&& rhs) noexcept ;
      DiskCalorimeter& operator=(const DiskCalorimeter& rhs) = delete;
      DiskCalorimeter& operator=(DiskCalorimeter&& rhs)      = delete;

      size_t            nDisks()            const override {return disks_.size();}
      const Disk&       disk(unsigned i)    const override {return disks_.at(i);}
      const Disks&      disks()             const override {return disks_;}
            Disks&      disks()                   override {return disks_;}

      size_t            nCrystals()         const override {return crystals_.size();}
      const Crystal&    crystal(unsigned i) const override {return *crystals_.at(i);}
      const Crystals&   crystals()          const override {return crystals_;}

      const CaloG4Info& G4Info()            const override {return G4Info_;}
      void              print(std::ostream &os = std::cout) const override;

      std::vector<unsigned> neighbors(unsigned globalId, unsigned level) const override;
      bool              isInsideAnyCrystal(const Hep3Vector& pos)    const override;
      bool              isInsideAnyDisk   (const Hep3Vector& pos)    const override;
      bool              isInsideSameDisk  (const Hep3Vector& front,
                                           const Hep3Vector& back)   const override;

      Hep3Vector        mu2eToDisk   (unsigned diskId,    const Hep3Vector& p) const override;
      Hep3Vector        diskToMu2e   (unsigned diskId,    const Hep3Vector& p) const override;
      Hep3Vector        mu2eToDiskFF (unsigned diskId,    const Hep3Vector& p) const override;
      Hep3Vector        diskFFToMu2e (unsigned diskId,    const Hep3Vector& p) const override;
      Hep3Vector        mu2eToCrystal(unsigned crystalId, const Hep3Vector& p) const override;
      Hep3Vector        crystalToMu2e(unsigned crystalId, const Hep3Vector& p) const override;

      // These should REALLY go to the tracker
      Hep3Vector        mu2eToTracker(const Hep3Vector& p) const override;
      Hep3Vector        trackerToMu2e(const Hep3Vector& p) const override;

    private:
      void rebuildCrystalPtrs();

      Disks             disks_;
      Crystals          crystals_; //non-owning cache for convenience
      CaloG4Info        G4Info_;
      CLHEP::Hep3Vector trackerCenter_;

  };
}
#endif
