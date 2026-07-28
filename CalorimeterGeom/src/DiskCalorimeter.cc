#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "Offline/CalorimeterGeom/inc/Disk.hh"
#include "Offline/Mu2eInterfaces/inc/ProditionsEntity.hh"

#include "CLHEP/Vector/ThreeVector.h"

#include <iostream>
#include <algorithm>


namespace mu2e {

  DiskCalorimeter::DiskCalorimeter() :
    ProditionsEntity(Calorimeter::cxname),
    disks_(),
    crystals_(),
    G4Info_()
  {
    rebuildCrystalPtrs();
  }

  DiskCalorimeter::DiskCalorimeter(const DiskCalorimeter& rhs) :
    ProditionsEntity(rhs),
    disks_(rhs.disks_),
    G4Info_(rhs.G4Info_)
  {
    rebuildCrystalPtrs();
  }

  DiskCalorimeter::DiskCalorimeter(DiskCalorimeter&& rhs) noexcept :
    ProditionsEntity(std::move(rhs)),
    disks_(std::move(rhs.disks_)),
    G4Info_(std::move(rhs.G4Info_))
  {
    rebuildCrystalPtrs();
  }

  void DiskCalorimeter::rebuildCrystalPtrs()
  {
    crystals_.clear();
    for (const auto& disk : disks_) {
      for (size_t i = 0; i < disk.nCrystals(); ++i) {
        crystals_.push_back(&disk.crystal(i));
      }
    }
  }

  std::vector<unsigned> DiskCalorimeter::neighbors(unsigned globalId, unsigned level) const
  {
    const Crystal& c = crystal(globalId);
    const Disk&    d = disk(c.diskID());
    auto ids = d.neighbors(c.localID(), level);
    for (auto& v : ids) v += d.crystalOffset();
    return ids;
  }

  bool DiskCalorimeter::isInsideAnyDisk(const CLHEP::Hep3Vector& pos) const
  {
    for (const auto& disk : disks_) {
      CLHEP::Hep3Vector posInDisk = mu2eToDisk(disk.id(),pos);
      if (disk.isInsideDisk(posInDisk)) return true;
    }
    return false;
  }

  bool DiskCalorimeter::isInsideAnyCrystal(const CLHEP::Hep3Vector& pos) const
  {
    for (const auto& disk : disks_) {
      CLHEP::Hep3Vector posInDisk = mu2eToDiskFF(disk.id(),pos);
      if (disk.isInsideCrystal(posInDisk)) return true;
    }
    return false;
  }

  bool DiskCalorimeter::isInsideSameDisk(const CLHEP::Hep3Vector& front,
                                         const CLHEP::Hep3Vector& back) const
  {
    for (const auto& disk : disks_) {
      CLHEP::Hep3Vector frontInDisk = mu2eToDisk(disk.id(),front);
      CLHEP::Hep3Vector backInDisk  = mu2eToDisk(disk.id(),back);
      if (disk.isInsideDisk(frontInDisk) && disk.isInsideDisk(backInDisk)) return true;
    }
    return false;
  }

  CLHEP::Hep3Vector DiskCalorimeter::mu2eToDisk  (unsigned diskId, const CLHEP::Hep3Vector& p) const {
    return disk(diskId).mu2eToDisk(p);
  }

  CLHEP::Hep3Vector DiskCalorimeter::diskToMu2e  (unsigned diskId, const CLHEP::Hep3Vector& p) const {
    return disk(diskId).diskToMu2e(p);
  }

  CLHEP::Hep3Vector DiskCalorimeter::mu2eToDiskFF(unsigned diskId, const CLHEP::Hep3Vector& p) const {
    return disk(diskId).mu2eToDiskFF(p);
  }

  CLHEP::Hep3Vector DiskCalorimeter::diskFFToMu2e(unsigned diskId, const CLHEP::Hep3Vector& p) const {
    return disk(diskId).diskFFToMu2e(p);
  }

  CLHEP::Hep3Vector DiskCalorimeter::mu2eToCrystal(unsigned crystalId, const CLHEP::Hep3Vector& p) const {
      const Crystal& c = crystal(crystalId);            // global crystal cache
      return disk(c.diskID()).mu2eToCrystal(c.localID(), p);
  }

  CLHEP::Hep3Vector DiskCalorimeter::crystalToMu2e(unsigned crystalId, const CLHEP::Hep3Vector& p) const {
      const Crystal& c = crystal(crystalId);
      return disk(c.diskID()).crystalToMu2e(c.localID(), p);
  }

  // MOVE TO TRACKER CODE
  //-----------------------------------------------------------------------------
  CLHEP::Hep3Vector DiskCalorimeter::mu2eToTracker(const CLHEP::Hep3Vector& pos) const
  {
    return pos - trackerCenter_;
  }

  //-----------------------------------------------------------------------------
  CLHEP::Hep3Vector DiskCalorimeter::trackerToMu2e(const CLHEP::Hep3Vector& pos) const
  {
    return pos + trackerCenter_;
  }



  void DiskCalorimeter::print(std::ostream &os) const
  {
     os<<"Disk calorimeter "<<std::endl;
     os<<"Number of disks :"<< disks_.size()<<std::endl;
     for (size_t idisk=0;idisk<disks_.size();++idisk) disk(idisk).print(os);
  }
}
