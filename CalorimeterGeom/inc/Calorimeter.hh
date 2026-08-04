#ifndef CalorimeterGeom_Calorimeter_hh
#define CalorimeterGeom_Calorimeter_hh
//
// Interface to the disk calorimeter
// Original author B. Echenard
//
#include "Offline/Mu2eInterfaces/inc/Detector.hh"
#include "Offline/CalorimeterGeom/inc/CaloG4Info.hh"
#include "Offline/CalorimeterGeom/inc/Disk.hh"
#include "Offline/CalorimeterGeom/inc/Crystal.hh"

#include "CLHEP/Vector/ThreeVector.h"
#include <vector>
#include <iostream>

namespace mu2e {

  class Calorimeter: virtual public Detector  {

    public:
      using Disks      = std::vector<Disk>;
      using Crystals   = std::vector<const Crystal*>;
      using Hep3Vector = CLHEP::Hep3Vector;

      constexpr static const char* cxname = {"Calorimeter"};

      //no constructor for this interface
      virtual ~Calorimeter() = default;

      virtual size_t            nDisks()            const = 0;
      virtual const Disk&       disk(unsigned i)    const = 0;
      virtual const Disks&      disks()             const = 0;

      virtual size_t            nCrystals()         const = 0;
      virtual const Crystal&    crystal(unsigned i) const = 0;
      virtual const Crystals&   crystals()          const = 0;

      virtual const CaloG4Info& G4Info()            const = 0;
      virtual void              print(std::ostream &os = std::cout) const = 0;

      virtual std::vector<unsigned> neighbors(unsigned globalId, unsigned level) const = 0;
      virtual bool              isInsideAnyCrystal(const Hep3Vector& pos)    const = 0;
      virtual bool              isInsideAnyDisk   (const Hep3Vector& pos)    const = 0;
      virtual bool              isInsideSameDisk  (const Hep3Vector& front,
                                                   const Hep3Vector& back)   const = 0;

      virtual Hep3Vector        mu2eToDisk   (unsigned diskId,    const Hep3Vector& p) const = 0;
      virtual Hep3Vector        diskToMu2e   (unsigned diskId,    const Hep3Vector& p) const = 0;
      virtual Hep3Vector        mu2eToDiskFF (unsigned diskId,    const Hep3Vector& p) const = 0;
      virtual Hep3Vector        diskFFToMu2e (unsigned diskId,    const Hep3Vector& p) const = 0;
      virtual Hep3Vector        mu2eToCrystal(unsigned crystalId, const Hep3Vector& p) const = 0;
      virtual Hep3Vector        crystalToMu2e(unsigned crystalId, const Hep3Vector& p) const = 0;

      // These should REALLY go to the tracker
      virtual Hep3Vector        mu2eToTracker(const Hep3Vector& p) const = 0;
      virtual Hep3Vector        trackerToMu2e(const Hep3Vector& p) const = 0;

  };
}
#endif
