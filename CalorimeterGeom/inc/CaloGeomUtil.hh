//
// Contains geometry utilities: coordinates transformations and checks if you are inside disks
// FF is the abrevation of front face
//
// Original author B. Echenard
//

#ifndef CalorimeterGeom_CaloGeomUtil_hh
#define CalorimeterGeom_CaloGeomUtil_hh

#include "Offline/CalorimeterGeom/inc/CaloGeomUtil.hh"
#include "Offline/CalorimeterGeom/inc/Disk.hh"

#include "CLHEP/Vector/ThreeVector.h"
#include "cetlib_except/exception.h"
#include <vector>
#include <memory>


namespace mu2e {


    class CaloGeomUtil {

       public:
          CaloGeomUtil(const std::vector<std::shared_ptr<Disk>>& disks,
                       const std::vector<const Crystal*>& fullCrystalList);

          void origin        (const CLHEP::Hep3Vector& vec) {origin_ = vec;}
          void trackerCenter (const CLHEP::Hep3Vector& vec) {trackerCenter_ = vec;}
          void crystalZLength(double value)                 {crystalZLength_ = value;}

          const CLHEP::Hep3Vector& origin()          const {return origin_;}
          const CLHEP::Hep3Vector& trackerCenter()   const {return trackerCenter_;}
          const double             crystalZLength()  const {return crystalZLength_;}
          const Disk&              disk(int i)       const {return *disks_.at(i);}


          CLHEP::Hep3Vector mu2eToCrystal(int crystalId, const CLHEP::Hep3Vector& pos) const;
          CLHEP::Hep3Vector mu2eToDisk   (int diskId,    const CLHEP::Hep3Vector& pos) const;
          CLHEP::Hep3Vector mu2eToDiskFF (int diskId,    const CLHEP::Hep3Vector& pos) const;
          CLHEP::Hep3Vector mu2eToTracker(const CLHEP::Hep3Vector& pos)                const;

          CLHEP::Hep3Vector crystalToMu2e(int crystalId, const CLHEP::Hep3Vector& pos) const;
          CLHEP::Hep3Vector diskToMu2e   (int diskId,    const CLHEP::Hep3Vector& pos) const;
          CLHEP::Hep3Vector diskFFToMu2e (int diskId,    const CLHEP::Hep3Vector& pos) const;
          CLHEP::Hep3Vector trackerToMu2e(const CLHEP::Hep3Vector& pos)                const;

          bool isInsideCalorimeter       (const CLHEP::Hep3Vector& pos)                                  const;
          bool isInsideSection           (unsigned iDisk, const CLHEP::Hep3Vector& pos)                  const;
          bool isContainedSection        (const CLHEP::Hep3Vector& front, const CLHEP::Hep3Vector& back) const;


       private:
          const std::vector<std::shared_ptr<Disk>>& disks_;
          const std::vector<const Crystal*>&        fullCrystalList_; //non-owning pointers
          CLHEP::Hep3Vector                         origin_;
          CLHEP::Hep3Vector                         trackerCenter_;
          double                                    crystalZLength_;
     };
}

#endif
