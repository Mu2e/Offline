#ifndef CalorimeterGeom_CaloGeomUtil_hh
#define CalorimeterGeom_CaloGeomUtil_hh
//
// Contains geometry utilities: coordinates transformations and checks if you are inside disks
// FF is the abrevation of front face
//
// Original author B. Echenard
//

#include "CalorimeterGeom/inc/CaloGeomUtil.hh"
#include "CalorimeterGeom/inc/CaloInfo.hh"
#include "CalorimeterGeom/inc/CaloGeomInfo.hh"
#include "CalorimeterGeom/inc/Disk.hh"

#include "CLHEP/Vector/ThreeVector.h"

#include <vector>
#include <memory>

namespace mu2e {

    

    class CaloGeomUtil {

       public:
          
          CaloGeomUtil(const std::vector<std::shared_ptr<Disk>>& disks,
                       const CaloInfo& diskInfo, 
                       const CaloGeomInfo& geomInfo,
                       const std::vector<Crystal const*>& fullCrystalList);
          
          ~CaloGeomUtil() {};

          CLHEP::Hep3Vector mu2eToCrystal(  int crystalId, const CLHEP::Hep3Vector& pos) const;
	  CLHEP::Hep3Vector mu2eToDisk(     int diskId, const CLHEP::Hep3Vector& pos) const;	    
	  CLHEP::Hep3Vector mu2eToDiskFF(   int diskId, const CLHEP::Hep3Vector& pos) const;	    
	  CLHEP::Hep3Vector mu2eToTracker(  const CLHEP::Hep3Vector& pos)                const;

	  CLHEP::Hep3Vector crystalToMu2e(  int crystalId, const CLHEP::Hep3Vector& pos) const;
	  CLHEP::Hep3Vector diskToMu2e(     int diskId, const CLHEP::Hep3Vector& pos) const;
	  CLHEP::Hep3Vector diskFFToMu2e(   int diskId, const CLHEP::Hep3Vector& pos) const;
	  CLHEP::Hep3Vector trackerToMu2e(  const CLHEP::Hep3Vector& pos)                const;

	  bool isInsideCalorimeter(const CLHEP::Hep3Vector& pos)                                 const;       	 	 
          bool isInsideSection(int iDisk, const CLHEP::Hep3Vector& pos)                       const;
	  bool isContainedSection(const CLHEP::Hep3Vector& front, const CLHEP::Hep3Vector& back) const;
	  
          const Disk&  disk(int i) const  {return *disks_.at(i);}



       private:
          
          const std::vector<std::shared_ptr<Disk>>& disks_;
       	  const CaloInfo&                           caloInfo_;
       	  const CaloGeomInfo&                       geomInfo_;
          const std::vector<Crystal const*>&        fullCrystalList_;
     };

}    

#endif 
