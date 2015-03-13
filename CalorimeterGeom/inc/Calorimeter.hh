#ifndef CalorimeterGeom_Calorimeter_hh
#define CalorimeterGeom_Calorimeter_hh
//
// $Id: Calorimeter.hh,v 1.27 2014/08/01 21:49:38 echenard Exp $
// $Author: echenard $
// $Date: 2014/08/01 21:49:38 $
//
// Original author B. Echenard
//

// C++ includes
#include <vector>

// Mu2e includes
#include "Mu2eInterfaces/inc/Detector.hh"
#include "CalorimeterGeom/inc/BaseCalorimeterInfoGeom.hh"
#include "CalorimeterGeom/inc/BaseCalorimeterInfoGeom.hh"
#include "CalorimeterGeom/inc/Crystal.hh"
#include "CalorimeterGeom/inc/CaloSection.hh"

//CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"


namespace mu2e {
    
    class Calorimeter: virtual public Detector{


	public:

	  enum CaloType {disk,vane};


	  //no constructor for this interface
	  virtual ~Calorimeter(){}

	  //accessor to geometry data
	  virtual CaloType const& caloType() const = 0;
	  
	  //accessor to geometry data
	  virtual BaseCalorimeterInfoGeom const& caloGeomInfo() const = 0;

	  
	  // coordinate position and transformation - origin refers to the Mu2e frame
	  virtual CLHEP::Hep3Vector const& origin()                                                           const = 0; 	  
	  
	  virtual CLHEP::Hep3Vector        toCrystalFrame(int crystalId, CLHEP::Hep3Vector const& pos)        const = 0;
	  virtual CLHEP::Hep3Vector        toSectionFrame(int sectionId, CLHEP::Hep3Vector const& pos)        const = 0;
	  virtual CLHEP::Hep3Vector        toSectionFrameFF(int sectionId, CLHEP::Hep3Vector const& pos)      const = 0;
	  virtual CLHEP::Hep3Vector        toTrackerFrame(CLHEP::Hep3Vector const& pos)                       const = 0;
	  
	  virtual CLHEP::Hep3Vector        fromCrystalFrame(int crystalId, CLHEP::Hep3Vector const& pos)      const = 0;
	  virtual CLHEP::Hep3Vector        fromSectionFrame(int sectionId, CLHEP::Hep3Vector const& pos)      const = 0;
	  virtual CLHEP::Hep3Vector        fromSectionFrameFF(int sectionId, CLHEP::Hep3Vector const& pos)    const = 0;
	  virtual CLHEP::Hep3Vector        fromTrackerFrame(CLHEP::Hep3Vector const& pos)                     const = 0;
          
	  virtual CLHEP::Hep3Vector        crystalOrigin(int crystalId)                                       const = 0;	  
	  virtual CLHEP::Hep3Vector        crystalOriginInSection(int crystalId)                              const = 0;
	  virtual int                      crystalIdxFromPosition(CLHEP::Hep3Vector const& pos)               const = 0;

	  virtual bool                     isInsideCalorimeter(CLHEP::Hep3Vector const& pos)                  const = 0;
	  virtual bool                     isInsideSection(int iSection, CLHEP::Hep3Vector const& pos)        const = 0;
	  virtual void                     print()                                                            const = 0;

	  

	  //calo sections
	  virtual unsigned int          nSection()                 const = 0;
	  virtual CaloSection const&    section(int i)             const = 0;


	  //crystal / readout section
	  virtual unsigned int          nRO()                      const = 0;
	  virtual unsigned int          nCrystal()                 const = 0;
          virtual Crystal  const&       crystal(int i)             const = 0;

	  virtual          int          crystalByRO(int roid)      const = 0;
	  virtual          int          ROBaseByCrystal(int id)    const = 0;


	  // neighbors, indexing 
          virtual std::vector<int>  const& neighbors(int crystalId)            const = 0;	  
          virtual std::vector<int>  const& nextNeighbors(int crystalId)        const = 0;	  
          virtual std::vector<int>  neighborsByLevel(int crystalId, int level) const = 0;	  


    };

}
#endif /* CalorimeterGeom_Calorimeter_hh */
