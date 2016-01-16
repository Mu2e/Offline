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

	   enum CaloType {none,disk,vane};


	   //no constructor for this interface
	   virtual ~Calorimeter(){}


	   virtual const CaloType& caloType()                    const = 0;
	   virtual void  print(std::ostream &os = std::cout)     const = 0;


	   //accessor to geometry data
	   virtual const BaseCalorimeterInfoGeom& caloGeomInfo() const = 0;


	   // Coordinate position and transformation  - The REFERENCE is _Mu2e_ frame
	   // "from" methods: from XXX to Mu2e 
           // "to" methods  : from Mu2e  to XXX     
	   virtual const CLHEP::Hep3Vector& origin() const = 0; 	  

	   virtual CLHEP::Hep3Vector        toCrystalFrame(int crystalId, const CLHEP::Hep3Vector &pos)            const = 0;
	   virtual CLHEP::Hep3Vector        toSectionFrame(int sectionId, const CLHEP::Hep3Vector &pos)            const = 0;
	   virtual CLHEP::Hep3Vector        toSectionFrameFF(int sectionId, const CLHEP::Hep3Vector &pos)          const = 0;
	   virtual CLHEP::Hep3Vector        toTrackerFrame(const CLHEP::Hep3Vector &pos)                           const = 0;

	   virtual CLHEP::Hep3Vector        fromCrystalFrame(int crystalId, const CLHEP::Hep3Vector &pos)          const = 0;
	   virtual CLHEP::Hep3Vector        fromSectionFrame(int sectionId, const CLHEP::Hep3Vector &pos)          const = 0;
	   virtual CLHEP::Hep3Vector        fromSectionFrameFF(int sectionId, const CLHEP::Hep3Vector &pos)        const = 0;
	   virtual CLHEP::Hep3Vector        fromTrackerFrame(const CLHEP::Hep3Vector &pos)                         const = 0;



	   //position checks
	   virtual int                      crystalIdxFromPosition(const CLHEP::Hep3Vector &pos)                   const = 0;
	   virtual bool                     isInsideCalorimeter(const CLHEP::Hep3Vector &pos)                      const = 0;
	   virtual bool                     isInsideSection(int iSection, const CLHEP::Hep3Vector &pos)            const = 0;
	   virtual bool                     isContainedSection(CLHEP::Hep3Vector const&, CLHEP::Hep3Vector const&) const = 0;



	   //access to calo sections
	   virtual unsigned int             nSection()                                 const = 0;
	   virtual const CaloSection&       section(int i)                             const = 0;


	   //access to crystal / readout section
	   virtual unsigned int             nCrystal()                                 const = 0;
           virtual const Crystal&           crystal(int i)                             const = 0;

	   virtual unsigned int             nRO()                                      const = 0;
	   virtual          int             crystalByRO(int roid)                      const = 0;
	   virtual          int             ROBaseByCrystal(int id)                    const = 0;

	   // access to neighbors, indexing 
           virtual const std::vector<int>&  neighbors(int crystalId)                   const = 0;	  
           virtual const std::vector<int>&  nextNeighbors(int crystalId)               const = 0;	  
           virtual std::vector<int>         neighborsByLevel(int crystalId, int level) const = 0; //use only for level >2	  


    };

}

#endif 
