#ifndef CalorimeterGeom_Calorimeter_hh
#define CalorimeterGeom_Calorimeter_hh
//
// Hold all geometry and identifier information about
// a Calorimeter.  In order to insulate this class from
// knowledge of databases etc, this class must not know
// how to make itself.
//
// $Id: Calorimeter.hh,v 1.26 2014/08/01 20:57:44 echenard Exp $
// $Author: echenard $
// $Date: 2014/08/01 20:57:44 $
//
// Original author R. Bernstein and Rob Kutschke
//

// C++ includes
#include <vector>

// Mu2e includes
#include "Mu2eInterfaces/inc/Detector.hh"
#include "BaseCalorimeterData.hh"

//CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"


namespace mu2e {
    
    class Calorimeter: virtual public Detector{


	public:
	  //no constructor for this interface
	  virtual ~Calorimeter(){}


	  //accessor to geometry data
	  virtual BaseCalorimeterData const& caloGeomInfo() const = 0;

	  // coordinate position and transformation - origin refers to the Mu2e frame
	  virtual CLHEP::Hep3Vector const& origin()                                                      const = 0; 	  
	  virtual CLHEP::Hep3Vector        toCrystalFrame(int crystalId, CLHEP::Hep3Vector const& pos)   const = 0;
	  virtual CLHEP::Hep3Vector        toSectionFrame(int sectionId, CLHEP::Hep3Vector const& pos)   const = 0;
	  virtual CLHEP::Hep3Vector        fromCrystalFrame(int crystalId, CLHEP::Hep3Vector const& pos) const = 0;
	  virtual CLHEP::Hep3Vector        fromSectionFrame(int sectionId, CLHEP::Hep3Vector const& pos) const = 0;
          virtual CLHEP::Hep3Vector        crystalOrigin(int crystalId)                                  const = 0;
	  virtual CLHEP::Hep3Vector        crystalOriginInSection(int crystalId)                         const = 0;
	  virtual bool                     isInsideCalorimeter(CLHEP::Hep3Vector const& pos)             const = 0  ;
	  virtual double                   crystalLongPos(int crystalId, CLHEP::Hep3Vector const& pos)   const = 0;
	  
          
	  // crystal neighbors, indexing 
          virtual std::vector<int>  const& neighbors(int crystalId)            const = 0;	  
          virtual std::vector<int>  const& nextNeighbors(int crystalId)        const = 0;	  
          virtual std::vector<int>  const& nextNextNeighbors(int crystalId)    const = 0;	  
          virtual std::vector<int>  neighborsByLevel(int crystalId, int level) const = 0;	  


	  // crystal / readout section
	  virtual unsigned int nCrystal()                   const = 0;
	  virtual unsigned int nRO()                        const = 0;
	  virtual          int crystalByRO(int roid)        const = 0;
	  virtual          int ROBaseByCrystal(int id)      const = 0;
 	  virtual          int caloSectionId(int crystalId) const = 0;




    };

}
#endif /* CalorimeterGeom_Calorimeter_hh */
