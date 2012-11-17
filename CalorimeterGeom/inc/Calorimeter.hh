#ifndef CalorimeterGeom_Calorimeter_hh
#define CalorimeterGeom_Calorimeter_hh
//
// Hold all geometry and identifier information about
// a Calorimeter.  In order to insulate this class from
// knowledge of databases etc, this class must not know
// how to make itself.
//
// $Id: Calorimeter.hh,v 1.21 2012/11/17 00:06:25 echenard Exp $
// $Author: echenard $
// $Date: 2012/11/17 00:06:25 $
//
// Original author R. Bernstein and Rob Kutschke
//

// C++ includes
#include <vector>

// Mu2e includes
#include "Mu2eInterfaces/inc/Detector.hh"

//CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"


namespace mu2e {
    
    class Calorimeter: virtual public Detector{


	public:
	  //no constructor for this interface
	  virtual ~Calorimeter(){}


	  // coordinate position and transformation
	  virtual CLHEP::Hep3Vector const& getOrigin(void) const = 0; 
	  
	  virtual CLHEP::Hep3Vector toCrystalFrame(int crystalId, CLHEP::Hep3Vector const& pos) const = 0;
	  virtual CLHEP::Hep3Vector toSectionFrame(int sectionId, CLHEP::Hep3Vector const& pos) const = 0;
	  virtual CLHEP::Hep3Vector fromSectionFrame(int sectionId, CLHEP::Hep3Vector const& pos) const = 0;

          virtual CLHEP::Hep3Vector getCrystalOrigin(int crystalId) const =0;
	  virtual CLHEP::Hep3Vector getLocalCrystalOrigin(int crystalId) const = 0;
          virtual std::vector<int>  getNeighbors(int crystalId, int level=1) const = 0;


	  //crystal / readout section
	  virtual unsigned int nROPerCrystal(void) const = 0;
	  virtual unsigned int nCrystal(void) const  = 0;
	  virtual unsigned int nRO(void) const  = 0;
	  virtual double       crystalVolume(void) const = 0;

	  //crystal - readout id
	  virtual int getCrystalByRO(int roid) const = 0;
	  virtual int getROBaseByCrystal(int id) const = 0;
	  virtual int getCaloSectionId(int crystalId) const = 0;


//bunch of accessors that will go away in later, keep them now to compile the code

	  //crystal characteristics      
	  virtual double getNonuniformity() const = 0;
	  virtual double getTimeGap() const = 0;
	  virtual double getElectronEdep() const = 0;
	  virtual double getElectronEmin() const = 0;
	  virtual double roHalfSize() const = 0;
	  virtual double roHalfThickness() const = 0;

	  //wrapper characteristics
	  virtual double wrapperThickness() const = 0;

 
    };

}
#endif /* CalorimeterGeom_Calorimeter_hh */
