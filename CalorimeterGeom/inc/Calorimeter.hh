#ifndef CalorimeterGeom_Calorimeter_hh
#define CalorimeterGeom_Calorimeter_hh
//
// Hold all geometry and identifier information about
// a Calorimeter.  In order to insulate this class from
// knowledge of databases etc, this class must not know
// how to make itself.
//
// $Id: Calorimeter.hh,v 1.22 2013/03/05 20:33:25 aluca Exp $
// $Author: aluca $
// $Date: 2013/03/05 20:33:25 $
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
	  virtual CLHEP::Hep3Vector const& origin(void) const = 0; 
	  
	  virtual CLHEP::Hep3Vector toCrystalFrame(int crystalId, CLHEP::Hep3Vector const& pos) const = 0;
	  virtual CLHEP::Hep3Vector toSectionFrame(int sectionId, CLHEP::Hep3Vector const& pos) const = 0;
	  virtual CLHEP::Hep3Vector fromSectionFrame(int sectionId, CLHEP::Hep3Vector const& pos) const = 0;

      virtual CLHEP::Hep3Vector crystalOrigin(int crystalId) const =0;
	  virtual CLHEP::Hep3Vector localCrystalOrigin(int crystalId) const = 0;
      virtual std::vector<int>  neighbors(int crystalId, int level=1) const = 0;


	  //crystal / readout section
	  virtual unsigned int nROPerCrystal(void) const = 0;
	  virtual unsigned int nCrystal(void) const  = 0;
	  virtual unsigned int nRO(void) const  = 0;
	  virtual double       crystalVolume(void) const = 0;
      virtual double       crystalHalfLength(void) const = 0;

	  //crystal - readout id
	  virtual int crystalByRO(int roid) const = 0;
	  virtual int ROBaseByCrystal(int id) const = 0;
	  virtual int caloSectionId(int crystalId) const = 0;

	  //calorimeter envelope
	  virtual double       envelopeRmin(void) const = 0;
	  virtual double       envelopeRmax(void) const = 0;
	  virtual double       envelopeHalfLength(void) const = 0;

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
