#ifndef CalorimeterGeom_Calorimeter_hh
#define CalorimeterGeom_Calorimeter_hh
//
// Hold all geometry and identifier information about
// a Calorimeter.  In order to insulate this class from
// knowledge of databases etc, this class must not know
// how to make itself.
//
// $Id: Calorimeter.hh,v 1.19 2012/09/08 02:24:25 echenard Exp $
// $Author: echenard $
// $Date: 2012/09/08 02:24:25 $
//
// Original author R. Bernstein and Rob Kutschke
//



// Mu2e includes
#include "Mu2eInterfaces/inc/Detector.hh"

//CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"



namespace mu2e {
    
    class Calorimeter: virtual public Detector{


	public:
	  //no constructor for this interface
	  ~Calorimeter(){}


	  //crystal characteristics      
	  virtual double getNonuniformity() const = 0;
	  virtual double getTimeGap() const = 0;
	  virtual double getElectronEdep() const = 0;
	  virtual double getElectronEmin() const = 0;

	  //readout section
	  virtual unsigned int nRO() const  = 0;
	  virtual unsigned int nROPerCrystal() const = 0;
	  virtual double roHalfSize() const = 0;
	  virtual double roHalfThickness() const = 0;

	  //wrapper characteristics
	  virtual double wrapperThickness() const = 0;

	  //crystal - readout id
	  virtual int getCrystalByRO(int roid) const = 0;
	  virtual int getROBaseByCrystal(int id) const = 0;

	  // coordinate position and transformation
	  virtual CLHEP::Hep3Vector const& getOrigin() const = 0;
	  virtual CLHEP::Hep3Vector toCrystalFrame(int, CLHEP::Hep3Vector const&) const = 0;

 
    };

} //namespace mu2e

#endif /* CalorimeterGeom_Calorimeter_hh */
