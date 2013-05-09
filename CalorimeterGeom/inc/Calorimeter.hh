#ifndef CalorimeterGeom_Calorimeter_hh
#define CalorimeterGeom_Calorimeter_hh
//
// Hold all geometry and identifier information about
// a Calorimeter.  In order to insulate this class from
// knowledge of databases etc, this class must not know
// how to make itself.
//
// $Id: Calorimeter.hh,v 1.24 2013/05/09 23:14:14 echenard Exp $
// $Author: echenard $
// $Date: 2013/05/09 23:14:14 $
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
	  virtual CLHEP::Hep3Vector const& origin() const = 0; 	  
	  virtual CLHEP::Hep3Vector        toCrystalFrame(int crystalId, CLHEP::Hep3Vector const& pos) const = 0;
	  virtual CLHEP::Hep3Vector        toSectionFrame(int sectionId, CLHEP::Hep3Vector const& pos) const = 0;
	  virtual CLHEP::Hep3Vector        fromSectionFrame(int sectionId, CLHEP::Hep3Vector const& pos) const = 0;
          virtual CLHEP::Hep3Vector        crystalOrigin(int crystalId) const =0;
	  virtual CLHEP::Hep3Vector        localCrystalOrigin(int crystalId) const = 0;
	  
          virtual std::vector<int>  neighbors(int crystalId, int level=1) const = 0;
	  
	  virtual bool isInsideCalorimeter(CLHEP::Hep3Vector const& pos) const =0  ;
	  virtual int  crystalIdxFromPosition(CLHEP::Hep3Vector const& pos) const =0;


	  //crystal / readout section
	  virtual unsigned int nROPerCrystal() const = 0;
	  virtual unsigned int nCrystal() const  = 0;
	  virtual unsigned int nRO() const  = 0;


	  //crystal - readout id
	  virtual int caloSectionId(int crystalId) const = 0;
	  virtual int crystalByRO(int roid) const = 0;
	  virtual int ROBaseByCrystal(int id) const = 0;

	  //crystal characteristics      
	  virtual double crystalVolume() const = 0;


          //a few accessors we keep for convenience
	  virtual double crystalHalfTrans() const = 0;
          virtual double crystalHalfLength() const = 0;
	  virtual double roHalfSize() const = 0;
	  virtual double roHalfThickness()   const = 0;
	  virtual double wrapperThickness() const = 0;
	  virtual double caseThickness() const = 0;

	  //bunch of accessors that will go away in later, keep them now to compile the code
	  virtual double getNonuniformity() const = 0;
	  virtual double getTimeGap() const = 0;
	  virtual double getElectronEdep() const = 0;
	  virtual double getElectronEmin() const = 0;
          virtual double apdMeanNoise()   const = 0;
	  virtual double apdSigmaNoise()  const = 0;
	  virtual double lysoLightYield() const = 0;
	  virtual double apdQuantumEff()  const = 0;
	  virtual double apdCollectEff()  const = 0;


    };

}
#endif /* CalorimeterGeom_Calorimeter_hh */
