//
// Make a Calorimeter.
//
// $Id: CalorimeterMaker.cc,v 1.6 2010/04/29 18:23:18 rhbob Exp $
// $Author: rhbob $
// $Date: 2010/04/29 18:23:18 $

// original authors Julie Managan and Robert Bernstein

//
// C++ includes

//
// Mu2e includes
#include "CalorimeterGeom/inc/CalorimeterMaker.hh"
#include "Mu2eUtilities/inc/hep3VectorFromStdVector.hh"
#include "CalorimeterGeom/inc/CrystalDetail.hh"
#include "CalorimeterGeom/inc/Crystal.hh"
#include "CalorimeterGeom/inc/Vane.hh"
#include "CalorimeterGeom/inc/RSlice.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/CrystalId.hh"

// Framework include files
#include "FWCore/Utilities/interface/Exception.h"



//
// other includes
#include "CLHEP/Vector/RotationX.h"
#include "CLHEP/Vector/RotationY.h"
#include "CLHEP/Vector/RotationZ.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h"

using namespace std;
using namespace CLHEP;



//
// naming conventions:  
// a) a vane is one of the sets of crystals; 
// b) a "z-slice" is one set of crystals in a vane with long side aligned, at constant
//    z in an ideal array;
// c) an "r-slice" is one set of crystals in a vane at the same distance
//    from the z-axis


namespace mu2e{

  namespace calorimeter {
    CalorimeterMaker::CalorimeterMaker( SimpleConfig const& config)
    {
      numberOfVanes            = config.getInt   ("calorimeter.numberOfVanes");
      crystalHalfTrans         = config.getDouble("calorimeter.crystalHalfTrans");
      crystalHalfLong          = config.getDouble("calorimeter.crystalHalfLong");
      nCrystalRSlices          = config.getInt   ("calorimeter.nCrystalRSlices");   
      nCrystalZSlices          = config.getInt   ("calorimeter.nCrystalZSlices");
      calorimeterCenter        = config.getHep3Vector("calorimeter.calorimeterCenter");
      calorimeterCenterOffset  = config.getHep3Vector("calorimeter.calorimeterCenterOffset");
      rInscribed               = config.getDouble("calorimeter.rInscribed");
      phi0                     = config.getDouble("calorimeter.phi0");
      theta0                   = config.getDouble("calorimeter.theta0");

      //
      //rotations of vanes (from system with theta0 = 0, phi0 = 0) in Goldstein convention
      config.getVectorDouble("calorimeter.calorimeterVaneRotationsX",calorimeterVaneRotationsX,numberOfVanes);
      config.getVectorDouble("calorimeter.calorimeterVaneRotationsY",calorimeterVaneRotationsY,numberOfVanes);
      config.getVectorDouble("calorimeter.calorimeterVaneRotationsZ",calorimeterVaneRotationsZ,numberOfVanes);

      //
      //offsets of vanes(wrt calorimeterCenter vector)
      config.getVectorDouble("calorimeter.calorimeterVaneOffsetsX",calorimeterVaneOffsetsX,numberOfVanes);
      config.getVectorDouble("calorimeter.calorimeterVaneOffsetsY",calorimeterVaneOffsetsY,numberOfVanes);
      config.getVectorDouble("calorimeter.calorimeterVaneOffsetsZ",calorimeterVaneOffsetsZ,numberOfVanes);
      
      //
      // crystal materials
      crystalMaterial              = config.getString("calorimeter.crystalMaterial");
      crystalWrapper               = config.getString("calorimeter.crystalWrapper");
      crystalWrapperHalfThickness  = config.getDouble("calorimeter.crystalWrapperHalfThickness");

      //
      //make sure above information is consistent
      CheckIt();

      //
      //do the work
      BuildIt();
    }


    CalorimeterMaker::~CalorimeterMaker() {}



    bool CalorimeterMaker::CheckIt(){ return true;}
    bool CalorimeterMaker::CheckVaneConsistency(){ return true;}

    void CalorimeterMaker::BuildIt()
    {
      _calorimeter = auto_ptr<Calorimeter>(new Calorimeter());

      //
      // here's the model:  I make crystals without regard
      // for orientation. Right now all crystals are identical.

      // I take a vane and orient it as a whole.  Every crystal 
      // has a vector pointing down its long axis from the inside to the 
      // outside, and that vector gets rotated or translated as the vane is 
      // put into place


      MakeGenericCrystals();  // if you make more than one type, not generic
      MakeVanes();
      MakeCalorimeter();

      return;
    }

    void CalorimeterMaker::MakeGenericCrystals(){

      //
      // make crystals and assign the details to them.  
      // for now, all crystals have identical details.  If we
      // want more than one type, this is a natural place to do it


      // 
      // the idea here is to create crystals without regard for their geometric
      // location or orientation.  

      _calorimeter->_crystalDetail.push_back
	(CrystalDetail(
		       crystalHalfTrans,
		       crystalHalfLong,
		       crystalMaterial,
		       crystalWrapper,
      		       crystalWrapperHalfThickness));
      return;
    }

    void CalorimeterMaker::MakeVanes()
    {
      
      deque  <Crystal>& allCrystals = _calorimeter->_allCrystals;
      vector <Vane>&    allVanes    = _calorimeter->_vanes;
      
      //
      // loop over vanes

      //
      // set initial direction for vector along crystal long axis in this vane
      Hep3Vector initialLongAxis(1.,0.,0.); // vane at 6 o'clock if four vanes
      for (uint32_t ithVane = 0; ithVane < numberOfVanes ; ++ithVane)
	{
	  allVanes.push_back(Vane(ithVane));
	  Vane& currentVane = allVanes[ithVane]; 

	  //
	  //start making RSlices and ZSlices in this vane
	  vector<ZSlice>& zsliceCurrentVane = currentVane._zslices;

	  // 
	  //all zslices in a vane have the same phi. 
	  double phiZSlice = ithVane*CLHEP::twopi/numberOfVanes + phi0;
	  //
	  // rotate "ideal vane" into its phi position, then adjust for angular offset of vane
	  Hep3Vector currentLongAxis = HepRotation(calorimeterVaneRotationsX[ithVane],
						   calorimeterVaneRotationsY[ithVane],
						   calorimeterVaneRotationsZ[ithVane])
	    *HepRotationZ(phiZSlice)*initialLongAxis;
	  cout <<"wire of vane number " << ithVane << " is " << currentLongAxis << endl;
	  //	  assert(2==1);
	}
      return;}

    void CalorimeterMaker::FillNearestNeighbours(){return;}
    void CalorimeterMaker::FillPointersAndIndices(){return;}
   
    uint32_t  CalorimeterMaker::NumberOfCrystalsPerVane() 
    {return nCrystalRSlices*nCrystalZSlices;}
    uint32_t  CalorimeterMaker::TotalNumberOfCrystals()   
    {return nCrystalRSlices*nCrystalZSlices*numberOfVanes;}
    void CalorimeterMaker::MakeCalorimeter()
    { return;}
  } //namespace calorimeter
} //namespace mu2e

