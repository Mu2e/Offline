//
// Make a Calorimeter.
//
// $Id: CalorimeterMaker.cc,v 1.14 2010/05/25 17:36:21 rhbob Exp $
// $Author: rhbob $
// $Date: 2010/05/25 17:36:21 $

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
#include "CLHEP/Matrix/DiagMatrix.h"
using namespace std;




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
      config.getVectorDouble("calorimeter.calorimeterVaneRotationsPhi",calorimeterVaneRotationsPhi,numberOfVanes);
      config.getVectorDouble("calorimeter.calorimeterVaneRotationsTheta",calorimeterVaneRotationsTheta,numberOfVanes);
      config.getVectorDouble("calorimeter.calorimeterVaneRotationsPsi",calorimeterVaneRotationsPsi,numberOfVanes);

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



    bool CalorimeterMaker::CheckIt(){

      //
      // just make sure everything is specified consistently
      return true;}
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
      FillNearestNeighbours();
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
      
      vector  <Crystal>& allCrystals = _calorimeter->_allCrystals;
      vector <Vane>&    allVanes    = _calorimeter->_vanes;

      allVanes.reserve(numberOfVanes);
      //
      // set initial direction for vector along crystal long axis
      // vane at 6 o'clock if four vanes
      CLHEP::Hep3Vector initialLongAxis(1.,0.,0.); 
      //
      // other useful numbers
      double angleOfRotation = CLHEP::twopi/numberOfVanes;
      double distanceBetweenCrystals = 2.*crystalHalfTrans;
      //
      // loop over vanes
      for (uint32_t ithVane = 0; ithVane < numberOfVanes ; ++ithVane)
        {
          allVanes.push_back(Vane(ithVane));
          Vane& currentVane = allVanes[ithVane]; 
          VaneId vid = ithVane;

          //
          //start making RSlices and ZSlices in this vane
          vector<ZSlice>& zSlicesCurrentVane = currentVane._zslices;

          // 
          //all zslices in a vane have the same phi. 
          double phiZSlice = ithVane*angleOfRotation + phi0;
          //
          // rotate "ideal vane" into its phi position, then adjust for angular offset of vane
          //
          // this assumes eulerian system
          CLHEP::HepRotation masterVaneRotation = 
            CLHEP::HepRotation(calorimeterVaneRotationsPhi[ithVane],
			       calorimeterVaneRotationsTheta[ithVane],
			       calorimeterVaneRotationsPsi[ithVane])  *  CLHEP::HepRotationZ(phiZSlice);
          //cout << "rotation matrix " << "\n" << masterVaneRotation << endl;
          CLHEP::Hep3Vector currentLongAxis = masterVaneRotation*initialLongAxis;
          //cout <<"wire of vane number " << ithVane << " is " << currentLongAxis << endl;

          //
          // create ZSlices
          zSlicesCurrentVane.reserve(nCrystalZSlices);
          for (uint32_t ithZSlice = 0; ithZSlice < nCrystalZSlices; ++ithZSlice)
            {
              zSlicesCurrentVane.push_back(ZSlice(ZSliceId(ithVane,ithZSlice)));
              //
              // get the address of this last CLHEP::piece so I can put things directly into the vector;
              // otherwise I need to fill the vector and then do a copy later
              ZSlice& currentZSlice = zSlicesCurrentVane.back();


              //
              // populate r slices at a fixed z in this vane
              double z0 = crystalHalfTrans + distanceBetweenCrystals*ithZSlice;
              for (uint32_t ithRSlice = 0; ithRSlice < nCrystalRSlices; ++ithRSlice)
                {
                  //
                  // calculate various useful things in format that which will be convenient for
                  // G4 instead of human beings; names, etc are for vanes parallel to z but we still
                  // tweak them with offsets

                  //
                  // center of crystal in this rslice again with first one at 6PM for 4-vane
                  double x0 = 0.;
                  double y0 = -(crystalHalfTrans + distanceBetweenCrystals*ithRSlice + rInscribed); 
                  CLHEP::Hep3Vector origin(x0,y0,0.);
                  CLHEP::Hep3Vector delta(0.,0.,distanceBetweenCrystals);
                  RSliceId rid = RSliceId(vid,ithZSlice,ithRSlice);
                  RSlice currentRSlice = RSlice(rid,
						nCrystalRSlices,origin,delta);
                  currentZSlice._rslices.push_back(currentRSlice);
                  //
                  // so now for this vane, we've made a ZSlice populated by RSlices.
                  // Next we need to relate them to crystals.
                  //
                  // the below seems (and probably is) stupid.  But for each ZSlice there are 
                  // nCrystalRSlices, but each of them has exactly one crystal.  This structure
                  // was set up for the Trackers, which are geometrically more complicated
                 
                  vector<CrystalIndex> indices = currentRSlice._indices;
                  indices.reserve(1);
                  vector<const Crystal*>& crystals = currentRSlice._crystals;
                  crystals.reserve(1);
                  //
                  // and drop the crystal where it's supposed to be
                  CLHEP::Hep3Vector finalPositionWithinZSlice = masterVaneRotation*origin;
                  //cout << "origin = " << origin << endl;

                  //
                  CLHEP::Hep3Vector finalLongAxis = masterVaneRotation*initialLongAxis;

                  CrystalIndex index = CrystalIndex(allCrystals.size());
                  //
                  // if we ever have more than one sort of crystal
                  CrystalDetail* standardCrystal = &_calorimeter->_crystalDetail[0]; 

                  allCrystals.push_back(Crystal(CrystalId(rid),
                                                index,finalPositionWithinZSlice,standardCrystal,finalLongAxis));
                  indices.push_back(index);

                  /*
		    cout << "dumpola: " << "\n" <<
                    "size of crystal array" << allCrystals.size() << "\n" <<
                    "rid " << rid << " " << "\n"
                    "crystalId "<< CrystalId(rid,1) << " \n" <<
                    "index " << index << "\n" <<
                    "finalPositionWithinZSlice " << finalPositionWithinZSlice << "\n" <<
                    "finalLongAxis" << finalLongAxis << endl;
		  */
                }// close rslice
            } //close zslice
        }//close vane
      return;}




    void CalorimeterMaker::FillNearestNeighbours(){
      int32_t ifoo = -1;
      // Build the nearest neighbour info for each crystal.
      for ( vector<Crystal>::iterator i=_calorimeter->_allCrystals.begin(), e=_calorimeter->_allCrystals.end();
	    i != e; ++i){
	++ifoo;
	//cout << "crystal number" << ifoo << endl;
	//
	// for readability
	Crystal& crystal = *i;

	//
	// pull out information about this crystal's Vane, RSlice and ZSlice
	const VaneId&   crystalVaneId      = crystal.Id().getVaneId();
	const ZSliceId& crystalZSliceId    = crystal.Id().getZSliceId();
	const RSliceId& crystalRSliceId    = crystal.Id().getRSliceId();
	const ZSlice&   crystalZSlice      = _calorimeter->getZSlice(crystalZSliceId);
	const RSlice&   crystalRSlice      = _calorimeter->getRSlice(crystalRSliceId);

	int jv0 = crystal.Id().getVane();
	int jz0 = crystal.Id().getZSlice();
	int jr0 = crystal.Id().getRSlice();

	int jC = ifoo; 

	//cout << "for this crystal the vane, zslice, rslice and crystal number are " << jv0 << " " << jz0 << " " << jr0 << " " << jC << endl;

	//
	//innermost and outermost RSlices:

	//
	//innermost rslice:
	// first zslice has nothing upstream
	if (jr0 == 0){
	  if (jz0 == 0){
	    crystal._nearestById.push_back(CrystalId(crystalVaneId,jz0+1,jr0));
	    crystal._nearestById.push_back(CrystalId(crystalVaneId,jz0,jr0+1));
	    crystal._nearestById.push_back(CrystalId(crystalVaneId,jz0+1,jr0+1));
	  } 
	  // last zslice has nothing downstream
	  else if (jz0 == nCrystalZSlices-1){
	    crystal._nearestById.push_back(CrystalId(crystalVaneId,jz0-1,jr0));
	    crystal._nearestById.push_back(CrystalId(crystalVaneId,jz0,jr0+1));
	    crystal._nearestById.push_back(CrystalId(crystalVaneId,jz0-1,jr0+1));
	  }
	  else{
	    //
	    //all other zslices in innermost rslice
	    crystal._nearestById.push_back(CrystalId(crystalVaneId,jz0-1,jr0));
	    crystal._nearestById.push_back(CrystalId(crystalVaneId,jz0+1,jr0));
	    crystal._nearestById.push_back(CrystalId(crystalVaneId,jz0-1,jr0+1));
	    crystal._nearestById.push_back(CrystalId(crystalVaneId,jz0  ,jr0+1));
	    crystal._nearestById.push_back(CrystalId(crystalVaneId,jz0+1,jr0+1));

	  }
	}


	//
	//outermost rslice
	// first zslice has nothing upstream
	if (jr0 == nCrystalRSlices-1){
	  if (jz0 == 0){
	    crystal._nearestById.push_back(CrystalId(crystalVaneId,jz0+1,jr0));
	    crystal._nearestById.push_back(CrystalId(crystalVaneId,jz0  ,jr0-1));
	    crystal._nearestById.push_back(CrystalId(crystalVaneId,jz0+1,jr0-1));
	  } 
	  // last zslice has nothing downstream
	  else if (jz0 == nCrystalZSlices-1){
	    crystal._nearestById.push_back(CrystalId(crystalVaneId,jz0-1,jr0));
	    crystal._nearestById.push_back(CrystalId(crystalVaneId,jz0-1,jr0-1));
	    crystal._nearestById.push_back(CrystalId(crystalVaneId,jz0  ,jr0-1));
	  }
	  else{
	    //
	    //all other zslices in outermost rslice
	    crystal._nearestById.push_back(CrystalId(crystalVaneId,jz0-1,jr0));
	    crystal._nearestById.push_back(CrystalId(crystalVaneId,jz0+1,jr0));
	    crystal._nearestById.push_back(CrystalId(crystalVaneId,jz0-1,jr0+1));
	    crystal._nearestById.push_back(CrystalId(crystalVaneId,jz0  ,jr0+1));
	    crystal._nearestById.push_back(CrystalId(crystalVaneId,jz0+1,jr0+1));

	  }
	}



	//
	//innermost and outermost ZSlices:

	//
	//innermost zslice; corners done
	if (jz0 == 0 && (jr0 !=0 && jr0 != nCrystalRSlices -1)) {
	  crystal._nearestById.push_back(CrystalId(crystalVaneId,jz0   ,jr0+1));
	  crystal._nearestById.push_back(CrystalId(crystalVaneId,jz0   ,jr0-1));
	  crystal._nearestById.push_back(CrystalId(crystalVaneId,jz0+1 ,jr0-1));
	  crystal._nearestById.push_back(CrystalId(crystalVaneId,jz0+1 ,jr0  ));
	  crystal._nearestById.push_back(CrystalId(crystalVaneId,jz0+1 ,jr0+1));
	}


	//
	//outermost zslice
	if (jz0 == nCrystalZSlices-1 && (jr0!=0 && jr0 != nCrystalZSlices-1)){
	  crystal._nearestById.push_back(CrystalId(crystalVaneId,jz0   ,jr0+1));
	  crystal._nearestById.push_back(CrystalId(crystalVaneId,jz0   ,jr0-1));
	  crystal._nearestById.push_back(CrystalId(crystalVaneId,jz0-1 ,jr0-1));
	  crystal._nearestById.push_back(CrystalId(crystalVaneId,jz0-1 ,jr0  ));
	  crystal._nearestById.push_back(CrystalId(crystalVaneId,jz0-1 ,jr0+1));
	} 








	//
	//everyone else ( yes, I know I could have written these special cases as a giant if else clause but it only gets executed once/job
	//and it's just easier for me to follow.  and I could have unrolled this loop too.  )
	if (jr0 != 0 && jr0 != nCrystalRSlices-1 && jz0 != 0 && jz0 != nCrystalZSlices-1){
	  for (uint32_t iz0 = jz0-1; iz0 <= jz0 + 1; ++iz0){
	    for (uint32_t ir0 = jr0 -1; ir0 <= jr0 + 1; ++ir0){
	      //   cout << "inside standard loop" << endl;
	      if (!(iz0 == jz0 && ir0 == jr0)){crystal._nearestById.push_back(CrystalId(crystalVaneId,iz0,ir0));}
	    }
	  }
	}

	// print out results
	for (uint32_t idx = 0; idx < crystal._nearestById.size(); ++idx)
	  {
	    // cout << "nearest neighbor " << idx+1 << crystal._nearestById[idx] << endl;
	  }
      }
    

      //      assert(2==1);
      return;
    }



    void CalorimeterMaker::FillPointersAndIndices(){

      // Fill the main link
      for ( vector<Vane>::iterator iVane = _calorimeter->_vanes.begin(),
	      eVane = _calorimeter->_vanes.end(); 
	    iVane != eVane;  ++iVane ){
	for ( vector<ZSlice>::iterator iZSlice = iVane->_zslices.begin(),
		eZSlice = iVane->_zslices.end();
	      iZSlice != eZSlice; ++iZSlice ){
	  for ( vector<RSlice>::iterator iRSlice = iZSlice->_rslices.begin(),
		  eRSlice = iZSlice->_rslices.end();
		iRSlice != eRSlice; ++iRSlice ){
	    iRSlice->_crystals.clear();

	    for ( vector<CrystalIndex>::iterator iCrys = iRSlice->_indices.begin(),
		    eCrys = iRSlice->_indices.end();
		  iCrys != eCrys ; ++iCrys ){
	      const Crystal& cryst = _calorimeter->_allCrystals[(*iCrys).asInt()];
	      iRSlice->_crystals.push_back( &cryst );
	    }
	  }
	}
      }
      return;
    }
   

  void CalorimeterMaker::FillPointersAndIndicesByNN(){

    // Fill nearest neighbour indices and pointers from the nearest neighbor Ids.
    for ( vector<Crystal>::iterator iCrys=_calorimeter->_allCrystals.begin(), 
            eCrys=_calorimeter->_allCrystals.end();
          iCrys != eCrys; 
          ++iCrys){
      vector<CrystalId>& byId = iCrys->_nearestById;
      vector<CrystalIndex>& byIndex = iCrys->_nearestByIndex;

      byIndex.clear();

      for ( vector<CrystalId>::iterator jCrys=byId.begin(), jCrysEnd=byId.end();
            jCrys != jCrysEnd; ++jCrys){
        const CrystalId& id = *jCrys;
        const Crystal& cryst = _calorimeter->getCrystal(id);
        byIndex.push_back( cryst.Index() );
      }
    }
    return;
  }

  uint32_t  CalorimeterMaker::NumberOfCrystalsPerVane() {
    return nCrystalRSlices*nCrystalZSlices;
  }

  uint32_t  CalorimeterMaker::TotalNumberOfCrystals(){
    return nCrystalRSlices*nCrystalZSlices*numberOfVanes;
  }
  void CalorimeterMaker::MakeCalorimeter() {
    return;
  }
} //namespace calorimeter
} //namespace mu2e

