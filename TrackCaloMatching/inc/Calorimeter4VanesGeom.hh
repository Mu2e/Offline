//
// $Id: Calorimeter4VanesGeom.hh,v 1.17 2014/08/01 20:57:45 echenard Exp $
// $Author: echenard $
// $Date: 2014/08/01 20:57:45 $
//
// Original author G. Pezzullo & G. Tassielli
//


#ifndef Calorimeter4VanesGeom_HH
#define Calorimeter4VanesGeom_HH


// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Principal/Handle.h"

//CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Geometry/HepPoint.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Matrix/SymMatrix.h"

#include "KalmanTests/inc/KalRepCollection.hh"

// From the art tool-chain
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//tracker includes
#include "TrkBase/TrkRep.hh"
#include "KalmanTrack/KalRep.hh"
#include "KalmanTests/inc/TrkFitDirection.hh"
// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TTrackerGeom/inc/TTracker.hh"

// BaBar
#include "BaBar/BaBar.hh"
#include "BaBar/include/DetectorModel/DetVolumeElem.hh"
#include "BaBar/include/DetectorModel/DetVolumeType.hh"
#include "BaBar/CLHEP/include/Geometry/Transformation.h"
#include "BaBar/include/DetectorModel/DetElem.hh"
#include "BaBar/include/TrkBase/TrkDifTraj.hh"
#include "BaBar/include/TrkBase/HelixTraj.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"
#include "KalmanTests/inc/KalRepCollection.hh"

//calorimeter include
#include "CalorimeterGeom/inc/VaneCalorimeter.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "TrackCaloMatching/inc/CaloVolumeElem.hh"
#include "TrackCaloMatching/inc/CaloVolumeType.hh"
#include "TrackCaloMatching/inc/CaloSurface.hh"
#include "CaloCluster/inc/CaloClusterer.hh"
//#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"

// Other includes.
#include "cetlib/exception.h"

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"

// From the art tool-chain
#include <cmath>
#include <deque>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <functional>

#include "CLHEP/Geometry/HepPoint.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Geometry/Transformation.h"

class HepPoint;
class HepTransformation;



namespace mu2e{

typedef std::vector<std::pair<double, double> > Length;


class Calorimeter4VanesGeom {
public :

    struct IntersectData_t {
      int    fVane;
      int    fRC;			// return code, 0=success, <0: failure, details TBD
      double fSEntr;
      double fSExit;
    };
    
        //Construct from a transform.
        Calorimeter4VanesGeom();

        //Destructor
        ~Calorimeter4VanesGeom ();

        CLHEP::Hep3Vector const& norm(){return _norm;}

        int nVanes() const{ return _nVanes;}

        //Half value of the radial dimension of the vane
        double dR() const{ return _dR;}

        //Half of the longitudinal dimension of the vane
        double dZ() const{ return _dZ;}

        //Half value of the crystal+wrapping+shell dimension
        double dU() const{ return _dU;}

        //Half value of the read-out device thickness
        double dAPD() const{ return _dAPD;}

        //Half value of the vane thickness
        double vaneHalfThickness() const{ return _vaneHalfThickness;}

        double solenoidOffSetX() const{ return _solenoidOffSetX;}

        double solenoidOffSetZ() const{ return _solenoidOffSetZ;}

        double innerRadius() const{ return _innerRadius;}

        double outherRadius() const{ return _outherRadius;}

        double ZfrontFaceCalo() const{ return _ZfrontFaceCalo;}

        double ZbackFaceCalo() const{ return _ZbackFaceCalo;}

        CaloVolumeElem* vane(int& i);
  
  CLHEP::Hep3Vector fromTrkToMu2eFrame(CLHEP::Hep3Vector  &vec);

        bool behindVane(double& posX, double& posY, int& vane);//{

        bool behindVane(HepPoint pos, int& vane);//{

        void caloExtrapol(int&          diagLevel,
			  int           evtNumber,
			  TrkFitDirection  fdir,
			  TrkRep const* trep,
			  double&       lowrange, 
			  double&       highrange,
			  HelixTraj     &trkHel, 
			  int           &res0, 
			  int&             NIntersections,
			  IntersectData_t* Intersections);

        void caloExtrapol(TrkRep const* trep, 
			  double& lowrange, 
			  double& highrange, 
			  HelixTraj &trkHel, 
			  int &res0, 
			  DetIntersection &intersec0, 
			  Length *pathLengths);

        void minimumPathLength(Length *length, int& vane, double& lowrange, double& highrange);
        
  void print(void){
    std::cout<<"Calorimeter 4 vanes info : "<<std::endl
	     <<"dR = "<< _dR<<std::endl
	     <<"dZ = "<< _dZ<<std::endl
	     <<"dU = "<< _dU<<std::endl
	     <<"dAPD = "<< _dAPD<<std::endl
	     <<"vaneHalfThickness = "<< _vaneHalfThickness<<std::endl
	     <<"solenoidOffSetX  = "<< _solenoidOffSetX<<std::endl
	     <<"solenoidOffSetZ  = "<< _solenoidOffSetZ<<std::endl
	     <<" innerRadius = "<< _innerRadius<<std::endl
	     <<"outherRadius = "<< _outherRadius<<std::endl
	     <<" ZfrontFaceCalo = "<< _ZfrontFaceCalo<<std::endl
	     <<"ZbackFaceCalo = "<<_ZbackFaceCalo<<std::endl
	     <<"nVanes  = "<<_nVanes<<std::endl;
  }
private :
        Hep3Vector _norm;
        double _dR;
        double _dZ;
        double _dU;
        double _dAPD;
        double _vaneHalfThickness;
        double _solenoidOffSetX;
        double _solenoidOffSetZ;
        double _innerRadius;
        double _outherRadius;
        double _ZfrontFaceCalo;
        double _ZbackFaceCalo;
        int _nVanes;


};
}


#endif
