//
// $Id: Calorimeter4VanesGeom.hh,v 1.2 2012/07/10 04:54:49 gianipez Exp $
// $Author: gianipez $
// $Date: 2012/07/10 04:54:49 $
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

#include "KalmanTests/inc/TrkRecoTrkCollection.hh"

// From the art tool-chain
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//tracker includes
#include "TrkBase/TrkRecoTrk.hh"
#include "TrkBase/TrkRep.hh"
#include "KalmanTrack/KalRep.hh"
#include "KalmanTests/inc/KalFitMC.hh"
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
#include "BaBar/include/TrkBase/TrkExchangePar.hh"
#include "KalmanTests/inc/TrkDef.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"
#include "KalmanTests/inc/KalFit.hh"
#include "KalmanTests/inc/KalFitMC.hh"
#include "KalmanTests/inc/TrkRecoTrkCollection.hh"
#include "TrkPatRec/inc/TrkHitFilter.hh"
#include "TrkPatRec/inc/TrkHelixFit.hh"
#include "TrkBase/TrkPoca.hh"

//calorimeter includes
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "TrackCaloMatching/inc/CaloVolumeElem.hh"
#include "TrackCaloMatching/inc/CaloVolumeType.hh"
#include "TrackCaloMatching/inc/CaloSurface.hh"
#include "CaloCluster/inc/CaloClusterer.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
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
        //Construct from a transform.
        Calorimeter4VanesGeom(){
                GeomHandle<Calorimeter> cg;
                art::ServiceHandle<GeometryService> geom;
                _norm.setX(0.0);
                _norm.setY(1.0);
                _norm.setZ(0.0);
                _dR = cg->nCrystalR()*cg->crystalHalfSize();
                _dZ = cg->nCrystalZ()*cg->crystalHalfSize();
                _dU = cg->crystalHalfLength();
                _dAPD = cg->roHalfThickness();
                _vaneHalfThickness = _dU + _dAPD;
                _solenoidOffSetX = geom->config().getDouble("mu2e.solenoidOffset");//3904.;//[mm]
                _solenoidOffSetZ = -geom->config().getDouble("mu2e.detectorSystemZ0");//-10200.;
                _innerRadius = cg->innerRaidus();
                _outherRadius = cg->outherRadius();
                _ZfrontFaceCalo = cg->getOrigin().z() + _solenoidOffSetZ - (cg->nCrystalZ() + 1.0)*cg->crystalHalfSize();
                _ZbackFaceCalo = cg->getOrigin().z() + _solenoidOffSetZ + (cg->nCrystalZ() + 1.0)*cg->crystalHalfSize();
                _nVanes = cg->nVane();
        }


        //Destructor
        ~Calorimeter4VanesGeom ();

        CLHEP::Hep3Vector const& norm(){return _norm;}

        int nVanes() const{ return _nVanes;}

        double dR() const{ return _dR;}

        double dZ() const{ return _dZ;}

        double dU() const{ return _dU;}

        double dAPD() const{ return _dAPD;}

        double vaneHalfThickness() const{ return _vaneHalfThickness;}

        double solenoidOffSetX() const{ return _solenoidOffSetX;}

        double solenoidOffSetZ() const{ return _solenoidOffSetZ;}

        double innerRadius() const{ return _innerRadius;}

        double outherRadius() const{ return _outherRadius;}

        double ZfrontFaceCalo() const{ return _ZfrontFaceCalo;}

        double ZbackFaceCalo() const{ return _ZbackFaceCalo;}

        CaloVolumeElem* vane(int& i);//{

        bool behindVane(double& posX, double& posY, int& vane);//{

        bool behindVane(HepPoint pos, int& vane);//{

        void caloExtrapol(int& diagLevel,int evtNumber,/*TrkDifTraj const& traj*/ TrkRep const* trep,double& lowrange, double& highrange,
                        HelixTraj &trkHel, int &res0, DetIntersection &intersec0, Length *pathLengths/*, int& maxNumberExtrapolPoints*/);

        void caloExtrapol(TrkRep const* trep,double& lowrange, double& highrange, HelixTraj &trkHel, int &res0, DetIntersection &intersec0, Length *pathLengths);

        void minimumPathLength(Length *length, int& vane, double& lowrange, double& highrange);

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
