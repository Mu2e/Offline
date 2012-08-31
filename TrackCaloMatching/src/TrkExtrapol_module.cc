//
//
//
// $Id: TrkExtrapol_module.cc,v 1.4 2012/08/31 22:34:53 brownd Exp $
// $Author: brownd $
// $Date: 2012/08/31 22:34:53 $
//
// Original author G. Pezzullo
//

// Framework includes.
#include "art/Framework/Core/EDProducer.h"
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
#include "BaBar/BaBar.hh"
#include "BaBar/Constants.hh"
#include "TrkBase/HelixParams.hh"
#include "TrkBase/TrkRep.hh"
#include "KalmanTrack/KalRep.hh"
#include "KalmanTests/inc/KalFitMC.hh"
// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TTrackerGeom/inc/TTracker.hh"
// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruth.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "TrackCaloMatching/inc/TrkToCaloExtrapolCollection.hh"


//calorimeter includes
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "TrackCaloMatching/inc/CaloVolumeElem.hh"
#include "TrackCaloMatching/inc/CaloVolumeType.hh"
#include "TrackCaloMatching/inc/CaloSurface.hh"
#include "CaloCluster/inc/CaloClusterer.hh"
#include "TrackCaloMatching/inc/Calorimeter4VanesGeom.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"


// Other includes.
#include "cetlib/exception.h"


// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/VisibleGenElTrack.hh"
#include "MCDataProducts/inc/VisibleGenElTrackCollection.hh"

//root includes
#include "TFile.h"
#include "TDirectory.h"
#include "TNtuple.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TGraph.h"
#include "TApplication.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TMath.h"

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


using namespace std;

namespace mu2e {

class TrkExtrapol : public art::EDProducer {
public:
        explicit TrkExtrapol(fhicl::ParameterSet const& pset):
        _fitterModuleLabel(pset.get<string>("fitterModuleLabel")),
        _kfitmc(pset.get<fhicl::ParameterSet>("KalFitMC")),
        _diagLevel(pset.get<int>("diagLevel",0)),
        _maxNumberStoresPoints(pset.get<int>("maxNumberStoresPoints",2)),
        _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),
        _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel", "g4run")),
        _caloReadoutModuleLabel(pset.get<std::string>("caloReadoutModuleLabel", "CaloReadoutHitsMaker")),
        _caloCrystalModuleLabel(pset.get<std::string>("caloCrystalModuleLabel", "CaloCrystalHitsMaker")),
        _application(0),
        _directory(0),
        _firstEvent(true){
                // Tell the framework what we make.
                produces<TrkToCaloExtrapolCollection>();
        }
        virtual ~TrkExtrapol() {
        }
        void beginJob() {}
        void endJob() {}

        void produce(art::Event & e );

private:

        void doExtrapolation(art::Event & evt, bool skip);
        // Module label of the module that performed the fits.
        std::string _fitterModuleLabel;
        // diagnostic of Kalman fit
        KalFitMC _kfitmc;

        // Diagnostic level
        int _diagLevel;

        int _maxNumberStoresPoints;

        // Label of the generator.
        std::string _generatorModuleLabel;

        // Label of the G4 module
        std::string _g4ModuleLabel;

        // Label of the calo readout hits maker
        std::string _caloReadoutModuleLabel;

        // Label of the calo crystal hists maker
        std::string _caloCrystalModuleLabel;

        std::auto_ptr<MCCaloUtilities> CaloManager;

        bool _skipEvent;

        // The job needs exactly one instance of TApplication.  See note 1.
        auto_ptr<TApplication> _application;

        // Save directory from beginJob so that we can go there in endJob. See note 3.
        TDirectory* _directory;
        bool _firstEvent;


};

void TrkExtrapol::produce(art::Event & evt ) {

        doExtrapolation(evt, _skipEvent);

} // end of analyze


void TrkExtrapol::doExtrapolation(art::Event & evt, bool skip){

        //create output
        auto_ptr<TrkToCaloExtrapolCollection> extrapolatedTracks(new TrkToCaloExtrapolCollection );
        TrkToCaloExtrapolCollection tmpExtrapolatedTracks;

        //Get handle to calorimeter
        art::ServiceHandle<GeometryService> geom;
        if(! geom->hasElement<Calorimeter>() ) return;
        GeomHandle<Calorimeter> cg;

        Calorimeter4VanesGeom *CaloVanes = new Calorimeter4VanesGeom();

        art::Handle<KalRepCollection> trksHandle;
        evt.getByLabel(_fitterModuleLabel,trksHandle);
        KalRepCollection const& trks = *trksHandle;

        if(_diagLevel>2){
                cout<<endl<<"Event Number : "<< evt.event()<< endl;
                cout<<"\n start TrkExtrapol..."<<endl;
                cout<<"trks.size() = "<< trks.size() <<endl;
        }

        double circleRadius = 0.0, centerCircleX=0.0, centerCircleY = 0.0, angle = 0.0;
        int res0 = -1;

        for ( size_t i=0; i< trks.size(); ++i ){

                KalRep const* trep = trks[i];
                if ( !trep ) continue;
                TrkDifTraj const& traj = trep->traj();
                double pos = 0.0;
                res0 = -1;

                double endTrk = trep->endFoundRange();
                HelixTraj trkHel(trep->helix(endTrk).params(),trep->helix(endTrk).covariance());//starting from the end of he tracker!!!FIXME

                angle = Constants::pi*0.5 + trkHel.phi0();

                circleRadius = 1.0/trkHel.omega();
                centerCircleX = trkHel.d0() + circleRadius;
                centerCircleY = centerCircleX*sin(angle);
                centerCircleX *= cos(angle);

                double lowrange = trkHel.zFlight(/*1740*/CaloVanes->ZfrontFaceCalo() ), highrange = trkHel.zFlight(/*3500*/ CaloVanes->ZbackFaceCalo() );

                DetIntersection intersec0;
                if(_diagLevel>2){

                        cout<<endl<<"Event Number : "<< evt.event()<< endl;
                        cout<<"------ trk number : "<<i<<" ------"<<endl;
                        cout<<"found traj, point of traj at "<<traj.position(pos)<<endl;
                        cout<<"*************** lowRange =  "<<trkHel.zFlight(1740)<<", highRange = "<<trkHel.zFlight(3500)<<endl;
                        cout<< " is the particle in the 0 quadrant?"<<endl<<
                                        ", traj.position(lowrange).x() = "<< traj.position(lowrange).x()<<
                                        ", traj.position(lowrange).y() = "<< traj.position(lowrange).y()<<
                                        ", traj.position(lowrange).z() = "<< traj.position(lowrange).z()<<endl;
                        cout<< ", circleRadius = "<< circleRadius<<
                                        ", centerCircleX = "<< centerCircleX<<
                                        ", centerCircleY = "<< centerCircleY<<endl;

                        cout<< " angle phi0 = " << angle<< endl;
                }

                const int nVanes = cg->nVane();
                Length length[nVanes];
                CaloVanes->caloExtrapol(_diagLevel, (int)evt.event(),trep, lowrange, highrange, trkHel,  res0, intersec0, length);

                cout<<"Event Number : "<< evt.event()<< endl;
                //}
                if(res0 != 1){

                        cout<< "ALLERT, intersection not found"<<endl;

                }else{
                        //if( _maxNumberStoresPoints>1){
                                for(int jVane=0; jVane <nVanes; ++jVane){

                                        for(Length::iterator it = length[jVane].begin(); it != length[jVane].end(); ++it){

                                                lowrange = it->first;
                                                highrange = it->second;

                                                KalRepPtr tmpRecTrk(trksHandle, i);
                                                tmpExtrapolatedTracks.push_back(
                                                                TrkToCaloExtrapol( jVane,tmpRecTrk,lowrange, highrange)
                                                );

                                                cout<< "Intersection found..."<<endl;
                                                //if(evt.event()%10==0){
                                                if(_diagLevel>2){
                                                        cout<<"vane intersected = "<< jVane <<
                                                                        ", entrance pathLength = "<<lowrange<<
                                                                        ", exit pathLength = "<<highrange<<endl;
                                                        cout<<"point of traj at entrance : "<<traj.position(lowrange)<<endl;
                                                        cout<<"errPoint of traj at entrance : "<<trep->positionErr(lowrange)<<endl;
                                                        cout<<"point of traj at the exit : "<<traj.position(highrange)<<endl;
                                                        cout<<"errPoint of traj at the exit : "<<trep->positionErr(highrange)<<endl;
                                                }

                                        }//end loop on pathLengths[]
                                }//end loop on vanes

                }
                std::sort(tmpExtrapolatedTracks.begin(), tmpExtrapolatedTracks.end());
                for(TrkToCaloExtrapolCollection::iterator it = tmpExtrapolatedTracks.begin(); it != tmpExtrapolatedTracks.end(); ++it){
                        extrapolatedTracks->push_back(*it);
                }

        }//end loop on recoTrj

        evt.put(extrapolatedTracks);

}

}

using mu2e::TrkExtrapol;
DEFINE_ART_MODULE(TrkExtrapol);
