//
//
//
// $Id: TrkExtrapol_module.cc,v 1.9 2013/03/27 18:37:24 murat Exp $
// $Author: murat $
// $Date: 2013/03/27 18:37:24 $
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
	_tpart((TrkParticle::type)(pset.get<int>("fitparticle",TrkParticle::e_minus))),
	_fdir((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection",TrkFitDirection::downstream))),
        _kfitmc(pset.get<fhicl::ParameterSet>("KalFitMC")),
        _diagLevel(pset.get<int>("diagLevel",0)),
        _maxNumberStoresPoints(pset.get<int>("maxNumberStoresPoints",2)),
        _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),
        _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel", "g4run")),
        _caloReadoutModuleLabel(pset.get<std::string>("caloReadoutModuleLabel", "CaloReadoutHitsMaker")),
        _caloCrystalModuleLabel(pset.get<std::string>("caloCrystalModuleLabel", "CaloCrystalHitsMaker")),
        CaloVanes(0),
        _application(nullptr),
        _directory(0),
        _firstEvent(true),
	_trkdiag(0){
	  // Tell the framework what we make.
	  produces<TrkToCaloExtrapolCollection>();
	  
	  // construct the data product instance name
	  _iname = _fdir.name() + _tpart.name();
	}
  
        virtual ~TrkExtrapol() {
                if(CaloVanes != 0){
                        delete CaloVanes;
                }
        }
        void beginJob(){

        }
        void endJob() {}

        void produce(art::Event & e );

private:

        void doExtrapolation(art::Event & evt, bool skip);
        // Module label of the module that performed the fits.
        std::string _fitterModuleLabel;
        
        TrkParticle _tpart;
        
        TrkFitDirection _fdir;
    
        std::string _iname;
        
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

        std::unique_ptr<MCCaloUtilities> CaloManager;

        Calorimeter4VanesGeom *CaloVanes;// = new Calorimeter4VanesGeom();

        bool _skipEvent;

        // The job needs exactly one instance of TApplication.  See note 1.
        unique_ptr<TApplication> _application;

        // Save directory from beginJob so that we can go there in endJob. See note 3.
        TDirectory* _directory;
        bool _firstEvent;

        TTree* _trkdiag;


};

void TrkExtrapol::produce(art::Event & evt ) {
        if(_firstEvent){
                CaloVanes =  new Calorimeter4VanesGeom();
		
		if(_diagLevel>0){
		  CaloVanes->print();
		}
		
                _firstEvent = false;
        }
        doExtrapolation(evt, _skipEvent);

} // end of analyze


//-----------------------------------------------------------------------------
void TrkExtrapol::doExtrapolation(art::Event & evt, bool skip){
  const char* oname = "TrkExtrapol::doExtrapolation";
  double lowrange, highrange, zmin, zmax;

  //create output
  auto_ptr<TrkToCaloExtrapolCollection> extrapolatedTracks(new TrkToCaloExtrapolCollection );
  TrkToCaloExtrapolCollection tmpExtrapolatedTracks;
  
  //Get handle to calorimeter
  art::ServiceHandle<GeometryService> geom;
  if(! geom->hasElement<VaneCalorimeter>() ) return;
  GeomHandle<VaneCalorimeter> cg;

  Calorimeter4VanesGeom *CaloVanes = new Calorimeter4VanesGeom();

  art::Handle<KalRepCollection> trksHandle;
  evt.getByLabel(_fitterModuleLabel,"DownstreameMinus",trksHandle);
  KalRepCollection const& trks = *trksHandle;
  
  if(_diagLevel>2){
    cout<<endl<<"Event Number : "<< evt.event()<< endl;
    cout<<"\n start TrkExtrapol..."<<endl;
    cout<<"trks.size() = "<< trks.size() <<endl;
  }

  double circleRadius = 0.0, centerCircleX=0.0, centerCircleY = 0.0, angle = 0.0;
  int res0 = -1;

  for ( size_t itrk=0; itrk< trks.size(); ++itrk ){

    KalRep const* trep = trks[itrk];
    if ( !trep ) continue;
    TrkDifTraj const& traj = trep->traj();
    double pos = 0.0;
    res0 = -1;
    
    double endTrk = trep->endFoundRange();
					// starting from the end of the tracker!!!FIXME
    HelixTraj trkHel(trep->helix(endTrk).params(),trep->helix(endTrk).covariance());

    angle = Constants::pi*0.5 + trkHel.phi0();

    circleRadius  = 1.0/trkHel.omega();
    centerCircleX = trkHel.d0() + circleRadius;
    centerCircleY = centerCircleX*sin(angle);
    centerCircleX *= cos(angle);

    // 2013-03-22 P.Murat : add 1cm tolerance
    zmin = CaloVanes->ZfrontFaceCalo()-10.;
    zmax = CaloVanes->ZbackFaceCalo ()+10.;
    
    lowrange  = trkHel.zFlight(zmin);  /*1740*/
    highrange = trkHel.zFlight(zmax); /*3500*/ 
    
    if (_diagLevel>2) {
      
      cout<<endl<<"Event Number : "<< evt.event()<< endl;
      cout<<"------ trk number : "<<itrk<<" ------"<<endl;
      cout<<"found traj, point of traj at "<<traj.position(pos)<<endl;
      cout<<"*************** lowRange =  "<< lowrange <<", highRange = "<< highrange << endl;
      cout<< " is the particle in the 0 quadrant?"<<endl<<
	", traj.position(lowrange).x() = "<< traj.position(lowrange).x()<<
	", traj.position(lowrange).y() = "<< traj.position(lowrange).y()<<
	", traj.position(lowrange).z() = "<< traj.position(lowrange).z()<<endl;

      printf("circle: R = %10.3f X0 = %10.3f  Y0 = %10.3f phi0 = %10.3f\n",
	     circleRadius,centerCircleX,centerCircleY,angle);
    }

    //    const int nVanes = cg->nVane();

    Calorimeter4VanesGeom::IntersectData_t  intersection[100];
    int                                     nint(0);

    CaloVanes->caloExtrapol(_diagLevel, (int)evt.event(), 
			    trep, lowrange, highrange, 
			    trkHel,  
			    res0, 
			    nint,
			    intersection);
    
    if (nint == 0) {
      printf("%s ERROR: intersection not found, res0 = %i\n",oname,res0);
    }

    for (int i=0; i<nint; i++) {
      KalRepPtr tmpRecTrk(trksHandle,itrk);
      tmpExtrapolatedTracks.push_back(
				      TrkToCaloExtrapol(intersection[i].fVane,
							tmpRecTrk,
							intersection[i].fSEntr,
							intersection[i].fSExit)
				      );
    }

    // P.Murat: why would one need to sort at this point?

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
