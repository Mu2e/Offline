//
// a producer module that creates a track-calorimeter association with a payload
// of information belonging to that combination.  The association looks for 
// clusters within a configurable time window of the track tzero.  That allows
// us to see other clusters associated witht he process that produced the track.

    //
    // normally we would think  a track-calo is one-to-one, but we want a track 
    // to be associated with all clusters that might be related to the underlying event,
    // and sometime you fnd multiple tracks and want to associate a cluster with more than
    // one track.  So many-to-many.


//
// Original author Robert Bernstein
//
#include "RecoDataProducts/inc/TrackCaloAssns.hh"

#include "CLHEP/Units/SystemOfUnits.h"

#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "GlobalConstantsService/inc/unknownPDGIdName.hh"
#include "art/Framework/Core/EDProducer.h"

#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/CalorimeterPhysicalConstants.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"


#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "CaloCluster/inc/CaloContentMC.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/VirtualDetector.hh"

#include "BTrk/TrkBase/HelixTraj.hh"
#include "BTrk/TrkBase/HelixParams.hh"
#include "RecoDataProducts/inc/KalRepCollection.hh"
#include "RecoDataProducts/inc/TrkFitDirection.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
//
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"

#include "MCDataProducts/inc/CaloHitMCTruthCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/GenId.hh"
#include "MCDataProducts/inc/CaloHitSimPartMCCollection.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"

#include "Mu2eUtilities/inc/CaloHitMCNavigator.hh"

#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/TrkCaloIntersectCollection.hh"
#include "RecoDataProducts/inc/TrkCaloMatchCollection.hh"

#include "RecoDataProducts/inc/CaloHit.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"


// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/TrackerHitTimeClusterCollection.hh"

#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Provenance.h"
#include "canvas/Persistency/Common/Assns.h"
#include "cetlib/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Utilities/InputTag.h"

#include "TDirectory.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TProfile.h"

#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <vector>


#include "GeneralUtilities/inc/sqrtOrThrow.hh"

using namespace std;
using CLHEP::Hep3Vector;
using CLHEP::HepVector;
using CLHEP::keV;



namespace mu2e{
  class TrackCaloAssns : public art::EDProducer {
  public:

    explicit TrackCaloAssns(fhicl::ParameterSet const& pset);
    virtual ~TrackCaloAssns() { }

    virtual void beginJob() override;
    virtual void endJob() override;
    virtual void produce(art::Event& event) override;

  private:

    typedef std::vector< art::Handle<StepPointMCCollection> > HandleVector;
    typedef art::Ptr< CaloCrystalHit> CaloCrystalHitPtr;
    typedef art::Ptr<SimParticle> SimParticlePtr;
    int _diagLevel;
    int _nProcess;
    std::string _caloClusterModuleLabel;
    std::string      _trkCaloMatchModuleLabel;
    std::string      _trkIntersectModuleLabel;
    std::string      _trkFitterModuleLabel;
    std::string _trkPatRecModuleLabel;
    TrkParticle _tpart;
    TrkFitDirection _fdir;
    std::string _shLabel;
    std::string _shpLabel; 
    std::string _shfLabel;
    std::string _bkfLabel;
    std::string _trkfitInstanceName; 
    std::string _caloClusterAlgorithm;
    std::string _caloClusterSeeding;
    const std::string _producerName;
    std::string _instanceName;

    double _kinetic;
    double _kFrac;
    double _eOverP;
 
  



    int _timeDiff;

    TTree* _Ntup;

  };


  TrackCaloAssns::TrackCaloAssns(fhicl::ParameterSet const& pset) :
    _diagLevel(pset.get<int>("diagLevel",0)),
    _nProcess(0),
    _caloClusterModuleLabel(pset.get<std::string>("caloClusterModuleLabel")),
    _trkCaloMatchModuleLabel(pset.get<std::string>("trkCaloMatchModuleLabel")),
    _trkPatRecModuleLabel(pset.get<string>("trkPatRecModuleLabel")),
    _tpart((TrkParticle::type)(pset.get<int>("fitparticle",TrkParticle::e_minus))),
    _fdir((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection",TrkFitDirection::downstream))),
    _timeDiff(pset.get<int>("timeDifference",20)),
    _Ntup(0)
  {
    _instanceName = _fdir.name() + _tpart.name();
    _trkfitInstanceName = _fdir.name() + _tpart.name();
    produces<TrackCaloMatchAssns>(); // this must match the make_unique, it is the name of the data product, not the class...
   }

  void TrackCaloAssns::beginJob(){

    art::ServiceHandle<art::TFileService> tfs;

    _Ntup  = tfs->make<TTree>("TrackCaloAssns", "TrackCaloAssns");
  }
 
  void TrackCaloAssns::endJob(){}




  void TrackCaloAssns::produce(art::Event& event) {
 

    auto outputTrackCaloMatch = std::make_unique<TrackCaloMatchAssns>();

    ++_nProcess;
    if (_nProcess%10==0 && _diagLevel > 0) std::cout<<"Processing event from TrackCaloAssns =  "<<_nProcess << " with instance name " << _instanceName <<std::endl;

    if (_diagLevel > 0){std::cout << "******************new event*******************" << std::endl;}

    //get clusters
    art::Handle<CaloClusterCollection> caloClustersHandle;
    event.getByLabel(_caloClusterModuleLabel, caloClustersHandle);
    CaloClusterCollection const& caloClusters(*caloClustersHandle);


    //get tracks
    art::Handle<KalRepPtrCollection> krepsHandle;
    event.getByLabel(_trkPatRecModuleLabel,_instanceName,krepsHandle);
    KalRepPtrCollection const& kreps = *krepsHandle;  

    //get matches
    art::Handle<TrkCaloMatchCollection>  trkCaloMatchHandle;
    event.getByLabel(_trkCaloMatchModuleLabel, trkCaloMatchHandle);
    TrkCaloMatchCollection const& trkCaloMatches(*trkCaloMatchHandle);

 
	for (KalRepPtrCollection::const_iterator ithTrk=kreps.begin(); ithTrk != kreps.end(); ++ithTrk)
	  {
	    KalRep const& trk = **ithTrk;
	    //
	    // now KalRepCollections do not have art::Ptrs to the elements (legacy BaBar); I want them, so Rob K tells me to do this:

	    size_t trackIndex = ithTrk - kreps.begin();

	    //
	    //what is the t0 for this track? since this is a very loose cut we can use t0 at tracker center
	    //aka (0) in arrivalTime, pathlength = 0
	    double _trkTime = trk.arrivalTime(0);
	    double _trkMomentum = trk.momentum(0).mag();
	    //
	    // and loop through clusters, find time for each
	    for (CaloClusterCollection::const_iterator ithCluster = caloClusters.begin(); ithCluster != caloClusters.end(); ++ithCluster)
	      {
		size_t caloIndex = ithCluster - caloClusters.begin();
		double _cluTime = ithCluster->time();
		double _cluEnergy = ithCluster->energyDep();

		// 
		// are they close in time?  If so write out assn
		if (abs (_cluTime - _trkTime) < _timeDiff)
		  {

		    //
		    // see below, need this match info for payload
		    double _chi2TimeMatch(-1.);
		    double _chi2PosMatch(-1.);
		    double _chi2Match(-1.);

		    //
		    // calculate some useful things for payload
		    _eOverP = _cluEnergy/_trkMomentum ;
		    // 
		    // get kinetic energy; which particle?
		    double _kinetic = sqrt(_trkMomentum*_trkMomentum + _tpart.mass()*_tpart.mass()) 
		      - _tpart.mass();
		    _kFrac = _cluEnergy/_kinetic;

		    //in order to get chi2 to put into payload I need to pull out the matches
		    for (auto const& trkCaloMatch: trkCaloMatches)
		      {
			//
			// match this track only
			if (trkCaloMatch.trkId() == static_cast<int>(trackIndex) )
			  {
			    _chi2TimeMatch = trkCaloMatch.chi2Time();
			    _chi2PosMatch = trkCaloMatch.chi2Pos();
			    _chi2Match = trkCaloMatch.chi2();
			  }
		      }

		    //
		    // create ptrs to track and to cluster and to intersection objects
		    art::Ptr<KalRepPtr> _trkPtr(krepsHandle,trackIndex); 
		    art::Ptr<CaloCluster> _caloPtr(caloClustersHandle,caloIndex);
		    auto trkCaloAssnInfo = TrackCaloMatchInfo(_chi2Match,_chi2PosMatch,_chi2TimeMatch,_kFrac,_kinetic,_eOverP);
		    outputTrackCaloMatch->addSingle(_trkPtr,_caloPtr,trkCaloAssnInfo);
		  }
	      }
	  }

    event.put(std::move(outputTrackCaloMatch));
  }
}// end namespace mu2e
DEFINE_ART_MODULE(mu2e::TrackCaloAssns);


