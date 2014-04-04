//
// A module to create simple stereo hits out of StrawHits.  This can work
// with either tracker.  StrawHit selection is done by flagging in an upstream module
//
// $Id: MakeStrawHitsNew_module.cc,v 1.3 2014/04/04 21:23:34 murat Exp $
// $Author: murat $
// $Date: 2014/04/04 21:23:34 $
// 
//  Original Author: David Brown, LBNL
//  

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruth.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"

#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "CaloCluster/inc/CaloClusterTools.hh"
#include "CaloCluster/inc/CaloClusterUtilities.hh"


// art includes.
#include "art/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileService.h"

// From the art tool-chain
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// root
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLine.h"
#include "TMarker.h"
#include "TList.h"
#include "TLegend.h"

// C++ includes.
#include <iostream>
#include <float.h>

using namespace std;

namespace mu2e {

  class MakeStrawHitsNew : public art::EDProducer {

  public:
    explicit MakeStrawHitsNew(fhicl::ParameterSet const& pset);
    // Accept compiler written d'tor.
    virtual ~MakeStrawHitsNew()
    {
      delete _cluTool;
    }
    void produce( art::Event& e);
    void beginJob();

  private:

    // Diagnostics level.
    int _diagLevel;
    
    // Name of the StrawHit collection
    string _shLabel;
    
    // Label of the calo clusters  maker
    std::string  _caloClusterModuleLabel;
    std::string  _caloClusterAlgorithm;
    std::string  _caloClusterSeeding;
    std::string  _producerName;

    int _nbins;
    double _tmin, _tmax;

 // Parameters
    double _maxEemc;
    double _t1, _t2; // maximum time separation between hits
    double _pitchAngle;

    
    //histograms for diag level
    TH1F* _oldHits;
    TH1F* _newHits;
    CaloClusterTools *_cluTool;
    
    double _solenoidOffSetZ;
    double _zFirstDisk;
    double _zSecondDisk;
    
    void doTimeDistance(double zcal, double zstraw, double &distace);
  };

  void MakeStrawHitsNew::doTimeDistance(double zcal, double zstraw, double &distance)
  {
    //printf("\nzcal = %.1f, \n zstraw = %.1f", zcal, zstraw);
    distance = (zcal-zstraw)*CLHEP::mm;
    //printf("\n zcal-zstraw = %.1f", distance);
    distance /= CLHEP::m / CLHEP::mm;
    //printf("\nconverted in to m = %.1f", distance);
    distance /= std::sin(_pitchAngle);
    //printf("\n /sinpitchAngle = %.1f", distance);
    distance /= CLHEP::c_light;
    //printf("\n /c = %.1e", distance);
    distance /= CLHEP::ms / CLHEP::s;
    //printf("\n converted in to ns = %.1e\n", distance);
  }

  MakeStrawHitsNew::MakeStrawHitsNew(fhicl::ParameterSet const& pset) :

    // Parameters
    _diagLevel(pset.get<int>("diagLevel",0)),
    _shLabel(pset.get<string>("StrawHitCollectionLabel","makeSH")),
    _caloClusterModuleLabel(pset.get<std::string>("caloClusterModuleLabel","CaloClusterMakerNew")),
    _caloClusterAlgorithm(pset.get<std::string>("caloClusterAlgorithm", "closest")),
    _caloClusterSeeding(pset.get<std::string>("caloClusterSeeding", "energy")),
    _producerName("Algo"+mu2e::TOUpper(_caloClusterAlgorithm)+"SeededBy"+mu2e::TOUpper(_caloClusterSeeding)),
    _nbins(pset.get<int>("histbins", 500)),
    _tmin(pset.get<double>("tmin", -1000)),
    _tmax(pset.get<double>("tmax", 1000)),
    _maxEemc(pset.get<double>("maxEemc",70.0)), // nsec
    _t1(pset.get<double>("t1",-70.0)), // nse
    _t2(pset.get<double>("t2",20.0)), // nsec
    _pitchAngle(pset.get<double>("pitchAngle", 0.67)),
    _oldHits(0),_newHits(0)
  {
    // Tell the framework what we make.
    produces<StrawHitCollection>();
    //03 - 02 - 2014 Gianipez added the following change for alligning the code with Dave changes
    //    produces<StrawHitMCTruthCollection>();
    produces<PtrStepPointMCVectorCollection>("StrawHitMCPtr");
  }

   void MakeStrawHitsNew::beginJob(){
    // create diagnostics if requested
     if(_diagLevel > 0)
       {
	 art::ServiceHandle<art::TFileService> tfs;
	 _oldHits = tfs->make<TH1F>("oldHits","Time distall hits;t_{EMC}-t_{straw}", _nbins, _tmin, _tmax);
	 _newHits = tfs->make<TH1F>("newHits","Time dist selected hits;t_{EMC}-t_{straw}", _nbins, _tmin, _tmax);
       }

     

  }

  
  void
  MakeStrawHitsNew::produce(art::Event& event) {

    if ( _diagLevel > 0 && (event.id().event() % 100 ==0)  ) cout << "MakeStrawHitsNew: produce() begin; event " << event.id().event() << endl;

    //data about clusters
    art::Handle<CaloClusterCollection> caloClustersHandle;
    event.getByLabel(_caloClusterModuleLabel,_producerName, caloClustersHandle);
    CaloClusterCollection const& caloClusters(*caloClustersHandle);
    
    art::Handle<mu2e::StrawHitCollection> strawhitsH; 
    const StrawHitCollection* strawhits(0);
    if(event.getByLabel(_shLabel,strawhitsH))
      strawhits = strawhitsH.product();
    if(strawhits == 0){
      throw cet::exception("RECO")<<"mu2e::MakeStrawHitsNew: No StrawHit collection found for label " <<  _shLabel << endl;
    } 

    //03 - 02 - 2014 Gianipez added the following change for alligning the code with Dave changes
    // art::Handle<StrawHitMCTruthCollection> truthHandle;
//     event.getByLabel(_shLabel, truthHandle);
//     const StrawHitMCTruthCollection* hits_truth = truthHandle.product();

    art::Handle<PtrStepPointMCVectorCollection> mcptrHandle;
    event.getByLabel(_shLabel,"StrawHitMCPtr",mcptrHandle);
    const PtrStepPointMCVectorCollection* hits_mcptr = mcptrHandle.product();


    art::ServiceHandle<GeometryService> geom;
    _solenoidOffSetZ = -geom->config().getDouble("mu2e.detectorSystemZ0");
    
    GeomHandle<DiskCalorimeter> cgDisks;
    _zFirstDisk = cgDisks->origin().z() + cgDisks->diskSeparation(0) + _solenoidOffSetZ;
    _zSecondDisk  = cgDisks->origin().z() + cgDisks->diskSeparation(1) + _solenoidOffSetZ;
    if(_diagLevel > 1)
      {
	printf("\nSolenoid offset z = %.2f, \n z first disk = %.2f, \n z first disk = %.2f\n", _solenoidOffSetZ,_zFirstDisk,_zSecondDisk);
      }
    
    //loop over the reconstructed emc-clusters for finding the most energetic one
    // and get the time from that for defining the time window
    CaloCluster cluster;
    double EemcMax( _maxEemc), timeEMC(-9999.0), tmpEnergy(0.0);
    double zcal;
    size_t nEMCclusters(caloClusters.size() );


    for(size_t i=0; i < nEMCclusters ; ++i)
      {
	cluster = caloClusters.at(i);
	_cluTool = new CaloClusterTools(cluster);
	tmpEnergy = cluster.energyDep();
	if( tmpEnergy >= EemcMax ) 
	  {
	    timeEMC = _cluTool->timeFasterCrystal();
	    EemcMax = tmpEnergy;
	    if(cluster.vaneId() == 0)
	      {
		zcal = _zFirstDisk;
	      } else
	      {
		zcal = _zSecondDisk;
	      }
	  }
      }
    


    // create a collection of StrawHitPosition, and intialize them using the time division
    size_t nsh = strawhits->size();

    unique_ptr<StrawHitCollection> shcol(new StrawHitCollection);
    //03 - 02 - 2014 Gianipez added the following change for alligning the code with Dave changes
    //    unique_ptr<StrawHitMCTruthCollection>      truthHits(new StrawHitMCTruthCollection);
    unique_ptr<PtrStepPointMCVectorCollection> mcptrHits(new PtrStepPointMCVectorCollection);
    

    const Tracker& tracker = getTrackerOrThrow();
    double timeToCal(-9999.);
    double deltaTime(-9999);
    //    if(timeEMC < 0.0) goto NEXT_EVENT;

     
    for(size_t ish=0;ish<nsh;++ish)
      {
	StrawHit const& hit = strawhits->at(ish);
	Straw const& straw = tracker.getStraw(hit.strawIndex());


	doTimeDistance(zcal, straw.getMidPoint().z(), timeToCal);
	deltaTime = timeEMC - hit.time() - timeToCal;
	
	if(_diagLevel > 0)
	  {
	    _oldHits->Fill(deltaTime);
	  }

	if ((timeEMC < 0) || ( (deltaTime >= _t1) && (deltaTime <= _t2))) 
	  {
	    shcol->push_back(StrawHit(hit.strawIndex(), 
				      hit.time(),
				      hit.dt(),
				      hit.energyDep()) );
    //03 - 02 - 2014 Gianipez added the following change for alligning the code with Dave changes    
	    // truthHits->push_back(StrawHitMCTruth(hits_truth->at(ish).driftTime(),
// 						 hits_truth->at(ish).driftDistance(),
// 						 hits_truth->at(ish).distanceToMid()) );
	    mcptrHits->push_back(PtrStepPointMCVector(hits_mcptr->at(ish)) );
	    if(_diagLevel > 0)
	      {
		_newHits->Fill(deltaTime);
	      }
	  }
      }
    
    // NEXT_EVENT:;
    event.put(std::move(shcol));
    //03 - 02 - 2014 Gianipez added the following change for alligning the code with Dave changes
    //event.put(std::move(truthHits)); 
    event.put(std::move(mcptrHits),"StrawHitMCPtr");
  } // end MakeStrawHitsNew::produce.

} // end namespace mu2e

using mu2e::MakeStrawHitsNew;
DEFINE_ART_MODULE(MakeStrawHitsNew)

