//
// This module transforms StrawDigi objects into StrawHit objects
// It also builds the truth match map (if MC truth info for the StrawDigis exists)
//
// $Id: PreloadDigiProducts_module.cc,v 1.12 2014/03/25 22:14:39 brownd Exp $
// $Author: brownd $ 
// $Date: 2014/03/25 22:14:39 $
//
// Original author David Brown, LBNL
//
// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "GeometryService/inc/GeomHandle.hh"
#include "art/Framework/Core/EDProducer.h"
#include "GeometryService/inc/DetectorSystem.hh"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "TrackerConditions/inc/StrawElectronics.hh"
#include "TrackerConditions/inc/StrawPhysics.hh"
// helpers
#include "TrkChargeReco/inc/PeakFit.hh"
#include "TrkChargeReco/inc/PeakFitRoot.hh"
#include "TrkChargeReco/inc/PeakFitFunction.hh"
#include "TrkChargeReco/inc/ComboPeakFitRoot.hh"
#include "TrackerMC/inc/SHInfo.hh"
//CLHEP
#include "CLHEP/Vector/ThreeVector.h"
// root 
#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TMarker.h"
#include "TTree.h"
// data
#include "RecoDataProducts/inc/StrawDigiCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/StrawDigiMCCollection.hh"

//#include "PerfLib/inc/perflib.hh"
//perf::PerfStats g_perf("PreloadDigiProducts 100") ;


using namespace std;
using CLHEP::Hep3Vector;
namespace mu2e {

  class PreloadDigiProducts : public art::EDProducer {

  public:
    explicit PreloadDigiProducts(fhicl::ParameterSet const& pset);
    virtual ~PreloadDigiProducts(); 
    virtual void beginJob();
    virtual void beginRun( art::Run& run );
    virtual void produce( art::Event& e);

  private:
    // Name of the StrawDigi collection
    string _strawDigis;
  };

  PreloadDigiProducts::PreloadDigiProducts(fhicl::ParameterSet const& pset) :
    _strawDigis(pset.get<string>("StrawDigis","makeSD"))
  {
   // produces<StrawHitCollection>();
   // produces<PtrStepPointMCVectorCollection>();
   // produces<StrawDigiMCCollection>();
  }

  PreloadDigiProducts::~PreloadDigiProducts() {}

  void PreloadDigiProducts::beginJob(){}

  void PreloadDigiProducts::beginRun( art::Run& run ){}

  void PreloadDigiProducts::produce(art::Event& event) { 
    // find the digis
 //   g_perf.read_begin_counters_inlined();

    art::Handle<mu2e::StrawDigiCollection> strawdigisH; 
    auto res [[gnu::unused]] = event.getByLabel(_strawDigis,strawdigisH);
    
    art::Handle<PtrStepPointMCVectorCollection> mcptrdigiH;    
    res=event.getByLabel(_strawDigis,mcptrdigiH);
    
    art::Handle<StrawDigiMCCollection> mcdigiH;
    res=event.getByLabel(_strawDigis,mcdigiH);

    
    //event.put(move(strawHits));

    //if(mcptrHits != 0)event.put(move(mcptrHits));
    //if(mchits != 0)event.put(move(mchits));
 //   g_perf.read_end_counters_inlined();
    
  }

}
using mu2e::PreloadDigiProducts;
DEFINE_ART_MODULE(PreloadDigiProducts);

