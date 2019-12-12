//Author: S Middleton
//Date: Dec 2019
//Purpose: Filter for Counter of Cosmics
#include "CLHEP/Units/SystemOfUnits.h"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"

#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"


#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"

#include "RecoDataProducts/inc/TrkFitFlag.hh"
#include "RecoDataProducts/inc/TriggerInfo.hh"
#include "RecoDataProducts/inc/CosmicTrackSeed.hh"
// #include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Provenance.h"

#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Utilities/InputTag.h"

#include "TDirectory.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TFile.h"

#include <cmath>
#include <string>
#include <vector>


using namespace std;

namespace mu2e {


  class CosmicCounter : public art::EDFilter {
     
  public:
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    struct Config {
      fhicl::Atom<int> diag{ Name("diagLevel"),
	  Comment("Diagnostic Level"), 0};
      fhicl::Atom<art::InputTag> CSTag{ Name("CosmicTrackSeedModuleLabel"),
	  Comment("CosmicTrackSeedModuleLabel producer")};
      fhicl::Atom<std::string> trgPath{ Name("triggerPath"),
	  Comment("label of the given trigger-path")};
    };

    virtual ~CaloClusterCounter() { }

    virtual void beginJob();
    virtual void endJob  ();
    virtual bool filter  (art::Event& event) override;
    virtual bool endRun( art::Run& run ) override;

    using Parameters = art::EDFilter::Table<Config>;
    explicit CaloClusterCounter(const Parameters& conf);

  private:
       
    typedef art::Ptr<CosmicTrackSeed> CosmicTrackSeedPtr;

    int                     _diagLevel;
    int                     _nProcess;
    int                     _nPass;
    art::InputTag           _csTag;
    std::string             _trigPath;
    
  };


  CosmicCounter::CaloClusterCounter(const Parameters& config):
    art::EDFilter{config},
    _diagLevel                   (config().diag()), 
    _nProcess                    (0),		     
    _nPass                       (0),		     
    _csTag                       (config().CSTag()),
    _trigPath                    (config().trgPath()){
      
      produces<TriggerInfo>();
    }
  
  void CosmicCounter::beginJob(){ }

  void CosmicCounter::endJob(){}

  bool CosmicCounter::endRun( art::Run& run ) {
    if(_diagLevel > 0 && _nProcess > 0){
      cout << "CosmicCounter" << " passed " <<  _nPass << " events out of " << _nProcess << " for a ratio of " << float(_nPass)/float(_nProcess) << endl;
    }
    return true;
  }
  
  //--------------------------------------------------------------------------------
  // Follow the body of the Filter logic
  //--------------------------------------------------------------------------------
  bool CosmicCounter::filter(art::Event& event) {

    ++_nProcess;
    if (_nProcess%10==0 && _diagLevel > 0) std::cout<<"Processing event from CosmicCounter =  "<<_nProcess  <<std::endl;
   
    unique_ptr<TriggerInfo> triginfo(new TriggerInfo);
    bool   retval(false);

    auto  cosH = event.getValidHandle<CosmicTrackSeedCollection>(_clTag);
    const CosmicTrackSeedCollection*  cosmics = cosH.product();

    int    nPass(0);

    for(auto icos = cosmics->begin();icos != cosmics->end(); ++icos){
      auto const& cosmic = *icos;
      //TODO - what to cut on here?
      ++nPass;
    }
      
    if (nPass>_npass) retval = true;  

    event.put(std::move(triginfo));
    return retval;
  }
 
}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::CosmicCounter);


