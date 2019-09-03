//
// A Filter module aimed to select events using a Likelihood defined with calorimeter cluster info
//
// $Id: $
// $Author: $
// $Date: $
//
// Original author G. Pezzullo
//

#include "CLHEP/Units/SystemOfUnits.h"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"

#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"

#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"

#include "CaloCluster/inc/ClusterMoments.hh"

#include "RecoDataProducts/inc/CrvCoincidenceClusterCollection.hh"
#include "RecoDataProducts/inc/TrkFitFlag.hh"
#include "RecoDataProducts/inc/TriggerInfo.hh"

// #include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
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


  class CrvCoincidenceClusterFilter : public art::EDFilter {
     
  public:
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    struct Config {
      fhicl::Atom<int> diag{ Name("diagLevel"),
	Comment("Diagnostic Level"), 0};
      fhicl::Atom<art::InputTag> CRVCCF{ Name("CrvCoincidenceClusterFinder"),
	Comment("CrvCoincidenceClusterFinder producer")};
      fhicl::Atom<size_t> minNCC{ Name("MinNCC"),
	Comment("Minimum number of cluster of coincidences")};
      fhicl::Atom<std::string> trgPath{ Name("triggerPath"),
	Comment("label of the given trigger-path")};
    };

    virtual ~CrvCoincidenceClusterFilter() { }

    virtual void beginJob();
    virtual void endJob  ();
    virtual bool filter  (art::Event& event) override;
    virtual bool endRun( art::Run& run ) override;

    using Parameters = art::EDFilter::Table<Config>;
    explicit CrvCoincidenceClusterFilter(const Parameters& conf);

  private:
       
    typedef art::Ptr< CaloCrystalHit> CaloCrystalHitPtr;

    int                     _diagLevel;
    int                     _nProcess;
    int                     _nPass;
    art::InputTag           _clTag;
    size_t                  _minNCC;
    std::string             _trigPath;
    
  };


  CrvCoincidenceClusterFilter::CrvCoincidenceClusterFilter(const Parameters& config) :
    art::EDFilter{config},
    _diagLevel                   (config().diag()),  
    _nProcess                    (0),		     
    _nPass                       (0),		     
    _clTag                       (config().CRVCCF()),
    _minNCC                      (config().minNCC()),
    _trigPath                    (config().trgPath()){

    produces<TriggerInfo>();
  }

  void CrvCoincidenceClusterFilter::beginJob(){ }

  void CrvCoincidenceClusterFilter::endJob(){}

  bool CrvCoincidenceClusterFilter::endRun( art::Run& run ) {
    if(_diagLevel > 0 && _nProcess > 0){
      cout << "CrvCoincidenceClusterFilter" << " passed " <<  _nPass << " events out of " << _nProcess << " for a ratio of " << float(_nPass)/float(_nProcess) << endl;
    }
    return true;
  }
  
  //--------------------------------------------------------------------------------
  // Follow the body of the Filter logic
  //--------------------------------------------------------------------------------
  bool CrvCoincidenceClusterFilter::filter(art::Event& event) {

    ++_nProcess;
    if (_nProcess%10==0 && _diagLevel > 0) std::cout<<"Processing event from CrvCoincidenceClusterFilter =  "<<_nProcess  <<std::endl;
   
    unique_ptr<TriggerInfo> triginfo(new TriggerInfo);
    bool   retval(false);

    //Get calo cluster collection
    auto  clH = event.getValidHandle<CrvCoincidenceClusterCollection>(_clTag);
    const CrvCoincidenceClusterCollection*  crvCoincidenceClusters = clH.product();

    if (crvCoincidenceClusters->size() > _minNCC) retval = true;
    
    event.put(std::move(triginfo));
    return retval;
  }
 
}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::CrvCoincidenceClusterFilter);


