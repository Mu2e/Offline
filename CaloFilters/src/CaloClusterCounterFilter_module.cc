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

#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloHit.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "RecoDataProducts/inc/TrkFitFlag.hh"
#include "RecoDataProducts/inc/TriggerInfo.hh"

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


  class CaloClusterCounter : public art::EDFilter {
     
  public:
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    struct Config {
      fhicl::Atom<int> diag{ Name("diagLevel"),
	  Comment("Diagnostic Level"), 0};
      fhicl::Atom<art::InputTag> CCTag{ Name("CaloClusterModuleLabel"),
	  Comment("CaloClusterModuleLabel producer")};
      fhicl::Atom<int> minNClE{ Name("MinClusterEnergy"),
	  Comment("Minimum energy deposit in a calo cluster")};
      fhicl::Atom<int> minNCl{ Name("MinNCl"),
	  Comment("Minimum number of calo cluster")};
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
       
    typedef art::Ptr< CaloCrystalHit> CaloCrystalHitPtr;

    int                     _diagLevel;
    int                     _nProcess;
    int                     _nPass;
    art::InputTag           _clTag;
    double                  _minClEnergy;
    int                     _minNCl;
    std::string             _trigPath;
    
  };


  CaloClusterCounter::CaloClusterCounter(const Parameters& config):
    art::EDFilter{config},
    _diagLevel                   (config().diag()), 
    _nProcess                    (0),		     
    _nPass                       (0),		     
    _clTag                       (config().CCTag()),
    _minClEnergy                 (config().minNClE()),
    _minNCl                      (config().minNCl()),
    _trigPath                    (config().trgPath()){
      
      produces<TriggerInfo>();
    }
  
  void CaloClusterCounter::beginJob(){ }

  void CaloClusterCounter::endJob(){}

  bool CaloClusterCounter::endRun( art::Run& run ) {
    if(_diagLevel > 0 && _nProcess > 0){
      cout << "CaloClusterCounter" << " passed " <<  _nPass << " events out of " << _nProcess << " for a ratio of " << float(_nPass)/float(_nProcess) << endl;
    }
    return true;
  }
  
  //--------------------------------------------------------------------------------
  // Follow the body of the Filter logic
  //--------------------------------------------------------------------------------
  bool CaloClusterCounter::filter(art::Event& event) {

    ++_nProcess;
    if (_nProcess%10==0 && _diagLevel > 0) std::cout<<"Processing event from CaloClusterCounter =  "<<_nProcess  <<std::endl;
   
    unique_ptr<TriggerInfo> triginfo(new TriggerInfo);
    bool   retval(false);

    //Get calo cluster collection
    auto  clH = event.getValidHandle<CaloClusterCollection>(_clTag);
    const CaloClusterCollection*  caloClusters = clH.product();

    int    nClusterAboveThreshold(0);

    //for loop over the clusters in the calorimeter
    for(auto icl = caloClusters->begin();icl != caloClusters->end(); ++icl){
      auto const& cluster = *icl;
      double                   clEnergy = cluster.energyDep();
      if ( clEnergy < _minClEnergy)                       continue;
      ++nClusterAboveThreshold;
    }
      
    if (nClusterAboveThreshold>_minNCl) retval = true;  

    event.put(std::move(triginfo));
    return retval;
  }
 
}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::CaloClusterCounter);


