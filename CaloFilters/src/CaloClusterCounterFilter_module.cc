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
#include "art/Framework/Services/Optional/TFileService.h"
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
    
    enum {
      kN1DVar    = 10,
      kN2DVar    = 10,
      kNCorHist  = 10
    };

    virtual ~CaloClusterCounter() { }

    virtual void beginJob();
    virtual void endJob  ();
    virtual bool filter  (art::Event& event) override;
    virtual bool endRun( art::Run& run ) override;

    explicit CaloClusterCounter(const fhicl::ParameterSet& PSet);

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


  CaloClusterCounter::CaloClusterCounter(const fhicl::ParameterSet & pset) :
    _diagLevel                   (pset.get<int>("diagLevel",0)),
    _nProcess                    (0),
    _nPass                       (0),
    _clTag                       (pset.get<art::InputTag> ("CaloClusterModuleLabel")),
    _minClEnergy                 (pset.get<double>        ("MinClusterEnergy"     ,   50.)),   // MeV
    _minNCl                      (pset.get<int>           ("MinNClsuter"          ,   2)),   // 
    _trigPath                    (pset.get<std::string>("triggerPath")){

    produces<TriggerInfo>();
  }

  void CaloClusterCounter::beginJob(){ }

  void CaloClusterCounter::endJob(){}

  bool CaloClusterCounter::endRun( art::Run& run ) {
    if(_diagLevel > 0 && _nProcess > 0){
      cout << *currentContext()->moduleLabel() << " passed " <<  _nPass << " events out of " << _nProcess << " for a ratio of " << float(_nPass)/float(_nProcess) << endl;
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


