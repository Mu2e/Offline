//
// Filter events whit killed tracks.
// $Id: EcalTrigger_module.cc,v 1.2 2013/03/11 23:18:11 brownd Exp $
// $Author: brownd $
// $Date: 2013/03/11 23:18:11 $
//

// Mu2e includes.

#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"

#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"

#include "RecoDataProducts/inc/KalRepCollection.hh"
#include "RecoDataProducts/inc/TrkFitDirection.hh"
#include "RecoDataProducts/inc/TrkQual.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/TrkBase/TrkParticle.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/TrkBase/HelixParams.hh"
#include "BTrk/TrkBase/HelixTraj.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"

#include "RecoDataProducts/inc/KalRepPtrCollection.hh"

#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/GenId.hh"

#include "MCDataProducts/inc/StepPointMCCollection.hh"

#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"

#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "RecoDataProducts/inc/CaloRecoDigi.hh"
#include "RecoDataProducts/inc/CaloRecoDigiCollection.hh"

// prefetching Digi
#include "RecoDataProducts/inc/CaloDigi.hh"
#include "RecoDataProducts/inc/CaloDigiCollection.hh"

#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitIndex.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"

// Framework includes.
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "RecoDataProducts/inc/TrackClusterMatch.hh"

// C++ includes
#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <memory>
#include <vector>

using namespace std;

namespace mu2e {

  class EcalTriggerPreselect : public art::EDFilter {
  public:
    explicit EcalTriggerPreselect(fhicl::ParameterSet const& pset);
    virtual ~EcalTriggerPreselect() { }

    bool filter( art::Event& event);

    //virtual bool beginRun(art::Run &run) override;

    virtual void beginJob() override;
    virtual void endJob() override;

  private:

    int _diagLevel;
    std::string _g4ModuleLabel;
    std::string _virtualDetectorLabel;
    std::string _caloDigiModuleLabel;
    art::InputTag _shTag;   
    art::InputTag _shpTag;


    std::string _trkPatRecModuleLabel;
    art::InputTag _tqTag;
    std::string _caloClusterModuleLabel;
    std::string _trackClusterMatchModuleLabel;
    std::string _instanceName;
    TrkParticle _tpart;
    TrkFitDirection _fdir;
    SimParticleTimeOffset _toff;  // time offset smearing
    int nch; // number of combo hits
    const ComboHit*     comboHit;
    const ComboHitCollection* _chcol;
    const TrkQualCollection* _tqcol;
    int _virtualhit_good;
    int _trk_good;
    int _match_good;
    float _D0MIN;
    float _D0MAX;
    float _RMAXMIN;
    float _RMAXMAX;
    float _T0MIN;
    float _T0MAX;
    float _TANDIPMIN;
    float _TANDIPMAX;
    float _MVAMIN;
    float _VDPID;
    float _PMIN;
    float _PMATCHMIN;
    float _ECLUMIN;
    float _DTMIN;
    float _DTMAX;
    float _CHI2MATCHMIN;

    int                        _nProcessed;
    // Clusters
    int   _nCluster;
    float _cluEnergys[16384];
    float _cluene;
  };

  EcalTriggerPreselect::EcalTriggerPreselect(fhicl::ParameterSet const& pset):
    art::EDFilter{pset},
    _diagLevel(pset.get<int>("diagLevel",0)),
    _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel","g4run")),
    _virtualDetectorLabel(pset.get<std::string>("virtualDetectorName","virtualdetector")),
    _caloDigiModuleLabel      (pset.get<std::string>("caloDigiModuleLabel","CaloDigiFromShower")), 
    _shTag(pset.get<art::InputTag>("StrawHitCollectionTag","makePH")),
    _trkPatRecModuleLabel(pset.get<std::string>("trkPatRecModuleLabel")),
    _tqTag(pset.get<art::InputTag>("TrkQualTag","KFFDeM")),
    _caloClusterModuleLabel(pset.get<std::string>("caloClusterModuleLabel")), 
    _trackClusterMatchModuleLabel(pset.get<std::string>("trackClusterMatchModuleLabel")), 
    _tpart((TrkParticle::type)(pset.get<int>("fitparticle",TrkParticle::e_minus))),
    _fdir((TrkFitDirection::FitDirection)(pset.get<int>("fitdirection",TrkFitDirection::downstream))),
    _toff(pset.get<fhicl::ParameterSet>("TimeOffsets", fhicl::ParameterSet())),
    _virtualhit_good(pset.get<int>("virtualhit",1)),
    _trk_good(pset.get<int>("trkgood",1)),
    _match_good(pset.get<int>("matchgood",1)),
    _D0MIN                    (pset.get<float>("D0MIN",-80.)),
    _D0MAX                    (pset.get<float>("D0MAX",105.)),
    _RMAXMIN                  (pset.get<float>("RMAXMIN",450.)),
    _RMAXMAX                  (pset.get<float>("RMAXMAX",680.)),
    _T0MIN                    (pset.get<float>("T0MIN",700.)),
    _T0MAX                    (pset.get<float>("T0MAX",1720.)),
    _TANDIPMIN                (pset.get<float>("TANDIPMIN",0.57735)),
    _TANDIPMAX                (pset.get<float>("TANDIPMAX",1.)),
    _MVAMIN                   (pset.get<float>("MVAMIN",0.4)),
    _VDPID                    (pset.get<int>("VDPID",11)),
    _PMIN                     (pset.get<float>("PMIN",90.)),
    _PMATCHMIN                (pset.get<float>("PMATCHMIN",100.)),
    _ECLUMIN                  (pset.get<float>("ECLUMIN",10.)),
    _DTMIN                    (pset.get<float>("DTMIN",-12.)),
    _DTMAX                    (pset.get<float>("DTMAX",0.)),
    _CHI2MATCHMIN             (pset.get<float>("CHI2MATCHMIN",100.)),
    _nProcessed(0){
  }

  void EcalTriggerPreselect::beginJob(){
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory tfdir = tfs->mkdir("diag");
    
  }

  bool EcalTriggerPreselect::filter(art::Event& event) {

    
    ++_nProcessed;
    if (_nProcessed%100==0 && _diagLevel > 0) std::cout<<"Processing event from EcalTriggerPreselect =  "<<_nProcessed <<std::endl;

    //Handle to the calorimeter
    art::ServiceHandle<GeometryService> geom;
    if( ! geom->hasElement<Calorimeter>() ) return false;

    //load the timeoffset
    ConditionsHandle<AcceleratorParams> accPar("ignored");
    double _mbtime = accPar->deBuncherPeriod;
    _toff.updateMap(event);


    //--------------------------  Prefetch Calo Digis  --------------------------------    
    auto caloDigiFlag = event.getValidHandle<mu2e::CaloDigiCollection>(_caloDigiModuleLabel);
    if ( caloDigiFlag.product() == 0){
      std::cout << "EcalTriggerPreselect: CaloDigiHandle NOT VALID" << std::endl;
      return false;
    }
    const CaloDigiCollection& caloDigis = *caloDigiFlag.product();

    int roId;
    for (const auto& caloDigi : caloDigis){
      roId=caloDigi.roId();
      if (roId==0) continue;
    }

    // ------------------------- Prefetch Straw Hits and Straw Hit Positions -------------------------------
   
    auto chcolH = event.getValidHandle<mu2e::ComboHitCollection>(_shTag);

    if (chcolH.product() != 0){
      _chcol= chcolH.product();
      nch=_chcol->size();
      for(int istr=0; istr<nch;++istr) {
	comboHit=&_chcol->at(istr);
      }
    }
    else{
      std::cout << "EcalTriggerPreselect: ComboHitCollection Handle NOT VALID" << std::endl;
      return false;
    }
    if (_diagLevel > 1) std::cout << "Combo Hits Number: "  <<  nch <<endl;
    //------------------------------------------------------------------------------------
    // Handle ECAL virtual hits
    art::Handle<StepPointMCCollection> vdhits;
    event.getByLabel(_g4ModuleLabel,_virtualDetectorLabel,vdhits);
        
    // Handle tracks
    art::Handle<KalRepPtrCollection> trksHandle;
    event.getByLabel(_trkPatRecModuleLabel, trksHandle);
    const KalRepPtrCollection& trks = *trksHandle;
    //----- Track quality (produced by TrkPatRec/src/KalFinalFit_module.cc)
    // from TrkDiag/src/TrkRecoDiag_module.cc
    auto tqH = event.getValidHandle<TrkQualCollection>(_tqTag);
    _tqcol = tqH.product();

    // Handle clusters 
    art::Handle<CaloClusterCollection> caloClustersHandle;
    event.getByLabel(_caloClusterModuleLabel, caloClustersHandle);
    CaloClusterCollection const& caloClusters(*caloClustersHandle);

    _nCluster = 0;
    for (CaloClusterCollection::const_iterator clusterIt = caloClusters.begin(); clusterIt != caloClusters.end(); ++clusterIt){
      _cluEnergys[_nCluster] = clusterIt->energyDep();
      ++_nCluster;
    }

    // Check virtual hits

    bool virtualhit_good=false; 
    if (_virtualhit_good ){
      if (vdhits.isValid()){
	for (auto iter=vdhits->begin(), ie=vdhits->end(); iter!=ie; ++iter){
	  // the virtual hit has the same structure of the Step Point
	  const StepPointMC& hit = *iter;
	  // HITS in CALORIMETER VIRTUAL DETECTORS               
	  if (hit.volumeId()<VirtualDetectorId::EMC_Disk_0_SurfIn || hit.volumeId()>VirtualDetectorId::EMC_Disk_1_EdgeOut) continue;
	
	  double hitTimeUnfolded = _toff.timeWithOffsetsApplied(hit);
	  double hitTime         = fmod(hitTimeUnfolded,_mbtime);

	  if (hit.simParticle()->pdgId()==_VDPID && hitTime>(double) _T0MIN ){
	    virtualhit_good=true;
	    break;
	  }
	}
      }
      if (! virtualhit_good ) return false;
    }

    // Check on track cluster matching
    art::Handle<TrackClusterMatchCollection>  trackClusterHandle;
    event.getByLabel(_trackClusterMatchModuleLabel, trackClusterHandle);
    TrackClusterMatchCollection const& trackClusterMatches(*trackClusterHandle);
       
    bool trk_good=false; 
    if (_trk_good ){
      for ( size_t itrk=0; itrk< trks.size(); ++itrk ){
	KalRep const* trk = trks.at(itrk).get();
	CLHEP::Hep3Vector mom = trk->momentum(0);
	if ( (trk->helix(0).d0()                    > _D0MIN       ) &&
	     (trk->helix(0).d0()                    < _D0MAX       ) &&
	     (trk->helix(0).d0()+2./trk->helix(0).omega() > _RMAXMIN ) &&
	     (trk->helix(0).d0()+2./trk->helix(0).omega() < _RMAXMAX ) &&
	     (trk->t0().t0()                        > _T0MIN ) &&
	     (trk->t0().t0()                        < _T0MAX ) &&
	     (trk->helix(0).tanDip()                > _TANDIPMIN   ) &&
	     (trk->helix(0).tanDip()                < _TANDIPMAX   ) &&
	     (_tqcol->at(itrk).MVAOutput()          > _MVAMIN ) &&
	     (mom.mag()                             > _PMIN   )){
	  trk_good=true;
	}
      }
      if (_diagLevel > 0) std::cout << "Tracks Number: "  <<  trks.size() << " good?" << trk_good << endl;
      if ( ! trk_good ) return false;
    }

    bool match_good=false;  
    if (_match_good ){
      for (TrackClusterMatchCollection::const_iterator matchIt = trackClusterMatches.begin(); matchIt != trackClusterMatches.end(); ++matchIt){
	// good matching
	if (matchIt->iex()<(int)trks.size() && matchIt->iex()>=0 ){
	  KalRep const* trk = trks.at(matchIt->iex()).get();
	  CLHEP::Hep3Vector mom = trk->momentum(0);
	  if (matchIt->icl()<_nCluster && matchIt->icl()>=0 ){
	    _cluene= _cluEnergys[matchIt->icl()];
	    if (_diagLevel > 0) std::cout << "Cluster energy=" << _cluene << " p=" << mom.mag() << " matchIt->dt()=" << matchIt->dt() << " _tqcol=" <<_tqcol->at(matchIt->iex()).MVAOutput() <<  endl;
	    if (mom.mag()>_PMATCHMIN &&
		_cluene>_ECLUMIN &&
		matchIt->dt()>_DTMIN &&  matchIt->dt()<_DTMAX &&
		matchIt->chi2()<_CHI2MATCHMIN &&
		_tqcol->at(matchIt->iex()).MVAOutput()>_MVAMIN){
	      match_good=true;
	      break;
	    }
	  }
	}
      }
      if (_diagLevel > 0) std::cout << "Matches Number: "  << trackClusterMatches.size() << " good?" << match_good << endl;
      if (! match_good ) return false;
    }
    return true;
      
  }
  void EcalTriggerPreselect::endJob(){
    cout << "EcalTriggerPreselect filter end job:" << _nProcessed << " events processed" << endl;
  }
  
}

using mu2e::EcalTriggerPreselect;
DEFINE_ART_MODULE(EcalTriggerPreselect);
