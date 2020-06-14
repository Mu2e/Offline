// Author : S. middleton
// Date: June 2020
// Purpose: RPC Combine Track Analysis
#include "CLHEP/Units/SystemOfUnits.h"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "GlobalConstantsService/inc/unknownPDGIdName.hh"

#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "TrackerGeom/inc/Tracker.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/VirtualDetector.hh"
#include "GeometryService/inc/DetectorSystem.hh"

#include "MCDataProducts/inc/CaloHitMCTruthCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"

#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/GenId.hh"
#include "MCDataProducts/inc/CaloHitSimPartMCCollection.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"

#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloHit.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "BTrkData/inc/TrkStrawHit.hh"
#include "RecoDataProducts/inc/TrkCaloIntersectCollection.hh"
#include "RecoDataProducts/inc/TrackClusterMatch.hh"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Provenance.h"
#include "cetlib_except/exception.h"
#include "GeneralUtilities/inc/ParameterSetHelpers.hh"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes.
#include "RecoDataProducts/inc/KalRepCollection.hh"
#include "RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "RecoDataProducts/inc/TrkFitDirection.hh"

// BaBar Kalman filter includes
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/TrkBase/HelixTraj.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/BbrGeom/TrkLineTraj.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrk/KalmanTrack/KalHit.hh"
#include "BTrk/TrkBase/HelixParams.hh"
#include "BTrk/BaBar/BaBar.hh"

// ROOT incldues
#include "TLegend.h"
#include "TLatex.h"
#include "TTree.h"
#include "TH2D.h"
#include "TF1.h"

#include "Rtypes.h"
#include "TApplication.h"
#include "TArc.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TNtuple.h"

#include "TStyle.h"
#include "TText.h"
#include "TRotMatrix.h"
#include "TColor.h"
#include "TLorentzVector.h"

#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <memory>
#include <vector>
#include <algorithm> 
using namespace std;
using CLHEP::Hep3Vector;
using CLHEP::keV;

#define me 0.511 //MeV

namespace mu2e {

  class TrkCombine : public art::EDAnalyzer {
    public:
      struct Config {
        using Name=fhicl::Name;
        using Comment=fhicl::Comment;
        fhicl::Sequence<art::InputTag> inputs{Name("KalRepPtrCollection"),Comment("outcome of Kalman filter (for tracker momentum info)")};
        
    };
    typedef art::EDAnalyzer::Table<Config> Parameters;

    explicit TrkCombine(const Parameters& conf);
    virtual ~TrkCombine() {};


    virtual void beginJob();
    virtual void analyze(const art::Event& e) override;

  private:
		
		Config _conf;
    typedef std::vector<art::InputTag> InputTags;
    InputTags inputs_;
    Int_t _evt, _run, _nEvts=0, _nNeg=0, _nPos=0;
    TTree* _Ntup;
   
    Int_t  _nTracks;
    Float_t _TrackT0, _TrackMom,  _TrackBackTime ,
    _TrackBackOmega ,_TrackBackD0 , _TrackBackZ0, _TrackBackPhi0,
    _TrackBackTanDip, _TrackChi2, _TrackChi2DOF, _TrackNHits, _Charge;
    Float_t _TrackT0Neg, _TrackMomNeg, _ChargeNeg,
    _TrackBackOmegaNeg ,_TrackBackD0Neg , _TrackBackZ0Neg, _TrackBackPhi0Neg,
    _TrackBackTanDipNeg, _TrackChi2Neg, _TrackChi2DOFNeg, _TrackNHitsNeg;
    Float_t _TrackT0Pos, _TrackMomPos,
    _TrackBackOmegaPos ,_TrackBackD0Pos , _TrackBackZ0Pos, _TrackBackPhi0Pos,
    _TrackBackTanDipPos, _TrackChi2Pos, _TrackChi2DOFPos, _TrackNHitsPos, _ChargePos;
	};

  TrkCombine::TrkCombine(const Parameters& conf):
    art::EDAnalyzer(conf)
    {
      for(const auto& i : conf().inputs()) {
        inputs_.emplace_back(i);
      }
    }

  void TrkCombine::beginJob(){
    art::ServiceHandle<art::TFileService> tfs;
    _Ntup  = tfs->make<TTree>("TrkCombineAna", "TrackCombineAna");
    _Ntup->Branch("evt",          	&_evt ,        "evt/I");
    _Ntup->Branch("run",          	&_run ,        "run/I");
    _Ntup->Branch("nTracks",	&_nTracks,		"nTracks/I");
    _Ntup->Branch("TrackT0", 	&_TrackT0, 		"TrackT0/F");
    _Ntup->Branch("TrackBackTime", 	&_TrackBackTime,   	"TrackBackTime/F");
    _Ntup->Branch("TrackBackOmega", &_TrackBackOmega,	"TrackBackOmega/F");
    _Ntup->Branch("TrackBackD0",  	&_TrackBackD0, 		"TrackBackD0/F");
    _Ntup->Branch("TrackBackZ0", 	&_TrackBackZ0, 		"TrackBackZ0/F");
    _Ntup->Branch("TrackBackPhi0",	&_TrackBackPhi0,	"TrackBackPhi0/F");
    _Ntup->Branch("TrackBackTanDip",&_TrackBackTanDip,	"TrackBackTanDip/F");
    _Ntup->Branch("TrackChi2",	&_TrackChi2,		"TrackChi2/F");
    _Ntup->Branch("TrackChi2DOF", 	&_TrackChi2DOF, 	"TrackCho2DOF/F");
    _Ntup->Branch("TrackMom",	&_TrackMom, 		"TrackMom/F");
    _Ntup->Branch("TrackNHits",		&_TrackNHits,			"TrackNHits/F");
    _Ntup->Branch("Charge",		&_Charge,			"Charge/F");

    _Ntup->Branch("TrackT0Neg", 	&_TrackT0Neg, 		"TrackT0Neg/F");
    _Ntup->Branch("TrackBackOmegaNeg", &_TrackBackOmegaNeg,	"TrackBackOmegaNeg/F");
    _Ntup->Branch("TrackBackD0Neg",  	&_TrackBackD0Neg, 		"TrackBackD0Neg/F");
    _Ntup->Branch("TrackBackZ0Neg", 	&_TrackBackZ0Neg, 		"TrackBackZ0Neg/F");
    _Ntup->Branch("TrackBackPhi0Neg",	&_TrackBackPhi0Neg,	"TrackBackPhi0Neg/F");
    _Ntup->Branch("TrackBackTanDipNeg",&_TrackBackTanDipNeg,	"TrackBackTanDipNeg/F");
    _Ntup->Branch("TrackChi2Neg",	&_TrackChi2Neg,		"TrackChi2Neg/F");
    _Ntup->Branch("TrackChi2DOFNeg", 	&_TrackChi2DOFNeg, 	"TrackCho2DOFNeg/F");
    _Ntup->Branch("TrackMomNeg",	&_TrackMomNeg, 		"TrackMomNeg/F");
    _Ntup->Branch("TrackNHitsNeg",		&_TrackNHitsNeg,			"TrackNHitsNeg/F");
    _Ntup->Branch("ChargeNeg",		&_ChargeNeg,			"ChargeNeg/F");
    
    _Ntup->Branch("TrackT0Pos", 	&_TrackT0Pos, 		"TrackT0Pos/F");
    _Ntup->Branch("TrackBackOmegaPos", &_TrackBackOmegaPos,	"TrackBackOmegaPos/F");
    _Ntup->Branch("TrackBackD0Pos",  	&_TrackBackD0Pos, 		"TrackBackD0Pos/F");
    _Ntup->Branch("TrackBackZ0Pos", 	&_TrackBackZ0Pos, 		"TrackBackZ0Pos/F");
    _Ntup->Branch("TrackBackPhi0Pos",	&_TrackBackPhi0Pos,	"TrackBackPhi0Pos/F");
    _Ntup->Branch("TrackBackTanDipPos",&_TrackBackTanDipPos,	"TrackBackTanDipPos/F");
    _Ntup->Branch("TrackChi2Pos",	&_TrackChi2Pos,		"TrackChi2Pos/F");
    _Ntup->Branch("TrackChi2DOFPos", 	&_TrackChi2DOFPos, 	"TrackCho2DOFPos/F");
    _Ntup->Branch("TrackMomPos",	&_TrackMomPos, 		"TrackMomPos/F");
    _Ntup->Branch("TrackNHitsPos",		&_TrackNHitsPos,			"TrackNHitsPos/F");
    _Ntup->Branch("ChargePos",		&_ChargePos,			"ChargePos/F");
  }


  void TrkCombine::analyze(const art::Event& event) {
    _evt = event.id().event();
    _run = event.run();
    art::ServiceHandle<mu2e::GeometryService>   geom;
    mu2e::GeomHandle<mu2e::DetectorSystem>      ds;
    mu2e::GeomHandle<mu2e::VirtualDetector>     vdet;
    Hep3Vector vd_tt_back = ds->toDetector(vdet->getGlobal(mu2e::VirtualDetectorId::TT_Back));
    double     Z      = vd_tt_back.z();
    _nTracks = 0;
    _nEvts++;
    for(const auto& cn : inputs_) {
      auto ih = event.getValidHandle<KalRepPtrCollection>(cn);
      if(ih->empty()) continue;
      art::Ptr<KalRep> const& ptr = ih->at(0);
      const KalRep* TrackKrep = ptr.get();
      double  ds(10.), s0, s1, s2, z0, z1, z2, dzds, sz, sz1, z01;
      const TrkHitVector* hots = &TrackKrep->hitVector();
      int nh = hots->size();
      const TrkHit *first(nullptr), *last(nullptr);

      for (int ih=0; ih<nh; ++ih) {
        const TrkHit* hit = hots->at(ih);
        if (hit  != nullptr) {
          if (first == nullptr) first = hit;
            last = hit;
          }
      }

      s1 = first->fltLen();
      s2 = last ->fltLen();

      z1     = TrackKrep->position(s1).z();
      z2     = TrackKrep->position(s2).z();

      dzds   = (z2-z1)/(s2-s1);

      if (fabs(Z-z1) > fabs(Z-z2)) {
        z0 = z2;
        s0 = s2;
      }
      else {
        z0 = z1;
        s0 = s1;
      }

      sz    = s0+(Z-z0)/dzds;

      z0     = TrackKrep->position(sz).z();     
      z01    = TrackKrep->position(sz+ds).z();

      dzds   = (z01-z0)/ds;
      sz1    = sz+(Z-z0)/dzds;	          
      HelixParams helx  = TrackKrep->helix(sz1);
      _TrackT0 = TrackKrep->t0().t0();
      _TrackMom = TrackKrep->momentum(sz1).mag();
      _TrackBackTime =   TrackKrep->arrivalTime(sz1);
      _TrackBackOmega    = helx.omega();
      _TrackBackD0       = helx.d0();
      _TrackBackZ0       = helx.z0();
      _TrackBackPhi0     = helx.phi0();
      _TrackBackTanDip  = helx.tanDip();
      _TrackChi2 =TrackKrep->chisq();
      _TrackChi2DOF= TrackKrep->chisq()/TrackKrep->nActive();
      _TrackNHits = TrackKrep->nActive();
      _Charge = TrackKrep->charge();
      if(_Charge == -1){
        _nNeg++;
        _TrackT0Neg = TrackKrep->t0().t0();
        _TrackMomNeg = TrackKrep->momentum(sz1).mag();
        _TrackBackOmegaNeg       = helx.omega();
        _TrackBackD0Neg       = helx.d0();
        _TrackBackZ0Neg       = helx.z0();
        _TrackBackPhi0Neg     = helx.phi0();
        _TrackBackTanDipNeg   = helx.tanDip();
        _TrackChi2Neg =TrackKrep->chisq();
        _TrackChi2DOFNeg= TrackKrep->chisq()/TrackKrep->nActive();
        _TrackNHitsNeg = TrackKrep->nActive();
        _ChargeNeg = _Charge;
        std::cout<<"Event "<<_evt<<" Run "<<_run<<" Track "<<_nTracks<<" Charge "<<_Charge<<" Neg D0 "<<_TrackBackD0Neg<<" phi0 "<<_TrackBackPhi0Neg<<" Chi2 dof "<<_TrackChi2DOFNeg<<" Mom "<<_TrackMomNeg<<" T0 "<<_TrackT0Neg<<std::endl;
      }
      
      if(_Charge == 1){
        _nPos++;
        _TrackT0Pos = TrackKrep->t0().t0();
        _TrackMomPos = TrackKrep->momentum(sz1).mag();
        _TrackBackOmegaPos       = helx.omega();
        _TrackBackD0Pos       = helx.d0();
        _TrackBackZ0Pos       = helx.z0();
        _TrackBackPhi0Pos     = helx.phi0();
        _TrackBackTanDipPos   = helx.tanDip();
        _TrackChi2Pos =TrackKrep->chisq();
        _TrackChi2DOFPos= TrackKrep->chisq()/TrackKrep->nActive();
        _TrackNHitsPos = TrackKrep->nActive();
        _ChargePos = _Charge;
        std::cout<<"Event "<<_evt<<" Run "<<_run<<" Track "<<_nTracks<<" Charge "<<_Charge<<" Neg D0 "<<_TrackBackD0Pos<<" phi0 "<<_TrackBackPhi0Pos<<" Chi2 dof "<<_TrackChi2DOFPos<<" Mom "<<_TrackMomPos<<" T0 "<<_TrackT0Pos<<std::endl;
      }

     
      _nTracks ++;
   
      }

      _Ntup->Fill();
 
  
 }
}

DEFINE_ART_MODULE(mu2e::TrkCombine);


