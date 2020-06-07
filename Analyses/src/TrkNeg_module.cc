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

  class TrkNeg : public art::EDAnalyzer {
    public:
      struct Config {
        using Name=fhicl::Name;
        using Comment=fhicl::Comment;
        fhicl::Sequence<art::InputTag> inputs{Name("KalRepPtrCollection"),Comment("outcome of Kalman filter (for tracker momentum info)")};
        
    };
    typedef art::EDAnalyzer::Table<Config> Parameters;

    explicit TrkNeg(const Parameters& conf);
    virtual ~TrkNeg() {};


    virtual void beginJob();
    virtual void analyze(const art::Event& e) override;

  private:
		
		Config _conf;
    typedef std::vector<art::InputTag> InputTags;
    InputTags inputs_;
    Int_t _evt, _run, _nEvts=0;
    TTree *_newNeg;
   
    Int_t  _nTracks;
   
    Float_t _TrackT0Neg, _TrackMomNeg,
    _TrackBackOmegaNeg ,_TrackBackD0Neg , _TrackBackZ0Neg, _TrackBackPhi0Neg,
    _TrackBackTanDipNeg, _TrackChi2Neg, _TrackChi2DOFNeg, _TrackNHitsNeg;
	};

  TrkNeg::TrkNeg(const Parameters& conf):
    art::EDAnalyzer(conf)
    {
      for(const auto& i : conf().inputs()) {
        inputs_.emplace_back(i);
      }
    }

  void TrkNeg::beginJob(){
    art::ServiceHandle<art::TFileService> tfs;
    _newNeg  = tfs->make<TTree>("TrkNegAnaNeg", "TrkNegAnaNeg");
    _newNeg->Branch("TrackT0Neg", 	&_TrackT0Neg, 		"TrackT0Neg/F");
    _newNeg->Branch("TrackBackOmegaNeg", &_TrackBackOmegaNeg,	"TrackBackOmegaNeg/F");
    _newNeg->Branch("TrackBackD0Neg",  	&_TrackBackD0Neg, 		"TrackBackD0Neg/F");
    _newNeg->Branch("TrackBackZ0Neg", 	&_TrackBackZ0Neg, 		"TrackBackZ0Neg/F");
    _newNeg->Branch("TrackBackPhi0Neg",	&_TrackBackPhi0Neg,	"TrackBackPhi0Neg/F");
    _newNeg->Branch("TrackBackTanDipNeg",&_TrackBackTanDipNeg,	"TrackBackTanDipNeg/F");
    _newNeg->Branch("TrackChi2Neg",	&_TrackChi2Neg,		"TrackChi2Neg/F");
    _newNeg->Branch("TrackChi2DOFNeg", 	&_TrackChi2DOFNeg, 	"TrackCho2DOFNeg/F");
    _newNeg->Branch("TrackMomNeg",	&_TrackMomNeg, 		"TrackMomNeg/F");
    _newNeg->Branch("TrackNHitsNeg",		&_TrackNHitsNeg,			"TrackNHitsNeg/F");

  }


  void TrkNeg::analyze(const art::Event& event) {
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
    
      
      if( TrackKrep->charge() == -1){
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

      }

      _nTracks ++;
   
      }

      _newNeg->Fill();    
  
 }
}

DEFINE_ART_MODULE(mu2e::TrkNeg);


