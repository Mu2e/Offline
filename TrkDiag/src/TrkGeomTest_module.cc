//
// This module analyzes the tracker alignment
//
// Original author David Brown, 11/16/2020 LBNL
//
// framework
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "GeometryService/inc/GeomHandle.hh"
#include "ProditionsService/inc/ProditionsHandle.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "TTree.h"
#include <iostream>

namespace mu2e {

  class TrkGeomTest : public art::EDAnalyzer
  {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      struct Config {
	fhicl::Atom<int> printLevel{ Name("printLevel"), Comment("Print level" ),0 };
	fhicl::Atom<int> diagLevel{ Name("diagLevel"), Comment("Diagnostic level" ),0 };
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit TrkGeomTest(Parameters const& config);
      virtual ~TrkGeomTest();
      void analyze( art::Event const & e) override;
      void beginJob ( ) override;
    private:
      ProditionsHandle<Tracker> _alignedTracker_h;
      TTree* trkalign_;
      int print_, diag_;
      int plane_, panel_, straw_;
      float nomx_, nomy_, nomz_;
      float nomdx_, nomdy_, nomdz_;
      float rnom_, phinom_, dphinom_;
      float deltax_, deltay_, deltaz_;
      bool first_;
  };

  TrkGeomTest::TrkGeomTest(Parameters const& config) :     art::EDAnalyzer(config),
  trkalign_(0),
  print_(config().printLevel()),
  diag_(config().diagLevel()),
  first_(true) {}

  TrkGeomTest::~TrkGeomTest() {}

  void TrkGeomTest::beginJob(){ 
    if(diag_ > 0){
      art::ServiceHandle<art::TFileService> tfs;
      trkalign_=tfs->make<TTree>("trkalign","tracker alignment");
      trkalign_->Branch("plane",&plane_,"plane/I");
      trkalign_->Branch("panel",&panel_,"panel/I");
      trkalign_->Branch("straw",&straw_,"straw/I");
      trkalign_->Branch("nomx",&nomx_,"nomx/F");
      trkalign_->Branch("nomy",&nomy_,"nomy/F");
      trkalign_->Branch("nomz",&nomz_,"nomz/F");
      trkalign_->Branch("nomdx",&nomdx_,"nomdx/F");
      trkalign_->Branch("nomdy",&nomdy_,"nomdy/F");
      trkalign_->Branch("nomdz",&nomdz_,"nomdz/F");
      trkalign_->Branch("rnom",&rnom_,"rnom/F");
      trkalign_->Branch("phinom",&phinom_,"phinom/F");
      trkalign_->Branch("dphinom",&dphinom_,"dphinom/F");
      trkalign_->Branch("deltax",&deltax_,"deltax/F");
      trkalign_->Branch("deltay",&deltay_,"deltay/F");
      trkalign_->Branch("deltaz",&deltaz_,"deltaz/F");
    }
  }

  void TrkGeomTest::analyze(art::Event const& event) {
    if(first_){
      // fetch aligned and nominal tracker objects
      GeomHandle<Tracker> nominalTracker_h;
      auto const& ntracker = *nominalTracker_h;
      auto const& atracker = _alignedTracker_h.get(event.id());
      if(ntracker.straws().size() != atracker.straws().size()){
	std::cout << "Trackers don't match" << std::endl;
      } else {
	for(size_t istr = 0; istr < ntracker.straws().size(); istr++){
	  auto const& nstraw = ntracker.straws()[istr];
	  auto const& astraw = atracker.straws()[istr];
	  if(nstraw.id() != astraw.id())
	    std::cout << "StrawIds don't match: nominal " << nstraw.id() << " aligned " << astraw.id() << std::endl;
	  plane_ = nstraw.id().plane();
	  panel_ = nstraw.id().panel();
	  straw_ = nstraw.id().straw();
	  auto npos = nstraw.origin();
	  auto ndir = nstraw.wireDirection();
	  auto apos = astraw.origin();
	  auto delta = apos-npos;
	  nomx_ = npos.x(); nomy_ = npos.y(); nomz_ = npos.z();
	  nomdx_ = ndir.x(); nomdy_ = ndir.y(); nomdz_ = ndir.z();
	  rnom_ = npos.rho();
	  phinom_ = npos.phi();
	  dphinom_ = ndir.phi();
	  deltax_ = delta.x(); deltay_ = delta.y(); deltaz_ = delta.z();
	  trkalign_->Fill();
	}
      }
    }
    first_ = false;
  }
}

using mu2e::TrkGeomTest;
DEFINE_ART_MODULE(TrkGeomTest);
