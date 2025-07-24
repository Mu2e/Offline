//
// This module analyzes the tracker alignment
//
// Original author David Brown, 11/16/2020 LBNL
//
// framework
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/TrackerConditions/inc/TrackerStatus.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
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
      TTree* strawtest_, *paneltest_;
      int print_, diag_;
      int splane_, spanel_, straw_;
      int pplane_, panel_;
      float nomx_, nomy_, nomz_;
      float nomdx_, nomdy_, nomdz_;
      float rnom_, phinom_, dphinom_;
      float deltax_, deltay_, deltaz_;
      bool first_;
      float uphi_, vphi_, wcost_;
      float oz_, or_, ophi_;
  };

  TrkGeomTest::TrkGeomTest(Parameters const& config) :     art::EDAnalyzer(config),
  strawtest_(0),
  print_(config().printLevel()),
  diag_(config().diagLevel()),
  first_(true) {}

  TrkGeomTest::~TrkGeomTest() {}

  void TrkGeomTest::beginJob(){
    if(diag_ > 0){
      art::ServiceHandle<art::TFileService> tfs;
      strawtest_=tfs->make<TTree>("strawtest","strawtest");
      strawtest_->Branch("plane",&splane_,"plane/I");
      strawtest_->Branch("panel",&spanel_,"panel/I");
      strawtest_->Branch("straw",&straw_,"straw/I");
      strawtest_->Branch("nomx",&nomx_,"nomx/F");
      strawtest_->Branch("nomy",&nomy_,"nomy/F");
      strawtest_->Branch("nomz",&nomz_,"nomz/F");
      strawtest_->Branch("nomdx",&nomdx_,"nomdx/F");
      strawtest_->Branch("nomdy",&nomdy_,"nomdy/F");
      strawtest_->Branch("nomdz",&nomdz_,"nomdz/F");
      strawtest_->Branch("rnom",&rnom_,"rnom/F");
      strawtest_->Branch("phinom",&phinom_,"phinom/F");
      strawtest_->Branch("dphinom",&dphinom_,"dphinom/F");
      strawtest_->Branch("deltax",&deltax_,"deltax/F");
      strawtest_->Branch("deltay",&deltay_,"deltay/F");
      strawtest_->Branch("deltaz",&deltaz_,"deltaz/F");
      paneltest_=tfs->make<TTree>("paneltest","paneltest");
      paneltest_->Branch("plane",&pplane_,"plane/I");
      paneltest_->Branch("panel",&panel_,"panel/I");
      paneltest_->Branch("uphi",&uphi_,"uphi/F");
      paneltest_->Branch("vphi",&vphi_,"vphi/F");
      paneltest_->Branch("wcost",&wcost_,"wcost/F");
      paneltest_->Branch("oz",&oz_,"oz/F");
      paneltest_->Branch("or",&or_,"or/F");
      paneltest_->Branch("ophi",&ophi_,"ophi/F");
    }
  }

  void TrkGeomTest::analyze(art::Event const& event) {
    // this is a test module, so only a single event is processed.
    if(first_){
      // fetch tracker status
      ProditionsHandle<TrackerStatus> trackerstatus_h;
      auto const& trkstatus = trackerstatus_h.get(event.id());
      if(print_ > 1){
        std::cout << "TrackerStatus " << std::endl;
        trkstatus.print(std::cout);
      }
      // fetch aligned and nominal tracker objects
      GeomHandle<Tracker> nominalTracker_h;
      auto const& ntracker = *nominalTracker_h;
      auto const& atracker = _alignedTracker_h.get(event.id());
      double totallength(0.0);
      if(ntracker.straws().size() != atracker.straws().size()){
        std::cout << "Trackers don't match" << std::endl;
      } else {
        std::set<double> zpos;
        for(size_t istr = 0; istr < ntracker.straws().size(); istr++){
          auto const& nstraw = ntracker.straws()[istr];
          auto const& astraw = atracker.straws()[istr];
          if(nstraw.id() != astraw.id())
            std::cout << "StrawIds don't match: nominal " << nstraw.id() << " aligned " << astraw.id() << std::endl;
          totallength += 2.0*nstraw.halfLength();
          splane_ = nstraw.id().plane();
          spanel_ = nstraw.id().panel();
          straw_ = nstraw.id().straw();
          auto npos = nstraw.origin();
          auto ndir = nstraw.wireDirection();
          auto adir = astraw.wireDirection();
          // test
          if(ndir.dot(adir)<0.0)std::cout << "Straw directions don't match: nominal " << ndir << " aligned " << adir << std::endl;
          auto apos = astraw.origin();
          auto delta = apos-npos;
          nomx_ = npos.x(); nomy_ = npos.y(); nomz_ = npos.z();
          nomdx_ = ndir.x(); nomdy_ = ndir.y(); nomdz_ = ndir.z();
          rnom_ = npos.rho();
          phinom_ = npos.phi();
          dphinom_ = ndir.phi();
          deltax_ = delta.x(); deltay_ = delta.y(); deltaz_ = delta.z();
          strawtest_->Fill();
          // compare Z positions for straws 1 and 0 in panels 0 and 1
          if(print_>0){
            if(splane_ == 9){
              if(spanel_ == 0 || spanel_ == 1){
                if(straw_ == 0 || straw_ == 1){
                  std::cout << std::setw(8) << "Straw " << straw_ << " Panel " << spanel_ << " Z " << npos.z() << std::endl;
                  zpos.insert(npos.z());
                }
              }
            }
          }
        }
        //
        // panel test
        //
        for(auto const& panel : ntracker.panels()){
          pplane_ = panel.id().plane();
          panel_ = panel.id().panel();
          uphi_ = panel.uDirection().phi();
          vphi_ = panel.vDirection().phi();
          wcost_ = cos(panel.wDirection().theta());
          oz_ = panel.origin().z();
          or_ = panel.origin().rho();
          ophi_ = panel.origin().phi();
          paneltest_->Fill();
        }
        if(print_ > 0){
          double zavg = 0.0;
          for(auto iz = zpos.begin(); iz != zpos.end(); ++iz) zavg += *iz;
          zavg /= 4.0;
          for(auto iz = zpos.begin(); iz != zpos.end(); ++iz) std::cout << std::setw(8) << "DZ =  " << *iz - zavg << std::endl;
        }

        first_ = false;
        std::cout << "Total # Straws " << ntracker.straws().size() << " Sum length = " << totallength
          << " volume = " << totallength*M_PI*ntracker.strawProperties().strawInnerRadius() << std::endl;
      }
    }
  }
}

using mu2e::TrkGeomTest;
DEFINE_ART_MODULE(TrkGeomTest)
