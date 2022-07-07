//
// This module analyzes StrawResposne
//
// Original author David Brown, 11/16/2020 LBNL
//
// framework
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "TDirectory.h"
#include "TTree.h"
#include "TGraph.h"
#include <iostream>

namespace mu2e {

  class StrawResponseTest : public art::EDAnalyzer
  {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      struct Config {
        fhicl::Atom<int> printLevel{ Name("printLevel"), Comment("Print level" ),0 };
        fhicl::Atom<int> diagLevel{ Name("diagLevel"), Comment("Diagnostic level" ),0 };
        fhicl::Atom<int> nbins{ Name("nBins"), Comment("Number of bins for TGraph objects" ),200 };
        fhicl::Atom<double> rmax{ Name("rmax"), Comment("Maximum radius to test" ),3.0 };
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit StrawResponseTest(Parameters const& config);
      virtual ~StrawResponseTest();
      void analyze( art::Event const & e) override;
      void beginJob ( ) override;
    private:
      ProditionsHandle<Tracker> _alignedTracker_h;
      TTree*  srtest_;
      float dist_, time_, tvar_, vinst_, vavg_;
      TGraph* d2t_, *d2tvar_, *d2v_;
      int plane_, panel_, straw_;
      int print_, diag_;
      int nbins_;
      double rmax_;
      bool first_;
      ProditionsHandle<StrawResponse> strawResponse_h_;
  };

  StrawResponseTest::StrawResponseTest(Parameters const& config) :     art::EDAnalyzer(config),
  srtest_(0),d2t_(0), d2tvar_(0), d2v_(0),
  print_(config().printLevel()),
  diag_(config().diagLevel()),
  nbins_(config().nbins()),
  rmax_(config().rmax()),
  first_(true) {}

  StrawResponseTest::~StrawResponseTest() {}

  void StrawResponseTest::beginJob(){
    if(diag_ > 0){
      art::ServiceHandle<art::TFileService> tfs;
      srtest_=tfs->make<TTree>("strawtest","strawtest");
//      srtest_->Branch("plane",&splane_,"plane/I");  Eventually add straw dependence TODO
//      srtest_->Branch("panel",&spanel_,"panel/I");
//      srtest_->Branch("straw",&straw_,"straw/I");
      srtest_->Branch("dist",&dist_,"dist/F");
//      srtest_->Branch("phi",&phi_,"phi/F"); // add phi dependence TODO
      srtest_->Branch("time",&time_,"time/F");
      srtest_->Branch("tvar",&tvar_,"tvar/F");
      srtest_->Branch("vinst",&vinst_,"vinst/F");
      srtest_->Branch("vavg",&vavg_,"vavg/F");
      d2t_ = tfs->make<TGraph>(nbins_);
      d2t_->SetName("D2T");
      d2t_->SetTitle("Distance to Time;Distance to Wire (mm);Mean Drift Time (ns)");
      gDirectory->Append(d2t_);
      d2tvar_ = tfs->make<TGraph>(nbins_);
      d2tvar_->SetName("D2TVar");
      d2tvar_->SetTitle("Drift Time Variance;Distance to Wire (mm);Drift Time Variance (ns^2)");
      gDirectory->Append(d2tvar_);
      d2v_ = tfs->make<TGraph>(nbins_);
      d2v_->SetName("D2V");
      d2v_->SetTitle("Instantaneous Drift Velocity;Distance to Wire (mm);Drift Velocity (mm/ns)");
      gDirectory->Append(d2v_);
    }
  }

  void StrawResponseTest::analyze(art::Event const& event) {
    // this is a test module, so only a single event is processed.
    StrawId sid;

    if(first_){
      // fetch StrawResposne
      auto const& sresponse = strawResponse_h_.getPtr(event.id());
      vavg_ = sresponse->driftConstantSpeed();
      double dstep = rmax_/(nbins_-1);
      for (int ibin=0;ibin<nbins_; ++ibin){
        double dist = dstep*ibin;
        auto dinfo = sresponse->driftInfoAtDistance(sid,dist,0.0);
        dist_ = dist; time_ = dinfo.time; vinst_ = dinfo.invSpeed; tvar_ = dinfo.variance;
        srtest_->Fill();
        d2t_->SetPoint(ibin,dist,dinfo.time);
        d2tvar_->SetPoint(ibin,dist,dinfo.variance);
        d2v_->SetPoint(ibin,dist,dinfo.invSpeed);
      }
    }
  }
}

using mu2e::StrawResponseTest;
DEFINE_ART_MODULE(StrawResponseTest);
