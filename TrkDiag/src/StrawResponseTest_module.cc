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
#include "TCanvas.h"
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
        fhicl::Atom<double> tmin{ Name("tmin"), Comment("Minimum time to test" ),-5.0 };
        fhicl::Atom<double> tmax{ Name("tmax"), Comment("Maximum time to test" ),45.0 };
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit StrawResponseTest(Parameters const& config);
      virtual ~StrawResponseTest();
      void analyze( art::Event const & e) override;
      void beginJob ( ) override;
      void endJob ( ) override;
    private:
      ProditionsHandle<Tracker> _alignedTracker_h;
      TTree*  srtest_;
      TCanvas* srtcan_;
      float dtime_, vinst_;
      float rderr_, rdrift_, cdrift_, rnerr_;
      TGraph* t2d_, *t2derr_, *t2nerr_, *t2v_;
      int plane_, panel_, straw_;
      int print_, diag_;
      int nbins_;
      double tmin_, tmax_;
      bool first_;
      ProditionsHandle<StrawResponse> strawResponse_h_;
  };

  StrawResponseTest::StrawResponseTest(Parameters const& config) :     art::EDAnalyzer(config),
  srtest_(0),t2d_(0), t2derr_(0), t2nerr_(0), t2v_(0),
  print_(config().printLevel()),
  diag_(config().diagLevel()),
  nbins_(config().nbins()),
  tmin_(config().tmin()),
  tmax_(config().tmax()),
  first_(true) {}

  StrawResponseTest::~StrawResponseTest() {}

  void StrawResponseTest::beginJob(){
    if(diag_ > 0){
      art::ServiceHandle<art::TFileService> tfs;
      srtest_=tfs->make<TTree>("SRT","SRT");
//      srtest_->Branch("plane",&splane_,"plane/I");  Eventually add straw dependence TODO
//      srtest_->Branch("panel",&spanel_,"panel/I");
//      srtest_->Branch("straw",&straw_,"straw/I");
      srtest_->Branch("rdrift",&rdrift_,"rdrift/F");
      srtest_->Branch("cdrift",&cdrift_,"cdrift/F");
      srtest_->Branch("dtime",&dtime_,"dtime/F");
      srtest_->Branch("vinst",&vinst_,"vinst/F");
      srtest_->Branch("rderr",&rderr_,"rderr/F");
      srtest_->Branch("rnerr",&rnerr_,"rnerr/F");
      t2d_ = tfs->make<TGraph>(nbins_);
      t2d_->SetName("T2D");
      t2d_->SetTitle("Time to Distance;Drift Time (ns);Distance to Wire (mm)");
      gDirectory->Append(t2d_);
      t2derr_ = tfs->make<TGraph>(nbins_);
      t2derr_->SetName("T2SDErr");
      t2derr_->SetTitle("Signed Drift Error;Drift Time (ns);Drift Error (mm)");
      gDirectory->Append(t2derr_);
      t2nerr_ = tfs->make<TGraph>(nbins_);
      t2nerr_->SetName("T2UDErr");
      t2nerr_->SetTitle("Unsigned Drift Error;Drift Time (ns);Drift Error (mm)");
      gDirectory->Append(t2nerr_);
      t2v_ = tfs->make<TGraph>(nbins_);
      t2v_->SetName("T2V");
      t2v_->SetTitle("Instantaneous Drift Velocity;Drift Time (ns);Instantaneous Drift Velocity (ns/mm)");
      gDirectory->Append(t2v_);
      srtcan_=tfs->make<TCanvas>("SRTcan","SRTcan",800,800);
      gDirectory->Append(srtcan_);
    }
  }

  void StrawResponseTest::analyze(art::Event const& event) {
    // this is a test module, so only a single event is processed.
    StrawId sid;

    if(first_){
      // fetch StrawResposne
      auto const& sresponse = strawResponse_h_.getPtr(event.id());
      double tstep = (tmax_ - tmin_)/(nbins_-1);
      for (int ibin=0;ibin<nbins_; ++ibin){
        dtime_ = tmin_ + tstep*ibin;
        DriftInfo dinfo = sresponse->driftInfo(sid,dtime_,0.0);
        rdrift_ = dinfo.rDrift_;
        cdrift_ = dinfo.cDrift_;
        rderr_ = dinfo.signedDriftError_;
        rnerr_ = dinfo.unsignedDriftError_;
        vinst_ = dinfo.driftVelocity_;
        srtest_->Fill();
        t2d_->SetPoint(ibin,dtime_,rdrift_);
        t2derr_->SetPoint(ibin,dtime_,rderr_);
        t2nerr_->SetPoint(ibin,dtime_,rnerr_);
        t2v_->SetPoint(ibin,dtime_,vinst_);
      }
    }
  }
  void StrawResponseTest::endJob(){
    if(diag_ > 0){
       srtcan_->Divide(2,2);
       srtcan_->cd(1);
       t2d_->Draw("APL");
       srtcan_->cd(2);
       t2derr_->Draw("APL");
       srtcan_->cd(3);
       t2v_->Draw("APL");
   }
  }
}

using mu2e::StrawResponseTest;
DEFINE_ART_MODULE(StrawResponseTest)
