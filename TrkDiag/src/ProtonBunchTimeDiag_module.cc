
//
// This module analyzes the proton bunch time reconstruction, comparing
// it to the MC true EventWindow offset
//
// Original author David Brown, 10/1/2020 LBNL
//
// framework
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "RecoDataProducts/inc/ProtonBunchTime.hh"
#include "DataProducts/inc/EventWindowMarker.hh"
#include "TTree.h"

namespace mu2e {

  class ProtonBunchTimeDiag : public art::EDAnalyzer
  {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      struct Config {
	fhicl::Atom<unsigned> minnhits{ Name("minNHits"), Comment("Minimum number of hits needed to use result"), 1};
	fhicl::Atom<art::InputTag> ewMarkerTag{ Name("EventWindowMarkerTag"), Comment("EventWindowMarker producer"),"EWMProducer" };
	fhicl::Atom<art::InputTag> pbtimeTag{ Name("ProtonBunchTimeTag"), Comment("ProtonBunchTime producer"),"PBTFSD" };
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit ProtonBunchTimeDiag(Parameters const& config);
      virtual ~ProtonBunchTimeDiag();
      void analyze( art::Event const & e) override;
      void beginJob ( ) override;
    private:
      art::InputTag ewmtag_, pbttag_;
      TTree* pbttest_;
      int iev_, nhits_;
      float tmean_, tvariance_, ewmoffset_;
  };

  ProtonBunchTimeDiag::ProtonBunchTimeDiag(Parameters const& config) :
    art::EDAnalyzer(config),
    ewmtag_(config().ewMarkerTag()),
    pbttag_(config().pbtimeTag()),
    pbttest_(0)
  {
    consumes<EventWindowMarker>(ewmtag_);
    consumes<ProtonBunchTime>(pbttag_);
  }

  ProtonBunchTimeDiag::~ProtonBunchTimeDiag() {}

  void ProtonBunchTimeDiag::beginJob(){ 
    art::ServiceHandle<art::TFileService> tfs;
    pbttest_=tfs->make<TTree>("pbtdiag","proton bunch time diagnostics");
    pbttest_->Branch("iev",&iev_,"iev/I");
    pbttest_->Branch("nhits",&nhits_,"nhits/I");
    pbttest_->Branch("tmean",&tmean_,"tmean/F");
    pbttest_->Branch("tvariance",&tvariance_,"tvariance/F");
    pbttest_->Branch("ewmoffset",&ewmoffset_,"ewmoffset/F");
  }

  void ProtonBunchTimeDiag::analyze(art::Event const& event) {        
    auto pbtH = event.getValidHandle<ProtonBunchTime>(pbttag_);
    auto const& pbt(*pbtH);
    auto ewmH = event.getValidHandle<EventWindowMarker>(ewmtag_);
    auto const& ewm(*ewmH);

    iev_ = event.id().event();
    nhits_ = pbt.nhits_; 
    tmean_= pbt.pbtime_; 
    tvariance_ = pbt.variance_; 
    ewmoffset_ = ewm._timeOffset;
    pbttest_->Fill();
  }
}

using mu2e::ProtonBunchTimeDiag;
DEFINE_ART_MODULE(ProtonBunchTimeDiag);

