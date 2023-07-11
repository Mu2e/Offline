
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
#include "Offline/RecoDataProducts/inc/ProtonBunchTime.hh"
#include "Offline/MCDataProducts/inc/ProtonBunchTimeMC.hh"
#include "Offline/DataProducts/inc/EventWindowMarker.hh"
#include "TTree.h"

namespace mu2e {

  class ProtonBunchTimeDiag : public art::EDAnalyzer
  {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      struct Config {
        fhicl::Atom<art::InputTag> ewMarkerTag{ Name("EventWindowMarkerTag"), Comment("EventWindowMarker producer"),"EWMProducer" };
        fhicl::Atom<art::InputTag> pbtmcTag{ Name("ProtonBunchTimeMCTag"), Comment("ProtonBunchTimeMC producer"),"EWMProducer" };
        fhicl::Atom<art::InputTag> pbtimeTag{ Name("ProtonBunchTimeTag"), Comment("ProtonBunchTime producer"),"PBTFSD" };
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit ProtonBunchTimeDiag(Parameters const& config);
      virtual ~ProtonBunchTimeDiag();
      void analyze( art::Event const & e) override;
      void beginJob ( ) override;
    private:
      art::InputTag ewmtag_, pbtmctag_, pbttag_;
      TTree* pbttest_;
      int iev_, nhits_, ewlength_;
      float tmean_, terr_, pbtimemc_;
  };

  ProtonBunchTimeDiag::ProtonBunchTimeDiag(Parameters const& config) :
    art::EDAnalyzer(config),
    ewmtag_(config().ewMarkerTag()),
    pbtmctag_(config().pbtmcTag()),
    pbttag_(config().pbtimeTag()),
    pbttest_(0)
  {
    consumes<EventWindowMarker>(ewmtag_);
    consumes<ProtonBunchTimeMC>(pbtmctag_);
    consumes<ProtonBunchTime>(pbttag_);
  }

  ProtonBunchTimeDiag::~ProtonBunchTimeDiag() {}

  void ProtonBunchTimeDiag::beginJob(){
    art::ServiceHandle<art::TFileService> tfs;
    pbttest_=tfs->make<TTree>("pbtdiag","proton bunch time diagnostics");
    pbttest_->Branch("iev",&iev_,"iev/I");
    pbttest_->Branch("tmean",&tmean_,"tmean/F");
    pbttest_->Branch("terr",&terr_,"terr/F");
    pbttest_->Branch("pbtimemc",&pbtimemc_,"pbtimemc/F");
    pbttest_->Branch("ewlength",&ewlength_,"ewlength/I");
  }

  void ProtonBunchTimeDiag::analyze(art::Event const& event) {
    auto pbtH = event.getValidHandle<ProtonBunchTime>(pbttag_);
    auto const& pbt(*pbtH);
    auto pbtmcH = event.getValidHandle<ProtonBunchTimeMC>(pbtmctag_);
    auto const& pbtmc(*pbtmcH);
    auto ewmH = event.getValidHandle<EventWindowMarker>(ewmtag_);
    auto const& ewm(*ewmH);

    iev_ = event.id().event();
    tmean_= pbt.pbtime_;
    terr_ = pbt.pbterr_;
    pbtimemc_ = pbtmc.pbtime_;
    ewlength_ = ewm.eventLength();
    pbttest_->Fill();
  }
}

using mu2e::ProtonBunchTimeDiag;
DEFINE_ART_MODULE(ProtonBunchTimeDiag)

