// Andrei Gaponenko, 2013

#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include "cetlib_except/exception.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Provenance.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "MCDataProducts/inc/GenEventCount.hh"

#include "TH1D.h"

namespace mu2e {

  //================================================================
  class GenEventCountReader : public art::EDAnalyzer {
    GenEventCount::count_t numEvents_;
    unsigned numSubRuns_;
    bool makeHistograms_;
  public:

    struct Config {
      fhicl::Atom<bool> makeHistograms{
        fhicl::Name("makeHistograms"),
          fhicl::Comment("Write out number of events and subruns as histograms, in addition to printing them out. "),
          true
          };
    };

    using Parameters = art::EDAnalyzer::Table<Config>;
    explicit GenEventCountReader(const Parameters& conf);

    virtual void analyze(const art::Event&) override {}
    virtual void endSubRun(const art::SubRun& sr) override;
    virtual void endJob() override;
  };

  //================================================================
  GenEventCountReader::GenEventCountReader(const Parameters& conf)
    : art::EDAnalyzer(conf)
    , numEvents_(0)
    , numSubRuns_(0)
    , makeHistograms_(conf().makeHistograms())
  {}

  //================================================================
  void GenEventCountReader::endSubRun(const art::SubRun& sr) {

    // We expect exactly one object of type GenEventCount per SubRun.
    std::vector<art::Handle<GenEventCount> > hh = sr.getMany<GenEventCount>();
    if(hh.size() > 1) {
      std::ostringstream os;
      os<<"GenEventCountReader: multiple GenEventCount objects found in "
        <<sr.id()<<":\n";
      for(const auto& h : hh) {
        os<<"    moduleLabel = "<<h.provenance()->moduleLabel()
          <<", instance = "<<h.provenance()->productInstanceName()
          <<", process = "<<h.provenance()->processName()
          <<"\n";
      }
      os<<"\n";
      throw cet::exception("BADCONFIG")<<os.str();
    }
    else if(hh.empty()) {
      throw cet::exception("BADCONFIG")
        <<"GenEventCountReader: no GenEventCount record in "<<sr.id()<<"\n";
    }

    mf::LogInfo("INFO")<<"GenEventCount: "
                       <<hh.front()->count()<<" events in "<<sr.id()
                       <<"\n";

    ++numSubRuns_;
    numEvents_ += hh.front()->count();
  }

  //================================================================
  void GenEventCountReader::endJob() {
    mf::LogInfo("Summary")<<"GenEventCount total: "
                <<numEvents_<<" events in "
                <<numSubRuns_<<" SubRuns"
                <<"\n";

    if(makeHistograms_) {
      art::ServiceHandle<art::TFileService> tfs;
      TH1* hEvents = tfs->make<TH1D>("numEvents", "numEvents", 1, -0.5, 0.5);
      hEvents->Fill(0., numEvents_);
      TH1* hSubRuns = tfs->make<TH1D>("numSubRuns", "numSubRuns", 1, -0.5, 0.5);
      hSubRuns->Fill(0., numSubRuns_);
    }
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::GenEventCountReader);
