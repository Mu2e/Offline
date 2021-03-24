// Pass events with StatusG4 <= maxAcceptedStatus.
// See MCDataProducts/inc/StatusG4.hh for the meaning of different status values.
// By default maxAcceptedStatus=0 and only perfect events are passed.
//
// Andrei Gaponenko, 2013

#include <string>
#include <map>
#include <sstream>

// art includes.
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "MCDataProducts/inc/StatusG4.hh"

namespace mu2e {

  //================================================================
  class FilterStatusG4 : public art::EDFilter {
    art::ProductToken<StatusG4> const token_;
    int maxAcceptedStatus_;
    typedef std::map<int,int> StatMap;
    StatMap stats_;
  public:

    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Atom<art::InputTag> input {
        Name("input"),
          Comment("The StatusG4 object to use.")
          };

      fhicl::Atom<int> maxAcceptedStatus {
        Name("maxAcceptedStatus"),
          Comment("The filter passes events with StatusG4 <= maxAcceptedStatus.\n"
                  "See MCDataProducts/inc/StatusG4.hh for the meaning of different status values.\n"
                  "By default maxAcceptedStatus=0 and only perfect events are passed."
                  ),
          0
          };
    };

    using Parameters = art::EDFilter::Table<Config>;
    explicit FilterStatusG4(const Parameters& conf);

    virtual bool filter(art::Event& event) override;
    virtual void endJob() override;
  };

  //================================================================
  FilterStatusG4::FilterStatusG4(const Parameters& conf)
    : art::EDFilter{conf}
    , token_{consumes<StatusG4>(conf().input())}
    , maxAcceptedStatus_(conf().maxAcceptedStatus())
  {}

  //================================================================
  bool FilterStatusG4::filter(art::Event& event) {
    auto ih = event.getValidHandle(token_);
    ++stats_[ih->status()];
    return (ih->status() <= maxAcceptedStatus_);
  }

  //================================================================
  void FilterStatusG4::endJob() {
    std::ostringstream os;
    os<< "FilterStatusG4_module summary:\n";
    for(const auto &i : stats_) {
      os<<"    "<< i.second<< " events with status "<<i.first<<"\n";
    }

    mf::LogInfo("Summary")<<os.str();
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::FilterStatusG4);
