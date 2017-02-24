// A deliberate memory leak to use to test the art --memcheck option.
//
// Rob Kutschke, 2015

#include <iostream>
#include <vector>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"

#include "cetlib_except/exception.h"

namespace mu2e {

  class Leak : public art::EDAnalyzer {
  public:
    explicit Leak(fhicl::ParameterSet const& pset);
    void analyze(const art::Event& evt) override;
  private:
    size_t everyN_       = 1;  // Leak every N events; when event counter mod everyN_=0.
    size_t nwords_       = 1;  // Number of 4 byte words to leak
    size_t sleepSeconds_ = 0;  // Number of seconds to sleep after leaking
                               //  - Use this to watch with top or another monitoring tool.
    size_t eventCounter_ = 0;  // Count number of events we have seen so far.

  };

  Leak::Leak(const fhicl::ParameterSet& pset)
    : art::EDAnalyzer(pset)
    , everyN_(pset.get<size_t>("everyN"))
    , nwords_(pset.get<size_t>("nwords"))
    , sleepSeconds_(pset.get<size_t>("sleepSeconds"))
  {
    std::cout << "Leak: will leak "
              << nwords_ << " words (of 4 bytes) every "
              << everyN_ << " events."
              << std::endl;
    if ( everyN_ == 0 ){
      throw cet::exception("CONFIG") << "Value of everyN must be positive;  it's value is: " << everyN_ << "\n";
    }
  }

  void Leak::analyze(const art::Event& event) {
    if ( (++eventCounter_)%everyN_ == 0 ){
      new std::vector<int>(nwords_,0);
    }
    if ( sleepSeconds_ > 0 ){
      sleep(sleepSeconds_);
    }
  }

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::Leak);
