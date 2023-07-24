///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////

// C++ includes
#include <iostream>
#include <stdexcept>
#include <string>

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

#include "Offline/CRVConditions/inc/CRVCalib.hh"
#include "Offline/CRVConditions/inc/CRVOrdinal.hh"
#include "Offline/CRVConditions/inc/CRVStatus.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"

namespace mu2e {

class ProditionsTest : public art::EDAnalyzer {
 public:
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Atom<int> verbose{Name("verbose"), Comment("verbose flag, 0 to 10"),
                             1};
  };

  // this line is required by art to allow the command line help print
  typedef art::EDAnalyzer::Table<Config> Parameters;

  explicit ProditionsTest(const Parameters& conf) :
      art::EDAnalyzer(conf), _conf(conf()) {}

  ~ProditionsTest() {}
  void analyze(const art::Event& event) override;

 private:
  Config _conf;

  ProditionsHandle<CRVOrdinal> _ordinal_h;
  ProditionsHandle<CRVStatus> _status_h;
  ProditionsHandle<CRVCalib> _calib_h;
};

//-----------------------------------------------------------------------------
void ProditionsTest::analyze(const art::Event& event) {
  std::cout << "ProditionsTest::analyze  " << event.id() << std::endl;

  auto const& status = _status_h.get(event.id());
  auto const& ordinal = _ordinal_h.get(event.id());
  auto const& calib = _calib_h.get(event.id());

  std::cout << "Test Ordinal\n";
  CRVROC rr = ordinal.online(1010);
  std::cout << " channel " << 1010 << " online: " << rr.ROC() << " " << rr.FEB()
            << " " << rr.FEBchannel() << "\n";
  std::cout << "online 2 2 2   offline " << ordinal.offline(CRVROC(2, 2, 2))
            << "\n";

  std::cout << "Test Status\n";
  for (auto const& ss : status.map()) {
    std::cout << ss.first << " " << ss.second << "\n";
  }

  std::cout << "Test Calib\n";
  std::size_t channel;
  channel = 11;
  std::cout << "channel " << channel << " ped: " << calib.pedestal(channel)
            << "   height " << calib.pulseHeight(channel)
            << "   area "   << calib.pulseArea(channel)
            << "   timeOffset " << calib.timeOffset(channel) << "\n";
  channel = 20111;
  std::cout << "channel " << channel << " ped: " << calib.pedestal(channel)
            << "   height " << calib.pulseHeight(channel)
            << "   area " << calib.pulseArea(channel)
            << "   timeOffset " << calib.timeOffset(channel) << "\n";
}

}  // namespace mu2e

DEFINE_ART_MODULE(mu2e::ProditionsTest)
