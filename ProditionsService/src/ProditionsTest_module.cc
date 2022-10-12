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

#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/STMConditions/inc/STMEnergyCalib.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"

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

  // ProditionsHandle<StrawResponse> _testh;
  ProditionsHandle<STMEnergyCalib> _testh;

};

//-----------------------------------------------------------------------------
void ProditionsTest::analyze(const art::Event& event) {
  std::cout << "ProditionsTest::analyze  " << event.id() << std::endl;

  auto const& stm = _testh.get(event.id());
  auto cl = STMChannel(STMChannel::enum_type::LaBr);
  auto ch = STMChannel(STMChannel::enum_type::HPGe);

  auto r = stm.calib(cl);
  std::cout << cl.name() << " " << r.p0 << " " << r.p1 << " " << r.p2 << "\n";
  r = stm.calib(ch);
  std::cout << ch.name() << " " << r.p0 << " " << r.p1 << " " << r.p2 << "\n";

  float pl = stm.pedestal(cl);
  float ph = stm.pedestal(ch);
  std::cout << "pedestals l,h "<< pl << " " <<ph << "\n";

}

}  // namespace mu2e

DEFINE_ART_MODULE(mu2e::ProditionsTest);
