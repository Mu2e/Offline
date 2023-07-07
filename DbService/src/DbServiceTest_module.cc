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
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
//#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

#include "Offline/DbService/inc/DbHandle.hh"
#include "Offline/DbService/inc/DbService.hh"
#include "Offline/DbTables/inc/TstCalib1.hh"
#include "Offline/DbTables/inc/TstCalib2.hh"

namespace mu2e {

class DbServiceTest : public art::EDAnalyzer {
 public:
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Atom<int> verbose{Name("verbose"), Comment("verbose flag, 0 to 10"),
                             1};

    fhicl::Sequence<std::string> tableList{
        fhicl::Name("tableList"), fhicl::Comment("Sequence of table names"),
        std::vector<std::string>()};
  };

  // this line is required by art to allow the command line help print
  typedef art::EDAnalyzer::Table<Config> Parameters;

  explicit DbServiceTest(const Parameters& conf) :
      art::EDAnalyzer(conf), _conf(conf()) {}

  ~DbServiceTest() {}
  void beginJob() override;
  void beginSubRun(const art::SubRun& subrun) override;
  void analyze(const art::Event& event) override;

 private:
  Config _conf;

  mu2e::DbHandle<mu2e::TstCalib1> _testCalib1_h;
  mu2e::DbHandle<mu2e::TstCalib2> _testCalib2_h;

  art::ServiceHandle<DbService> _dbh;
};

//-----------------------------------------------------------------------------

void DbServiceTest::beginJob() {}

void DbServiceTest::beginSubRun(const art::SubRun& subrun) {}

//-----------------------------------------------------------------------------
void DbServiceTest::analyze(const art::Event& event) {
  std::cout << "DbServiceTest::analyze  " << event.id() << std::endl;

  for (auto const& tname : _conf.tableList()) {
    if (tname == "TstCalib1") {
      auto const& tc1 = _testCalib1_h.get(event.id());
      std::cout << tc1.name() << " got cid=" << _testCalib1_h.cid()
                << "  rows=" << tc1.nrow() << "\n";
      std::cout << tc1.csv();
    } else if (tname == "TstCalib2") {
      auto const& tc2 = _testCalib2_h.get(event.id());
      std::cout << tc2.name() << " got cid=" << _testCalib2_h.cid()
                << "  rows=" << tc2.nrow() << "\n";
      std::cout << tc2.csv();
    } else {
      art::EventID eid = event.id();
      int tid = _dbh->engine().tidByName(tname);
      auto const& liveTable = _dbh->engine().update(tid, uint32_t(eid.run()),
                                                    uint32_t(eid.subRun()));
      auto const& tt = liveTable.table();
      std::cout << tt.name() << " got cid=" << liveTable.cid()
                << "  rows=" << tt.nrow() << "\n";
    }
  }
}
}  // namespace mu2e

DEFINE_ART_MODULE(mu2e::DbServiceTest)
