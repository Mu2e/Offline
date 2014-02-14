//
// Write a text file that can be used as input for the Source00 module.
//
//   Original author Rob Kutschke
//
// The format of the file is one event per line:
//    runNumber subRunNumber eventNumber dataProduct
//
// where dataProduct is a trivial data product: it is a single
// integer whose value is:
//
//    10000*runNumber + 100*subRunNumber + eventNumber
//

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"

#include <iostream>
#include <fstream>

using namespace std;

namespace mu2e {

  class MakeTestInputFile : public art::EDAnalyzer {
  public:

    explicit MakeTestInputFile(fhicl::ParameterSet const& pset);

    void analyze    ( art::Event  const&  event) override;

  private:

    ofstream ofile_;

  };

  MakeTestInputFile::MakeTestInputFile(fhicl::ParameterSet const& pset):
    EDAnalyzer(pset),
    ofile_( pset.get<std::string>("filename")){
  }

  void MakeTestInputFile::analyze(art::Event const& event) {

    double dataProduct = 10000*event.id().run() + 100*event.id().subRun() + event.id().event();
    ofile_ << event.id().run()    <<  " "
           << event.id().subRun() <<  " "
           << event.id().event()  <<  " "
           << dataProduct << endl;
  }

}  // end namespace mu2e

using mu2e::MakeTestInputFile;
DEFINE_ART_MODULE(MakeTestInputFile)
