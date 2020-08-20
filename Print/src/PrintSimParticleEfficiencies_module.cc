//
//  A module to print SimParticleEfficiencies
//  -- we want to be able to print to a file and not print extraneous information
//

#include <vector>
#include <ostream>
#include <iomanip>
#include <memory>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "Print/inc/ProductPrinter.hh"
#include "Print/inc/SimStageEfficiencyPrinter.hh"

using namespace std;

namespace mu2e {

  class PrintSimParticleEfficiencies : public art::EDAnalyzer {

  public:

    struct Config {
      fhicl::Table<SimStageEfficiencyPrinter::Config> simStageEfficiencyPrinter { 
	fhicl::Name("simStageEfficiencyPrinter") }; 
      fhicl::Atom<std::string> outFileName {
	fhicl::Name("outFileName")};
    };

    // this line is required by art to allow the command line help print
    typedef art::EDAnalyzer::Table<Config> Parameters;

    explicit PrintSimParticleEfficiencies(const Parameters& conf);
    void analyze  ( art::Event const&  event  ) override;
    void beginSubRun( art::SubRun const& subrun) override;
    void beginRun( art::Run const& run) override;
    void endRun( art::Run const& run) override;

  private:

    // set by fcl, turn all prints on (1) or off (0)
    int _verbose;
    // each of these object prints a different product
    vector< unique_ptr<mu2e::ProductPrinter> > _printers;

    std::string _outfilename;
    std::ofstream _outfile;
  };

}


mu2e::PrintSimParticleEfficiencies::PrintSimParticleEfficiencies(const Parameters& conf):
  art::EDAnalyzer(conf),
  _outfilename(conf().outFileName()),
  _outfile(_outfilename.c_str())
{
  //cout << "start main pset\n"<< pset.to_string() << "\n end main pset"<< endl;
  _outfile << "TABLE SimEfficiencies" << std::endl;

  _printers.push_back( make_unique<SimStageEfficiencyPrinter>( conf().simStageEfficiencyPrinter() ) );
}


void mu2e::PrintSimParticleEfficiencies::analyze(art::Event const& event) {
  for(auto& prod_printer: _printers) prod_printer->Print(event, _outfile);
}

void mu2e::PrintSimParticleEfficiencies::beginSubRun(art::SubRun const& subrun) {
  for(auto& prod_printer: _printers) prod_printer->PrintSubRun(subrun, _outfile);
}

void mu2e::PrintSimParticleEfficiencies::beginRun(art::Run const& run) {
  for(auto& prod_printer: _printers) prod_printer->PrintRun(run, _outfile);
}

void mu2e::PrintSimParticleEfficiencies::endRun(art::Run const& run) {
}


DEFINE_ART_MODULE(mu2e::PrintSimParticleEfficiencies)
