//
//  A module to print products in an event
//

#include <vector>
#include <ostream>
#include <iomanip>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "Print/inc/ProductPrinter.hh"
#include "Print/inc/GenParticlePrinter.hh"
#include "Print/inc/SimParticlePrinter.hh"
#include "Print/inc/StepPointMCPrinter.hh"
#include "Print/inc/MCTrajectoryPrinter.hh"
#include "Print/inc/CaloHitPrinter.hh"
#include "Print/inc/CaloDigiPrinter.hh"
#include "Print/inc/CaloRecoDigiPrinter.hh"
#include "Print/inc/CaloCrystalHitPrinter.hh"
#include "Print/inc/CaloClusterPrinter.hh"
#include "Print/inc/StrawDigiPrinter.hh"
#include "Print/inc/StrawHitPrinter.hh"
#include "Print/inc/TrackClusterMatchPrinter.hh"
#include "Print/inc/TrkCaloIntersectPrinter.hh"
#include "Print/inc/TrackSummaryPrinter.hh"
#include "Print/inc/KalRepPrinter.hh"

namespace mu2e {

  class PrintModule : public art::EDAnalyzer {

  public:

    explicit PrintModule(fhicl::ParameterSet const& );
    void analyze  ( art::Event const&  event  ) override;

  private:

    // set by fcl, turn all prints on (1) or off (0)
    int _verbose;
    // each of these object prints a different product
    std::vector<mu2e::ProductPrinter*> _printers;
  };

}


mu2e::PrintModule::PrintModule(fhicl::ParameterSet const& pset ):
  art::EDAnalyzer(pset) {
  //std::cout << "start main pset\n"<< pset.to_string() << "\n end main pset"<< std::endl;

  _printers.push_back( new GenParticlePrinter(pset) );
  _printers.push_back( new SimParticlePrinter(pset) );
  _printers.push_back( new StepPointMCPrinter(pset) );
  _printers.push_back( new MCTrajectoryPrinter(pset) );
  _printers.push_back( new CaloHitPrinter(pset) );
  _printers.push_back( new CaloDigiPrinter(pset) );
  _printers.push_back( new CaloRecoDigiPrinter(pset) );
  _printers.push_back( new CaloCrystalHitPrinter(pset) );
  _printers.push_back( new CaloClusterPrinter(pset) );
  _printers.push_back( new StrawDigiPrinter(pset) );
  _printers.push_back( new StrawHitPrinter(pset) );
  _printers.push_back( new TrackClusterMatchPrinter(pset) );
  _printers.push_back( new TrkCaloIntersectPrinter(pset) );
  _printers.push_back( new TrackSummaryPrinter(pset) );
  _printers.push_back( new KalRepPrinter(pset) );

}


void
mu2e::PrintModule::analyze(art::Event const& event) {
  std::cout 
    << "\n"
    << " ###############  PrintModule Run/Subrun/Event " 
    << std::setw(9) << event.run()
    << std::setw(9) << event.subRun()
    << std::setw(9) << event.event()
    << std::endl;

  for(auto& prod_printer: _printers) prod_printer->Print(event);

  std::cout << std::endl;

}


DEFINE_ART_MODULE(mu2e::PrintModule)
