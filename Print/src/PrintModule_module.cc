//
//  A module to print products in an event
//

#include <vector>
#include <ostream>
#include <iomanip>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "Print/inc/ProductPrinter.hh"
#include "Print/inc/StatusG4Printer.hh"
#include "Print/inc/GenParticlePrinter.hh"
#include "Print/inc/SimParticlePrinter.hh"
#include "Print/inc/SimParticlePtrPrinter.hh"
#include "Print/inc/StepPointMCPrinter.hh"
#include "Print/inc/MCTrajectoryPrinter.hh"
#include "Print/inc/CaloShowerStepPrinter.hh"
#include "Print/inc/CaloHitPrinter.hh"
#include "Print/inc/CaloDigiPrinter.hh"
#include "Print/inc/CaloRecoDigiPrinter.hh"
#include "Print/inc/CaloCrystalHitPrinter.hh"
#include "Print/inc/CaloClusterPrinter.hh"
#include "Print/inc/StrawDigiPrinter.hh"
#include "Print/inc/StrawDigiMCPrinter.hh"
#include "Print/inc/StrawHitPrinter.hh"
#include "Print/inc/StrawHitFlagPrinter.hh"
#include "Print/inc/BkgClusterPrinter.hh"
#include "Print/inc/BkgQualPrinter.hh"
#include "Print/inc/TrackClusterMatchPrinter.hh"
#include "Print/inc/TrkCaloIntersectPrinter.hh"
#include "Print/inc/TrackSummaryPrinter.hh"
#include "Print/inc/KalRepPrinter.hh"
#include "Print/inc/SimParticleTimeMapPrinter.hh"
#include "Print/inc/ComboHitPrinter.hh"
#include "Print/inc/TimeClusterPrinter.hh"
#include "Print/inc/KalSeedPrinter.hh"
#include "Print/inc/PhysicalVolumePrinter.hh"

namespace mu2e {

  class PrintModule : public art::EDAnalyzer {

  public:

    explicit PrintModule(fhicl::ParameterSet const& );
    void analyze  ( art::Event const&  event  ) override;
    void beginSubRun( art::SubRun const& subrun) override;

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

  _printers.push_back( new StatusG4Printer(pset) );
  _printers.push_back( new GenParticlePrinter(pset) );
  _printers.push_back( new SimParticlePrinter(pset) );
  _printers.push_back( new SimParticlePtrPrinter(pset) );
  _printers.push_back( new StepPointMCPrinter(pset) );
  _printers.push_back( new MCTrajectoryPrinter(pset) );
  _printers.push_back( new CaloShowerStepPrinter(pset) );
  _printers.push_back( new CaloHitPrinter(pset) );
  _printers.push_back( new CaloDigiPrinter(pset) );
  _printers.push_back( new CaloRecoDigiPrinter(pset) );
  _printers.push_back( new CaloCrystalHitPrinter(pset) );
  _printers.push_back( new CaloClusterPrinter(pset) );
  _printers.push_back( new StrawDigiPrinter(pset) );
  _printers.push_back( new StrawDigiMCPrinter(pset) );
  _printers.push_back( new StrawHitPrinter(pset) );
  _printers.push_back( new BkgClusterPrinter(pset) );
  _printers.push_back( new BkgQualPrinter(pset) );
  _printers.push_back( new TrackClusterMatchPrinter(pset) );
  _printers.push_back( new TrkCaloIntersectPrinter(pset) );
  _printers.push_back( new TrackSummaryPrinter(pset) );
  _printers.push_back( new KalRepPrinter(pset) );
  _printers.push_back( new SimParticleTimeMapPrinter(pset) );
  _printers.push_back( new StrawHitFlagPrinter(pset) );
  _printers.push_back( new ComboHitPrinter(pset) );
  _printers.push_back( new TimeClusterPrinter(pset) );
  _printers.push_back( new KalSeedPrinter(pset) );
  _printers.push_back( new PhysicalVolumePrinter(pset) );
}


void mu2e::PrintModule::analyze(art::Event const& event) {
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

void mu2e::PrintModule::beginSubRun(art::SubRun const& subrun) {
  std::cout 
    << "\n"
    << " ###############  PrintModule Run/Subrun " 
    << std::setw(9) << subrun.run()
    << std::setw(9) << subrun.subRun()
    << std::endl;

  for(auto& prod_printer: _printers) prod_printer->PrintSubRun(subrun);

  std::cout << std::endl;

}


DEFINE_ART_MODULE(mu2e::PrintModule)
