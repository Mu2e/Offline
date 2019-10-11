//
//  A module to print products in an event
//

#include <vector>
#include <ostream>
#include <iomanip>
#include <memory>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
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
#include "Print/inc/CrvDigiPrinter.hh"
#include "Print/inc/CrvDigiMCPrinter.hh"
#include "Print/inc/CrvRecoPulsePrinter.hh"
#include "Print/inc/CrvCoincidenceClusterPrinter.hh"
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
#include "Print/inc/TriggerResultsPrinter.hh"

using namespace std;

namespace mu2e {

  class PrintModule : public art::EDAnalyzer {

  public:

    struct Config {
      fhicl::Table<ProductPrinter::Config> statusG4Printer { 
	fhicl::Name("statusG4Printer") }; 
      fhicl::Table<ProductPrinter::Config> genParticlePrinter { 
	fhicl::Name("genParticlePrinter") }; 
      fhicl::Table<SimParticlePrinter::Config> simParticlePrinter { 
	fhicl::Name("simParticlePrinter") }; 
      fhicl::Table<ProductPrinter::Config> simParticlePtrPrinter { 
	fhicl::Name("simParticlePtrPrinter") }; 
      fhicl::Table<StepPointMCPrinter::Config> stepPointMCPrinter { 
	fhicl::Name("stepPointMCPrinter") }; 
      fhicl::Table<ProductPrinter::Config> mcTrajectoryPrinter { 
	fhicl::Name("mcTrajectoryPrinter") }; 
      fhicl::Table<ProductPrinter::ConfigE> caloShowerStepPrinter { 
	fhicl::Name("caloShowerStepPrinter") }; 
      fhicl::Table<ProductPrinter::ConfigE> caloHitPrinter { 
	fhicl::Name("caloHitPrinter") }; 
      fhicl::Table<ProductPrinter::Config> caloDigiPrinter { 
	fhicl::Name("caloDigiPrinter") }; 
      fhicl::Table<ProductPrinter::ConfigE> caloRecoDigiPrinter { 
	fhicl::Name("caloRecoDigiPrinter") }; 
      fhicl::Table<ProductPrinter::ConfigE> caloCrystalHitPrinter { 
	fhicl::Name("caloCrystalHitPrinter") }; 
      fhicl::Table<ProductPrinter::ConfigE> caloClusterPrinter { 
	fhicl::Name("caloClusterPrinter") }; 
      fhicl::Table<ProductPrinter::Config> crvDigiPrinter { 
	fhicl::Name("crvDigiPrinter") }; 
      fhicl::Table<ProductPrinter::Config> crvDigiMCPrinter { 
	fhicl::Name("crvDigiMCPrinter") }; 
      fhicl::Table<ProductPrinter::Config> crvRecoPulsePrinter { 
	fhicl::Name("crvRecoPulsePrinter") }; 
      fhicl::Table<ProductPrinter::Config> crvCoincidenceClusterPrinter { 
	fhicl::Name("crvCoincidenceClusterPrinter") }; 
      fhicl::Table<ProductPrinter::Config> strawDigiPrinter { 
      	fhicl::Name("strawDigiPrinter") }; 
      fhicl::Table<ProductPrinter::Config> strawDigiMCPrinter { 
      	fhicl::Name("strawDigiMCPrinter") }; 
      fhicl::Table<ProductPrinter::ConfigE> strawHitPrinter { 
      	fhicl::Name("strawHitPrinter") }; 
      fhicl::Table<ProductPrinter::Config> strawHitFlagPrinter { 
      	fhicl::Name("strawHitFlagPrinter") }; 
      fhicl::Table<ProductPrinter::Config> bkgClusterPrinter { 
	fhicl::Name("bkgClusterPrinter") }; 
      fhicl::Table<ProductPrinter::ConfigE> bkgQualPrinter { 
	fhicl::Name("bkgQualPrinter") }; 


      fhicl::Table<ProductPrinter::Config> trackClusterMatchPrinter { 
	fhicl::Name("trackClusterMatchPrinter") }; 
      fhicl::Table<ProductPrinter::Config> trkCaloIntersectPrinter { 
	fhicl::Name("trkCaloIntersectPrinter") }; 
      fhicl::Table<ProductPrinter::Config> trackSummaryPrinter { 
	fhicl::Name("trackSummaryPrinter") }; 
      fhicl::Table<ProductPrinter::Config> kalRepPrinter { 
	fhicl::Name("kalRepPrinter") }; 
      fhicl::Table<ProductPrinter::Config> simParticleTimeMapPrinter { 
	fhicl::Name("simParticleTimeMapPrinter") }; 
      fhicl::Table<ProductPrinter::Config> comboHitPrinter { 
	fhicl::Name("comboHitPrinter") }; 
      fhicl::Table<ProductPrinter::Config> timeClusterPrinter { 
	fhicl::Name("timeClusterPrinter") }; 
      fhicl::Table<ProductPrinter::Config> kalSeedPrinter { 
	fhicl::Name("kalSeedPrinter") }; 
      fhicl::Table<ProductPrinter::Config> physicalVolumePrinter { 
	fhicl::Name("physicalVolumePrinter") }; 
      fhicl::Table<ProductPrinter::Config> triggerResultsPrinter { 
	fhicl::Name("triggerResultsPrinter") }; 

    };

    // this line is required by art to allow the command line help print
    typedef art::EDAnalyzer::Table<Config> Parameters;

    explicit PrintModule(const Parameters& conf);
    void analyze  ( art::Event const&  event  ) override;
    void beginSubRun( art::SubRun const& subrun) override;

  private:

    // set by fcl, turn all prints on (1) or off (0)
    int _verbose;
    // each of these object prints a different product
    vector< unique_ptr<mu2e::ProductPrinter> > _printers;
  };

}


mu2e::PrintModule::PrintModule(const Parameters& conf):
  art::EDAnalyzer(conf) {
  //cout << "start main pset\n"<< pset.to_string() << "\n end main pset"<< endl;

  _printers.push_back( make_unique<StatusG4Printer>( conf().statusG4Printer() ) );
  _printers.push_back( make_unique<GenParticlePrinter>( conf().genParticlePrinter() ) );
  _printers.push_back( make_unique<SimParticlePrinter>( conf().simParticlePrinter() ) );
  _printers.push_back( make_unique<SimParticlePtrPrinter>( conf().simParticlePtrPrinter() ) );
  _printers.push_back( make_unique<StepPointMCPrinter>( conf().stepPointMCPrinter() ) );
  _printers.push_back( make_unique<MCTrajectoryPrinter>( conf().mcTrajectoryPrinter() ) );
  _printers.push_back( make_unique<CaloShowerStepPrinter>( conf().caloShowerStepPrinter() ) );
  _printers.push_back( make_unique<CaloHitPrinter>( conf().caloHitPrinter() ) );
  _printers.push_back( make_unique<CaloDigiPrinter>( conf().caloDigiPrinter() ) );
  _printers.push_back( make_unique<CaloRecoDigiPrinter>( conf().caloRecoDigiPrinter() ) );
  _printers.push_back( make_unique<CaloCrystalHitPrinter>( conf().caloCrystalHitPrinter() ) );
  _printers.push_back( make_unique<CaloClusterPrinter>( conf().caloClusterPrinter() ) );
  _printers.push_back( make_unique<CrvDigiPrinter>( conf().crvDigiPrinter() ) );
  _printers.push_back( make_unique<CrvDigiMCPrinter>( conf().crvDigiMCPrinter() ) );
  _printers.push_back( make_unique<CrvRecoPulsePrinter>( conf().crvRecoPulsePrinter() ) );
  _printers.push_back( make_unique<CrvCoincidenceClusterPrinter>( conf().crvCoincidenceClusterPrinter() ) );
  _printers.push_back( make_unique<StrawDigiPrinter>( conf().strawDigiPrinter() ) );
  _printers.push_back( make_unique<StrawDigiMCPrinter>( conf().strawDigiMCPrinter() ) );
  _printers.push_back( make_unique<StrawHitPrinter>( conf().strawHitPrinter() ) );
  _printers.push_back( make_unique<StrawHitFlagPrinter>( conf().strawHitFlagPrinter() ) );

  _printers.push_back( make_unique<BkgClusterPrinter>( conf().bkgClusterPrinter() ) );
  _printers.push_back( make_unique<BkgQualPrinter>( conf().bkgQualPrinter() ) );
  _printers.push_back( make_unique<TrackClusterMatchPrinter>( conf().trackClusterMatchPrinter() ) );
  _printers.push_back( make_unique<TrkCaloIntersectPrinter>( conf().trkCaloIntersectPrinter() ) );
  _printers.push_back( make_unique<TrackSummaryPrinter>( conf().trackSummaryPrinter() ) );
  _printers.push_back( make_unique<KalRepPrinter>( conf().kalRepPrinter() ) );
  _printers.push_back( make_unique<SimParticleTimeMapPrinter>( conf().simParticleTimeMapPrinter() ) );
  _printers.push_back( make_unique<ComboHitPrinter>( conf().comboHitPrinter() ) );
  _printers.push_back( make_unique<TimeClusterPrinter>( conf().timeClusterPrinter() ) );
  _printers.push_back( make_unique<KalSeedPrinter>( conf().kalSeedPrinter() ) );
  _printers.push_back( make_unique<PhysicalVolumePrinter>( conf().physicalVolumePrinter() ) );
  _printers.push_back( make_unique<TriggerResultsPrinter>( conf().triggerResultsPrinter() ) );
}


void mu2e::PrintModule::analyze(art::Event const& event) {
  cout 
    << "\n"
    << " ###############  PrintModule Run/Subrun/Event " 
    << setw(9) << event.run()
    << setw(9) << event.subRun()
    << setw(9) << event.event()
    << endl;

  for(auto& prod_printer: _printers) prod_printer->Print(event);

  cout << endl;

}

void mu2e::PrintModule::beginSubRun(art::SubRun const& subrun) {
  cout 
    << "\n"
    << " ###############  PrintModule Run/Subrun " 
    << setw(9) << subrun.run()
    << setw(9) << subrun.subRun()
    << endl;

  for(auto& prod_printer: _printers) prod_printer->PrintSubRun(subrun);

  cout << endl;

}


DEFINE_ART_MODULE(mu2e::PrintModule)
