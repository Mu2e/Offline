//
//  A module to print products in an event
//

#include <iomanip>
#include <memory>
#include <ostream>
#include <vector>

#include "Offline/Print/inc/BkgClusterPrinter.hh"
#include "Offline/Print/inc/BkgQualPrinter.hh"
#include "Offline/Print/inc/CaloClusterMCPrinter.hh"
#include "Offline/Print/inc/CaloClusterPrinter.hh"
#include "Offline/Print/inc/CaloDigiPrinter.hh"
#include "Offline/Print/inc/CaloHitMCPrinter.hh"
#include "Offline/Print/inc/CaloHitPrinter.hh"
#include "Offline/Print/inc/CaloRecoDigiPrinter.hh"
#include "Offline/Print/inc/CaloShowerStepPrinter.hh"
#include "Offline/Print/inc/ComboHitPrinter.hh"
#include "Offline/Print/inc/CosmicLivetimePrinter.hh"
#include "Offline/Print/inc/CrvCoincidenceClusterPrinter.hh"
#include "Offline/Print/inc/CrvDigiMCPrinter.hh"
#include "Offline/Print/inc/CrvDigiPrinter.hh"
#include "Offline/Print/inc/CrvRecoPulsePrinter.hh"
#include "Offline/Print/inc/CrvStepPrinter.hh"
#include "Offline/Print/inc/EventWindowMarkerPrinter.hh"
#include "Offline/Print/inc/GenParticlePrinter.hh"
#include "Offline/Print/inc/KalSeedPrinter.hh"
#include "Offline/Print/inc/MCTrajectoryPrinter.hh"
#include "Offline/Print/inc/PhysicalVolumePrinter.hh"
#include "Offline/Print/inc/PrimaryParticlePrinter.hh"
#include "Offline/Print/inc/ProductPrinter.hh"
#include "Offline/Print/inc/ProtonBunchIntensityPrinter.hh"
#include "Offline/Print/inc/ProtonBunchTimeMCPrinter.hh"
#include "Offline/Print/inc/ProtonBunchTimePrinter.hh"
#include "Offline/Print/inc/SimParticlePrinter.hh"
#include "Offline/Print/inc/SimParticlePtrPrinter.hh"
#include "Offline/Print/inc/StatusG4Printer.hh"
#include "Offline/Print/inc/StepPointMCPrinter.hh"
#include "Offline/Print/inc/STMWaveformDigiPrinter.hh"
#include "Offline/Print/inc/StrawDigiADCWaveformPrinter.hh"
#include "Offline/Print/inc/StrawDigiMCPrinter.hh"
#include "Offline/Print/inc/StrawDigiPrinter.hh"
#include "Offline/Print/inc/StrawGasStepPrinter.hh"
#include "Offline/Print/inc/SurfaceStepPrinter.hh"
#include "Offline/Print/inc/StrawHitFlagPrinter.hh"
#include "Offline/Print/inc/StrawHitPrinter.hh"
#include "Offline/Print/inc/TimeClusterPrinter.hh"
#include "Offline/Print/inc/HelixSeedPrinter.hh"
#include "Offline/Print/inc/CosmicTrackSeedPrinter.hh"
#include "Offline/Print/inc/TrackSummaryPrinter.hh"
#include "Offline/Print/inc/TriggerInfoPrinter.hh"
#include "Offline/Print/inc/TriggerResultsPrinter.hh"
#include "art/Framework/Core/EDAnalyzer.h"

using namespace std;

namespace mu2e {

class PrintModule : public art::EDAnalyzer {
 public:
  struct Config {
    fhicl::Atom<bool> printEvent { fhicl::Name("PrintEvent"), fhicl::Comment("Print Event Products"),true};
    fhicl::Atom<bool> printSubRun { fhicl::Name("PrintSubRun"), fhicl::Comment("Print SubRun Products"),true};
    fhicl::Atom<int> verbose{fhicl::Name("verbose"),fhicl::Comment("verbose flag, 0 to 1 or 2"), 1};

    fhicl::Table<ProductPrinter::Config> statusG4Printer{
        fhicl::Name("statusG4Printer")};
    fhicl::Table<ProductPrinter::Config> ProtonBunchTimePrinter{
        fhicl::Name("ProtonBunchTimePrinter")};
    fhicl::Table<ProductPrinter::Config> ProtonBunchTimeMCPrinter{
        fhicl::Name("ProtonBunchTimeMCPrinter")};
    fhicl::Table<ProductPrinter::Config> ProtonBunchIntensityPrinter{
        fhicl::Name("ProtonBunchIntensityPrinter")};
    fhicl::Table<ProductPrinter::Config> CosmicLivetimePrinter{
        fhicl::Name("CosmicLivetimePrinter")};
    fhicl::Table<ProductPrinter::Config> EventWindowMarkerPrinter{
        fhicl::Name("EventWindowMarkerPrinter")};
    fhicl::Table<ProductPrinter::Config> genParticlePrinter{
        fhicl::Name("genParticlePrinter")};
    fhicl::Table<SimParticlePrinter::Config> simParticlePrinter{
        fhicl::Name("simParticlePrinter")};
    fhicl::Table<ProductPrinter::Config> simParticlePtrPrinter{
        fhicl::Name("simParticlePtrPrinter")};
    fhicl::Table<StepPointMCPrinter::Config> stepPointMCPrinter{
        fhicl::Name("stepPointMCPrinter")};
    fhicl::Table<ProductPrinter::Config> mcTrajectoryPrinter{
        fhicl::Name("mcTrajectoryPrinter")};
    fhicl::Table<ProductPrinter::ConfigE> caloShowerStepPrinter{
        fhicl::Name("caloShowerStepPrinter")};
    fhicl::Table<ProductPrinter::Config> caloDigiPrinter{
        fhicl::Name("caloDigiPrinter")};
    fhicl::Table<ProductPrinter::ConfigE> caloRecoDigiPrinter{
        fhicl::Name("caloRecoDigiPrinter")};
    fhicl::Table<ProductPrinter::ConfigE> CaloHitPrinter{
        fhicl::Name("caloHitPrinter")};
    fhicl::Table<ProductPrinter::ConfigE> CaloHitMCPrinter{
        fhicl::Name("caloHitMCPrinter")};
    fhicl::Table<ProductPrinter::ConfigE> caloClusterPrinter{
        fhicl::Name("caloClusterPrinter")};
    fhicl::Table<ProductPrinter::ConfigE> caloClusterMCPrinter{
        fhicl::Name("caloClusterMCPrinter")};
    fhicl::Table<ProductPrinter::ConfigE> crvStepPrinter{
        fhicl::Name("crvStepPrinter")};
    fhicl::Table<ProductPrinter::Config> crvDigiPrinter{
        fhicl::Name("crvDigiPrinter")};
    fhicl::Table<ProductPrinter::Config> crvDigiMCPrinter{
        fhicl::Name("crvDigiMCPrinter")};
    fhicl::Table<ProductPrinter::Config> crvRecoPulsePrinter{
        fhicl::Name("crvRecoPulsePrinter")};
    fhicl::Table<ProductPrinter::Config> crvCoincidenceClusterPrinter{
        fhicl::Name("crvCoincidenceClusterPrinter")};
    fhicl::Table<ProductPrinter::Config> strawGasStepPrinter{
        fhicl::Name("strawGasStepPrinter")};
    fhicl::Table<ProductPrinter::Config> surfaceStepPrinter{
        fhicl::Name("surfaceStepPrinter")};
    fhicl::Table<ProductPrinter::Config> strawDigiPrinter{
        fhicl::Name("strawDigiPrinter")};
    fhicl::Table<ProductPrinter::Config> strawDigiADCWaveformPrinter{
        fhicl::Name("strawDigiADCWaveformPrinter")};
    fhicl::Table<ProductPrinter::Config> strawDigiMCPrinter{
        fhicl::Name("strawDigiMCPrinter")};
    fhicl::Table<ProductPrinter::ConfigE> strawHitPrinter{
        fhicl::Name("strawHitPrinter")};
    fhicl::Table<ProductPrinter::Config> strawHitFlagPrinter{
        fhicl::Name("strawHitFlagPrinter")};
    fhicl::Table<ProductPrinter::Config> bkgClusterPrinter{
        fhicl::Name("bkgClusterPrinter")};
    fhicl::Table<ProductPrinter::ConfigE> bkgQualPrinter{
        fhicl::Name("bkgQualPrinter")};
    fhicl::Table<ProductPrinter::Config> trkCaloIntersectPrinter{
        fhicl::Name("trkCaloIntersectPrinter")};
    fhicl::Table<ProductPrinter::Config> trackSummaryPrinter{
        fhicl::Name("trackSummaryPrinter")};
    fhicl::Table<ProductPrinter::Config> comboHitPrinter{
        fhicl::Name("comboHitPrinter")};
    fhicl::Table<ProductPrinter::Config> timeClusterPrinter{
        fhicl::Name("timeClusterPrinter")};
    fhicl::Table<ProductPrinter::Config> helixSeedPrinter{
        fhicl::Name("helixSeedPrinter")};
    fhicl::Table<ProductPrinter::Config> cosmicTrackSeedPrinter{
        fhicl::Name("cosmicTrackSeedPrinter")};
    fhicl::Table<ProductPrinter::Config> kalSeedPrinter{
        fhicl::Name("kalSeedPrinter")};
    fhicl::Table<ProductPrinter::Config> stmWaveformDigiPrinter{
        fhicl::Name("stmWaveformDigiPrinter")};
    fhicl::Table<ProductPrinter::Config> physicalVolumePrinter{
        fhicl::Name("physicalVolumePrinter")};
    fhicl::Table<ProductPrinter::Config> triggerResultsPrinter{
        fhicl::Name("triggerResultsPrinter")};
    fhicl::Table<ProductPrinter::Config> triggerInfoPrinter{
        fhicl::Name("triggerInfoPrinter")};
    fhicl::Table<ProductPrinter::Config> primaryParticlePrinter{
        fhicl::Name("primaryParticlePrinter")};
  };

  // this line is required by art to allow the command line help print
  typedef art::EDAnalyzer::Table<Config> Parameters;

  explicit PrintModule(const Parameters& conf);
  void analyze(art::Event const& event) override;
  void beginSubRun(art::SubRun const& subrun) override;
  void endJob() override;

 private:
  bool _printevent, _printsubrun;
  int _verbose;
  // each of these object prints a different product
  vector<unique_ptr<mu2e::ProductPrinter> > _printers;
};

}  // namespace mu2e

mu2e::PrintModule::PrintModule(const Parameters& conf) : art::EDAnalyzer(conf),
  _printevent(conf().printEvent()), _printsubrun(conf().printSubRun()), _verbose(conf().verbose()){

  _printers.push_back(make_unique<StatusG4Printer>(conf().statusG4Printer()));
  _printers.push_back(
      make_unique<ProtonBunchTimePrinter>(conf().ProtonBunchTimePrinter()));
  _printers.push_back(
      make_unique<ProtonBunchTimeMCPrinter>(conf().ProtonBunchTimeMCPrinter()));
  _printers.push_back(make_unique<ProtonBunchIntensityPrinter>(
      conf().ProtonBunchIntensityPrinter()));
  _printers.push_back(
      make_unique<CosmicLivetimePrinter>(conf().CosmicLivetimePrinter()));
  _printers.push_back(
      make_unique<EventWindowMarkerPrinter>(conf().EventWindowMarkerPrinter()));
  _printers.push_back(
      make_unique<GenParticlePrinter>(conf().genParticlePrinter()));
  _printers.push_back(
      make_unique<SimParticlePrinter>(conf().simParticlePrinter()));
  _printers.push_back(
      make_unique<SimParticlePtrPrinter>(conf().simParticlePtrPrinter()));
  _printers.push_back(
      make_unique<StepPointMCPrinter>(conf().stepPointMCPrinter()));
  _printers.push_back(
      make_unique<MCTrajectoryPrinter>(conf().mcTrajectoryPrinter()));
  _printers.push_back(
      make_unique<CaloShowerStepPrinter>(conf().caloShowerStepPrinter()));
  _printers.push_back(make_unique<CaloDigiPrinter>(conf().caloDigiPrinter()));
  _printers.push_back(
      make_unique<CaloRecoDigiPrinter>(conf().caloRecoDigiPrinter()));
  _printers.push_back(make_unique<CaloHitPrinter>(conf().CaloHitPrinter()));
  _printers.push_back(make_unique<CaloHitMCPrinter>(conf().CaloHitMCPrinter()));
  _printers.push_back(
      make_unique<CaloClusterPrinter>(conf().caloClusterPrinter()));
  _printers.push_back(
      make_unique<CaloClusterMCPrinter>(conf().caloClusterMCPrinter()));
  _printers.push_back(make_unique<CrvStepPrinter>(conf().crvStepPrinter()));
  _printers.push_back(make_unique<CrvDigiPrinter>(conf().crvDigiPrinter()));
  _printers.push_back(make_unique<CrvDigiMCPrinter>(conf().crvDigiMCPrinter()));
  _printers.push_back(
      make_unique<CrvRecoPulsePrinter>(conf().crvRecoPulsePrinter()));
  _printers.push_back(make_unique<CrvCoincidenceClusterPrinter>(
      conf().crvCoincidenceClusterPrinter()));
  _printers.push_back(
      make_unique<StrawGasStepPrinter>(conf().strawGasStepPrinter()));
  _printers.push_back(
      make_unique<SurfaceStepPrinter>(conf().surfaceStepPrinter()));
  _printers.push_back(make_unique<StrawDigiPrinter>(conf().strawDigiPrinter()));
  _printers.push_back(make_unique<StrawDigiADCWaveformPrinter>(
      conf().strawDigiADCWaveformPrinter()));
  _printers.push_back(
      make_unique<StrawDigiMCPrinter>(conf().strawDigiMCPrinter()));
  _printers.push_back(make_unique<StrawHitPrinter>(conf().strawHitPrinter()));
  _printers.push_back(
      make_unique<StrawHitFlagPrinter>(conf().strawHitFlagPrinter()));
  _printers.push_back(
      make_unique<BkgClusterPrinter>(conf().bkgClusterPrinter()));
  _printers.push_back(make_unique<BkgQualPrinter>(conf().bkgQualPrinter()));
  _printers.push_back(
      make_unique<TrackSummaryPrinter>(conf().trackSummaryPrinter()));
  _printers.push_back(make_unique<ComboHitPrinter>(conf().comboHitPrinter()));
  _printers.push_back(make_unique<TimeClusterPrinter>(conf().timeClusterPrinter()));
  _printers.push_back(make_unique<HelixSeedPrinter>(conf().helixSeedPrinter()));
  _printers.push_back(make_unique<CosmicTrackSeedPrinter>(conf().cosmicTrackSeedPrinter()));
  _printers.push_back(make_unique<KalSeedPrinter>(conf().kalSeedPrinter()));
  _printers.push_back(make_unique<STMWaveformDigiPrinter>(
      conf().stmWaveformDigiPrinter()));
  _printers.push_back(
      make_unique<PhysicalVolumePrinter>(conf().physicalVolumePrinter()));
  _printers.push_back(
      make_unique<TriggerResultsPrinter>(conf().triggerResultsPrinter()));
  _printers.push_back(
      make_unique<TriggerInfoPrinter>(conf().triggerInfoPrinter()));
  _printers.push_back(
      make_unique<PrimaryParticlePrinter>(conf().primaryParticlePrinter()));
}

void mu2e::PrintModule::analyze(art::Event const& event) {
  if(_printevent) {
    if(_verbose > 0){
      cout << "\n"
        << " ###############  PrintModule Run/Subrun/Event " << setw(9)
        << event.run() << setw(9) << event.subRun() << setw(9) << event.event()
        << endl;
    }

    for (auto& prod_printer : _printers) prod_printer->Print(event);

    if(_verbose > 0)cout << endl;
  }
}

void mu2e::PrintModule::beginSubRun(art::SubRun const& subrun) {
  if(_printsubrun) {
    if(_verbose > 0){
      cout << "\n"
        << " ###############  PrintModule Run/Subrun " << setw(9) << subrun.run()
        << setw(9) << subrun.subRun() << endl;
    }

    for (auto& prod_printer : _printers) prod_printer->PrintSubRun(subrun);

    if(_verbose > 0)cout << endl;
  }
}

void mu2e::PrintModule::endJob() {
  for (auto& prod_printer : _printers) prod_printer->PrintEndJob();
  cout << endl;
}

DEFINE_ART_MODULE(mu2e::PrintModule)
