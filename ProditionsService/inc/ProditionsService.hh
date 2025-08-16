#ifndef ProditionsService_ProditionsService_hh
#define ProditionsService_ProditionsService_hh

//
// Service to hold and deliver time-dependent and
// database-backed conditions quantities
//

#include "Offline/AnalysisConfig/inc/MVACatalogConfig.hh"
#include "Offline/CRVConfig/inc/CRVCalibConfig.hh"
#include "Offline/CRVConfig/inc/CRVOrdinalConfig.hh"
#include "Offline/CRVConfig/inc/CRVPhotonYieldConfig.hh"
#include "Offline/CRVConfig/inc/CRVStatusConfig.hh"
#include "Offline/CaloConfig/inc/CaloDAQMapConfig.hh"
#include "Offline/DAQConfig/inc/EventTimingConfig.hh"
#include "Offline/Mu2eInterfaces/inc/ProditionsCache.hh"
#include "Offline/Mu2eInterfaces/inc/ProditionsEntity.hh"
#include "Offline/STMConfig/inc/STMEnergyCalibConfig.hh"
#include "Offline/SimulationConfig/inc/SimBookkeeperConfig.hh"
#include "Offline/TrackerConfig/inc/AlignedTrackerConfig.hh"
#include "Offline/TrackerConfig/inc/FullReadoutStrawConfig.hh"
#include "Offline/TrackerConfig/inc/Mu2eDetectorConfig.hh"
#include "Offline/TrackerConfig/inc/Mu2eMaterialConfig.hh"
#include "Offline/TrackerConfig/inc/StrawDriftConfig.hh"
#include "Offline/TrackerConfig/inc/StrawElectronicsConfig.hh"
#include "Offline/TrackerConfig/inc/StrawPhysicsConfig.hh"
#include "Offline/TrackerConfig/inc/StrawResponseConfig.hh"
#include "Offline/TrackerConfig/inc/TrackerStatusConfig.hh"
#include "Offline/TrackerConfig/inc/TrackerPanelMapConfig.hh"

#include "Offline/AnalysisConfig/inc/MVACatalogConfig.hh"
#include "Offline/SimulationConfig/inc/SimBookkeeperConfig.hh"
#include "Offline/CaloConfig/inc/CalCalibConfig.hh"

#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include <string>


namespace mu2e {

class ProditionsService {
 public:
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    fhicl::Atom<int> verbose{Name("verbose"), Comment("verbosity 0 or 1"), 0};
    fhicl::Table<CRVOrdinalConfig> crvOrdinal{
        Name("crvOrdinal"),
        Comment("CRV online-offline numbering configuration")};
    fhicl::Table<CRVPhotonYieldConfig> crvPhotonYield{
        Name("crvPhotonYield"),
        Comment("CRV photon yield deviation map for each channel")};
    fhicl::Table<CRVStatusConfig> crvStatus{
        Name("crvStatus"), Comment("CRV bad channels configuration")};
    fhicl::Table<CRVCalibConfig> crvCalib{
        Name("crvCalib"), Comment("CRV SiPM calibration configuration")};
    fhicl::Table<EventTimingConfig> eventTiming{
        Name("eventTiming"), Comment("Event timing configuration")};
    fhicl::Table<STMEnergyCalibConfig> stmEnergyCalib{
        Name("stmEnergyCalib"),
        Comment("STM energy calibration configuration")};
    fhicl::Table<FullReadoutStrawConfig> fullReadoutStraw{
        Name("fullReadoutStraw"),
        Comment("Straws with no time window in readout")};
    fhicl::Table<TrackerStatusConfig> trackerStatus{
        Name("trackerStatus"),
        Comment("Status of tracker elements (straws, panels, planes, ...)")};
    fhicl::Table<TrackerPanelMapConfig> trackerPanelMap {
      Name("trackerPanelMap"),Comment("TrackerPanelMap config")};
    fhicl::Table<StrawDriftConfig> strawDrift{
        Name("strawDrift"), Comment("Straw drift model function and binning")};
    fhicl::Table<StrawPhysicsConfig> strawPhysics{
        Name("strawPhysics"), Comment("Straw physics model")};
    fhicl::Table<StrawElectronicsConfig> strawElectronics{
        Name("strawElectronics"), Comment("Straw electronics model")};
    fhicl::Table<StrawResponseConfig> strawResponse{
        Name("strawResponse"), Comment("Straw response model")};
    fhicl::Table<AlignedTrackerConfig> alignedTracker{
        Name("alignedTracker"), Comment("Tracker alignment in reco code")};
    fhicl::Table<AlignedTrackerConfig> alignedTrackerSim{
        Name("alignedTrackerSim"), Comment("Tracker alignment in sim code")};
    fhicl::Table<Mu2eMaterialConfig> mu2eMaterial{
        Name("mu2eMaterial"), Comment("Mu2e material for BTrk")};
    fhicl::Table<Mu2eDetectorConfig> mu2eDetector{
        Name("mu2eDetector"), Comment("Mu2e detector model for BTrk")};
    fhicl::Table<CaloDAQMapConfig> caloDAQConditions{
        Name("caloDAQConditions"), Comment("DAQ channel maps for calorimeter")};
    fhicl::Table<MVACatalogConfig> trkQualCatalog{
        Name("trkQualCatalog"), Comment("Catalog of TrkQual trainings")};
    fhicl::Table<SimBookkeeperConfig> simbookkeeper{
        Name("simbookkeeper"), Comment("simulation bookkeeping")};
    fhicl::Table<CalCalibConfig> calCalib{
          Name("calCalib"),
          Comment("calorimeter energy calib") };
          
  };

  // this line is required by art to allow the command line help print
  typedef art::ServiceTable<Config> Parameters;

  ProditionsService(Parameters const& config, art::ActivityRegistry& iRegistry);

  ProditionsCache::ptr getCache(std::string name) {
    if (_caches.count(name) == 0) return ProditionsCache::ptr();
    return _caches[name];
  }

  // void postBeginJob();

 private:
  // This is not copyable or assignable - private and unimplemented.
  ProditionsService const& operator=(ProditionsService const& rhs);
  ProditionsService(ProditionsService const& rhs);

  Config _config;
  std::map<std::string, ProditionsCache::ptr> _caches;
};

}  // namespace mu2e

DECLARE_ART_SERVICE(mu2e::ProditionsService, SHARED)
#endif /* ProditionsService_ProditionsService_hh */
