#ifndef DbExample_ProditionsService_hh
#define DbExample_ProditionsService_hh

//
// Service to hold and deliver time-dependent and 
// database-backed conditions quantities
//

#include <string>

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "cetlib_except/exception.h"

#include "Mu2eInterfaces/inc/ProditionsEntity.hh"
#include "Mu2eInterfaces/inc/ProditionsCache.hh"

#include "TrackerConfig/inc/FullReadoutStrawConfig.hh"
#include "TrackerConfig/inc/DeadStrawConfig.hh"
#include "TrackerConfig/inc/StrawDriftConfig.hh"
#include "TrackerConfig/inc/StrawPhysicsConfig.hh"
#include "TrackerConfig/inc/StrawElectronicsConfig.hh"
#include "TrackerConfig/inc/StrawResponseConfig.hh"
#include "TrackerConfig/inc/AlignedTrackerConfig.hh"
#include "TrackerConfig/inc/Mu2eMaterialConfig.hh"
#include "TrackerConfig/inc/Mu2eDetectorConfig.hh"


namespace mu2e {

  class ProditionsService {

  public:

    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<int> verbose{Name("verbose"),
	  Comment("verbosity 0 or 1"),0};
      fhicl::Table<FullReadoutStrawConfig> fullReadoutStraw{
	  Name("fullReadoutStraw"), 
	  Comment("Straws with no time window in readout") };
      fhicl::Table<DeadStrawConfig> deadStraw{
	  Name("deadStraw"), 
	  Comment("Dead Straw List by Plane, Panel and Straw") };
      fhicl::Table<StrawDriftConfig> strawDrift{
	  Name("strawDrift"), 
	  Comment("Straw drift model function and binning") };
      fhicl::Table<StrawPhysicsConfig> strawPhysics{
	  Name("strawPhysics"), 
	  Comment("Straw physics model") };
      fhicl::Table<StrawElectronicsConfig> strawElectronics{
	  Name("strawElectronics"), 
	  Comment("Straw electronics model") };
      fhicl::Table<StrawResponseConfig> strawResponse{
	  Name("strawResponse"), 
	  Comment("Straw response model") };
      fhicl::Table<AlignedTrackerConfig> alignedTracker{
	  Name("alignedTracker"), 
	  Comment("Tracker alignment in reco code") };
      fhicl::Table<Mu2eMaterialConfig> mu2eMaterial{
	  Name("mu2eMaterial"), 
	  Comment("Mu2e material for BTrk") };
      fhicl::Table<Mu2eDetectorConfig> mu2eDetector{
	  Name("mu2eDetector"), 
	  Comment("Mu2e detector model for BTrk") };

    };

    // this line is required by art to allow the command line help print
    typedef art::ServiceTable<Config> Parameters;

    ProditionsService(Parameters const& config, 
		       art::ActivityRegistry& iRegistry);

    ProditionsCache::ptr getCache(std::string name) { 
      return _caches[name];
    }
    //void postBeginJob();

  private:

    // This is not copyable or assignable - private and unimplemented.
    ProditionsService const& operator=(ProditionsService const& rhs);
    ProditionsService(ProditionsService const& rhs);

    Config _config;
    std::map<std::string,ProditionsCache::ptr> _caches;
  };

}

DECLARE_ART_SERVICE(mu2e::ProditionsService, SHARED)
#endif /* ProditionsService_ProditionsService_hh */
