//
// The btrk code needs some information and tools that are supplied externally.
// This service is the agent that supplies the that information and those tools.
//
// Original author Rob Kutschke
//

#include "BTrkHelper/inc/BTrkHelper.hh"

#include "BaBar/BaBar/include/ExternalInfo.hh"

mu2e::BTrkHelper::BTrkHelper(fhicl::ParameterSet const& pset,
                             art::ActivityRegistry&     registry):
  particleInfo_(),
  fileFinder_(pset.get<fhicl::ParameterSet>("dictionaries"))
{

  // Register callbacks.
  registry.sPreBeginRun.watch   ( this, &BTrkHelper::beginRun    );
  registry.sPreBeginSubRun.watch( this, &BTrkHelper::beginSubRun );

  // These are available at c'tor time - so push them to the BaBar code now.
  ExternalInfo::set( &fileFinder_);
  ExternalInfo::set( &particleInfo_);
}

void
mu2e::BTrkHelper::beginRun(art::Run const &) {
  // When we have information that needs to be pushed at beginSubRun time, then this is where to do it.
  // Follow the same model for any information that needs to be pushed at beginRun time.
}

void
mu2e::BTrkHelper::beginSubRun(art::SubRun const &) {
  // When we have information that needs to be pushed at beginSubRun time, then this is where to do it.
  // Follow the same model for any information that needs to be pushed at beginRun time.
}

DEFINE_ART_SERVICE(mu2e::BTrkHelper);
