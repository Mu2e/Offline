// Pixel digitization: create ExtMonFNALRawClusters from raw hits.
//
//
// Original author Andrei Gaponenko
//

#include <string>
#include <cmath>
#include <memory>
#include <iostream>
#include <iterator>
#include <algorithm>


#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/Ptr.h"

#include "RecoDataProducts/inc/ExtMonFNALRawHit.hh"
#include "RecoDataProducts/inc/ExtMonFNALRawHitCollection.hh"
#include "RecoDataProducts/inc/ExtMonFNALRawCluster.hh"
#include "RecoDataProducts/inc/ExtMonFNALRawClusterCollection.hh"

#include "GeometryService/inc/GeomHandle.hh"

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALModule.hh"

#include "ExtinctionMonitorFNAL/Geometry/inc/PixelNeighbors.hh"
#include "ExtinctionMonitorFNAL/Reconstruction/inc/PixelHitLookup.hh"

namespace mu2e {

  //================================================================
  class ExtMonFNALRawClusterization : public art::EDProducer {

  public:
    explicit ExtMonFNALRawClusterization(fhicl::ParameterSet const& pset)
      : EDProducer{pset}
      , verbosityLevel_(pset.get<int>("verbosityLevel", 0))
      , inputModuleLabel_(pset.get<std::string>("inputModuleLabel"))
      , inputInstanceName_(pset.get<std::string>("inputInstanceName", ""))
      , geomModuleLabel_(pset.get<std::string>("geomModuleLabel"))
      , geomInstanceName_(pset.get<std::string>("geomInstanceName", ""))
      , extmon_(0)
    {
      produces<ExtMonFNALRawClusterCollection>();
    }

    virtual void produce(art::Event& evt);
    virtual void beginRun(art::Run& run);

  private:
    int verbosityLevel_;
    std::string inputModuleLabel_;
    std::string inputInstanceName_;
    std::string geomModuleLabel_;
    std::string geomInstanceName_;

    // Non-owning pointers to the geometry and conditions objects. The
    // current Mu2e infrastructure does not allow the use of a Handle
    // as a class member.
    const ExtMonFNAL::ExtMon *extmon_;

    void formClusters(ExtMonFNALRawClusterCollection *clusters,
                      const art::Handle<ExtMonFNALRawHitCollection>& hits,
                      const art::EDProductGetter* hitsGetter
                      );

  };

  //================================================================
  void ExtMonFNALRawClusterization::beginRun(art::Run& run) {
    if(!geomModuleLabel_.empty()) {
      if(verbosityLevel_ > 0) {
        std::cout<<"ExtMonFNALRawClusterization: using recorded geometry: "
                 <<"("<<geomModuleLabel_<<", "<<geomInstanceName_<<")"
                 <<std::endl;
      }
      art::Handle<ExtMonFNAL::ExtMon> emf;
      run.getByLabel(geomModuleLabel_, geomInstanceName_, emf);
      extmon_ = &*emf;
    }
    else {
      if(verbosityLevel_ > 0) {
        std::cout<<"ExtMonFNALRawClusterization: using GeometryService"<<std::endl;
      }
      GeomHandle<ExtMonFNAL::ExtMon> emf;
      extmon_ = &*emf;
    }
  }

  //================================================================
  void ExtMonFNALRawClusterization::produce(art::Event& event) {

    std::unique_ptr<ExtMonFNALRawClusterCollection> clusters(new ExtMonFNALRawClusterCollection());

    art::Handle<ExtMonFNALRawHitCollection> hits;
    event.getByLabel(inputModuleLabel_, inputInstanceName_, hits);

    const art::EDProductGetter *hitsGetter = event.productGetter(hits.id());
    formClusters(&*clusters, hits, hitsGetter);

    event.put(std::move(clusters));
  }

  //================================================================
  void ExtMonFNALRawClusterization::formClusters(ExtMonFNALRawClusterCollection *clusters,
                                              const art::Handle<ExtMonFNALRawHitCollection>& hitsHandle,
                                              const art::EDProductGetter* hitsGetter
                                              )
  {
    const ExtMonFNALRawHitCollection& hits(*hitsHandle);

    PixelNeighbors pn(extmon_->module(), extmon_->chip());
    std::set<std::size_t> used;

    PixelHitLookup pixmap(hits);

    // Go through all hits
    for(unsigned iseed = 0; iseed < hits.size(); ++iseed) {

      // If the hit is not yet used in a cluster
      if(used.insert(iseed).second) {

        // start a new cluster
        ExtMonFNALRawCluster::Hits clusterHits;
        clusterHits.push_back(art::Ptr<ExtMonFNALRawHit>(hitsHandle.id(), iseed, hitsGetter));

        unsigned previousClusterSize = 0;

        // Iteratively grow the cluster, merging neighbour hits with
        // the correct time.  We define clusters in a commutative way:
        // a clusters with b iff b clusters with a, therefore it does
        // not matter what pixel served as a seed - the final set of
        // clusters is invariant w.r.t. the input collection ordering.

        while(clusterHits.size() != previousClusterSize) {
          previousClusterSize = clusterHits.size();

          // Go through all pixels already in the cluster and try to
          // find any neighbors that should be joined.

          for(unsigned i=0; i < clusterHits.size(); ++i) {

            const PixelNeighbors::Collection nb(pn.neighbors(clusterHits[i]->pixelId()));

            for(PixelNeighbors::Collection::const_iterator nid = nb.begin(); nid != nb.end(); ++nid) {

              typedef PixelHitLookup::HitIndex HitIndex;

              // We only cluster hits in the same time bin.
              // This may need to change if timewalk is important.
              //
              // Note: clock reference here is from the original seed.
              HitIndex candidate = pixmap.findHit(*nid, hits[iseed].clock());

              if( (candidate != HitIndex(-1)) && used.insert(candidate).second) {
                clusterHits.push_back(art::Ptr<ExtMonFNALRawHit>(hitsHandle.id(), candidate, hitsGetter));
              }

            } // for(all neighboring pixels)
          } // for(all pixels already in the cluster)
        }

        // Compute parameters and add cluster to the output
        clusters->push_back(ExtMonFNALRawCluster(clusterHits));

      } // if(hit is not yet used)
    } // for(all input hits)

  } // formClusters()

  //================================================================

} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::ExtMonFNALRawClusterization)
