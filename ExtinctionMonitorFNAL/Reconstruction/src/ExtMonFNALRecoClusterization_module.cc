// Compute calibrated pixel clusters from raw ones.
//
//
// Original author Andrei Gaponenko
//

#include <string>
#include <iostream>
#include <cmath>
#include <memory>
#include <iterator>
#include <algorithm>
#include <cassert>

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
#include "canvas/Persistency/Common/Ptr.h"

#include "RecoDataProducts/inc/ExtMonFNALRawHit.hh"
#include "RecoDataProducts/inc/ExtMonFNALRawCluster.hh"
#include "RecoDataProducts/inc/ExtMonFNALRawClusterCollection.hh"
#include "RecoDataProducts/inc/ExtMonFNALRecoCluster.hh"
#include "RecoDataProducts/inc/ExtMonFNALRecoClusterCollection.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"

namespace mu2e {

  //================================================================
  class ExtMonFNALRecoClusterization : public art::EDProducer {

  public:
    explicit ExtMonFNALRecoClusterization(fhicl::ParameterSet const& pset)
      : EDProducer{pset}
      , verbosityLevel_(pset.get<int>("verbosityLevel", 0))
      , inputModuleLabel_(pset.get<std::string>("inputModuleLabel"))
      , inputInstanceName_(pset.get<std::string>("inputInstanceName", ""))
      , geomModuleLabel_(pset.get<std::string>("geomModuleLabel"))
      , geomInstanceName_(pset.get<std::string>("geomInstanceName", ""))
      , extmon_(0)
    {
      produces<ExtMonFNALRecoClusterCollection>();
    }

    virtual void produce(art::Event& evt);
    virtual void beginRun(art::Run& run);

  private:
    int verbosityLevel_;
    std::string inputModuleLabel_;
    std::string inputInstanceName_;
    std::string geomModuleLabel_;
    std::string geomInstanceName_;

    // Non-owning pointers to the geometry object. The current Mu2e
    // infrastructure does not allow the use of a Handle as a class
    // member.
    const ExtMonFNAL::ExtMon *extmon_;

    ExtMonFNALRecoCluster computeCluster(const art::Ptr<ExtMonFNALRawCluster>& raw);

    void computeClusters(ExtMonFNALRecoClusterCollection *clusters,
                         const art::Handle<ExtMonFNALRawClusterCollection>& raw,
                         const art::EDProductGetter* rawGetter
                         );

  };

  //================================================================
  void ExtMonFNALRecoClusterization::beginRun(art::Run& run) {
    if(!geomModuleLabel_.empty()) {
      if(verbosityLevel_ > 0) {
        std::cout<<"ExtMonFNALRecoClusterization: using recorded geometry: "
                 <<"("<<geomModuleLabel_<<", "<<geomInstanceName_<<")"
                 <<std::endl;
      }
      art::Handle<ExtMonFNAL::ExtMon> emf;
      run.getByLabel(geomModuleLabel_, geomInstanceName_, emf);
      extmon_ = &*emf;
    }
    else {
      if(verbosityLevel_ > 0) {
        std::cout<<"ExtMonFNALRecoClusterization: using GeometryService"<<std::endl;
      }
      GeomHandle<ExtMonFNAL::ExtMon> emf;
      extmon_ = &*emf;
    }
  }

  //================================================================
  void ExtMonFNALRecoClusterization::produce(art::Event& event) {

    std::unique_ptr<ExtMonFNALRecoClusterCollection> reco
      (new ExtMonFNALRecoClusterCollection(extmon_->up().nplanes()+extmon_->dn().nplanes()));

    art::Handle<ExtMonFNALRawClusterCollection> raw;
    event.getByLabel(inputModuleLabel_, inputInstanceName_, raw);

    const art::EDProductGetter *rawGetter = event.productGetter(raw.id());
    computeClusters(&*reco, raw, rawGetter);

    event.put(std::move(reco));
  }

  //================================================================
  void ExtMonFNALRecoClusterization::computeClusters(ExtMonFNALRecoClusterCollection *clusters,
                                                     const art::Handle<ExtMonFNALRawClusterCollection>& rawh,
                                                     const art::EDProductGetter* rawGetter
                                                     )
  {
    for(unsigned i = 0, max=rawh->size(); i < max; ++i) {
      art::Ptr<ExtMonFNALRawCluster> raw(rawh.id(), i, rawGetter);
      clusters->insert(computeCluster(raw));
    }
  }

  //================================================================
  struct CmpPixelX {
    bool operator()(const art::Ptr<ExtMonFNALRawHit>& a, const art::Ptr<ExtMonFNALRawHit>& b) {
      return (a->pixelId().chip().chipCol() < b->pixelId().chip().chipCol()) ||
        (a->pixelId().chip().chipCol() == b->pixelId().chip().chipCol() && a->pixelId().col() < b->pixelId().col());
    }
  };

  struct CmpPixelY {
    bool operator()(const art::Ptr<ExtMonFNALRawHit>& a, const art::Ptr<ExtMonFNALRawHit>& b) {
      return (a->pixelId().chip().chipRow() < b->pixelId().chip().chipRow()) ||
        (a->pixelId().chip().chipRow() == b->pixelId().chip().chipRow() && a->pixelId().row() < b->pixelId().row());
    }
  };

  ExtMonFNALRecoCluster ExtMonFNALRecoClusterization::computeCluster(const art::Ptr<ExtMonFNALRawCluster>& raw) {
    typedef ExtMonFNALRawCluster::Hits::const_iterator Iter;

    Iter pxmax = std::max_element(raw->hits().begin(), raw->hits().end(), CmpPixelX() );
    assert(pxmax != raw->hits().end());
    Iter pxmin = std::min_element(raw->hits().begin(), raw->hits().end(), CmpPixelX() );
    assert(pxmin != raw->hits().end());
    Iter pymax = std::max_element(raw->hits().begin(), raw->hits().end(), CmpPixelY() );
    assert(pymax != raw->hits().end());
    Iter pymin = std::min_element(raw->hits().begin(), raw->hits().end(), CmpPixelY() );
    assert(pymin != raw->hits().end());

    const int xWidth =
      extmon_->chip().nColumns() * ((*pxmax)->pixelId().chip().chipCol() -
                                    (*pxmin)->pixelId().chip().chipCol())
      + (*pxmax)->pixelId().col() - (*pxmin)->pixelId().col() + 1;

    const int yWidth =
      extmon_->chip().nRows() * ((*pymax)->pixelId().chip().chipRow() -
                                 (*pymin)->pixelId().chip().chipRow())
      + (*pymax)->pixelId().row() - (*pymin)->pixelId().row() + 1;

    double sumw(0), sumClock(0);
    CLHEP::Hep3Vector sumPos;
    for(Iter i = raw->hits().begin(); i != raw->hits().end(); ++i) {
      const double weight = 1.; // could weight by a function of ToT
      sumw     += weight;
      sumClock += weight * (*i)->clock();
      sumPos   += weight * extmon_->pixelPositionInPlaneStack((*i)->pixelId());
    }

    unsigned int plane = raw->hits().front()->pixelId().chip().module().plane();
    return ExtMonFNALRecoCluster(raw, plane, sumPos/sumw, xWidth, yWidth, int(sumClock/sumw));
  }

  //================================================================

} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::ExtMonFNALRecoClusterization)
