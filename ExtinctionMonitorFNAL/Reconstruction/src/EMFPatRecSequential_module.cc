// Fast search for tracks in a narrow range of angles around ExtMonFNAL axis.
// The algorithm requires that track has hits on all pixel planes.
//
// $Id: EMFPatRecSequential_module.cc,v 1.2 2012/11/01 23:36:27 gandr Exp $
// $Author: gandr $
// $Date: 2012/11/01 23:36:27 $
//
// Original author Andrei Gaponenko
//


#include <string>
#include <iostream>
#include <cmath>
#include <memory>
#include <iterator>
#include <algorithm>
#include <limits>
#include <cassert>

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "cetlib/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/Assns.h"

#include "RecoDataProducts/inc/ExtMonFNALRecoCluster.hh"
#include "RecoDataProducts/inc/ExtMonFNALRecoClusterCollection.hh"
#include "RecoDataProducts/inc/ExtMonFNALTrkParam.hh"
#include "RecoDataProducts/inc/ExtMonFNALTrkParamCollection.hh"
#include "RecoDataProducts/inc/ExtMonFNALPatRecTrackAssns.hh"


#include "GeometryService/inc/GeomHandle.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "ExtinctionMonitorFNAL/Reconstruction/inc/TrackExtrapolator.hh"

#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "TH1D.h"
#include "TH2D.h"

//#define AGDEBUG(stuff) do { std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<", func "<<__func__<<": "<<stuff<<std::endl; } while(0)
#define AGDEBUG(stuff)

namespace mu2e {
  namespace ExtMonFNAL {

    //================================================================
    struct TrackSeed {

      const ExtMonFNALRecoClusterCollection *coll; // non-owning

      std::vector<unsigned> ihits; // global indexes into ExtMonFNALRecoClusterCollection, one entry per plane

      // pars at the extrapolation head (the last hit)
      float xmin;
      float xmax;
      float ymin;
      float ymax;
      float slopexmin;
      float slopexmax;
      float slopeymin;
      float slopeymax;
      float rinvmin;
      float rinvmax;

      TrackSeed(const ExtMonFNALRecoClusterCollection& c,
                unsigned icluster,
                float posx,
                float eposx,
                float posy,
                float eposy,
                float slopex,
                float eslopex,
                float slopey,
                float eslopey,
                float rinv1,
                float rinv2
                )
        : coll(&c)
        , xmin(posx - eposx), xmax(posx + eposx)
        , ymin(posy - eposy), ymax(posy + eposy)
        , slopexmin(slopex - eslopex), slopexmax(slopex + eslopex)
        , slopeymin(slopey - eslopey), slopeymax(slopey + eslopey)
        , rinvmin(rinv1), rinvmax(rinv2)
      {
        ihits.push_back(icluster);
      }
    };

    typedef std::list<TrackSeed> TrackSeeds;

    std::ostream& operator<<(std::ostream& os, const TrackSeed& ts) {
      return os<<"TrackSeed(nHits="<<ts.ihits.size()
               <<", xmin="<<ts.xmin
               <<", xmax="<<ts.xmax
               <<", ymin="<<ts.ymin
               <<", ymax="<<ts.ymax
               <<", sxmin="<<ts.slopexmin
               <<", sxmax="<<ts.slopexmax
               <<", symin="<<ts.slopeymin
               <<", symax="<<ts.slopeymax
               <<", rinvmin="<<ts.rinvmin
               <<", rinvmax="<<ts.rinvmax
               <<" )";
    }

    struct ExtrapolationWindow {
      int clockmin;
      int clockmax;
      float xmin;
      float xmax;
      float ymin;
      float ymax;

      // a hack: store here value from track extrapolation trhough magnet
      // that are not needed for cluster-seed compatibility evaluation
      // but are useful later
      float symin; // slope in y
      float symax;

      ExtrapolationWindow(int c1, int c2, float x1, float x2, float y1, float y2, float sy1, float sy2)
        : clockmin(c1), clockmax(c2), xmin(x1), xmax(x2), ymin(y1), ymax(y2), symin(sy1), symax(sy2)
      {}
    };

    std::ostream& operator<<(std::ostream& os, const ExtrapolationWindow& win) {
      return os<<"ExtrapolationWindow(cmin="<<win.clockmin
               <<", cmax="<<win.clockmax
               <<", xmin="<<win.xmin
               <<", xmax="<<win.xmax
               <<", ymin="<<win.ymin
               <<", ymax="<<win.ymax
               <<", symin="<<win.symin
               <<", symax="<<win.symax
               <<" )";
    }

    //================================================================
    class EMFPatRecSequential : public art::EDProducer {

    public:
      explicit EMFPatRecSequential(fhicl::ParameterSet const& pset)
        : verbosityLevel_(pset.get<int>("verbosityLevel", 0))
        , inputModuleLabel_(pset.get<std::string>("inputModuleLabel"))
        , inputInstanceName_(pset.get<std::string>("inputInstanceName", ""))
        , geomModuleLabel_(pset.get<std::string>("geomModuleLabel"))
        , geomInstanceName_(pset.get<std::string>("geomInstanceName", ""))
        , extmon_(0)
        , extrapolator_(extmon_)
        , slopexmax_(pset.get<float>("maxTrackSlope"))
        , slopeymax_(slopexmax_)
        , pmin_(pset.get<float>("pmin"))
        , pmax_(pset.get<float>("pmax"))
        , clusterClockTolerance_(pset.get<unsigned>("clusterClockTolerance"))
        , perPlaneScatterTolerance_(pset.get<float>("perPlaneScatterTolerance"))
        , alignmentToleranceX_(pset.get<float>("alignmentToleranceX"))
        , alignmentToleranceY_(pset.get<float>("alignmentToleranceY"))
        , hSeedMultiplicity_(0)
      {
        produces<ExtMonFNALTrkParamCollection>();
        produces<ExtMonFNALPatRecTrackAssns>();
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
      const ExtMon *extmon_;
      TrackExtrapolator extrapolator_;

      //----------------------------------------------------------------
      // Track search parameters

      // We look for tracks that are parallel to the detector axis within
      float slopexmax_;
      float slopeymax_;
      // and have momenta within the range
      float pmin_;
      float pmax_;

      int clusterClockTolerance_;
      float perPlaneScatterTolerance_;
      float alignmentToleranceX_;
      float alignmentToleranceY_;

      //----------------------------------------------------------------
      TH2D *hSeedMultiplicity_;

      //----------------------------------------------------------------
      void findTracks(ExtMonFNALPatRecTrackAssns *tracks,
                      ExtMonFNALTrkParamCollection *params,
                      const art::ProductID& paramsPID,
                      const art::EDProductGetter *paramsGetter,
                      const art::Handle<ExtMonFNALRecoClusterCollection>& clusters);

      void addPlaneClusters(TrackSeeds* seeds,
                            const ExtMonFNALRecoClusterCollection& coll,
                            unsigned plane);

      ExtrapolationWindow extrapolate(const TrackSeed& ts, unsigned toPlane);

      bool compatible(const ExtrapolationWindow& win, const ExtMonFNALRecoCluster& cl);

      TrackSeed createTrackSeed(const TrackSeed& parent,
                                unsigned nextClusterIndex,
                                const ExtrapolationWindow& win);

      void modifyInPlace(TrackSeed *ts,
                         unsigned nextClusterIndex,
                         const ExtrapolationWindow& win);

      bool sameSensorStack(unsigned plane1, unsigned plane2) {
        bool dn1 = (plane1 < extmon_->dn().nplanes());
        bool dn2 = (plane2 < extmon_->dn().nplanes());
        return !(dn1^dn2);
      }
    };

    //================================================================
    void EMFPatRecSequential::beginRun(art::Run& run) {
      // This is a workaround for geometry not being available at beginJob()
      if(!extmon_) {
        if(!geomModuleLabel_.empty()) {
          if(verbosityLevel_ > 0) {
            std::cout<<"EMFPatRecSequential: using recorded geometry: "
                     <<"("<<geomModuleLabel_<<", "<<geomInstanceName_<<")"
                     <<std::endl;
          }
          art::Handle<ExtMonFNAL::ExtMon> emf;
          run.getByLabel(geomModuleLabel_, geomInstanceName_, emf);
          extmon_ = &*emf;
          extrapolator_ = TrackExtrapolator(extmon_);
        }
        else {
          if(verbosityLevel_ > 0) {
            std::cout<<"EMFPatRecSequential: using GeometryService"<<std::endl;
          }
          GeomHandle<ExtMonFNAL::ExtMon> emf;
          extmon_ = &*emf;
          extrapolator_ = TrackExtrapolator(extmon_);
        }

        //----------------
        const unsigned nplanes = extmon_->up().nplanes() + extmon_->dn().nplanes();
        art::ServiceHandle<art::TFileService> tfs;
        hSeedMultiplicity_ = tfs->make<TH2D>("seedMultiplicity", "Plane vs track seed multiplicity at the plane",
                                             200, -0.5, 999.5, nplanes, -0.5, nplanes-0.5);

        hSeedMultiplicity_->SetOption("colz");
        hSeedMultiplicity_->GetXaxis()->SetTitle("track seed multiplicity");
        hSeedMultiplicity_->GetYaxis()->SetTitle("plane number");
      }
    }

    //================================================================
    void EMFPatRecSequential::produce(art::Event& event) {
      art::Handle<ExtMonFNALRecoClusterCollection> clusters;
      event.getByLabel(inputModuleLabel_, inputInstanceName_, clusters);

      std::auto_ptr<ExtMonFNALTrkParamCollection> params(new ExtMonFNALTrkParamCollection);
      std::auto_ptr<ExtMonFNALPatRecTrackAssns> tracks(new ExtMonFNALPatRecTrackAssns);

      const art::ProductID paramsPID = getProductID<ExtMonFNALTrkParamCollection>(event);
      const art::EDProductGetter *paramsGetter = event.productGetter(paramsPID);
      findTracks(&*tracks, &*params, paramsPID, paramsGetter, clusters);

      event.put(params);
      event.put(tracks);
    }

    //================================================================
    void EMFPatRecSequential::findTracks(ExtMonFNALPatRecTrackAssns *tracks,
                                         ExtMonFNALTrkParamCollection *params,
                                         const art::ProductID& paramsPID,
                                         const art::EDProductGetter *paramsGetter,
                                         const art::Handle<ExtMonFNALRecoClusterCollection>& coll)
    {
      TrackSeeds seeds;

      const unsigned lastPlane = extmon_->dn().nplanes() + extmon_->up().nplanes() - 1;
      typedef ExtMonFNALRecoClusterCollection::PlaneClusters PlaneClusters;

      const PlaneClusters& plane = coll->clusters(lastPlane);
      for(unsigned i=0; i<plane.size(); ++i) {
        const ExtMonFNALRecoCluster cl = plane[i];

        const float ex = 0.5*cl.xWidth()*extmon_->chip().xPitch();
        const float ey = 0.5*cl.yWidth()*extmon_->chip().yPitch();

        const ExtMonFNALMagnet& mag = extmon_->spectrometerMagnet();
        const float rinvmin = 1./mag.trackBendRadius(pmax_);
        const float rinvmax = 1./mag.trackBendRadius(pmin_);

        seeds.insert(seeds.begin(),
                     TrackSeed(*coll,
                               coll->globalIndex(lastPlane, i),
                               cl.position().x(),
                               ex,
                               cl.position().y(),
                               ey,
                               0., // parallel to the axis
                               slopexmax_,
                               0., // parallel to the axis
                               slopeymax_,
                               rinvmin,
                               rinvmax
                               ));
      }

      hSeedMultiplicity_->Fill(seeds.size(), lastPlane);
      AGDEBUG("Got "<<seeds.size()<<" track seeds from the last plane");

      for(unsigned plane=lastPlane; plane > 0; ) {
        --plane;

        AGDEBUG("Adding clusters from plane "<<plane);
        addPlaneClusters(&seeds, *coll, plane);
        hSeedMultiplicity_->Fill(seeds.size(), plane);
      }

      AGDEBUG("Here: seeds.size() = "<<seeds.size());

      for(TrackSeeds::const_iterator iseed=seeds.begin(); iseed != seeds.end(); ++iseed) {
        ExtMonFNALTrkParam par;
        const ExtMonFNALRecoCluster& lastCluster = (*iseed->coll)[iseed->ihits.back()];
        par.setz0(lastCluster.position().z());
        par.setposx(0.5*(iseed->xmin + iseed->xmax));
        par.setslopex(0.5*(iseed->slopexmin + iseed->slopexmax));
        par.setposy(0.5*(iseed->ymin + iseed->ymax));
        par.setslopey(0.5*(iseed->slopeymin + iseed->slopeymax));
        par.setrinv(iseed->rinvmin);

        params->push_back(par);
        art::Ptr<ExtMonFNALTrkParam> ppar(paramsPID, params->size()-1, paramsGetter);
        for(unsigned ihit=0; ihit < iseed->ihits.size(); ++ihit) {
          art::Ptr<ExtMonFNALRecoCluster> phit(coll, iseed->ihits[ihit]);
          tracks->addSingle(phit, ppar);
        }
      }

    }

    //================================================================
    void EMFPatRecSequential::addPlaneClusters(TrackSeeds* seeds,
                                               const ExtMonFNALRecoClusterCollection& coll,
                                               unsigned plane
                                               )
    {
      AGDEBUG("Here");

      const ExtMonFNALRecoClusterCollection::PlaneClusters& clusters = coll.clusters(plane);

      AGDEBUG("Here: seeds->size() = "<<seeds->size());

      for(TrackSeeds::iterator i = seeds->begin(); i != seeds->end(); ) {

        AGDEBUG("Here");

        // find all clusters compatible with seed i
        std::vector<unsigned> compatibleClusters;
        ExtrapolationWindow win(extrapolate(*i, plane));
        AGDEBUG(win);

        for(unsigned ic=0; ic<clusters.size(); ++ic) {
          if(compatible(win, clusters[ic])) {
            compatibleClusters.push_back(coll.globalIndex(plane, ic));
          }
        }

        AGDEBUG("compatibleClusters.size() = "<<compatibleClusters.size());

        if(compatibleClusters.empty()) {
          // We don't allow missing space points.  Kill the seed
          seeds->erase(i++);
        }
        else {
          // loop for the case numCompatibleClusters > 1
          for(unsigned count = 1; count < compatibleClusters.size(); ++count) {
            // insert BEFORE the current point
            seeds->insert(i, createTrackSeed(*i, compatibleClusters[count], win));
          }
          // One of the clusters can be added to the current seed
          // in place
          modifyInPlace(&*i, compatibleClusters[0], win);

          ++i;
        }

      }
    }

    //================================================================
    ExtrapolationWindow EMFPatRecSequential::extrapolate(const TrackSeed& ts, unsigned toPlane)
    {
      const ExtMonFNALRecoCluster& lastCluster = (*ts.coll)[ts.ihits.back()];

      // Times are always compared to the time of the initial hit so that
      // there is no "drift"
      int clock0 = (*ts.coll)[ts.ihits[0]].clock();

      float xmin(0.), xmax(0.), ymin(0.), ymax(0.);

      int clockmin = clock0 - clusterClockTolerance_;
      int clockmax = clock0 + clusterClockTolerance_;

      float symin(std::numeric_limits<float>::max());
      float symax(std::numeric_limits<float>::min());

      AGDEBUG("Input seed: "<<ts);

      if(sameSensorStack(lastCluster.plane(), toPlane)) {

        const bool tgtDN(toPlane < extmon_->dn().nplanes());
        const ExtMonFNALSensorStack& stack = tgtDN ? extmon_->dn() : extmon_->up();
        const unsigned stackPlane = tgtDN ? toPlane : toPlane - extmon_->dn().nplanes();
        const float toZ = stack.sensor_zoffset()[stackPlane];

        AGDEBUG("toZ="<<toZ<<", z0="<<lastCluster.position().z());

        // we go from larger to smaller z, therefore slope min/max are reversed
        xmin = ts.xmin + ts.slopexmax * (toZ - lastCluster.position().z());
        xmax = ts.xmax + ts.slopexmin * (toZ - lastCluster.position().z());
        ymin = ts.ymin + ts.slopeymax * (toZ - lastCluster.position().z());
        ymax = ts.ymax + ts.slopeymin * (toZ - lastCluster.position().z());
        symin = ts.slopeymin;
        symax = ts.slopeymax;
      }
      else { // account for the magnet and coordinate system change

        std::vector<ExtMonFNALTrkParam> tests;
        tests.reserve(4);

        // xpos and ypos can be added after the extrapolation, but
        // input slopey and rinv affect the result in a non-trivial way.
        // try different combinations of input min/max values to figure out
        // target min/max.  FIXME: this can be done more efficiently.
        {
          ExtMonFNALTrkParam test1;
          test1.setz0(lastCluster.position().z());
          test1.setslopex(1.); // measure the path length
          test1.setslopey(ts.slopeymin);
          test1.setrinv(ts.rinvmin);
          tests.push_back(test1);
        }

        {
          ExtMonFNALTrkParam test2;
          test2.setz0(lastCluster.position().z());
          test2.setslopex(1.); // measure the path length
          test2.setslopey(ts.slopeymax);
          test2.setrinv(ts.rinvmin);
          tests.push_back(test2);
        }

        {
          ExtMonFNALTrkParam test3;
          test3.setz0(lastCluster.position().z());
          test3.setslopex(1.); // measure the path length
          test3.setslopey(ts.slopeymax);
          test3.setrinv(ts.rinvmax);
          tests.push_back(test3);
        }

        {
          ExtMonFNALTrkParam test4;
          test4.setz0(lastCluster.position().z());
          test4.setslopex(1.); // measure the path length
          test4.setslopey(ts.slopeymin);
          test4.setrinv(ts.rinvmax);
          tests.push_back(test4);
        }

        float testymin(std::numeric_limits<float>::max());
        float testymax(std::numeric_limits<float>::min());
        float pathmin(std::numeric_limits<float>::max());
        float pathmax(std::numeric_limits<float>::min());
        bool anyPassed(false);
        for(unsigned i=0; i<tests.size(); ++i) {
          bool res = extrapolator_.extrapolateToPlane(toPlane, &tests[i]);
          anyPassed |= res;
          if(res) {
            testymin = std::min<float>(testymin, tests[i].posy());
            testymax = std::max<float>(testymax, tests[i].posy());
            pathmin = std::min<float>(pathmin, std::abs(tests[i].posx()));
            pathmax = std::max<float>(pathmax, std::abs(tests[i].posx()));
            symin = std::min<float>(symin, tests[i].slopey());
            symax = std::max<float>(symax, tests[i].slopey());
          }
          else {
            // a track in range may hit the bottom of the plane
            testymin = extmon_->dn().sensor_yoffset().back() - extmon_->sensor().halfSize()[1] - ts.ymin;
            // a rough extimate
            pathmax =  extmon_->up().sensor_zoffset()[0] - extmon_->dn().sensor_zoffset().back()
              + (M_PI - 2.)*extmon_->spectrometerMagnet().outerHalfSize()[2];

            symax = std::numeric_limits<float>::max();
          }
        }
        AGDEBUG("anyPassed = "<<anyPassed<<", testymin = "<<testymin<<", testymax = "<<testymax
                <<", pathmin = "<<pathmin<<", pathmax = "<<pathmax);

        if(anyPassed) {
          ymin = ts.ymin + testymin;
          ymax = ts.ymax + testymax;
          xmin = ts.xmin - (ts.slopexmax > 0 ? ts.slopexmax*pathmax : ts.slopexmax * pathmin);
          xmax = ts.xmax - (ts.slopexmin < 0 ? ts.slopexmin*pathmax : ts.slopexmin * pathmin);
        }
        else { // No tracks in the current range would make it through
          clockmax = -1;
          clockmin =  0;
        }
      }

      return ExtrapolationWindow(clockmin, clockmax, xmin, xmax, ymin, ymax, symin, symax);
    }

    //================================================================
    bool EMFPatRecSequential::compatible(const ExtrapolationWindow& win, const ExtMonFNALRecoCluster& cl) {
      // AGDEBUG("clock: "<<win.clockmin<<", "<<cl.clock()<<", "<<win.clockmax);
      // AGDEBUG("x:     "<<win.xmin<<", "<<cl.position().x()<<", "<<win.xmax);
      // AGDEBUG("y:     "<<win.ymin<<", "<<cl.position().y()<<", "<<win.ymax);

      return
        (win.clockmin <= cl.clock()) && (cl.clock() <= win.clockmax) &&
        (win.xmin <= cl.position().x()) && (cl.position().x() <= win.xmax) &&
        (win.ymin <= cl.position().y()) && (cl.position().y() <= win.ymax)
        ;
    }

    //================================================================
    TrackSeed EMFPatRecSequential::createTrackSeed(const TrackSeed& parent,
                                                   unsigned nextClusterIndex,
                                                   const ExtrapolationWindow& win)
    {
      TrackSeed res(parent);
      modifyInPlace(&res, nextClusterIndex, win);
      return res;
    }

    //================================================================
    void EMFPatRecSequential::modifyInPlace(TrackSeed *ts,
                                            unsigned nextClusterIndex,
                                            const ExtrapolationWindow& win)
    {
      // Add the new cluster
      // and propagate parameters to the new plane

      const ExtMonFNALRecoCluster& lastCluster = (*ts->coll)[ts->ihits.back()];
      const ExtMonFNALRecoCluster& newCluster = (*ts->coll)[nextClusterIndex];

      ts->ihits.push_back(nextClusterIndex);

      const float newToleranceX = alignmentToleranceX_ + 0.5*newCluster.xWidth()*extmon_->chip().xPitch();
      const float newToleranceY = alignmentToleranceY_ + 0.5*newCluster.yWidth()*extmon_->chip().yPitch();

      const float newxmin = newCluster.position().x() - newToleranceX;
      const float newxmax = newCluster.position().x() + newToleranceX;
      const float newymin = newCluster.position().y() - newToleranceY;
      const float newymax = newCluster.position().y() + newToleranceY;

      if(sameSensorStack(lastCluster.plane(), newCluster.plane())) {

        const float dz = newCluster.position().z() - lastCluster.position().z();

        // compute the slope as measured from the last two clusters
        // note that dz<0 changes the meaning of min/max
        const float measslopexmax = (newxmin - ts->xmax)/dz;
        const float measslopexmin = (newxmax - ts->xmin)/dz;
        const float measslopeymax = (newymin - ts->ymax)/dz;
        const float measslopeymin = (newymax - ts->ymin)/dz;

        // pick the more constraining of the new and pre-existing slope measurements
        // account for scattering in the new plane
        // and limit to the max search slope
        const float newslopexmin = std::max(-slopexmax_, std::max(measslopexmin, ts->slopexmin) - perPlaneScatterTolerance_);
        const float newslopexmax = std::min(+slopexmax_, std::min(measslopexmax, ts->slopexmax) + perPlaneScatterTolerance_);
        const float newslopeymin = std::max(-slopeymax_, std::max(measslopeymin, ts->slopeymin) - perPlaneScatterTolerance_);
        const float newslopeymax = std::min(+slopeymax_, std::min(measslopeymax, ts->slopeymax) + perPlaneScatterTolerance_);

        ts->slopexmin = newslopexmin;
        ts->slopexmax = newslopexmax;
        ts->slopeymin = newslopeymin;
        ts->slopeymax = newslopeymax;
      }
      else {
        // Don't try to improve slope in x at this stage
        // Slopey needs fixing because of track bending in the magnet
        ts->slopeymin = win.symin;
        ts->slopeymax = win.symax;
      }

      ts->xmin = newxmin;
      ts->xmax = newxmax;
      ts->ymin = newymin;
      ts->ymax = newymax;
    }

    //================================================================

  } // end namespace ExtMonFNAL
} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::ExtMonFNAL::EMFPatRecSequential)
