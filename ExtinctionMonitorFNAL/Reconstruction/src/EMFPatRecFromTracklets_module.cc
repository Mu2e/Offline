// Fast search for tracks in a narrow range of angles around
// ExtMonFNAL axis.  The algorithm starts by finding straight line
// tracklets in upstream and downstream stacks separately.  A pair of
// tracklets (one up, one donwstream) measures bend angle in the
// magnet.  The upstream tracklet is assigned the corresponding
// momentum and extrapolated donwstream where the consistency with the
// donwstream tracklet is checked.
//
// $Id: EMFPatRecFromTracklets_module.cc,v 1.2 2012/09/19 03:58:29 gandr Exp $
// $Author: gandr $
// $Date: 2012/09/19 03:58:29 $
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

#include "CLHEP/Vector/TwoVector.h"
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
#include "ExtinctionMonitorFNAL/Reconstruction/inc/PixelRecoUtils.hh"

#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "TH1D.h"
#include "TH2D.h"

//#define AGDEBUG(stuff) do { std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<", func "<<__func__<<": "<<stuff<<std::endl; } while(0)
#define AGDEBUG(stuff)

namespace mu2e {
  namespace ExtMonFNAL {

    //================================================================
    struct Tracklet {

      art::Ptr<ExtMonFNALRecoCluster> firstCluster;
      art::Ptr<ExtMonFNALRecoCluster> lastCluster;

      std::vector<art::Ptr<ExtMonFNALRecoCluster> > middleClusters;

      Tracklet(const art::Ptr<ExtMonFNALRecoCluster>& fc,
               const art::Ptr<ExtMonFNALRecoCluster>& lc)
        : firstCluster(fc)
        , lastCluster(lc)
      {}
    };

    typedef std::list<Tracklet> Tracklets;

    std::ostream& operator<<(std::ostream& os, const Tracklet& tl) {
      return os<<"Tracklet(fc="<<tl.firstCluster
               <<", lc="<<tl.lastCluster
               <<", nmiddle="<<tl.middleClusters.size()
               <<" )";
    }

    //================================================================
    class EMFPatRecFromTracklets : public art::EDProducer {

    public:
      explicit EMFPatRecFromTracklets(fhicl::ParameterSet const& pset)
        : verbosityLevel_(pset.get<int>("verbosityLevel", 0))
        , inputModuleLabel_(pset.get<std::string>("inputModuleLabel"))
        , inputInstanceName_(pset.get<std::string>("inputInstanceName", ""))
        , geomModuleLabel_(pset.get<std::string>("geomModuleLabel"))
        , geomInstanceName_(pset.get<std::string>("geomInstanceName", ""))
        , extmon_(0)
        , extrapolator_(extmon_)

        , clusterClockTolerance_(pset.get<unsigned>("clusterClockTolerance"))
        , slopexmax_(pset.get<double>("maxTrackSlope"))
        , slopeymax_(slopexmax_)
        , stackScatterAngleTolerance_(pset.get<double>("stackScatterAngleTolerance"))

        , trackletMatchSlopeXTolerance_(pset.get<double>("trackletMatchSlopeXTolerance"))
        , trackletMatchXTolerance_(pset.get<double>("trackletMatchXTolerance"))
        , trackletMatchYTolerance_(pset.get<double>("trackletMatchYTolerance"))

        , alignmentToleranceX_(pset.get<double>("alignmentToleranceX"))
        , alignmentToleranceY_(pset.get<double>("alignmentToleranceY"))

        , hTrackletMultiplicity_()
        , hTrackletMatchXY_()
        , hTrackletMatchSlopeX_()
        , hClusterAddXY_()
        , hClusterAddXXMax_()
        , hClusterAddYYMax_()
        , hClockDiffTrackletSeedClusters_()
        , hClockDiffClusterTracklet_()

        , singleParticleMode_(pset.get<bool>("singleParticleMode", false))
     {
        produces<ExtMonFNALTrkParamCollection>();
        produces<ExtMonFNALPatRecTrackAssns>();

        if(singleParticleMode_) {
          std::cout<<"EMFPatRecFromTracklets: working in the single particle mode"<<std::endl;
        }
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

      int clusterClockTolerance_;

      // We look for tracks that are parallel to the detector axis within
      double slopexmax_;
      double slopeymax_;
      double stackScatterAngleTolerance_;
      double trackletMatchSlopeXTolerance_;
      double trackletMatchXTolerance_;
      double trackletMatchYTolerance_;
      double alignmentToleranceX_;
      double alignmentToleranceY_;

      //----------------------------------------------------------------
      TH2D *hTrackletMultiplicity_;
      TH2D *hTrackletMatchXY_;
      TH1D *hTrackletMatchSlopeX_;
      std::vector<TH2D*> hClusterAddXY_;
      std::vector<TH2D*> hClusterAddXXMax_;
      std::vector<TH2D*> hClusterAddYYMax_;
      TH1D *hClockDiffTrackletSeedClusters_;
      TH1D *hClockDiffClusterTracklet_;

      //----------------------------------------------------------------
      // debugging and tuning aid
      bool singleParticleMode_;

      bool acceptSingleParticleEvent(const art::Event& event);

      //----------------------------------------------------------------
      Tracklets formTracklets(const ExtMonFNALSensorStack& stack,
                              const art::Handle<ExtMonFNALRecoClusterCollection>& clusters);

      void findTracks(art::Event& event,
                      ExtMonFNALPatRecTrackAssns *tracks,
                      ExtMonFNALTrkParamCollection *params,
                      const art::Handle<ExtMonFNALRecoClusterCollection>& clusters);

      void addPlaneClusters(Tracklets* seeds,
                            const art::Handle<ExtMonFNALRecoClusterCollection>& coll,
                            const ExtMonFNALSensorStack& stack,
                            unsigned stackPlane,
                            unsigned anchorPlane1,
                            unsigned anchorPlane2
                            );

      void recordTrackClusters(ExtMonFNALPatRecTrackAssns *outAssns,
                               const art::Ptr<ExtMonFNALTrkParam>& trackPar,
                               const Tracklet& tl);

      Tracklet createTracklet(const Tracklet& orig, const art::Ptr<ExtMonFNALRecoCluster>& cl);

      void modifyInPlace(Tracklet *inout, const art::Ptr<ExtMonFNALRecoCluster>& cl);

      double yangle(const Tracklet& tl);
      double slopex(const Tracklet& tl);

      bool inTime(const Tracklet& tl, const ExtMonFNALRecoCluster& cl);
      bool inTime(const Tracklet& tl1, const Tracklet& tl2);
    };

    //================================================================
    void EMFPatRecFromTracklets::beginRun(art::Run& run) {
      if(!geomModuleLabel_.empty()) {
        if(verbosityLevel_ > 0) {
          std::cout<<"EMFPatRecFromTracklets: using recorded geometry: "
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
          std::cout<<"EMFPatRecFromTracklets: using GeometryService"<<std::endl;
        }
        GeomHandle<ExtMonFNAL::ExtMon> emf;
        extmon_ = &*emf;
        extrapolator_ = TrackExtrapolator(extmon_);
      }

      //----------------
      art::ServiceHandle<art::TFileService> tfs;
      hTrackletMultiplicity_ = tfs->make<TH2D>("trackletMultiplicity", "Tracklet multiplicity down- vs upstream",
                                               200, -0.5, 199.5, 200, -0.5, 199.5);
      hTrackletMultiplicity_->SetOption("colz");

      hTrackletMatchXY_ = tfs->make<TH2D>("trackletMatchXY", "Tracklet match (extrapolation - cluster) xy",
                                          200, -5., 5., 200, -5., 5.);

      hTrackletMatchXY_->SetOption("colz");

      hTrackletMatchSlopeX_ = tfs->make<TH1D>("trackletMatchSlopeX", "Tracklet match dslopex", 200, -15.e-3, 15.e-3);

      hClusterAddXY_.resize(extmon_->nplanes());
      hClusterAddXXMax_.resize(extmon_->nplanes());
      hClusterAddYYMax_.resize(extmon_->nplanes());
      for(unsigned plane=0; plane < extmon_->nplanes(); ++plane) {
        std::ostringstream xyname;
        xyname<<"clusterAddXY"<<plane;
        std::ostringstream xytitle;
        xytitle<<"dy/dymax vs dx/dxmax for plane "<<plane;
        hClusterAddXY_[plane] = tfs->make<TH2D>(xyname.str().c_str(), xytitle.str().c_str(), 200, -5., 5., 200, -5., 5.);
        hClusterAddXY_[plane]->SetOption("colz");

        std::ostringstream xxname;
        xxname<<"clusterAddXXMax"<<plane;
        std::ostringstream xxtitle;
        xxtitle<<"dx vs dxmax for "<<plane;
        hClusterAddXXMax_[plane] = tfs->make<TH2D>(xxname.str().c_str(), xxtitle.str().c_str(), 100, 0., 5., 200, -5., 5.);
        hClusterAddXXMax_[plane]->SetOption("colz");

        std::ostringstream yyname;
        yyname<<"clusterAddYYMax"<<plane;
        std::ostringstream yytitle;
        yytitle<<"dy vs dymax for "<<plane;
        hClusterAddYYMax_[plane] = tfs->make<TH2D>(yyname.str().c_str(), yytitle.str().c_str(), 100, 0., 5., 200, -5., 5.);
        hClusterAddYYMax_[plane]->SetOption("colz");
      }

      hClockDiffTrackletSeedClusters_ = tfs->make<TH1D>("clockDiffTrackletSeedClusters", "Tracklet seed cluster time diff", 120, -30.5, 89.5);
      hClockDiffClusterTracklet_ = tfs->make<TH1D>("clockDiffClusterTracklet", "Cluster - tracklet time clock", 120, -30.5, 89.5);
    }

    //================================================================
    void EMFPatRecFromTracklets::produce(art::Event& event) {
      art::Handle<ExtMonFNALRecoClusterCollection> clusters;
      event.getByLabel(inputModuleLabel_, inputInstanceName_, clusters);

      std::auto_ptr<ExtMonFNALTrkParamCollection> params(new ExtMonFNALTrkParamCollection);
      std::auto_ptr<ExtMonFNALPatRecTrackAssns> tracks(new ExtMonFNALPatRecTrackAssns);

      if(!singleParticleMode_ || acceptSingleParticleEvent(event)) {
        findTracks(event, &*tracks, &*params, clusters);
      }

      event.put(params);
      event.put(tracks);
    }

    //================================================================
    void EMFPatRecFromTracklets::findTracks(art::Event& event,
                                            ExtMonFNALPatRecTrackAssns *tracks,
                                            ExtMonFNALTrkParamCollection *params,
                                            const art::Handle<ExtMonFNALRecoClusterCollection>& coll)
    {
      Tracklets tup = formTracklets(extmon_->up(), coll);
      Tracklets tdn = formTracklets(extmon_->dn(), coll);
      hTrackletMultiplicity_->Fill(tdn.size(), tup.size());

      // merge compatible tracklet pairs into tracks
      const double nominalBendHalfAngle = extmon_->spectrometerMagnet().nominalBendHalfAngle();
      const double twoOverL(1./extmon_->spectrometerMagnet().outerHalfSize()[2]);
      for(Tracklets::const_iterator iup = tup.begin(); iup != tup.end(); ++iup) {
        for(Tracklets::const_iterator idn = tdn.begin(); idn != tdn.end(); ++idn) {

          const double dslopex = slopex(*iup) - slopex(*idn);
          hTrackletMatchSlopeX_->Fill(dslopex);

          if(inTime(*iup, *idn) && (std::abs(dslopex) < trackletMatchSlopeXTolerance_)) {
            const double halfBend = 0.5*(yangle(*idn) - yangle(*iup)) + nominalBendHalfAngle;
            const double rinv = twoOverL * sin(halfBend);
            const CLHEP::Hep3Vector& startPos = iup->firstCluster->position();
            const CLHEP::Hep3Vector startDir(iup->lastCluster->position() - iup->firstCluster->position());
            ExtMonFNALTrkParam testpar;
            testpar.setz0(startPos.z());
            testpar.setposx(startPos.x());
            testpar.setposy(startPos.y());
            testpar.setslopex(startDir.x()/startDir.z());
            testpar.setslopey(startDir.y()/startDir.z());
            testpar.setrinv(rinv);

            const ExtMonFNALRecoCluster& targetCluster = *idn->lastCluster;
            if(extrapolator_.extrapolateToPlane(targetCluster.plane(), &testpar)) {
              AGDEBUG("targetCluster = "<<targetCluster<<", extrapolation result = "<<testpar);

              const double missx = testpar.posx() - targetCluster.position().x();
              const double missy = testpar.posy() - targetCluster.position().y();
              hTrackletMatchXY_->Fill(missx, missy);

              if((std::abs(missx) < trackletMatchXTolerance_) && (std::abs(missy) < trackletMatchYTolerance_)) {
                // Accept the tracklet pair.
                // Will store track parameter estimate at the signal particle entrance - that is, the last
                // plane of the upstream stack.

                const CLHEP::Hep3Vector& trackPos(iup->lastCluster->position());
                const CLHEP::Hep3Vector& trackDir(startDir);
                ExtMonFNALTrkParam trackPar;
                trackPar.setz0(trackPos.z());
                trackPar.setposx(trackPos.x());
                trackPar.setposy(trackPos.y());
                trackPar.setslopex(trackDir.x()/trackDir.z());
                trackPar.setslopey(trackDir.y()/trackDir.z());
                trackPar.setrinv(testpar.rinv());

                // Add to the output collection
                params->push_back(trackPar);

                // Ptr to the newly added track parameters
                const art::ProductID paramsPID = getProductID<ExtMonFNALTrkParamCollection>(event);
                const art::EDProductGetter *paramsGetter = event.productGetter(paramsPID);
                art::Ptr<ExtMonFNALTrkParam> ppar(paramsPID, params->size()-1, paramsGetter);

                // Associate clusters to the track
                recordTrackClusters(tracks, ppar, *idn);
                recordTrackClusters(tracks, ppar, *iup);
              }
            }
          }
        }
      }

    }

    //================================================================
    Tracklets EMFPatRecFromTracklets::formTracklets(const ExtMonFNALSensorStack& stack,
                                                    const art::Handle<ExtMonFNALRecoClusterCollection>& coll)
    {
      Tracklets res;

      // Start with pair combinations of clusters on the outside stack planes
      const unsigned plane1 = stack.planeNumberOffset() + 0;
      const unsigned plane2 = stack.planeNumberOffset() + stack.nplanes() - 1;
      AGDEBUG("plane1 = "<<plane1<<", plane2 = "<<plane2);

      typedef ExtMonFNALRecoClusterCollection::PlaneClusters PC;
      const PC pc1 = coll->clusters(plane1);
      const PC pc2 = coll->clusters(plane2);

      const double dz = stack.sensor_zoffset().back() - stack.sensor_zoffset().front();
      const double dxmax = std::abs(dz * slopexmax_);
      const double dymax = std::abs(dz * slopeymax_);

      for(unsigned i1 = 0; i1 < pc1.size(); ++i1) {
        art::Ptr<ExtMonFNALRecoCluster> c1(coll, coll->globalIndex(plane1, i1));

        for(unsigned i2 = 0; i2 < pc2.size(); ++i2) {
          art::Ptr<ExtMonFNALRecoCluster> c2(coll, coll->globalIndex(plane2, i2));

          hClockDiffTrackletSeedClusters_->Fill(c2->clock() - c1->clock());

          if( (std::abs(c1->clock() - c2->clock()) <= clusterClockTolerance_) &&
              (std::abs(c1->position().x() - c2->position().x()) < dxmax) &&
              (std::abs(c1->position().y() - c2->position().y()) < dymax) )
            {
              // filter out or add candidates using information from middle planes
              Tracklets tmp;
              tmp.insert(tmp.begin(), Tracklet(c1,c2));
              for(unsigned stackPlane=0; stackPlane < stack.nplanes(); ++stackPlane) {
                if(stackPlane == plane1 - stack.planeNumberOffset()) continue;
                if(stackPlane == plane2 - stack.planeNumberOffset()) continue;
                addPlaneClusters(&tmp, coll, stack, stackPlane,
                                 plane1 -stack.planeNumberOffset(),
                                 plane2 - stack.planeNumberOffset());
              }

              res.insert(res.begin(), tmp.begin(), tmp.end());
            }
        }
      }

      return res;
    }

    //================================================================
    void EMFPatRecFromTracklets::addPlaneClusters(Tracklets* seeds,
                                                  const art::Handle<ExtMonFNALRecoClusterCollection>& coll,
                                                  const ExtMonFNALSensorStack& stack,
                                                  unsigned stackPlane,
                                                  unsigned anchorPlane1,
                                                  unsigned anchorPlane2
                                                  )
    {
      const double planeZ = stack.sensor_zoffset()[stackPlane];
      AGDEBUG("Here: stackPlane="<<stackPlane<<", planeZ = "<<planeZ<<", seeds->size() = "<<seeds->size()
              <<", anchorPlane1 = "<<anchorPlane1<<", anchorPlane2 = "<<anchorPlane2);

      const double dxtol1 =  alignmentToleranceX_ + stackScatterAngleTolerance_ * std::abs(planeZ - stack.sensor_zoffset()[anchorPlane1]);
      const double dxtol2 =  alignmentToleranceX_ + stackScatterAngleTolerance_ * std::abs(planeZ - stack.sensor_zoffset()[anchorPlane2]);

      const double dytol1 =  alignmentToleranceY_ + stackScatterAngleTolerance_ * std::abs(planeZ - stack.sensor_zoffset()[anchorPlane1]);
      const double dytol2 =  alignmentToleranceY_ + stackScatterAngleTolerance_ * std::abs(planeZ - stack.sensor_zoffset()[anchorPlane2]);
      AGDEBUG("dxtol1 ="<<dxtol1<<", dxtol2 ="<<dxtol2<<", dytol1 = "<<dytol1<<", dytol2 = "<<dytol2);

      for(Tracklets::iterator i = seeds->begin(); i != seeds->end(); ) {

        const CLHEP::Hep3Vector& pos1 = i->firstCluster->position();
        const CLHEP::Hep3Vector& pos2 = i->lastCluster->position();
        AGDEBUG("pos1 = "<<pos1<<", pos2 = "<<pos2);

        const double z1(pos1.z()), z2(pos2.z());
        const double dz(z2-z1);

        const CLHEP::Hep3Vector interpolated = pos2 * (planeZ - z1)/dz + pos1 * (z2 - planeZ)/dz;

        const double dxmax1 = dxtol1 + 0.5*i->firstCluster->xWidth()*extmon_->chip().xPitch();
        const double dxmax2 = dxtol2 + 0.5*i->lastCluster->xWidth()*extmon_->chip().xPitch();
        // satisfy both constraints
        const double dxmax(std::min(dxmax1, dxmax2));

        const double dymax1 = dytol1 + 0.5*i->firstCluster->yWidth()*extmon_->chip().yPitch();
        const double dymax2 = dytol2 + 0.5*i->lastCluster->yWidth()*extmon_->chip().yPitch();
        // satisfy both constraints
        const double dymax(std::min(dymax1, dymax2));
        AGDEBUG("dxmax="<<dxmax<<", dymax="<<dymax<<", interpolated="<<interpolated);

        // find all clusters compatible with seed i
        std::vector<art::Ptr<ExtMonFNALRecoCluster> > compatibleClusters;
        const unsigned int globalPlane = stackPlane + stack.planeNumberOffset();
        const ExtMonFNALRecoClusterCollection::PlaneClusters& clusters = coll->clusters(globalPlane);
        for(unsigned ic=0; ic<clusters.size(); ++ic) {
          const ExtMonFNALRecoCluster& cl = clusters[ic];
          if(inTime(*i, cl)) {
            const double dx = interpolated.x() - cl.position().x();
            const double dy = interpolated.y() - cl.position().y();

            hClusterAddXY_[globalPlane]->Fill(dx/dxmax, dy/dymax);
            hClusterAddXXMax_[globalPlane]->Fill(dxmax, dx);
            hClusterAddYYMax_[globalPlane]->Fill(dymax, dy);

            if((std::abs(dx) < dxmax) && (std::abs(dy) < dymax) ) {
              compatibleClusters.push_back(art::Ptr<ExtMonFNALRecoCluster>(coll, coll->globalIndex(globalPlane, ic)));
            }

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
            seeds->insert(i, createTracklet(*i, compatibleClusters[count]));
          }
          // One of the clusters can be added to the current seed
          // in place
          modifyInPlace(&*i, compatibleClusters[0]);

          ++i;
        } // One or more clusters added

      } // for(tracklets)
    }

    //================================================================
    Tracklet EMFPatRecFromTracklets::createTracklet(const Tracklet& orig, const art::Ptr<ExtMonFNALRecoCluster>& cl) {
      Tracklet res(orig);
      modifyInPlace(&res, cl);
      return res;
    }

    //================================================================
    void EMFPatRecFromTracklets::modifyInPlace(Tracklet *inout, const art::Ptr<ExtMonFNALRecoCluster>& cl) {
      inout->middleClusters.push_back(cl);
    }

    //================================================================
    double EMFPatRecFromTracklets::yangle(const Tracklet& tl) {
      const CLHEP::Hep3Vector dir(tl.lastCluster->position() - tl.firstCluster->position());
      return atan(dir.y()/dir.z());
    }

    //================================================================
    double EMFPatRecFromTracklets::slopex(const Tracklet& tl) {
      const CLHEP::Hep3Vector dir(tl.lastCluster->position() - tl.firstCluster->position());
      return dir.x()/dir.z();
    }

    //================================================================
    bool EMFPatRecFromTracklets::inTime(const Tracklet& tl, const ExtMonFNALRecoCluster& cl) {
      const int dtFirst = cl.clock() - tl.firstCluster->clock();
      const int dtLast  = cl.clock() - tl.lastCluster->clock();

      hClockDiffClusterTracklet_->Fill(dtFirst);
      hClockDiffClusterTracklet_->Fill(dtLast);

      return
        (std::abs(dtFirst) <= clusterClockTolerance_) &&
        (std::abs(dtLast) <= clusterClockTolerance_);
    }

    //================================================================
    bool EMFPatRecFromTracklets::inTime(const Tracklet& tl1, const Tracklet& tl2) {
      return inTime(tl1, *tl2.firstCluster) && inTime(tl1, *tl2.lastCluster);
    }

    //================================================================
    void EMFPatRecFromTracklets::recordTrackClusters(ExtMonFNALPatRecTrackAssns *outAssns,
                                                     const art::Ptr<ExtMonFNALTrkParam>& trackPar,
                                                     const Tracklet& tl)
    {
      outAssns->addSingle(tl.firstCluster, trackPar);
      for(unsigned i=0; i<tl.middleClusters.size(); ++i) {
        outAssns->addSingle(tl.middleClusters[i], trackPar);
      }
      outAssns->addSingle(tl.lastCluster, trackPar);
    }

    //================================================================
    bool EMFPatRecFromTracklets::acceptSingleParticleEvent(const art::Event& event) {
      art::Handle<ExtMonFNALRecoClusterCollection> coll;
      event.getByLabel(inputModuleLabel_, inputInstanceName_, coll);
      return perfectSingleParticleEvent(*coll, extmon_->nplanes());
    }

    //================================================================

  } // end namespace ExtMonFNAL
} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::ExtMonFNAL::EMFPatRecFromTracklets)
