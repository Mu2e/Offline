// Fast search for tracks in a narrow range of angles around
// ExtMonFNAL axis.  The algorithm starts by finding straight line
// tracklets in upstream and downstream stacks separately.  A pair of
// tracklets (one up, one donwstream) measures bend angle in the
// magnet.  The upstream tracklet is assigned the corresponding
// momentum and extrapolated donwstream where the consistency with the
// donwstream tracklet is checked.
//
// $Id: EMFPatRecFromTracklets_module.cc,v 1.9 2013/03/15 15:52:04 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/15 15:52:04 $
//
// Original author Andrei Gaponenko
//


#include "ExtinctionMonitorFNAL/Reconstruction/inc/Tracklet.hh"
#include <string>
#include <iostream>
#include <cmath>
#include <memory>
#include <iterator>
#include <algorithm>
#include <limits>
#include <cassert>

#include "boost/noncopyable.hpp"

#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/GenericFunctions/CumulativeChiSquare.hh"

#include "cetlib/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"

#include "RecoDataProducts/inc/ExtMonFNALRecoCluster.hh"
#include "RecoDataProducts/inc/ExtMonFNALRecoClusterCollection.hh"
#include "RecoDataProducts/inc/ExtMonFNALTrkParam.hh"
#include "RecoDataProducts/inc/ExtMonFNALTrkFitQuality.hh"
#include "RecoDataProducts/inc/ExtMonFNALTrkClusterResiduals.hh"
#include "RecoDataProducts/inc/ExtMonFNALTrkFit.hh"
#include "RecoDataProducts/inc/ExtMonFNALTrkFitCollection.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "ExtinctionMonitorFNAL/Reconstruction/inc/TrackExtrapolator.hh"
#include "ExtinctionMonitorFNAL/Reconstruction/inc/PixelRecoUtils.hh"
#include "ExtinctionMonitorFNAL/Reconstruction/inc/ClusterOnTrackPrecisionTool.hh"
#include "ExtinctionMonitorFNAL/Reconstruction/inc/LinearRegression.hh"

#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "TH1D.h"
#include "TH2D.h"

//#define AGDEBUG(stuff) do { std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<", func "<<__func__<<": "<<stuff<<std::endl; } while(0)
#define AGDEBUG(stuff)

namespace mu2e {
  namespace ExtMonFNAL {

    //================================================================
    class HistTracklet : private boost::noncopyable {
    public:
      // Book histograms in the subdirectory, given by the relativePath; that path is
      // relative to the root TFileDirectory for the current module.
      // The default is to use the current module's TFileDirectory
      void book(const ExtMon& extmon, const std::string& relativePath="");

      // Book histograms in the specified TFileDirectory.
      void book(const ExtMon& extmon, art::TFileDirectory& tfdir);

      void fill(const Tracklet& tl);
      void fill(const Tracklets& tls);

    private:
      TH1D* hclock_;
      TH2D* hPosBegin_;
      TH2D* hPosEnd_;
      TH2D* hSlopes_;
    };

    // Book histograms in the subdirectory, given by the relativePath; that path is
    // relative to the root TFileDirectory for the current module.
    void HistTracklet::book(const ExtMonFNAL::ExtMon& extmon, const std::string& relativePath)
    {
      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory tfdir = relativePath.empty() ? *tfs : tfs->mkdir(relativePath.c_str());
      book (extmon, tfdir);
    }

    // Book the histograms.
    void HistTracklet::book(const ExtMonFNAL::ExtMon& extmon, art::TFileDirectory& tfdir) {

      hclock_ = tfdir.make<TH1D>("clock", "Tracklet clock", 100, -20.5, 79.5);

      hPosBegin_ = tfdir.make<TH2D>("posBegin", "Tracklet start position",
                                    400, -extmon.plane(4).halfSize()[0], +extmon.plane(4).halfSize()[0],
                                    400, -extmon.plane(4).halfSize()[1], +extmon.plane(4).halfSize()[1]
                                    );
      hPosBegin_->SetOption("colz");

      hPosEnd_ = tfdir.make<TH2D>("posEnd", "Tracklet end position",
                                  400, -extmon.plane(4).halfSize()[0], +extmon.plane(4).halfSize()[0],
                                  400, -extmon.plane(4).halfSize()[1], +extmon.plane(4).halfSize()[1]
                                  );
      hPosEnd_->SetOption("colz");

      hSlopes_ = tfdir.make<TH2D>("slopes", "Tracklet slope Y vs X", 300, -0.15, +0.15, 300, -0.15, +0.15);
      hSlopes_->SetOption("colz");

    } // end HistTracklet::book()

    void HistTracklet::fill(const Tracklet& tl) {
      hclock_->Fill(0.5*(tl.firstCluster->clock() + tl.lastCluster->clock()));
      hPosBegin_->Fill(tl.firstCluster->position().x(), tl.firstCluster->position().y());
      hPosEnd_->Fill(tl.lastCluster->position().x(), tl.lastCluster->position().y());

      const CLHEP::Hep3Vector dir(tl.lastCluster->position() - tl.firstCluster->position());
      hSlopes_->Fill(dir.x()/dir.z(), dir.y()/dir.z());

    } // end HistTracklet::fill()

    void HistTracklet::fill(const Tracklets& tls) {
      for(Tracklets::const_iterator i=tls.begin(); i!=tls.end(); ++i) {
        fill(*i);
      }
    }

    //================================================================
    class HistTrkMatch : private boost::noncopyable {
    public:
      // Book histograms in the subdirectory, given by the relativePath; that path is
      // relative to the root TFileDirectory for the current module.
      // The default is to use the current module's TFileDirectory
      void book(const ExtMon& extmon, const std::string& relativePath="");

      // Book histograms in the specified TFileDirectory.
      void book(const ExtMon& extmon, art::TFileDirectory& tfdir);

      void fill(const Tracklet& up, const Tracklet& dn);
      void fill(const Tracklets& ups, const Tracklets& dns);

    private:
      TH2D* hclockMatch_;
      TH1D* hslopexMatch_;
    };

    // Book histograms in the subdirectory, given by the relativePath; that path is
    // relative to the root TFileDirectory for the current module.
    void HistTrkMatch::book(const ExtMonFNAL::ExtMon& extmon, const std::string& relativePath)
    {
      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory tfdir = relativePath.empty() ? *tfs : tfs->mkdir(relativePath.c_str());
      book (extmon, tfdir);
    }

    // Book the histograms.
    void HistTrkMatch::book(const ExtMonFNAL::ExtMon& extmon, art::TFileDirectory& tfdir) {

      hclockMatch_ = tfdir.make<TH2D>("clockMatch", "Tracklet clock dn vs up",
                                      100, -20.5, 79.5, 100, -20.5, 79.5);
      hclockMatch_->SetOption("colz");

      hslopexMatch_ = tfdir.make<TH1D>("slopexMatch", "Tracklet slope X dn vs up", 300, -0.015, +0.015);

    } // end HistTrkMatch::book()

    void HistTrkMatch::fill(const Tracklet& up, const Tracklet& dn) {

      hclockMatch_->Fill(0.5*(up.firstCluster->clock() + up.lastCluster->clock()),
                         0.5*(dn.firstCluster->clock() + dn.lastCluster->clock()));

      const CLHEP::Hep3Vector diru(up.lastCluster->position() - up.firstCluster->position());
      const CLHEP::Hep3Vector dird(dn.lastCluster->position() - dn.firstCluster->position());
      hslopexMatch_->Fill(dird.x()/dird.z() - diru.x()/diru.z());
    } // end HistTrkMatch::fill()

    void HistTrkMatch::fill(const Tracklets& ups, const Tracklets& dns) {
      for(Tracklets::const_iterator u=ups.begin(); u!=ups.end(); ++u) {
        for(Tracklets::const_iterator d=dns.begin(); d!=dns.end(); ++d) {
          fill(*u, *d);
        }
      }
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
        , lr_()

        , clusterClockTolerance_(pset.get<unsigned>("clusterClockTolerance"))
        , slopexmax_(pset.get<double>("maxTrackSlope"))
        , slopeymax_(slopexmax_)
        , stackScatterAngleTolerance_(pset.get<double>("stackScatterAngleTolerance"))

        , trackletMatchSlopeXTolerance_(pset.get<double>("trackletMatchSlopeXTolerance"))

        , alignmentToleranceX_(pset.get<double>("alignmentToleranceX"))
        , alignmentToleranceY_(pset.get<double>("alignmentToleranceY"))
        , thetaScatterOnePlane_(pset.get<double>("thetaScatterOnePlane"))
        , cutMinTrackProb_(pset.get<double>("cutMinTrackProb"))

        , hTrackletMultiplicity_()
        , hTrackletMatchSlopeX_()
        , hClusterAddXY_()
        , hClusterAddXXMax_()
        , hClusterAddYYMax_()
        , hClockDiffTrackletSeedClusters_()
        , hClockDiffClusterTracklet_()

        , singleParticleMode_(pset.get<bool>("singleParticleMode", false))
      {
        produces<ExtMonFNALTrkFitCollection>();

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
      ClusterOnTrackPrecisionTool clTool_;
      LinearRegression lr_;

      //----------------------------------------------------------------
      // Track search parameters

      int clusterClockTolerance_;

      // We look for tracks that are parallel to the detector axis within
      double slopexmax_;
      double slopeymax_;
      double stackScatterAngleTolerance_;
      double trackletMatchSlopeXTolerance_;
      double alignmentToleranceX_;
      double alignmentToleranceY_;
      double thetaScatterOnePlane_;
      double cutMinTrackProb_;

      //----------------------------------------------------------------
      HistTracklet htup_;
      HistTracklet htdn_;
      HistTrkMatch hudm_;

      TH2D *hTrackletMultiplicity_;
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
      void formTrackletsWithPlanes(const ExtMonFNALPlaneStack& stack, Tracklets& res,
			      const unsigned plane1, const unsigned plane2,
                              const art::Handle<ExtMonFNALRecoClusterCollection>& clusters);

      Tracklets formTracklets(const ExtMonFNALPlaneStack& stack,

                              const art::Handle<ExtMonFNALRecoClusterCollection>& clusters);

      void findTracks(art::Event& event,
                      ExtMonFNALTrkFitCollection *tracks,
                      const art::Handle<ExtMonFNALRecoClusterCollection>& clusters);

      void addPlaneClusters(Tracklets* seeds,
                            const art::Handle<ExtMonFNALRecoClusterCollection>& coll,
                            const ExtMonFNALPlaneStack& stack,
                            unsigned stackPlane,
                            unsigned anchorPlane1,
                            unsigned anchorPlane2
                            );

      Tracklet createTracklet(const Tracklet& orig, const art::Ptr<ExtMonFNALRecoCluster>& cl);

      void modifyInPlace(Tracklet *inout, const art::Ptr<ExtMonFNALRecoCluster>& cl);

      double slopex(const Tracklet& tl);

      bool inTime(const Tracklet& tl, const ExtMonFNALRecoCluster& cl);
      bool inTime(const Tracklet& tl1, const Tracklet& tl2);

      void addToClusters(std::vector<art::Ptr<ExtMonFNALRecoCluster> > *clusters, const Tracklet& tl);

      ExtMonFNALTrkFitQuality evaluateFit(std::vector<ExtMonFNALTrkClusterResiduals> *residuals,
                                          const ExtMonFNALTrkParam& pars,
                                          const std::vector<art::Ptr<ExtMonFNALRecoCluster> >& clusters);
    };

    //================================================================
    void EMFPatRecFromTracklets::beginRun(art::Run& run) {
      // This is a workaround for geometry not being available at beginJob()
      if(!extmon_) {

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
          clTool_ = ClusterOnTrackPrecisionTool(*extmon_, thetaScatterOnePlane_);
          lr_ = LinearRegression(extmon_, clTool_);
        }
        else {
          if(verbosityLevel_ > 0) {
            std::cout<<"EMFPatRecFromTracklets: using GeometryService"<<std::endl;
          }
          GeomHandle<ExtMonFNAL::ExtMon> emf;
          extmon_ = &*emf;
          extrapolator_ = TrackExtrapolator(extmon_);
          clTool_ = ClusterOnTrackPrecisionTool(*extmon_, thetaScatterOnePlane_);
          lr_ = LinearRegression(extmon_, clTool_);
        }

        std::cout<<"EMFPatRecFromTracklets: spectrometer nominalMomentum = "
                 <<extmon_->spectrometerMagnet().nominalMomentum()
                 <<std::endl;

        //----------------
        htup_.book(*extmon_, "TrackletsUp");
        htdn_.book(*extmon_, "TrackletsDn");
        hudm_.book(*extmon_, "TrackletsMatch");

        art::ServiceHandle<art::TFileService> tfs;
        hTrackletMultiplicity_ = tfs->make<TH2D>("trackletMultiplicity", "Tracklet multiplicity down- vs upstream",
                                                 200, -0.5, 199.5, 200, -0.5, 199.5);
        hTrackletMultiplicity_->SetOption("colz");

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
    }

    //================================================================
    void EMFPatRecFromTracklets::produce(art::Event& event) {
      art::Handle<ExtMonFNALRecoClusterCollection> clusters;
      event.getByLabel(inputModuleLabel_, inputInstanceName_, clusters);

      std::unique_ptr<ExtMonFNALTrkFitCollection> tracks(new ExtMonFNALTrkFitCollection);

      if(!singleParticleMode_ || acceptSingleParticleEvent(event)) {
        findTracks(event, &*tracks, clusters);
      }

      event.put(std::move(tracks));
    }

    //================================================================
    void EMFPatRecFromTracklets::findTracks(art::Event& event,
                                            ExtMonFNALTrkFitCollection *tracks,
                                            const art::Handle<ExtMonFNALRecoClusterCollection>& coll)
    {
      Tracklets tup = formTracklets(extmon_->up(), coll);
      Tracklets tdn = formTracklets(extmon_->dn(), coll);
      hTrackletMultiplicity_->Fill(tdn.size(), tup.size());
      htup_.fill(tup);
      htdn_.fill(tdn);
      hudm_.fill(tup, tdn);

      // merge compatible tracklet pairs into tracks
      for(Tracklets::const_iterator iup = tup.begin(); iup != tup.end(); ++iup) {
        for(Tracklets::const_iterator idn = tdn.begin(); idn != tdn.end(); ++idn) {

          const double dslopex = slopex(*iup) - slopex(*idn);
          hTrackletMatchSlopeX_->Fill(dslopex);

          if(inTime(*iup, *idn) && (std::abs(dslopex) < trackletMatchSlopeXTolerance_)) {

            // Compute track parameter estimates
            ExtMonFNALTrkParam trkpar = lr_.estimatePars(*iup, *idn);

            // Compute fit quality and residuals
            std::vector<art::Ptr<ExtMonFNALRecoCluster> > clusters;
            addToClusters(&clusters, *idn);
            addToClusters(&clusters, *iup);

            std::vector<ExtMonFNALTrkClusterResiduals> residuals;
            ExtMonFNALTrkFitQuality quality =
              evaluateFit(&residuals, trkpar, clusters);

            Genfun::CumulativeChiSquare pf(quality.ndf());
            const double prob = 1. - pf(quality.chi2());
            if(cutMinTrackProb_ <= prob) { // Accept the track
              tracks->push_back(ExtMonFNALTrkFit(trkpar, quality, clusters, residuals));
            } // if(fit quality)

          } // if(inTime and slope match)
        }
      }

    }

    //================================================================
    void EMFPatRecFromTracklets::formTrackletsWithPlanes(const ExtMonFNALPlaneStack& stack, Tracklets& res,
						    const unsigned plane1, const unsigned plane2,
                                                    const art::Handle<ExtMonFNALRecoClusterCollection>& coll)
    {
      const double dz = stack.plane_zoffset().back() - stack.plane_zoffset().front();
      const double dxmax = std::abs(dz * slopexmax_);
      const double dymax = std::abs(dz * slopeymax_);

      typedef ExtMonFNALRecoClusterCollection::PlaneClusters PC;
      const PC pc1 = coll->clusters(plane1);
      const PC pc2 = coll->clusters(plane2);

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


    }

    //================================================================
    Tracklets EMFPatRecFromTracklets::formTracklets(const ExtMonFNALPlaneStack& stack,
                                                    const art::Handle<ExtMonFNALRecoClusterCollection>& coll)
    {
      Tracklets res;

      // find tracklets starting with pair combinations of clusters on the even stack planes
      const unsigned plane1 = stack.planeNumberOffset() + 0;
      const unsigned plane2 = stack.planeNumberOffset() + stack.nplanes() - 2;
      AGDEBUG("plane1 = "<<plane1<<", plane2 = "<<plane2);

      formTrackletsWithPlanes(stack,res,plane1,plane2,coll);
 
     // add tracklets starting with pair combinations of clusters on the odd stack planes
      const unsigned plane3 = stack.planeNumberOffset() + 1;
      const unsigned plane4 = stack.planeNumberOffset() + stack.nplanes() - 1;
      AGDEBUG("plane3 = "<<plane3<<", plane4 = "<<plane4);

      formTrackletsWithPlanes(stack,res,plane3,plane4,coll);

      return res;
    }

    //================================================================
    void EMFPatRecFromTracklets::addPlaneClusters(Tracklets* seeds,
                                                  const art::Handle<ExtMonFNALRecoClusterCollection>& coll,
                                                  const ExtMonFNALPlaneStack& stack,
                                                  unsigned stackPlane,
                                                  unsigned anchorPlane1,
                                                  unsigned anchorPlane2
                                                  )
    {
      const double planeZ = stack.plane_zoffset()[stackPlane];
      AGDEBUG("Here: stackPlane="<<stackPlane<<", planeZ = "<<planeZ<<", seeds->size() = "<<seeds->size()
              <<", anchorPlane1 = "<<anchorPlane1<<", anchorPlane2 = "<<anchorPlane2);

      const double dxtol1 =  alignmentToleranceX_ + stackScatterAngleTolerance_ * std::abs(planeZ - stack.plane_zoffset()[anchorPlane1]);
      const double dxtol2 =  alignmentToleranceX_ + stackScatterAngleTolerance_ * std::abs(planeZ - stack.plane_zoffset()[anchorPlane2]);

      const double dytol1 =  alignmentToleranceY_ + stackScatterAngleTolerance_ * std::abs(planeZ - stack.plane_zoffset()[anchorPlane1]);
      const double dytol2 =  alignmentToleranceY_ + stackScatterAngleTolerance_ * std::abs(planeZ - stack.plane_zoffset()[anchorPlane2]);
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
          // Must be at least one match.  Kill the seed
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
    bool EMFPatRecFromTracklets::acceptSingleParticleEvent(const art::Event& event) {
      art::Handle<ExtMonFNALRecoClusterCollection> coll;
      event.getByLabel(inputModuleLabel_, inputInstanceName_, coll);
      return perfectSingleParticleEvent(*coll, extmon_->nplanes());
    }

    //================================================================
    void EMFPatRecFromTracklets::
    addToClusters(std::vector<art::Ptr<ExtMonFNALRecoCluster> > *clusters, const Tracklet& tl) {
      clusters->push_back(tl.firstCluster);
      for(unsigned i=0; i<tl.middleClusters.size(); ++i) {
        clusters->push_back(tl.middleClusters[i]);
      }
      clusters->push_back(tl.lastCluster);
    }

    //================================================================
    ExtMonFNALTrkFitQuality EMFPatRecFromTracklets::
    evaluateFit(std::vector<ExtMonFNALTrkClusterResiduals> *residuals,
                const ExtMonFNALTrkParam& pars,
                const std::vector<art::Ptr<ExtMonFNALRecoCluster> >& clusters) {

      residuals->resize(clusters.size());
      double chi2(0);

      ExtMonFNALTrkParam tmp(pars); // updated on each loop iteration

      for(int i=clusters.size()-1; i>=0; --i) {

        const ExtMonFNALRecoCluster& cl(*clusters[i]);
        if(extrapolator_.extrapolateToPlane(cl.plane(), &tmp)) {
          const double dx = cl.position().x() - tmp.posx();
          const double dy = cl.position().y() - tmp.posy();

          const ClusterOnTrackPrecision cp = clTool_.clusterPrecision(cl);
          const double ux = dx/sqrt(cp.sigma2x);
          const double uy = dy/sqrt(cp.sigma2y);

          (*residuals)[i] = ExtMonFNALTrkClusterResiduals(dx, ux, dy, uy);

          chi2 += ux*ux + uy*uy;
        }
        else {
          throw cet::exception("RECO")<<"EMFPatRecFromTracklets::evaluateFit(): error in track extrapolation\n";
        }

      } // for(clusters)

      // Each cluster is two measurements.
      const int ndf = 2*clusters.size() - 5;

      return ExtMonFNALTrkFitQuality(ndf, chi2);

    } //evaluateFit()

    //================================================================

  } // end namespace ExtMonFNAL
} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::ExtMonFNAL::EMFPatRecFromTracklets)
