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

#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"

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

#include "art_root_io/TFileService.h"
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
      hclock_->Fill(0.5*(tl.firstSeedCluster->clock() + tl.secondSeedCluster->clock()));
      hPosBegin_->Fill(tl.firstSeedCluster->position().x(), tl.firstSeedCluster->position().y());
      hPosEnd_->Fill(tl.secondSeedCluster->position().x(), tl.secondSeedCluster->position().y());

      const CLHEP::Hep3Vector dir(tl.secondSeedCluster->position() - tl.firstSeedCluster->position());
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

      hclockMatch_->Fill(0.5*(up.firstSeedCluster->clock() + up.secondSeedCluster->clock()),
                         0.5*(dn.firstSeedCluster->clock() + dn.secondSeedCluster->clock()));

      const CLHEP::Hep3Vector diru(up.secondSeedCluster->position() - up.firstSeedCluster->position());
      const CLHEP::Hep3Vector dird(dn.secondSeedCluster->position() - dn.firstSeedCluster->position());
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
        : EDProducer{pset}
        , verbosityLevel_(pset.get<int>("verbosityLevel", 0))
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

        , maxMissedHits_(pset.get<unsigned>("maxMissedHits"))

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

      unsigned maxMissedHits_;

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
      void findTracks(art::Event& event,
                      ExtMonFNALTrkFitCollection *tracks,
                      const art::Handle<ExtMonFNALRecoClusterCollection>& clusters);

      Tracklets formTracklets(const ExtMonFNALPlaneStack& stack,
                              const art::Handle<ExtMonFNALRecoClusterCollection>& clusters);

      // finds all clusters in stackPlane that are compatible with a straigth line track going through c1 and c2
      std::vector<unsigned> findCompatibleClusterIndices(
							 const ExtMonFNALRecoCluster& c1, 
							 const ExtMonFNALRecoCluster& c2,
							 const art::Handle<ExtMonFNALRecoClusterCollection>& coll,
							 const ExtMonFNALPlaneStack& stack,
							 unsigned stackPlane
							 );

      Tracklet createTracklet(const Tracklet& orig, const art::Ptr<ExtMonFNALRecoCluster>& cl);

      void modifyInPlace(Tracklet *inout, const art::Ptr<ExtMonFNALRecoCluster>& cl);

      double slopex(const Tracklet& tl);

      bool inTime(const ExtMonFNALRecoCluster& c1, const ExtMonFNALRecoCluster& c2);
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
    Tracklets EMFPatRecFromTracklets::formTracklets(const ExtMonFNALPlaneStack& stack,
                                                    const art::Handle<ExtMonFNALRecoClusterCollection>& coll)
    {
      Tracklets res;
      const double dz = stack.plane_zoffset().back() - stack.plane_zoffset().front();
      const double dxmax = std::abs(dz * slopexmax_);
      const double dymax = std::abs(dz * slopeymax_);
      typedef ExtMonFNALRecoClusterCollection::PlaneClusters PC;
      // each vector is nplanes long and contains cluster indexes or -1u
      // This is the list of cluster combinations found using previous seed planes
      std::list< std::vector<unsigned> > foundCombinations;
      for(unsigned seedPlane=0; seedPlane <= maxMissedHits_; ++seedPlane) {
        unsigned backPlane = stack.nplanes() - maxMissedHits_  + seedPlane - 1;
	const unsigned plane1 = stack.planeNumberOffset() + seedPlane;
        const unsigned plane2 = stack.planeNumberOffset() + backPlane;


	// This vector is used to simplify loooping over the non-seed planes
	std::vector< unsigned > additionalPlanes ( stack.nplanes() -2, 0 );
	for(unsigned count = 0, stackPlane=0; stackPlane < stack.nplanes(); ++stackPlane) {
          if( stackPlane == seedPlane ) continue;
          if( stackPlane == backPlane ) continue;
	  additionalPlanes[count] = stackPlane; ++count;
	}

        const PC pc1 = coll->clusters(plane1);
        const PC pc2 = coll->clusters(plane2);

        // This will hold the list of cluster combinations found using this seed plane
	std::list< std::vector<unsigned> > newCombinations;

        // find seed pairs
        for(unsigned i1 = 0; i1 < pc1.size(); ++i1) {
          art::Ptr<ExtMonFNALRecoCluster> c1(coll, coll->globalIndex(plane1, i1));
          for(unsigned i2 = 0; i2 < pc2.size(); ++i2) {
            art::Ptr<ExtMonFNALRecoCluster> c2(coll, coll->globalIndex(plane2, i2));

            hClockDiffTrackletSeedClusters_->Fill(c2->clock() - c1->clock());

            if( (std::abs(c1->clock() - c2->clock()) <= clusterClockTolerance_) &&
		(std::abs(c1->position().x() - c2->position().x()) < dxmax) &&
		(std::abs(c1->position().y() - c2->position().y()) < dymax) )
              {

              // This will hold the list of cluster combinations found using this seed pair
	      std::list< std::vector<unsigned> > newSeedCombinations;
	
	      for(unsigned p=0; p < stack.nplanes() - 2; ++p) {
	        unsigned stackPlane = additionalPlanes[p];
	        std::vector<unsigned> clusters = findCompatibleClusterIndices(*c1, *c2, coll, stack, stackPlane);

                // This will hold the list of extra cluster combinations found using this plane
	        std::list< std::vector<unsigned> > newPlaneCombinations;

	        for( unsigned i = 0; i < clusters.size(); i++ ) {

	          // check if this combination was found using an earlier seed plane
	          bool notfound = true;
	          for(auto comb = foundCombinations.begin(); comb != foundCombinations.end(); ++comb ) {
	            if( (*comb)[seedPlane] == i1 &&
	        	(*comb)[backPlane] == i2 &&
	        	(*comb)[stackPlane] == i )
	              {
	        	notfound = false; break;
	              }
	          }

	          if( notfound ) {
	            if( newSeedCombinations.size() == 0 ) {
		    // This is the first new combination using this seed pair
	        	std::vector<unsigned> planes( stack.nplanes(), -1u );
	        	planes[seedPlane] = i1;
	        	planes[backPlane] = i2;
	        	newSeedCombinations.push_back(planes);
		    }
		    for(auto comb = newSeedCombinations.begin(); comb != newSeedCombinations.end(); ++comb ) {
		    // Form all possible coombinations of this plane's hits with those of earlier planes
		    //  the 1st one can be included into the existing combinations will go into copies of them 
	              if( i == 0 ) {
	        	(*comb)[stackPlane] = i;
	              } else {
	        	std::vector<unsigned> planes( *comb );
	        	planes[stackPlane] = i;
	        	newPlaneCombinations.push_back(planes);
	              }
	            }

	          } // if combination found
	        } // additional plane clusters

                newSeedCombinations.splice( newSeedCombinations.begin(), newPlaneCombinations );

	      } // additional planes

	      newCombinations.splice( newCombinations.begin(), newSeedCombinations );

	    } // if seed ok
          } // back plane clusters
        } // seed plane clusters

	// create new tracklets
	for(auto comb = newCombinations.begin(); comb != newCombinations.end(); ++comb ) {
	  // this piece of code can be replaced with an appropriate Tracklet constructor
	  std::vector<art::Ptr<ExtMonFNALRecoCluster> > additionalClusters;
          art::Ptr<ExtMonFNALRecoCluster> c1(coll, coll->globalIndex(plane1, (*comb)[seedPlane]));
          art::Ptr<ExtMonFNALRecoCluster> c2(coll, coll->globalIndex(plane2, (*comb)[backPlane]));
	  for(unsigned p=0; p < stack.nplanes() - 2; ++p) {
	    unsigned stackPlane = additionalPlanes[p];
	    if( (*comb)[stackPlane] == -1u ) continue;
	    unsigned globalPlane = stack.planeNumberOffset() + stackPlane;
	    additionalClusters.push_back(art::Ptr<ExtMonFNALRecoCluster>(coll, coll->globalIndex(globalPlane, (*comb)[stackPlane])));
	  }
	  if( additionalClusters.size() >= stack.nplanes() - maxMissedHits_  - 2 ) {
	    Tracklet track(c1,c2);
	    for( unsigned ic = 0; ic < additionalClusters.size(); ic++ ) {
	      modifyInPlace(&track, additionalClusters[ic]);
	    }
	    res.push_back( track );
	  }
	} // new combinations

	foundCombinations.splice( foundCombinations.begin(), newCombinations );

      } // seed plane

      return res;
    }

    //================================================================
    std::vector<unsigned> EMFPatRecFromTracklets::findCompatibleClusterIndices(const ExtMonFNALRecoCluster& c1,
									       const ExtMonFNALRecoCluster& c2,
									       const art::Handle<ExtMonFNALRecoClusterCollection>& coll,
									       const ExtMonFNALPlaneStack& stack,
									       unsigned stackPlane
									       )
    {
      const unsigned anchorPlane1 = c1.plane();
      const unsigned anchorPlane2 = c2.plane();

      const double planeZ = stack.plane_zoffset()[stackPlane];
      AGDEBUG("Here: stackPlane="<<stackPlane<<", planeZ = "<<planeZ<<", seeds->size() = "<<seeds->size()
              <<", anchorPlane1 = "<<anchorPlane1<<", anchorPlane2 = "<<anchorPlane2);

      const double dxtol1 =  alignmentToleranceX_ + stackScatterAngleTolerance_ * std::abs(planeZ - stack.plane_zoffset()[anchorPlane1]);
      const double dxtol2 =  alignmentToleranceX_ + stackScatterAngleTolerance_ * std::abs(planeZ - stack.plane_zoffset()[anchorPlane2]);

      const double dytol1 =  alignmentToleranceY_ + stackScatterAngleTolerance_ * std::abs(planeZ - stack.plane_zoffset()[anchorPlane1]);
      const double dytol2 =  alignmentToleranceY_ + stackScatterAngleTolerance_ * std::abs(planeZ - stack.plane_zoffset()[anchorPlane2]);
      AGDEBUG("dxtol1 ="<<dxtol1<<", dxtol2 ="<<dxtol2<<", dytol1 = "<<dytol1<<", dytol2 = "<<dytol2);

      std::vector<unsigned> res;

      const CLHEP::Hep3Vector& pos1 = c1.position();
      const CLHEP::Hep3Vector& pos2 = c2.position();
      AGDEBUG("pos1 = "<<pos1<<", pos2 = "<<pos2);

      const double z1(pos1.z()), z2(pos2.z());
      const double dz(z2-z1);

      const CLHEP::Hep3Vector interpolated = pos2 * (planeZ - z1)/dz + pos1 * (z2 - planeZ)/dz;

      const double dxmax1 = dxtol1 + 0.5*c1.xWidth()*extmon_->chip().xPitch();
      const double dxmax2 = dxtol2 + 0.5*c2.xWidth()*extmon_->chip().xPitch();
      // satisfy both constraints
      const double dxmax(std::min(dxmax1, dxmax2));

      const double dymax1 = dytol1 + 0.5*c1.yWidth()*extmon_->chip().yPitch();
      const double dymax2 = dytol2 + 0.5*c2.yWidth()*extmon_->chip().yPitch();
      // satisfy both constraints
      const double dymax(std::min(dymax1, dymax2));
      AGDEBUG("dxmax="<<dxmax<<", dymax="<<dymax<<", interpolated="<<interpolated);

      // find all clusters compatible with seed i
      const unsigned int globalPlane = stackPlane + stack.planeNumberOffset();
      const ExtMonFNALRecoClusterCollection::PlaneClusters& clusters = coll->clusters(globalPlane);
      for(unsigned ic=0; ic<clusters.size(); ++ic) {
	const ExtMonFNALRecoCluster& cl = clusters[ic];
	if(inTime(c1, cl) || inTime(c2, cl)) {
	  const double dx = interpolated.x() - cl.position().x();
	  const double dy = interpolated.y() - cl.position().y();

	  hClusterAddXY_[globalPlane]->Fill(dx/dxmax, dy/dymax);
	  hClusterAddXXMax_[globalPlane]->Fill(dxmax, dx);
	  hClusterAddYYMax_[globalPlane]->Fill(dymax, dy);

	  if((std::abs(dx) < dxmax) && (std::abs(dy) < dymax) ) { res.push_back(ic); }
	}
      }
      return res;
    }

    //================================================================
    Tracklet EMFPatRecFromTracklets::createTracklet(const Tracklet& orig, const art::Ptr<ExtMonFNALRecoCluster>& cl) {
      Tracklet res(orig);
      modifyInPlace(&res, cl);
      return res;
    }

    //================================================================
    void EMFPatRecFromTracklets::modifyInPlace(Tracklet *inout, const art::Ptr<ExtMonFNALRecoCluster>& cl) {
      inout->addedClusters.push_back(cl);
    }

    //================================================================
    double EMFPatRecFromTracklets::slopex(const Tracklet& tl) {
      const CLHEP::Hep3Vector dir(tl.secondSeedCluster->position() - tl.firstSeedCluster->position());
      return dir.x()/dir.z();
    }

    //================================================================
    bool EMFPatRecFromTracklets::inTime(const ExtMonFNALRecoCluster& c1, const ExtMonFNALRecoCluster& c2) {
      const int dt = c1.clock() - c2.clock();
      hClockDiffClusterTracklet_->Fill(dt);
      return (std::abs(dt) <= clusterClockTolerance_);
    }

    //================================================================
    bool EMFPatRecFromTracklets::inTime(const Tracklet& tl, const ExtMonFNALRecoCluster& cl) {
      const int dtFirst = cl.clock() - tl.firstSeedCluster->clock();
      const int dtLast  = cl.clock() - tl.secondSeedCluster->clock();

      hClockDiffClusterTracklet_->Fill(dtFirst);
      hClockDiffClusterTracklet_->Fill(dtLast);

      return
        (std::abs(dtFirst) <= clusterClockTolerance_) &&
        (std::abs(dtLast) <= clusterClockTolerance_);
    }

    //================================================================
    bool EMFPatRecFromTracklets::inTime(const Tracklet& tl1, const Tracklet& tl2) {
      return inTime(tl1, *tl2.firstSeedCluster) && inTime(tl1, *tl2.secondSeedCluster);
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
      clusters->push_back(tl.firstSeedCluster);
      for(unsigned i=0; i<tl.addedClusters.size(); ++i) {
        clusters->push_back(tl.addedClusters[i]);
      }
      clusters->push_back(tl.secondSeedCluster);
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
