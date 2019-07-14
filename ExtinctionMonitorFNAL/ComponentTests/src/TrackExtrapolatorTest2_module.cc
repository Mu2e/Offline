// Compare extrapolator result with truth info for G4 simulated single particles.
//
// Andrei Gaponenko, 2012

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <sstream>

#include "boost/noncopyable.hpp"

#include "cetlib_except/exception.h"

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Utilities/InputTag.h"

#include "RecoDataProducts/inc/ExtMonFNALRecoCluster.hh"
#include "RecoDataProducts/inc/ExtMonFNALRecoClusterCollection.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "RecoDataProducts/inc/ExtMonFNALTrkParam.hh"

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALMagnet.hh"
#include "GeometryService/inc/GeomHandle.hh"

#include "ExtinctionMonitorFNAL/Reconstruction/inc/TrackExtrapolator.hh"

#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "TH1D.h"
#include "TH2D.h"

namespace mu2e {
  namespace ExtMonFNAL {

    //================================================================
    namespace {
      class TestHistos // : private boost::noncopyable
      {

        TH2D *hxcorr_;
        TH2D *hycorr_;
        TH2D *hresiduals_;

      public:

        TestHistos() : hxcorr_(), hycorr_(), hresiduals_() {}

        // Book histograms in the subdirectory, given by the relativePath; that path is
        // relative to the root TFileDirectory for the current module.
        // The default is to use the current module's TFileDirectory
        void book(const ExtMonFNAL::ExtMon& extmon, const std::string& relativePath="");

        // Book histograms in the specified TFileDirectory.
        void book(const ExtMonFNAL::ExtMon& extmon, art::TFileDirectory& tfdir);

        void fill(const ExtMonFNALTrkParam& extrapolated, const ExtMonFNALRecoCluster& cluster);
      };

      void TestHistos::book(const ExtMon& extmon, const std::string& relativePath)
      {
        art::ServiceHandle<art::TFileService> tfs;
        art::TFileDirectory tfdir = relativePath.empty() ? *tfs : tfs->mkdir(relativePath.c_str());
        book (extmon, tfdir);
      }

      // Book the histograms.
      void TestHistos::book(const ExtMon& extmon, art::TFileDirectory& tfdir) {
        // Bin cluster position histograms according to pixel size
        const unsigned nx = extmon.sensor().nxChips() * extmon.chip().nColumns();
        const double   lx = nx * extmon.chip().xPitch();
        const unsigned ny = extmon.sensor().nyChips() * extmon.chip().nRows();
        const double   ly = ny * extmon.chip().yPitch();

        hxcorr_ = tfdir.make<TH2D>("xcorr", "x extrapolated vs hit", nx, -lx/2, +lx/2, nx, -lx/2, +lx/2);
        hxcorr_->SetOption("colz");

        hycorr_ = tfdir.make<TH2D>("ycorr", "y extrapolated vs hit", ny, -ly/2, +ly/2, ny, -ly/2, +ly/2);
        hycorr_->SetOption("colz");


        // Smaller range for residuals
        const double approximateHalfRange = 1.; // mm

        const double halfXbinsEstimate = approximateHalfRange / extmon.chip().xPitch();
        const int nXbins =  1 + 2*halfXbinsEstimate;
        const double halfXRange = 0.5 * nXbins * extmon.chip().xPitch();

        const double halfYbinsEstimate = approximateHalfRange / extmon.chip().yPitch();
        const int nYbins =  1 + 2*halfYbinsEstimate;
        const double halfYRange = 0.5 * nYbins * extmon.chip().yPitch();

        hresiduals_ = tfdir.make<TH2D>("residuals", "(extrapolator-cluster) position y vs x",
                                       nXbins, -halfXRange, + halfXRange,
                                       nYbins, -halfYRange, + halfYRange);

        hresiduals_->SetOption("colz");
        hresiduals_->GetXaxis()->SetTitle("dx [mm]");
        hresiduals_->GetYaxis()->SetTitle("dy [mm]");

      } // end TestHistos::book()

      void TestHistos::fill(const ExtMonFNALTrkParam& extrapolated,
                            const ExtMonFNALRecoCluster& cluster)
      {
        hxcorr_->Fill(cluster.position().x(), extrapolated.posx());

        hycorr_->Fill(cluster.position().y(), extrapolated.posy());

        hresiduals_->Fill(extrapolated.posx() - cluster.position().x(),
                          extrapolated.posy() - cluster.position().y());
      }

    } // namespace {}


    //================================================================
    class TrackExtrapolatorTest2 : public art::EDAnalyzer {
      std::string particleModuleLabel_;
      std::string particleInstanceName_;
      std::string recoClusterModuleLabel_;
      std::string recoClusterInstanceName_;
      std::string geomModuleLabel_; // emtpy to take info from Run
      std::string geomInstanceName_;

      const ExtMonFNAL::ExtMon *extmon_;

      // by plane
      std::map<unsigned, TestHistos> hh_;

    public:
      explicit TrackExtrapolatorTest2(const fhicl::ParameterSet& pset);
      virtual void beginRun(const art::Run& run);
      virtual void analyze(const art::Event& event);
    };

    //================================================================
    TrackExtrapolatorTest2::TrackExtrapolatorTest2(const fhicl::ParameterSet& pset)
      : art::EDAnalyzer(pset)
      , particleModuleLabel_(pset.get<std::string>("particleModuleLabel"))
      , particleInstanceName_(pset.get<std::string>("particleInstanceName", ""))
      , recoClusterModuleLabel_(pset.get<std::string>("recoClusterModuleLabel"))
      , recoClusterInstanceName_(pset.get<std::string>("recoClusterInstanceName", ""))
      , geomModuleLabel_(pset.get<std::string>("geomModuleLabel"))
      , geomInstanceName_(pset.get<std::string>("geomInstanceName", ""))

      , extmon_()

      , hh_()
    {}

    //================================================================
    void TrackExtrapolatorTest2::beginRun(const art::Run& run) {
      if(!geomModuleLabel_.empty()) {
        art::Handle<ExtMonFNAL::ExtMon> emf;
        run.getByLabel(geomModuleLabel_, geomInstanceName_, emf);
        extmon_ = &*emf;
      }
      else {
        GeomHandle<ExtMonFNAL::ExtMon> emf;
        extmon_ = &*emf;
      }

      //----------------------------------------------------------------
      const unsigned nplanes = extmon_->dn().nplanes() + extmon_->up().nplanes();
      for(unsigned plane=0; plane < nplanes; ++plane) {
        std::ostringstream osname;
        osname<<"plane"<<plane;
        hh_[plane].book(*extmon_, osname.str());
      }
    }

    //================================================================
    void TrackExtrapolatorTest2::analyze(const art::Event& event) {

      art::Handle<SimParticleCollection> ih;
      event.getByLabel(particleModuleLabel_, particleInstanceName_, ih);
      const SimParticleCollection& particles(*ih);

      if(particles.size() == 1) { // look only at unambiguous one-to-one cases

        art::Handle<ExtMonFNALRecoClusterCollection> ch;
        event.getByLabel(recoClusterModuleLabel_, recoClusterInstanceName_, ch);
        const ExtMonFNALRecoClusterCollection& cc(*ch);

        bool allSingleHits(true);
        const unsigned nplanes = extmon_->dn().nplanes() + extmon_->up().nplanes();
        for(unsigned i=0; i<nplanes; ++i) {
          const ExtMonFNALRecoClusterCollection::PlaneClusters& pc = cc.clusters(i);
          allSingleHits &= (pc.size() == 1);
        }

        if(allSingleHits) { // look only at unambiguous one-to-one cases

          const ExtMonFNALMagnet& mag = extmon_->spectrometerMagnet();
          TrackExtrapolator ex(&*extmon_);

          const CLHEP::Hep3Vector& startPos = extmon_->mu2eToExtMon_position(particles.front().second.startPosition());
          const CLHEP::Hep3Vector& startMom = extmon_->mu2eToExtMon_momentum(particles.front().second.startMomentum());

          ExtMonFNALTrkParam start;
          start.setposx(startPos.x());
          start.setposy(startPos.y());
          start.setz0(startPos.z());
          start.setslopex(startMom.x()/startMom.z());
          start.setslopey(startMom.y()/startMom.z());
          start.setrinv(1./mag.trackBendRadius(startMom.mag()));

          std::cout<<"AG: start = "<<start<<" sp ="<<startPos<<", sm = "<<startMom<<std::endl;

          for(unsigned i = 0; i < nplanes; ++i) {
            ExtMonFNALTrkParam current(start);
            if(ex.extrapolateToPlane(i, &current)) {
              std::cout<<"AG: current = "<<current<<" for plane "<<i<<std::endl;
              hh_[i].fill(current, cc.clusters(i)[0]);
            }
            else {
              throw cet::exception("EXTMON")<<": Could not extrapolate "<<start<<" to plane 0\n";
            }
          }
        }
      }
    } // analyze()

    //================================================================
  } // namespace ExtMonFNAL
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::ExtMonFNAL::TrackExtrapolatorTest2);
