// ExtMonFNAL PatRec efficiency/fake rate analysis
//
// Andrei Gaponenko, 2012

#include <iostream>
#include <string>
#include <vector>
#include <set>

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Utilities/InputTag.h"

#include "RecoDataProducts/inc/ExtMonFNALTrkFit.hh"
#include "RecoDataProducts/inc/ExtMonFNALTrkFitCollection.hh"
#include "RecoDataProducts/inc/ExtMonFNALRecoClusterCollection.hh"

#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/ExtMonFNALRecoClusterTruthAssn.hh"
#include "MCDataProducts/inc/ExtMonFNALPatRecTruthAssns.hh"

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "ExtinctionMonitorFNAL/Reconstruction/inc/PixelRecoUtils.hh"
#include "ExtinctionMonitorFNAL/Analyses/inc/EMFPatRecEffHistograms.hh"
#include "ExtinctionMonitorFNAL/Analyses/inc/EMFPatRecFakeHistograms.hh"

#include "ExtinctionMonitorFNAL/Reconstruction/inc/TrackExtrapolator.hh"

#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"

#define AGDEBUG(stuff) do { std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<", func "<<__func__<<": "<<stuff<<std::endl; } while(0)
//#define AGDEBUG(stuff)

namespace mu2e {
  namespace ExtMonFNAL {

    //================================================================
    class EMFDetMCHistPatRec : public art::EDAnalyzer {
      std::string patRecModuleLabel_;
      std::string patRecInstanceName_;
      std::string trkTruthModuleLabel_;
      std::string trkTruthInstanceName_;
      std::string clusterTruthModuleLabel_;
      std::string clusterTruthInstanceName_;
      std::string particleModuleLabel_;
      std::string particleInstanceName_;

      std::string geomModuleLabel_; // emtpy to take info from Run
      std::string geomInstanceName_;

      unsigned cutParticleMinClusters_; //

      // signal particle coordinates at first and last plane must be within the limits
      double cutHitXmax_;
      double cutHitYmax_;

      const ExtMonFNAL::ExtMon *extmon_;
      TrackExtrapolator extrapolator_;

      //----------------
      TH2D *hMultiplicitySignal_;
      TH2D *hCommonClusters_;

      TProfile *hRecoTruthRatio_;
      TProfile *hTruthRecoRatio_;
      TH2D * hRecoVsTruth_;

      EMFPatRecEffHistograms effPhysics_;
      EMFPatRecEffHistograms effSoftware_;
      EMFPatRecFakeHistograms fakes_;

      bool signalParticlePhysics(const SimParticle& particle);
      bool inAcceptance(const ExtMonFNALTrkParam& par);

      bool signalParticleSofware(const art::FindMany<ExtMonFNALRecoCluster,ExtMonFNALRecoClusterTruthBits>& clusterFinder,
                                 unsigned iParticle);

      //----------------------------------------------------------------
      // cuts tuning: single particle mode
      bool singleParticleMode_;
      std::string clusterModuleLabel_; // used only in single particle mode
      std::string clusterInstanceName_; // used only in single particle mode
      bool printEventsWithFakes_;
      bool printInefficiencies_;
      bool acceptSingleParticleEvent(const art::Event& event);

      void fillRTR(const art::Event& event);

    public:
      explicit EMFDetMCHistPatRec(const fhicl::ParameterSet& pset);
      virtual void beginRun(const art::Run& run);
      virtual void analyze(const art::Event& event);
    };

    //================================================================
    EMFDetMCHistPatRec::EMFDetMCHistPatRec(const fhicl::ParameterSet& pset)
      : art::EDAnalyzer(pset)
      , patRecModuleLabel_(pset.get<std::string>("patRecModuleLabel"))
      , patRecInstanceName_(pset.get<std::string>("patRecInstanceName", ""))
      , trkTruthModuleLabel_(pset.get<std::string>("trkTruthModuleLabel"))
      , trkTruthInstanceName_(pset.get<std::string>("trkTruthInstanceName", ""))
      , clusterTruthModuleLabel_(pset.get<std::string>("clusterTruthModuleLabel"))
      , clusterTruthInstanceName_(pset.get<std::string>("clusterTruthInstanceName", ""))
      , particleModuleLabel_(pset.get<std::string>("particleModuleLabel"))
      , particleInstanceName_(pset.get<std::string>("particleInstanceName", ""))

      , geomModuleLabel_(pset.get<std::string>("geomModuleLabel"))
      , geomInstanceName_(pset.get<std::string>("geomInstanceName", ""))

        // Signal particle cuts
      , cutParticleMinClusters_(pset.get<unsigned>("cutParticleMinClusters"))

      , cutHitXmax_(pset.get<double>("cutHitXmax"))
      , cutHitYmax_(pset.get<double>("cutHitYmax"))

      , extmon_()
      , extrapolator_(extmon_)

      , hMultiplicitySignal_()
      , hCommonClusters_()
      , hRecoTruthRatio_()
      , hTruthRecoRatio_()
      , hRecoVsTruth_()
      , effPhysics_(pset.get<unsigned>("cutMinCommonClusters"))
      , effSoftware_(pset.get<unsigned>("cutMinCommonClusters"))
      , fakes_(pset.get<unsigned>("cutMinCommonClusters"))


      , singleParticleMode_(pset.get<bool>("singleParticleMode", false))
      , clusterModuleLabel_(singleParticleMode_ ? pset.get<std::string>("singleParticleClusterModuleLabel") : "")
      , clusterInstanceName_(pset.get<std::string>("singleParticleClusterInstanceName", ""))
      , printEventsWithFakes_(pset.get<bool>("printEventsWithFakes", false))
      , printInefficiencies_(pset.get<bool>("printInefficiencies", false))
    {
      std::cout<<"EMFDetMCHistPatRec: working in the singleParticleMode = "<<singleParticleMode_
               <<", printEventsWithFakes = "<<printEventsWithFakes_
               <<", printInefficiencies = "<<printInefficiencies_
               <<std::endl;
    }

    //================================================================
    void EMFDetMCHistPatRec::beginRun(const art::Run& run) {
      // This is a workaround for geometry not being available at beginJob()
      if(!extmon_) {
        if(!geomModuleLabel_.empty()) {
          art::Handle<ExtMonFNAL::ExtMon> emf;
          run.getByLabel(geomModuleLabel_, geomInstanceName_, emf);
          extmon_ = &*emf;
        }
        else {
          GeomHandle<ExtMonFNAL::ExtMon> emf;
          extmon_ = &*emf;
        }

        extrapolator_ = TrackExtrapolator(extmon_);

        //----------------------------------------------------------------
        art::ServiceHandle<art::TFileService> tfs;

        hMultiplicitySignal_ = tfs->make<TH2D>("multiplicitySignal", "Num PatRec tracks vs num signal SimParticles",
                                               500, -0.5, 499.5, 500, -0.5, 499.5);

        hMultiplicitySignal_->SetOption("colz");
        hMultiplicitySignal_->GetXaxis()->SetTitle("num signal particles");
        hMultiplicitySignal_->GetYaxis()->SetTitle("num PatRec tracks");


        hCommonClusters_ = tfs->make<TH2D>("commonClusters", "Track&SimParticle vs SimParticle clusters for best track",
                                           10, -0.5, 9.5, 10, -0.5, 9.5);

        hCommonClusters_->SetOption("colz");
        hCommonClusters_->GetXaxis()->SetTitle("particle clusters");
        hCommonClusters_->GetYaxis()->SetTitle("common clusters");

        hRecoTruthRatio_ = tfs->make<TProfile>("rtr", "nReco/nTruth vs nTruth", 500, 0.5, 500.5);
        hTruthRecoRatio_ = tfs->make<TProfile>("trt", "nTruth/nReco vs nReco",  500, 0.5, 500.5);

        hRecoVsTruth_ =  tfs->make<TH2D>("recoVsTruth", "nReco vs nTruth", 500, 0.5, 500.5, 500, 0.5, 500.5);
        hRecoVsTruth_->SetOption("colz");

        effPhysics_.book(*extmon_, "effPhysics");
        effSoftware_.book(*extmon_, "effSoftware");
        fakes_.book(*extmon_, "fakes");
      }
    }

    //================================================================
    void EMFDetMCHistPatRec::analyze(const art::Event& event) {

      fillRTR(event);

      if(singleParticleMode_ && !acceptSingleParticleEvent(event)) {
        return;
      }

      art::Handle<SimParticleCollection> ih;
      event.getByLabel(particleModuleLabel_, particleInstanceName_, ih);

      // FIXME: need to create sequence of SimParticles by hand
      // because map_vector does not work with FindMany
      // https://cdcvs.fnal.gov/redmine/issues/2967
      std::vector<art::Ptr<SimParticle> > particles;
      for(SimParticleCollection::const_iterator i = ih->begin(), iend = ih->end(); i != iend; ++i) {
        particles.push_back(art::Ptr<SimParticle>(ih, i->first.asUint()));
      }

      art::FindMany<ExtMonFNALTrkFit,ExtMonFNALTrkMatchInfo>
        trackFinder(particles, event, art::InputTag(trkTruthModuleLabel_, trkTruthInstanceName_));

      art::FindMany<ExtMonFNALRecoCluster,ExtMonFNALRecoClusterTruthBits>
        clusterFinder(particles, event, art::InputTag(clusterTruthModuleLabel_, clusterTruthInstanceName_));

      // different denominator definitions
      std::set<unsigned> signalPhysics;
      std::set<unsigned> signalSW;

      for(unsigned ip=0; ip<particles.size(); ++ip) {
        if(signalParticlePhysics(*particles[ip])) {
          signalPhysics.insert(ip);

          //----------------
          if(true) { // fill an extra histogam
            const std::vector<const ExtMonFNALTrkMatchInfo*>& matchInfo = trackFinder.data(ip);
            // Figure out the best match
            unsigned maxCommonClusters(0), numClustersOnParticle(0);
            for(unsigned itrack = 0; itrack < matchInfo.size(); ++itrack) {
              if(maxCommonClusters < matchInfo[itrack]->nCommonClusters()) {
                maxCommonClusters =  matchInfo[itrack]->nCommonClusters();
                numClustersOnParticle = matchInfo[itrack]->nParticleClusters();
              }
            }
            hCommonClusters_->Fill(numClustersOnParticle, maxCommonClusters);
          }

          //----------------
          if(signalParticleSofware(clusterFinder, ip)) {
            signalSW.insert(ip);
          }

        } // if(physics)
      } // for(ip)

      art::Handle<ExtMonFNALTrkFitCollection> tracks;
      event.getByLabel(patRecModuleLabel_, patRecInstanceName_, tracks);
      hMultiplicitySignal_->Fill(signalPhysics.size(), tracks->size());

      //----------------------------------------------------------------
      EMFPatRecEffHistograms::Fillable phys =
        effPhysics_.fillable(particles,
                             event,
                             art::InputTag(trkTruthModuleLabel_, trkTruthInstanceName_),
                             signalPhysics.size());

      for(std::set<unsigned>::const_iterator i=signalPhysics.begin(); i != signalPhysics.end(); ++i) {
        phys.fill(*i);
      }

      //----------------------------------------------------------------
      EMFPatRecEffHistograms::Fillable sw =
        effSoftware_.fillable(particles,
                              event,
                              art::InputTag(trkTruthModuleLabel_, trkTruthInstanceName_),
                              // Use the same X axis for both physics and SW plots
                              signalPhysics.size());

      for(std::set<unsigned>::const_iterator i=signalSW.begin(); i != signalSW.end(); ++i) {
        if(!sw.fill(*i)) {
          if(printInefficiencies_) {
            std::cout<<"Event "<<event.id()
                     <<" got inefficiency: particle id = "
                     <<particles[*i]->id()
                     <<", pdgId = "<<particles[*i]->pdgId()
                     <<", startPos = "<<particles[*i]->startPosition()
                     <<", startMom = "<<particles[*i]->startMomentum()
                     <<std::endl;
          }
        }
      }

      //----------------------------------------------------------------
      EMFPatRecFakeHistograms::Fillable fk =
        fakes_.fillable(tracks,
                        event,
                        art::InputTag(trkTruthModuleLabel_, trkTruthInstanceName_),
                        // Use the same X axis for both physics and SW plots
                        signalPhysics.size());

      for(unsigned i=0; i < tracks->size(); ++i) {
        if(fk.fill(i)) {
          if(printEventsWithFakes_) {
            std::cout<<"EMFPatRecFakeHistograms flagged event "<<event.id()<<std::endl;
          }
        }
      }
    } // analyze()

    //================================================================
    bool EMFDetMCHistPatRec::signalParticlePhysics(const SimParticle& particle) {
      bool res = false;

      const CLHEP::Hep3Vector& startPos = extmon_->mu2eToExtMon_position(particle.startPosition());
      const CLHEP::Hep3Vector& startMom = extmon_->mu2eToExtMon_momentum(particle.startMomentum());

      // Make sure tha particle starts at a place where it can go through the whole detector
      if((extmon_->up().plane_zoffset().back() < startPos.z()) && /* FIXME: plane_zoffset plus a sensor offset */
         (startMom.z() < -1.) // <= FIXME: zero pz problematic for extrapolation
         ) {

        const double rTrack = extmon_->spectrometerMagnet().trackBendRadius(startMom.mag());
        ExtMonFNALTrkParam mcpar;
        mcpar.setz0(startPos.z());
        mcpar.setposx(startPos.x());
        mcpar.setposy(startPos.y());
        mcpar.setslopex(startMom.x()/startMom.z());
        mcpar.setslopey(startMom.y()/startMom.z());
        mcpar.setrinv(1./rTrack);

        res =  extrapolator_.extrapolateToPlane(extmon_->nplanes()-1, &mcpar) &&
          inAcceptance(mcpar) &&
          extrapolator_.extrapolateToPlane(0, &mcpar) &&
          inAcceptance(mcpar)
          ;
      }

      return res;
    }

    //================================================================
    bool EMFDetMCHistPatRec::inAcceptance(const ExtMonFNALTrkParam& par) {
      return
        (std::abs(par.posx()) < cutHitXmax_) &&
        (std::abs(par.posy()) < cutHitYmax_);
    }


    //================================================================
    bool EMFDetMCHistPatRec::signalParticleSofware(const art::FindMany<ExtMonFNALRecoCluster,ExtMonFNALRecoClusterTruthBits>& clusterFinder,
                                                   unsigned iParticle)
    {
      std::set<unsigned> hitPlanes;
      typedef std::vector<const ExtMonFNALRecoCluster*> Clusters;
      const Clusters& clusters = clusterFinder.at(iParticle);
      for(Clusters::const_iterator i=clusters.begin(); i!=clusters.end(); ++i) {
        hitPlanes.insert( (*i)->plane());
      }
      return hitPlanes.size() == extmon_->nplanes();
    }


    //================================================================
    bool EMFDetMCHistPatRec::acceptSingleParticleEvent(const art::Event& event) {
      art::Handle<ExtMonFNALRecoClusterCollection> coll;
      event.getByLabel(clusterModuleLabel_, clusterInstanceName_, coll);
      return perfectSingleParticleEvent(*coll, extmon_->nplanes());
    }


    //================================================================
    void EMFDetMCHistPatRec::fillRTR(const art::Event& event) {

      // Count the number of primaries
      int nPrimaries(0);
      art::Handle<SimParticleCollection> mch;
      event.getByLabel(particleModuleLabel_, particleInstanceName_, mch);
      for(SimParticleCollection::const_iterator i=mch->begin(), iend = mch->end(); i != iend; ++i) {
        if(!i->second.genParticle().isNull()) {
          ++nPrimaries;
        }
      }

      // Get the number of tracks
      art::Handle<ExtMonFNALTrkFitCollection> tracks;
      event.getByLabel(patRecModuleLabel_, patRecInstanceName_, tracks);
      const int nTracks = tracks->size();

      // Fill the histos
      hRecoVsTruth_->Fill(nPrimaries, nTracks);

      if(nPrimaries > 0) {
        hRecoTruthRatio_->Fill(nPrimaries, double(nTracks)/double(nPrimaries));
      }

      if(nTracks>0) {
        hTruthRecoRatio_->Fill(nTracks, double(nPrimaries)/double(nTracks));
      }

    } // fillRTR()

    //================================================================
  } // namespace ExtMonFNAL
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::ExtMonFNAL::EMFDetMCHistPatRec);
