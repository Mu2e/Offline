// ExtMonFNAL PatRec efficiency/fake rate analysis
//
// Andrei Gaponenko, 2012

#include <iostream>
#include <string>
#include <vector>

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/FindMany.h"
#include "art/Utilities/InputTag.h"

#include "RecoDataProducts/inc/ExtMonFNALTrkParam.hh"
#include "RecoDataProducts/inc/ExtMonFNALTrkParamCollection.hh"
#include "RecoDataProducts/inc/ExtMonFNALPatRecTrackAssns.hh"

#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/ExtMonFNALRecoClusterTruthAssn.hh"
#include "MCDataProducts/inc/ExtMonFNALPatRecTruthAssns.hh"

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "GeometryService/inc/GeomHandle.hh"

#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "TH1D.h"
#include "TH2D.h"

namespace mu2e {

  //================================================================
  class EMFDetHistPatRec : public art::EDAnalyzer {
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

    unsigned cutParticleMinClusters_;
    double cutParticleMaxAngle_;

    const ExtMonFNAL::ExtMon *extmon_;

    TH2D *hMultiplicityAll_;
    TH2D *hMultiplicitySignal_;

    TH2D *hCommonClusters_;

  public:
    explicit EMFDetHistPatRec(const fhicl::ParameterSet& pset);
    virtual void beginRun(const art::Run& run);
    virtual void analyze(const art::Event& event);
  };

  //================================================================
  EMFDetHistPatRec::EMFDetHistPatRec(const fhicl::ParameterSet& pset)
    : patRecModuleLabel_(pset.get<std::string>("patRecModuleLabel"))
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
    , cutParticleMaxAngle_(pset.get<double>("cutParticleMaxAngle"))

    , extmon_()

    , hMultiplicityAll_()
    , hMultiplicitySignal_()
    , hCommonClusters_()
  {}

  //================================================================
  void EMFDetHistPatRec::beginRun(const art::Run& run) {
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
    art::ServiceHandle<art::TFileService> tfs;

    hMultiplicityAll_ = tfs->make<TH2D>("multiplicityAll", "Num PatRec tracks vs num SimParticles with hits",
                                        200, -0.5, 199.5, 200, -0.5, 199.5);

    hMultiplicityAll_->SetOption("colz");
    hMultiplicityAll_->GetXaxis()->SetTitle("num particles");
    hMultiplicityAll_->GetYaxis()->SetTitle("num PatRec tracks");


    hMultiplicitySignal_ = tfs->make<TH2D>("multiplicitySignal", "Num PatRec tracks vs num signal SimParticles",
                                           200, -0.5, 199.5, 200, -0.5, 199.5);

    hMultiplicitySignal_->SetOption("colz");
    hMultiplicitySignal_->GetXaxis()->SetTitle("num signal particles");
    hMultiplicitySignal_->GetYaxis()->SetTitle("num PatRec tracks");


    hCommonClusters_ = tfs->make<TH2D>("commonClusters", "Track&SimParticle vs SimParticle clusters for best track",
                                       10, -0.5, 9.5, 10, -0.5, 9.5);

    hCommonClusters_->SetOption("colz");
    hCommonClusters_->GetXaxis()->SetTitle("particle clusters");
    hCommonClusters_->GetYaxis()->SetTitle("common clusters");
  }

  //================================================================
  void EMFDetHistPatRec::analyze(const art::Event& event) {


    art::Handle<SimParticleCollection> ih;
    event.getByLabel(particleModuleLabel_, particleInstanceName_, ih);

    // FIXME: need to create sequence of SimParticles by hand
    // because map_vector does not work with FindMany
    // https://cdcvs.fnal.gov/redmine/issues/2967
    std::vector<art::Ptr<SimParticle> > particles;
    for(SimParticleCollection::const_iterator i = ih->begin(), iend = ih->end(); i != iend; ++i) {
      particles.push_back(art::Ptr<SimParticle>(ih, i->first.asUint()));
    }

    art::FindMany<ExtMonFNALTrkParam,ExtMonFNALTrkMatchInfo>
      trackFinder(particles, event, art::InputTag(trkTruthModuleLabel_, trkTruthInstanceName_));

    art::FindMany<ExtMonFNALRecoCluster,ExtMonFNALRecoClusterTruthBits>
      clusterFinder(particles, event, art::InputTag(clusterTruthModuleLabel_, clusterTruthInstanceName_));

    unsigned numParticlesWithHits(0), numSignalParticles(0);
    for(unsigned ip=0; ip<particles.size(); ++ip) {

      const unsigned numClustersOnParticle = clusterFinder.at(ip).size();

      if(numClustersOnParticle > 0) {
        ++numParticlesWithHits;
      }

      if(numClustersOnParticle >= cutParticleMinClusters_) {

        const CLHEP::Hep3Vector mom = extmon_->mu2eToExtMon_momentum(particles[ip]->startMomentum());

        if( (-mom).theta() < cutParticleMaxAngle_) {

          ++numSignalParticles;

          const std::vector<const ExtMonFNALTrkMatchInfo*>& matchInfo = trackFinder.data(ip);

          // Figure out the best match
          unsigned maxCommonClusters(0);
          for(unsigned itrack = 0; itrack < matchInfo.size(); ++itrack) {
            maxCommonClusters = std::max(maxCommonClusters, matchInfo[itrack]->nCommonClusters());
          }

          hCommonClusters_->Fill(numClustersOnParticle, maxCommonClusters);
        }
      }
    }

    art::Handle<ExtMonFNALTrkParamCollection> tracks;
    event.getByLabel(patRecModuleLabel_, patRecInstanceName_, tracks);

    hMultiplicityAll_->Fill(numParticlesWithHits, tracks->size());
    hMultiplicitySignal_->Fill(numSignalParticles, tracks->size());

  } // analyze()

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::EMFDetHistPatRec);
