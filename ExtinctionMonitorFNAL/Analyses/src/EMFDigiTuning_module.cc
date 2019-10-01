// This module evaluates pixel sensor efficiency and other parameters
// for digitization tuning.  It requres a specially prepared single
// particle file, with protons shot at ExtMonFNAL that should be
// guaranteed to hit the most upstream pixel sensor somewhere not too
// near to the sensor edge.
//
// Andrei Gaponenko, 2012

#include <iostream>
#include <string>
#include <vector>

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/ParameterSetRegistry.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"

#include "RecoDataProducts/inc/ExtMonFNALRawHit.hh"
#include "RecoDataProducts/inc/ExtMonFNALRawHitCollection.hh"
#include "RecoDataProducts/inc/ExtMonFNALRecoCluster.hh"
#include "RecoDataProducts/inc/ExtMonFNALRecoClusterCollection.hh"

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "GeometryService/inc/GeomHandle.hh"

#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "TH1.h"
#include "TH2.h"
#include "TEfficiency.h"

namespace mu2e {
  namespace ExtMonFNAL {

    //================================================================
    class EMFDigiTuning : public art::EDAnalyzer {

      std::string hitsModuleLabel_;
      std::string hitsInstanceName_;

      std::string clustersModuleLabel_;
      std::string clustersInstanceName_;

      std::string geomModuleLabel_;
      std::string geomInstanceName_;

      const ExtMon *extmon_;

      // an artifact of booking in initialization
      art::TFileService *tfs_;

      TEfficiency *sensorEffVsThreshold_;
      TEfficiency *sensorEffVsCurrent_;
      TEfficiency *sensorEffVsNumClusters_;

      TEfficiency *sensorIneffMap_;

      TH1 *numRawHits_;
      TH1 *hitToT_;
      TH1 *hitClockDiff_; // (ti - tfirst) for numRawHits>1

      TH1 *numRecoClusters_;
      TH1 *clusterNumPixels_;
      TH1 *clusterNx_;
      TH1 *clusterNy_;

    public:
      explicit EMFDigiTuning(const fhicl::ParameterSet& pset);
      virtual void beginRun(const art::Run& run);
      virtual void analyze(const art::Event& event);
    };

    //================================================================
    EMFDigiTuning::EMFDigiTuning(const fhicl::ParameterSet& pset)
      : art::EDAnalyzer(pset)
      , hitsModuleLabel_(pset.get<std::string>("hitsModuleLabel"))
      , hitsInstanceName_(pset.get<std::string>("hitsInstanceName", ""))
      , clustersModuleLabel_(pset.get<std::string>("clustersModuleLabel"))
      , clustersInstanceName_(pset.get<std::string>("clustersInstanceName", ""))
      , geomModuleLabel_(pset.get<std::string>("geomModuleLabel"))
      , geomInstanceName_(pset.get<std::string>("geomInstanceName", ""))

      , tfs_(&*art::ServiceHandle<art::TFileService>())

      , sensorEffVsThreshold_(tfs_->make<TEfficiency>("sensorEffVsTh", "Sensor efficiency vs threshold", 20, 250., 10250.))
      , sensorEffVsCurrent_(tfs_->make<TEfficiency>("sensorEffVsI", "Sensor efficiency vs log2(current/1000)", 15, -0.5, 14.5))
      , sensorEffVsNumClusters_(tfs_->make<TEfficiency>("sensorEffVsNc", "Sensor efficiency vs log2(numDigiClusters)", 15, -0.5, 14.5))

      , sensorIneffMap_(tfs_->make<TEfficiency>("sensorIneffMap", "Sensor inefficiency: log2(current/1000) vs threshold",
                                                20, 250., 10250., 15, -0.5, 14.5))

      , numRawHits_(tfs_->make<TH1D>("numRawHits", "num raw hits", 40, -0.5, 39.5))
      , hitToT_(tfs_->make<TH1D>("hitToT", "hit ToT", 15, -0.5, 14.5))
      , hitClockDiff_(tfs_->make<TH1D>("hitClockDiff", "raw hit clock(i)-clock(first)", 10, -0.5, 9.5))
      , numRecoClusters_(tfs_->make<TH1D>("numRecoClusters", "num reco clusters", 20, -0.5, 19.5))
      , clusterNumPixels_(tfs_->make<TH1D>("clusterNumPixels", "cluster num pixels", 40, 0.5, 40.5))
      , clusterNx_(tfs_->make<TH1D>("clusterNx", "cluster ny", 20, 0.5, 20.5))
      , clusterNy_(tfs_->make<TH1D>("clusterNy", "cluster ny", 20, 0.5, 20.5))
    {}

    //================================================================
    void EMFDigiTuning::beginRun(const art::Run& run) {
      if(!geomModuleLabel_.empty()) {
        art::Handle<ExtMonFNAL::ExtMon> emf;
        run.getByLabel(geomModuleLabel_, geomInstanceName_, emf);
        extmon_ = &*emf;
      }
      else {
        GeomHandle<ExtMonFNAL::ExtMon> emf;
        extmon_ = &*emf;
      }
    }

    //================================================================
    void EMFDigiTuning::analyze(const art::Event& event) {
      // only look at this single sensor plane
      const unsigned iplane = extmon_->nplanes() - 1;

      art::Handle<ExtMonFNALRawHitCollection> hitsh;
      event.getByLabel(hitsModuleLabel_, hitsInstanceName_, hitsh);
      const ExtMonFNALRawHitCollection& allHits(*hitsh);

      typedef std::vector<const ExtMonFNALRawHit*> SensorHits;
      SensorHits hits;
      for(ExtMonFNALRawHitCollection::const_iterator i=allHits.begin(); i!=allHits.end(); ++i) {
        if(i->pixelId().chip().module().plane() == iplane) {
          hits.push_back(&*i);
        }
      }

      art::Handle<ExtMonFNALRecoClusterCollection> clustersh;
      event.getByLabel(clustersModuleLabel_, clustersInstanceName_, clustersh);
      typedef ExtMonFNALRecoClusterCollection::PlaneClusters PlaneClusters;
      const PlaneClusters& clusters(clustersh->clusters(iplane));

      if(clusters.empty() && !hits.empty()) {
        throw cet::exception("BUG")<<"No reco clusters for non-empty hits";
      }

      //----------------
      // Data points for the "scan" plots

      {
        //// this needs new art:
        //const fhicl::ParameterSet& digipset(hitsh.provenance()->parameterSet());

        if(hitsh.provenance()->psetIDs().size() != 1) {
          throw cet::exception("BADINPUTS")
            <<"Raw hits provenance: more than on psetIDs"
            <<" --- did you forget to change name in the writer configuration?\n";
        }
        const fhicl::ParameterSet& digiConfig =
          fhicl::ParameterSetRegistry::get(*hitsh.provenance()->psetIDs().begin());

        const double threshold = digiConfig.get<double>("discriminatorThreshold");
        const double qCalib = digiConfig.get<double>("qCalib");
        const int totCalib = digiConfig.get<int>("totCalib");
        const double k = (qCalib - threshold)/(totCalib);
        //std::cout<<"EMFDigiTuning: got digi config."<<" discriminatorThreshold = "<<discriminatorThreshold<<std::endl;
        const double logI = log2(k/1000.);

        const bool eff = !hits.empty();
        sensorEffVsThreshold_->Fill(eff, threshold);
        sensorEffVsCurrent_->Fill(eff, logI);
        sensorEffVsNumClusters_->Fill(eff, log2(digiConfig.get<int>("numClustersPerHit")));
        sensorIneffMap_->Fill(!eff, threshold, logI);
      }

      //----------------
      numRawHits_->Fill(hits.size());

      std::vector<int> hitTimes;
      for(SensorHits::const_iterator i=hits.begin(); i!=hits.end(); ++i) {
        hitToT_->Fill((*i)->tot());
        hitTimes.push_back((*i)->clock());
      }
      std::sort(hitTimes.begin(), hitTimes.end());
      for(unsigned i=1; i<hitTimes.size(); ++i) {
        hitClockDiff_->Fill(hitTimes[i]-hitTimes[0]);
      }

      numRecoClusters_->Fill(clusters.size());
      for(unsigned i = 0; i < clusters.size(); ++i) {
        const ExtMonFNALRecoCluster& cl = clusters[i];
        clusterNumPixels_->Fill(cl.raw()->hits().size());
        clusterNx_->Fill(cl.xWidth());
        clusterNy_->Fill(cl.yWidth());
      }
      //----------------
    }

    //================================================================
  } // namespace ExtMonFNAL
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::ExtMonFNAL::EMFDigiTuning);
