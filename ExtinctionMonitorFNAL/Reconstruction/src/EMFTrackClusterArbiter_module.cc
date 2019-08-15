// Eliminates some tracks from a track fit collection to satisty a
// limit on the number of shared clusters on a track.
//
// $Id: EMFTrackClusterArbiter_module.cc,v 1.3 2013/03/15 15:52:04 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/15 15:52:04 $
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
#include <map>
#include <cassert>

#include "boost/noncopyable.hpp"

//#include "CLHEP/GenericFunctions/CumulativeChiSquare.hh"

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

#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "TH1D.h"
#include "TH2D.h"

//#define AGDEBUG(stuff) do { std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<", func "<<__func__<<": "<<stuff<<std::endl; } while(0)
#define AGDEBUG(stuff)

namespace mu2e {
  namespace ExtMonFNAL {

    //================================================================
    namespace {

      typedef
      std::map<art::Ptr<ExtMonFNALRecoCluster>, std::vector<const ExtMonFNALTrkFit*> >
      ClusterTrackMap;

      //----------------------------------------------------------------
      struct TrkSummary {
        const ExtMonFNALTrkFit *track;
        unsigned numSharedClusters;
        TrkSummary(const ExtMonFNALTrkFit *t, int n) : track(t), numSharedClusters(n) {}
      };

      typedef std::vector<TrkSummary> TrkSummaries;

      //----------------------------------------------------------------
      // Define which of two tracks is "better"
      struct TrkSorter {
        bool operator()(const TrkSummary& a, const TrkSummary& b) {
          return
            (a.numSharedClusters < b.numSharedClusters) ||
            ((a.numSharedClusters == b.numSharedClusters) &&
             // ndf is the same for all tracks, it's sufficient to compare chi2
             (a.track->quality().chi2() < b.track->quality().chi2()));
        }
      };

      //----------------------------------------------------------------
      class TrkSummaryProvider {
        ClusterTrackMap ctm_;
        TrkSummaries summaries_;

      public:

        template<class TrkIter> TrkSummaryProvider(TrkIter begin, TrkIter end);

        const ClusterTrackMap& ctm() const { return ctm_; }

        // Ordered using TrkSorter
        const TrkSummaries& summaries() const { return summaries_; }
      };

      template<class Iter> TrkSummaryProvider::TrkSummaryProvider(Iter begin, Iter end) {

        // fill the cluster map
        for(Iter it = begin; it != end; ++it) {
          for(ExtMonFNALTrkFit::Clusters::const_iterator ic=(*it)->clusters().begin();
              ic != (*it)->clusters().end(); ++ic)
            {
              ctm_[*ic].push_back(&**it);
            }
        }

        // now can compute track summaries
        for(Iter it=begin; it != end; ++it) {

          int numShared = 0;

          for(ExtMonFNALTrkFit::Clusters::const_iterator ic=(*it)->clusters().begin();
              ic != (*it)->clusters().end(); ++ic)
            {
              ClusterTrackMap::const_iterator clusterTracks = ctm_.find(*ic);
              assert(clusterTracks != ctm_.end());
              if(clusterTracks->second.size() > 1) {
                ++numShared;
              }
            }

          summaries_.push_back(TrkSummary(&**it, numShared));
        }

        // Order the tracks
        std::sort(summaries_.begin(), summaries_.end(), TrkSorter());
      }

    }

    //================================================================
    class EMFTrackClusterArbiter : public art::EDProducer {

      int verbosityLevel_;
      std::string tracksModuleLabel_;
      std::string tracksInstanceName_;
      std::string clustersModuleLabel_;
      std::string clustersInstanceName_;

      unsigned cutMaxSharedClustersOnTrack_;
      //// could also have
      //unsigned cutMaxTracksPerCluster_;

      TH1* hTracksPerClusterIn_;
      TH1* hNumSharedPerTrackIn_;
      TH1* hTracksPerClusterOut_;
      TH1* hNumSharedPerTrackOut_;

      void fillTracksPerCluster(TH1 *hh, const art::Handle<ExtMonFNALRecoClusterCollection>& clusters, const ClusterTrackMap& cm);
      void fillNumShared(TH1 *hh, const TrkSummaries& trackSummaries);

    public:
      explicit EMFTrackClusterArbiter(fhicl::ParameterSet const& pset)
        : EDProducer{pset}
        , verbosityLevel_(pset.get<int>("verbosityLevel", 0))
        , tracksModuleLabel_(pset.get<std::string>("tracksModuleLabel"))
        , tracksInstanceName_(pset.get<std::string>("tracksInstanceName", ""))
        , clustersModuleLabel_(pset.get<std::string>("clustersModuleLabel"))
        , clustersInstanceName_(pset.get<std::string>("clustersInstanceName", ""))
        , cutMaxSharedClustersOnTrack_(pset.get<unsigned>("cutMaxSharedClustersOnTrack"))
        , hTracksPerClusterIn_()
        , hNumSharedPerTrackIn_()
        , hTracksPerClusterOut_()
        , hNumSharedPerTrackOut_()
      {
        produces<ExtMonFNALTrkFitCollection>();

        if(clustersModuleLabel_.empty()) {
          std::cout<<"EMFTrackClusterArbiter: clustersModuleLabel not given. Not filling per-cluster histograms."<<std::endl;
        }
      }

      virtual void beginJob();
      virtual void produce(art::Event& evt);
    };

    //================================================================
    void EMFTrackClusterArbiter::beginJob() {
      art::ServiceHandle<art::TFileService> tfs;

      if(!clustersModuleLabel_.empty()) {
        hTracksPerClusterIn_ = tfs->make<TH1D>("tracksPerClusterIn", "Num tracks per cluster, in", 20, -0.5, 19.5);
        hTracksPerClusterOut_ = tfs->make<TH1D>("tracksPerClusterOut", "Num tracks per cluster, out", 20, -0.5, 19.5);
      }

      hNumSharedPerTrackIn_ = tfs->make<TH1D>("sharedClustersPerTrackIn", "Num shared clusters on track, in", 7, -0.5, 6.5);
      hNumSharedPerTrackOut_ = tfs->make<TH1D>("sharedClustersPerTrackOut", "Num shared clusters on track, out", 7, -0.5, 6.5);
    }

    //================================================================
    void EMFTrackClusterArbiter::produce(art::Event& event) {
      art::Handle<ExtMonFNALTrkFitCollection> intracksh;
      event.getByLabel(tracksModuleLabel_, tracksInstanceName_, intracksh);

      art::Handle<ExtMonFNALRecoClusterCollection> clustersh;
      if(!clustersModuleLabel_.empty()) {
        event.getByLabel(clustersModuleLabel_, clustersInstanceName_, clustersh);
      }

      typedef std::set<const ExtMonFNALTrkFit*> TrkSet;
      TrkSet remainingTracks;
      for(ExtMonFNALTrkFitCollection::const_iterator i = intracksh->begin(), iend = intracksh->end(); i != iend; ++i) {
        remainingTracks.insert(&*i);
      }

      TrkSummaryProvider tsp(remainingTracks.begin(), remainingTracks.end());

      // fill the "before" histos
      fillNumShared(hNumSharedPerTrackIn_, tsp.summaries());
      if(hTracksPerClusterIn_) {
        fillTracksPerCluster(hTracksPerClusterIn_, clustersh, tsp.ctm());
      }
      // filter the tracks
      if(!tsp.summaries().empty()) {
        while(tsp.summaries().back().numSharedClusters > cutMaxSharedClustersOnTrack_) {
          const auto status [[gnu::unused]] = remainingTracks.erase(tsp.summaries().back().track);
          assert(status);
          tsp = TrkSummaryProvider(remainingTracks.begin(), remainingTracks.end());
        }
      }

      // Fill the "after" histos
      fillNumShared(hNumSharedPerTrackOut_, tsp.summaries());
      if(hTracksPerClusterOut_) {
        fillTracksPerCluster(hTracksPerClusterOut_, clustersh, tsp.ctm());
      }

      // Prepare the output collection
      std::unique_ptr<ExtMonFNALTrkFitCollection> outtracks(new ExtMonFNALTrkFitCollection);
      for(TrkSet::const_iterator i = remainingTracks.begin(); i!=remainingTracks.end(); ++i) {
        outtracks->push_back(**i);
      }
      event.put(std::move(outtracks));
    }

    //================================================================
    void EMFTrackClusterArbiter::fillTracksPerCluster(TH1* hh,
                                                      const art::Handle<ExtMonFNALRecoClusterCollection>& clusters,
                                                      const ClusterTrackMap& cm) {

      for(unsigned i=0; i<clusters->size(); ++i) {
        art::Ptr<ExtMonFNALRecoCluster> cl(clusters, i);
        ClusterTrackMap::const_iterator entry = cm.find(cl);
        const unsigned numShared = (entry == cm.end()) ? 0 : entry->second.size();
        hh->Fill(numShared);
      }
    }

    //================================================================
    void  EMFTrackClusterArbiter::fillNumShared(TH1 *hh, const TrkSummaries& tsms) {
      for(TrkSummaries::const_iterator is=tsms.begin(); is != tsms.end(); ++is) {
        hh->Fill(is->numSharedClusters);
      }
    }

    //================================================================
  } // end namespace ExtMonFNAL
} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::ExtMonFNAL::EMFTrackClusterArbiter)
