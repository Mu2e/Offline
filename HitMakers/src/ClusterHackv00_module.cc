//
// A hack at makeing a cluster finder driven from CrudeTrackerHit objects.
//
// $Id: ClusterHackv00_module.cc,v 1.3 2011/05/18 02:27:16 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:16 $
//
// Original author Rob Kutschke
//

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Persistency/Common/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Root includes.
#include "TFile.h"
#include "TH1F.h"

// Mu2e includes.
#include "ToyDP/inc/ProtoStrawCluster.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "LTrackerGeom/inc/LTracker.hh"
#include "LTrackerGeom/inc/CrudeStrawHitCollection.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"
#include "Mu2eUtilities/inc/resolveDPIndices.hh"
#include "HitMakers/inc/growCluster.hh"

using namespace std;

namespace mu2e {

  //--------------------------------------------------------------------
  //
  //
  class ClusterHackv00 : public art::EDAnalyzer {
  public:
    explicit ClusterHackv00(fhicl::ParameterSet const& pset):
      _diagLevel(pset.get<int>("diagLevel",0)),
      _maxFullPrint(pset.get<int>("maxFullPrint",5)),
      _hReHit(0),
      _hNClusters(0),
      _hClusterSize(0),
      _hNCombo(0){
    }
    virtual ~ClusterHackv00() { }

    virtual void beginJob();

    void analyze( art::Event const& e);

  private:

    // Diagnostics level.
    int _diagLevel;

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;

    // Some diagnostic histograms.

    // Counts how often one straw is hit multiple times in one event.
    TH1F* _hReHit;

    // Number of clusters per event.
    TH1F* _hNClusters;

    // Cluster size.
    TH1F* _hClusterSize;

    // Number of triplet combinations per event.
    TH1F* _hNCombo;

  };

  void ClusterHackv00::beginJob(){

    // Create histograms if diagnostics are enabled.
    if ( _diagLevel > 0 ){

      art::ServiceHandle<art::TFileService> tfs;

      _hReHit       = tfs->make<TH1F>( "hReHit",       "Number of multiply Hit Straws in one Event;", 100,  0.,   3.  );
      _hNClusters   = tfs->make<TH1F>( "hNClusters",   "Number of clusters per Event;;",  40,  0.,   40.  );
      _hClusterSize = tfs->make<TH1F>( "hClusterSize", "Cluster size;;",                  20,  0.,   20.  );
      _hNCombo      = tfs->make<TH1F>( "hNCombo",      "Number of triplets per Event;;", 100,  0.,  100.  );

    }

  }

  void
  ClusterHackv00::analyze(art::Event const& evt) {

    static int ncalls(0);
    ++ncalls;

    // Instance name of the module that created the hits of interest;
    static const string creatorName("makeCSH");

    // Geometry info for the LTracker.
    GeomHandle<LTracker> ltrackerHandle;
    LTracker const& ltracker(*ltrackerHandle);

    // Get the persistent data about the CrudeStrawHits.
    art::Handle<CrudeStrawHitPData> pdataHandle;
    evt.getByLabel(creatorName,pdataHandle);

    // Form a fully functional collection of these hits.
    CrudeStrawHitCollection crudeHits(evt, pdataHandle);

    // A convenient alias.  This checks for the validity of the handle.
    CrudeStrawHitPData const& pdata(*pdataHandle);

    // Number of hits.
    int const nhits(pdata.size());

    // Number of straws in the LTracker.
    int const nstraws(ltracker.getAllStraws().size());

    // Set to 1 when this hit has been used in a cluster.
    vector<int> used(nhits,0);

    // The list of clusters.
    vector<ProtoStrawCluster> clusterList;

    // Loop over all hits.
    for ( int i=0; i<nhits; ++i){

      // Skip hits that are already used in a cluster.
      if ( used[i] == 1 ) continue;

      // Aliases for readability.
      CrudeStrawHit const&  hit(crudeHits.get(i));
      Straw const&          straw(ltracker.getStraw(hit.strawIndex));

      // Start a new cluster containing just this straw.
      clusterList.push_back( ProtoStrawCluster( straw.Id().getSectorId(), i ) );
      ProtoStrawCluster& cluster = clusterList.back();

      // Mark this hit as used.
      used[i] = 1;

      // Iteratively grow the cluster until no new hits are added.
      int last(-1);
      while (true){

        int start(last+1);
        last = cluster.size()-1;
        int const startingHit(i+1);

        int nadded = growCluster( cluster, start, startingHit, pdataHandle, used, ltracker );

        if ( nadded == 0 ) break;

      }

    } // end of main loop over hits.


    // Number of found clusters.
    int nClusters(clusterList.size());

    // Diagnostics
    if ( _diagLevel > 0 ){

      _hNClusters->Fill( nClusters);

      if ( _diagLevel > 1 ){
        cout << "\nStarting diagnsotics for evevt: " << evt.id() << endl;
      }

      // Sum of cluster sizes.
      int sum(0);

      for ( int i=0; i<nClusters; ++i){
        _hClusterSize->Fill( clusterList[i].size() );
        sum += clusterList[i].size();

        if ( _diagLevel > 1 ){
          ProtoStrawCluster const& cluster = clusterList[i];
          cout << "Cluster contents: "
               << i << " : ";
          for ( int j=0; j<cluster.size(); ++j ){
            CrudeStrawHit const& hit(pdata[cluster.at(j)]);
            Straw const&  straw(ltracker.getStraw(hit.strawIndex));
            cout << " " << straw.Id();
          }
          cout << endl;
        }
      }

      if ( _diagLevel > 1 ){
        cout << "Check size: "
             << sum << " "
             << pdata.size()
             << endl;
      }

    }

    // Minimum cluster size for forming triplets.
    static int const minClusterSize(3);

    // Number of triplet combinations of clusters in this event.
    int nCombo(0);

    // Form all combinations of three clusters which:
    //  - have a minimum size
    //  - are from 3 different sectors
    for ( int i=0; i<nClusters-2; ++i){
      ProtoStrawCluster const& ci( clusterList.at(i) );
      if ( ci.size() < minClusterSize ) continue;

      for ( int j=i+1; j<nClusters-1; ++j){
        ProtoStrawCluster const& cj( clusterList.at(j) );
        if ( cj.size() < minClusterSize ) continue;

        // Skip a combination if two clusters are in the same sector.
        if ( ci.id == cj.id ) continue;

        for ( int k=j+1; k<nClusters; ++k){
          ProtoStrawCluster const& ck( clusterList.at(k) );
          if ( ck.size() < minClusterSize ) continue;

          // Skip a combination if two clusters are in the same sector.
          if ( ci.id == ck.id  || cj.id == ck.id ) continue;

          // We like this combination.
          ++nCombo;

          if ( _diagLevel > 1 ){
            cout << "Combo: "
                 << evt.id().event() << " "
                 << nCombo << " : ("
                 << ci.id << "),  ("
                 << cj.id << "),  ("
                 << ck.id << ") "
                 << endl;
          }

        } // end of k loop
      }   // end of j loop
    }     // end of i loop

    if ( _diagLevel > 0 ) {
      _hNCombo->Fill(nCombo);

      // Count how often each straw appears in this hit list.
      vector<int> strawHit(nstraws,0);
      for ( int i=0; i<nhits; ++i){
        CrudeStrawHit const&  hit(crudeHits.get(i));
        ++strawHit[hit.strawIndex.asInt()];
      }

      for ( int i=0; i<nstraws; ++i ){
        if ( strawHit[i] > 1 ){
          _hReHit->Fill(strawHit[i]);
        }
      }
    }

  } // end of ::analyze.

}


using mu2e::ClusterHackv00;
DEFINE_ART_MODULE(ClusterHackv00);
