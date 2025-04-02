///////////////////////////////////////////////////////////////////////////////
// takes inputs from two helix finding algorithms, produces one helix collection
// on output to be used for the track seed-fit
//
// also use to merge collections of negative and positive helices. In this case,
// expect algorithm ID colelctions to be present
//
// Original author P. Murat
///////////////////////////////////////////////////////////////////////////////
// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
// BaBar
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/KalmanTrack/KalHit.hh"
#include "BTrk/BbrGeom/TrkLineTraj.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrk/TrkBase/HelixParams.hh"

#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"

#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"

#include "TVector2.h"

#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/RecoDataProducts/inc/HelixHit.hh"

#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"

// #include "Offline/BTrkData/inc/Doublet.hh"
// #include "Offline/TrkReco/inc/DoubletAmbigResolver.hh"

#include "Offline/Mu2eUtilities/inc/LsqSums4.hh"
#include "Offline/CalPatRec/inc/ObjectDumpUtils.hh"

#include "Offline/RecoDataProducts/inc/AlgorithmID.hh"

//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Vector/ThreeVector.h"
// root
#include "TMath.h"
// C++
#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <functional>
#include <float.h>
#include <vector>
#include <set>
#include <map>
using namespace std;
using CLHEP::Hep3Vector;

namespace mu2e {
  class MergeHelixFinder : public art::EDProducer {
  public:
    explicit MergeHelixFinder(fhicl::ParameterSet const&);
    virtual ~MergeHelixFinder();
    virtual void beginJob();
    virtual void beginRun(art::Run&);
    virtual void produce(art::Event& event );
    void endJob();

    void saveBest(HelixSeedCollection* HelixColl, const HelixSeed* Helix,
                  art::Event& AnEvent, const art::Provenance* Prov,
                  int Index, AlgorithmIDCollection* AidColl);

    void checkPresenceOfDuplicates(HelixSeedCollection* HelixColl, std::vector<int> &IndexToSkip);

  private:
    unsigned         _iev;
                                        // configuration parameters
    int              _diag;
    int              _debugLevel;
    float            _minTprChi2;
    float            _minCprChi2;
    int              _printfreq;
    bool             _addhits;
                                        // event object labels
    std::string      _tprHelixCollTag;
    std::string      _cprHelixCollTag;
  };

  MergeHelixFinder::MergeHelixFinder(fhicl::ParameterSet const& pset) :
    EDProducer{pset},
    _diag            (pset.get<int>        ("diagLevel"      )),
    _debugLevel      (pset.get<int>        ("debugLevel"     )),
    _tprHelixCollTag (pset.get<std::string>("tprHelixCollTag")),
    _cprHelixCollTag (pset.get<std::string>("cprHelixCollTag"))
  {

    produces<AlgorithmIDCollection>  ();
    produces<HelixSeedCollection>    ();

  }

  MergeHelixFinder::~MergeHelixFinder() {
  }

  void MergeHelixFinder::beginJob() {
  }

  void MergeHelixFinder::beginRun(art::Run& ) {
  }

//-----------------------------------------------------------------------------
// search within the collection the presence of duplicates
//-----------------------------------------------------------------------------
  void MergeHelixFinder::checkPresenceOfDuplicates(HelixSeedCollection* HelixColl, std::vector<int> &IndexToSkip){
    int    nHel = HelixColl->size();
    int    nh1(0), nh2(0), natc(0);

    const HelixSeed          *h1(0), *h2(0);
    const ComboHitCollection *list1(0), *list2(0);
    const mu2e::HelixHit     *hitt, *hitc;

    for (int i=0; i<nHel; ++i){
      h1    = &HelixColl->at(i);
      list1 = &h1->hits();
      nh1   = list1->size();
      natc  = 0;

      for (int j=i+1; j<nHel; ++j){
        h2    = &HelixColl->at(j);
        list2 = &h2->hits();
        nh2   = list2->size();

        for (int k=0; k<nh1; ++k){
          hitt = &list1->at(k);
          for (int l=0; l<nh2; l++){
            hitc = &list2->at(l);
            if (hitt->index() == hitc->index()) {
              natc += 1;
              break;
            }
          }
        }
//-----------------------------------------------------------------------------
// if > 50% of the helix hits are common, it unlikely to be an "independent" object
// logic of the choice:
// 1. take the track which has more hits
// 2. if two helices have the same number of hits, pending future studies,
//    the choose CalPatRec one
//-----------------------------------------------------------------------------
        if ((natc > nh1/2.) || (natc > nh2/2.)) {
          //-----------------------------------------------------------------------------
          // for one of the two helices, the number of shared hits > 50%
          //-----------------------------------------------------------------------------
          if (nh2 > nh1) {
            //-----------------------------------------------------------------------------
            // h2 is a winner, no need to save h1
            //-----------------------------------------------------------------------------
            IndexToSkip.push_back(i);
            break;
          }
          else {
            //-----------------------------------------------------------------------------
            // h1 is a winner, mark h2 in hope that it will be OK, continue looping
            //-----------------------------------------------------------------------------
            IndexToSkip.push_back(j);
          }
        }

      }//end loop over j
    }//end loop over i
  }

//-----------------------------------------------------------------------------
// algorithm collections are unnamed
//-----------------------------------------------------------------------------
  void MergeHelixFinder::saveBest(HelixSeedCollection* HelixColl, const HelixSeed* Helix,
                                  art::Event& AnEvent, const art::Provenance* Prov,
                                  int Index, AlgorithmIDCollection* AidColl) {

    AlgorithmID                              aid;
    short                                    best(0),  mask(-1);
    art::Handle<mu2e::AlgorithmIDCollection> algH;

    HelixColl->push_back(*Helix);

    AnEvent.getByLabel(Prov->moduleLabel(),"",Prov->processName(),algH);

    if (algH.isValid()) {
      mask = algH->at(Index).AlgMask();
      best = algH->at(Index).BestID();
    }
    else                 {
      if      (_cprHelixCollTag.find("HelixFinder:"  ) == 0) {
        mask = 1 << AlgorithmID::TrkPatRecBit;
        best = AlgorithmID::TrkPatRecBit;
      }
      else if (_cprHelixCollTag.find("CalHelixFinder") == 0) {
        mask = 1 << AlgorithmID::CalPatRecBit;
        best = AlgorithmID::CalPatRecBit;
      }
    }

    aid.Set(best,mask);
    AidColl->push_back(aid);
  }


//-----------------------------------------------------------------------------
  void MergeHelixFinder::produce(art::Event& AnEvent) {

                                        // assume less than 100 tracks

    int nhel[2] {0,0} ;

    art::Handle<mu2e::HelixSeedCollection>    hcH[2];

    mu2e::HelixSeedCollection                 *hcoll[2] {nullptr,nullptr};

    // mu2e::GeomHandle<mu2e::DetectorSystem>    ds;
    // mu2e::GeomHandle<mu2e::VirtualDetector>   vdet;

    unique_ptr<AlgorithmIDCollection>         algs     (new AlgorithmIDCollection );
    unique_ptr<HelixSeedCollection>           helixPtrs(new HelixSeedCollection   );

    if (_debugLevel > 0) ObjectDumpUtils::printEventHeader(&AnEvent,"MergeHelixFinder::produce");
//-----------------------------------------------------------------------------
// internal function
//-----------------------------------------------------------------------------
    art::Handle<mu2e::AlgorithmIDCollection> algH;

    const art::Provenance  *prov[2] {nullptr,nullptr};

    AnEvent.getByLabel(_tprHelixCollTag,hcH[0]);
    AnEvent.getByLabel(_cprHelixCollTag,hcH[1]);

    for (int i=0; i<2; i++) {
      if (hcH[i].isValid()) {
        prov [i] = hcH[i].provenance();
        hcoll[i] = (mu2e::HelixSeedCollection*) &(*hcH[i]);
        nhel [i] = hcoll[i]->size();
      }
    }

    int mnhel = max(nhel[0],nhel[1]);

    vector<int> flag[2];
    flag[0].reserve(mnhel);
    flag[1].reserve(mnhel);

    for (int i=0; i<mnhel; i++) {
      flag[0][i] = 1;
      flag[1][i] = 1;
    }

    std::vector<int>         helixToSkip0, helixToSkip1;

    checkPresenceOfDuplicates(hcoll[0], helixToSkip0);
    checkPresenceOfDuplicates(hcoll[1], helixToSkip1);

    const HelixSeed          *h1, *h2;
    const ComboHitCollection *tlist, *clist;
    int                       nh1, nh2, natc;
    const mu2e::HelixHit     *hitt, *hitc;

    for (int i1=0; i1<nhel[0]; i1++) {
      std::vector<int>::iterator it1;
      it1 = find(helixToSkip0.begin(), helixToSkip0.end(), i1);
      if (it1 != helixToSkip0.end())                  continue;

      h1     = &hcoll[0]->at(i1);
//------------------------------------------------------------------------------
// check if an AlgorithmID collection has been created by the process
//-----------------------------------------------------------------------------
      tlist        = &h1->hits();
      nh1          = tlist->size();
      natc         = 0;

      for (int i2=0; i2<nhel[1]; i2++) {
        std::vector<int>::iterator it2;
        it2 = find(helixToSkip1.begin(), helixToSkip1.end(), i2);
        if (it2 != helixToSkip1.end())                  continue;

        h2 = &hcoll[1]->at(i2);
//-----------------------------------------------------------------------------
// at Mu2e, 2 helices with different helicity could be duplicates of each other
//-----------------------------------------------------------------------------
        clist        = &h2->hits();
        nh2          = clist->size();
//-----------------------------------------------------------------------------
// check the number of common hits: do we need to check also if they have
// close momentum?
//-----------------------------------------------------------------------------
        for (int k=0; k<nh1; ++k){
          hitt = &tlist->at(k);
          for (int l=0; l<nh2; l++){
            hitc = &clist->at(l);
            if (hitt->index() == hitc->index()) {
              natc += 1;
              break;
            }
          }
        }
//-----------------------------------------------------------------------------
// if > 50% of the helix hits are common, it unlikely to be an "independent" object
// logic of the choice:
// 1. take the track which has more hits
// 2. if two helices have the same number of hits, pending future studies,
//    the choose CalPatRec one
//-----------------------------------------------------------------------------
        if ((natc > nh1/2.) || (natc > nh2/2.)) {
//-----------------------------------------------------------------------------
// for one of the two helices, the number of shared hits > 50%
//-----------------------------------------------------------------------------
          if (nh2 > nh1) {
//-----------------------------------------------------------------------------
// h2 is a winner, no need to save h1
//-----------------------------------------------------------------------------
            saveBest(helixPtrs.get(),h2,AnEvent,prov[1],i2,algs.get());

            flag[0][i1] = 0;
            flag[1][i2] = 0;
            break;
          }
          else {
//-----------------------------------------------------------------------------
// h1 is a winner, mark h2 in hope that it will be OK, continue looping
//-----------------------------------------------------------------------------
            flag[1][i2] = 0;
          }
        }
      }

      if (flag[0][i1] == 1) {
//-----------------------------------------------------------------------------
// looped over all helices from the second coll, nothing looking like h1
//-----------------------------------------------------------------------------
        saveBest(helixPtrs.get(),h1,AnEvent,prov[0],i1,algs.get());
      }
    }
//-----------------------------------------------------------------------------
// account for presence of multiple helices - loop over the 2nd collection
// and pass the unique helices to output
//-----------------------------------------------------------------------------
    for (int i=0; i<nhel[1]; i++) {
      if (flag[1][i] == 1) {
        h2 = &hcoll[1]->at(i);
        saveBest(helixPtrs.get(),h2,AnEvent,prov[1],i,algs.get());
      }
    }

    AnEvent.put(std::move(helixPtrs));
    AnEvent.put(std::move(algs     ));
  }
//-----------------------------------------------------------------------------
// end job :
//-----------------------------------------------------------------------------
  void MergeHelixFinder::endJob() {
  }

}

using mu2e::MergeHelixFinder;
DEFINE_ART_MODULE(MergeHelixFinder)
