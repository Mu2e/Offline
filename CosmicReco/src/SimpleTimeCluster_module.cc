//Author: R Bonventre 
//Date: Feb 2020
//Purpose: To improve TimeClustering for Cosmics

#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
// Mu2e
#include "GeneralUtilities/inc/Angles.hh"
#include "Mu2eUtilities/inc/MVATools.hh"
#include "Mu2eUtilities/inc/polyAtan2.hh"
// data
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
// tracking
#include "TrkReco/inc/TrkUtilities.hh"
#include "TrkReco/inc/TrkTimeCalculator.hh"
// C++
#include <memory>
#include <algorithm>
#include <utility>
using namespace std;

namespace mu2e {
  class SimpleTimeCluster : public art::EDProducer {
    public:
      enum Mode{flag=0,filter};
      explicit SimpleTimeCluster(fhicl::ParameterSet const& pset);

      void produce(art::Event& e) override;

    private:
      int               _iev;
      int               _debug;
      int               _printfreq;
      bool		_testflag;
      art::ProductToken<ComboHitCollection> const _chToken;
      art::ProductToken<StrawHitFlagCollection> const _shfToken;
      const StrawHitFlagCollection *_shfcol;
      const ComboHitCollection *_chcol;
      StrawHitFlag      _hsel, _hbkg;


      bool goodHit(const StrawHitFlag& flag) const;
  };

  SimpleTimeCluster::SimpleTimeCluster(fhicl::ParameterSet const& pset) :
    art::EDProducer{pset},
    _debug             (pset.get<int>("debugLevel",0)),
    _printfreq         (pset.get<int>("printFrequency",101)),
    _testflag(pset.get<bool>("TestFlag")),
    _chToken{consumes<ComboHitCollection>(pset.get<art::InputTag>("ComboHitCollection"))},
    _shfToken{mayConsume<StrawHitFlagCollection>(pset.get<art::InputTag>("StrawHitFlagCollection"))},
    _hsel              (pset.get<std::vector<std::string> >("HitSelectionBits",vector<string>{"EnergySelection","TimeSelection"})),
    _hbkg              (pset.get<vector<string> >("HitBackgroundBits",vector<string>{}))
    {
      produces<TimeClusterCollection>();
    }


  //--------------------------------------------------------------------------------------------------------------
  void SimpleTimeCluster::produce(art::Event & event ){
    _iev = event.id().event();

    if (_debug > 0 && (_iev%_printfreq)==0) std::cout<<"SimpleTimeCluster: event="<<_iev<<std::endl;

    auto const& chH = event.getValidHandle(_chToken);
    _chcol = chH.product();

    if(_testflag){
      auto shfH = event.getValidHandle(_shfToken);
      _shfcol = shfH.product();
      if(_shfcol->size() != _chcol->size())
        throw cet::exception("RECO")<<"SimpleTimeCluster: inconsistent flag collection length " << endl;
    }

    std::unique_ptr<TimeClusterCollection> tccol(new TimeClusterCollection);

    TimeCluster tclust;
    tclust._nsh = 0;
    double min = 9e9;
    double max = -9e9;
    double avg = 0;
    for (size_t j=0;j<_chcol->size();j++){
      if (_testflag && !goodHit((*_shfcol)[j]))
        continue;
      if (_chcol->at(j).time() < min)
        min = _chcol->at(j).time();
      if (_chcol->at(j).time() > max)
        max = _chcol->at(j).time();
      avg += _chcol->at(j).time();

      tclust._strawHitIdxs.push_back(j);
      tclust._nsh += _chcol->at(j).nStrawHits();
    }
    tclust._t0 = TrkT0(avg/tclust._strawHitIdxs.size(),(max-min)/2.);
    (*tccol).push_back(tclust);


    if (_debug > 0) std::cout << "Found " << tccol->size() << " Time Clusters " << std::endl;

    if (_debug > 1){
      for(auto const& tc : *tccol) {
        std::cout << "Time Cluster time = " << tc.t0().t0() << " +- " << tc.t0().t0Err()
          << " position = " << tc._pos << std::endl;
        if(_debug > 3){
          for (auto shi : tc._strawHitIdxs ) {
            std::cout << "Time Cluster hit at index " << shi << std::endl;
          }
        }
      }
    }
    event.put(std::move(tccol));
  }

  bool SimpleTimeCluster::goodHit(const StrawHitFlag& flag) const
  {
    return flag.hasAllProperties(_hsel) && !flag.hasAnyProperty(_hbkg);
  }

}

using mu2e::SimpleTimeCluster;
DEFINE_ART_MODULE(SimpleTimeCluster);

