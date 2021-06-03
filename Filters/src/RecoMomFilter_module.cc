
// ======================================================================
//
// RecoMomFilter_module: allows filtering on the momentum of tracks
//   in KalSeed
//
// ======================================================================

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Principal/Handle.h"

// Mu2e includes.
#include "RecoDataProducts/inc/KalSeed.hh"



#include <iostream>
#include <string>

using namespace std;
namespace mu2e {

  class RecoMomFilter : public art::EDFilter {
    public:
      explicit RecoMomFilter(fhicl::ParameterSet const& pset);

    private:
      bool beginRun(art::Run& run) override;
      bool endRun(art::Run& run) override;
      bool filter(art::Event& event) override;
      
      std::string _moduleRoot;
      std::vector<std::string> _trkTags;
      std::vector<double> _momCutoff;
  };

  RecoMomFilter::RecoMomFilter(fhicl::ParameterSet const& pset):
    art::EDFilter{pset},
    _moduleRoot(pset.get<std::string>("KalFinalTagRoot")),
    _trkTags(pset.get<std::vector<std::string> >("TrkTags")),
    _momCutoff(pset.get<std::vector<double> >("MomentumCutoff"))
    {
      if (_trkTags.size() != _momCutoff.size())
	throw cet::exception("CONFIG")
	  << "RecoMomFilter module: TrkTags length must be the same as MomentumCutoff length.\n";
    }

  bool RecoMomFilter::beginRun(art::Run& run) {
    return true;
  }
  bool RecoMomFilter::endRun(art::Run& run) {
    return true;
  }

  bool RecoMomFilter::filter(art::Event& event) {

    std::vector<art::Handle<KalSeedCollection> > colls = event.getMany<KalSeedCollection>();

    std::vector<std::string> moduleNames;
    for (size_t i=0;i<_trkTags.size();i++){
      moduleNames.push_back(_moduleRoot + _trkTags[i]);
    }

    bool pass = false;
    for(const auto& ic : colls) {
      std::string moduleLabel = ic.provenance()->moduleLabel();
      auto mniter = std::find(moduleNames.begin(),moduleNames.end(),moduleLabel);
      if (mniter == moduleNames.end())
        continue;
      for(const auto& ks : *ic) {
        if (ks.segments().begin()->mom() > _momCutoff[mniter-moduleNames.begin()]){
          pass = true;
          break;
        }
      }
    }

    if (pass)
      return true;

    return false;
  }
}

using mu2e::RecoMomFilter;
DEFINE_ART_MODULE(RecoMomFilter);
