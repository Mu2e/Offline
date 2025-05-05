
// ======================================================================
//
// WeightSamplingFilterProducer_module:  Allows filtering based on eventweights
//   in order to create new subset of events with correct distribution
//
// ======================================================================

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "art_root_io/TFileService.h"

#include "CLHEP/Random/RandFlat.h"

// ROOT includes
#include "TFile.h"
#include "TH1.h"

// Mu2e includes.
#include "Offline/SeedService/inc/SeedService.hh"
#include "Offline/MCDataProducts/inc/EventWeight.hh"

// c++ includes
#include <iostream>

using namespace std;
namespace mu2e {

  class WeightSamplingFilter : public art::EDFilter {
    public:
      struct Config {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;
        fhicl::Atom<art::InputTag> input{Name("inputTag"), Comment("Module tag creating the event weights")};
        fhicl::Atom<double>        maxWeight{Name("maximumWeight"), Comment("Maximum weight value assumed in accept/reject algorithm"), 1.};
        fhicl::Atom<int>           debug{Name("debugLevel"), Comment("Debugging level"), 0};
        fhicl::Atom<bool>          hists{Name("makeHistograms"), Comment("Make debug histograms"), false};
      };
      using Parameters = art::EDFilter::Table<Config>;
      explicit WeightSamplingFilter(const Parameters& conf);

    private:
      void endJob() override;
      bool filter(art::Event& event) override;

      art::InputTag _evtWtModule;
      double _maxWeight; //accept/reject algorithm needs to know the maximum weight in the input data when sampling
      int _debug;
      bool _hists;
      art::RandomNumberGenerator::base_engine_t& _engine;
      CLHEP::RandFlat _randflat;
      unsigned      _nevt, _npass;
      double        _sumwt;

    struct Hist_t {
      TH1* hwt_ = nullptr;
      TH1* hlogwt_ = nullptr;
      TH1* hpass_ = nullptr;
    };
    enum {kHistSets = 3};
    Hist_t _Hists[kHistSets];
  };

  WeightSamplingFilter::WeightSamplingFilter(const Parameters& config):
    art::EDFilter{config},
    _evtWtModule(config().input()),
    _maxWeight(config().maxWeight()),
    _debug(config().debug()),
    _hists(config().hists()),
    _engine(createEngine( art::ServiceHandle<SeedService>()->getSeed())),
    _randflat( _engine ),
    _nevt(0), _npass(0), _sumwt(0.) {

    // Produces event weights to potentially correct the accept/reject algorithm in cases where the maximum weight is wrong
    produces<mu2e::EventWeight>();

    // Book histograms if requested
    if(_hists) {
      art::ServiceHandle<art::TFileService> tfs;
      {
        art::TFileDirectory tfdir = tfs->mkdir("all_events");
        _Hists[0].hwt_    = tfdir.make<TH1F>("hwt", "Event weight", 200,0.,(_maxWeight > 0.) ? 2.*_maxWeight : 2.);
        _Hists[0].hlogwt_ = tfdir.make<TH1F>("hlogwt", "log10(Event weight)", 200,-10.,1.);
        _Hists[0].hpass_  = tfdir.make<TH1D>("hpass", "Filter result", 2,-0.5,1.5);
      }
      {
        art::TFileDirectory tfdir = tfs->mkdir("weighted_events");
        _Hists[1].hwt_    = tfdir.make<TH1F>("hwt", "Event weight", 200,0.,(_maxWeight > 0.) ? 2.*_maxWeight : 2.);
        _Hists[1].hlogwt_ = tfdir.make<TH1F>("hlogwt", "log10(Event weight)", 200,-10.,1.);
        _Hists[1].hpass_  = tfdir.make<TH1D>("hpass", "Filter result", 2,-0.5,1.5);
      }
      {
        art::TFileDirectory tfdir = tfs->mkdir("accepted_events");
        _Hists[2].hwt_    = tfdir.make<TH1F>("hwt", "Event weight", 200,0.,(_maxWeight > 0.) ? 2.*_maxWeight : 2.);
        _Hists[2].hlogwt_ = tfdir.make<TH1F>("hlogwt", "log10(Event weight)", 200,-10.,1.);
        _Hists[2].hpass_  = tfdir.make<TH1D>("hpass", "Filter result", 2,-0.5,1.5);
      }
    }
  }

  void WeightSamplingFilter::endJob() {
    // Report the filtering efficiency unless printout is actively suppressed
    if(_debug > -1) printf("[WeightSamplingFilter::%s::%s] Saw %u events and accepted %u --> Efficiency = %.3g, sum of weights = %.3g\n",
                           __func__, moduleDescription().moduleLabel().c_str(), _nevt, _npass, (_nevt > 0) ? _npass*1./_nevt : 0., _sumwt);

  }

  bool WeightSamplingFilter::filter(art::Event& event) {
    ++_nevt;
    // Retrieve the event weight
    const double evtwt = event.getValidHandle<EventWeight>( _evtWtModule )->weight();
    if (_debug > 0 && evtwt > _maxWeight){
      printf("[WeightSamplingFilter::%s::%s] Event %i:%i:%u: Event weight is greater than maximum weight (weight = %.3g, max weight = %.3g)\n",
             __func__, moduleDescription().moduleLabel().c_str(), event.run(), event.subRun(), event.event(), evtwt, _maxWeight);
    }
    _sumwt += evtwt;

    // Accept/reject algorithm: Fire a uniform random number and accept if weight is greater than the value
    const double rand = _maxWeight*_randflat.fire();
    const bool pass = evtwt > rand;
    if(pass) ++_npass;
    if(_debug > 1) printf("[WeightSamplingFilter::%s::%s] Event %i:%i:%u: Scaled event weight = %.3g, rand = %.3g, pass = %o\n",
                          __func__, moduleDescription().moduleLabel().c_str(), event.run(), event.subRun(), event.event(), evtwt, rand, pass);

    // Add the output event weight
    const double out_weight = std::max(_maxWeight, evtwt); // if event weight exceeds this maximum, need weights to account for it
    std::unique_ptr<EventWeight> pw(new EventWeight(out_weight));
    event.put(std::move(pw));

    // Fill histograms if requested
    if(_hists) {
      const double logwt = (evtwt > 0.) ? std::log10(evtwt) : -1e10;
      _Hists[0].hwt_   ->Fill(evtwt);
      _Hists[0].hlogwt_->Fill(logwt);
      _Hists[0].hpass_ ->Fill(pass );
      _Hists[1].hwt_   ->Fill(evtwt, evtwt); //weighted to compare to the accepted shapes
      _Hists[1].hlogwt_->Fill(logwt, evtwt);
      _Hists[1].hpass_ ->Fill(pass , evtwt);
      if(pass) {
        _Hists[2].hwt_   ->Fill(evtwt, out_weight);
        _Hists[2].hlogwt_->Fill(logwt, out_weight);
        _Hists[2].hpass_ ->Fill(pass , out_weight);
      }
    }
    return pass;
  }
}

using mu2e::WeightSamplingFilter;
DEFINE_ART_MODULE(WeightSamplingFilter)
