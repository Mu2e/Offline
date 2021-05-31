// Select hits from an StepPointMCCollection-s based on momentum 
// and position and copy selected hits to output collections, 
// preserving product instance names.
//
// Zhengyun You, 2013

#include <string>
#include <vector>
#include <set>
#include <map>
#include <algorithm>

// art includes.
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "DataProducts/inc/PDGCode.hh"

namespace mu2e {

  //================================================================
  class FilterStepPointPositionMomentum : public art::EDFilter {
    typedef std::vector<art::InputTag> InputTags;
    InputTags inputs_;

    // Output instance names.
    typedef std::set<std::string> OutputNames;
    OutputNames outNames_;

    // Particle type lists are typically short, so a vector works
    // better than a set.
    std::vector<PDGCode::type> pdgToDrop_;
    std::vector<PDGCode::type> pdgToKeep_;

    double xMin_;
    double xMax_;
    double yMin_;
    double yMax_;
    double zMin_;
    double zMax_;

    double pxMin_;
    double pxMax_;
    double pyMin_;
    double pyMax_;
    double pzMin_;
    double pzMax_;  
    double pMin_;
    double pMax_;

    // statistics
    unsigned numInputHits_;
    unsigned numOutputHits_;

    unsigned numInputEvents_;
    unsigned numPassedEvents_;

  public:
    explicit FilterStepPointPositionMomentum(const fhicl::ParameterSet& pset);
    virtual bool filter(art::Event& event) override;
    virtual void endJob() override;
  };

  //================================================================
  FilterStepPointPositionMomentum::FilterStepPointPositionMomentum(const fhicl::ParameterSet& pset)
    : art::EDFilter{pset}
    , numInputHits_()
    , numOutputHits_()
    , numInputEvents_()
    , numPassedEvents_()
  {
    const auto tags(pset.get<std::vector<std::string> >("inputs"));
    for(const auto& i : tags) {
      inputs_.emplace_back(i);
      // Coalesce same instance names from multiple input modules/processes.
      outNames_.insert(inputs_.back().instance());
    }
    for(const auto& i : outNames_) {
      produces<StepPointMCCollection>(i);
    }

    const auto drop(pset.get<std::vector<int> >("pdgToDrop", std::vector<int>()));
    for(const auto i : drop) {
      pdgToDrop_.emplace_back(PDGCode::type(i));
    }

    const auto keep(pset.get<std::vector<int> >("pdgToKeep", std::vector<int>()));
    for(const auto i : keep) {
      pdgToKeep_.emplace_back(PDGCode::type(i));
    }

    if(drop.empty() && keep.empty()) {
      throw cet::exception("BADCONFIG")
        <<"FilterStepPointPositionMomentum: either pdgToDrop or pdgToKeep must be specified.\n";
    }
    if(!drop.empty() && ! keep.empty()) {
      throw cet::exception("BADCONFIG")
        <<"FilterStepPointPositionMomentum: either pdgToDrop or pdgToKeep,"
        <<" but not both, can be used at a time.\n";
    }

    const double kMin = -9e99;
    const double kMax =  9e99;

    xMin_ = pset.get<double>("xMin", kMin);
    xMax_ = pset.get<double>("xMax", kMax);
    yMin_ = pset.get<double>("yMin", kMin);
    yMax_ = pset.get<double>("yMax", kMax);
    zMin_ = pset.get<double>("zMin", kMin);
    zMax_ = pset.get<double>("zMax", kMax);

    pxMin_ = pset.get<double>("pxMin", kMin);
    pxMax_ = pset.get<double>("pxMax", kMax);
    pyMin_ = pset.get<double>("pyMin", kMin);
    pyMax_ = pset.get<double>("pyMax", kMax);    
    pzMin_ = pset.get<double>("pzMin", kMin);
    pzMax_ = pset.get<double>("pzMax", kMax);
    pMin_  = pset.get<double>("pMin",  kMin);
    pMax_  = pset.get<double>("pMax",  kMax);
  }

  //================================================================
  bool FilterStepPointPositionMomentum::filter(art::Event& event) {
    bool passed = false;

    typedef std::map<std::string, std::unique_ptr<StepPointMCCollection> > OutMap;
    OutMap outHits;
    for(const auto& i : outNames_) {
      std::unique_ptr<StepPointMCCollection> p(new StepPointMCCollection());
      outHits.insert(std::move(std::make_pair(i, std::move(p))));
    }


    for(const auto& tag : inputs_) {

      auto ih = event.getValidHandle<StepPointMCCollection>(tag);
      StepPointMCCollection& out = *outHits[tag.instance()];

      for(const auto& hit : *ih) {

        const double x = hit.position().x();
        const double y = hit.position().y();
        const double z = hit.position().z();

        const double px = hit.momentum().x();
        const double py = hit.momentum().y();
        const double pz = hit.momentum().z();
        const double p  = hit.momentum().mag();

        if ( x <  xMin_ ||  x >  xMax_ ||
             y <  yMin_ ||  y >  yMax_ ||
             z <  zMin_ ||  z >  zMax_ ||
            px < pxMin_ || px > pxMax_ ||
            py < pyMin_ || py > pyMax_ ||
            pz < pzMin_ || pz > pzMax_ || 
            p  < pMin_  || p  > pMax_  ) {
          continue;
        }

        const PDGCode::type pdgId = hit.simParticle()->pdgId();

        if(!pdgToDrop_.empty() &&
           (std::find(pdgToDrop_.begin(), pdgToDrop_.end(), pdgId) == pdgToDrop_.end()))
          {
            out.emplace_back(hit);
            passed |= true;
          }

        if(!pdgToKeep_.empty() &&
           (std::find(pdgToKeep_.begin(), pdgToKeep_.end(), pdgId) != pdgToKeep_.end()))
          {
            out.emplace_back(hit);
            passed |= true;
          }

      }

      numInputHits_ += ih->size();
      numOutputHits_ += out.size();
    }

    for(const auto& i : outNames_) {
      event.put(std::move(outHits[i]), i);
    }

    ++numInputEvents_;
    if(passed) {
      ++numPassedEvents_;
    }

    return passed;
  }

  //================================================================
  void FilterStepPointPositionMomentum::endJob() {
    mf::LogInfo("Summary")<<"FilterStepPointPositionMomentum: passed "
                          <<numOutputHits_ <<" / "<<numInputHits_
                          <<" StepPointMCs, "
                          <<numPassedEvents_<<" / "<<numInputEvents_
                          <<" events\n";
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::FilterStepPointPositionMomentum);
