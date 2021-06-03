// Select hits from an StepPointMCCollection-s based on its pz 
// (in target frame) and copy selected hits to output collections, 
// preserving product instance names.
//
// Zhengyun You, 2013

#include <string>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include "CLHEP/Units/SystemOfUnits.h"

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
  class FilterStepPointPzVsTarget : public art::EDFilter {
    typedef std::vector<art::InputTag> InputTags;
    InputTags inputs_;

    // Output instance names.
    typedef std::set<std::string> OutputNames;
    OutputNames outNames_;

    // Particle type lists are typically short, so a vector works
    // better than a set.
    std::vector<PDGCode::type> pdgToDrop_;
    std::vector<PDGCode::type> pdgToKeep_;

    double targetAngle_;
    double targetDirx_;
    double targetDiry_;
    double targetDirz_;    

    double pzMin_;
    double pzMax_;

    // statistics
    unsigned numInputHits_;
    unsigned numOutputHits_;

    unsigned numInputEvents_;
    unsigned numPassedEvents_;

  public:
    explicit FilterStepPointPzVsTarget(const fhicl::ParameterSet& pset);
    virtual bool filter(art::Event& event) override;
    virtual void endJob() override;
  };

  //================================================================
  FilterStepPointPzVsTarget::FilterStepPointPzVsTarget(const fhicl::ParameterSet& pset)
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
        <<"FilterStepPointPzVsTarget: either pdgToDrop or pdgToKeep must be specified.\n";
    }
    if(!drop.empty() && ! keep.empty()) {
      throw cet::exception("BADCONFIG")
        <<"FilterStepPointPzVsTarget: either pdgToDrop or pdgToKeep,"
        <<" but not both, can be used at a time.\n";
    }

    targetAngle_ = pset.get<double>("targetAngle", 14.);
    pzMin_ = pset.get<double>("pzMin", -9e99);
    pzMax_ = pset.get<double>("pzMax",  9e99);

    targetDirx_ = -1.0*sin(targetAngle_*CLHEP::deg);
    targetDiry_ = 0.0;
    targetDirz_ = -1.0*cos(targetAngle_*CLHEP::deg);
    //std::cout << "target direction (" << targetDirx_ << ", " << targetDiry_ << ", " << targetDirz_ << ")" << std::endl; 
  }

  //================================================================
  bool FilterStepPointPzVsTarget::filter(art::Event& event) {
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

        //const double x = hit.position().x();
        //const double y = hit.position().y();
        //const double z = hit.position().z();

        const double px = hit.momentum().x();
        const double py = hit.momentum().y();
        const double pz = hit.momentum().z();
        const double p = sqrt(px*px+py*py+pz*pz); 
        if (p == 0) continue;  

        // rotate around y by targetAngle
        double rotateAngle = targetAngle_*CLHEP::deg;
        //std::cout << "rotateAngle " << rotateAngle << std::endl;

        //double pxp = px*cos(rotateAngle) - pz*sin(rotateAngle);
        //double pyp = py;
        double pzp = px*sin(rotateAngle) + pz*cos(rotateAngle);
        //std::cout << "p  (" << px  << ", " << py  << ", " << pz  << ")" << std::endl;
        //std::cout << "pp (" << pxp << ", " << pyp << ", " << pzp << ")" << std::endl;

        if ( pzp < pzMin_ || pzp > pzMax_ ) {
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
  void FilterStepPointPzVsTarget::endJob() {
    mf::LogInfo("Summary")<<"FilterStepPointPzVsTarget: passed "
                          <<numOutputHits_ <<" / "<<numInputHits_
                          <<" StepPointMCs, "
                          <<numPassedEvents_<<" / "<<numInputEvents_
                          <<" events\n";
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::FilterStepPointPzVsTarget);
