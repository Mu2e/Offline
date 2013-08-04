// Given a set of "interesting" hits, write out all SimParticles
// making the hits and all their parents.  Also write out all hits
// made by the selected particles in a specified set of "extra"
// collections, preserving the original collection instance names.
//
//   E.g:
//         "interesting" = TS3Vacuum, DS2Vacuum, DS3Vacuum, CRV
//         "extra" = virtualdetector
//
//   Outputs:  TS3Vacuum, DS2Vacuum, DS3Vacuum, CRV, virtualdetector, SimParticleCollection
//
// $Id: FilterStepPointMCs_module.cc,v 1.1 2013/08/04 14:28:27 gandr Exp $
// $Author: gandr $
// $Date: 2013/08/04 14:28:27 $
//
// Andrei Gaponenko, 2013

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <iterator>

// art includes.
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SelectorBase.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"

#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"

#include "MCDataProducts/inc/GenParticleSPMHistory.hh"
#include "MCDataProducts/inc/GenSimParticleLink.hh"

#include "Mu2eUtilities/inc/compressSimParticleCollection.hh"
#include "Mu2eUtilities/inc/SimParticleParentGetter.hh"

namespace mu2e {

  namespace {

    // Use art::Ptr instead of bare ptr to get the desired sorting
    typedef std::set<art::Ptr<SimParticle> > SPSet;
    // Keep particles from different collections separately
    struct SPSuperSet : public std::map<art::ProductID, SPSet> {
      bool isKept(const art::Ptr<SimParticle>& part) const {
        bool res = false;
        const auto si = this->find(part.id());
        if(si != this->end()) {
          const auto sp = si->second.find(part);
          res = (sp != si->second.end());
        }
        return res;
      }
    };

    // Adapter for compressSimParticleCollection()
    class ParticleSelector {
    public:
      ParticleSelector(const SPSet& m) {
        for(SPSet::const_iterator i = m.begin(); i!=m.end(); ++i) {
          m_keys.insert((*i)->id());
        }
      }

      bool operator[]( cet::map_vector_key key ) const {
        return m_keys.find(key) != m_keys.end();
      }

    private:
      std::set<cet::map_vector_key> m_keys;
    };

    //----------------
    class ProductIDSelector : public art::SelectorBase {
      const art::Event *evt_;
      art::ProductID pid_;
      virtual bool doMatch(const art::BranchDescription& p) const override {
        return evt_->branchIDToProductID(p.branchID()) == pid_;
      }
    public:
      ProductIDSelector(const art::Event& evt, const art::ProductID& pid) : evt_(&evt), pid_(pid) {}
      virtual ProductIDSelector* clone() const { return new ProductIDSelector(*this); }
    };

  } // anonymous namespace

  //================================================================
  class FilterStepPointMCs : public art::EDFilter {

    typedef std::vector<art::InputTag> InputTags;
    InputTags mainHitInputs_; // define a set of SimParticles
    InputTags extraHitInputs_;  // other StepPointMCCollections to keep

    InputTags gpStepLinkInputs_;
    InputTags gpSimLinkInputs_;

    InputTags particleInputs_;

    // Output instance names.
    typedef std::set<std::string> OutputNames;
    OutputNames mainOutputNames_;
    OutputNames extraOutputNames_;
    OutputNames simPartOutNames;

    // statistics counters
    unsigned numInputEvents_;
    unsigned numPassedEvents_;

    unsigned numMainHits_;
    unsigned numInputExtraHits_;
    unsigned numPassedExtraHits_;

    unsigned numInputParticles_;
    unsigned numPassedParticles_;

  public:
    explicit FilterStepPointMCs(const fhicl::ParameterSet& pset);
    virtual bool filter(art::Event& event) override;
    virtual void endJob() override;
  };

  //================================================================
  FilterStepPointMCs::FilterStepPointMCs(const fhicl::ParameterSet& pset)
    : numInputEvents_(), numPassedEvents_()
    , numMainHits_(), numInputExtraHits_(), numPassedExtraHits_()
    , numInputParticles_(), numPassedParticles_()
  {
    typedef std::vector<std::string> VS;

    const VS mainStrings(pset.get<VS>("mainHitInputs"));
    for(const auto& i : mainStrings) {
      mainHitInputs_.emplace_back(i);
      // Coalesce same instance names from multiple input modules/processes.
      mainOutputNames_.insert(mainHitInputs_.back().instance());
    }
    for(const auto& i : mainOutputNames_) {
      produces<StepPointMCCollection>(i);
    }

    const VS extraStrings(pset.get<VS>("extraHitInputs", VS()));
    for(const auto& i : extraStrings) {
      extraHitInputs_.emplace_back(i);
      // Coalesce same instance names from multiple input modules/processes.
      extraOutputNames_.insert(extraHitInputs_.back().instance());
    }
    for(const auto& i : extraOutputNames_) {
      produces<StepPointMCCollection>(i);
    }

    // We can't merge different SimParticle collections (unlike the hits)
    // Need to have a separate output collection for every input one.
    // How many SimParticleCollection-s do we produce?
    const unsigned numSimPartOuts(pset.get<unsigned>("numSimParticleCollections"));
    for(unsigned i = 0; i < numSimPartOuts; ++i) {
      std::ostringstream os;
      os<<"s"<<i;
      simPartOutNames.insert(os.str());
      produces<SimParticleCollection>(os.str());
    }

    // FIXME:
    // Need to update hits pointers in GenParticleSPMHistory objects.
    const VS gpStepLinks(pset.get<VS>("gpStepLinks", VS()));
    for(const auto& i : gpStepLinks) {
      gpStepLinkInputs_.emplace_back(i);
    }
    if(!gpStepLinkInputs_.empty()) {
      produces<GenParticleSPMHistory>();
    }

    // Need to update hits pointers in GenParticleSPMHistory objects.
    const VS gpSimLinks(pset.get<VS>("gpSimLinks", VS()));
    for(const auto& i : gpSimLinks) {
      gpSimLinkInputs_.emplace_back(i);
    }
    if(!gpSimLinkInputs_.empty()) {
      produces<GenSimParticleLink>();
    }
  }

  //================================================================
  // Return true is any hits are passed
  bool FilterStepPointMCs::filter(art::Event& event) {
    bool passed = false;
    typedef std::map<std::string, std::unique_ptr<StepPointMCCollection> > OutMap;

    //----------------------------------------------------------------
    // Build list of all the SimParticles we want to preserve:
    SPSuperSet toBeKept;

    // These are all particles with hits in the "main" collections
    for(const auto& i : mainHitInputs_) {
      auto ih = event.getValidHandle<StepPointMCCollection>(i);
      for(StepPointMCCollection::const_iterator i=ih->begin(); i!=ih->end(); ++i) {
        const art::Ptr<SimParticle>& particle(i->simParticle());
        toBeKept[particle.id()].insert(i->simParticle());
      }
    }

    // and their parents
    SimParticleParentGetter pg(event);
    for(const auto& iset : toBeKept) {
      for(const auto& ipart : iset.second) {
        art::Ptr<SimParticle> next = pg.parent(ipart);
        while(next) {
          // Insertion into a set does not invalidate the iterator
          toBeKept[next.id()].insert(next);
          next = pg.parent(next);
        }
      }
    }

    //----------------------------------------------------------------
    // Prepare output SimParticle collections

    // old => new collection
    typedef std::map<art::ProductID, art::ProductID> PIDMap;
    PIDMap partCollMap;

    if(toBeKept.size() > simPartOutNames.size()) {
      throw cet::exception("BADCONFIG")
        <<"FilterStepPointMCs: configured numSimParticleCollections = "<<simPartOutNames.size()
        <<" but used "<<toBeKept.size()<<" collections in the event\n";
    }

    // We map input to output collections "randomly".
    // is there a reason to worry about the output names?
    // We always output the specified number of collections
    // event if they are empty.

    SPSuperSet::const_iterator iss = toBeKept.begin();
    for(const auto& outInstance : simPartOutNames) {

      std::unique_ptr<SimParticleCollection> outparts(new SimParticleCollection());
      art::ProductID newParticlesPID(getProductID<SimParticleCollection>(event, outInstance));
      const art::EDProductGetter *newParticlesGetter(event.productGetter(newParticlesPID));

      // Is there anything to copy into this output?
      if(iss != toBeKept.end()) {
        ProductIDSelector csel(event, iss->first);
        art::Handle<SimParticleCollection> inputParticles;
        event.get(csel, inputParticles);

        // old => new particle collection ID
        partCollMap[iss->first] = newParticlesPID;

        ParticleSelector psel(iss->second);
        compressSimParticleCollection(newParticlesPID,
                                      newParticlesGetter,
                                      *inputParticles,
                                      psel,
                                      *outparts);

        numInputParticles_ += inputParticles->size();
        numPassedParticles_ += outparts->size();

        ++iss;
      }

      passed = passed || !outparts->empty();
      event.put(std::move(outparts), outInstance);
    }

    //----------------------------------------------------------------
    OutMap outMain;
    for(const auto& i : mainOutputNames_) {
      std::unique_ptr<StepPointMCCollection> p(new StepPointMCCollection());
      outMain.insert(std::move(std::make_pair(i, std::move(p))));
    }

    for(const auto& inTag : mainHitInputs_) {
      auto ih = event.getValidHandle<StepPointMCCollection>(inTag);
      numMainHits_ += ih->size();

      for(StepPointMCCollection::const_iterator i=ih->begin(); i!=ih->end(); ++i) {

        art::Ptr<SimParticle> oldPtr(i->simParticle());
        assert(toBeKept.isKept(oldPtr));

        art::ProductID newParticlesPID = partCollMap[oldPtr.id()];
        const art::EDProductGetter *newParticlesGetter(event.productGetter(newParticlesPID));

        StepPointMCCollection& output = *outMain[inTag.instance()];

        output.emplace_back(*i);
        output.back().simParticle() = art::Ptr<SimParticle>(newParticlesPID, oldPtr->id().asUint(), newParticlesGetter);
      }
    }

    for(const auto& i : mainOutputNames_) {
      event.put(std::move(outMain[i]), i);
    }

    //----------------------------------------------------------------
    OutMap outExtra;
    for(const auto& i : extraOutputNames_) {
      std::unique_ptr<StepPointMCCollection> p(new StepPointMCCollection());
      outExtra.insert(std::move(std::make_pair(i, std::move(p))));
    }

    // Need to filter extra hits. Keep only those with a preserved SimParticle.
    for(const auto& inTag : extraHitInputs_) {
      auto ih = event.getValidHandle<StepPointMCCollection>(inTag);
      numInputExtraHits_ += ih->size();

      for(StepPointMCCollection::const_iterator i=ih->begin(); i!=ih->end(); ++i) {

        art::Ptr<SimParticle> oldPtr(i->simParticle());
        if(toBeKept.isKept(oldPtr)) {

          art::ProductID newParticlesPID = partCollMap[oldPtr.id()];
          const art::EDProductGetter *newParticlesGetter(event.productGetter(newParticlesPID));

          StepPointMCCollection& output = *outExtra[inTag.instance()];

          output.emplace_back(*i);
          output.back().simParticle() = art::Ptr<SimParticle>(newParticlesPID, oldPtr->id().asUint(), newParticlesGetter);
        }
      }
    }

    for(const auto& i : extraOutputNames_) {
      numPassedExtraHits_ += outExtra[i]->size();
      event.put(std::move(outExtra[i]), i);
    }

    //----------------------------------------------------------------
    ++numInputEvents_;
    if(passed) ++numPassedEvents_;

    return passed;
  }

  //================================================================
  void FilterStepPointMCs::endJob() {
    mf::LogInfo("Summary")
      << "FilterStepPointMCs_module stats: passed "
      << numPassedEvents_ <<" / "<<numInputEvents_<<" events, "
      << numMainHits_ <<" main hits, "
      << numPassedExtraHits_ <<" / "<<numInputExtraHits_<<" extra hits, "
      << numPassedParticles_ <<" / "<<numInputParticles_<<"+ particles, "
      << "\n";
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::FilterStepPointMCs);
