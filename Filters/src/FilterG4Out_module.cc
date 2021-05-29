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
// An optional vetoDaughters input specifies a set of particles whose
// daughters (and  further descendants) should be excluded from output,
// both the particle collection and the hits.    The intent is to
// veto daughters of particles stopped in the stopping target, because
// they will be simulated with a foil generator in other jobs.
//
// vetoParticles works similar, but also vetoes the particle listed, not
// just the daughters.
//
// The other use mode is to specify a SimParticlePtrCollection of stuff to keep.
// Intended to write out framework files of stopped muons.
//
//
// Andrei Gaponenko, 2013

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <iterator>

// art includes.
#include "fhiclcpp/types/OptionalSequence.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SelectorBase.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/MCTrajectory.hh"
#include "MCDataProducts/inc/MCTrajectoryCollection.hh"
#include "Mu2eUtilities/inc/compressSimParticleCollection.hh"
#include "Mu2eUtilities/inc/SimParticleParentGetter.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/SimParticlePtrCollection.hh"
#include "MCDataProducts/inc/SimParticleRemapping.hh"
#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"

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
      art::ProductID pid_;
      virtual bool doMatch(const art::BranchDescription& p) const override {
        return p.productID() == pid_;
      }
    public:
      ProductIDSelector( const art::ProductID& pid) : pid_(pid) {}
      std::string doPrint(std::string const& indent) const override {
        std::ostringstream os;
        os << pid_;
        return indent + "art::ProductID: " + os.str();
      }
    };

    //----------------
    void addDaughterTree(SPSet *res, const art::Ptr<SimParticle>& p) {
      for(const auto& d : p->daughters()) {
        res->insert(d);
        addDaughterTree(res, d);
      }
    }

  } // anonymous namespace

  //================================================================
  class FilterG4Out : public art::EDFilter {
  public:

    struct IOMapEntry {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Atom<art::InputTag> in {
        Name("in"),
          Comment("The input collection")
          };

      fhicl::Atom<std::string> out {
        Name("out"),
          Comment("Map the input to this instance name on the output")
          };
    };

    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Sequence<art::InputTag> mainHitInputs {
        Name("mainHitInputs"),
          Comment("Particles and StepPointMCs mentioned in thise collections will be preserved.")
          };

      fhicl::Sequence<art::InputTag> mainSPPtrInputs {
        Name("mainSPPtrInputs"),
          Comment("A list of SimParticlePtrCollections. The pointed to SimParticles will be preserved."),
          std::vector<art::InputTag>()
          };

      fhicl::Sequence<art::InputTag> extraHitInputs {
        Name("extraHitInputs"),
          Comment("If a hit in one of these collection is produced by a preserved SimParticle, the hit will be kept."),
          std::vector<art::InputTag>()
          };

      fhicl::Sequence<art::InputTag> mcTrajectoryInputs {
        Name("mcTrajectoryInputs"),
          Comment("Trajectories to preserve if the corresponding SimParticle is preserved."),
          std::vector<art::InputTag>()
          };

      fhicl::Sequence<art::InputTag> vetoDaughters {
        Name("vetoDaughters"),
          Comment("SimParticlePtrCollections; the daughters and further descendants of the pointed-to particles will be dropped."),
          std::vector<art::InputTag>()
          };

      fhicl::Sequence<art::InputTag> vetoParticles {
        Name("vetoParticles"),
          Comment("SimParticlePtrCollections; the particles and all their descendants will be dropped.."),
          std::vector<art::InputTag>()
          };

      fhicl::OptionalSequence<fhicl::Table<IOMapEntry> > simParticleIOMap {
        Name("simParticleIOMap"),
          Comment("Specify an output instance name for each SimParticleCollection referenced by the inputs.\n"
                  "This should only be used when inputs depend on multiple SimParticle collections (e.g. in mixed events)."
                  )
          };



      // It looks like the GenParticle compression code can write out multiple copies
      // of the same particle.  The setting is not used by any fcl files in the
      // current master branch of Offline.  Comment it out untile the issue is
      // verified and resolved.

      //FIXME:fhicl::Atom<bool> compressGenParticles {
      //FIXME:  Name("compressGenParticles"),
      //FIXME:    Comment("By default GenParticles are not touched, and the output SimParticles will point\n"
      //FIXME:            "to the original GenParticleCollection. If GenParticle compression is requested,\n"
      //FIXME:            "a new GenParticleCollection will be produced, and output SimParticles will point to it."
      //FIXME:            ),
      //FIXME:    false
      //FIXME:    };



    };

    using Parameters = art::EDFilter::Table<Config>;
    explicit FilterG4Out(const Parameters& conf);

    virtual bool filter(art::Event& event) override;
    virtual void endJob() override;

  private:
    typedef std::vector<art::InputTag> InputTags;
    InputTags mainHitInputs_; // These hits, and SimParticles referenced here, will be kept
    InputTags mainSPPtrInputs_; // SimParticles referenced here, will be kept

    InputTags extraHitInputs_;  // other StepPointMCCollections to keep

    InputTags trajectoryInputs_; // Optionally, filter and keep MCTrajectoryCollection

    InputTags vetoDaughtersInputs_;
    InputTags vetoParticlesInputs_;

    bool compressGenParticles_;

    // Output instance names.
    typedef std::set<std::string> OutputNames;
    OutputNames mainOutputNames_;
    OutputNames extraOutputNames_;

    typedef std::vector<IOMapEntry> VPM;
    VPM simParticleIOVec_;
    OutputNames simPartOutNames;

    // statistics counters
    unsigned numInputEvents_;
    unsigned numPassedEvents_;

    unsigned numMainHits_;
    unsigned numInputExtraHits_;
    unsigned numPassedExtraHits_;

    unsigned numInputParticles_;
    unsigned numPassedParticles_;

    unsigned numVetoedParticles_;
    unsigned numVetoedHits_;
  };

  //================================================================
  FilterG4Out::FilterG4Out(const Parameters& conf)
    : art::EDFilter{conf}
    , compressGenParticles_(false /*FIXME: compressGenParticles_(conf().compressGenParticles() */)
    , numInputEvents_(), numPassedEvents_()
    , numMainHits_(), numInputExtraHits_(), numPassedExtraHits_()
    , numInputParticles_(), numPassedParticles_()
    , numVetoedParticles_(), numVetoedHits_()
  {
    for(const auto& i : conf().mainHitInputs()) {
      mainHitInputs_.emplace_back(i);
      // Coalesce same instance names from multiple input modules/processes.
      mainOutputNames_.insert(mainHitInputs_.back().instance());
    }
    for(const auto& i : mainOutputNames_) {
      produces<StepPointMCCollection>(i);
    }

    for(const auto& i : conf().mainSPPtrInputs()) {
      mainSPPtrInputs_.emplace_back(i);
    }

    for(const auto& i : conf().extraHitInputs()) {
      extraHitInputs_.emplace_back(i);
      // Coalesce same instance names from multiple input modules/processes.
      extraOutputNames_.insert(extraHitInputs_.back().instance());
    }
    for(const auto& i : extraOutputNames_) {
      produces<StepPointMCCollection>(i);
    }

    for(const auto& i : conf().mcTrajectoryInputs()) {
      trajectoryInputs_.emplace_back(i);
    }
    if(!trajectoryInputs_.empty()) {
      produces<MCTrajectoryCollection>();
    }

    for(const auto& i : conf().vetoDaughters()) {
      vetoDaughtersInputs_.emplace_back(i);
    }

    for(const auto& i : conf().vetoParticles()) {
      vetoParticlesInputs_.emplace_back(i);
    }

    // We can't merge different SimParticle collections (unlike the hits)
    // Need to have a separate output collection for every input one.
    conf().simParticleIOMap(simParticleIOVec_);

    if(!simParticleIOVec_.empty()) {

      // FIXME: the check for duplicate InputTags
      // needs https://cdcvs.fnal.gov/redmine/issues/23320
      //
      //redmine#23320:std::set<art::InputTag> inset;

      for(const auto & i : simParticleIOVec_) {
        //redmine#23320:if(!inset.insert(i.in()).second) {
        //redmine#23320:  throw cet::exception("BADCONFIG")
        //redmine#23320:    <<"FilterG4Out: duplicate input tag in simParticleIOMap: in: "<<i.in()<<"\n";
        //redmine#23320:}

        if(!simPartOutNames.insert(i.out()).second) {
          throw cet::exception("BADCONFIG")
            <<"FilterG4Out: duplicate out name in simParticleIOMap: "<<" out: "<<i.out()<<"\n";
        }

        produces<SimParticleCollection>(i.out());
      }
    }
    else {
      const std::string defaultInstance;
      simPartOutNames.insert(defaultInstance);
      produces<SimParticleCollection>(defaultInstance);
    }

    produces<SimParticleRemapping>();

    if(compressGenParticles_) {
      produces<mu2e::GenParticleCollection>();
    }
  }

  //================================================================
  // Return true is any hits are passed
  bool FilterG4Out::filter(art::Event& event) {
    bool passed = false;
    typedef std::map<std::string, std::unique_ptr<StepPointMCCollection> > OutMap;

    //----------------------------------------------------------------
    // Build a full list of the vetoed particles

    SPSet vetoedParticles;

    for(const auto& tag : vetoDaughtersInputs_) {
      auto ih = event.getValidHandle<SimParticlePtrCollection>(tag);
      for(const auto& p : *ih) {
        addDaughterTree(&vetoedParticles, p);
      }
    }

    for(const auto& tag : vetoParticlesInputs_) {
      auto ih = event.getValidHandle<SimParticlePtrCollection>(tag);
      for(const auto& p : *ih) {
        vetoedParticles.insert(p);
        addDaughterTree(&vetoedParticles, p);
      }
    }

    numVetoedParticles_ += vetoedParticles.size();

    //----------------------------------------------------------------
    // Build list of all the SimParticles we want to preserve:
    SPSuperSet toBeKept;

    // Non-vetoed particles with hits in the "main" collections
    for(const auto& tag : mainHitInputs_) {
      auto ih = event.getValidHandle<StepPointMCCollection>(tag);
      for(const auto& i : *ih) {
        const art::Ptr<SimParticle>& particle(i.simParticle());
        if(vetoedParticles.find(particle) == vetoedParticles.end()) {
          toBeKept[particle.id()].insert(particle);
        }
      }
    }

    // Non-vetoed particles explicitly requiested
    for(const auto& tag : mainSPPtrInputs_) {
      auto ih = event.getValidHandle<SimParticlePtrCollection>(tag);
      for(const auto& particle : *ih) {
        if(vetoedParticles.find(particle) == vetoedParticles.end()) {
          toBeKept[particle.id()].insert(particle);
        }
      }
    }

    // and their parents (can not be vetoed since the daughter is in)
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

    // input PID => output instance
    typedef std::map<art::ProductID, std::string> SimParticleInstanceMap;
    SimParticleInstanceMap spim;

    if(!simParticleIOVec_.empty()) {
      for(const auto& i : simParticleIOVec_) {
        auto ih = event.getValidHandle<SimParticleCollection>(i.in());
        if(!spim.insert(std::make_pair(ih.id(), i.out())).second) {
          throw cet::exception("BADCONFIG")
            <<"FilterG4Out: Different entries of simParticleIOMap resolved to the same product!\n"
            <<"The current one: { in : "<<i.in()<<" out: "<<i.out()<<" }\n";
        }
      }
    }
    else { // abritrary assignment out output instance names
      if(toBeKept.size() > simPartOutNames.size()) {
        throw cet::exception("BADCONFIG")
          <<"FilterG4Out: configured numSimParticleCollections = "<<simPartOutNames.size()
          <<" but used "<<toBeKept.size()<<" collections in the event\n";
      }
      auto instance = simPartOutNames.begin();
      for(const auto& i : toBeKept) {
        spim.insert(std::make_pair(i.first, *instance++));
      }

      // We always output the specified number of collections
      // event if they are empty.
      while(instance != simPartOutNames.end()) {
        std::unique_ptr<SimParticleCollection> outparts(new SimParticleCollection());
        event.put(std::move(outparts), *instance++);
      }
    }

    // old => new collection
    typedef std::map<art::ProductID, art::ProductID> PIDMap;
    PIDMap partCollMap;

    std::unique_ptr<GenParticleCollection> genParts(new GenParticleCollection());
    art::ProductID newGenPID(compressGenParticles_ ? event.getProductID<GenParticleCollection>() : art::ProductID());
    const art::EDProductGetter *newGenGetter(compressGenParticles_ ? event.productGetter(newGenPID) : nullptr);

    for(const auto& iopair : spim) {
      SPSuperSet::const_iterator iss = toBeKept.find(iopair.first);
      const auto& outInstance = iopair.second;

      std::unique_ptr<SimParticleCollection> outparts(new SimParticleCollection());
      art::ProductID newParticlesPID(event.getProductID<SimParticleCollection>(outInstance));
      const art::EDProductGetter *newParticlesGetter(event.productGetter(newParticlesPID));

      // Is there anything to copy into this output?
      if(iss != toBeKept.end()) {
        ProductIDSelector csel(iss->first);
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


        if(compressGenParticles_) {
          for(auto& i : *outparts) {
            mu2e::SimParticle& newsim = i.second;
            if(!newsim.genParticle().isNull()) { // will crash if not resolvable
              // Copy GenParticle to the new collection
              genParts->emplace_back(*newsim.genParticle());
              newsim.genParticle() = art::Ptr<GenParticle>(newGenPID, genParts->size()-1, newGenGetter);
            }
          }
        }

        numInputParticles_ += inputParticles->size();
        numPassedParticles_ += outparts->size();
      }

      passed = passed || !outparts->empty();
      event.put(std::move(outparts), outInstance);
    }

    if(compressGenParticles_) {
      event.put(std::move(genParts));
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
        if(toBeKept.isKept(oldPtr)) {

          art::ProductID newParticlesPID = partCollMap[oldPtr.id()];
          const art::EDProductGetter *newParticlesGetter(event.productGetter(newParticlesPID));

          StepPointMCCollection& output = *outMain[inTag.instance()];

          art::Ptr<SimParticle> newParticle(newParticlesPID, oldPtr->id().asUint(), newParticlesGetter);
          output.emplace_back(*i);
          output.back().simParticle() = newParticle;
        }
        else {
          ++numVetoedHits_;
        }
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

    for(const auto& inTag : extraHitInputs_) {
      auto ih = event.getValidHandle<StepPointMCCollection>(inTag);
      numInputExtraHits_ += ih->size();

      for(StepPointMCCollection::const_iterator i=ih->begin(); i!=ih->end(); ++i) {

        art::Ptr<SimParticle> oldPtr(i->simParticle());
        if(toBeKept.isKept(oldPtr)) {
          ++numPassedExtraHits_;

          art::ProductID newParticlesPID = partCollMap[oldPtr.id()];
          const art::EDProductGetter *newParticlesGetter(event.productGetter(newParticlesPID));

          StepPointMCCollection& output = *outExtra[inTag.instance()];

          art::Ptr<SimParticle> newParticle(newParticlesPID, oldPtr->id().asUint(), newParticlesGetter);
          output.emplace_back(*i);
          output.back().simParticle() = newParticle;
        }
        else {
          if(vetoedParticles.find(oldPtr) != vetoedParticles.end()) {
            ++numVetoedHits_;
          }
        }
      }
    }

    for(const auto& i : extraOutputNames_) {
      numPassedExtraHits_ += outExtra[i]->size();
      event.put(std::move(outExtra[i]), i);
    }

    //----------------------------------------------------------------
    if(!trajectoryInputs_.empty()) {

      std::unique_ptr<MCTrajectoryCollection> outTrajectory(new MCTrajectoryCollection());

      for(const auto& in : trajectoryInputs_) {

        auto ih = event.getValidHandle<MCTrajectoryCollection>(in);

        for(MCTrajectoryCollection::const_iterator i=ih->begin(); i!=ih->end(); ++i) {

          art::Ptr<SimParticle> oldPtr(i->first);

          if(toBeKept.isKept(oldPtr)) {

            art::ProductID newParticlesPID = partCollMap[oldPtr.id()];
            const art::EDProductGetter *newParticlesGetter(event.productGetter(newParticlesPID));
            art::Ptr<SimParticle> newParticle(newParticlesPID, oldPtr->id().asUint(), newParticlesGetter);

            (*outTrajectory)[newParticle] = i->second;
            (*outTrajectory)[newParticle].sim() = newParticle;
          }
        }
      }

      event.put(std::move(outTrajectory));
    }

    //----------------------------------------------------------------
    std::unique_ptr<SimParticleRemapping> remap(new SimParticleRemapping());

    for(const auto& x : toBeKept) {
      for(const auto& oldPtr : x.second) {
        art::ProductID newParticlesPID = partCollMap[oldPtr.id()];
        const art::EDProductGetter *newParticlesGetter(event.productGetter(newParticlesPID));
        (*remap)[oldPtr] = art::Ptr<SimParticle>(newParticlesPID, oldPtr->id().asUint(), newParticlesGetter);
      }
    }

    event.put(std::move(remap));

    //----------------------------------------------------------------
    ++numInputEvents_;
    if(passed) ++numPassedEvents_;

    return passed;
  }

  //================================================================
  void FilterG4Out::endJob() {
    mf::LogInfo("Summary")
      << "FilterG4Out_module stats: passed "
      << numPassedEvents_ <<" / "<<numInputEvents_<<" events, "
      << numMainHits_ <<" main hits, "
      << numPassedExtraHits_ <<" / "<<numInputExtraHits_<<" extra hits, "
      << numPassedParticles_ <<" / "<<numInputParticles_<<"+ particles, "
      << "vetoed " << numVetoedParticles_ << " particles and "
      << numVetoedHits_ <<" hits."
      << "\n";
  }

  //================================================================
  std::ostream& operator<<(std::ostream& out, const FilterG4Out::IOMapEntry& a) {
    return out<<"IOMapEntry{ in: "<<a.in()<<", out: "<<a.out()<<" }";
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::FilterG4Out);
