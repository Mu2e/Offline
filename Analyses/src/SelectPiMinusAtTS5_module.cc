// Andrei Gaponenko, 2011

#include <string>
#include <vector>
#include <algorithm>

// Mu2e includes.
#include "MCDataProducts/inc/StatusG4.hh"
#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "DataProducts/inc/PDGCode.hh"

#include "Mu2eUtilities/inc/compressSimParticleCollection.hh"

// art includes.
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"

#define AGDEBUG(stuff) std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<": "<<stuff<<std::endl;
//#define AGDEBUG(stuff)

namespace mu2e {

  // Use art::Ptr instead of bare ptr to get the desired sorting
  typedef std::set<art::Ptr<SimParticle> > TrackSet;

  // Adapter for compressSimParticleCollection()
  class ParticleSelector {
  public:
    ParticleSelector(const TrackSet& m) {
      for(TrackSet::const_iterator i = m.begin(); i!=m.end(); ++i) {
        m_keys.insert((*i)->id());
      }
    }

    bool operator[]( cet::map_vector_key key ) const {
      return m_keys.find(key) != m_keys.end();
    }

  private:
    std::set<cet::map_vector_key> m_keys;
  };

  //================================================================
  class SelectPiMinusAtTS5 : public art::EDFilter {
    std::string _outInstanceName;
    std::string _inModuleLabel;
    std::string _inInstanceName;

    typedef StepPointMC::VolumeId_type VolumeId;

    // The anticipated use case is to select hits for one or a few
    // VDs.  Therefore we use a vector: the log-vs-linear efficiency
    // gain from using a set is probably negated by the overheads
    // of dealing with a more complicated data structure that is not
    // local in memory.
    std::vector<VolumeId> _vids;  // preserve hits in just those volumes

    bool _positionCut;
    std::vector<double> _positionCenter;
    std::vector<double> _positionHalfLength;

    bool _storeParents;

    bool _storeExtraHits; // store all VD hits for saved particles, not just in _vids.

  public:
    explicit SelectPiMinusAtTS5(const fhicl::ParameterSet& pset);
    virtual bool filter(art::Event& event);
    virtual bool beginRun(art::Run& run) { std::cout<<"AG: beginRun() called"<<std::endl; return true; }
    virtual bool endRun(art::Run& run) { std::cout<<"AG: endRun() called"<<std::endl; return true; }
  };

  //================================================================
  SelectPiMinusAtTS5::SelectPiMinusAtTS5(const fhicl::ParameterSet& pset)
    : EDFilter{pset}
    , _inModuleLabel(pset.get<std::string>("inputModuleLabel"))
    , _inInstanceName(pset.get<std::string>("inputInstanceName"))
    , _vids(pset.get<std::vector<VolumeId> >("acceptedVids"))
    , _positionCut(pset.get<bool>("positionCut", false))
    , _positionCenter(std::vector<double>(3, 0.0))
    , _positionHalfLength(std::vector<double>(3, 0.0))
    , _storeParents(pset.get<bool>("storeParents"))
      // default to false for compatibility with existing .fcl files.
    , _storeExtraHits(pset.get<bool>("storeExtraHits", false))
  {
    if (_positionCut) {
      _positionCenter = pset.get<std::vector<double> >("positionCenter");
      _positionHalfLength = pset.get<std::vector<double> >("positionHalfLength");
    }

    std::cout<<"SelectPiMinusAtTS5(): storeParents = "<<_storeParents<<std::endl;
    std::cout<<"SelectPiMinusAtTS5(): storeExtraVDs = "<<_storeExtraHits<<std::endl;

    if(_storeExtraHits && !_storeParents) {
      throw cet::exception("BADCONFIG")
        <<"SelectPiMinusAtTS5: attempting to storeExtraHits without storeParents probably does not make sense.";
    }

    produces<StepPointMCCollection>();
    produces<SimParticleCollection>();
    if(_storeExtraHits) {
      produces<StepPointMCCollection>("extraHits");
    }
  }

  //================================================================
  bool SelectPiMinusAtTS5::filter(art::Event& event) {
    AGDEBUG("SelectPiMinusAtTS5 begin event "<<event.id());

    // Use art::Ptr instead of bare ptr to get the desired sorting
    TrackSet particlesWithHits;

    art::Handle<StepPointMCCollection> ih;
    event.getByLabel(_inModuleLabel, _inInstanceName, ih);
    const StepPointMCCollection& inhits(*ih);

    std::unique_ptr<StepPointMCCollection> outhits(new StepPointMCCollection());
    std::unique_ptr<StepPointMCCollection> extrahits(new StepPointMCCollection());


    for(StepPointMCCollection::const_iterator i=inhits.begin(); i!=inhits.end(); ++i) {

      if(std::find(_vids.begin(), _vids.end(), i->volumeId()) != _vids.end()) {
        std::cout << "particle type " << ((*i).simParticle())->pdgId() << std::endl;
        if ( ((*i).simParticle())->pdgId() != PDGCode::pi_minus ) continue;

        if ( _positionCut &&
             ( fabs(i->position().x() - _positionCenter[0]) > _positionHalfLength[0] ||
               fabs(i->position().y() - _positionCenter[1]) > _positionHalfLength[1] ||
               fabs(i->position().z() - _positionCenter[2]) > _positionHalfLength[2] ) )
          continue;

        AGDEBUG("here");
        const art::Ptr<SimParticle>& particle = outhits->back().simParticle();
        AGDEBUG("here");

        if (particle->pdgId() != PDGCode::pi_minus) continue;

          AGDEBUG("here: found a piminus!! "<<particle<<" (internal id = "<< particle->id()<< " pdgId = " << particle->pdgId() <<")"<<" for hit "<<*i);


        if(!particle.get()) {
          throw cet::exception("MISSINGINFO")
            <<"NULL particle pointer for StepPointMC = "<<*i
            <<" in event "<<event.id()
            ;
        }
        else {
          AGDEBUG("here: particle = "<<particle<<" (internal id = "<< particle->id()<< " pdgId = " << particle->pdgId() <<")"<<" for hit "<<*i);
          particlesWithHits.insert(particle);
        }
      }
    }
    AGDEBUG("here");

    if(_storeParents) {
      // For each particle hitting our VD also save
      // all the parents in the chain up to the primary.
      for(TrackSet::const_iterator i=particlesWithHits.begin(); i!=particlesWithHits.end(); ++i) {
        art::Ptr<SimParticle> current = *i;
        while(current->hasParent()) {
          current = current->parent();
          // Insertion into the set does not invalidate the iterator (i)
          particlesWithHits.insert(current);
        }
      }
    }

    // The case when parents are stored is supported by compressSimParticleCollection()
    // otherwise we need to prepare the output SimParticleCollection by hand
    std::unique_ptr<SimParticleCollection> outparts(new SimParticleCollection());
    art::ProductID newProductId(event.getProductID<SimParticleCollection>());
    const art::EDProductGetter *newProductGetter(event.productGetter(newProductId));
    if(!_storeParents) {

      // Need to save SimParticles that produced the hits Particles
      // must come in order, and there may be no one-to-one
      // correspondence with saved hits, so this must be a separate
      // loop.
      //
      // Note that the art::Ptr comparison operator sorts the set in
      // the correct order for our use.
      for(TrackSet::const_iterator i=particlesWithHits.begin(); i!=particlesWithHits.end(); ++i) {

        AGDEBUG("here");
        cet::map_vector_key key = (*i)->id();
        AGDEBUG("here: key = "<<key);
        // Copy the original particle to the output

        // FIXME: why does push_back() screw up the given key?
        // outparts->push_back(std::make_pair(key, **i));

        (*outparts)[key] = **i;

        SimParticle& particle(outparts->getOrThrow(key));

        if (particle.pdgId() != PDGCode::pi_minus) continue;


        // Zero internal pointers: intermediate particles are not preserved to reduce data size
        AGDEBUG("here");
        particle.setDaughterPtrs(std::vector<art::Ptr<SimParticle> >());
        AGDEBUG("here");
        particle.parent() = art::Ptr<SimParticle>();

        AGDEBUG("after p/d reset: particle id = "<<particle.id());
      }
    } // if(_storeParents)
    else {
      art::Handle<SimParticleCollection> inparts;
      event.getByLabel(_inModuleLabel, "", inparts);

      ParticleSelector selector(particlesWithHits);
      compressSimParticleCollection(newProductId, newProductGetter, *inparts, selector, *outparts);
    } // else(_storeParents)

    //----------------------------------------------------------------
    // We have a set of saved particles.  Some of their hits are
    // already in the output collection, _storeExtraHits requrests to
    // add any remaining hits made by those particles.
    if(_storeExtraHits) {

      for(StepPointMCCollection::const_iterator i=inhits.begin(); i!=inhits.end(); ++i) {

        // Hits in _vids are already stored.  Just store the complement.
        if(std::find(_vids.begin(), _vids.end(), i->volumeId()) == _vids.end() ) {

          const art::Ptr<SimParticle>& particle = i->simParticle();

        if (particle->pdgId() != PDGCode::pi_minus) continue;

          if(!particle.get()) {
            throw cet::exception("MISSINGINFO")
              <<"storeExtraHits: NULL particle pointer for StepPointMC = "<<*i
              <<" in event "<<event.id()
              ;
          }
          else {
            AGDEBUG("here: particle = "<<particle<<" (internal id = "<< particle->id()<<")"<<" for hit "<<*i);
            if(outparts->has(particle->id())) {
              extrahits->push_back(*i);
            }
          }
        }
      }
    } // if(_storeExtraHits)

    //----------------------------------------------------------------

    AGDEBUG("here");
    event.put(std::move(outparts));

    // Update pointers in the hit collection
    AGDEBUG("here");
    for(StepPointMCCollection::iterator i=outhits->begin(); i!=outhits->end(); ++i) {
      AGDEBUG("here");
      art::Ptr<SimParticle> oldPtr(i->simParticle());
      AGDEBUG("here: settting id = "<< oldPtr->id()<<", newProductId = "<<newProductId <<" for hit "<<*i);
      i->simParticle() = art::Ptr<SimParticle>(newProductId, oldPtr->id().asUint(), newProductGetter);
    }
    AGDEBUG("here");

    if(_storeExtraHits) {
      // Update pointers in the extra hits collection
      AGDEBUG("here");
      for(StepPointMCCollection::iterator i=extrahits->begin(); i!=extrahits->end(); ++i) {
        AGDEBUG("here");
        art::Ptr<SimParticle> oldPtr(i->simParticle());
        AGDEBUG("here: settting id = "<< oldPtr->id()<<", newProductId = "<<newProductId <<" for hit "<<*i);
        i->simParticle() = art::Ptr<SimParticle>(newProductId, oldPtr->id().asUint(), newProductGetter);
      }
      AGDEBUG("here");

      event.put(std::move(extrahits), "extraHits");
    }

    //----------------------------------------------------------------

    bool nonEmpty = !outhits->empty();
    event.put(std::move(outhits));

    std::cout << " event nonempty? " << nonEmpty << std::endl;

    return nonEmpty;
  }

  //================================================================

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::SelectPiMinusAtTS5);
