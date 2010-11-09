#ifndef CONVELECUTILITIES_HH
#define CONVELECUTILITIES_HH

// C++ includes
#include <ostream>
#include <string>
#include <vector>
#include <map>
#include <cstdlib>

// Framework includes
#include "FWCore/Framework/interface/Event.h"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"

// Mu2e includes
#include "ToyDP/inc/StepPointMCCollection.hh"
#include "ToyDP/inc/ToyGenParticleCollection.hh"
#include "ToyDP/inc/SimParticleCollection.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"
#include "TrackerGeom/inc/StrawIndex.hh"



namespace mu2e {

  class ConvElecUtilities{
    
    typedef SimParticleCollection::key_type key_type;

  public:
    
    /**
     * Constructor.  Takes an event as an argument
     * plus g4module name (string) and the trackerStepPoints name (string)
     */
    ConvElecUtilities(const edm::Event & event,
                      std::string g4ModuleLabel,
                      std::string trackerStepPoints);

    ~ConvElecUtilities();
    
  public:

    int nConvElec() const;
    std::vector<size_t> convElecHitsIdx();
    std::vector<StrawIndex> convElecStrawIdx();
    size_t hasStepPointMC() const;
    StepPointMC const& firstHit();
    StrawIndex earliestStrawIndex() const;
    const SimParticle& simConvElec();
    const ToyGenParticle& genConvElec();

  private:
    
    void checkConvElec(const edm::Event & event);
    void lookAtHits(const edm::Event & event);
    int _nconv;
    edm::Handle<SimParticleCollection> _simParticles;
    edm::Handle<ToyGenParticleCollection> _genParticles;
    edm::Handle<StepPointMCCollection> hits;
    std::vector<size_t> _convElecHits; 
    std::vector<StrawIndex> _convElecStrawIdx;
    std::string _g4ModuleLabel, _trackerStepPoints;
    key_type _convTrackId;
    size_t _earliestidx;

  }; //end of class ConvElecUtilities

} // end namespace mu2e
#endif
