//
// $Id: ConvElecUtilities.hh,v 1.6 2011/05/17 15:36:00 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:36:00 $
//
// Original author Gianni Onorato
//

#ifndef CONVELECUTILITIES_HH
#define CONVELECUTILITIES_HH

// C++ includes
#include <ostream>
#include <string>
#include <vector>
#include <map>
#include <cstdlib>
#include <memory>

// Framework includes
#include "art/Framework/Core/Event.h"

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
    ConvElecUtilities(const art::Event & event,
                      std::string g4ModuleLabel,
                      std::string trackerStepPoints);

    ~ConvElecUtilities();
    
  public:

    StepPointMC const& firstHit();
    StrawIndex earliestStrawIndex() const;
    const SimParticle& simConvElec() const;
    const ToyGenParticle& genConvElec();


    //Trivial accessors defined here

    //Number of Conversion Electrons. Used as a check
    int nConvElec() const {
      return _nconv;
    }
    
    //Returns how many hits come from the convElectron
    size_t hasStepPointMC() const {
      return _convElecHits.size();
    }
        
    //return a vector of index related to the stepPointMCCollection
    //identifying hits of the conversion electron
    const std::vector<size_t> & convElecHitsIdx() const {
      return _convElecHits;
    } //maybe it could be transformed in a vector of references to event hits 
    
    //return a vector of StrawIndex related to hits of the conversion electron
    const std::vector<StrawIndex> & convElecStrawIdx() const {
      return _convElecStrawIdx;
    }

  private:
    
    void checkConvElec(const art::Event & event);
    void lookAtHits(const art::Event & event);
    int _nconv;
    art::Handle<SimParticleCollection> _simParticles;
    art::Handle<ToyGenParticleCollection> _genParticles;
    art::Handle<StepPointMCCollection> hits;
    std::vector<size_t> _convElecHits; 
    std::vector<StrawIndex> _convElecStrawIdx;
    std::string _g4ModuleLabel, _trackerStepPoints;
    key_type _convTrackId;
    size_t _earliestidx;
    std::auto_ptr<SimParticle> _simParticle;


  }; //end of class ConvElecUtilities

} // end namespace mu2e
#endif
