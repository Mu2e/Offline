//
// $Id: ConvElecUtilities.hh,v 1.5 2013/03/15 15:52:05 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/15 15:52:05 $
//
// Original author Gianni Onorato
//

#ifndef Mu2eG4_ConvElecUtilities_hh
#define Mu2eG4_ConvElecUtilities_hh

// C++ includes
#include <cstdlib>
#include <map>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

// Framework includes
#include "art/Framework/Principal/Event.h"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"

// Mu2e includes
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"

namespace mu2e {

  class ConvElecUtilities{

    typedef SimParticleCollection::key_type key_type;

  public:

    /**
     * Constructor.  Takes an event as an argument
     * plus g4module name (string) and the trackerStepPoints name (string)
     */
    ConvElecUtilities(const art::Event & event,
                      std::string const &generatorModuleLabel,
                      std::string const &g4ModuleLabel,
                      std::string const &trackerStepPoints,
                      std::string const &caloROModuleLabel);

    ~ConvElecUtilities();

  public:

    StepPointMC const& firstHit();
    StepPointMC const& firstCaloHit();
    StrawIndex earliestStrawIndex() const;
    const SimParticle& simConvElec() const;
    const GenParticle& genConvElec();

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

    double totEDep() const {
      return _totEDep;
    }

    bool gotCaloHit();

    private:



    std::string _generatorModuleLabel, _g4ModuleLabel, _trackerStepPoints, 
      _caloROlabel, _caloCrylabel;
    double _totEDep;
    void checkConvElec(const art::Event & event);
    void lookAtHits(const art::Event & event);
    void lookAtCalo(const art::Event & event);
    int _nconv;
    art::Handle<SimParticleCollection> _simParticles;
    art::Handle<GenParticleCollection> _genParticles;
    art::Handle<StepPointMCCollection> hits;
    art::Handle<StepPointMCCollection> calohits;
    std::vector<size_t> _convElecHits;
    std::vector<StrawIndex> _convElecStrawIdx;
    key_type _convTrackId;
    size_t _earliestidx, _earliestcry, _earliestvector;
    std::unique_ptr<SimParticle> _simParticle;
    bool _stepincalo;
    StepPointMC _earliestSPMC;

    //return a vector of index related to the stepPointMCCollection


  }; //end of class ConvElecUtilities

} // end namespace mu2e
#endif /* Mu2eG4_ConvElecUtilities_hh */
