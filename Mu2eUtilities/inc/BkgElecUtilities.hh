//
// $Id: BkgElecUtilities.hh,v 1.3 2013/03/15 15:52:05 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/15 15:52:05 $
//
// Original author Gianni Onorato
//

#ifndef Mu2eG4_BkgElecUtilities_hh
#define Mu2eG4_BkgElecUtilities_hh

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

  class BkgElecUtilities{

    typedef SimParticleCollection::key_type key_type;

  public:

    /**
     * Constructor.  Takes an event as an argument
     * plus g4module name (string) and the trackerStepPoints name (string)
     */
    BkgElecUtilities(const art::Event & event,
                      std::string const &generatorModuleLabel,
                      std::string const &g4ModuleLabel,
                      std::string const &trackerStepPoints,
                      std::string const &caloROModuleLabel);

    ~BkgElecUtilities();

  public:

    StepPointMC const& firstHit();
    StepPointMC const& firstCaloHit();
    StrawIndex earliestStrawIndex() const;
    const SimParticle& simBkgElec() const;
    const GenParticle& genBkgElec();

    //Trivial accessors defined here

    //Number of Conversion Electrons. Used as a check
    int nBkgElec() const {
      return _nbkg;
    }

    //Returns how many hits come from the bkgElectron
    size_t hasStepPointMC() const {
      return _bkgElecHits.size();
    }
   
    //return a vector of index related to the stepPointMCCollection
    //identifying hits of the bkg electron
    const std::vector<size_t> & bkgElecHitsIdx() const {
      return _bkgElecHits;
    } //maybe it could be transformed in a vector of references to event hits

    //return a vector of StrawIndex related to hits of the bkg electron
    const std::vector<StrawIndex> & bkgElecStrawIdx() const {
      return _bkgElecStrawIdx;
    }

    double totEDep() const {
      return _totEDep;
    }

    bool gotCaloHit();

    private:



    std::string _generatorModuleLabel, _g4ModuleLabel, _trackerStepPoints, 
      _caloROlabel, _caloCrylabel;
    double _totEDep;
    void checkBkgElec(const art::Event & event);
    void lookAtHits(const art::Event & event);
    void lookAtCalo(const art::Event & event);
    int _nbkg;
    art::Handle<SimParticleCollection> _simParticles;
    art::Handle<GenParticleCollection> _genParticles;
    art::Handle<StepPointMCCollection> hits;
    art::Handle<StepPointMCCollection> calohits;
    std::vector<size_t> _bkgElecHits;
    std::vector<StrawIndex> _bkgElecStrawIdx;
    key_type _bkgTrackId;
    size_t _earliestidx, _earliestcry, _earliestvector;
    std::unique_ptr<SimParticle> _simParticle;
    bool _stepincalo;
    StepPointMC _earliestSPMC;

    //return a vector of index related to the stepPointMCCollection


  }; //end of class BkgElecUtilities

} // end namespace mu2e
#endif /* Mu2eG4_BkgElecUtilities_hh */
