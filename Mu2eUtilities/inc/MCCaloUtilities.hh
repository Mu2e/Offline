//
// $Id: MCCaloUtilities.hh,v 1.3 2011/05/17 15:41:36 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:36 $
//
// Original author Gianni Onorato
//

#ifndef Mu2eUtilities_MCCaloUtilities_hh
#define Mu2eUtilities_MCCaloUtilities_hh

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/ParameterSet/FileInPath.h"
#include "art/Utilities/Exception.h"
#include "art/Framework/Core/Event.h"

// Mu2e includes

#include "ToyDP/inc/SimParticleCollection.hh"

namespace mu2e {

  class MCCaloUtilities {

  public:

    MCCaloUtilities();

    ~MCCaloUtilities();

    void setTrackAndRO(const art::Event & event,
                       SimParticleCollection::key_type track,
                       uint32_t RO);

    void printOutCaloInfo();
    
    bool fromOutside();
    
    bool primary();
    
    bool generated();
    
    int startingVane();

    int getStartingVane(CLHEP::Hep3Vector origin);

    int localVane();

  private:

    uint32_t _localRO;
    uint32_t _localCrystal;
    uint32_t _localVane;
    int _startingVane;
    bool _fromOutside, _primary, _generated;

  };

}

#endif /* Mu2eUtilities_MCCaloUtilities_hh */
