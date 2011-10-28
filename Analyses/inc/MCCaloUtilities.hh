//
// $Id: MCCaloUtilities.hh,v 1.4 2011/10/28 18:47:05 greenc Exp $
// $Author: greenc $
// $Date: 2011/10/28 18:47:05 $
//
// Original author Gianni Onorato
//

#ifndef Mu2eUtilities_MCCaloUtilities_hh
#define Mu2eUtilities_MCCaloUtilities_hh

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Utilities/Exception.h"
#include "art/Framework/Principal/Event.h"

// Mu2e includes

#include "MCDataProducts/inc/SimParticleCollection.hh"

#include <string>

namespace mu2e {

  class MCCaloUtilities {

  public:

    MCCaloUtilities();

    ~MCCaloUtilities();

    void setTrackAndRO(const art::Event & event,
                       std::string const &_g4ModuleLabel,
                       SimParticleCollection::key_type track,
                       unsigned RO);

    void printOutCaloInfo();

    bool fromOutside();

    bool primary();

    bool generated();

    int startingVane();

    int getStartingVane(CLHEP::Hep3Vector origin);

    int localVane();

  private:

    unsigned _localRO;
    unsigned _localCrystal;
    unsigned _localVane;
    int _startingVane;
    bool _fromOutside, _primary, _generated;

  };

}

#endif /* Mu2eUtilities_MCCaloUtilities_hh */
